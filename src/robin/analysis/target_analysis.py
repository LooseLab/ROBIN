#!/usr/bin/env python3
"""
Target Analysis Module for robin.

Requires Python 3.12+. Automated target analysis; integrates with robin's
workflow system and processes files for target-specific analysis.

Classes: TargetMetadata, TargetAnalysis.
"""

from __future__ import annotations

import sys
if sys.version_info < (3, 12):
    raise RuntimeError("robin target_analysis requires Python 3.12 or newer")

import os
import tempfile
import logging
import time
import json
import subprocess
import shutil
import glob
import fcntl
import uuid
import re
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple
from dataclasses import dataclass
from io import StringIO
import numpy as np
import pandas as pd
import pysam
from robin.logging_config import get_job_logger
from robin.analysis.snp_processing import build_snp_display_data

# Optional import for Docker functionality
try:
    import docker
except ImportError:
    docker = None

# BAM flag: supplementary (0x800) | secondary (0x100) — alignments to skip for "primary-only" logic
_BAM_NON_PRIMARY_MASK = 0x800 | 0x100

# Fast Parquet I/O (matches fusion_work / accumulation pipelines)
_PARQUET_WRITE_KWARGS = {"index": False, "engine": "pyarrow", "compression": "snappy"}
_PARQUET_READ_KWARGS = {"engine": "pyarrow"}


def is_docker_available_for_snp_analysis() -> tuple[bool, str]:
    """
    Check if Docker is available for SNP analysis (Clair3 pipeline).
    Returns (True, "") if Docker is ready, (False, "error message") otherwise.
    """
    if docker is None:
        return False, "Docker Python package is not installed. Install it with: pip install docker"
    try:
        client = docker.from_env()
        client.ping()
        return True, ""
    except Exception as e:
        return False, f"Docker daemon is not available: {e}"


def json_serializable(obj):
    """Convert numpy/pandas types to JSON-serializable (3.12: type patterns)."""
    match obj:
        case int():
            return int(obj)
        case float():
            return float(obj)
        case _ if isinstance(obj, np.ndarray):
            return obj.tolist()
        case _ if isinstance(obj, pd.Series):
            return obj.tolist()
        case _ if isinstance(obj, pd.DataFrame):
            return obj.to_dict("records")
        case _:
            return obj


def _is_non_empty_file(path: Path) -> bool:
    """Return True when a path exists, is a file, and is non-empty."""
    return path.exists() and path.is_file() and path.stat().st_size > 0


def _resolve_clinvar_db_for_snpsift(logger: logging.Logger) -> Optional[str]:
    """
    Resolve ClinVar DB for SnpSift using bgzip + tabix index.

    Preferred and required format is `clinvar.vcf.gz` with a sibling `.tbi` index.
    If the index is missing, try to generate it in place.
    """
    try:
        # Use an alias to avoid shadowing module-level `resources`.
        from robin import resources as robin_resources

        resources_dir = Path(os.path.dirname(os.path.abspath(robin_resources.__file__)))
        clinvar_gz = resources_dir / "clinvar.vcf.gz"
        clinvar_tbi = resources_dir / "clinvar.vcf.gz.tbi"

        if not _is_non_empty_file(clinvar_gz):
            logger.warning(
                "ClinVar bgzipped VCF not found or empty at %s", clinvar_gz
            )
            return None

        if not _is_non_empty_file(clinvar_tbi):
            logger.info(
                "ClinVar tabix index missing at %s; attempting to generate it.",
                clinvar_tbi,
            )
            try:
                pysam.tabix_index(
                    str(clinvar_gz),
                    preset="vcf",
                    force=True,
                    keep_original=True,
                )
            except Exception as index_exc:
                logger.warning(
                    "Failed to create ClinVar tabix index for %s: %s",
                    clinvar_gz,
                    index_exc,
                )
                return None

        if not _is_non_empty_file(clinvar_tbi):
            logger.warning(
                "ClinVar tabix index still missing/empty after generation attempt: %s",
                clinvar_tbi,
            )
            return None

        logger.info(
            "Using ClinVar DB for SnpSift: %s (index: %s)",
            clinvar_gz,
            clinvar_tbi,
        )
        return str(clinvar_gz)
    except Exception as exc:
        logger.warning("Could not resolve ClinVar DB for SnpSift: %s", exc)
        return None


class FileLock:
    """
    Simple file-based lock for coordinating access across processes/threads.
    Uses fcntl for POSIX systems.
    """
    
    def __init__(self, lock_file: str, timeout: float = 30.0):
        self.lock_file = lock_file
        self.timeout = timeout
        self.fd = None
    
    def __enter__(self):
        self.acquire()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.release()
    
    def acquire(self):
        """Acquire the lock with timeout"""
        os.makedirs(os.path.dirname(self.lock_file), exist_ok=True)
        self.fd = open(self.lock_file, 'w')
        
        start_time = time.time()
        while True:
            try:
                fcntl.flock(self.fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                return
            except IOError:
                if time.time() - start_time > self.timeout:
                    raise TimeoutError(f"Could not acquire lock {self.lock_file} within {self.timeout}s")
                time.sleep(0.1)
    
    def release(self):
        """Release the lock"""
        if self.fd:
            fcntl.flock(self.fd, fcntl.LOCK_UN)
            self.fd.close()
            self.fd = None


# Import robin resources for BED files
try:
    from robin import resources
except ImportError:
    resources = None


def _bam_has_any_alignment(bam_path: str) -> bool:
    """
    True if the BAM contains at least one alignment record (stops after first read).
    Faster than AlignmentFile.count(until_eof=True) for emptiness checks.
    """
    try:
        if not os.path.exists(bam_path) or os.path.getsize(bam_path) == 0:
            return False
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for _ in bam.fetch(until_eof=True):
                return True
    except OSError:
        return False
    return False


def _load_bed_regions(bedfile: str) -> List[Tuple[str, int, int]]:
    """Load BED regions into a list of (chrom, start, end) tuples."""
    regions: List[Tuple[str, int, int]] = []
    with open(bedfile, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                regions.append((chrom, start, end))
    return regions


def run_bedtools(bamfile, bedfile, tempbamfile, regions: Optional[List[Tuple[str, int, int]]] = None):
    """
    Extract target regions from BAM file, keeping all mappings (primary, secondary, supplementary)
    for reads that overlap the target regions.
    
    This function uses a two-step process:
    1. Extract read names from primary alignments that overlap target regions
    2. Extract ALL alignments (primary, secondary, supplementary) for those read names
    
    Parameters
    ----------
    bamfile : str
        Path to the input BAM file
    bedfile : str
        Path to the BED file defining regions
    tempbamfile : str
        Path where the output BAM file should be written
    regions : list of tuples, optional
        Pre-loaded BED regions to avoid repeated parsing
    """
    logger = logging.getLogger("robin.target")

    try:
        # Step 1: Extract read names from primary alignments overlapping target regions
        # Use bedtools to find primary alignments (filtering out supplementary and secondary)
        # Flag 2304 = 0x800 (supplementary) | 0x100 (secondary)
        # We use -F 2304 to exclude supplementary and secondary, keeping only primary alignments
        logger.debug(f"Step 1: Extracting read names from primary alignments overlapping {bedfile}")
        
        # Read BED file to get regions (unless preloaded)
        if regions is None:
            regions = _load_bed_regions(bedfile)
        
        if not regions:
            logger.warning(f"No valid regions found in BED file: {bedfile}")
            # Create empty BAM file
            with pysam.AlignmentFile(bamfile, "rb") as in_bam:
                header = in_bam.header.copy()
                with pysam.AlignmentFile(tempbamfile, "wb", header=header) as out_bam:
                    pass
            pysam.index(tempbamfile)
            return
        
        # Open input BAM and collect read names from primary alignments overlapping regions
        read_names = set()
        with pysam.AlignmentFile(bamfile, "rb") as in_bam:
            # Require BAM index for performance and correctness
            index_file = f"{bamfile}.bai"
            if not os.path.exists(index_file):
                raise FileNotFoundError(f"BAM index (.bai) not found for {bamfile}")
            
            # Use indexed access (faster)
            for chrom, start, end in regions:
                try:
                    # Fetch reads overlapping this region
                    for read in in_bam.fetch(chrom, start, end):
                        # Only primary alignments (not supplementary or secondary)
                        if (read.flag & _BAM_NON_PRIMARY_MASK) == 0:
                            read_names.add(read.query_name)
                except ValueError:
                    # Chromosome not found in BAM, skip
                    logger.debug(f"Chromosome {chrom} not found in BAM file, skipping")
                    continue
        
        logger.debug(f"Found {len(read_names)} unique read names overlapping target regions")
        
        if not read_names:
            logger.warning("No reads found overlapping target regions")
            # Create empty BAM file with same header
            with pysam.AlignmentFile(bamfile, "rb") as in_bam:
                header = in_bam.header.copy()
                with pysam.AlignmentFile(tempbamfile, "wb", header=header) as out_bam:
                    pass
            pysam.index(tempbamfile)
            return
        
        # Step 2: Extract ALL alignments (primary, secondary, supplementary) for those read names
        logger.debug(f"Step 2: Extracting all alignments for {len(read_names)} read names")
        
        reads_written = 0
        names = read_names
        with pysam.AlignmentFile(bamfile, "rb") as in_bam:
            header = in_bam.header.copy()
            with pysam.AlignmentFile(tempbamfile, "wb", header=header) as out_bam:
                # Single pass: emit all alignments for reads that hit target regions in step 1
                for read in in_bam.fetch(until_eof=True):
                    if read.query_name in names:
                        out_bam.write(read)
                        reads_written += 1
                
                # Ensure all data is written to disk before closing
                out_bam.flush()
                logger.debug(f"Wrote {reads_written} alignments (including secondary/supplementary) to output BAM")
        
        # Verify the file was written successfully before indexing
        if not os.path.exists(tempbamfile):
            raise RuntimeError(f"Output BAM file was not created: {tempbamfile}")
        
        file_size = os.path.getsize(tempbamfile)
        if file_size == 0:
            logger.warning(f"Output BAM file is empty: {tempbamfile}")
        else:
            logger.debug(f"Output BAM file size: {file_size} bytes")
        
        # Index the output BAM (only if file has content)
        if file_size > 0:
            try:
                pysam.index(tempbamfile)
                logger.info(f"Successfully extracted target regions to {tempbamfile} ({reads_written} alignments)")
            except Exception as e:
                logger.error(f"Failed to index BAM file {tempbamfile}: {e}")
                # Try to verify if the BAM file is valid
                try:
                    with pysam.AlignmentFile(tempbamfile, "rb") as test_bam:
                        test_count = test_bam.count(until_eof=True)
                        logger.info(f"BAM file is readable, contains {test_count} reads")
                except Exception as verify_error:
                    logger.error(f"BAM file appears corrupted: {verify_error}")
                    raise
        else:
            logger.warning(f"Skipping indexing for empty BAM file: {tempbamfile}")
        
    except Exception as e:
        logger.error(f"Error in run_bedtools: {e}")
        import traceback
        logger.error(traceback.format_exc())


def get_covdfs(bamfile, bedfile=None):
    """
    Extract coverage information from a BAM file.

    This function calculates coverage statistics for both the entire genome
    and specific target regions defined in a BED file.

    Parameters
    ----------
    bamfile : str
        Path to the input BAM file.
    bedfile : str, optional
        Path to the target BED file. If not provided, will use unique_genes.bed.

    Returns
    -------
    tuple
        A tuple containing two pandas DataFrames:
        - Coverage statistics for the entire genome
        - Coverage statistics for target regions

    Notes
    -----
    Uses pysam for efficient BAM file processing.
    """
    logger = logging.getLogger("robin.target")

    try:
        # Get genome-wide coverage using pysam.coverage
        coverage_output = pysam.coverage(f"{bamfile}")

        newcovdf = pd.read_csv(StringIO(coverage_output), sep="\t")

        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(f"Raw pysam.coverage columns: {list(newcovdf.columns)}")
            logger.debug(
                "Sample raw coverage data: %s",
                newcovdf.head(2).to_dict("records"),
            )

        newcovdf.drop(
            columns=["coverage", "meanbaseq", "meanmapq"],
            inplace=True,
        )

        if logger.isEnabledFor(logging.DEBUG):
            logger.debug(f"After dropping columns: {list(newcovdf.columns)}")

        # Find target BED file - use provided bedfile or fallback to unique_genes.bed
        target_bed = bedfile
        if target_bed is None:
            # Fallback to unique_genes.bed for backward compatibility
            if resources is not None:
                try:
                    target_bed = os.path.join(
                        os.path.dirname(os.path.abspath(resources.__file__)),
                        "unique_genes.bed",
                    )
                    if not os.path.exists(target_bed):
                        target_bed = None
                except Exception:
                    pass

            # Fallback paths for unique_genes.bed
            if target_bed is None:
                possible_paths = [
                    "unique_genes.bed",
                    "data/unique_genes.bed",
                    "/usr/local/share/unique_genes.bed",
                ]
                for path in possible_paths:
                    if os.path.exists(path):
                        target_bed = path
                        break

        if target_bed is None:
            logger.warning("Target BED file not found, skipping bedcov analysis")
            bedcovdf = pd.DataFrame(
                columns=["chrom", "startpos", "endpos", "name", "bases"]
            )
        else:
            # Get target region coverage using pysam.bedcov
            bedcovdf = pd.read_csv(
                StringIO(
                    pysam.bedcov(
                        target_bed,
                        f"{bamfile}",
                    )
                ),
                names=["chrom", "startpos", "endpos", "name", "bases"],
                sep="\t",
            )

        logger.info(f"Successfully extracted coverage data from {bamfile}")
        logger.info(f"Genome coverage: {len(newcovdf)} regions")
        logger.info(f"Target coverage: {len(bedcovdf)} regions")

        return newcovdf, bedcovdf

    except Exception as e:
        logger.error(f"Error in get_covdfs: {str(e)}")
        return None, None


def get_read_counts_per_target(bamfile, bedfile):
    """
    Count reads overlapping each target region in a BED file.
    
    Parameters
    ----------
    bamfile : str
        Path to the input BAM file
    bedfile : str
        Path to the BED file defining target regions
    
    Returns
    -------
    pd.DataFrame
        DataFrame with columns: chrom, startpos, endpos, name, reads
        Returns empty DataFrame if extraction fails
    """
    logger = logging.getLogger("robin.target")
    
    try:
        # Read BED file to get regions
        bed_regions = []
        with open(bedfile, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 4:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    name = parts[3] if parts[3].strip() else f"{chrom}:{start}-{end}"
                    bed_regions.append((chrom, start, end, name))
        
        if not bed_regions:
            logger.warning(f"No valid regions found in BED file: {bedfile}")
            return pd.DataFrame(columns=["chrom", "startpos", "endpos", "name", "reads"])
        
        # Count reads per region
        read_counts = []
        with pysam.AlignmentFile(bamfile, "rb") as bam:
            index_file = f"{bamfile}.bai"
            if not os.path.exists(index_file):
                logger.debug(
                    "BAM index (.bai) not found for %s; pysam may work if BAM is coordinate-sorted",
                    bamfile,
                )
            for chrom, start, end, name in bed_regions:
                try:
                    read_count = 0
                    for read in bam.fetch(chrom, start, end):
                        if (read.flag & _BAM_NON_PRIMARY_MASK) == 0:
                            read_count += 1
                    
                    read_counts.append({
                        'chrom': chrom,
                        'startpos': start,
                        'endpos': end,
                        'name': name,
                        'reads': read_count
                    })
                except ValueError:
                    logger.debug(f"Chromosome {chrom} not found in BAM file, skipping")
                    read_counts.append({
                        'chrom': chrom,
                        'startpos': start,
                        'endpos': end,
                        'name': name,
                        'reads': 0
                    })
                    continue
        
        df = pd.DataFrame(read_counts)
        logger.debug(f"Extracted read counts for {len(df)} target regions from {bamfile}")
        return df
        
    except Exception as e:
        logger.error(f"Error extracting read counts per target: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return pd.DataFrame(columns=["chrom", "startpos", "endpos", "name", "reads"])


def run_bedmerge(newcovdf, cov_df_main, bedcovdf, bedcov_df_main):
    """
    Merge coverage dataframes for incremental processing.

    Parameters
    ----------
    newcovdf : pd.DataFrame
        New genome coverage data
    cov_df_main : pd.DataFrame
        Existing genome coverage data
    bedcovdf : pd.DataFrame
        New target coverage data
    bedcov_df_main : pd.DataFrame
        Existing target coverage data

    Returns
    -------
    tuple
        Merged coverage dataframes
    """
    logger = logging.getLogger("robin.target")

    logger.debug("Starting coverage merge:")
    logger.debug(
        f"  New genome coverage: {len(newcovdf) if newcovdf is not None else 0} regions"
    )
    logger.debug(
        f"  Existing genome coverage: {len(cov_df_main) if cov_df_main is not None else 0} regions"
    )
    logger.debug(
        f"  New target coverage: {len(bedcovdf) if bedcovdf is not None else 0} targets"
    )
    logger.debug(
        f"  Existing target coverage: {len(bedcov_df_main) if bedcov_df_main is not None else 0} targets"
    )

    # Merge genome coverage data
    if cov_df_main is not None and not cov_df_main.empty:
        merged_df = pd.merge(
            newcovdf,
            cov_df_main,
            on=["#rname", "startpos", "endpos"],
            suffixes=("_df1", "_df2"),
        )
        merged_df["numreads"] = merged_df["numreads_df1"] + merged_df["numreads_df2"]
        merged_df["covbases"] = merged_df["covbases_df1"] + merged_df["covbases_df2"]
        merged_df["meandepth"] = merged_df["meandepth_df1"] + merged_df["meandepth_df2"]

        merged_df.drop(
            columns=[
                "numreads_df1",
                "numreads_df2",
                "meandepth_df1",
                "meandepth_df2",
                "covbases_df1",
                "covbases_df2",
            ],
            inplace=True,
        )

        logger.debug(
            f"Merged genome coverage: {len(cov_df_main)} + {len(newcovdf)} -> {len(merged_df)} regions"
        )
        logger.debug(
            f"Sample merged genome data: {merged_df.head(3).to_dict('records')}"
        )
    else:
        merged_df = newcovdf
        logger.debug(
            f"No existing genome coverage, using new data: {len(merged_df)} regions"
        )

    # Merge target coverage data
    if bedcov_df_main is not None and not bedcov_df_main.empty:
        merged_bed_df = pd.merge(
            bedcovdf,
            bedcov_df_main,
            on=["chrom", "startpos", "endpos", "name"],
            suffixes=("_df1", "_df2"),
        )
        merged_bed_df["bases"] = merged_bed_df["bases_df1"] + merged_bed_df["bases_df2"]
        merged_bed_df.drop(columns=["bases_df1", "bases_df2"], inplace=True)

        logger.debug(
            f"Merged target coverage: {len(bedcov_df_main)} + {len(bedcovdf)} -> {len(merged_bed_df)} targets"
        )
        logger.debug(
            f"Sample merged target data: {merged_bed_df.head(3).to_dict('records')}"
        )
    else:
        merged_bed_df = bedcovdf
        logger.debug(
            f"No existing target coverage, using new data: {len(merged_bed_df)} targets"
        )

    return merged_df, merged_bed_df


@dataclass(slots=True)
class TargetMetadata:
    """Container for target analysis metadata and results."""

    sample_id: str
    file_path: str
    analysis_timestamp: float
    target_data_path: Optional[str] = None
    target_plot_path: Optional[str] = None
    analysis_results: Optional[Dict] = None
    processing_steps: List[str] = None
    error_message: Optional[str] = None
    coverage_data: Optional[Dict] = None
    target_bam_path: Optional[str] = None
    coverage_over_time: Optional[np.ndarray] = None

    def __post_init__(self):
        if self.processing_steps is None:
            self.processing_steps = []
        if self.analysis_results is None:
            self.analysis_results = {}
        if self.coverage_data is None:
            self.coverage_data = {}


class TargetAnalysis:
    """Target analysis worker"""

    def __init__(self, work_dir=None, config_path=None, threads=4, target_panel=None, 
                 batch_size=10, use_staging=True):
        logger = logging.getLogger("robin.target")

        self.work_dir = work_dir or os.getcwd()
        self.config_path = config_path
        self.threads = threads
        self.target_panel = target_panel
        self.batch_size = batch_size  # Trigger accumulation after N files
        self.use_staging = use_staging  # Enable/disable staging optimization

        # Initialize counter for incremental file naming
        self.file_counter = 1

        # Load configuration if provided
        self.config = self._load_config()

        # Note: State is now persisted to files instead of in-memory dictionaries

        # Find BED file for target extraction
        self.bedfile = self._find_target_bed(target_panel)

        # Configuration parameters
        self.callthreshold = self.config.get("call_threshold", 0.1)
        self.simtime = self.config.get("simtime", False)
        self.reference = self.config.get("reference", None)
        self.snp_calling = self.config.get("snp_calling", False)

        # Check environment variable for reference genome
        if not self.reference:
            self.reference = os.environ.get("robin_REFERENCE")

        # Log reference genome status
        if self.reference:
            logger.info(f"Reference genome configured: {self.reference}")
        else:
            logger.info(
                "No reference genome configured - SNP calling will not be available"
            )

        logger.info(f"Target Analysis initialized (staging={'enabled' if use_staging else 'disabled'}, batch_size={batch_size})")

    def _get_master_bed_path(self, sample_id: str) -> Optional[str]:
        """
        Get the path to the master BED file for a sample if it exists.
        Master BED includes target panel + CNV breakpoints + fusion breakpoints + master BED breakpoints.
        
        Args:
            sample_id: Sample ID
            
        Returns:
            Path to master BED file, or None if not found
        """
        logger = logging.getLogger("robin.target")
        try:
            import glob
            from robin.analysis.fusion_work import _load_analysis_counter
            
            sample_dir = os.path.join(self.work_dir, sample_id)
            bed_dir = os.path.join(sample_dir, "bed_files")
            
            if not os.path.exists(bed_dir):
                return None
            
            # Try to get analysis counter
            analysis_counter = _load_analysis_counter(sample_id, self.work_dir)
            master_bed_path = os.path.join(bed_dir, f"master_{analysis_counter:03d}.bed")
            
            if os.path.exists(master_bed_path):
                return master_bed_path
            
            # Try to find the latest master BED file if counter-based doesn't exist
            master_bed_files = glob.glob(os.path.join(bed_dir, "master_*.bed"))
            if master_bed_files:
                # Sort by modification time and return the latest
                latest = max(master_bed_files, key=os.path.getmtime)
                logger.debug(f"Using latest master BED file: {latest}")
                return latest
        except Exception as e:
            logger.debug(f"Error finding master BED file: {e}")
        
        return None

    def _find_target_bed(self, target_panel: str) -> str:
        """Find the target BED file from robin resources based on panel type"""
        logger = logging.getLogger("robin.target")
        
        # Determine the correct BED file name based on panel type
        bed_filename = None
        if target_panel == "rCNS2":
            bed_filename = "rCNS2_panel_name_uniq.bed"
        elif target_panel == "AML":
            bed_filename = "AML_panel_name_uniq.bed"
        else:
            # Check for custom panel
            bed_filename = f"{target_panel}_panel_name_uniq.bed"
            logger.info(f"Using custom panel: {target_panel}")
        
        if resources is not None:
            try:
                bed_path = os.path.join(
                    os.path.dirname(os.path.abspath(resources.__file__)),
                    bed_filename,
                )
                if os.path.exists(bed_path):
                    logger.info(f"Found {target_panel} BED file: {bed_path}")
                    return bed_path
            except Exception:
                pass

        # Fallback paths
        possible_paths = [
            bed_filename,
            f"data/{bed_filename}",
            f"/usr/local/share/{bed_filename}",
        ]

        for path in possible_paths:
            if os.path.exists(path):
                logger.info(f"Found {target_panel} BED file: {path}")
                return path

        # If not found, create a placeholder (this will cause an error later)
        logger.warning(f"Target BED file '{bed_filename}' not found for panel '{target_panel}', will use placeholder")
        return bed_filename

    def _load_config(self) -> Dict[str, Any]:
        """Load configuration from file or use defaults"""
        if self.config_path and os.path.exists(self.config_path):
            try:
                with open(self.config_path, "r") as f:
                    return json.load(f)
            except Exception as e:
                logging.warning(f"Error loading config from {self.config_path}: {e}")

        # Default configuration
        return {
            "analysis_type": "standard",
            "output_format": "json",
            "include_plots": True,
            "threshold": 0.5,
            "call_threshold": 0.1,
            "simtime": False,
            "snp_calling": False,
        }

    def _get_next_file_number(self) -> int:
        """Get the next file number for incremental processing"""
        return self.file_counter

    def _check_and_create_folder(self, base_dir: str, sample_id: str) -> str:
        """Create sample-specific directory and return path"""
        sample_dir = os.path.join(base_dir, sample_id)
        os.makedirs(sample_dir, exist_ok=True)
        return sample_dir
    
    def _get_staging_dir(self, sample_id: str) -> str:
        """Get staging directory for temporary per-file results"""
        staging_dir = os.path.join(self.work_dir, sample_id, "_staging")
        os.makedirs(staging_dir, exist_ok=True)
        return staging_dir
    
    def _get_lock_file(self, sample_id: str, lock_type: str = "counter") -> str:
        """Get lock file path for coordinating concurrent access"""
        lock_dir = os.path.join(self.work_dir, sample_id, "_locks")
        os.makedirs(lock_dir, exist_ok=True)
        return os.path.join(lock_dir, f"{lock_type}.lock")
    
    def _get_pending_count(self, sample_id: str) -> int:
        """Get count of complete staging sets pending accumulation (thread-safe)."""
        staging_dir = self._get_staging_dir(sample_id)
        complete_sets, _, _ = self._collect_complete_staging_sets(staging_dir)
        return len(complete_sets)

    def _collect_complete_staging_sets(
        self, staging_dir: str
    ) -> Tuple[List[Tuple[str, str, str, str]], Dict[str, int], List[str]]:
        """
        Build aligned staging tuples using the shared numeric suffix.

        Returns:
            - complete_sets: list of (coverage, bedcov, timestamp, source_bam)
            - counts: raw file counts per staging type
            - incomplete_ids: staged suffix IDs that are missing one or more companion files
        """
        pattern_map = {
            "coverage": ("coverage", "parquet"),
            "bedcov": ("bedcov", "parquet"),
            "timestamp": ("timestamp", "txt"),
            "source_bam": ("source_bam", "txt"),
        }

        indexed_paths: Dict[str, Dict[str, str]] = {k: {} for k in pattern_map}
        counts: Dict[str, int] = {}

        for key, (prefix, ext) in pattern_map.items():
            files = sorted(glob.glob(os.path.join(staging_dir, f"{prefix}_*.{ext}")))
            counts[key] = len(files)
            suffix_re = re.compile(rf"^{re.escape(prefix)}_(\d+)\.{re.escape(ext)}$")
            for file_path in files:
                name = os.path.basename(file_path)
                match = suffix_re.match(name)
                if not match:
                    continue
                indexed_paths[key][match.group(1)] = file_path

        id_sets = [set(v.keys()) for v in indexed_paths.values()]
        shared_ids = sorted(set.intersection(*id_sets)) if id_sets else []
        all_ids = set.union(*id_sets) if id_sets else set()
        incomplete_ids = sorted(all_ids - set(shared_ids))

        complete_sets = [
            (
                indexed_paths["coverage"][suffix_id],
                indexed_paths["bedcov"][suffix_id],
                indexed_paths["timestamp"][suffix_id],
                indexed_paths["source_bam"][suffix_id],
            )
            for suffix_id in shared_ids
        ]

        return complete_sets, counts, incomplete_ids
    
    def _atomic_counter_increment(self, sample_id: str) -> int:
        """
        Atomically increment and return the file counter for a sample.
        Uses file locking to prevent race conditions.
        
        Returns:
            The counter value to use for this file
        """
        lock_file = self._get_lock_file(sample_id, "counter")
        counter_file = os.path.join(self.work_dir, sample_id, "target_analysis_counter.txt")
        
        with FileLock(lock_file, timeout=30.0):
            # Read current counter
            if os.path.exists(counter_file):
                try:
                    with open(counter_file, "r") as f:
                        counter = int(f.read().strip())
                except (ValueError, IOError):
                    counter = 0
            else:
                counter = 0
            
            # Write incremented counter
            os.makedirs(os.path.dirname(counter_file), exist_ok=True)
            with open(counter_file, "w") as f:
                f.write(str(counter + 1))
            
            return counter
    
    def process_file_with_staging(
        self,
        file_path: str,
        metadata: Dict[str, Any],
        timestamp: Optional[float] = None,
    ) -> Tuple[TargetMetadata, bool]:
        """
        Fast per-file processing that saves results to staging area.
        Does NOT merge with accumulated data - much faster for large datasets.
        
        Args:
            file_path: Path to the input file
            metadata: File metadata from preprocessing
            timestamp: Optional timestamp for coverage tracking
        
        Returns:
            Tuple of (TargetMetadata, should_accumulate)
            - TargetMetadata: Results from this file
            - should_accumulate: True if batch accumulation should run now
        """
        logger = logging.getLogger("robin.target")
        
        logger.info(f"Processing file with staging: {file_path}")
        start_time = time.time()
        
        # Extract sample ID from metadata
        sample_id = metadata.get("sample_id", "unknown")
        logger.debug(f"Extracted sample_id: {sample_id}")
        
        target_result = TargetMetadata(
            sample_id=sample_id, file_path=file_path, analysis_timestamp=start_time
        )
        
        try:
            # Step 1: Validate input file
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"Input file not found: {file_path}")
            
            target_result.processing_steps.append("file_validation")
            
            # Step 2: Create sample-specific output directory
            sample_output_dir = self._check_and_create_folder(self.work_dir, sample_id)
            logger.info(f"Sample output directory: {sample_output_dir}")
            
            # Step 3: Get atomic counter (thread-safe)
            analysis_counter = self._atomic_counter_increment(sample_id)
            logger.info(f"Assigned file counter: {analysis_counter}")
            
            target_result.processing_steps.append("counter_assigned")
            
            # Step 4: Extract coverage data (no loading of accumulated data)
            logger.info("Extracting coverage data...")
            newcovdf, bedcovdf = get_covdfs(file_path, self.bedfile)
            
            if newcovdf is None or bedcovdf is None:
                raise RuntimeError("Failed to extract coverage data from BAM file")
            
            target_result.processing_steps.append("coverage_extracted")
            logger.info(
                f"Coverage data extracted: genome={newcovdf.shape}, targets={bedcovdf.shape}"
            )
            
            # Step 5: Save to staging (Parquet is ~5-10x faster than CSV)
            staging_dir = self._get_staging_dir(sample_id)
            
            coverage_staging = os.path.join(
                staging_dir, f"coverage_{analysis_counter:06d}.parquet"
            )
            bedcov_staging = os.path.join(
                staging_dir, f"bedcov_{analysis_counter:06d}.parquet"
            )
            timestamp_staging = os.path.join(
                staging_dir, f"timestamp_{analysis_counter:06d}.txt"
            )
            
            # Save coverage data to staging
            newcovdf.to_parquet(coverage_staging, **_PARQUET_WRITE_KWARGS)
            bedcovdf.to_parquet(bedcov_staging, **_PARQUET_WRITE_KWARGS)
            
            # Save timestamp for coverage tracking
            current_timestamp = timestamp * 1000 if (self.simtime and timestamp) else time.time() * 1000
            with open(timestamp_staging, "w") as f:
                f.write(str(current_timestamp))
            
            # Save source BAM file path for target.bam creation during accumulation
            source_bam_staging = os.path.join(
                staging_dir, f"source_bam_{analysis_counter:06d}.txt"
            )
            with open(source_bam_staging, "w") as f:
                f.write(file_path)
            
            target_result.processing_steps.append("saved_to_staging")
            logger.info(f"Saved to staging: {coverage_staging}")
            
            # Step 6: Store minimal coverage data in metadata (for logging)
            target_result.coverage_data = {
                "genome_coverage_shape": newcovdf.shape,
                "target_coverage_shape": bedcovdf.shape,
                "staging_file": coverage_staging,
            }
            
            # Step 7: Check if accumulation should run
            pending_count = self._get_pending_count(sample_id)
            should_accumulate = pending_count >= self.batch_size
            
            logger.info(
                f"File staged successfully. Pending files: {pending_count}/{self.batch_size}"
            )
            
            if should_accumulate:
                logger.info(
                    f"Accumulation threshold reached ({pending_count} >= {self.batch_size})"
                )
            
            target_result.processing_steps.append("staging_complete")
            elapsed = time.time() - start_time
            logger.info(
                f"Staging complete for {sample_id} in {elapsed:.2f}s (vs ~{elapsed*10:.1f}s without staging)"
            )
            
            return target_result, should_accumulate
            
        except Exception as e:
            error_details = f"Error in staging for {sample_id}: {str(e)}"
            logger.error(error_details)
            target_result.error_message = error_details
            target_result.processing_steps.append("staging_failed")
            return target_result, False

    def process_file(
        self,
        file_path: str,
        metadata: Dict[str, Any],
        timestamp: Optional[float] = None,
    ) -> TargetMetadata:
        """
        Process a file for target analysis.

        Args:
            file_path: Path to the input file
            metadata: File metadata from preprocessing
            timestamp: Optional timestamp for coverage tracking

        Returns:
            TargetMetadata object with analysis results
        """
        logger = logging.getLogger("robin.target")

        logger.info(f"Processing file: {file_path} for target analysis")
        start_time = time.time()

        # Extract sample ID from metadata
        sample_id = metadata.get("sample_id", "unknown")
        logger.debug(f"Extracted sample_id: {sample_id}")

        target_result = TargetMetadata(
            sample_id=sample_id, file_path=file_path, analysis_timestamp=start_time
        )
        logger.debug(f"Created TargetMetadata object: {target_result}")

        logger.info(f"Starting target analysis for sample: {sample_id}")
        logger.debug(f"File: {os.path.basename(file_path)}")
        logger.debug(f"Work directory: {self.work_dir}")

        try:
            # Step 1: Validate input file
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"Input file not found: {file_path}")

            target_result.processing_steps.append("file_validation")

            # Step 2: Create sample-specific output directory
            sample_output_dir = self._check_and_create_folder(self.work_dir, sample_id)
            logger.info(f"Created sample output directory: {sample_output_dir}")

            target_result.processing_steps.append("directory_created")

            # Step 3: Load analysis counter from disk
            analysis_counter = self._load_analysis_counter(
                sample_id, self.work_dir, logger
            )
            logger.info(f"Loaded analysis counter for {sample_id}: {analysis_counter}")

            # Step 4: Load existing accumulated data from files
            existing_covdf = self._load_existing_coverage_data(
                sample_output_dir, "coverage_main.csv", logger
            )
            existing_bedcovdf = self._load_existing_coverage_data(
                sample_output_dir, "bed_coverage_main.csv", logger
            )
            existing_coverage_over_time = self._load_existing_coverage_over_time(
                sample_output_dir, logger
            )
            existing_target_bam_path = os.path.join(sample_output_dir, "target.bam")

            # Step 5: Create temporary directory for processing
            with tempfile.TemporaryDirectory(dir=sample_output_dir) as temp_dir:
                logger.info(f"Created temporary directory: {temp_dir}")

                # Step 6: Extract coverage data using get_covdfs
                logger.info("Extracting coverage data...")
                newcovdf, bedcovdf = get_covdfs(file_path, self.bedfile)

                if newcovdf is None or bedcovdf is None:
                    raise RuntimeError("Failed to extract coverage data from BAM file")

                # Store coverage data in metadata
                target_result.coverage_data = {
                    "genome_coverage_shape": newcovdf.shape,
                    "target_coverage_shape": bedcovdf.shape,
                }

                target_result.processing_steps.append("coverage_extracted")
                logger.info(
                    f"Coverage data extracted: genome={newcovdf.shape}, targets={bedcovdf.shape}"
                )

                # Step 7: Extract target regions using bedtools
                logger.info("Extracting target regions using bedtools...")
                tempbamfile = tempfile.NamedTemporaryFile(
                    dir=sample_output_dir, suffix=".bam"
                )

                # Prefer master BED file if available, otherwise use original target panel BED
                # Master BED includes target panel + CNV + fusion + master BED breakpoints
                sample_id = metadata.get("sample_id", "unknown")
                targets_bed = self._get_master_bed_path(sample_id) or self.bedfile
                
                # Run bedtools intersection
                run_bedtools(file_path, targets_bed, tempbamfile.name)

                # Check if target BAM has reads (first record only, not full count)
                if _bam_has_any_alignment(tempbamfile.name):
                    logger.info(f"Target regions found in {sample_id}")

                    # Handle target BAM file accumulation
                    if os.path.exists(existing_target_bam_path):
                        # Merge with existing target BAM
                        tempbamholder = tempfile.NamedTemporaryFile(
                            dir=sample_output_dir, suffix=".bam"
                        )
                        pysam.cat(
                            "-o",
                            tempbamholder.name,
                            existing_target_bam_path,
                            tempbamfile.name,
                        )
                        shutil.copy2(tempbamholder.name, existing_target_bam_path)
                        target_result.target_bam_path = existing_target_bam_path
                        logger.info(
                            f"Updated target BAM file: {existing_target_bam_path}"
                        )

                        # Clean up temporary files
                        try:
                            os.remove(f"{tempbamholder.name}.bai")
                        except FileNotFoundError:
                            pass
                    else:
                        # Create new target BAM file
                        shutil.copy2(tempbamfile.name, existing_target_bam_path)
                        target_result.target_bam_path = existing_target_bam_path
                        logger.info(
                            f"Created target BAM file: {existing_target_bam_path}"
                        )
                else:
                    logger.info(f"No target regions found in {sample_id}")

                # Clean up temporary BAM index
                try:
                    os.remove(f"{tempbamfile.name}.bai")
                except FileNotFoundError:
                    pass

                target_result.processing_steps.append("target_regions_extracted")

                # Step 8: Update coverage data
                if existing_covdf is not None:
                    updated_covdf, updated_bedcovdf = run_bedmerge(
                        newcovdf, existing_covdf, bedcovdf, existing_bedcovdf
                    )
                else:
                    updated_covdf, updated_bedcovdf = newcovdf, bedcovdf

                target_result.processing_steps.append("coverage_updated")

                # Step 9: Calculate coverage statistics
                bases = updated_covdf["covbases"].sum()
                genome = updated_covdf["endpos"].sum()
                coverage = bases / genome if genome > 0 else 0.0

                logger.info(
                    f"Coverage calculated: {coverage:.4f} ({bases} bases / {genome} genome)"
                )
                logger.debug(f"Coverage data shape: {updated_covdf.shape}")
                logger.debug(f"Coverage columns: {list(updated_covdf.columns)}")
                logger.debug("Sample coverage data:")
                logger.debug(
                    f"  First few rows: {updated_covdf.head(3).to_dict('records')}"
                )
                logger.debug(f"  covbases sum: {bases}")
                logger.debug(f"  endpos sum: {genome}")

                # Step 10: Track coverage over time
                if self.simtime and timestamp:
                    currenttime = timestamp * 1000
                else:
                    currenttime = time.time() * 1000

                if existing_coverage_over_time is not None:
                    updated_coverage_over_time = np.vstack(
                        [existing_coverage_over_time, [(currenttime, coverage)]]
                    )
                else:
                    updated_coverage_over_time = np.array([[currenttime, coverage]])

                target_result.coverage_over_time = updated_coverage_over_time
                target_result.processing_steps.append("coverage_tracked")

                # Step 11: Save updated coverage data to files
                np.save(
                    os.path.join(sample_output_dir, "coverage_time_chart.npy"),
                    updated_coverage_over_time,
                )

                updated_covdf.to_csv(
                    os.path.join(sample_output_dir, "coverage_main.csv"),
                    index=False,
                )

                # Calculate length and coverage for bed_coverage_main.csv (matching reference format)
                bed_coverage_main_df = updated_bedcovdf.copy()
                bed_coverage_main_df["length"] = (
                    bed_coverage_main_df["endpos"]
                    - bed_coverage_main_df["startpos"]
                    + 1
                )

                bed_coverage_main_df["coverage"] = (
                    bed_coverage_main_df["bases"] / bed_coverage_main_df["length"]
                )

                # Reorder columns to match reference format: chrom,startpos,endpos,name,length,coverage,bases
                bed_coverage_main_df = bed_coverage_main_df[
                    [
                        "chrom",
                        "startpos",
                        "endpos",
                        "name",
                        "length",
                        "coverage",
                        "bases",
                    ]
                ]

                bed_coverage_main_df.to_csv(
                    os.path.join(sample_output_dir, "bed_coverage_main.csv"),
                    index=False,
                )

                target_result.processing_steps.append("coverage_saved")

                # Step 12: Calculate target coverage
                target_coverage_df = updated_bedcovdf.copy()
                target_coverage_df["length"] = (
                    target_coverage_df["endpos"] - target_coverage_df["startpos"] + 1
                )

                target_coverage_df["coverage"] = (
                    target_coverage_df["bases"] / target_coverage_df["length"]
                )

                # Reorder columns to match reference format: chrom,startpos,endpos,name,length,coverage,bases
                target_coverage_df = target_coverage_df[
                    [
                        "chrom",
                        "startpos",
                        "endpos",
                        "name",
                        "length",
                        "coverage",
                        "bases",
                    ]
                ]

                target_coverage_df.to_csv(
                    os.path.join(sample_output_dir, "target_coverage.csv"),
                    index=False,
                )

                target_result.processing_steps.append("target_coverage_calculated")

                # Step 13b: Create target coverage with timestamp and read counts (optimized)
                logger.info("Calculating read counts per target...")
                new_read_counts_df = get_read_counts_per_target(file_path, self.bedfile)
                
                if not new_read_counts_df.empty:
                    # Use optimized approach: store latest cumulative reads in separate Parquet file for fast access
                    time_coverage_file = os.path.join(sample_output_dir, "target_coverage_time.csv")
                    latest_reads_cache = os.path.join(sample_output_dir, "_target_coverage_latest_reads.parquet")
                    
                    # Load previous cumulative reads from cache (much faster than reading entire CSV)
                    previous_cumulative_reads = None
                    if os.path.exists(latest_reads_cache):
                        try:
                            previous_cumulative_reads = pd.read_parquet(
                                latest_reads_cache, **_PARQUET_READ_KWARGS
                            )
                            previous_cumulative_reads.rename(columns={'reads': 'previous_reads'}, inplace=True)
                        except Exception as e:
                            logger.debug(f"Could not load cached latest reads, will try CSV: {e}")
                            # Fallback to CSV if cache doesn't exist
                            if os.path.exists(time_coverage_file):
                                try:
                                    # Only read last chunk for efficiency
                                    existing_time_df = pd.read_csv(time_coverage_file)
                                    if not existing_time_df.empty:
                                        latest_timestamp = existing_time_df['timestamp'].max()
                                        previous_cumulative_reads = existing_time_df[
                                            existing_time_df['timestamp'] == latest_timestamp
                                        ][['chrom', 'startpos', 'endpos', 'name', 'reads']].copy()
                                        previous_cumulative_reads.rename(columns={'reads': 'previous_reads'}, inplace=True)
                                except Exception as e2:
                                    logger.warning(f"Error loading existing target_coverage_time.csv: {e2}")
                    
                    # Merge new read counts with target coverage data
                    target_coverage_with_reads = target_coverage_df.merge(
                        new_read_counts_df[['chrom', 'startpos', 'endpos', 'name', 'reads']],
                        on=['chrom', 'startpos', 'endpos', 'name'],
                        how='left'
                    )
                    # Fill missing reads with 0
                    target_coverage_with_reads['reads'] = target_coverage_with_reads['reads'].fillna(0).astype(int)
                    
                    # Accumulate with previous cumulative reads if available
                    if previous_cumulative_reads is not None:
                        target_coverage_with_reads = target_coverage_with_reads.merge(
                            previous_cumulative_reads,
                            on=['chrom', 'startpos', 'endpos', 'name'],
                            how='left'
                        )
                        target_coverage_with_reads['previous_reads'] = target_coverage_with_reads['previous_reads'].fillna(0).astype(int)
                        # Add new reads to previous cumulative reads
                        target_coverage_with_reads['reads'] = (
                            target_coverage_with_reads['reads'] + target_coverage_with_reads['previous_reads']
                        )
                        target_coverage_with_reads.drop(columns=['previous_reads'], inplace=True)
                    
                    # Calculate normalized reads (reads per length)
                    target_coverage_with_reads['reads_per_length'] = (
                        target_coverage_with_reads['reads'] / target_coverage_with_reads['length']
                    )
                    
                    # Add timestamp
                    if self.simtime and timestamp:
                        current_timestamp = timestamp * 1000
                    else:
                        current_timestamp = time.time() * 1000
                    target_coverage_with_reads['timestamp'] = current_timestamp
                    
                    # Reorder columns: chrom, startpos, endpos, name, length, coverage, bases, timestamp, reads, reads_per_length
                    target_coverage_with_reads = target_coverage_with_reads[
                        ['chrom', 'startpos', 'endpos', 'name', 'length', 'coverage', 'bases',
                         'timestamp', 'reads', 'reads_per_length']
                    ]
                    
                    # Save latest cumulative reads to cache for next time (fast access)
                    try:
                        target_coverage_with_reads[
                            ['chrom', 'startpos', 'endpos', 'name', 'reads']
                        ].to_parquet(latest_reads_cache, **_PARQUET_WRITE_KWARGS)
                    except Exception as e:
                        logger.debug(f"Could not save latest reads cache: {e}")
                    
                    # Append to CSV using append mode (much faster than reading entire file)
                    try:
                        # Check if file exists to determine if we need header
                        file_exists = os.path.exists(time_coverage_file)
                        target_coverage_with_reads.to_csv(
                            time_coverage_file, 
                            mode='a', 
                            header=not file_exists,
                            index=False
                        )
                    except Exception as e:
                        logger.warning(f"Error appending to target_coverage_time.csv: {e}")
                        # Fallback: read and concat (slower but works)
                        if os.path.exists(time_coverage_file):
                            try:
                                existing_time_df = pd.read_csv(time_coverage_file)
                                target_coverage_with_reads = pd.concat(
                                    [existing_time_df, target_coverage_with_reads],
                                    ignore_index=True
                                )
                                target_coverage_with_reads.to_csv(time_coverage_file, index=False)
                            except Exception as e2:
                                logger.error(f"Error in fallback CSV write: {e2}")
                    
                    logger.info(f"Saved target coverage with timestamp and cumulative read counts: {time_coverage_file}")
                    target_result.processing_steps.append("target_coverage_time_saved")
                else:
                    logger.warning("No read counts extracted, skipping target_coverage_time.csv")

                # Step 13: Identify targets exceeding threshold
                run_list = target_coverage_df[
                    target_coverage_df["coverage"].ge(self.callthreshold)
                ]

                # Load existing targets exceeding threshold count
                targets_exceeding_file = os.path.join(
                    sample_output_dir, "targets_exceeding_threshold_count.txt"
                )
                existing_targets_count = 0
                if os.path.exists(targets_exceeding_file):
                    try:
                        with open(targets_exceeding_file, "r") as f:
                            existing_targets_count = int(f.read().strip())
                    except (ValueError, IOError):
                        existing_targets_count = 0

                # Always generate targets_exceeding_threshold.bed for potential SNP analysis
                if len(run_list) > 0:
                    logger.info(f"Found {len(run_list)} regions exceeding threshold")

                    # Save updated count
                    with open(targets_exceeding_file, "w") as f:
                        f.write(str(len(run_list)))

                    # Generate BED file for regions exceeding threshold
                    run_list[["chrom", "startpos", "endpos"]].to_csv(
                        os.path.join(
                            sample_output_dir, "targets_exceeding_threshold.bed"
                        ),
                        sep="\t",
                        header=None,
                        index=None,
                    )

                    logger.info(
                        f"Targets exceeding threshold saved: {len(run_list)} regions"
                    )

                    # Log the results for debugging (no storage to avoid memory leaks)
                    logger.info(
                        f"Generated BED file with {len(run_list)} regions exceeding threshold"
                    )
                else:
                    logger.info("No regions exceed coverage threshold")

                    # Create empty BED file to avoid errors
                    empty_bed_file = os.path.join(
                        sample_output_dir, "targets_exceeding_threshold.bed"
                    )
                    with open(empty_bed_file, "w") as f:
                        # Write minimal BED header
                        pass

                    logger.info("Created empty targets_exceeding_threshold.bed file")

                    # Log the empty file creation (no storage to avoid memory leaks)
                    logger.info("Created empty BED file - no regions exceed threshold")

                # Log reference genome status for SNP calling
                if self.reference:
                    logger.info(
                        f"Reference genome found: {self.reference} - SNP calling will be available"
                    )
                else:
                    logger.info(
                        "No reference genome provided - SNP calling will not be available"
                    )

                target_result.processing_steps.append("threshold_analysis_complete")

                # Step 13: Perform target analysis
                logger.info("Performing target analysis...")
                analysis_results = self._perform_target_analysis(
                    file_path, temp_dir, metadata, newcovdf, bedcovdf
                )
                target_result.analysis_results = analysis_results
                target_result.processing_steps.append("analysis_completed")

                # Step 14: Generate output files
                logger.info("Generating output files...")
                output_files = self._generate_output_files(
                    sample_output_dir,
                    analysis_counter,
                    analysis_results,
                    newcovdf,
                    bedcovdf,
                    logger,
                    metadata,
                )
                # No JSON files created - data is accumulated in CSV files
                target_result.target_data_path = None
                target_result.target_plot_path = None
                target_result.processing_steps.append("output_generated")

                # Step 15: Update analysis counter
                analysis_counter += 1
                self._save_analysis_counter(
                    sample_id, analysis_counter, self.work_dir, logger
                )
                target_result.processing_steps.append("counter_updated")

                target_result.processing_steps.append("target_analysis_complete")
                logger.info(f"Target analysis complete for {sample_id}")
                logger.info(f"   Results: {len(analysis_results)} items")
                logger.info(f"   Output files: {len(output_files)} files")
                logger.info(
                    f"   Targets exceeding threshold: {len(run_list) if 'run_list' in locals() else 0}"
                )

                return target_result

        except Exception as e:
            error_details = f"Error in target analysis for {sample_id}: {str(e)}"
            logger.error(error_details)
            target_result.error_message = error_details
            target_result.processing_steps.append("analysis_failed")
            return target_result

    def _perform_target_analysis(
        self,
        file_path: str,
        temp_dir: str,
        metadata: Dict[str, Any],
        newcovdf: pd.DataFrame,
        bedcovdf: pd.DataFrame,
    ) -> Dict[str, Any]:
        """
        Perform the actual target analysis.

        Args:
            file_path: Path to input file
            temp_dir: Temporary directory for processing
            metadata: File metadata
            newcovdf: Genome coverage DataFrame
            bedcovdf: Target coverage DataFrame

        Returns:
            Dictionary containing analysis results
        """
        logger = logging.getLogger("robin.target")

        # This is where the actual target analysis logic would go
        # For now, we'll create a placeholder implementation that uses the coverage data

        results = {
            "analysis_type": self.config.get("analysis_type", "standard"),
            "input_file": os.path.basename(file_path),
            "sample_id": metadata.get("sample_id", "unknown"),
            "timestamp": time.time(),
            "targets_found": 0,
            "targets_analyzed": [],
            "coverage_summary": {
                "genome_regions": len(newcovdf) if newcovdf is not None else 0,
                "target_regions": len(bedcovdf) if bedcovdf is not None else 0,
                "total_bases_covered": 0,
                "average_coverage": 0.0,
            },
            "summary_stats": {
                "total_targets": 0,
                "significant_targets": 0,
                "average_score": 0.0,
            },
        }

        # Calculate coverage statistics
        if newcovdf is not None and not newcovdf.empty:
            results["coverage_summary"]["total_bases_covered"] = newcovdf.get(
                "numreads", pd.Series([0])
            ).sum()
            results["coverage_summary"]["average_coverage"] = newcovdf.get(
                "numreads", pd.Series([0])
            ).mean()

        # Analyze target regions
        if bedcovdf is not None and not bedcovdf.empty:
            results["targets_found"] = len(bedcovdf)
            results["summary_stats"]["total_targets"] = len(bedcovdf)

            # Calculate significance based on coverage
            significant_threshold = self.config.get("threshold", 0.5)
            significant_targets = bedcovdf[
                bedcovdf["bases"] > significant_threshold * bedcovdf["bases"].max()
            ]
            results["summary_stats"]["significant_targets"] = len(significant_targets)

            # Add target analysis results (list comp inlined in 3.12)
            max_bases = bedcovdf["bases"].max()
            results["targets_analyzed"] = [
                {
                    "target_id": row.get("name", f"target_{idx+1:03d}"),
                    "position": f"{row['chrom']}:{row['startpos']}-{row['endpos']}",
                    "coverage": row["bases"],
                    "significance": (
                        "high" if row["bases"] > significant_threshold * max_bases else "medium"
                    ),
                }
                for idx, row in bedcovdf.iterrows()
            ]

            # Calculate average score
            if results["targets_analyzed"]:
                scores = [t["coverage"] for t in results["targets_analyzed"]]
                results["summary_stats"]["average_score"] = sum(scores) / len(scores)

        logger.info(f"Analysis completed: {results['targets_found']} targets found")
        logger.info(f"Coverage summary: {results['coverage_summary']}")
        return results

    def _generate_output_files(
        self,
        sample_dir: str,
        analysis_counter: int,
        analysis_results: Dict[str, Any],
        newcovdf: pd.DataFrame,
        bedcovdf: pd.DataFrame,
        logger,
        metadata: Dict[str, Any],
    ) -> Dict[str, str]:
        """
        Generate minimal output files for the analysis.

        Args:
            sample_dir: Sample output directory
            analysis_counter: Current analysis counter
            analysis_results: Analysis results
            newcovdf: Genome coverage DataFrame
            bedcovdf: Target coverage DataFrame
            logger: Logger instance
            metadata: File metadata

        Returns:
            Dictionary with paths to generated files
        """
        output_files = {}

        try:
            # No JSON files needed - accumulated data is already saved in CSV files
            # The coverage data is persisted in:
            # - coverage_main.csv (genome coverage)
            # - bed_coverage_main.csv (target coverage)
            # - coverage_time_chart.npy (coverage over time)
            # - target_coverage.csv (calculated target coverage)

            logger.info(
                "No additional output files needed - data accumulated in CSV files"
            )
            return output_files

        except Exception as e:
            logger.error(f"Error in output file generation: {e}")
            return {}

    def _load_analysis_counter(self, sample_id: str, work_dir: str, logger) -> int:
        """Load the analysis counter for a sample from disk"""
        counter_file = os.path.join(work_dir, sample_id, "target_analysis_counter.txt")
        if os.path.exists(counter_file):
            try:
                with open(counter_file, "r") as f:
                    return int(f.read().strip())
            except (ValueError, IOError) as e:
                logger.warning(f"Error loading counter for {sample_id}: {e}")
        return 0

    def _save_analysis_counter(
        self, sample_id: str, counter: int, work_dir: str, logger
    ) -> None:
        """Save the analysis counter for a sample to disk"""
        counter_file = os.path.join(work_dir, sample_id, "target_analysis_counter.txt")
        try:
            # Create sample directory if it doesn't exist
            os.makedirs(os.path.dirname(counter_file), exist_ok=True)
            with open(counter_file, "w") as f:
                f.write(str(counter))
        except IOError as e:
            logger.error(f"Error saving counter for {sample_id}: {e}")

    def _load_existing_coverage_data(
        self, sample_dir: str, filename: str, logger
    ) -> Optional[pd.DataFrame]:
        """Load existing coverage data from CSV file if it exists"""
        file_path = os.path.join(sample_dir, filename)
        if os.path.exists(file_path):
            try:
                df = pd.read_csv(file_path)
                logger.debug(
                    f"Loaded existing coverage data from {filename}: {df.shape}"
                )

                # For bed_coverage_main.csv, extract only the raw data columns needed for merging
                if (
                    filename == "bed_coverage_main.csv"
                    and "length" in df.columns
                    and "coverage" in df.columns
                ):
                    # Extract only the raw data columns: chrom,startpos,endpos,name,bases
                    df = df[["chrom", "startpos", "endpos", "name", "bases"]]
                    logger.debug(
                        f"Extracted raw data columns from {filename}: {df.shape}"
                    )

                return df
            except Exception as e:
                logger.warning(
                    f"Error loading existing coverage data from {filename}: {e}"
                )
                return None
        return None

    def _load_existing_coverage_over_time(
        self, sample_dir: str, logger
    ) -> Optional[np.ndarray]:
        """Load existing coverage over time data from numpy file if it exists"""
        file_path = os.path.join(sample_dir, "coverage_time_chart.npy")
        if os.path.exists(file_path):
            try:
                data = np.load(file_path)
                logger.debug(f"Loaded existing coverage over time data: {data.shape}")
                return data
            except Exception as e:
                logger.warning(f"Error loading existing coverage over time data: {e}")
                return None
        return None
    
    def accumulate_staged_files(self, sample_id: str, force: bool = False) -> Dict[str, Any]:
        """
        Batch accumulation of staged files with minimal lock hold time.
        Lock is held only to: (1) claim staging files, (2) load existing + merge, (3) write outputs + cleanup.
        Heavy work (load parquets, BAM filtering, read counts) runs outside the lock.
        """
        logger = logging.getLogger("robin.target")
        lock_file = self._get_lock_file(sample_id, "accumulation")
        staging_dir = self._get_staging_dir(sample_id)
        sample_output_dir = os.path.join(self.work_dir, sample_id)
        start_time = time.time()

        try:
            # --- Critical section 1: claim staging files (short) ---
            batch_dir = None
            num_claimed = 0
            with FileLock(lock_file, timeout=60.0):
                complete_sets, staging_counts, incomplete_ids = self._collect_complete_staging_sets(
                    staging_dir
                )
                if not complete_sets:
                    if any(staging_counts.values()):
                        logger.warning(
                            "Skipping accumulation for %s - no complete staging sets "
                            "(counts=%s, incomplete_ids=%d)",
                            sample_id,
                            staging_counts,
                            len(incomplete_ids),
                        )
                    else:
                        logger.info(f"No staged files to accumulate for {sample_id}")
                    return {
                        "status": "no_files",
                        "files_processed": 0,
                        "files_pending": 0,
                        "incomplete_staging_sets": len(incomplete_ids),
                    }
                if incomplete_ids:
                    logger.warning(
                        "Found %d incomplete staging sets for %s; accumulating only complete sets",
                        len(incomplete_ids),
                        sample_id,
                    )
                complete_count = len(complete_sets)
                if not force and complete_count < self.batch_size:
                    logger.info(
                        f"Skipping accumulation - only {complete_count} complete files staged "
                        f"(threshold: {self.batch_size}, force={force})"
                    )
                    return {
                        "status": "below_threshold",
                        "files_pending": complete_count,
                        "incomplete_staging_sets": len(incomplete_ids),
                    }
                n_claim = complete_count if force else min(self.batch_size, complete_count)
                batch_id = str(uuid.uuid4())
                batch_dir = os.path.join(staging_dir, f"_batch_{batch_id}")
                os.makedirs(batch_dir, exist_ok=True)
                for coverage_src, bedcov_src, timestamp_src, source_bam_src in complete_sets[:n_claim]:
                    for src in (coverage_src, bedcov_src, timestamp_src, source_bam_src):
                        if os.path.exists(src):
                            shutil.move(src, os.path.join(batch_dir, os.path.basename(src)))
                num_claimed = n_claim
            # Lock released

            # --- No lock: load from batch dir ---
            batch_sets, batch_counts, batch_incomplete = self._collect_complete_staging_sets(
                batch_dir
            )
            if batch_incomplete:
                logger.warning(
                    "Batch %s has %d incomplete staging sets after claim (counts=%s)",
                    os.path.basename(batch_dir),
                    len(batch_incomplete),
                    batch_counts,
                )
            cov_frames: List[pd.DataFrame] = []
            bed_frames: List[pd.DataFrame] = []
            timestamps = []
            source_bam_paths = []
            loaded_files = 0
            for cov_file, bed_file, ts_file, source_bam_file in batch_sets:
                try:
                    cov_frames.append(
                        pd.read_parquet(cov_file, **_PARQUET_READ_KWARGS)
                    )
                    bed_frames.append(
                        pd.read_parquet(bed_file, **_PARQUET_READ_KWARGS)
                    )
                    with open(ts_file, "r") as f:
                        timestamps.append(float(f.read().strip()))
                    with open(source_bam_file, "r") as f:
                        source_bam_paths.append(f.read().strip())
                    loaded_files += 1
                except Exception as e:
                    logger.warning(f"Error loading staging file {cov_file}: {e}")
            if not cov_frames or not bed_frames:
                batch_covdf = None
                batch_bedcovdf = None
            else:
                batch_covdf = (
                    pd.concat(cov_frames, ignore_index=True)
                    .groupby(
                        ["#rname", "startpos", "endpos"],
                        as_index=False,
                        sort=False,
                    )
                    .agg(
                        {
                            "numreads": "sum",
                            "covbases": "sum",
                            "meandepth": "sum",
                        }
                    )
                )
                batch_bedcovdf = (
                    pd.concat(bed_frames, ignore_index=True)
                    .groupby(
                        ["chrom", "startpos", "endpos", "name"],
                        as_index=False,
                        sort=False,
                    )
                    .agg({"bases": "sum"})
                )
            if batch_covdf is None or batch_bedcovdf is None:
                logger.error("No staging files could be loaded successfully")
                if batch_dir and os.path.exists(batch_dir):
                    for f in glob.glob(os.path.join(batch_dir, "*")):
                        try:
                            os.remove(f)
                        except OSError:
                            pass
                    try:
                        os.rmdir(batch_dir)
                    except OSError:
                        pass
                return {"status": "load_failed", "files_attempted": num_claimed}
            logger.info(f"Loaded {loaded_files} staging files; merging with accumulated data...")

            # --- Critical section 2: load existing + merge (short) ---
            existing_covdf = None
            existing_bedcovdf = None
            existing_coverage_over_time = None
            updated_covdf = None
            updated_bedcovdf = None
            with FileLock(lock_file, timeout=60.0):
                existing_covdf = self._load_existing_coverage_data(
                    sample_output_dir, "coverage_main.csv", logger
                )
                existing_bedcovdf = self._load_existing_coverage_data(
                    sample_output_dir, "bed_coverage_main.csv", logger
                )
                existing_coverage_over_time = self._load_existing_coverage_over_time(
                    sample_output_dir, logger
                )
                if existing_covdf is not None:
                    updated_covdf, updated_bedcovdf = run_bedmerge(
                        batch_covdf, existing_covdf, batch_bedcovdf, existing_bedcovdf
                    )
                else:
                    updated_covdf, updated_bedcovdf = batch_covdf, batch_bedcovdf
            # Lock released

            logger.info(f"Final accumulated: genome={updated_covdf.shape}, targets={updated_bedcovdf.shape}")
            bases = updated_covdf["covbases"].sum()
            genome = updated_covdf["endpos"].sum()
            coverage = bases / genome if genome > 0 else 0.0
            logger.info(f"Coverage: {coverage:.4f} ({bases} bases / {genome} genome)")
            current_timestamp = max(timestamps) if timestamps else time.time() * 1000
            if existing_coverage_over_time is not None:
                updated_coverage_over_time = np.vstack(
                    [existing_coverage_over_time, [[current_timestamp, coverage]]]
                )
            else:
                updated_coverage_over_time = np.array([[current_timestamp, coverage]])

            target_coverage_df = updated_bedcovdf.copy()
            target_coverage_df["length"] = (
                target_coverage_df["endpos"] - target_coverage_df["startpos"] + 1
            )
            target_coverage_df["coverage"] = (
                target_coverage_df["bases"] / target_coverage_df["length"]
            )
            target_coverage_df = target_coverage_df[
                ["chrom", "startpos", "endpos", "name", "length", "coverage", "bases"]
            ]
            run_list = target_coverage_df[target_coverage_df["coverage"].ge(self.callthreshold)]

            # --- No lock: BAM read counts and filtered BAMs ---
            batch_read_counts = None
            if source_bam_paths:
                batch_read_counts_list = []
                for bam_path in source_bam_paths:
                    read_counts_df = get_read_counts_per_target(bam_path, self.bedfile)
                    if not read_counts_df.empty:
                        batch_read_counts_list.append(read_counts_df)
                if batch_read_counts_list:
                    batch_read_counts = pd.concat(batch_read_counts_list, ignore_index=True)
                    batch_read_counts = batch_read_counts.groupby(
                        ['chrom', 'startpos', 'endpos', 'name'], as_index=False
                    ).agg({'reads': 'sum'})

            if len(run_list) > 0:
                logger.info("Processing filtered BAM files...")
                try:
                    original_bam_files = list(set(source_bam_paths))
                    if original_bam_files:
                        targets_bed = self._get_master_bed_path(sample_id) or self.bedfile
                        if targets_bed != self.bedfile:
                            logger.info(f"Using master BED file for target.bam: {targets_bed}")
                        else:
                            logger.info(f"Using original target panel BED file for target.bam: {targets_bed}")
                        filtered_bams_dir = os.path.join(sample_output_dir, "_filtered_bams")
                        os.makedirs(filtered_bams_dir, exist_ok=True)
                        new_filtered_bams = []
                        bed_regions_cache = None
                        try:
                            bed_regions_cache = _load_bed_regions(targets_bed)
                        except Exception as e:
                            logger.warning(f"Could not preload BED regions from {targets_bed}: {e}")
                        for i, source_bam in enumerate(original_bam_files):
                            filtered_bam_name = f"filtered_{i:06d}_{os.path.basename(source_bam)}"
                            filtered_bam_path = os.path.join(filtered_bams_dir, filtered_bam_name)
                            logger.info(f"Filtering BAM {i+1}/{len(original_bam_files)}: {os.path.basename(source_bam)}")
                            run_bedtools(source_bam, targets_bed, filtered_bam_path, regions=bed_regions_cache)
                            if os.path.exists(filtered_bam_path):
                                new_filtered_bams.append(filtered_bam_path)
                        if new_filtered_bams:
                            batch_timestamp = int(time.time() * 1000)
                            batch_merged_bam = os.path.join(sample_output_dir, f"batch_{batch_timestamp}.bam")
                            if len(new_filtered_bams) > 1:
                                pysam.merge("-o", batch_merged_bam, *new_filtered_bams)
                            else:
                                shutil.copy2(new_filtered_bams[0], batch_merged_bam)
                            pysam.index(batch_merged_bam)
                            for filtered_bam in new_filtered_bams:
                                try:
                                    if os.path.exists(filtered_bam):
                                        os.remove(filtered_bam)
                                    if os.path.exists(f"{filtered_bam}.bai"):
                                        os.remove(f"{filtered_bam}.bai")
                                except OSError as e:
                                    logger.warning(f"Could not remove filtered BAM {filtered_bam}: {e}")
                            try:
                                if not os.listdir(filtered_bams_dir):
                                    os.rmdir(filtered_bams_dir)
                            except OSError:
                                pass
                except Exception as e:
                    logger.error(f"Error creating target.bam: {e}")
                    import traceback
                    logger.error(traceback.format_exc())

            # --- Critical section 3: write all outputs + cleanup (short) ---
            time_coverage_file = os.path.join(sample_output_dir, "target_coverage_time.csv")
            latest_reads_cache = os.path.join(sample_output_dir, "_target_coverage_latest_reads.parquet")
            with FileLock(lock_file, timeout=60.0):
                np.save(
                    os.path.join(sample_output_dir, "coverage_time_chart.npy"),
                    updated_coverage_over_time,
                )
                updated_covdf.to_csv(
                    os.path.join(sample_output_dir, "coverage_main.csv"), index=False
                )
                bed_coverage_main_df = updated_bedcovdf.copy()
                bed_coverage_main_df["length"] = (
                    bed_coverage_main_df["endpos"] - bed_coverage_main_df["startpos"] + 1
                )
                bed_coverage_main_df["coverage"] = (
                    bed_coverage_main_df["bases"] / bed_coverage_main_df["length"]
                )
                bed_coverage_main_df = bed_coverage_main_df[
                    ["chrom", "startpos", "endpos", "name", "length", "coverage", "bases"]
                ]
                bed_coverage_main_df.to_csv(
                    os.path.join(sample_output_dir, "bed_coverage_main.csv"), index=False
                )
                target_coverage_df.to_csv(
                    os.path.join(sample_output_dir, "target_coverage.csv"), index=False
                )
                if batch_read_counts is not None:
                    previous_cumulative_reads = None
                    if os.path.exists(latest_reads_cache):
                        try:
                            previous_cumulative_reads = pd.read_parquet(
                                latest_reads_cache, **_PARQUET_READ_KWARGS
                            )
                            previous_cumulative_reads.rename(columns={'reads': 'previous_reads'}, inplace=True)
                        except Exception:
                            pass
                    if previous_cumulative_reads is None and os.path.exists(time_coverage_file):
                        try:
                            existing_time_df = pd.read_csv(time_coverage_file)
                            if not existing_time_df.empty:
                                latest_ts = existing_time_df['timestamp'].max()
                                previous_cumulative_reads = existing_time_df[
                                    existing_time_df['timestamp'] == latest_ts
                                ][['chrom', 'startpos', 'endpos', 'name', 'reads']].copy()
                                previous_cumulative_reads.rename(columns={'reads': 'previous_reads'}, inplace=True)
                        except Exception:
                            pass
                    target_coverage_with_reads = target_coverage_df.merge(
                        batch_read_counts[['chrom', 'startpos', 'endpos', 'name', 'reads']],
                        on=['chrom', 'startpos', 'endpos', 'name'], how='left'
                    )
                    target_coverage_with_reads['reads'] = target_coverage_with_reads['reads'].fillna(0).astype(int)
                    if previous_cumulative_reads is not None:
                        target_coverage_with_reads = target_coverage_with_reads.merge(
                            previous_cumulative_reads,
                            on=['chrom', 'startpos', 'endpos', 'name'], how='left'
                        )
                        target_coverage_with_reads['previous_reads'] = target_coverage_with_reads['previous_reads'].fillna(0).astype(int)
                        target_coverage_with_reads['reads'] = (
                            target_coverage_with_reads['reads'] + target_coverage_with_reads['previous_reads']
                        )
                        target_coverage_with_reads.drop(columns=['previous_reads'], inplace=True)
                    target_coverage_with_reads['reads_per_length'] = (
                        target_coverage_with_reads['reads'] / target_coverage_with_reads['length']
                    )
                    target_coverage_with_reads['timestamp'] = current_timestamp
                    target_coverage_with_reads = target_coverage_with_reads[
                        ['chrom', 'startpos', 'endpos', 'name', 'length', 'coverage', 'bases',
                         'timestamp', 'reads', 'reads_per_length']
                    ]
                    try:
                        target_coverage_with_reads[
                            ['chrom', 'startpos', 'endpos', 'name', 'reads']
                        ].to_parquet(latest_reads_cache, **_PARQUET_WRITE_KWARGS)
                    except Exception:
                        pass
                    file_exists = os.path.exists(time_coverage_file)
                    target_coverage_with_reads.to_csv(
                        time_coverage_file, mode='a', header=not file_exists, index=False
                    )
                targets_exceeding_file = os.path.join(
                    sample_output_dir, "targets_exceeding_threshold_count.txt"
                )
                if len(run_list) > 0:
                    with open(targets_exceeding_file, "w") as f:
                        f.write(str(len(run_list)))
                    run_list[["chrom", "startpos", "endpos"]].to_csv(
                        os.path.join(sample_output_dir, "targets_exceeding_threshold.bed"),
                        sep="\t", header=None, index=None,
                    )
                else:
                    with open(
                        os.path.join(sample_output_dir, "targets_exceeding_threshold.bed"), "w"
                    ) as f:
                        pass
                if batch_dir and os.path.exists(batch_dir):
                    for f in glob.glob(os.path.join(batch_dir, "*")):
                        try:
                            os.remove(f)
                        except OSError:
                            pass
                    try:
                        os.rmdir(batch_dir)
                    except OSError:
                        pass

            elapsed = time.time() - start_time
            logger.info(
                f"Batch accumulation complete for {sample_id}: "
                f"{num_claimed} files in {elapsed:.2f}s"
            )
            return {
                "status": "success",
                "files_processed": num_claimed,
                "coverage": coverage,
                "targets_exceeding_threshold": len(run_list) if len(run_list) > 0 else 0,
                "elapsed_time": elapsed,
            }

        except TimeoutError as e:
            logger.error(f"Could not acquire accumulation lock for {sample_id}: {e}")
            return {"status": "lock_timeout", "error": str(e)}
        except Exception as e:
            logger.error(f"Error during batch accumulation for {sample_id}: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return {"status": "error", "error": str(e)}


def process_single_file(
    file_path: str,
    metadata: Dict[str, Any],
    work_dir: str,
    logger,
    reference: Optional[str] = None,
    target_panel: str = "rCNS2",
) -> Dict[str, Any]:
    """
    Process a single file for target analysis using the complete pipeline.

    Args:
        file_path: Path to input file
        metadata: File metadata
        work_dir: Working directory
        logger: Logger instance
        reference: Optional reference genome path
        target_panel: Target panel type (rCNS2, AML, or custom panel)

    Returns:
        Dictionary with target analysis results
    """
    sample_id = metadata.get("sample_id", "unknown")
    logger.info(f"Starting target analysis for sample: {sample_id}")

    # Create sample-specific output directory
    sample_output_dir = os.path.join(work_dir, sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)
    logger.info(f"Created sample output directory: {sample_output_dir}")

    analysis_result = {
        "sample_id": sample_id,
        "file_path": file_path,
        "analysis_timestamp": time.time(),
        "target_data_path": None,
        "target_plot_path": None,
        "analysis_results": {},
        "processing_steps": [],
        "error_message": None,
        "coverage_data": {},
        "target_bam_path": None,
        "coverage_over_time": None,
    }

    try:
        # Initialize target analysis
        target_analysis = TargetAnalysis(work_dir=work_dir, target_panel=target_panel)

        # Set reference genome if provided
        if reference:
            target_analysis.reference = reference
            logger.info(f"Using reference genome: {reference}")

        # Process the file
        target_metadata = target_analysis.process_file(file_path, metadata)

        # Update analysis result with target metadata
        analysis_result.update(
            {
                "target_data_path": target_metadata.target_data_path,
                "target_plot_path": target_metadata.target_plot_path,
                "analysis_results": target_metadata.analysis_results,
                "processing_steps": target_metadata.processing_steps,
                "error_message": target_metadata.error_message,
                "coverage_data": target_metadata.coverage_data,
                "target_bam_path": target_metadata.target_bam_path,
                "coverage_over_time": target_metadata.coverage_over_time,
            }
        )

        if target_metadata.error_message:
            logger.error(f"Target analysis failed: {target_metadata.error_message}")
        else:
            logger.info(f"Target analysis complete for {sample_id}")
            logger.info(f"   Results: {len(target_metadata.analysis_results)} items")
            logger.info(
                f"   Processing steps: {', '.join(target_metadata.processing_steps)}"
            )
            logger.info(
                f"   Coverage data: genome={target_metadata.coverage_data.get('genome_coverage_shape', (0,0))}, targets={target_metadata.coverage_data.get('target_coverage_shape', (0,0))}"
            )
            if target_metadata.target_bam_path:
                logger.info(f"   Target BAM: {target_metadata.target_bam_path}")

        return analysis_result

    except Exception as e:
        error_details = f"Error in target analysis for {sample_id}: {str(e)}"
        logger.error(error_details)
        analysis_result["error_message"] = error_details
        analysis_result["processing_steps"].append("analysis_failed")
        return analysis_result


def process_multiple_files(bam_paths, metadata_list, work_dir, logger, reference=None, target_panel=None):
    """
    Process multiple BAM files for target analysis using staged processing.
    
    This function processes multiple BAM files for the same sample using the
    existing staging infrastructure. Each file is processed individually and
    staged, then all staged files are accumulated in a single batch operation.

    Args:
        bam_paths: List of paths to BAM files
        metadata_list: List of metadata dictionaries (one per BAM file)
        work_dir: Working directory
        logger: Logger instance
        reference: Optional path to reference genome for SNP calling
        target_panel: Target panel type (rCNS2, AML, or custom panel)

    Returns:
        Dictionary with aggregated target analysis results
    """
    if not bam_paths or not metadata_list:
        raise ValueError("bam_paths and metadata_list must not be empty")
    
    if len(bam_paths) != len(metadata_list):
        raise ValueError("bam_paths and metadata_list must have the same length")
    
    # Get sample ID from first metadata (assuming all BAMs are from same sample)
    sample_id = metadata_list[0].get("sample_id", "unknown")
    
    logger.info(f"🎯 Starting multi-file target analysis for sample: {sample_id}")
    logger.info(f"Processing {len(bam_paths)} BAM files for sample {sample_id}")
    
    # Log essential metadata only
    for i, (bam_path, metadata) in enumerate(zip(bam_paths, metadata_list)):
        logger.debug(f"BAM file {i+1}: {os.path.basename(bam_path)}")

    analysis_result = {
        "sample_id": sample_id,
        "bam_paths": bam_paths,
        "analysis_timestamp": time.time(),
        "processing_steps": [],
        "error_message": None,
        "files_processed": 0,
        "total_files": len(bam_paths),
        "coverage_data": {},
        "target_bam_path": None,
        "coverage_over_time": None,
    }

    try:
        # Initialize target analysis with staging enabled
        target_analysis = TargetAnalysis(
            work_dir=work_dir, 
            target_panel=target_panel,
            batch_size=1,  # Force accumulation after each batch
            use_staging=True
        )
        
        # Set reference genome if provided
        if reference:
            target_analysis.reference = reference
            logger.info(f"Using reference genome: {reference}")

        # Process each BAM file individually using staging
        logger.info("Processing files with staging (fast path)")
        processed_files = 0
        
        for i, (bam_path, metadata) in enumerate(zip(bam_paths, metadata_list)):
            logger.info(f"Processing BAM file {i+1}/{len(bam_paths)}: {os.path.basename(bam_path)}")
            
            try:
                # Process file with staging
                target_metadata, should_accumulate = target_analysis.process_file_with_staging(
                    bam_path, metadata
                )
                
                if target_metadata.error_message:
                    logger.warning(f"Error processing {os.path.basename(bam_path)}: {target_metadata.error_message}")
                    continue
                
                processed_files += 1
                logger.debug(f"Successfully staged file {i+1}: {os.path.basename(bam_path)}")
                
            except Exception as e:
                logger.warning(f"Error processing {os.path.basename(bam_path)}: {e}")
                continue

        if processed_files == 0:
            analysis_result["error_message"] = "No files could be processed successfully"
            analysis_result["processing_steps"].append("no_files_processed")
            return analysis_result

        analysis_result["files_processed"] = processed_files
        analysis_result["processing_steps"].append("files_staged")

        # Force accumulation of all staged files
        logger.info(f"Accumulating {processed_files} staged files for sample {sample_id}")
        accumulation_result = target_analysis.accumulate_staged_files(
            sample_id, force=True
        )
        
        if accumulation_result.get("status") != "success":
            analysis_result["error_message"] = f"Accumulation failed: {accumulation_result.get('error', 'Unknown error')}"
            analysis_result["processing_steps"].append("accumulation_failed")
            return analysis_result

        analysis_result["processing_steps"].append("accumulation_complete")
        logger.info(f"Accumulation completed: {accumulation_result}")

        # Load final accumulated data for result metadata
        sample_output_dir = os.path.join(work_dir, sample_id)
        
        # Load final coverage data
        final_covdf = target_analysis._load_existing_coverage_data(
            sample_output_dir, "coverage_main.csv", logger
        )
        final_bedcovdf = target_analysis._load_existing_coverage_data(
            sample_output_dir, "bed_coverage_main.csv", logger
        )
        final_coverage_over_time = target_analysis._load_existing_coverage_over_time(
            sample_output_dir, logger
        )
        
        # Store final results
        analysis_result["coverage_data"] = {
            "genome_coverage_shape": final_covdf.shape if final_covdf is not None else (0, 0),
            "target_coverage_shape": final_bedcovdf.shape if final_bedcovdf is not None else (0, 0),
            "coverage": accumulation_result.get("coverage", 0.0),
            "targets_exceeding_threshold": accumulation_result.get("targets_exceeding_threshold", 0),
        }
        
        analysis_result["target_bam_path"] = os.path.join(sample_output_dir, "target.bam")
        analysis_result["coverage_over_time"] = final_coverage_over_time
        
        analysis_result["processing_steps"].append("analysis_complete")
        logger.info(f"Multi-file target analysis completed for {sample_id}")
        logger.info(f"Files successfully processed: {analysis_result['files_processed']}/{analysis_result['total_files']}")
        logger.info(f"Final coverage: {analysis_result['coverage_data']['coverage']:.4f}")
        logger.info(f"Targets exceeding threshold: {analysis_result['coverage_data']['targets_exceeding_threshold']}")
        
        return analysis_result

    except Exception as e:
        logger.error(f"Error in multi-file target analysis for {sample_id}: {e}")
        analysis_result["error_message"] = str(e)
        analysis_result["processing_steps"].append("analysis_failed")
        return analysis_result


def target_handler(job, work_dir=None, reference=None, target_panel=None):
    """
    Handler function for target analysis jobs.
    This function processes files for target-specific analysis.

    Args:
        job: The workflow job containing file and metadata
        work_dir: Optional base directory for output (defaults to file directory)
        reference: Optional path to reference genome for SNP calling
        target_panel: Target panel type (required)
    """
    # Validate required parameters
    if not target_panel:
        raise ValueError("target_panel is required for target analysis")
    
    # Get job-specific logger
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
    
    # Check if this is a batched job
    batched_job = job.context.metadata.get("_batched_job")
    if batched_job:
        batch_size = batched_job.get_file_count()
        sample_id = batched_job.get_sample_id()
        batch_id = batched_job.batch_id
        logger.info(f"Processing target analysis batch: {batch_size} files for sample '{sample_id}' (batch_id: {batch_id})")
        
        # Get all filepaths in the batch
        filepaths = batched_job.get_filepaths()
        
        # Log individual files in the batch
        for i, filepath in enumerate(filepaths):
            logger.info(f"  Batch file {i+1}/{batch_size}: {os.path.basename(filepath)}")
        
        # Prepare metadata list for all BAM files in the batch (list comp inlined in 3.12)
        def _batch_metadata(i: int) -> dict:
            ctx = batched_job.contexts[i]
            sid = ctx.get_sample_id()
            return {**ctx.metadata.get("bam_metadata", {}), "sample_id": sid if sid != "unknown" else sample_id}
        metadata_list = [_batch_metadata(i) for i in range(len(filepaths))]

        # Determine work directory for the batch
        if work_dir is None:
            # Default to first BAM file directory
            batch_work_dir = os.path.dirname(filepaths[0])
        else:
            # Use specified work directory, create if it doesn't exist
            os.makedirs(work_dir, exist_ok=True)
            batch_work_dir = work_dir
            logger.debug(f"Using specified work directory: {batch_work_dir}")
        
        # Log and validate target panel
        job_panel = batched_job.contexts[0].metadata.get("target_panel", target_panel)
        if job_panel != target_panel:
            logger.warning(f"Panel mismatch: job metadata has '{job_panel}' but handler received '{target_panel}'. Using '{job_panel}' from metadata.")
            target_panel = job_panel
        logger.info(f"Using target panel: {target_panel}")
        
        # Debug: Log reference genome status
        if reference:
            logger.info(f"Reference genome provided to target_handler: {reference}")
        else:
            logger.info("No reference genome provided to target_handler")

        # Also check job metadata for reference
        job_reference = batched_job.contexts[0].metadata.get("reference")
        if job_reference:
            logger.info(f"Reference genome found in job metadata: {job_reference}")
            # Use job metadata reference if not provided directly
            if not reference:
                reference = job_reference
                logger.info(f"Using reference from job metadata: {reference}")
        else:
            logger.info("No reference genome found in job metadata")
        
        # Process all BAM files in the batch using the new aggregated function
        logger.info(f"Processing {batch_size} BAM files as aggregated batch for sample '{sample_id}'")
        batch_result = process_multiple_files(
            bam_paths=filepaths,
            metadata_list=metadata_list,
            work_dir=batch_work_dir,
            logger=logger,
            reference=reference,
            target_panel=target_panel
        )
        
        # Store batch results in job context (maintain compatibility with existing structure)
        job.context.add_metadata("target_analysis", {
            "batch_result": batch_result,  # Single aggregated result
            "batch_size": batch_size,
            "sample_id": sample_id,
            "batch_id": batch_id,
            "files_processed": batch_result.get("files_processed", batch_size),
            "total_files": batch_result.get("total_files", batch_size)
        })
        
        logger.info(f"Completed target analysis batch processing: {batch_size} files for sample '{sample_id}'")
        logger.info(f"Files successfully processed: {batch_result.get('files_processed', batch_size)}/{batch_result.get('total_files', batch_size)}")
        
        if batch_result.get("error_message"):
            logger.error(f"Batch processing completed with errors: {batch_result['error_message']}")
            job.context.add_error("target_analysis", batch_result["error_message"])
        else:
            logger.info("Batch processing completed successfully with aggregated target analysis")
            job.context.add_result(
                "target_analysis",
                {
                    "status": "success",
                    "sample_id": sample_id,
                    "analysis_time": batch_result.get("analysis_timestamp", 0),
                    "targets_found": batch_result.get("coverage_data", {}).get("targets_exceeding_threshold", 0),
                    "processing_steps": batch_result.get("processing_steps", []),
                    "target_data_path": batch_result.get("target_data_path", ""),
                    "target_plot_path": batch_result.get("target_plot_path", ""),
                    "target_bam_path": batch_result.get("target_bam_path", ""),
                    "coverage_summary": batch_result.get("coverage_data", {}),
                    "files_processed": batch_result.get("files_processed", batch_size),
                    "total_files": batch_result.get("total_files", batch_size),
                },
            )
        
        # Update master.csv with panel information
        if sample_id and not batch_result.get("error_message"):
            try:
                from robin.analysis.master_csv_manager import MasterCSVManager
                
                # Update master.csv with panel information
                csv_manager = MasterCSVManager(batch_work_dir)
                csv_manager.update_analysis_panel(sample_id, target_panel)
                logger.info(f"Updated master.csv with panel '{target_panel}' for sample {sample_id}")
                
            except Exception as e:
                logger.warning(f"Could not update master.csv with panel info for {sample_id}: {e}")
        
        return
        
    else:
        # Single file processing (backward compatibility)
        try:
            file_path = job.context.filepath

            logger.info(f"Starting target analysis for: {os.path.basename(file_path)}")
            
            # Log and validate target panel
            job_panel = job.context.metadata.get("target_panel", target_panel)
            if job_panel != target_panel:
                logger.warning(f"Panel mismatch: job metadata has '{job_panel}' but handler received '{target_panel}'. Using '{job_panel}' from metadata.")
                target_panel = job_panel
            logger.info(f"Using target panel: {target_panel}")

            # Debug: Log reference genome status
            if reference:
                logger.info(f"Reference genome provided to target_handler: {reference}")
            else:
                logger.info("No reference genome provided to target_handler")

            # Also check job metadata for reference
            job_reference = job.context.metadata.get("reference")
            if job_reference:
                logger.info(f"Reference genome found in job metadata: {job_reference}")
                # Use job metadata reference if not provided directly
                if not reference:
                    reference = job_reference
                    logger.info(f"Using reference from job metadata: {reference}")
            else:
                logger.info("No reference genome found in job metadata")

            # Get metadata from preprocessing
            file_metadata = job.context.metadata.get("bam_metadata", {})

            # Determine work directory
            if work_dir is None:
                # Default to file directory
                work_dir = os.path.dirname(file_path)
            else:
                # Use specified work directory, create if it doesn't exist
                os.makedirs(work_dir, exist_ok=True)

            # Initialize target analysis with staging enabled
            target_analysis = TargetAnalysis(
                work_dir=work_dir, 
                target_panel=target_panel,
                batch_size=10,  # Accumulate every 10 files
                use_staging=True
            )
            
            # Set reference genome if provided
            if reference:
                target_analysis.reference = reference
            
            # Use fast staging-based processing
            logger.info("Using staging-based processing (fast path)")
            target_metadata, should_accumulate = target_analysis.process_file_with_staging(
                file_path, file_metadata
            )
            
            # Convert TargetMetadata to dict for storage
            result = {
                "sample_id": target_metadata.sample_id,
                "file_path": target_metadata.file_path,
                "analysis_timestamp": target_metadata.analysis_timestamp,
                "processing_steps": target_metadata.processing_steps,
                "error_message": target_metadata.error_message,
                "coverage_data": target_metadata.coverage_data,
            }
            
            # Store results in job context
            job.context.add_metadata("target_analysis", result)
            
            # Trigger accumulation if threshold reached
            if should_accumulate:
                logger.info("Accumulation threshold reached - running batch accumulation")
                accumulation_result = target_analysis.accumulate_staged_files(
                    target_metadata.sample_id, force=False
                )
                logger.info(f"Accumulation result: {accumulation_result}")
                job.context.add_metadata("accumulation_result", accumulation_result)
            
            # Store flag for potential end-of-queue accumulation
            job.context.add_metadata("needs_final_accumulation", True)

            # Update master.csv with panel information
            if result.get("sample_id") and not result.get("error_message"):
                try:
                    from robin.analysis.master_csv_manager import MasterCSVManager
                    
                    # Determine work directory
                    if work_dir is None:
                        work_dir = os.path.dirname(file_path)
                    
                    # Update master.csv with panel information
                    csv_manager = MasterCSVManager(work_dir)
                    csv_manager.update_analysis_panel(result["sample_id"], target_panel)
                    logger.info(f"Updated master.csv with panel '{target_panel}' for sample {result['sample_id']}")
                    
                except Exception as e:
                    logger.warning(f"Could not update master.csv with panel info for {result.get('sample_id', 'unknown')}: {e}")

            if result.get("error_message"):
                job.context.add_error("target_analysis", result["error_message"])
                logger.error(f"Target analysis failed: {result['error_message']}")
            else:
                job.context.add_result(
                    "target_analysis",
                    {
                        "status": "success",
                        "sample_id": result.get("sample_id", "unknown"),
                        "analysis_time": result.get("analysis_timestamp", 0),
                        "targets_found": result.get("analysis_results", {}).get(
                            "targets_found", 0
                        ),
                        "processing_steps": result.get("processing_steps", []),
                        "target_data_path": result.get("target_data_path", ""),
                        "target_plot_path": result.get("target_plot_path", ""),
                        "target_bam_path": result.get("target_bam_path", ""),
                        "coverage_summary": result.get("analysis_results", {}).get(
                            "coverage_summary", {}
                        ),
                    },
                )
                logger.info(f"Target analysis complete for {os.path.basename(file_path)}")
                logger.info(f"Sample ID: {result.get('sample_id', 'unknown')}")
                logger.info(
                    f"Targets found: {result.get('analysis_results', {}).get('targets_found', 0)}"
                )

        except Exception as e:
            error_details = f"Error in target analysis for {job.context.filepath}: {str(e)}"
            job.context.add_error("target_analysis", error_details)
            logger.error(error_details)


def finalize_accumulation_for_sample(
    sample_id: str, work_dir: str, target_panel: str, reference: Optional[str] = None
) -> Dict[str, Any]:
    """
    Force final accumulation of any remaining staged files for a sample
    and merge all batch BAMs into final target.bam.
    
    This should be called when:
    - All files for a sample have been processed
    - The workflow is completing
    - There are staged files that haven't been accumulated yet
    
    Args:
        sample_id: Sample identifier
        work_dir: Working directory containing sample data
        target_panel: Target panel type
    
    Returns:
        Dictionary with accumulation results
    """
    logger = logging.getLogger("robin.target")
    
    try:
        logger.info(f"Finalizing accumulation for sample {sample_id}")
        
        sample_output_dir = os.path.join(work_dir, sample_id)
        
        # Initialize target analysis
        target_analysis = TargetAnalysis(
            work_dir=work_dir,
            target_panel=target_panel,
            batch_size=1,  # Force accumulation regardless of count
            use_staging=True
        )
        if reference:
            target_analysis.reference = reference
            logger.info(f"Finalize path using reference genome: {reference}")
        
        # Check if there are pending files
        pending_count = target_analysis._get_pending_count(sample_id)
        
        if pending_count > 0:
            logger.info(f"Found {pending_count} pending files for {sample_id} - forcing accumulation")
            # Force accumulation of remaining files
            result = target_analysis.accumulate_staged_files(sample_id, force=True)
            logger.info(f"Final accumulation complete for {sample_id}: {result}")
        else:
            logger.info(f"No pending files for {sample_id} - skipping accumulation")
            result = {"status": "no_pending_files", "sample_id": sample_id}
        
        # Now merge all batch BAMs into target.bam
        logger.info(f"Merging all batch BAMs into final target.bam for {sample_id}")
        
        target_bam_path = os.path.join(sample_output_dir, "target.bam")
        
        # Check if target.bam already exists but batch files weren't cleaned up from a previous run
        if os.path.exists(target_bam_path) and os.path.exists(f"{target_bam_path}.bai"):
            existing_batch_bams = sorted(glob.glob(os.path.join(sample_output_dir, "batch_*.bam")))
            if existing_batch_bams:
                logger.info(f"Found {len(existing_batch_bams)} leftover batch BAM files - cleaning up since target.bam already exists")
                cleaned_count = 0
                for batch_bam in existing_batch_bams:
                    try:
                        if os.path.exists(batch_bam):
                            os.remove(batch_bam)
                            cleaned_count += 1
                        if os.path.exists(f"{batch_bam}.bai"):
                            os.remove(f"{batch_bam}.bai")
                    except OSError as e:
                        logger.warning(f"Could not remove leftover batch BAM {os.path.basename(batch_bam)}: {e}")
                logger.info(f"Cleaned up {cleaned_count} leftover batch files")
        
        # Find all batch BAM files
        batch_bams = sorted(glob.glob(os.path.join(sample_output_dir, "batch_*.bam")))
        
        if batch_bams:
            logger.info(f"Found {len(batch_bams)} batch BAM files to merge")
            
            # Create temp merged output
            temp_merged_bam = os.path.join(sample_output_dir, ".final_merged.bam.tmp")
            
            # Clean stale temp outputs from a prior interrupted finalize, then force overwrite.
            if os.path.exists(temp_merged_bam):
                try:
                    os.remove(temp_merged_bam)
                except OSError as e:
                    logger.warning(f"Could not remove stale temp merged BAM {temp_merged_bam}: {e}")
            if os.path.exists(f"{temp_merged_bam}.bai"):
                try:
                    os.remove(f"{temp_merged_bam}.bai")
                except OSError as e:
                    logger.warning(f"Could not remove stale temp merged BAI {temp_merged_bam}.bai: {e}")

            # Merge all batch BAMs (force overwrite in case file appears between checks)
            pysam.merge("-f", "-o", temp_merged_bam, *batch_bams)
            logger.info("Merged all batch BAMs into temporary file")
            
            # Index the merged BAM
            pysam.index(temp_merged_bam)
            logger.info("Indexed merged target.bam")

            existing_target_read_count = 0
            existing_target_has_reads = False
            if os.path.exists(target_bam_path) and os.path.exists(f"{target_bam_path}.bai"):
                try:
                    with pysam.AlignmentFile(target_bam_path, "rb") as bam_file:
                        existing_target_read_count = bam_file.count(until_eof=True)
                        existing_target_has_reads = existing_target_read_count > 0
                except Exception as e:
                    logger.warning(f"Could not read existing target.bam before replacement: {e}")

            merged_read_count = 0
            try:
                with pysam.AlignmentFile(temp_merged_bam, "rb") as bam_file:
                    merged_read_count = bam_file.count(until_eof=True)
            except Exception as e:
                logger.warning(f"Could not read merged temp BAM before replacement: {e}")

            replaced_target_bam = False

            # Safety guard: preserve known-good existing BAM if newly merged BAM is empty.
            if merged_read_count == 0 and existing_target_has_reads:
                logger.warning(
                    "Merged BAM for sample %s has 0 reads; preserving existing target.bam with %s reads.",
                    sample_id,
                    existing_target_read_count,
                )
                try:
                    if os.path.exists(temp_merged_bam):
                        os.remove(temp_merged_bam)
                    if os.path.exists(f"{temp_merged_bam}.bai"):
                        os.remove(f"{temp_merged_bam}.bai")
                except OSError as e:
                    logger.warning(f"Could not clean temp merged BAM after preserve decision: {e}")

                target_bam_exists = os.path.exists(target_bam_path) and os.path.exists(f"{target_bam_path}.bai")
                result["final_merge"] = "preserved_existing"
                result["warning"] = (
                    "Merged BAM was empty; preserved existing target.bam and continued."
                )
            else:
                # Atomically replace target.bam
                if os.path.exists(target_bam_path):
                    os.remove(target_bam_path)
                    if os.path.exists(f"{target_bam_path}.bai"):
                        os.remove(f"{target_bam_path}.bai")

                shutil.move(temp_merged_bam, target_bam_path)
                if os.path.exists(f"{temp_merged_bam}.bai"):
                    shutil.move(f"{temp_merged_bam}.bai", f"{target_bam_path}.bai")
                replaced_target_bam = True

                # Verify target.bam was created successfully
                target_bam_exists = os.path.exists(target_bam_path) and os.path.exists(f"{target_bam_path}.bai")

                if target_bam_exists:
                    try:
                        with pysam.AlignmentFile(target_bam_path, "rb") as bam_file:
                            read_count = bam_file.count(until_eof=True)
                            if read_count > 0:
                                logger.info(
                                    f"Successfully created target.bam with {read_count} reads from {len(batch_bams)} batch files"
                                )
                            else:
                                logger.warning("target.bam created but contains no reads")
                    except Exception as e:
                        logger.warning(f"Could not verify target.bam: {e}")
                else:
                    logger.error("Failed to create target.bam file")
                result["final_merge"] = "success" if target_bam_exists else "failed"
            
            # Clean up batch BAM files and their associated BAI files
            # Only clean up if target.bam was successfully created
            if target_bam_exists and replaced_target_bam:
                logger.info(f"Cleaning up {len(batch_bams)} batch BAM files and their index files after final merge")
                cleaned_count = 0
                failed_count = 0
                
                for batch_bam in batch_bams:
                    try:
                        # Remove the batch BAM file
                        if os.path.exists(batch_bam):
                            os.remove(batch_bam)
                            cleaned_count += 1
                            logger.debug(f"Removed batch BAM: {os.path.basename(batch_bam)}")
                        
                        # Remove the associated BAI file
                        batch_bai = f"{batch_bam}.bai"
                        if os.path.exists(batch_bai):
                            os.remove(batch_bai)
                            logger.debug(f"Removed batch BAI: {os.path.basename(batch_bai)}")
                    except OSError as e:
                        failed_count += 1
                        logger.warning(f"Could not remove batch BAM {os.path.basename(batch_bam)}: {e}")
                
                logger.info(f"Batch cleanup complete: {cleaned_count} batch files removed, {failed_count} failures")
                
                if failed_count > 0:
                    logger.warning(f"Failed to remove {failed_count} batch file(s) - they may need manual cleanup")
            elif target_bam_exists:
                logger.info("Preserved existing target.bam; keeping batch files for investigation/retry.")
            else:
                logger.warning("Skipping batch cleanup - target.bam was not successfully created")
            
            logger.info(f"Final merge complete for {sample_id}")
            result["batch_files_merged"] = len(batch_bams)
        else:
            logger.info(f"No batch BAM files found for {sample_id}")
            result["final_merge"] = "no_batch_files"
        
        return result
        
    except Exception as e:
        logger.error(f"Error during final accumulation for {sample_id}: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return {"status": "error", "error": str(e), "sample_id": sample_id}


def target_bam_finalize_handler(job, work_dir: Optional[str] = None) -> None:
    """
    Workflow handler to finalize target BAMs for a sample.
    Runs accumulation and merges batch BAMs into target.bam.

    Required metadata:
    - sample_id: Sample identifier (optional; defaults to context sample ID)
    - work_dir: Parent directory containing sample folder (optional)
    - target_panel: Target panel type (optional; resolved from master.csv if missing)
    """
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)

    try:
        sid = (
            job.context.get_sample_id()
            if hasattr(job.context, "get_sample_id")
            else "unknown"
        )
        metadata = job.context.metadata or {}
        sample_id = metadata.get("sample_id") or sid
        sample_dir = metadata.get("sample_dir")
        base = work_dir or metadata.get("work_dir")

        if not sample_dir:
            if base and sample_id and sample_id != "unknown":
                sample_dir = os.path.join(base, sample_id)
            elif os.path.isdir(job.context.filepath):
                sample_dir = job.context.filepath
            else:
                sample_dir = os.path.dirname(job.context.filepath)

        if not base and sample_dir:
            base = os.path.dirname(sample_dir)

        if not base or not sample_id:
            raise RuntimeError(
                f"Missing work_dir/sample_id for target BAM finalization (work_dir={base}, sample_id={sample_id})"
            )

        target_panel = metadata.get("target_panel")
        reference = metadata.get("reference")
        if not target_panel and sample_dir and os.path.isdir(sample_dir):
            try:
                import csv
                master_csv = os.path.join(sample_dir, "master.csv")
                if os.path.exists(master_csv):
                    with open(master_csv, "r", newline="") as fh:
                        reader = csv.DictReader(fh)
                        first_row = next(reader, None)
                        if first_row:
                            panel = first_row.get("analysis_panel", "").strip()
                            if panel:
                                target_panel = panel
            except Exception:
                pass

        if not target_panel:
            target_panel = "rCNS2"

        logger.info(
            f"Finalizing target BAM: sample_id={sample_id}, work_dir={base}, target_panel={target_panel}, reference={reference}"
        )

        result = finalize_accumulation_for_sample(
            sample_id=sample_id, work_dir=base, target_panel=target_panel, reference=reference
        )

        if result.get("status") == "error":
            job.context.add_error("target_bam_finalize", result.get("error", "Unknown error"))
        else:
            job.context.add_result("target_bam_finalize", result)

    except Exception as e:
        job.context.add_error("target_bam_finalize", str(e))
        logger.error(f"Error in target BAM finalization handler: {e}")


def ensure_sorted_igv_bam(
    sample_dir: str, threads: int = 4, force_regenerate: bool = False
) -> str:
    """
    Ensure an IGV-ready coordinate-sorted and indexed BAM exists for a sample.

    Preferred output: ``sample_dir/igv/igv_ready.bam`` (+ .bai).
    Backward compat: if an existing Clair3-sorted BAM is present, reuse it.
    If neither exists, sort and index ``sample_dir/target.bam`` into
    ``sample_dir/igv/igv_ready.bam``.

    Args:
        sample_dir: Path to the sample output directory.
        threads: Number of threads to pass to pysam.sort.
        force_regenerate: If True, force recreation even if output exists.

    Returns:
        The path to the sorted BAM file (empty string on failure).
    """
    logger = logging.getLogger("robin.target")
    try:
        sample_dir = os.path.abspath(sample_dir)
        igv_dir = os.path.join(sample_dir, "igv")
        clair_dir = os.path.join(sample_dir, "clair3")
        os.makedirs(igv_dir, exist_ok=True)

        # Preferred IGV-ready path
        out_bam = os.path.join(igv_dir, "igv_ready.bam")
        out_bai = f"{out_bam}.bai"

        # Check if files exist and force_regenerate is not set
        if not force_regenerate and os.path.exists(out_bam) and os.path.exists(out_bai):
            logger.info(f"IGV BAM already present: {out_bam}")
            return out_bam

        # Backward compatibility: reuse Clair3-sorted BAM if present
        legacy_bams = [
            os.path.join(clair_dir, "sorted_targets_exceeding_rerun.bam"),
            os.path.join(clair_dir, "sorted_targets_exceeding.bam"),
        ]
        for legacy in legacy_bams:
            if os.path.exists(legacy) and os.path.exists(f"{legacy}.bai"):
                logger.info(f"Reusing legacy Clair3 BAM for IGV: {legacy}")
                # Optionally copy to IGV location to standardize
                try:
                    import shutil

                    shutil.copy2(legacy, out_bam)
                    shutil.copy2(f"{legacy}.bai", out_bai)
                    return out_bam
                except Exception:
                    # Fallback: just return legacy path
                    return legacy

        in_bam = os.path.join(sample_dir, "target.bam")
        if not os.path.exists(in_bam):
            logger.warning(
                f"target.bam not found in {sample_dir}; cannot build IGV BAM"
            )
            return ""

        logger.info(f"Sorting BAM for IGV: {in_bam} -> {out_bam}")
        pysam.sort(f"-@{threads}", "-o", out_bam, in_bam)
        pysam.index(out_bam, out_bai)
        logger.info(f"IGV BAM ready: {out_bam}")
        return out_bam
    except Exception as e:
        logger.error(f"Failed to create IGV BAM in {sample_dir}: {e}")
        return ""


def igv_bam_handler(job, work_dir: Optional[str] = None) -> None:
    """
    Workflow handler to build an IGV-ready BAM for a sample.

    Expects a per-sample directory containing `target.bam`. Outputs to
    `sample_dir/igv/igv_ready.bam` (and .bai). On success, records the output
    path in the job context under results['igv_bam'].
    """
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
    try:
        # Determine sample dir: prefer provided work_dir + sample_id
        sid = (
            job.context.get_sample_id()
            if hasattr(job.context, "get_sample_id")
            else "unknown"
        )
        base = work_dir or job.context.metadata.get("work_dir")
        sample_dir = None

        logger.info(f"IGV BAM job: sample_id={sid}, work_dir={work_dir}, base={base}")

        if base and sid and sid != "unknown":
            sample_dir = os.path.join(base, sid)
            logger.info(f"Using work_dir + sample_id path: {sample_dir}")
        else:
            # Fallback: directory of the current file
            sample_dir = os.path.dirname(job.context.filepath)
            logger.info(f"Using fallback path: {sample_dir}")

        # Additional fallback: if the filepath is already a sample directory, use it directly
        if os.path.isdir(job.context.filepath):
            sample_dir = job.context.filepath
            logger.info(f"Using filepath directly as sample directory: {sample_dir}")

        if not os.path.isdir(sample_dir):
            raise RuntimeError(f"Sample directory not found: {sample_dir}")

        logger.info(f"Final sample directory: {sample_dir}")
        logger.info(
            f"Looking for target.bam in: {os.path.join(sample_dir, 'target.bam')}"
        )

        # Check if force_regenerate is requested
        force_regenerate = job.context.metadata.get("force_regenerate", False)
        if force_regenerate:
            logger.info(
                "Force regenerate requested - will recreate IGV BAM even if it exists"
            )

        out = ensure_sorted_igv_bam(sample_dir, force_regenerate=force_regenerate)
        if not out:
            raise RuntimeError("Failed to create IGV-ready BAM")

        # Record result
        try:
            job.context.add_result("igv_bam", {"bam": out, "bai": f"{out}.bai"})
        except Exception:
            pass
        logger.info(f"IGV-ready BAM created: {out}")
    except Exception as e:
        job.context.add_error("igv_bam", str(e))
        logger.error(f"igv_bam handler failed: {e}")


def run_snp_analysis(
    sample_dir: str,
    threads: int = 4,
    force_regenerate: bool = False,
    reference: Optional[str] = None,
    annotation_only: bool = False,
    annotation_verbose: bool = False,
    target_panel: Optional[str] = None,
) -> str:
    """
    Run SNP analysis using Clair3 for a sample.

    This function runs the complete SNP calling pipeline including:
    - BAM sorting and preparation
    - Clair3 variant calling
    - snpEff annotation
    - SnpSift annotation with ClinVar data
    - Optional annotation-only reruns with verbose logging

    Args:
        sample_dir: Path to the sample output directory
        threads: Number of threads to use for processing
        force_regenerate: If True, force recreation even if output exists
        annotation_only: If True, run only the snpEff/SnpSift annotation steps
        annotation_verbose: If True, run annotation steps with verbose logging enabled
        target_panel: Target panel name used to resolve full target BED file

    Returns:
        The path to the output directory containing SNP results (empty string on failure)
    """
    logger = logging.getLogger("robin.target")

    # CRITICAL: Immediate logging to see if we even get here
    logger.info("ENTERING run_snp_analysis function")
    logger.info(
        "Parameters: sample_dir=%s, threads=%s, force_regenerate=%s, "
        "reference=%s, annotation_only=%s, annotation_verbose=%s, target_panel=%s",
        sample_dir,
        threads,
        force_regenerate,
        reference,
        annotation_only,
        annotation_verbose,
        target_panel,
    )

    # Also log to stderr to make sure we see it

    try:
        logger.info("STEP 1: Setting up sample directory")
        sample_dir = os.path.abspath(sample_dir)
        logger.info(f"Absolute sample_dir: {sample_dir}")

        clair_dir = os.path.join(sample_dir, "clair3")
        logger.info(f"Clair3 output directory: {clair_dir}")

        os.makedirs(clair_dir, exist_ok=True)
        logger.info(f"Clair3 directory created/verified: {clair_dir}")

        logger.info("STEP 2: Checking for existing SNP analysis")

        # Check if SNP analysis already exists and force_regenerate is not set
        snp_output_files = [
            os.path.join(clair_dir, "snpsift_output.vcf"),
            os.path.join(clair_dir, "snpsift_indel_output.vcf"),
        ]

        if annotation_only:
            logger.info("Annotation-only mode enabled; skipping existing output shortcut.")
        elif not force_regenerate and all(os.path.exists(f) for f in snp_output_files):
            logger.info(f"SNP analysis already present in {clair_dir}")
            logger.info("Rebuilding SNP display JSON from existing snpsift output.")
            try:
                snp_display_vcf = Path(clair_dir) / "snpsift_output.vcf"
                snp_display_path = Path(clair_dir) / "snpsift_output_display.json"
                if snp_display_vcf.exists():
                    snp_display = build_snp_display_data(snp_display_vcf)
                    if snp_display is not None:
                        with snp_display_path.open("w", encoding="utf-8") as f_out:
                            json.dump(snp_display, f_out)
                        logger.info(
                            f"SNP display data refreshed at {snp_display_path}"
                        )
                    else:
                        logger.warning(
                            "Could not regenerate SNP display data from existing VCF."
                        )
                else:
                    logger.warning(
                        "snpsift_output.vcf not found while refreshing display JSON."
                    )
            except Exception as display_exc:
                logger.warning(f"Failed to refresh SNP display JSON: {display_exc}")
            return clair_dir

        logger.info("STEP 3: Checking for required input files")
        target_bam = os.path.join(sample_dir, "target.bam")
        threshold_targets_bed = os.path.join(sample_dir, "targets_exceeding_threshold.bed")
        targets_bed = threshold_targets_bed
        output_snv_vcf = os.path.join(clair_dir, "output_done.vcf.gz")
        output_indel_vcf = os.path.join(clair_dir, "output_indel_done.vcf.gz")

        def resolve_full_targets_bed() -> str:
            # 1) Prefer target panel BED from resources when panel is known.
            if target_panel:
                if target_panel == "rCNS2":
                    bed_filename = "rCNS2_panel_name_uniq.bed"
                elif target_panel == "AML":
                    bed_filename = "AML_panel_name_uniq.bed"
                else:
                    bed_filename = f"{target_panel}_panel_name_uniq.bed"

                panel_paths = []
                if resources is not None:
                    try:
                        panel_paths.append(
                            os.path.join(
                                os.path.dirname(os.path.abspath(resources.__file__)),
                                bed_filename,
                            )
                        )
                    except Exception:
                        pass
                panel_paths.extend(
                    [bed_filename, f"data/{bed_filename}", f"/usr/local/share/{bed_filename}"]
                )
                for path in panel_paths:
                    if os.path.exists(path):
                        return path

            # 2) Fall back to sample-specific master BED if available.
            try:
                import glob

                bed_dir = os.path.join(sample_dir, "bed_files")
                if os.path.isdir(bed_dir):
                    master_bed_files = glob.glob(os.path.join(bed_dir, "master_*.bed"))
                    if master_bed_files:
                        return max(master_bed_files, key=os.path.getmtime)
            except Exception as e:
                logger.debug(f"Error resolving master BED from sample bed_files/: {e}")

            # 3) Final fallback: legacy threshold BED.
            return threshold_targets_bed

        logger.info("CHECKING INPUT FILES FOR CLAIR3")
        logger.info(f"  Sample directory: {sample_dir}")
        logger.info(f"  Target BAM path: {target_bam}")
        targets_bed = resolve_full_targets_bed()
        logger.info(f"  Threshold BED path: {threshold_targets_bed}")
        logger.info(f"  Targets BED path selected for ClairS: {targets_bed}")
        logger.info(f"  Reference path: {reference}")

        if annotation_only:
            missing_annotation_inputs = [
                path for path in [output_snv_vcf, output_indel_vcf] if not os.path.exists(path)
            ]
            if missing_annotation_inputs:
                for missing in missing_annotation_inputs:
                    logger.error(f"Annotation input not found: {missing}")
                logger.error(
                    "Annotation-only mode requires existing Clair3 outputs. "
                    "Run the full SNP pipeline first."
                )
                return ""
        else:
            if not os.path.exists(target_bam):
                logger.error(f"target.bam not found: {target_bam}")
                return ""

            if not os.path.exists(targets_bed):
                logger.error(f"Target BED file not found: {targets_bed}")
                return ""

        # Check for reference genome
        if not reference:
            # Try to find reference from common locations
            possible_refs = [
                os.path.join(sample_dir, "reference.fasta"),
                os.path.join(sample_dir, "reference.fa"),
                os.path.join(os.path.dirname(sample_dir), "reference.fasta"),
                os.path.join(os.path.dirname(sample_dir), "reference.fa"),
                "/usr/local/share/reference.fasta",
                "/usr/local/share/reference.fa",
            ]

            for ref_path in possible_refs:
                if os.path.exists(ref_path):
                    reference = ref_path
                    break

        if not reference:
            if annotation_only:
                logger.warning(
                    "Reference genome not found, continuing because annotation-only mode is enabled."
                )
            else:
                logger.error(
                    "Reference genome not found. SNP calling requires a reference genome."
                )
                return ""

        logger.info(f"Starting SNP analysis for sample: {os.path.basename(sample_dir)}")
        logger.info(f"Reference: {reference}")
        logger.info(f"Target BAM: {target_bam}")
        logger.info(f"Targets BED: {targets_bed}")

        if annotation_only:
            logger.info("Annotation-only mode enabled; skipping Clair3 variant calling.")
            sorted_bam = os.path.join(clair_dir, "sorted_targets_exceeding.bam")
        else:
            logger.info("STEP 4: Sorting target BAM for Clair3")
            sorted_bam = os.path.join(clair_dir, "sorted_targets_exceeding.bam")
            logger.info(f"Sorting BAM: {target_bam} -> {sorted_bam}")

            try:
                # Check if sorted BAM already exists and is valid
                if os.path.exists(sorted_bam) and os.path.exists(sorted_bam + ".bai"):
                    logger.info(f"Sorted BAM already exists: {sorted_bam}")
                    # Verify the file is readable
                    try:
                        test_bam = pysam.AlignmentFile(sorted_bam, "rb")
                        test_bam.close()
                        logger.info("Sorted BAM file is valid and readable")
                    except Exception as e:
                        logger.warning(
                            f"Existing sorted BAM appears corrupted, regenerating: {e}"
                        )
                        os.remove(sorted_bam)
                        if os.path.exists(sorted_bam + ".bai"):
                            os.remove(sorted_bam + ".bai")
                        raise Exception("Corrupted BAM file")
                else:
                    logger.info("Creating new sorted BAM file...")

                # Sort the BAM file
                logger.info(f"Sorting BAM with {threads} threads...")
                pysam.sort(f"-@{threads}", "-o", sorted_bam, target_bam)

                # Verify the sorted file was created
                if not os.path.exists(sorted_bam):
                    raise Exception(f"Sorted BAM file was not created: {sorted_bam}")

                # Create index
                logger.info("Creating BAM index...")
                pysam.index(sorted_bam)

                # Verify index was created
                if not os.path.exists(sorted_bam + ".bai"):
                    raise Exception(f"BAM index was not created: {sorted_bam}.bai")

                # Final verification
                file_size = os.path.getsize(sorted_bam)
                logger.info(
                    f"BAM sorting completed successfully: {sorted_bam} ({file_size / (1024**2):.1f} MB)"
                )

            except Exception as e:
                logger.error(f"BAM sorting failed: {e}")
                # Clean up any partial files
                for file_path in [sorted_bam, sorted_bam + ".bai"]:
                    if os.path.exists(file_path):
                        os.remove(file_path)
                        logger.info(f"Cleaned up partial file: {file_path}")
                return ""

            logger.info("STEP 5: Running Clair3 pipeline")
            logger.info("Running Clair3 pipeline...")

            # Run Clair3 pipeline using Docker
            logger.info("CHECKING DOCKER AVAILABILITY")
            if docker is None:
                logger.error("Docker module not available. Cannot run Clair3 pipeline.")
                return ""
            else:
                logger.info("Docker module is available")

            client = docker.from_env()

            # Test Docker connectivity
            try:
                client.ping()
                logger.info("Docker daemon is accessible")
            except Exception as docker_error:
                logger.error(f"Docker daemon not accessible: {docker_error}")
                return ""

            # Define container-specific paths using the original working pattern
            container_bamfile = f"/data/output/{os.path.basename(sorted_bam)}"  # BAM goes to output directory
            container_bedfile = f"/data/bed/{os.path.basename(targets_bed)}"  # BED goes to bed directory
            container_reference = f"/data/reference/{os.path.basename(reference)}"  # Reference goes to reference directory
            container_output = "/data/output"

            # Check if Clair3 image exists
            try:
                client.images.get("hkubal/clairs-to:latest")
                logger.info("Clair3 image found: hkubal/clairs-to:latest")
            except Exception:
                logger.warning(
                    "Clair3 image not found, attempting to pull: hkubal/clairs-to:latest"
                )
                try:
                    client.images.pull("hkubal/clairs-to:latest")
                    logger.info("Clair3 image pulled successfully")
                except Exception as pull_error:
                    logger.error(f"Failed to pull Clair3 image: {pull_error}")
                    return ""

            # Clean up old debugging - no longer needed with the fixed approach

            # Use the original working pattern: bind each directory to its own container path
            volume_bindings = {
                os.path.abspath(os.path.dirname(sorted_bam)): {
                    "bind": "/data/bam",
                    "mode": "ro",
                },  # BAM directory (clair3)
                os.path.abspath(os.path.dirname(targets_bed)): {
                    "bind": "/data/bed",
                    "mode": "ro",
                },  # BED directory (Sample_103)
                os.path.abspath(os.path.dirname(reference)): {
                    "bind": "/data/reference",
                    "mode": "ro",
                },  # Reference directory
                os.path.abspath(clair_dir): {
                    "bind": "/data/output",
                    "mode": "rw",
                },  # Output directory (clair3)
            }

            logger.info("Verifying volume bindings...")
            for source_path, binding_info in volume_bindings.items():
                if not os.path.exists(source_path):
                    logger.error(f"Volume source path does not exist: {source_path}")
                    return ""
                if not os.access(source_path, os.R_OK):
                    logger.error(f"Volume source path not readable: {source_path}")
                    return ""
                logger.info(
                    f"Volume binding verified: {source_path} -> {binding_info['bind']}"
                )

            # Verify the specific files exist in their directories
            logger.info("Verifying input files in volume directories...")
            
            if not os.path.exists(sorted_bam):
                logger.error(f"Sorted BAM file not found in directory: {sorted_bam}")
                return ""
            if not os.path.exists(targets_bed):
                logger.error(
                    f"Targets BED file not found in directory: {targets_bed}"
                )
                return ""
            if not os.path.exists(reference):
                logger.error(f"Reference file not found in directory: {reference}")
                return ""

            logger.info("All input files verified in their directories")

            host_config = client.api.create_host_config(binds=volume_bindings)

            # Function to split BED file into manageable chunks.
            # Chunks can include multiple chromosomes, constrained by covered genomic span.
            def split_bed_into_chunks(bed_file, max_chunk_size=150000000):
                """Split BED entries into chunks up to max_chunk_size covered genomic span."""
                chunks = []
                try:
                    # Read all BED entries and sort them
                    bed_entries = []
                    with open(bed_file, "r") as f:
                        for line_num, line in enumerate(f, 1):
                            raw_line = line.rstrip("\n")
                            line = raw_line.strip()
                            if not line or line.startswith("#"):
                                continue

                            parts = line.split("\t")
                            if len(parts) < 3:
                                logger.warning(
                                    f"Skipping invalid BED line {line_num}: {line}"
                                )
                                continue

                            chrom = parts[0]
                            start = int(parts[1])
                            end = int(parts[2])
                            if end <= start:
                                continue
                            entry_len = end - start
                            bed_entries.append((chrom, start, end, raw_line, entry_len))

                    # Sort by chromosome and start position
                    bed_entries.sort(key=lambda x: (x[0], x[1]))

                    def covered_span(chrom_bounds):
                        # Sum per-chromosome spans so sparse targets cannot inflate a chunk indefinitely.
                        return sum((end - start) for start, end in chrom_bounds.values())

                    current_chunk = []
                    current_bounds = {}

                    for chrom, start, end, raw_line, entry_len in bed_entries:
                        projected_bounds = dict(current_bounds)
                        if chrom in projected_bounds:
                            prev_start, prev_end = projected_bounds[chrom]
                            projected_bounds[chrom] = (min(prev_start, start), max(prev_end, end))
                        else:
                            projected_bounds[chrom] = (start, end)

                        if current_chunk and covered_span(projected_bounds) > max_chunk_size:
                            chunks.append(current_chunk)
                            current_chunk = []
                            current_bounds = {}

                        if chrom in current_bounds:
                            prev_start, prev_end = current_bounds[chrom]
                            current_bounds[chrom] = (min(prev_start, start), max(prev_end, end))
                        else:
                            current_bounds[chrom] = (start, end)
                        current_chunk.append((chrom, start, end, raw_line, entry_len))

                    if current_chunk:
                        chunks.append(current_chunk)

                    logger.info(
                        f"Combined BED entries into {len(chunks)} chunks (max chunk size: {max_chunk_size:,} bases)"
                    )
                    return chunks
                except Exception as e:
                    logger.error(f"Error splitting BED file into chunks: {e}")
                    return []

            # Default to region-split mode to keep INDEL calling memory usage bounded.
            # Set ROBIN_CLAIRS_SPLIT_REGIONS=0 to force single-pass mode.
            split_regions_env = os.environ.get("ROBIN_CLAIRS_SPLIT_REGIONS")
            if split_regions_env is None:
                use_split_regions = True
            else:
                use_split_regions = split_regions_env.lower() in ("1", "true", "yes", "on")

            if use_split_regions:
                logger.info(
                    "Using region-split ClairS mode (recommended for INDEL memory stability)."
                )
                chunked_entries = split_bed_into_chunks(targets_bed, max_chunk_size=250000000)
                if not chunked_entries:
                    logger.error("Failed to split BED file into chunks")
                    return ""
                chunk_bed_dir = os.path.join(clair_dir, "chunk_beds")
                os.makedirs(chunk_bed_dir, exist_ok=True)

                regions = []
                for i, chunk in enumerate(chunked_entries, 1):
                    chunk_host_bed = os.path.join(chunk_bed_dir, f"region_{i}.bed")
                    with open(chunk_host_bed, "w") as f_out:
                        for _, _, _, raw_line, _ in chunk:
                            f_out.write(raw_line + "\n")

                    first_chrom, first_start = chunk[0][0], chunk[0][1]
                    last_chrom, last_end = chunk[-1][0], chunk[-1][2]
                    chrom_bounds = {}
                    for chrom, start, end, _, _ in chunk:
                        if chrom in chrom_bounds:
                            prev_start, prev_end = chrom_bounds[chrom]
                            chrom_bounds[chrom] = (min(prev_start, start), max(prev_end, end))
                        else:
                            chrom_bounds[chrom] = (start, end)
                    span_bases = sum((end - start) for start, end in chrom_bounds.values())
                    label = f"{first_chrom}:{first_start+1}-{last_end}"
                    if first_chrom != last_chrom:
                        label = f"{first_chrom}:{first_start+1}..{last_chrom}:{last_end}"

                    regions.append(
                        {
                            "label": label,
                            "bed_container": f"{container_output}/chunk_beds/region_{i}.bed",
                            "is_split": True,
                            "span_bases": span_bases,
                        }
                    )
                    logger.info(
                        f"Chunk {i}: {label} (entries: {len(chunk)}, span: {span_bases:,} bases)"
                    )
            else:
                logger.info(
                    "Running ClairS in single-pass mode over full BED."
                )
                logger.info(
                    "Set ROBIN_CLAIRS_SPLIT_REGIONS=1 to re-enable region-splitting mode."
                )
                regions = [
                    {
                        "label": "full_target_bed",
                        "bed_container": container_bedfile,
                        "is_split": False,
                        "span_bases": None,
                    }
                ]

            # Process each region separately (or single pass)
            all_snv_vcfs = []
            all_indel_vcfs = []

            logger.info(f"Container input BAM: {container_bamfile}")
            logger.info(f"Container input BED: {container_bedfile}")
            logger.info(f"Container output dir: {container_output}")
            logger.info(f"Container reference: {container_reference}")
            logger.info(f"Volume bindings: {volume_bindings}")

            # Function to merge VCF files
            def merge_vcf_files(vcf_files, output_file, variant_type):
                """Merge multiple VCF files into a single output file"""
                try:
                    if not vcf_files:
                        logger.error(f"No {variant_type} VCF files to merge")
                        return False

                    logger.info(
                        f"Merging {len(vcf_files)} {variant_type} VCF files into {output_file}"
                    )

                    # Read headers from first file
                    import gzip

                    headers = []
                    first_vcf = vcf_files[0]

                    with gzip.open(first_vcf, "rt") as f:
                        for line in f:
                            if line.startswith("#"):
                                headers.append(line)
                            else:
                                break

                    # Write merged output
                    with gzip.open(output_file, "wt") as out_f:
                        # Write headers
                        for header in headers:
                            out_f.write(header)

                        # Write variants from each file
                        for vcf_file in vcf_files:
                            if os.path.exists(vcf_file):
                                with gzip.open(vcf_file, "rt") as in_f:
                                    for line in in_f:
                                        if not line.startswith("#"):
                                            out_f.write(line)

                    logger.info(
                        f"Successfully merged {variant_type} VCF files into {output_file}"
                    )
                    return True

                except Exception as e:
                    logger.error(f"Error merging {variant_type} VCF files: {e}")
                    return False

            # Log system resources before starting
            try:
                import psutil

                memory = psutil.virtual_memory()
                logger.info(
                    f"System memory: {memory.available / (1024**3):.1f} GB available out of {memory.total / (1024**3):.1f} GB"
                )
                logger.info(f"System CPU cores: {psutil.cpu_count()}")
            except ImportError:
                logger.info("psutil not available, skipping resource logging")

            # Process each region
            for i, region_spec in enumerate(regions):
                region_label = region_spec["label"]
                print(
                    f"Processing region {i+1}/{len(regions)}: {region_label} - PRINT STATEMENT"
                )
                logger.info(f"Processing region {i+1}/{len(regions)}: {region_label}")

                # Create region-specific output directory
                if region_spec["is_split"]:
                    region_output = f"{container_output}/region_{i+1}"
                else:
                    region_output = container_output

                # Build Clair3 command for this region
                command = (
                    f"/opt/bin/run_clairs_to "
                    f"--tumor_bam_fn {container_bamfile} "
                    f"--ref_fn {container_reference} "
                    f"--threads {threads} "
                    f"--remove_intermediate_dir "
                    f"--platform ont_r10_guppy_hac_5khz "
                    f"--output_dir {region_output} "
                    f"-b {region_spec['bed_container']} "
                    f"--chunk_size 5000000 "
                    f"--disable_intermediate_phasing "
                )

                logger.info(f"Running Clair3 command for region {i+1}: {command}")

                # Create and start container for this region
                logger.info(f"Creating Clair3 container for region {i+1}...")
                container = client.api.create_container(
                    image="hkubal/clairs-to:latest",
                    command=command,
                    volumes=list(volume_bindings.keys()),
                    host_config=host_config,
                )
                logger.info(f"Container created with ID: {container.get('Id')}")

                # Start container
                client.api.start(container=container.get("Id"))

                # Stream container logs
                for line in client.api.logs(
                    container=container.get("Id"), stream=True, follow=True
                ):
                    logger.info(f"Clair3 Region {i+1}: {line.decode().strip()}")

                # Wait for completion
                result = client.api.wait(container=container.get("Id"))
                if result["StatusCode"] != 0:
                    status_code = result.get("StatusCode")
                    container_id = container.get("Id")

                    # Classify abnormal termination (OOM/SIGKILL/SIGTERM) for clearer diagnostics.
                    oom_killed = False
                    kill_reason = None
                    try:
                        inspect_data = client.api.inspect_container(container=container_id)
                        state = inspect_data.get("State", {}) if inspect_data else {}
                        oom_killed = bool(state.get("OOMKilled", False))
                    except Exception as inspect_exc:
                        logger.debug(
                            f"Could not inspect Clair3 container state for region {i+1}: {inspect_exc}"
                        )

                    if oom_killed:
                        kill_reason = "oom_killed"
                    elif status_code == 137:
                        kill_reason = "sigkill_or_oom"
                    elif status_code == 143:
                        kill_reason = "sigterm"

                    # Get detailed error logs before removing container
                    error_logs = client.api.logs(
                        container=container_id, stderr=True, stdout=False
                    )
                    error_details = error_logs.decode().strip()
                    if error_details:
                        logger.error(f"Clair3 region {i+1} error logs: {error_details}")

                    # Also get stdout logs for context
                    stdout_logs = client.api.logs(
                        container=container_id, stdout=True, stderr=False
                    )
                    stdout_details = stdout_logs.decode().strip()
                    if stdout_details:
                        logger.info(
                            f"Clair3 region {i+1} stdout logs: {stdout_details}"
                        )

                    if kill_reason is not None:
                        logger.error(
                            "Clair3 region %s terminated abnormally (%s, exit=%s). "
                            "This is often memory pressure (especially for oom_killed/sigkill_or_oom). "
                            "Consider lowering threads or reducing chunk span.",
                            i + 1,
                            kill_reason,
                            status_code,
                        )

                    logger.warning(
                        f"Clair3 region {i+1} failed with status code {status_code}, skipping this region"
                    )
                    continue

                # Clean up container
                client.api.remove_container(container=container.get("Id"))

                # Check if output files were created for this region/pass
                if region_spec["is_split"]:
                    region_snv = f"{clair_dir}/region_{i+1}/snv.vcf.gz"
                    region_indel = f"{clair_dir}/region_{i+1}/indel.vcf.gz"
                else:
                    region_snv = f"{clair_dir}/snv.vcf.gz"
                    region_indel = f"{clair_dir}/indel.vcf.gz"

                if os.path.exists(region_snv):
                    all_snv_vcfs.append(region_snv)
                    logger.info(f"Region {i+1} SNV output created: {region_snv}")
                else:
                    logger.warning(f"Region {i+1} SNV output not found: {region_snv}")

                if os.path.exists(region_indel):
                    all_indel_vcfs.append(region_indel)
                    logger.info(f"Region {i+1} INDEL output created: {region_indel}")
                else:
                    logger.warning(
                        f"Region {i+1} INDEL output not found: {region_indel}"
                    )

            if use_split_regions:
                # Merge all split-region outputs
                logger.info(
                    f"Merging {len(all_snv_vcfs)} SNV VCF files and {len(all_indel_vcfs)} INDEL VCF files..."
                )

                if all_snv_vcfs:
                    merged_snv = f"{clair_dir}/merged_snv.vcf.gz"
                    merge_vcf_files(all_snv_vcfs, merged_snv, "SNV")
                    shutil.copy2(merged_snv, f"{clair_dir}/output_done.vcf.gz")
                    logger.info(f"Merged SNV VCFs into: {clair_dir}/output_done.vcf.gz")
                else:
                    logger.error("No SNV VCF files were successfully created")
                    return ""

                if all_indel_vcfs:
                    merged_indel = f"{clair_dir}/merged_indel.vcf.gz"
                    merge_vcf_files(all_indel_vcfs, merged_indel, "INDEL")
                    shutil.copy2(merged_indel, f"{clair_dir}/output_indel_done.vcf.gz")
                    logger.info(
                        f"Merged INDEL VCFs into: {clair_dir}/output_indel_done.vcf.gz"
                    )
                else:
                    logger.warning("No INDEL VCF files were successfully created")
                logger.info("Clair3 pipeline completed successfully for all regions")
            else:
                # Single-pass output normalization (no merge step needed).
                single_snv = f"{clair_dir}/snv.vcf.gz"
                single_indel = f"{clair_dir}/indel.vcf.gz"
                if not os.path.exists(single_snv):
                    logger.error(f"Single-pass SNV output not found: {single_snv}")
                    return ""
                shutil.copy2(single_snv, f"{clair_dir}/output_done.vcf.gz")
                logger.info(f"Single-pass SNV output copied to: {clair_dir}/output_done.vcf.gz")

                if os.path.exists(single_indel):
                    shutil.copy2(single_indel, f"{clair_dir}/output_indel_done.vcf.gz")
                    logger.info(
                        f"Single-pass INDEL output copied to: {clair_dir}/output_indel_done.vcf.gz"
                    )
                else:
                    logger.warning(f"Single-pass INDEL output not found: {single_indel}")
                logger.info("Clair3 pipeline completed successfully in single-pass mode")


        if annotation_only:
            logger.info(
                "Proceeding directly to annotation pipeline using existing Clair3 outputs."
            )

        logger.info("STEP 6: Running annotation pipeline")
        logger.info("Running annotation pipeline...")

        try:
            # SNP annotation with snpEff
            logger.info("Starting snpEff annotation for SNPs...")
            logger.info(f"Input file: {output_snv_vcf}")
            logger.info(f"Output file: {clair_dir}/snpeff_output.vcf")

            # Define output file paths
            snpeff_out = f"{clair_dir}/snpeff_output.vcf"
            snpsift_out = f"{clair_dir}/snpsift_output.vcf"
            snpeff_indel_out = f"{clair_dir}/snpeff_indel_output.vcf"
            snpsift_indel_out = f"{clair_dir}/snpsift_indel_output.vcf"

            # Check if input file exists and has content
            if os.path.exists(output_snv_vcf):
                file_size = os.path.getsize(output_snv_vcf)
                logger.info(f"Input VCF file exists, size: {file_size} bytes")

                # Check first few lines to verify it's a valid VCF
                try:
                    import gzip

                    with gzip.open(output_snv_vcf, "rt") as f:
                        first_lines = [next(f) for _ in range(5)]
                    logger.info("First 5 lines of input VCF:")
                    for i, line in enumerate(first_lines):
                        logger.info(f"  Line {i+1}: {line.strip()}")
                except Exception as e:
                    logger.warning(f"Could not read input VCF file: {e}")
            else:
                logger.error(
                    f"Input VCF file does not exist: {output_snv_vcf}"
                )

            snpeff_cmd = ["snpEff"]
            snpeff_cmd.append("-v" if annotation_verbose else "-q")
            snpeff_cmd.extend(["hg38", output_snv_vcf])
            logger.info(f"Running snpEff command: {' '.join(snpeff_cmd)}")

            annotation_env = os.environ.copy()
            annotation_env["_JAVA_OPTIONS"] = "-Xmx16g"

            # Check if snpEff is available
            try:
                snpeff_version = subprocess.run(
                    ["snpEff", "-version"], capture_output=True, text=True
                )
                if snpeff_version.returncode == 0:
                    logger.info(f"snpEff version: {snpeff_version.stdout.strip()}")
                else:
                    logger.warning(
                        f"Could not get snpEff version: {snpeff_version.stderr}"
                    )
            except FileNotFoundError:
                logger.error("snpEff command not found in PATH")

            with open(snpeff_out, "w") as fout:
                logger.info("Starting snpEff execution...")
                result = subprocess.run(
                    snpeff_cmd,
                    stdout=fout,
                    stderr=subprocess.PIPE,
                    text=True,
                    env=annotation_env,
                )
                logger.info(
                    f"snpEff execution completed with return code: {result.returncode}"
                )
                if annotation_verbose and result.stderr:
                    logger.info(f"snpEff stderr output:\n{result.stderr}")

            if result.returncode != 0:
                logger.warning(f"snpEff failed with return code: {result.returncode}")
                logger.warning(f"snpEff stderr output: {result.stderr}")
                logger.info("Using raw Clair3 output for SNP annotation")
                # Decompress the gzipped file when copying to .vcf extension
                import gzip

                with gzip.open(output_snv_vcf, "rt") as f_in:
                    with open(snpeff_out, "w") as f_out:
                        f_out.write(f_in.read())
                logger.info("Created uncompressed VCF from gzipped input")
            else:
                logger.info("snpEff completed successfully!")
                # Check the output file
                if os.path.exists(snpeff_out):
                    output_size = os.path.getsize(snpeff_out)
                    logger.info(
                        f"snpEff output file created, size: {output_size} bytes"
                    )

                    # Check first few lines of output to verify annotation
                    try:
                        with open(snpeff_out, "r") as f:
                            first_lines = [next(f) for _ in range(5)]
                        logger.info("First 5 lines of snpEff output:")
                        for i, line in enumerate(first_lines):
                            logger.info(f"  Line {i+1}: {line.strip()}")
                    except Exception as e:
                        logger.warning(f"Could not read snpEff output: {e}")
                else:
                    logger.error("snpEff output file was not created!")

            # SnpSift annotation
            logger.info("Starting SnpSift annotation for SNPs...")
            logger.info(f"Input file: {snpeff_out}")
            logger.info(f"Output file: {snpsift_out}")

            # Check if snpEff output exists and has content
            if os.path.exists(snpeff_out):
                snpeff_size = os.path.getsize(snpeff_out)
                logger.info(f"snpEff output file exists, size: {snpeff_size} bytes")
            else:
                logger.error(f"snpEff output file does not exist: {snpeff_out}")

            clinvar_db_path = _resolve_clinvar_db_for_snpsift(logger)
            try:
                if clinvar_db_path:
                    clinvar_size = os.path.getsize(clinvar_db_path)
                    logger.info(
                        f"ClinVar DB found, size: {clinvar_size} bytes"
                    )

                    # Check if SnpSift is available
                    try:
                        snpsift_version = subprocess.run(
                            ["SnpSift"], capture_output=True, text=True
                        )
                        version_output = "\n".join(
                            part.strip()
                            for part in (
                                snpsift_version.stdout,
                                snpsift_version.stderr,
                            )
                            if part
                        )
                        version_line = ""
                        for line in version_output.splitlines():
                            if "SnpSift version" in line:
                                version_line = line.strip()
                                break
                        if version_line:
                            logger.info(f"SnpSift version: {version_line}")
                        elif snpsift_version.returncode == 0:
                            logger.info("SnpSift command executed successfully")
                        else:
                            logger.warning(
                                "Could not determine SnpSift version from command output"
                            )
                    except FileNotFoundError:
                        logger.error("SnpSift command not found in PATH")

                    snpsift_cmd = ["SnpSift", "annotate"]
                    if annotation_verbose:
                        snpsift_cmd.append("-v")
                    snpsift_cmd.extend([clinvar_db_path, snpeff_out])
                    logger.info(f"Running SnpSift command: {' '.join(snpsift_cmd)}")

                    snpsift_env = annotation_env.copy()

                    with open(snpsift_out, "w") as fout:
                        logger.info("Starting SnpSift execution...")
                        result = subprocess.run(
                            snpsift_cmd,
                            stdout=fout,
                            stderr=subprocess.PIPE,
                            text=True,
                            env=snpsift_env,
                        )
                        logger.info(
                            f"SnpSift execution completed with return code: {result.returncode}"
                        )
                        if annotation_verbose and result.stderr:
                            logger.info(f"SnpSift stderr output:\n{result.stderr}")

                    if result.returncode != 0:
                        logger.warning(
                            f"SnpSift failed with return code: {result.returncode}"
                        )
                        logger.warning(f"SnpSift stderr output: {result.stderr}")
                        logger.info("Using snpEff output for final SNP annotation")
                        shutil.copy2(snpeff_out, snpsift_out)
                        logger.info("Copied snpEff output as final SNP annotation")
                    else:
                        logger.info("SnpSift completed successfully!")
                        # Check the output file
                        snpsift_out = f"{clair_dir}/snpsift_output.vcf"
                        if os.path.exists(snpsift_out):
                            output_size = os.path.getsize(snpsift_out)
                            logger.info(
                                f"SnpSift output file created, size: {output_size} bytes"
                            )

                            # Check first few lines of output to verify annotation
                            try:
                                with open(snpsift_out, "r") as f:
                                    first_lines = [next(f) for _ in range(5)]
                                logger.info("First 5 lines of SnpSift output:")
                                for i, line in enumerate(first_lines):
                                    logger.info(f"  Line {i+1}: {line.strip()}")
                            except Exception as e:
                                logger.warning(f"Could not read SnpSift output: {e}")
                        else:
                            logger.error("SnpSift output file was not created!")
                else:
                    logger.warning(
                        "ClinVar bgzip/tabix database unavailable, using snpEff output"
                    )
                    shutil.copy2(snpeff_out, snpsift_out)
                    logger.info(
                        "Copied snpEff output as final SNP annotation (no ClinVar)"
                    )

            except Exception as e:
                logger.warning(f"SnpSift annotation failed with exception: {e}")
                logger.info("Using snpEff output for final SNP annotation")
                shutil.copy2(snpeff_out, snpsift_out)
                logger.info(
                    "Copied snpEff output as final SNP annotation due to exception"
                )

            # INDEL annotation
            logger.info("Starting snpEff annotation for INDELs...")
            logger.info(f"Input file: {output_indel_vcf}")
            logger.info(f"Output file: {snpeff_indel_out}")

            # Check if input file exists and has content
            if os.path.exists(output_indel_vcf):
                file_size = os.path.getsize(output_indel_vcf)
                logger.info(f"INDEL input VCF file exists, size: {file_size} bytes")

                # Check first few lines to verify it's a valid VCF
                try:
                    import gzip

                    with gzip.open(output_indel_vcf, "rt") as f:
                        first_lines = [next(f) for _ in range(5)]
                    logger.info("First 5 lines of INDEL input VCF:")
                    for i, line in enumerate(first_lines):
                        logger.info(f"  Line {i+1}: {line.strip()}")
                except Exception as e:
                    logger.warning(f"Could not read INDEL input VCF file: {e}")
            else:
                logger.error(
                    f"INDEL input VCF file does not exist: {output_indel_vcf}"
                )

            snpeff_indel_cmd = ["snpEff"]
            snpeff_indel_cmd.append("-v" if annotation_verbose else "-q")
            snpeff_indel_cmd.extend(["hg38", output_indel_vcf])
            logger.info(f"Running snpEff INDEL command: {' '.join(snpeff_indel_cmd)}")

            with open(snpeff_indel_out, "w") as fout:
                logger.info("Starting snpEff INDEL execution...")
                result = subprocess.run(
                    snpeff_indel_cmd,
                    stdout=fout,
                    stderr=subprocess.PIPE,
                    text=True,
                    env=annotation_env,
                )
                logger.info(
                    f"snpEff INDEL execution completed with return code: {result.returncode}"
                )
                if annotation_verbose and result.stderr:
                    logger.info(f"snpEff INDEL stderr output:\n{result.stderr}")

            if result.returncode != 0:
                logger.warning(
                    f"snpEff (INDEL) failed with return code: {result.returncode}"
                )
                logger.warning(f"snpEff (INDEL) stderr output: {result.stderr}")
                logger.info("Using raw Clair3 INDEL output for annotation")
                # Decompress the gzipped file when copying to .vcf extension
                import gzip

                with gzip.open(output_indel_vcf, "rt") as f_in:
                    with open(snpeff_indel_out, "w") as f_out:
                        f_out.write(f_in.read())
                logger.info("Created uncompressed INDEL VCF from gzipped input")
            else:
                logger.info("snpEff INDEL completed successfully!")
                # Check the output file
                if os.path.exists(snpeff_indel_out):
                    output_size = os.path.getsize(snpeff_indel_out)
                    logger.info(
                        f"snpEff INDEL output file created, size: {output_size} bytes"
                    )

                    # Check first few lines of output to verify annotation
                    try:
                        with open(snpeff_indel_out, "r") as f:
                            first_lines = [next(f) for _ in range(5)]
                        logger.info("First 5 lines of snpEff INDEL output:")
                        for i, line in enumerate(first_lines):
                            logger.info(f"  Line {i+1}: {line.strip()}")
                    except Exception as e:
                        logger.warning(f"Could not read snpEff INDEL output: {e}")
                else:
                    logger.error("snpEff INDEL output file was not created!")

            # SnpSift annotation for INDELs
            logger.info("Starting SnpSift annotation for INDELs...")
            logger.info(f"Input file: {snpeff_indel_out}")
            logger.info(f"Output file: {snpsift_indel_out}")

            # Check if snpEff INDEL output exists and has content
            if os.path.exists(snpeff_indel_out):
                snpeff_indel_size = os.path.getsize(snpeff_indel_out)
                logger.info(
                    f"snpEff INDEL output file exists, size: {snpeff_indel_size} bytes"
                )
            else:
                logger.error(
                    f"snpEff INDEL output file does not exist: {snpeff_indel_out}"
                )

            try:
                if clinvar_db_path:
                    logger.info(
                        f"Using ClinVar file for INDEL annotation: {clinvar_db_path}"
                    )

                    snpsift_indel_cmd = ["SnpSift", "annotate"]
                    if annotation_verbose:
                        snpsift_indel_cmd.append("-v")
                    snpsift_indel_cmd.extend([clinvar_db_path, snpeff_indel_out])
                    logger.info(
                        f"Running SnpSift INDEL command: {' '.join(snpsift_indel_cmd)}"
                    )

                    snpsift_indel_env = annotation_env.copy()

                    with open(snpsift_indel_out, "w") as fout:
                        logger.info("Starting SnpSift INDEL execution...")
                        result = subprocess.run(
                            snpsift_indel_cmd,
                            stdout=fout,
                            stderr=subprocess.PIPE,
                            text=True,
                            env=snpsift_indel_env,
                        )
                        logger.info(
                            f"SnpSift INDEL execution completed with return code: {result.returncode}"
                        )
                        if annotation_verbose and result.stderr:
                            logger.info(f"SnpSift INDEL stderr output:\n{result.stderr}")

                    if result.returncode != 0:
                        logger.warning(
                            f"SnpSift (INDEL) failed with return code: {result.returncode}"
                        )
                        logger.warning(
                            f"SnpSift (INDEL) stderr output: {result.stderr}"
                        )
                        logger.info(
                            "Using snpEff INDEL output for final INDEL annotation"
                        )
                        shutil.copy2(snpeff_indel_out, snpsift_indel_out)
                        logger.info(
                            "Copied snpEff INDEL output as final INDEL annotation"
                        )
                    else:
                        logger.info("SnpSift INDEL completed successfully!")
                        # Check the output file
                        if os.path.exists(snpsift_indel_out):
                            output_size = os.path.getsize(snpsift_indel_out)
                            logger.info(
                                f"SnpSift INDEL output file created, size: {output_size} bytes"
                            )

                            # Check first few lines of output to verify annotation
                            try:
                                with open(snpsift_indel_out, "r") as f:
                                    first_lines = [next(f) for _ in range(5)]
                                logger.info("First 5 lines of SnpSift INDEL output:")
                                for i, line in enumerate(first_lines):
                                    logger.info(f"  Line {i+1}: {line.strip()}")
                            except Exception as e:
                                logger.warning(
                                    f"Could not read SnpSift INDEL output: {e}"
                                )
                        else:
                            logger.error("SnpSift INDEL output file was not created!")
                else:
                    logger.warning(
                        "ClinVar file not found for INDELs, using snpEff INDEL output"
                    )
                    shutil.copy2(snpeff_indel_out, snpsift_indel_out)
                    logger.info(
                        "Copied snpEff INDEL output as final INDEL annotation (no ClinVar)"
                    )

            except Exception as e:
                logger.warning(f"SnpSift INDEL annotation failed with exception: {e}")
                logger.info("Using snpEff INDEL output for final INDEL annotation")
                shutil.copy2(snpeff_indel_out, snpsift_indel_out)
                logger.info(
                    "Copied snpEff INDEL output as final INDEL annotation due to exception"
                )

            logger.info("Annotation pipeline completed successfully")

        except Exception as e:
            logger.error(f"Error in annotation pipeline: {e}")
            return ""

        # Parse VCF files to CSV format
        logger.info("Parsing VCF files...")
        try:
            # Import parse_vcf function or implement basic parsing
            def basic_vcf_to_csv(vcf_file, csv_file):
                """Basic VCF to CSV conversion"""
                try:
                    with open(vcf_file, "r") as f:
                        lines = [line.strip() for line in f if not line.startswith("#")]

                    if not lines:
                        return

                    # Parse header from first non-comment line
                    header = lines[0].split("\t")

                    # Parse data
                    data = []
                    for line in lines[1:]:
                        fields = line.split("\t")
                        if len(fields) >= len(header):
                            data.append(fields[: len(header)])

                    # Create DataFrame and save
                    df = pd.DataFrame(data, columns=header)
                    df.to_csv(csv_file, index=False)
                    logger.info(f"VCF parsed to CSV: {csv_file}")

                except Exception as e:
                    logger.warning(f"Basic VCF parsing failed for {vcf_file}: {e}")

            # Parse both SNP and INDEL VCFs
            basic_vcf_to_csv(
                f"{clair_dir}/snpsift_output.vcf", f"{clair_dir}/snpsift_output.vcf.csv"
            )
            basic_vcf_to_csv(
                f"{clair_dir}/snpsift_indel_output.vcf",
                f"{clair_dir}/snpsift_indel_output.vcf.csv",
            )

            # Build pre-formatted SNP display data for the GUI
            try:
                snp_display_path = Path(clair_dir) / "snpsift_output_display.json"
                snp_display = build_snp_display_data(Path(clair_dir) / "snpsift_output.vcf")
                if snp_display is not None:
                    with snp_display_path.open("w", encoding="utf-8") as f_out:
                        json.dump(snp_display, f_out)
                    logger.info(f"SNP display data written to {snp_display_path}")
                else:
                    logger.warning("Could not generate SNP display data from VCF.")
            except Exception as display_exc:
                logger.warning(f"Failed to build SNP display data: {display_exc}")

            # Note: VCF files are now properly uncompressed when they have .vcf extensions
            # Compressed .vcf.gz files are only used as input to processing tools

        except Exception as e:
            logger.warning(f"VCF parsing failed: {e}")

        logger.info("FINAL SUCCESS: SNP analysis completed successfully!")
        logger.info(
            f"SNP analysis completed successfully for {os.path.basename(sample_dir)}"
        )
        logger.info(f"Output directory: {clair_dir}")
        return clair_dir

    except Exception as e:
        print(f"EXCEPTION in run_snp_analysis: {e} - PRINT STATEMENT")
        logger.error(f"Failed to run SNP analysis in {sample_dir}: {e}")
        import traceback

        print(f"TRACEBACK: {traceback.format_exc()} - PRINT STATEMENT")
        return ""


def snp_analysis_handler(job, work_dir: Optional[str] = None) -> None:
    """
    Workflow handler to run SNP analysis using Clair3 for a sample.

    Expects a per-sample directory containing `target.bam` and `targets_exceeding_threshold.bed`.
    Outputs to `sample_dir/clair3/` with annotated VCF files and CSV summaries.
    On success, records the output directory in the job context under results['snp_analysis'].

    Required metadata:
    - reference: Path to reference genome (optional, will auto-detect)
    - threads: Number of threads to use (default: 4)
    - force_regenerate: Whether to force regeneration of existing results (default: False)
    - annotation_only: Run only the annotation steps (default: False)
    - annotation_verbose: Enable verbose logging for snpEff/SnpSift (default: False)
    """
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)

    try:
        # Determine sample dir: prefer provided work_dir + sample_id
        sid = (
            job.context.get_sample_id()
            if hasattr(job.context, "get_sample_id")
            else "unknown"
        )
        base = work_dir or job.context.metadata.get("work_dir")
        sample_dir = None

        logger.info(
            f"SNP analysis job: sample_id={sid}, work_dir={work_dir}, base={base}"
        )

        if base and sid and sid != "unknown":
            sample_dir = os.path.join(base, sid)
            logger.info(f"Using work_dir + sample_id path: {sample_dir}")
        else:
            # Fallback: directory of the current file
            sample_dir = os.path.dirname(job.context.filepath)
            logger.info(f"Using fallback path: {sample_dir}")

        # Additional fallback: if the filepath is already a sample directory, use it directly
        if os.path.isdir(job.context.filepath):
            sample_dir = job.context.filepath
            logger.info(f"Using filepath directly as sample directory: {sample_dir}")

        if not os.path.isdir(sample_dir):
            raise RuntimeError(f"Sample directory not found: {sample_dir}")

        logger.info(f"Final sample directory: {sample_dir}")

        # Let's see what's actually in this directory
        try:
            if os.path.exists(sample_dir):
                dir_contents = os.listdir(sample_dir)
                logger.info(f"Directory contents of {sample_dir}:")
                for item in dir_contents:
                    item_path = os.path.join(sample_dir, item)
                    if os.path.isfile(item_path):
                        size = os.path.getsize(item_path)
                        logger.info(f"  {item} ({size / (1024**2):.1f} MB)")
                    else:
                        logger.info(f"  {item}/")
            else:
                logger.error(f"Sample directory does not exist: {sample_dir}")
        except Exception as e:
            logger.warning(f"Could not list directory contents: {e}")

        # Check for required input files
        target_bam = os.path.join(sample_dir, "target.bam")
        targets_bed = os.path.join(sample_dir, "targets_exceeding_threshold.bed")

        if not os.path.exists(target_bam):
            raise RuntimeError(f"target.bam not found: {target_bam}")

        if not os.path.exists(targets_bed):
            raise RuntimeError(
                f"targets_exceeding_threshold.bed not found: {targets_bed}"
            )

        # Get job parameters
        threads = job.context.metadata.get("threads", 4)
        force_regenerate = job.context.metadata.get("force_regenerate", False)
        reference = job.context.metadata.get("reference")
        annotation_only = job.context.metadata.get("annotation_only", False)
        annotation_verbose = job.context.metadata.get("annotation_verbose", False)
        target_panel = job.context.metadata.get("target_panel")
        if not target_panel:
            try:
                master_csv = os.path.join(sample_dir, "master.csv")
                if os.path.exists(master_csv):
                    import csv

                    with open(master_csv, "r", newline="") as fh:
                        reader = csv.DictReader(fh)
                        first_row = next(reader, None)
                        if first_row:
                            panel = first_row.get("analysis_panel", "").strip()
                            if panel:
                                target_panel = panel
            except Exception:
                pass

        # Debug logging for reference genome
        logger.info("=== SNP ANALYSIS HANDLER DEBUGGING ===")
        logger.info(f"Job metadata keys: {list(job.context.metadata.keys())}")
        logger.info(f"Job metadata: {job.context.metadata}")
        logger.info(f"Extracted reference: {reference}")
        logger.info(f"Annotation-only: {annotation_only}")
        logger.info(f"Annotation-verbose: {annotation_verbose}")
        logger.info(f"Target panel: {target_panel}")
        logger.info("=== END SNP ANALYSIS HANDLER DEBUGGING ===")

        if force_regenerate:
            logger.info(
                "Force regenerate requested - will recreate SNP analysis even if it exists"
            )

        if reference:
            logger.info(f"Using reference genome from job metadata: {reference}")
        else:
            logger.warning(
                "No reference genome in job metadata - will try to auto-detect"
            )

        # Run SNP analysis
        logger.info("ABOUT TO CALL run_snp_analysis function")
        logger.info(
            "Calling run_snp_analysis with: sample_dir=%s, threads=%s, "
            "force_regenerate=%s, reference=%s, annotation_only=%s, annotation_verbose=%s, target_panel=%s"
            % (
                sample_dir,
                threads,
                force_regenerate,
                reference,
                annotation_only,
                annotation_verbose,
                target_panel,
            )
        )

        output_dir = run_snp_analysis(
            sample_dir,
            threads=threads,
            force_regenerate=force_regenerate,
            reference=reference,
            annotation_only=annotation_only,
            annotation_verbose=annotation_verbose,
            target_panel=target_panel,
        )

        logger.info(f"run_snp_analysis returned: {output_dir}")

        if not output_dir:
            raise RuntimeError("Failed to complete SNP analysis")

        # Record result
        try:
            job.context.add_result(
                "snp_analysis",
                {
                    "output_dir": output_dir,
                    "snps_vcf": os.path.join(output_dir, "snpsift_output.vcf"),
                    "indels_vcf": os.path.join(output_dir, "snpsift_indel_output.vcf"),
                    "snps_csv": os.path.join(output_dir, "snpsift_output.vcf.csv"),
                    "indels_csv": os.path.join(
                        output_dir, "snpsift_indel_output.vcf.csv"
                    ),
                    "sample_id": sid,
                    "analysis_timestamp": time.time(),
                },
            )
        except Exception:
            pass

        logger.info(f"SNP analysis completed successfully for {sid}")
        logger.info(f"Output directory: {output_dir}")

    except Exception as e:
        job.context.add_error("snp_analysis", str(e))
        logger.error(f"snp_analysis handler failed: {e}")
