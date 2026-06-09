#!/usr/bin/env python3
"""
Master BED File Generator for ROBIN

This module generates merged master BED files from CNV and fusion breakpoint BED files.
The master BED file combines:
- CNV regions (new_file_{counter}.bed) - split into start/end breakpoints, padded by +/- 10kb
- CNV breakpoints (breakpoints_{counter}.bed) - merged as-is
- Fusion breakpoints (fusion_breakpoints_{counter}.bed) - merged as-is

All regions are intersected with the target gene panel and merged to remove overlaps.
"""

import os
import sys
import logging
import glob
import fcntl
import time
import argparse
import json
import threading
from typing import List, Tuple, Dict, Optional, Any
from pathlib import Path

import pandas as pd
import numpy as np

logger = logging.getLogger("robin.analysis.master_bed")


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
        try:
            while True:
                try:
                    fcntl.flock(self.fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    return
                except OSError:
                    if time.time() - start_time > self.timeout:
                        self.fd.close()
                        self.fd = None
                        raise TimeoutError(
                            f"Could not acquire lock {self.lock_file} within {self.timeout}s"
                        )
                    time.sleep(0.25)
        except Exception:
            if self.fd is not None:
                try:
                    self.fd.close()
                except Exception:
                    pass
                self.fd = None
            raise
    
    def release(self):
        """Release the lock"""
        if self.fd:
            fcntl.flock(self.fd, fcntl.LOCK_UN)
            self.fd.close()
            self.fd = None


def _get_latest_bed_file(bed_dir: str, pattern: str) -> Optional[str]:
    """
    Find the latest BED file matching a pattern based on counter.
    
    Args:
        bed_dir: Directory containing BED files
        pattern: Glob pattern to match (e.g., "new_file_*.bed")
        
    Returns:
        Path to the latest BED file, or None if not found
    """
    bed_files = glob.glob(os.path.join(bed_dir, pattern))
    if not bed_files:
        return None
    
    # Extract counter from filename and find the latest
    latest_file = None
    latest_counter = -1
    
    for bed_file in bed_files:
        try:
            # Extract counter from filename (e.g., "new_file_001.bed" -> 1)
            basename = os.path.basename(bed_file)
            # Pattern: name_counter.bed
            parts = basename.rsplit("_", 1)
            if len(parts) == 2:
                counter_str = parts[1].replace(".bed", "")
                counter = int(counter_str)
                if counter > latest_counter:
                    latest_counter = counter
                    latest_file = bed_file
        except (ValueError, IndexError):
            # If we can't parse the counter, use modification time as fallback
            if latest_file is None or os.path.getmtime(bed_file) > os.path.getmtime(latest_file):
                latest_file = bed_file
    
    return latest_file


def _load_bed_file(bed_path: str, require_bed6: bool = False) -> pd.DataFrame:
    """
    Load a BED file into a DataFrame, preserving strand information.
    
    Args:
        bed_path: Path to BED file
        require_bed6: If True, expects BED6 format (chrom, start, end, name, score, strand)
        
    Returns:
        DataFrame with columns: chrom, start, end, name, score, strand
    """
    if not os.path.exists(bed_path):
        return pd.DataFrame()
    
    try:
        # Read BED file - handle variable number of columns (BED format can have 3-12 columns)
        # First, read without column names to detect actual number of columns
        df_temp = pd.read_csv(
            bed_path,
            sep="\t",
            header=None,
            nrows=1,  # Just read first row to check column count
        )
        num_cols = len(df_temp.columns)
        
        # Standard BED6 column names
        bed6_cols = ["chrom", "start", "end", "name", "score", "strand"]
        
        # Read full file with appropriate column names based on actual column count
        if num_cols >= 6:
            # File has at least 6 columns (BED6 format)
            col_names = bed6_cols + [f"col{i}" for i in range(6, num_cols)]
            df = pd.read_csv(
                bed_path,
                sep="\t",
                header=None,
                names=col_names[:num_cols],
                dtype={"chrom": str, "start": int, "end": int, "name": str, "score": str, "strand": str},
                na_values=["."],
            )
            # Keep BED6 columns
            df = df[bed6_cols].copy()
        elif num_cols >= 4:
            # File has at least 4 columns (chrom, start, end, name)
            col_names = ["chrom", "start", "end", "name"] + [f"col{i}" for i in range(4, num_cols)]
            df = pd.read_csv(
                bed_path,
                sep="\t",
                header=None,
                names=col_names[:num_cols],
                dtype={"chrom": str, "start": int, "end": int, "name": str},
                na_values=["."],
            )
            # Add missing BED6 columns
            df["score"] = "0"
            df["strand"] = "."
            df = df[bed6_cols].copy()
        elif num_cols == 3:
            # File has only 3 columns (chrom, start, end)
            df = pd.read_csv(
                bed_path,
                sep="\t",
                header=None,
                names=["chrom", "start", "end"],
                dtype={"chrom": str, "start": int, "end": int},
            )
            # Add missing BED6 columns
            df["name"] = "."
            df["score"] = "0"
            df["strand"] = "."
            df = df[bed6_cols].copy()
        else:
            logger.warning(f"BED file has too few columns ({num_cols}): {bed_path}")
            return pd.DataFrame()
        
        # Fill missing values
        df["name"] = df["name"].fillna(".")
        df["score"] = df["score"].fillna("0")
        df["strand"] = df["strand"].fillna(".")
        
        # Normalize strand values
        df["strand"] = df["strand"].replace({"": ".", "?": "."})
        
        return df[bed6_cols]
    except Exception as e:
        logger.warning(f"Error loading BED file {bed_path}: {e}")
        return pd.DataFrame()


def _expand_cnv_regions(df: pd.DataFrame, expand_size: int = 10000) -> pd.DataFrame:
    """
    Convert CNV regions into breakpoint regions at start and end positions.
    Each breakpoint is padded by +/- expand_size and duplicated for both strands.
    
    Processing:
    1. Split each CNV region into two breakpoints: one at start, one at end
    2. Pad each breakpoint by +/- expand_size (creates 2*expand_size wide regions)
    3. Duplicate each breakpoint for both + and - strands
    
    Args:
        df: DataFrame with chrom, start, end, name, score, strand columns
        expand_size: Size to pad around each breakpoint in base pairs (default 10kb)
        
    Returns:
        DataFrame with breakpoint regions, duplicated for both strands
    """
    if df.empty:
        return df
    
    breakpoint_regions = []
    
    for _, row in df.iterrows():
        chrom = row["chrom"]
        region_start = row["start"]
        region_end = row["end"]
        name = row.get("name", ".")
        score = row.get("score", "0")
        
        # Create start breakpoint: pad around region_start
        start_bp_start = max(0, region_start - expand_size)
        start_bp_end = region_start + expand_size
        
        # Create end breakpoint: pad around region_end
        end_bp_start = max(0, region_end - expand_size)
        end_bp_end = region_end + expand_size
        
        # Add start breakpoint for both strands
        breakpoint_regions.append({
            "chrom": chrom,
            "start": start_bp_start,
            "end": start_bp_end,
            "name": name,
            "score": score,
            "strand": "+"
        })
        breakpoint_regions.append({
            "chrom": chrom,
            "start": start_bp_start,
            "end": start_bp_end,
            "name": name,
            "score": score,
            "strand": "-"
        })
        
        # Add end breakpoint for both strands
        breakpoint_regions.append({
            "chrom": chrom,
            "start": end_bp_start,
            "end": end_bp_end,
            "name": name,
            "score": score,
            "strand": "+"
        })
        breakpoint_regions.append({
            "chrom": chrom,
            "start": end_bp_start,
            "end": end_bp_end,
            "name": name,
            "score": score,
            "strand": "-"
        })
    
    result = pd.DataFrame(breakpoint_regions)
    return result


def _merge_overlapping_regions(df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge overlapping regions in a BED DataFrame, preserving strand information.
    Regions are merged separately by chromosome and strand.
    
    Args:
        df: DataFrame with chrom, start, end, name, score, strand columns
        
    Returns:
        DataFrame with merged regions
    """
    if df.empty:
        return df
    
    merged_regions = []
    
    # Sort by chromosome, strand, then start position
    df_sorted = df.sort_values(by=["chrom", "strand", "start", "end"]).copy()
    
    # Group by chromosome and strand
    for (chrom, strand), group in df_sorted.groupby(["chrom", "strand"], observed=True):
        if group.empty:
            continue
        
        # Merge overlapping regions within chromosome and strand
        current_start = None
        current_end = None
        current_name = None
        
        for _, row in group.iterrows():
            start = row["start"]
            end = row["end"]
            name = row.get("name", "merged")
            
            if current_start is None:
                # First region
                current_start = start
                current_end = end
                current_name = name
            elif start <= current_end:
                # Overlaps with current region - extend it
                current_end = max(current_end, end)
                # Combine names if different
                if current_name != name and name != "merged":
                    if current_name == "merged":
                        current_name = name
                    elif name not in current_name:
                        current_name = f"{current_name},{name}"
            else:
                # No overlap - save current region and start new one
                merged_regions.append({
                    "chrom": chrom,
                    "start": current_start,
                    "end": current_end,
                    "name": current_name if current_name else "merged",
                    "score": row.get("score", "0"),
                    "strand": strand
                })
                current_start = start
                current_end = end
                current_name = name
        
        # Don't forget the last region
        if current_start is not None:
            merged_regions.append({
                "chrom": chrom,
                "start": current_start,
                "end": current_end,
                "name": current_name if current_name else "merged",
                "score": "0",
                "strand": strand
            })
    
    if not merged_regions:
        return pd.DataFrame()
    
    return pd.DataFrame(merged_regions)


def _duplicate_unstranded_regions(df: pd.DataFrame) -> pd.DataFrame:
    """
    Duplicate regions with unstranded ('.') strand to both + and - strands.
    
    Args:
        df: DataFrame with chrom, start, end, name, score, strand columns
        
    Returns:
        DataFrame with unstranded regions duplicated for both strands
    """
    if df.empty:
        return df
    
    # Separate stranded and unstranded regions
    stranded = df[df["strand"].isin(["+", "-"])].copy()
    unstranded = df[df["strand"] == "."].copy()
    
    if unstranded.empty:
        return stranded
    
    # Duplicate unstranded regions for both strands
    plus_strand = unstranded.copy()
    plus_strand["strand"] = "+"
    
    minus_strand = unstranded.copy()
    minus_strand["strand"] = "-"
    
    # Combine stranded regions with duplicated unstranded regions
    result = pd.concat([stranded, plus_strand, minus_strand], ignore_index=True)
    
    return result


def _intersect_with_target_genes(
    regions_df: pd.DataFrame,
    target_bed_path: str
) -> pd.DataFrame:
    """
    Filter regions to only include those that overlap with target gene BED file.
    Keeps the original regions that overlap, preserving strand information.
    
    Args:
        regions_df: DataFrame with chrom, start, end, name, score, strand columns
        target_bed_path: Path to target gene BED file
        
    Returns:
        DataFrame with regions that overlap target genes
    """
    if regions_df.empty:
        return pd.DataFrame()
    
    if not os.path.exists(target_bed_path):
        logger.warning(f"Target BED file not found: {target_bed_path}, skipping intersection")
        return regions_df
    
    try:
        # Load target genes (BED6 format)
        target_df = _load_bed_file(target_bed_path, require_bed6=True)
        if target_df.empty:
            logger.warning(f"Target BED file is empty: {target_bed_path}")
            return pd.DataFrame()
        
        # Duplicate unstranded target regions for both strands
        target_df = _duplicate_unstranded_regions(target_df)
        
        # Filter regions to only those that overlap with target genes
        overlapping_regions = []
        
        for _, region_row in regions_df.iterrows():
            region_chrom = region_row["chrom"]
            region_start = region_row["start"]
            region_end = region_row["end"]
            region_strand = region_row.get("strand", ".")
            
            # Check if this region overlaps with any target gene (same chromosome)
            overlapping_targets = target_df[
                (target_df["chrom"] == region_chrom) &
                (target_df["start"] < region_end) &
                (target_df["end"] > region_start)
            ]
            
            if not overlapping_targets.empty:
                # Keep the original region (it overlaps with target genes)
                # Preserve the region's strand information
                overlapping_regions.append({
                    "chrom": region_chrom,
                    "start": region_start,
                    "end": region_end,
                    "name": region_row.get("name", "overlap"),
                    "score": region_row.get("score", "0"),
                    "strand": region_strand
                })
        
        if not overlapping_regions:
            return pd.DataFrame()
        
        return pd.DataFrame(overlapping_regions)
        
    except Exception as e:
        logger.error(f"Error intersecting with target genes: {e}")
        return regions_df


def _load_fai_file(fai_path: str) -> Dict[str, int]:
    """
    Load a FASTA index (.fai) file and return chromosome lengths.
    
    Args:
        fai_path: Path to .fai file
        
    Returns:
        Dictionary mapping chromosome names to lengths
    """
    chrom_lengths = {}
    
    if not os.path.exists(fai_path):
        logger.warning(f"FAI file not found: {fai_path}")
        return chrom_lengths
    
    try:
        with open(fai_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                # FAI format: chrom_name    length    offset    linebases    linewidth
                parts = line.split('\t')
                if len(parts) >= 2:
                    chrom_name = parts[0]
                    chrom_length = int(parts[1])
                    chrom_lengths[chrom_name] = chrom_length
    except Exception as e:
        logger.error(f"Error loading FAI file {fai_path}: {e}")
    
    return chrom_lengths


def _calculate_genome_coverage(
    bed_df: pd.DataFrame,
    fai_path: Optional[str] = None
) -> Optional[Dict[str, float]]:
    """
    Calculate the proportion of the genome covered by BED regions, accounting for strand.
    
    Args:
        bed_df: DataFrame with chrom, start, end, strand columns
        fai_path: Optional path to .fai file for genome sizes
        
    Returns:
        Dictionary with coverage statistics, or None if calculation fails
    """
    if bed_df.empty:
        return None
    
    if not fai_path or not os.path.exists(fai_path):
        logger.warning("FAI file not provided or not found - cannot calculate genome coverage")
        return None
    
    # Load chromosome lengths
    chrom_lengths = _load_fai_file(fai_path)
    if not chrom_lengths:
        logger.warning("No chromosome lengths loaded from FAI file")
        return None
    
    # Calculate total genome size (sum of all chromosome lengths)
    total_genome_size = sum(chrom_lengths.values())
    if total_genome_size == 0:
        logger.warning("Total genome size is 0")
        return None
    
    # Check for chromosomes in BED file that aren't in FAI file
    bed_chroms = set(bed_df["chrom"].unique())
    fai_chroms = set(chrom_lengths.keys())
    missing_chroms = bed_chroms - fai_chroms
    if missing_chroms:
        logger.warning(f"Chromosomes in BED file not found in FAI file: {sorted(missing_chroms)}")
    
    # Calculate coverage for each strand separately
    strand_coverage = {"+": 0, "-": 0, ".": 0}
    strand_total_genome = {"+": total_genome_size, "-": total_genome_size, ".": total_genome_size}
    
    # Group by strand
    for strand in ["+", "-", "."]:
        strand_df = bed_df[bed_df["strand"] == strand].copy()
        if strand_df.empty:
            continue
        
        # Calculate total covered bases for this strand
        # IMPORTANT: Regions must be merged before calling this function to avoid double-counting overlaps
        covered_bases = 0
        skipped_regions = 0
        for _, row in strand_df.iterrows():
            chrom = row["chrom"]
            start = row["start"]
            end = row["end"]
            
            # Only count if chromosome exists in FAI file
            if chrom in chrom_lengths:
                # Clip region to chromosome boundaries
                region_length = max(0, min(end, chrom_lengths[chrom]) - max(start, 0))
                covered_bases += region_length
            else:
                skipped_regions += 1
        
        if skipped_regions > 0:
            logger.debug(f"Skipped {skipped_regions} {strand} strand regions due to missing chromosomes in FAI file")
        
        strand_coverage[strand] = covered_bases
    
    # Calculate proportions
    results = {
        "total_genome_size": total_genome_size,
        "plus_strand_coverage": strand_coverage["+"],
        "minus_strand_coverage": strand_coverage["-"],
        "unstranded_coverage": strand_coverage["."],
        "plus_strand_proportion": strand_coverage["+"] / strand_total_genome["+"] if strand_total_genome["+"] > 0 else 0.0,
        "minus_strand_proportion": strand_coverage["-"] / strand_total_genome["-"] if strand_total_genome["-"] > 0 else 0.0,
        "unstranded_proportion": strand_coverage["."] / strand_total_genome["."] if strand_total_genome["."] > 0 else 0.0,
    }
    
    # Calculate combined coverage
    # + and - strands are treated separately (each can have different coverage)
    # Unstranded regions (if any) are counted once
    total_covered = strand_coverage["+"] + strand_coverage["-"] + strand_coverage["."]
    # Total possible coverage: 2 * genome_size (for + and - strands separately)
    # plus unstranded regions if they exist
    total_possible_coverage = 2 * total_genome_size
    results["total_covered_bases"] = total_covered
    results["total_proportion"] = total_covered / total_possible_coverage if total_possible_coverage > 0 else 0.0
    
    return results


def _sort_bed_regions(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sort BED regions by chromosome (numeric order), strand, then start, then end.
    
    Args:
        df: DataFrame with chrom, start, end, strand columns
        
    Returns:
        Sorted DataFrame
    """
    if df.empty:
        return df
    
    def chrom_sort_key(chrom):
        """Sort chromosomes: chr1-22, then chrX, chrY, chrM, then others"""
        if not chrom.startswith("chr"):
            return (300, chrom)  # Non-standard chromosomes at end
        
        chrom_suffix = chrom[3:]
        if chrom_suffix == "X":
            return (100, 0)
        elif chrom_suffix == "Y":
            return (101, 0)
        elif chrom_suffix == "M":
            return (102, 0)
        else:
            try:
                return (int(chrom_suffix), 0)
            except ValueError:
                return (300, chrom_suffix)
    
    def strand_sort_key(strand):
        """Sort strands: +, then -, then ."""
        if strand == "+":
            return 0
        elif strand == "-":
            return 1
        else:
            return 2
    
    df = df.copy()
    df["_chrom_sort_key"] = df["chrom"].apply(chrom_sort_key)
    df["_strand_sort_key"] = df.get("strand", ".").apply(strand_sort_key)
    df = df.sort_values(by=["_chrom_sort_key", "_strand_sort_key", "start", "end"])
    df = df.drop(columns=["_chrom_sort_key", "_strand_sort_key"])
    
    return df.reset_index(drop=True)


def _get_reference_path() -> Optional[str]:
    """
    Get the reference genome path from config or environment variable.
    
    Returns:
        Path to reference FASTA file (expanded), or None if not found
    """
    # Try to load from config file first
    try:
        import yaml  # type: ignore
        config_paths = [
            os.path.join(os.getcwd(), "config.yaml"),
            os.path.expanduser("~/.robin/config.yaml"),
            "/etc/robin/config.yaml",
        ]
        for config_path in config_paths:
            if os.path.exists(config_path):
                with open(config_path, 'r') as f:
                    config = yaml.safe_load(f)
                    if config and isinstance(config, dict):
                        reference = config.get("reference")
                        if reference:
                            # Expand user home directory if present
                            reference = os.path.expanduser(reference)
                            if os.path.exists(reference):
                                return reference
    except (ImportError, Exception):
        pass
    
    # Check environment variable
    reference_env = os.environ.get("robin_REFERENCE")
    if reference_env:
        # Expand user home directory if present
        reference_env = os.path.expanduser(reference_env)
        if os.path.exists(reference_env):
            return reference_env
    
    return None


def _find_fai_file(reference: Optional[str] = None) -> Optional[str]:
    """
    Find the FAI (FASTA index) file path from reference genome or environment variable.
    If reference exists but FAI doesn't, attempts to create it using samtools faidx.
    
    Args:
        reference: Optional path to reference FASTA file
        
    Returns:
        Path to FAI file, or None if not found
    """
    # Get reference if not provided
    if not reference:
        reference = _get_reference_path()
    
    # If reference is provided, expand user home directory and check for corresponding .fai file
    if reference:
        # Expand user home directory if present (handles ~/path)
        reference = os.path.expanduser(reference)
        fai_path = f"{reference}.fai"
        if os.path.exists(fai_path):
            return fai_path
        
        # If reference exists but FAI doesn't, try to create it
        if os.path.exists(reference):
            try:
                import subprocess
                logger.debug(f"FAI file not found for {reference}, attempting to create it with samtools faidx")
                result = subprocess.run(
                    ["samtools", "faidx", reference],
                    capture_output=True,
                    text=True,
                    timeout=300
                )
                if result.returncode == 0 and os.path.exists(fai_path):
                    logger.info(f"Successfully created FAI index: {fai_path}")
                    return fai_path
                else:
                    logger.warning(f"Failed to create FAI index for {reference}: {result.stderr}")
            except (subprocess.TimeoutExpired, FileNotFoundError, Exception) as e:
                logger.debug(f"Could not create FAI index for {reference}: {e}")
        else:
            logger.debug(f"Reference file does not exist: {reference}")
    
    # Check environment variable directly (fallback)
    reference_env = os.environ.get("robin_REFERENCE")
    if reference_env:
        # Expand user home directory if present
        reference_env = os.path.expanduser(reference_env)
        fai_path = f"{reference_env}.fai"
        if os.path.exists(fai_path):
            return fai_path
    
    # Common locations to check
    common_paths = [
        "/data/reference/genome.fa.fai",
        "/usr/local/share/reference/genome.fa.fai",
        os.path.expanduser("~/reference/genome.fa.fai"),
    ]
    
    for path in common_paths:
        if os.path.exists(path):
            return path
    
    return None


def _log_bed_coverage_data(
    sample_id: str,
    work_dir: str,
    analysis_counter: int,
    coverage_data: Dict[str, Dict[str, float]],
    log_file: str = "bed_coverage_log.json"
) -> None:
    """
    Log genome coverage data for each BED file type to a JSON file.
    Appends to existing log if it exists.
    
    Args:
        sample_id: Sample ID
        work_dir: Working directory
        analysis_counter: Analysis counter
        coverage_data: Dictionary mapping BED file type to coverage statistics
        log_file: Name of the log file (default: bed_coverage_log.json)
    """
    try:
        sample_dir = os.path.join(work_dir, sample_id)
        log_path = os.path.join(sample_dir, log_file)
        
        # Load existing log if it exists
        existing_data = []
        if os.path.exists(log_path):
            try:
                with open(log_path, 'r') as f:
                    existing_data = json.load(f)
                if not isinstance(existing_data, list):
                    # If file is corrupted or in wrong format, start fresh
                    existing_data = []
            except (json.JSONDecodeError, IOError) as e:
                logger.warning(f"Could not read existing coverage log: {e}")
                existing_data = []
        
        # Create new log entry
        log_entry = {
            "timestamp": time.time(),
            "analysis_counter": analysis_counter,
            "coverage": coverage_data
        }
        
        # Append new entry
        existing_data.append(log_entry)
        
        # Write updated log
        with open(log_path, 'w') as f:
            json.dump(existing_data, f, indent=2)
        
        logger.debug(f"Logged BED coverage data to {log_path}")
        
        # Generate pre-processed visualization data
        _generate_visualization_data(sample_dir, existing_data)
        
    except Exception as e:
        logger.warning(f"Could not log BED coverage data: {e}")


def _generate_visualization_data(sample_dir: str, log_entries: List[Dict[str, Any]]) -> None:
    """
    Generate pre-processed visualization data for the GUI.
    This avoids doing heavy DataFrame operations in the GUI for each viewer.
    
    Args:
        sample_dir: Sample directory path
        log_entries: List of log entry dictionaries
    """
    try:
        # pandas is already imported at module level
        
        if not log_entries:
            return
        
        # Color map for different BED types
        colors = {
            "cnv_regions": "#1f77b4",
            "cnv_breakpoints": "#ff7f0e",
            "fusion_breakpoints": "#2ca02c",
            "master_bed_breakpoints": "#d62728",
            "target_panel": "#9467bd",
            "master_bed": "#8c564b",
        }
        
        # Prepare DataFrame from log entries
        rows = []
        for entry in log_entries:
            timestamp = entry.get("timestamp", 0)
            counter = entry.get("analysis_counter", 0)
            coverage = entry.get("coverage", {})
            
            # Extract coverage for each BED file type
            for bed_type, bed_data in coverage.items():
                if isinstance(bed_data, dict):
                    total_prop = bed_data.get("total_proportion") if "total_proportion" in bed_data else None
                    rows.append({
                        "timestamp": timestamp,
                        "analysis_counter": counter,
                        "bed_type": bed_type,
                        "total_proportion": total_prop,
                    })
        
        if not rows:
            return
        
        df = pd.DataFrame(rows)
        
        # Format data for ECharts
        series = []
        for bed_type in sorted(df["bed_type"].unique()):
            bed_df = df[df["bed_type"] == bed_type].sort_values("timestamp")
            if bed_df.empty:
                continue
            
            # Format data for ECharts time series: [[timestamp_ms, value], ...]
            data = []
            for _, row in bed_df.iterrows():
                total_prop = row["total_proportion"]
                # Only include if value is not None and not NaN
                if total_prop is not None and not pd.isna(total_prop):
                    timestamp_ms = int(row["timestamp"] * 1000)  # Convert to milliseconds
                    value = float(total_prop * 100)  # Convert to percentage
                    data.append([timestamp_ms, value])
            
            # Only add series if we have data points to plot
            if not data:
                continue
            
            color = colors.get(bed_type, "#7f7f7f")
            series.append({
                "name": bed_type.replace("_", " ").title(),
                "type": "line",
                "smooth": True,
                "symbol": "circle",
                "symbolSize": 6,
                "data": data,
                "itemStyle": {"color": color},
                "lineStyle": {"width": 2, "color": color},
            })
        
        # Write pre-processed visualization data
        viz_path = os.path.join(sample_dir, "bed_coverage_viz.json")
        with open(viz_path, 'w') as f:
            json.dump({"series": series}, f, indent=2)
        
        logger.debug(f"Generated visualization data with {len(series)} series to {viz_path}")
        
    except Exception as e:
        logger.warning(f"Could not generate visualization data: {e}")


def _get_target_bed_path(target_panel: Optional[str] = None) -> Optional[str]:
    """
    Get the path to the target gene BED file.
    
    Args:
        target_panel: Target panel name (rCNS2, AML, custom, etc.)
        
    Returns:
        Path to target BED file, or None if not found
    """
    try:
        from robin import resources
        resources_dir = os.path.dirname(resources.__file__)
        
        if target_panel == "rCNS2":
            bed_path = os.path.join(resources_dir, "rCNS2_panel_name_uniq.bed")
        elif target_panel == "AML":
            bed_path = os.path.join(resources_dir, "AML_panel_name_uniq.bed")
        else:
            # Try custom panel
            bed_path = os.path.join(resources_dir, f"{target_panel}_panel_name_uniq.bed")
        
        if os.path.exists(bed_path):
            return bed_path
        
        # Fallback paths
        if target_panel == "rCNS2":
            fallback = "rCNS2_panel_name_uniq.bed"
        elif target_panel == "AML":
            fallback = "AML_panel_name_uniq.bed"
        else:
            fallback = f"{target_panel}_panel_name_uniq.bed"
        
        if os.path.exists(fallback):
            return fallback
            
    except ImportError:
        pass
    
    return None


def _try_get_target_panel_from_fusion_metadata(
    sample_id: str,
    work_dir: str,
) -> Optional[str]:
    """
    Try to get target_panel from fusion metadata if available.
    
    Args:
        sample_id: Sample ID
        work_dir: Working directory
        
    Returns:
        Target panel name if found, None otherwise
    """
    try:
        from robin.analysis.fusion_work import _load_fusion_metadata
        
        fusion_metadata = _load_fusion_metadata(work_dir, sample_id)
        if fusion_metadata and hasattr(fusion_metadata, 'target_panel'):
            return fusion_metadata.target_panel
    except Exception:
        pass
    
    return None


def _build_master_bed_data(
    sample_id: str,
    work_dir: str,
    analysis_counter: int,
    target_panel: Optional[str],
    reference: Optional[str],
    log: logging.Logger,
) -> Optional[Tuple[pd.DataFrame, Dict[str, Dict[str, float]]]]:
    """
    Build master BED DataFrame and coverage data from all BED sources (no lock).
    Caller must hold lock only when writing the result to disk.
    """
    sample_dir = os.path.join(work_dir, sample_id)
    bed_dir = os.path.join(sample_dir, "bed_files")
    if not os.path.exists(bed_dir):
        return None

    ref = reference
    if not ref:
        ref = _get_reference_path()
        if ref:
            ref = os.path.expanduser(ref)
            log.debug(f"Using reference genome from config/environment: {ref}")
    elif reference:
        ref = os.path.expanduser(reference)
    fai_path = _find_fai_file(ref)
    if fai_path:
        log.debug(f"Using FAI file for coverage calculations: {fai_path}")
    else:
        log.debug("FAI file not found - coverage proportions will not be calculated")
    coverage_data: Dict[str, Dict[str, float]] = {}
    all_regions: List[pd.DataFrame] = []

    # 1. CNV regions
    cnv_bed = _get_latest_bed_file(bed_dir, "new_file_*.bed")
    if cnv_bed:
        log.debug(f"Loading CNV regions from: {cnv_bed}")
        cnv_df = _load_bed_file(cnv_bed)
        if not cnv_df.empty:
            cnv_expanded_df = _expand_cnv_regions(cnv_df, expand_size=10000)
            all_regions.append(cnv_expanded_df)
            log.debug(f"Loaded {len(cnv_df)} CNV regions (converted to {len(cnv_expanded_df)} breakpoint regions)")
            cnv_df_merged = _merge_overlapping_regions(cnv_df)
            cnv_entry = {"region_count": len(cnv_df)}
            if fai_path:
                cnv_coverage = _calculate_genome_coverage(cnv_df_merged, fai_path)
                if cnv_coverage:
                    cnv_entry.update({
                        "total_proportion": cnv_coverage.get("total_proportion", 0.0),
                        "plus_strand_proportion": cnv_coverage.get("plus_strand_proportion", 0.0),
                        "minus_strand_proportion": cnv_coverage.get("minus_strand_proportion", 0.0),
                    })
            coverage_data["cnv_regions"] = cnv_entry

    # 2. CNV breakpoints
    breakpoints_bed = _get_latest_bed_file(bed_dir, "breakpoints_*.bed")
    if breakpoints_bed:
        log.debug(f"Loading CNV breakpoints from: {breakpoints_bed}")
        bp_df = _load_bed_file(breakpoints_bed)
        if not bp_df.empty:
            bp_df_processed = _duplicate_unstranded_regions(bp_df)
            all_regions.append(bp_df_processed)
            log.debug(f"Loaded {len(bp_df_processed)} CNV breakpoints (with both strands)")
            bp_df_merged = _merge_overlapping_regions(bp_df_processed)
            bp_entry = {"region_count": len(bp_df_processed)}
            if fai_path:
                bp_coverage = _calculate_genome_coverage(bp_df_merged, fai_path)
                if bp_coverage:
                    bp_entry.update({
                        "total_proportion": bp_coverage.get("total_proportion", 0.0),
                        "plus_strand_proportion": bp_coverage.get("plus_strand_proportion", 0.0),
                        "minus_strand_proportion": bp_coverage.get("minus_strand_proportion", 0.0),
                    })
            coverage_data["cnv_breakpoints"] = bp_entry

    # 3. Fusion breakpoints
    fusion_bed = _get_latest_bed_file(bed_dir, "fusion_breakpoints_*.bed")
    if fusion_bed:
        log.debug(f"Loading fusion breakpoints from: {fusion_bed}")
        fusion_df = _load_bed_file(fusion_bed)
        if not fusion_df.empty:
            fusion_df_processed = _duplicate_unstranded_regions(fusion_df)
            all_regions.append(fusion_df_processed)
            log.debug(f"Loaded {len(fusion_df_processed)} fusion breakpoints (preserving strand)")
            fusion_df_merged = _merge_overlapping_regions(fusion_df_processed)
            fusion_entry = {"region_count": len(fusion_df_processed)}
            if fai_path:
                fusion_coverage = _calculate_genome_coverage(fusion_df_merged, fai_path)
                if fusion_coverage:
                    fusion_entry.update({
                        "total_proportion": fusion_coverage.get("total_proportion", 0.0),
                        "plus_strand_proportion": fusion_coverage.get("plus_strand_proportion", 0.0),
                        "minus_strand_proportion": fusion_coverage.get("minus_strand_proportion", 0.0),
                    })
            coverage_data["fusion_breakpoints"] = fusion_entry

    # 4. Master BED breakpoints
    master_bed_bp = _get_latest_bed_file(bed_dir, "master_bed_breakpoints_*.bed")
    if master_bed_bp:
        log.debug(f"Loading master BED breakpoints from: {master_bed_bp}")
        master_bed_bp_df = _load_bed_file(master_bed_bp)
        if not master_bed_bp_df.empty:
            master_bed_bp_df_processed = _duplicate_unstranded_regions(master_bed_bp_df)
            all_regions.append(master_bed_bp_df_processed)
            log.debug(f"Loaded {len(master_bed_bp_df_processed)} master BED breakpoints (with both strands)")
            master_bed_bp_df_merged = _merge_overlapping_regions(master_bed_bp_df_processed)
            master_bed_bp_entry = {"region_count": len(master_bed_bp_df_processed)}
            if fai_path:
                master_bed_bp_coverage = _calculate_genome_coverage(master_bed_bp_df_merged, fai_path)
                if master_bed_bp_coverage:
                    master_bed_bp_entry.update({
                        "total_proportion": master_bed_bp_coverage.get("total_proportion", 0.0),
                        "plus_strand_proportion": master_bed_bp_coverage.get("plus_strand_proportion", 0.0),
                        "minus_strand_proportion": master_bed_bp_coverage.get("minus_strand_proportion", 0.0),
                    })
            coverage_data["master_bed_breakpoints"] = master_bed_bp_entry

    if not all_regions:
        log.debug("No BED regions found to merge")
        return None

    combined_df = pd.concat(all_regions, ignore_index=True)
    log.debug(f"Combined {len(combined_df)} total regions")
    merged_df = _merge_overlapping_regions(combined_df)
    log.debug(f"Merged to {len(merged_df)} regions")

    if target_panel:
        target_bed_path = _get_target_bed_path(target_panel)
        if target_bed_path:
            log.debug(f"Loading target panel regions: {target_panel}")
            target_df = _load_bed_file(target_bed_path, require_bed6=True)
            if not target_df.empty:
                target_df_processed = _duplicate_unstranded_regions(target_df)
                target_df_merged = _merge_overlapping_regions(target_df_processed)
                target_entry = {"region_count": len(target_df_processed)}
                if fai_path:
                    target_coverage = _calculate_genome_coverage(target_df_merged, fai_path)
                    if target_coverage:
                        target_entry.update({
                            "total_proportion": target_coverage.get("total_proportion", 0.0),
                            "plus_strand_proportion": target_coverage.get("plus_strand_proportion", 0.0),
                            "minus_strand_proportion": target_coverage.get("minus_strand_proportion", 0.0),
                        })
                coverage_data["target_panel"] = target_entry
                log.debug(f"Keeping all {len(merged_df)} breakpoint regions and adding {len(target_df_processed)} target gene regions")
                all_regions_with_target = pd.concat([merged_df, target_df_processed], ignore_index=True)
                merged_df = _merge_overlapping_regions(all_regions_with_target)
                log.debug(f"Final merged regions including target genes: {len(merged_df)}")
            else:
                log.warning(f"Target panel BED file is empty: {target_bed_path}")
        else:
            log.warning(f"Target panel BED file not found for {target_panel}")
    else:
        log.debug("No target panel specified - master BED will include all breakpoint regions")

    if merged_df.empty:
        log.debug("No regions remaining after processing")
        return None

    sorted_df = _sort_bed_regions(merged_df)
    bed6_cols = ["chrom", "start", "end", "name", "score", "strand"]
    for col in bed6_cols:
        if col not in sorted_df.columns:
            if col == "score":
                sorted_df[col] = "0"
            elif col == "strand":
                sorted_df[col] = "."
            else:
                sorted_df[col] = "."
    master_entry = {"region_count": len(sorted_df)}
    if fai_path:
        master_coverage = _calculate_genome_coverage(sorted_df, fai_path)
        if master_coverage:
            master_entry.update({
                "total_proportion": master_coverage.get("total_proportion", 0.0),
                "plus_strand_proportion": master_coverage.get("plus_strand_proportion", 0.0),
                "minus_strand_proportion": master_coverage.get("minus_strand_proportion", 0.0),
            })
    coverage_data["master_bed"] = master_entry
    return (sorted_df, coverage_data)


def generate_master_bed(
    sample_id: str,
    work_dir: str,
    analysis_counter: int,
    target_panel: Optional[str] = None,
    logger_instance: Optional[logging.Logger] = None,
    reference: Optional[str] = None,
) -> Optional[str]:
    """
    Generate master BED file from all available BED file types.
    Lock is held only for the initial existence check and for the final write + log.
    Heavy work (load/merge/coverage) is done outside the lock to avoid contention.
    """
    if logger_instance:
        log = logger_instance
    else:
        log = logger

    sample_dir = os.path.join(work_dir, sample_id)
    bed_dir = os.path.join(sample_dir, "bed_files")
    if not os.path.exists(bed_dir):
        log.debug(f"BED directory does not exist: {bed_dir}")
        return None

    master_bed_path = os.path.join(bed_dir, f"master_{analysis_counter:03d}.bed")
    lock_file = os.path.join(sample_dir, "_locks", "master_bed.lock")
    lock_timeout = 300.0

    try:
        # Critical section 1: quick existence check only
        with FileLock(lock_file, timeout=lock_timeout):
            if os.path.exists(master_bed_path):
                log.debug(f"Master BED file already exists: {master_bed_path}")
                return master_bed_path
        # Lock released here

        log.debug(f"Generating master BED file for sample {sample_id} (counter: {analysis_counter})")
        build_result = _build_master_bed_data(
            sample_id, work_dir, analysis_counter, target_panel, reference, log
        )
        if not build_result:
            return None
        sorted_df, coverage_data = build_result
        bed6_cols = ["chrom", "start", "end", "name", "score", "strand"]

        # Critical section 2: write BED and log coverage only
        with FileLock(lock_file, timeout=lock_timeout):
            if os.path.exists(master_bed_path):
                log.debug(f"Master BED was generated by another process: {master_bed_path}")
                return master_bed_path
            sorted_df[bed6_cols].to_csv(
                master_bed_path,
                sep="\t",
                header=False,
                index=False,
            )
            if coverage_data:
                _log_bed_coverage_data(
                    sample_id=sample_id,
                    work_dir=work_dir,
                    analysis_counter=analysis_counter,
                    coverage_data=coverage_data,
                )
            log.info(f"Generated master BED file: {master_bed_path} with {len(sorted_df)} regions")
            return master_bed_path

    except TimeoutError as e:
        if os.path.exists(master_bed_path):
            log.debug(
                f"Master BED was generated by another process while waiting for lock: {master_bed_path}"
            )
            return master_bed_path
        log.error(f"Could not acquire lock for master BED generation: {e}")
        return None
    except Exception as e:
        log.error(f"Error generating master BED file: {e}")
        import traceback
        log.debug(f"Traceback: {traceback.format_exc()}")
        return None


def generate_master_bed_async(
    sample_id: str,
    work_dir: str,
    analysis_counter: int,
    target_panel: Optional[str] = None,
    logger_instance: Optional[logging.Logger] = None,
    reference: Optional[str] = None,
) -> None:
    """
    Generate master BED file asynchronously in a background thread (non-blocking).
    
    This function checks if the master BED file already exists, and if not,
    spawns a background thread to generate it. Returns immediately without
    waiting for generation to complete.
    
    Args:
        sample_id: Sample ID
        work_dir: Working directory
        analysis_counter: Current analysis counter
        target_panel: Target panel name (for gene intersection)
        logger_instance: Optional logger instance
        reference: Optional reference genome path
        
    Returns:
        None (returns immediately, generation happens in background)
    """
    # Check if file already exists (idempotency)
    sample_dir = os.path.join(work_dir, sample_id)
    bed_dir = os.path.join(sample_dir, "bed_files")
    master_bed_path = os.path.join(bed_dir, f"master_{analysis_counter:03d}.bed")
    
    if os.path.exists(master_bed_path):
        if logger_instance:
            logger_instance.debug(f"Master BED file already exists: {master_bed_path} - skipping async generation")
        return
    
    # Spawn background thread for generation
    def _generate_in_background():
        try:
            generate_master_bed(
                sample_id=sample_id,
                work_dir=work_dir,
                analysis_counter=analysis_counter,
                target_panel=target_panel,
                logger_instance=logger_instance,
                reference=reference,
            )
        except Exception as e:
            log = logger_instance if logger_instance else logger
            log.error(f"Error in background master BED generation for {sample_id}: {e}")
            import traceback
            log.debug(f"Traceback: {traceback.format_exc()}")
    
    thread = threading.Thread(target=_generate_in_background, daemon=True)
    thread.start()
    
    if logger_instance:
        logger_instance.debug(f"Spawned background thread for master BED generation: {sample_id} (counter: {analysis_counter})")


def generate_master_bed_from_files(
    cnv_regions_file: Optional[str] = None,
    cnv_breakpoints_file: Optional[str] = None,
    fusion_breakpoints_file: Optional[str] = None,
    master_bed_breakpoints_file: Optional[str] = None,
    target_bed_file: Optional[str] = None,
    output_file: str = "master.bed",
    expand_cnv_size: int = 10000,
    logger_instance: Optional[logging.Logger] = None,
) -> Optional[str]:
    """
    Generate master BED file from explicit file paths (for CLI/testing).
    
    Args:
        cnv_regions_file: Path to CNV regions BED file (will be expanded by +/- expand_cnv_size)
        cnv_breakpoints_file: Path to CNV breakpoints BED file
        fusion_breakpoints_file: Path to fusion breakpoints BED file
        target_bed_file: Path to target gene panel BED file (optional)
        output_file: Path to output master BED file
        expand_cnv_size: Size to expand CNV regions in base pairs (default 10kb)
        logger_instance: Optional logger instance
        
    Returns:
        Path to generated master BED file, or None if generation failed
    """
    if logger_instance:
        log = logger_instance
    else:
        log = logger
    
    try:
        all_regions = []
        
        # 1. Load and convert CNV regions to breakpoints
        if cnv_regions_file and os.path.exists(cnv_regions_file):
            log.debug(f"Loading CNV regions from: {cnv_regions_file}")
            cnv_df = _load_bed_file(cnv_regions_file)
            if not cnv_df.empty:
                expanded_df = _expand_cnv_regions(cnv_df, expand_size=expand_cnv_size)
                all_regions.append(expanded_df)
                log.debug(f"Loaded {len(cnv_df)} CNV regions (converted to {len(expanded_df)} breakpoint regions)")
        
        # 2. Load CNV breakpoints
        if cnv_breakpoints_file and os.path.exists(cnv_breakpoints_file):
            log.debug(f"Loading CNV breakpoints from: {cnv_breakpoints_file}")
            bp_df = _load_bed_file(cnv_breakpoints_file)
            if not bp_df.empty:
                # Ensure CNV breakpoints have both strands
                bp_df = _duplicate_unstranded_regions(bp_df)
                all_regions.append(bp_df)
                log.debug(f"Loaded {len(bp_df)} CNV breakpoints (with both strands)")
        
        # 3. Load fusion breakpoints
        if fusion_breakpoints_file and os.path.exists(fusion_breakpoints_file):
            log.debug(f"Loading fusion breakpoints from: {fusion_breakpoints_file}")
            fusion_df = _load_bed_file(fusion_breakpoints_file)
            if not fusion_df.empty:
                # Fusion breakpoints preserve their strand
                fusion_df = _duplicate_unstranded_regions(fusion_df)
                all_regions.append(fusion_df)
                log.debug(f"Loaded {len(fusion_df)} fusion breakpoints (preserving strand)")
        
        # 4. Load master BED breakpoints (if provided as argument)
        # Note: In the main generate_master_bed function, these are loaded automatically
        # This parameter is for the CLI/testing function
        if master_bed_breakpoints_file and os.path.exists(master_bed_breakpoints_file):
            log.debug(f"Loading master BED breakpoints from: {master_bed_breakpoints_file}")
            master_bed_bp_df = _load_bed_file(master_bed_breakpoints_file)
            if not master_bed_bp_df.empty:
                # Master BED breakpoints should have both strands
                master_bed_bp_df = _duplicate_unstranded_regions(master_bed_bp_df)
                all_regions.append(master_bed_bp_df)
                log.debug(f"Loaded {len(master_bed_bp_df)} master BED breakpoints (with both strands)")
        
        if not all_regions:
            log.error("No BED regions found to merge")
            return None
        
        # Combine all regions
        combined_df = pd.concat(all_regions, ignore_index=True)
        log.debug(f"Combined {len(combined_df)} total regions")
        
        # Merge overlapping regions
        merged_df = _merge_overlapping_regions(combined_df)
        log.debug(f"Merged to {len(merged_df)} regions")
        
        # Add target gene panel regions if available
        if target_bed_file and os.path.exists(target_bed_file):
            log.debug(f"Loading target panel regions from: {target_bed_file}")
            target_df = _load_bed_file(target_bed_file, require_bed6=True)
            if not target_df.empty:
                # Duplicate unstranded target regions for both strands
                target_df = _duplicate_unstranded_regions(target_df)
                
                # Keep all breakpoint regions and add target gene regions
                log.debug(f"Keeping all {len(merged_df)} breakpoint regions and adding {len(target_df)} target gene regions")
                
                # Combine all breakpoints with target gene regions
                all_regions_with_target = pd.concat([merged_df, target_df], ignore_index=True)
                
                log.debug(f"Combined {len(merged_df)} breakpoints with {len(target_df)} target gene regions")
                
                # Merge overlapping regions
                merged_df = _merge_overlapping_regions(all_regions_with_target)
                log.debug(f"Final merged regions including target genes: {len(merged_df)}")
        
        if merged_df.empty:
            log.error("No regions remaining after processing")
            return None
        
        # Sort regions
        sorted_df = _sort_bed_regions(merged_df)
        
        # Ensure we have all BED6 columns
        bed6_cols = ["chrom", "start", "end", "name", "score", "strand"]
        for col in bed6_cols:
            if col not in sorted_df.columns:
                if col == "score":
                    sorted_df[col] = "0"
                elif col == "strand":
                    sorted_df[col] = "."
                else:
                    sorted_df[col] = "."
        
        # Write master BED file in BED6 format
        sorted_df[bed6_cols].to_csv(
            output_file,
            sep="\t",
            header=False,
            index=False,
        )
        
        log.info(f"Generated master BED file: {output_file} with {len(sorted_df)} regions")
        return output_file
        
    except Exception as e:
        log.error(f"Error generating master BED file: {e}")
        import traceback
        log.debug(f"Traceback: {traceback.format_exc()}")
        return None


def main():
    """CLI entry point for standalone testing."""
    parser = argparse.ArgumentParser(
        description="Generate master BED file from CNV and fusion breakpoint BED files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Combine CNV regions and breakpoints:
  python master_bed_generator.py --cnv-regions new_file_001.bed \\
                                  --cnv-breakpoints breakpoints_001.bed \\
                                  --output master.bed

  # Include fusion breakpoints and target panel:
  python master_bed_generator.py --cnv-regions new_file_001.bed \\
                                  --cnv-breakpoints breakpoints_001.bed \\
                                  --fusion-breakpoints fusion_breakpoints_001.bed \\
                                  --target-panel rCNS2_panel_name_uniq.bed \\
                                  --output master.bed

  # Include master BED breakpoints (new target regions from supplementary alignments):
  python master_bed_generator.py --cnv-regions new_file_001.bed \\
                                  --cnv-breakpoints breakpoints_001.bed \\
                                  --fusion-breakpoints fusion_breakpoints_001.bed \\
                                  --master-bed-breakpoints master_bed_breakpoints_001.bed \\
                                  --output master.bed

  # Custom CNV expansion size:
  python master_bed_generator.py --cnv-regions new_file_001.bed \\
                                  --expand-size 5000 \\
                                  --output master.bed
        """
    )
    
    parser.add_argument(
        "--cnv-regions",
        type=str,
        help=(
            "Path to CNV regions BED file (new_file_*.bed). "
            "Processing: Each CNV region is split into two breakpoints (start and end positions). "
            "Each breakpoint is padded by +/- expand-size (default 10kb) around the breakpoint position. "
            "Each breakpoint region is duplicated for both + and - strands. "
            "Regions are then merged if they overlap (separately by strand)."
        )
    )
    parser.add_argument(
        "--cnv-breakpoints",
        type=str,
        help=(
            "Path to CNV breakpoints BED file (breakpoints_*.bed). "
            "Processing: Regions are loaded as-is. If regions have unstranded ('.') strand, "
            "they are duplicated for both + and - strands. Regions are merged if they overlap (separately by strand)."
        )
    )
    parser.add_argument(
        "--fusion-breakpoints",
        type=str,
        help=(
            "Path to fusion breakpoints BED file (fusion_breakpoints_*.bed). "
            "Processing: Regions preserve their original strand information. "
            "If regions have unstranded ('.') strand, they are duplicated for both + and - strands. "
            "Regions are merged if they overlap (separately by strand)."
        )
    )
    parser.add_argument(
        "--master-bed-breakpoints",
        type=str,
        help=(
            "Path to master BED breakpoints BED file (master_bed_breakpoints_*.bed). "
            "These are new target regions identified from supplementary alignments of reads mapping to master BED regions. "
            "Processing: Regions are loaded as-is. If regions have unstranded ('.') strand, "
            "they are duplicated for both + and - strands (breaks can occur on either strand). "
            "Regions are merged if they overlap (separately by strand)."
        )
    )
    parser.add_argument(
        "--target-panel",
        type=str,
        help=(
            "Path to target gene panel BED file (optional, for intersection). "
            "Processing: Target regions are loaded and unstranded regions are duplicated for both + and - strands. "
            "All breakpoint regions are kept (not filtered). "
            "All target gene regions are then added to the final output. "
            "Final regions are merged if they overlap (separately by strand)."
        )
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to output master BED file"
    )
    parser.add_argument(
        "--expand-size",
        type=int,
        default=10000,
        help="Size to expand CNV regions in base pairs (default: 10000 = 10kb)"
    )
    parser.add_argument(
        "--fai-file",
        type=str,
        help=(
            "Path to FASTA index (.fai) file (optional). "
            "If provided, calculates the total proportion of the genome covered by the BED file, "
            "accounting for strand information (+ and - strands are calculated separately)."
        )
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    
    # Check that at least one input file is provided
    if not any([args.cnv_regions, args.cnv_breakpoints, args.fusion_breakpoints, args.master_bed_breakpoints]):
        parser.error(
            "At least one input BED file must be provided "
            "(--cnv-regions, --cnv-breakpoints, --fusion-breakpoints, or --master-bed-breakpoints)"
        )
    
    # Generate master BED file
    result = generate_master_bed_from_files(
        cnv_regions_file=args.cnv_regions,
        cnv_breakpoints_file=args.cnv_breakpoints,
        fusion_breakpoints_file=args.fusion_breakpoints,
        master_bed_breakpoints_file=args.master_bed_breakpoints,
        target_bed_file=args.target_panel,
        output_file=args.output,
        expand_cnv_size=args.expand_size,
    )
    
    if result:
        print(f"Successfully generated master BED file: {result}")
        
        # Calculate genome coverage if FAI file is provided
        if args.fai_file:
            try:
                # Load the generated BED file
                bed_df = _load_bed_file(result)
                if not bed_df.empty:
                    coverage_stats = _calculate_genome_coverage(bed_df, args.fai_file)
                    if coverage_stats:
                        print("\n" + "="*60)
                        print("Genome Coverage Statistics")
                        print("="*60)
                        print(f"Total genome size: {coverage_stats['total_genome_size']:,} bp")
                        print(f"\nCoverage by strand:")
                        print(f"  + strand: {coverage_stats['plus_strand_coverage']:,} bp "
                              f"({coverage_stats['plus_strand_proportion']*100:.4f}%)")
                        print(f"  - strand: {coverage_stats['minus_strand_coverage']:,} bp "
                              f"({coverage_stats['minus_strand_proportion']*100:.4f}%)")
                        if coverage_stats['unstranded_coverage'] > 0:
                            print(f"  . strand: {coverage_stats['unstranded_coverage']:,} bp "
                                  f"({coverage_stats['unstranded_proportion']*100:.4f}%)")
                        print(f"\nTotal covered bases: {coverage_stats['total_covered_bases']:,} bp")
                        print(f"Total proportion covered: {coverage_stats['total_proportion']*100:.4f}%")
                        print("="*60)
                    else:
                        print("Warning: Could not calculate genome coverage statistics", file=sys.stderr)
                else:
                    print("Warning: Generated BED file is empty, cannot calculate coverage", file=sys.stderr)
            except Exception as e:
                logger.error(f"Error calculating genome coverage: {e}")
                print(f"Warning: Error calculating genome coverage: {e}", file=sys.stderr)
        
        sys.exit(0)
    else:
        print("Failed to generate master BED file", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
