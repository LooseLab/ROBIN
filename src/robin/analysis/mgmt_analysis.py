#!/usr/bin/env python3
"""
MGMT (O6-methylguanine-DNA methyltransferase) promoter methylation analysis.

This module provides automated processing of BAM files to analyze methylation
patterns in the MGMT promoter region (chr10:129466536-129467536). The analysis
integrates with robin's preprocessing pipeline and provides rich metadata
output.

Features:
- Automated BAM processing with robin's preprocessing pipeline
- MGMT region extraction using bedtools
- Methylation analysis using matkit from robin package
- R-based prediction model for methylation status
- Methylation visualization using methylartist
- Comprehensive metadata extraction and logging
"""

import logging
import os
import shutil
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import pysam

from robin.logging_config import get_job_logger

# Try to find HV path from robin package, fallback to common locations
try:
    HVPATH = os.path.join(
        os.path.dirname(os.path.abspath(__import__("robin").__file__)),
        "submodules",
        "hv_rapidCNS2",
    )
except ImportError:
    # Fallback if robin package is not importable
    HVPATH = None


def has_reads(bam_file, chrom, start, end):
    """
    Quickly checks if any reads span a specific genomic locus.

    Args:
        bam_file (str): Path to your BAM file.
        chrom (str): Chromosome name (e.g., 'chr10').
        start (int): Start coordinate (0-based, inclusive).
        end (int): End coordinate (0-based, exclusive).

    Returns:
        bool: True if at least one read spans the region, False otherwise.
    """
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for _ in bam.fetch(chrom, start, end):
            return True  # Found at least one read
    return False  # No reads found


def safely_sort_and_index_bam(
    input_bam: str,
    output_bam: str,
    logger: logging.Logger,
    threads: int = 4,
    verify_readable: bool = True,
) -> bool:
    """
    Safely sort and index a BAM file with proper verification to prevent truncated file issues.

    This function ensures that:
    1. The input BAM file exists and is readable
    2. The sorted BAM file is fully written before indexing
    3. The indexed BAM file is verified to be readable
    4. Proper error handling and logging throughout

    Args:
        input_bam: Path to input BAM file
        output_bam: Path to output sorted BAM file
        logger: Logger instance for logging
        threads: Number of threads to use for sorting
        verify_readable: If True, verify the output BAM is readable after creation

    Returns:
        True if successful, False otherwise
    """
    try:
        # Verify input file exists
        if not os.path.exists(input_bam):
            logger.error(f"Input BAM file does not exist: {input_bam}")
            return False

        input_size = os.path.getsize(input_bam)
        if input_size == 0:
            logger.warning(f"Input BAM file is empty: {input_bam}")
            return False

        logger.debug(
            f"Sorting BAM file: {input_bam} -> {output_bam} (input size: {input_size} bytes)"
        )

        # Remove output file if it exists (to avoid issues with partial writes)
        if os.path.exists(output_bam):
            try:
                os.remove(output_bam)
                if os.path.exists(f"{output_bam}.bai"):
                    os.remove(f"{output_bam}.bai")
            except OSError as e:
                logger.warning(
                    f"Could not remove existing output file {output_bam}: {e}"
                )

        # Sort the BAM file
        try:
            if threads > 1:
                pysam.sort(f"-@{threads}", "-o", output_bam, input_bam)
            else:
                pysam.sort("-o", output_bam, input_bam)
        except Exception as e:
            logger.error(f"Failed to sort BAM file {input_bam}: {e}")
            return False

        # Verify sorted file was created and has content
        if not os.path.exists(output_bam):
            logger.error(f"Sorted BAM file was not created: {output_bam}")
            return False

        output_size = os.path.getsize(output_bam)
        if output_size == 0:
            logger.error(f"Sorted BAM file is empty: {output_bam}")
            return False

        logger.debug(
            f"Sorted BAM file created: {output_bam} (size: {output_size} bytes)"
        )

        # Verify the sorted BAM is readable before indexing
        if verify_readable:
            try:
                with pysam.AlignmentFile(output_bam, "rb") as test_bam:
                    # Try to read the header and first few reads to verify integrity
                    _ = test_bam.header
                    # Count reads to verify file is not truncated
                    read_count = 0
                    for _ in test_bam.fetch(until_eof=True):
                        read_count += 1
                        if read_count >= 10:  # Just verify first 10 reads
                            break
                    logger.debug(
                        f"Verified sorted BAM is readable (checked {read_count} reads)"
                    )
            except Exception as e:
                logger.error(f"Sorted BAM file appears corrupted or truncated: {e}")
                # Clean up corrupted file
                try:
                    os.remove(output_bam)
                except OSError:
                    pass
                return False

        # Index the sorted BAM file
        try:
            index_file = f"{output_bam}.bai"
            pysam.index(output_bam, index_file)

            # Verify index was created
            if not os.path.exists(index_file):
                logger.error(f"BAM index file was not created: {index_file}")
                return False

            logger.debug(f"BAM index created: {index_file}")
        except Exception as e:
            logger.error(f"Failed to index BAM file {output_bam}: {e}")
            # Clean up files if indexing failed
            try:
                if os.path.exists(output_bam):
                    os.remove(output_bam)
            except OSError:
                pass
            return False

        # Final verification: ensure both BAM and index are readable
        if verify_readable:
            try:
                with pysam.AlignmentFile(output_bam, "rb") as final_bam:
                    # Verify we can access the index
                    _ = final_bam.header
                    # Try a simple fetch to ensure index works
                    try:
                        # Get first chromosome from header
                        if final_bam.references:
                            first_chrom = final_bam.references[0]
                            list(
                                final_bam.fetch(first_chrom, 0, 1000)
                            )  # Try fetching a small region
                    except Exception:
                        # If fetch fails, that's okay - file might be empty or have no reads
                        pass
                logger.info(f"Successfully sorted and indexed BAM: {output_bam}")
            except Exception as e:
                logger.error(
                    f"Final verification failed for sorted BAM {output_bam}: {e}"
                )
                return False

        return True

    except Exception as e:
        logger.error(f"Unexpected error in safely_sort_and_index_bam: {e}")
        import traceback

        logger.error(traceback.format_exc())
        return False


def validate_methylation_data(bam_file: str, logger: logging.Logger) -> Dict[str, Any]:
    """
    Validate that BAM file has sufficient methylation data for methylartist.

    Args:
        bam_file (str): Path to BAM file
        logger (logging.Logger): Logger instance

    Returns:
        Dict containing validation results and metadata
    """
    validation_result = {
        "has_sufficient_data": False,
        "total_reads": 0,
        "reads_with_mm_tags": 0,
        "reads_with_ml_tags": 0,
        "modification_types": set(),
        "error_message": None,
    }

    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            read_count = 0
            mm_count = 0
            ml_count = 0

            for read in bam.fetch(until_eof=True):
                read_count += 1

                # Check for MM tags (modification calls)
                if read.has_tag("MM"):
                    mm_count += 1
                    mm_tag = read.get_tag("MM")
                    if mm_tag:
                        # Parse modification types from MM tag
                        mods = mm_tag.split(";")
                        for mod in mods:
                            if ":" in mod:
                                mod_type = mod.split(":")[0]
                                validation_result["modification_types"].add(mod_type)

                # Check for ML tags (modification likelihoods)
                if read.has_tag("ML"):
                    ml_count += 1

                # Stop counting after reasonable sample size to avoid long processing
                if read_count >= 10000:
                    break

            validation_result["total_reads"] = read_count
            validation_result["reads_with_mm_tags"] = mm_count
            validation_result["reads_with_ml_tags"] = ml_count

            # Determine if we have sufficient data
            # Need at least some reads with methylation data
            min_reads_with_mods = 10
            has_modifications = len(validation_result["modification_types"]) > 0

            if mm_count >= min_reads_with_mods and has_modifications:
                validation_result["has_sufficient_data"] = True
                logger.info(
                    f"Methylation data validation passed: {mm_count} reads with MM tags, mods: {validation_result['modification_types']}"
                )
            else:
                validation_result["error_message"] = (
                    f"Insufficient methylation data: {mm_count} reads with MM tags (min: {min_reads_with_mods}), modifications: {validation_result['modification_types']}"
                )
                logger.warning(validation_result["error_message"])

    except Exception as e:
        validation_result["error_message"] = (
            f"Error validating methylation data: {str(e)}"
        )
        logger.error(validation_result["error_message"])

    return validation_result


def extract_mgmt_site_rows_from_bed(bed_path: str) -> List[Dict[str, Any]]:
    """
    Extract MGMT-specific CpG site information from a bedmethyl file.

    Args:
        bed_path (str): Path to bedmethyl file

    Returns:
        List of dictionaries containing site information for annotation
    """
    try:
        import pandas as pd

        if not os.path.exists(bed_path):
            return []

        df = pd.read_csv(bed_path, sep="\t", header=None)

        # Check if column 10 contains space-separated values (old format)
        has_space_separated_col10 = False
        if df.shape[1] > 10 and len(df) > 0:
            sample_val = str(df.iloc[0, 9])
            has_space_separated_col10 = " " in sample_val or "\t" in sample_val

        # Check if this is the new bedmethyl format (separate columns) or old format
        if df.shape[1] >= 12 and not has_space_separated_col10:
            # New bedmethyl format with separate columns
            cols = [
                "Chromosome",
                "Start",
                "End",
                "Modified_Base_Code",
                "Score",
                "Strand",
                "Start2",
                "End2",
                "RGB",
                "Nvalid_cov",
                "Fraction_Modified",
                "Nmod",
            ]
            num_cols_to_read = min(len(cols), df.shape[1])
            df = df.iloc[:, :num_cols_to_read]
            df.columns = cols[:num_cols_to_read]

            df["Nvalid_cov"] = df["Nvalid_cov"].astype(float)
            df["Fraction_Modified"] = df["Fraction_Modified"].astype(float)
            # Some bedmethyl outputs store fraction as percent (0-100).
            # Normalize to 0-1 if needed.
            if len(df) > 0 and (df["Fraction_Modified"] > 1.0).any():
                df["Fraction_Modified"] = df["Fraction_Modified"] / 100.0
            if "Nmod" in df.columns:
                df["Nmod"] = df["Nmod"].astype(float)
            else:
                df["Nmod"] = df["Nvalid_cov"] * df["Fraction_Modified"]

            df["Start"] = df["Start"].astype(int)
            df["Coverage"] = df["Nvalid_cov"]
            df["Modified_Fraction"] = df["Fraction_Modified"] * 100.0
        elif df.shape[1] >= 10:
            # Old format with space-separated Coverage_Info in column 10
            cols = [
                "Chromosome",
                "Start",
                "End",
                "Name",
                "Score",
                "Strand",
                "Start2",
                "End2",
                "RGB",
                "Coverage_Info",
            ]
            df = df.iloc[:, : len(cols)]
            df.columns = cols

            cov_split = df["Coverage_Info"].astype(str).str.split()
            df["Coverage"] = cov_split.str[0].astype(float)
            fraction_val = cov_split.str[1].astype(float).fillna(0.0)
            is_percentage = (
                (fraction_val > 1.0).any() if len(fraction_val) > 0 else False
            )

            if is_percentage:
                df["Modified_Fraction"] = fraction_val
                df["Fraction_Modified"] = df["Modified_Fraction"] / 100.0
            else:
                df["Fraction_Modified"] = fraction_val
                df["Modified_Fraction"] = df["Fraction_Modified"] * 100.0

            df["Nvalid_cov"] = df["Coverage"]
            df["Nmod"] = df["Coverage"] * df["Fraction_Modified"]
            df["Start"] = df["Start"].astype(int)
        else:
            return []

        cpg_pairs = [
            (129467255, 129467256),
            (129467258, 129467259),
            (129467262, 129467263),
            (129467272, 129467273),
        ]
        rows: List[Dict[str, Any]] = []
        label_map = {
            "129467255/129467256": "1",
            "129467258/129467259": "2",
            "129467262/129467263": "3",
            "129467272/129467273": "4",
        }

        for p1, p2 in cpg_pairs:
            pos_key = f"{p1}/{p2}"
            site_label = label_map.get(pos_key, "Unknown")

            # Check forward strand reads at position p1
            fwd_p1 = df[
                (df["Chromosome"] == "chr10")
                & (df["Start"] == p1 - 1)
                & (df["Strand"] == "+")
            ]

            # Check reverse strand reads at position p2
            rev_p2 = df[
                (df["Chromosome"] == "chr10")
                & (df["Start"] == p2 - 1)
                & (df["Strand"] == "-")
            ]

            # Get forward strand data from p1
            if not fwd_p1.empty:
                cov_f = float(fwd_p1["Nvalid_cov"].iloc[0])
                mf = float(fwd_p1["Fraction_Modified"].iloc[0])
                nmod_f = float(fwd_p1["Nmod"].iloc[0])
                meth_fwd_count = int(round(nmod_f))
            else:
                cov_f = 0.0
                mf = 0.0
                meth_fwd_count = 0

            # Get reverse strand data from p2
            if not rev_p2.empty:
                cov_r = float(rev_p2["Nvalid_cov"].iloc[0])
                mr = float(rev_p2["Fraction_Modified"].iloc[0])
                nmod_r = float(rev_p2["Nmod"].iloc[0])
                meth_rev_count = int(round(nmod_r))
            else:
                cov_r = 0.0
                mr = 0.0
                meth_rev_count = 0

            # Only add row if we have data
            if cov_f > 0 or cov_r > 0:
                tot = cov_f + cov_r
                weighted = ((cov_f * mf) + (cov_r * mr)) / tot if tot > 0 else 0.0
                weighted_pct = weighted * 100.0

                rows.append(
                    {
                        "site": f"{site_label} (CpG {pos_key})",
                        "chr": "chr10",
                        "pos": pos_key,
                        "cov_fwd": int(cov_f),
                        "cov_rev": int(cov_r),
                        "cov_total": int(tot),
                        "meth": round(weighted_pct, 2),
                        "meth_fwd": int(meth_fwd_count),
                        "meth_rev": int(meth_rev_count),
                        "notes": "Combined methylation from both strands of CpG pair",
                    }
                )

        return rows
    except Exception as e:
        logging.getLogger("robin.mgmt").debug(
            f"Failed to extract site rows from {bed_path}: {e}"
        )
        return []


def generate_mgmt_visualization(
    bam_file: str,
    output_png: str,
    interval: str = "chr10:129466536-129467536",
    reference: Optional[str] = None,
    bed_file: Optional[str] = None,
    logger: Optional[logging.Logger] = None,
    extra_cli: Optional[List[str]] = None,
    use_fallback: bool = False,
) -> Dict[str, Any]:
    """
    Generate MGMT methylation visualization with both PNG and pickle outputs.

    This unified function uses the locus_figure approach to capture the matplotlib
    figure object, then saves both a PNG file (for reports) and a pickle file
    (for fast GUI loading). The figure includes MGMT CpG site annotations.

    Args:
        bam_file (str): Path to BAM file
        output_png (str): Output PNG file path
        interval (str): Genomic interval for analysis
        reference (Optional[str]): Reference genome path
        bed_file (Optional[str]): Path to bedmethyl file for extracting site annotations
        logger (Optional[logging.Logger]): Logger instance
        extra_cli (Optional[List[str]]): Extra CLI parameters for methylartist
        use_fallback (bool): Whether to use fallback parameters (more lenient)

    Returns:
        Dict containing execution results, metadata, and paths to generated files
    """
    if logger is None:
        logger = logging.getLogger("robin.mgmt")

    result = {
        "success": False,
        "error_message": None,
        "validation_passed": False,
        "fallback_used": use_fallback,
        "png_path": output_png,
        "pickle_path": None,
        "figure": None,
    }

    # Determine pickle path (same directory as PNG, with .pkl extension)
    pickle_path = os.path.splitext(output_png)[0] + ".pkl"
    result["pickle_path"] = pickle_path

    try:
        # Step 1: Validate methylation data
        logger.info(f"Validating methylation data in {os.path.basename(bam_file)}")
        validation = validate_methylation_data(bam_file, logger)

        if not validation["has_sufficient_data"]:
            result["error_message"] = (
                f"Methylation data validation failed: {validation['error_message']}"
            )
            logger.warning(result["error_message"])
            return result

        result["validation_passed"] = True
        logger.info("Methylation data validation passed, proceeding with methylartist")

        # Step 2: Extract site_rows from bed file if provided
        site_rows = None
        if bed_file and os.path.exists(bed_file):
            try:
                site_rows = extract_mgmt_site_rows_from_bed(bed_file)
                if site_rows:
                    logger.debug(
                        f"Extracted {len(site_rows)} site rows from {bed_file}"
                    )
            except Exception as e:
                logger.debug(f"Failed to extract site rows from {bed_file}: {e}")

        # Step 3: Build extra_cli parameters
        if extra_cli is None:
            extra_cli = []

        # Add figure size parameters if not already specified
        has_width = any("--width" in str(arg) for arg in extra_cli)
        has_height = any("--height" in str(arg) for arg in extra_cli)
        if not has_width:
            extra_cli.extend(["--width", "24"])
        if not has_height:
            extra_cli.extend(["--height", "12"])

        # Add quality/read parameters based on fallback mode
        if use_fallback:
            # More lenient parameters for fallback
            if not any("--minreads" in str(arg) for arg in extra_cli):
                extra_cli.extend(["--minreads", "1"])
            if not any("--minqual" in str(arg) for arg in extra_cli):
                extra_cli.extend(["--minqual", "5"])
            if not any("--smoothwindowsize" in str(arg) for arg in extra_cli):
                extra_cli.extend(["--smoothwindowsize", "1"])
            if not any("--reads" in str(arg) for arg in extra_cli):
                extra_cli.extend(["--reads", "100"])
        else:
            # Standard parameters
            if not any("--minreads" in str(arg) for arg in extra_cli):
                extra_cli.extend(["--minreads", "5"])
            if not any("--minqual" in str(arg) for arg in extra_cli):
                extra_cli.extend(["--minqual", "10"])
            if not any("--smoothwindowsize" in str(arg) for arg in extra_cli):
                extra_cli.extend(["--smoothwindowsize", "5"])

        # Step 4: Generate figure using locus_figure
        try:
            from robin.analysis.methylation_wrapper import (
                locus_figure,
                save_figure_pickle,
            )

            logger.info(f"Generating methylation visualization using locus_figure")
            fig = locus_figure(
                interval=interval,
                bam_path=bam_file,
                motif="CG",
                mods="m",
                extra_cli=extra_cli,
                site_rows=site_rows,
            )

            result["figure"] = fig

            # Step 5: Save PNG file
            try:
                fig.savefig(output_png, dpi=150, bbox_inches="tight")
                logger.info(f"Saved PNG visualization: {os.path.basename(output_png)}")
            except Exception as e:
                logger.warning(f"Failed to save PNG file: {e}")
                # Continue anyway - pickle is more important for GUI

            # Step 6: Save pickle file for fast GUI loading
            try:
                save_figure_pickle(fig, pickle_path)
                logger.info(
                    f"Saved pickle file for fast GUI loading: {os.path.basename(pickle_path)}"
                )
            except Exception as e:
                logger.warning(f"Failed to save pickle file: {e}")
                # Not fatal, but less optimal

            result["success"] = True
            logger.info("Methylation visualization generated successfully")

        except RuntimeError as e:
            # Check if this is a validation error we can handle
            error_msg = str(e)
            if "not indexed" in error_msg or "index" in error_msg.lower():
                result["error_message"] = f"BAM file indexing issue: {error_msg}"
            elif "did not produce a figure" in error_msg:
                result["error_message"] = (
                    f"Methylartist failed to produce figure: {error_msg}"
                )
            else:
                result["error_message"] = f"Methylartist error: {error_msg}"
            logger.error(result["error_message"])

            # Try fallback if not already using it
            if not use_fallback:
                logger.info("Attempting fallback with more lenient parameters")
                return generate_mgmt_visualization(
                    bam_file=bam_file,
                    output_png=output_png,
                    interval=interval,
                    reference=reference,
                    bed_file=bed_file,
                    logger=logger,
                    extra_cli=extra_cli,  # Keep any custom extra_cli
                    use_fallback=True,
                )
        except Exception as e:
            result["error_message"] = (
                f"Unexpected error generating visualization: {str(e)}"
            )
            logger.error(result["error_message"], exc_info=True)

            # Try fallback if not already using it
            if not use_fallback:
                logger.info("Attempting fallback with more lenient parameters")
                return generate_mgmt_visualization(
                    bam_file=bam_file,
                    output_png=output_png,
                    interval=interval,
                    reference=reference,
                    bed_file=bed_file,
                    logger=logger,
                    extra_cli=extra_cli,
                    use_fallback=True,
                )

    except Exception as e:
        result["error_message"] = (
            f"Unexpected error in visualization generation: {str(e)}"
        )
        logger.error(result["error_message"], exc_info=True)

    return result


def run_methylartist_safely(
    bam_file: str,
    output_file: str,
    interval: str = "chr10:129466536-129467536",
    reference: Optional[str] = None,
    logger: Optional[logging.Logger] = None,
    bed_file: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Run methylartist with robust error handling and fallback options.

    This function now uses the unified generate_mgmt_visualization function
    which saves both PNG and pickle files for optimal performance.

    Args:
        bam_file (str): Path to BAM file
        output_file (str): Output PNG file path
        interval (str): Genomic interval for analysis
        reference (Optional[str]): Reference genome path
        logger (Optional[logging.Logger]): Logger instance
        bed_file (Optional[str]): Path to bedmethyl file for extracting site annotations

    Returns:
        Dict containing execution results and metadata (backward compatible format)
    """
    # Use the unified visualization function
    viz_result = generate_mgmt_visualization(
        bam_file=bam_file,
        output_png=output_file,
        interval=interval,
        reference=reference,
        bed_file=bed_file,
        logger=logger,
        use_fallback=False,
    )

    # Convert to backward-compatible format
    result = {
        "success": viz_result["success"],
        "error_message": viz_result["error_message"],
        "command_used": None,  # Not applicable with locus_figure approach
        "validation_passed": viz_result["validation_passed"],
        "fallback_used": viz_result["fallback_used"],
    }

    return result


@dataclass
class MGMTMetadata:
    """Container for MGMT analysis metadata and results"""

    sample_id: str
    bam_path: str
    analysis_timestamp: float
    processing_steps: List[str] = field(default_factory=list)
    error_message: Optional[str] = None
    results: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        if self.processing_steps is None:
            self.processing_steps = []
        if self.results is None:
            self.results = {}


def _find_mgmt_bed(work_dir: str) -> str:
    """Find or create the MGMT BED file"""
    possible_paths = [
        "mgmt_hg38.bed",
        "bin/mgmt_hg38.bed",
        "data/mgmt_hg38.bed",
        "/usr/local/share/mgmt_hg38.bed",
    ]

    for path in possible_paths:
        if os.path.exists(path):
            return path

    # Create default MGMT BED file
    default_path = os.path.join(work_dir, "mgmt_hg38.bed")
    with open(default_path, "w") as f:
        f.write("chr10\t129466536\t129467536\tMGMT_promoter\t1000\t+\n")

    return default_path


def _find_hv_path(work_dir: str) -> str:
    """Find the HV rapidCNS2 path from robin package"""
    # Use the already defined HVPATH constant if available
    if HVPATH and os.path.exists(HVPATH):
        return HVPATH

    # Fallback: create a basic structure if HVPATH doesn't exist
    hv_path = os.path.join(work_dir, "hv_rapidCNS2")
    os.makedirs(hv_path, exist_ok=True)
    os.makedirs(os.path.join(hv_path, "bin"), exist_ok=True)

    logger = logging.getLogger("robin.mgmt")
    if HVPATH:
        logger.warning(
            f"HVPATH not found at {HVPATH}, created basic structure at: {hv_path}"
        )
    else:
        logger.warning(
            f"robin package not available, created basic HV structure at: {hv_path}"
        )
    return hv_path


def _get_next_file_number(sample_dir: str) -> int:
    """Get the next file number for incremental naming based on existing files"""
    # Look for existing numbered files in the sample directory
    if not os.path.exists(sample_dir):
        return 1

    existing_numbers = set()
    for filename in os.listdir(sample_dir):
        # Look for files with pattern like "1_mgmt.bed", "2_mgmt.csv", etc.
        if filename.endswith((".bed", ".csv", ".png")):
            parts = filename.split("_")
            if len(parts) >= 2 and parts[0].isdigit():
                existing_numbers.add(int(parts[0]))

    # Return the next available number
    if existing_numbers:
        return max(existing_numbers) + 1
    else:
        return 1


def convert_to_mixed_delim(input_file: str, output_file: str) -> bool:
    """
    Converts standard BEDMethyl format to mixed-delim format for R script compatibility.

    The R script expects column V10 to contain a space-separated string like "10 0.5"
    instead of separate columns for coverage and fraction.
    """
    try:
        with open(input_file, "r") as infile, open(output_file, "w") as outfile:
            for line in infile:
                if line.strip():
                    fields = line.strip().split("\t")
                    if len(fields) >= 18:
                        # Standard BEDMethyl format has:
                        # V10: coverage (column 10, 0-indexed)
                        # V11: fraction modified (column 11, 0-indexed)

                        # Create mixed-delim format by combining V10 and V11 into V10
                        coverage = fields[9]  # V10 (0-indexed)
                        fraction = fields[10]  # V11 (0-indexed)

                        # Create the mixed-delim format: "coverage fraction"
                        mixed_delim = f"{coverage} {fraction}"

                        # Replace V10 with the mixed-delim string
                        fields[9] = mixed_delim

                        # Write the converted line
                        outfile.write("\t".join(fields) + "\n")
                    else:
                        # If the line doesn't have enough columns, write it as-is
                        outfile.write(line)

        logger = logging.getLogger("robin.mgmt")
        logger.info(f"Converted BEDMethyl format: {input_file} → {output_file}")
        return True

    except Exception as e:
        logger = logging.getLogger("robin.mgmt")
        logger.error(f"Error converting BEDMethyl format: {e}")
        return False


def run_matkit(
    temp_dir: str, mgmt_bamfile: str, hv_path: str, sample_dir: str, file_number: int
) -> bool:
    """Processes BAM file with matkit and runs R script for MGMT prediction."""
    temp_files = []  # Track temporary files for cleanup

    logger = logging.getLogger("robin.mgmt")
    try:
        logger.info("Running matkit for methylation analysis")

        # Sort and index the BAM file using safe function to prevent truncated file issues
        sorted_bam = os.path.join(sample_dir, f"mgmt_sorted_{file_number}.bam")
        temp_files.append(sorted_bam)
        temp_files.append(f"{sorted_bam}.bai")

        if not safely_sort_and_index_bam(
            input_bam=mgmt_bamfile,
            output_bam=sorted_bam,
            logger=logger,
            threads=4,
            verify_readable=True,
        ):
            logger.error(f"Failed to sort/index BAM file for matkit: {mgmt_bamfile}")
            return False

        # Use matkit from robin package
        try:
            from robin.analysis.utilities.matkit import run_matkit as run_matkit_util

            temp_bed_file = os.path.join(sample_dir, f"temp_mgmt_{file_number}.bed")
            temp_files.append(temp_bed_file)
            run_matkit_util(sorted_bam, temp_bed_file)

            # Convert to mixed delimiter format for R script
            output_mgmt_bed = os.path.join(sample_dir, f"mgmt_mixed_{file_number}.bed")
            temp_files.append(output_mgmt_bed)
            convert_to_mixed_delim(temp_bed_file, output_mgmt_bed)

            # Copy BED file to sample directory with incremental naming
            final_bed_file = os.path.join(sample_dir, f"{file_number}_mgmt.bed")
            shutil.copy2(output_mgmt_bed, final_bed_file)

            # Run R script for MGMT prediction
            r_script_path = os.path.join(hv_path, "bin", "mgmt_pred_v0.3.R")
            if os.path.exists(r_script_path):
                cmd = f"Rscript {r_script_path} --input={output_mgmt_bed} --out_dir={sample_dir} --probes={hv_path}/bin/mgmt_probes.Rdata --model={hv_path}/bin/mgmt_137sites_mean_model.Rdata --sample=live_analysis"
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

                if result.returncode == 0:
                    # Rename the output CSV file to use incremental naming
                    default_csv = os.path.join(
                        sample_dir, "live_analysis_mgmt_status.csv"
                    )
                    final_csv = os.path.join(sample_dir, f"{file_number}_mgmt.csv")
                    if os.path.exists(default_csv):
                        shutil.move(default_csv, final_csv)
                else:
                    logger.warning(f"R script failed: {result.stderr}")
                    # Continue anyway as we still have the BED file
            else:
                logger.warning(f"R script not found at {r_script_path}")

            logger.info("Matkit processing completed")
            return True

        except ImportError:
            logger.warning("Robin package not available, skipping matkit processing")
            return False

    except Exception as e:
        logger.error(f"Error running matkit: {e}")
        return False
    finally:
        # Clean up temporary files
        for temp_file in temp_files:
            try:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
            except Exception as e:
                logger.warning(f"Could not remove temporary file {temp_file}: {e}")


def process_bam_file(
    bam_path: str,
    metadata: Dict[str, Any],
    work_dir: str,
    threads: int = 4,
    reference: Optional[str] = None,
) -> MGMTMetadata:
    """Process a single BAM file for MGMT analysis"""
    # print(f"Processing MGMT for BAM file: {bam_path}\n\n")
    logger = logging.getLogger("robin.mgmt")

    logger.info(f"Processing BAM file: {bam_path} in MGMT analysis")
    sample_id = metadata.get("sample_id", "unknown")
    start_time = time.time()

    # Create sample-specific directory
    sample_dir = os.path.join(work_dir, sample_id)
    os.makedirs(sample_dir, exist_ok=True)

    # Get next file number for incremental naming based on existing files
    file_number = _get_next_file_number(sample_dir)

    # Find required files and paths
    mgmt_bed_path = _find_mgmt_bed(work_dir)
    hv_path = _find_hv_path(work_dir)

    mgmt_result = MGMTMetadata(
        sample_id=sample_id, bam_path=bam_path, analysis_timestamp=start_time
    )

    logger.info(f"Starting MGMT analysis for sample: {sample_id}")
    logger.debug(f"BAM file: {os.path.basename(bam_path)}")
    logger.debug(f"Work directory: {work_dir}")
    logger.debug(f"MGMT BED file: {mgmt_bed_path}")
    logger.debug(f"HV path: {hv_path}")
    logger.debug(f"Threads: {threads}")

    try:
        # Step 1: Check if BAM file exists and is readable
        if not os.path.exists(bam_path):
            raise FileNotFoundError(f"BAM file not found: {bam_path}")

        mgmt_result.processing_steps.append("file_validation")

        # Step 2: Check if BAM file contains MGMT spanning reads using preprocessing data
        # This saves time by using the MGMT read detection that was already performed during BAM preprocessing
        has_mgmt_reads = metadata.get("has_mgmt_reads", False)
        mgmt_read_count = metadata.get("mgmt_read_count", 0)

        # Log the MGMT read information from preprocessing
        logger.info(
            f"MGMT read information from preprocessing: has_mgmt_reads={has_mgmt_reads}, count={mgmt_read_count}"
        )

        # Fallback: if preprocessing data is not available, run our own check
        if "has_mgmt_reads" not in metadata:
            logger.warning(
                "MGMT read information not found in preprocessing metadata, running own check"
            )
            has_mgmt_reads = has_reads(bam_path, "chr10", 129467242, 129467244)
            mgmt_read_count = 0  # We don't have the count from the fallback check

        if not has_mgmt_reads:
            mgmt_result.processing_steps.append("no_mgmt_reads")
            mgmt_result.error_message = "No MGMT spanning reads found"
            return mgmt_result

        mgmt_result.processing_steps.append("mgmt_reads_found")

        # Step 3: Extract MGMT region using bedtools
        MGMT_BED: str = mgmt_bed_path

        # Create temporary BAM file for this extraction
        temp_bamfile = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)
        temp_bamfile.close()

        # Extract MGMT region to temporary file
        bedtools_cmd = ["bedtools", "intersect", "-a", bam_path, "-b", MGMT_BED]

        try:
            with open(temp_bamfile.name, "wb") as f:
                result = subprocess.run(
                    bedtools_cmd,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    timeout=300,  # 5 minute timeout
                )

            if result.returncode != 0:
                # Check if the output file was created and has content
                if (
                    not os.path.exists(temp_bamfile.name)
                    or os.path.getsize(temp_bamfile.name) == 0
                ):
                    raise RuntimeError("Bedtools failed to create valid BAM file")

            # Create index for the BAM file
            try:
                pysam.index(temp_bamfile.name, f"{temp_bamfile.name}.bai")
            except Exception as e:
                logger.warning(f"Failed to create BAM index: {e}")
                # Continue without index if it fails

            # Check if the extracted BAM has reads
            if pysam.AlignmentFile(temp_bamfile.name, "rb").count(until_eof=True) > 0:
                # Accumulate reads in the sample's mgmt.bam file
                mgmt_bam_output = os.path.join(sample_dir, "mgmt.bam")

                if os.path.exists(mgmt_bam_output):
                    # Concatenate with existing mgmt.bam
                    temp_holder = tempfile.NamedTemporaryFile(
                        suffix=".bam", delete=False
                    )
                    temp_holder.close()

                    pysam.cat(
                        "-o", temp_holder.name, mgmt_bam_output, temp_bamfile.name
                    )
                    shutil.copy2(temp_holder.name, mgmt_bam_output)
                    try:
                        os.remove(f"{temp_holder.name}.bai")
                        os.remove(f"{temp_bamfile.name}.bai")
                    except FileNotFoundError:
                        pass
                else:
                    # Copy to create new mgmt.bam
                    shutil.copy2(temp_bamfile.name, mgmt_bam_output)
                    try:
                        os.remove(f"{temp_bamfile.name}.bai")
                    except FileNotFoundError:
                        pass

                mgmt_result.processing_steps.append("mgmt_region_extraction")
            else:
                mgmt_result.processing_steps.append("no_mgmt_region_reads")
                # mgmt_result.error_message = "No reads found in extracted MGMT region"
                return mgmt_result

        except subprocess.TimeoutExpired:
            raise RuntimeError("Bedtools extraction timed out")
        except FileNotFoundError:
            raise RuntimeError("bedtools not found in PATH")

        # Step 4: Generate visualization (if methylartist is available)
        temp_files = []  # Track temporary files for cleanup

        try:
            # Sort and index the accumulated mgmt.bam file for methylartist
            sorted_mgmt_bam = os.path.join(sample_dir, "mgmt_sorted.bam")
            # temp_files.append(sorted_mgmt_bam)
            # temp_files.append(f"{sorted_mgmt_bam}.bai")

            # Use safe sort and index function to prevent truncated file issues
            if safely_sort_and_index_bam(
                input_bam=mgmt_bam_output,
                output_bam=sorted_mgmt_bam,
                logger=logger,
                threads=4,
                verify_readable=True,
            ):
                logger.info(
                    f"Sorted and indexed accumulated MGMT BAM: {sorted_mgmt_bam}"
                )
                mgmt_bam_for_plot = sorted_mgmt_bam
            else:
                logger.warning(
                    f"Failed to sort/index accumulated MGMT BAM, using unsorted BAM"
                )
                mgmt_bam_for_plot = mgmt_bam_output

            plot_out = os.path.join(sample_dir, f"{file_number}_mgmt.png")

            # Look for corresponding bed file for site annotations
            bed_file = os.path.join(sample_dir, f"{file_number}_mgmt.bed")
            if not os.path.exists(bed_file):
                # Try alternative naming
                alt_bed = os.path.join(sample_dir, f"{file_number}_mgmt_mgmt.bed")
                bed_file = alt_bed if os.path.exists(alt_bed) else None

            # Use the new safe methylartist wrapper
            methylartist_result = run_methylartist_safely(
                bam_file=mgmt_bam_for_plot,
                output_file=plot_out,
                interval="chr10:129466536-129467536",
                reference=reference,
                logger=logger,
                bed_file=bed_file if bed_file and os.path.exists(bed_file) else None,
            )

            if methylartist_result["success"]:
                mgmt_result.processing_steps.append("visualization")
                logger.info("Visualization complete")
                if methylartist_result["fallback_used"]:
                    logger.info("Visualization completed using fallback parameters")
            else:
                logger.warning(
                    f"Methylartist visualization failed: {methylartist_result['error_message']}"
                )
                # Don't treat this as a fatal error - analysis can continue without visualization

        except Exception as e:
            logger.warning(f"Error during visualization step: {e}")
        finally:
            # Clean up temporary files
            for temp_file in temp_files:
                try:
                    if os.path.exists(temp_file):
                        os.remove(temp_file)
                except Exception as e:
                    logger.warning(f"Could not remove temporary file {temp_file}: {e}")

        # Step 5: Run matkit analysis on accumulated data (if available)
        try:
            logger.info("Running matkit")
            # Use the proper run_matkit function with robin package
            if run_matkit(work_dir, mgmt_bam_output, hv_path, sample_dir, file_number):
                mgmt_result.processing_steps.append("matkit_analysis")
                logger.info("Matkit analysis complete")
            else:
                logger.warning("Matkit analysis failed or robin package not available")

        except Exception as e:
            logger.error(f"Matkit analysis error: {e}")

        # Clean up temporary files after all analysis is complete
        temp_files_to_clean = [temp_bamfile.name]
        if "temp_holder" in locals():
            temp_files_to_clean.append(temp_holder.name)

        for temp_file in temp_files_to_clean:
            try:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
            except FileNotFoundError:
                pass

        # Step 6: Collect results and metadata
        analysis_time = time.time() - start_time
        mgmt_result.results = {
            "analysis_time": analysis_time,
            "mgmt_bam_file": mgmt_bam_output,
            "processing_steps": mgmt_result.processing_steps.copy(),
            "tools_available": {
                "bedtools": shutil.which("bedtools") is not None,
                "robin_package": HVPATH is not None,
                "matkit": HVPATH is not None,  # matkit comes from robin package
                "r": shutil.which("Rscript") is not None,
                "methylartist": shutil.which("methylartist") is not None,
            },
            "hv_path": hv_path,
            "hv_path_source": (
                "robin_module" if HVPATH and hv_path == HVPATH else "fallback"
            ),
            "mgmt_read_count_from_preprocessing": mgmt_read_count,
        }

        mgmt_result.processing_steps.append("analysis_complete")
        logger.info(f"MGMT analysis completed for {sample_id} in {analysis_time:.2f}s")

    except Exception as e:
        error_msg = f"MGMT analysis failed for {sample_id}: {str(e)}"
        logger.error(f"{error_msg}")
        mgmt_result.error_message = error_msg
        mgmt_result.processing_steps.append("analysis_failed")

    return mgmt_result


def run_final_combined_analysis(sample_id: str, work_dir: str) -> bool:
    """Run final combined analysis on all accumulated data for a sample."""
    sample_dir = os.path.join(work_dir, sample_id)
    mgmt_bam_output = os.path.join(sample_dir, "mgmt.bam")

    logger = logging.getLogger("robin.mgmt")
    if not os.path.exists(mgmt_bam_output) or os.path.getsize(mgmt_bam_output) == 0:
        logger.warning(f"No accumulated MGMT BAM file found for sample {sample_id}")
        return False

    logger.info(f"Running final combined analysis for sample {sample_id}")
    try:
        # Run matkit on the accumulated data with "final" naming
        hv_path = _find_hv_path(work_dir)
        if run_matkit(work_dir, mgmt_bam_output, hv_path, sample_dir, "final"):
            logger.info(f"Final combined analysis completed for sample {sample_id}")
            return True
        else:
            logger.warning(f"Final combined analysis failed for sample {sample_id}")
            return False
    except Exception as e:
        logger.error(f"Error in final combined analysis for sample {sample_id}: {e}")
        return False


def process_multiple_bams(bam_paths, metadata_list, work_dir, logger, reference=None):
    """
    Process multiple BAM files for MGMT analysis using aggregated MGMT data.

    This function processes multiple BAM files for the same sample, accumulating
    MGMT reads across all files before performing downstream analysis. This is
    more efficient than processing each BAM file individually and then trying
    to combine results.

    Args:
        bam_paths: List of paths to BAM files
        metadata_list: List of metadata dictionaries (one per BAM file)
        work_dir: Working directory
        logger: Logger instance
        reference: Optional path to reference genome for methylartist

    Returns:
        Dictionary with aggregated MGMT analysis results
    """
    if not bam_paths or not metadata_list:
        raise ValueError("bam_paths and metadata_list must not be empty")

    if len(bam_paths) != len(metadata_list):
        raise ValueError("bam_paths and metadata_list must have the same length")

    # Get sample ID from first metadata (assuming all BAMs are from same sample)
    sample_id = metadata_list[0].get("sample_id", "unknown")

    logger.info(f"🧬 Starting multi-BAM MGMT analysis for sample: {sample_id}")
    logger.info(f"Processing {len(bam_paths)} BAM files for sample {sample_id}")

    # Log essential metadata only
    for i, (bam_path, metadata) in enumerate(zip(bam_paths, metadata_list)):
        logger.debug(f"BAM file {i+1}: {os.path.basename(bam_path)}")

    # Create sample-specific output directory
    sample_dir = os.path.join(work_dir, sample_id)
    os.makedirs(sample_dir, exist_ok=True)

    analysis_result = {
        "sample_id": sample_id,
        "bam_paths": bam_paths,
        "analysis_timestamp": time.time(),
        "mgmt_bam_file": None,
        "processing_steps": [],
        "error_message": None,
        "files_processed": 0,
        "total_files": len(bam_paths),
        "tools_available": {},
        "mgmt_read_count_from_preprocessing": 0,
    }

    try:
        # Find required files and paths
        mgmt_bed_path = _find_mgmt_bed(work_dir)
        hv_path = _find_hv_path(work_dir)

        logger.info(f"Starting MGMT analysis for sample: {sample_id}")
        logger.debug(f"Work directory: {work_dir}")
        logger.debug(f"MGMT BED file: {mgmt_bed_path}")
        logger.debug(f"HV path: {hv_path}")

        # Check if any BAM has MGMT reads
        valid_bam_paths = []
        valid_metadata_list = []
        total_mgmt_read_count = 0

        for bam_path, metadata in zip(bam_paths, metadata_list):
            # Check MGMT reads from preprocessing metadata
            has_mgmt_reads = metadata.get("has_mgmt_reads", False)
            mgmt_read_count = metadata.get("mgmt_read_count", 0)

            if has_mgmt_reads:
                valid_bam_paths.append(bam_path)
                valid_metadata_list.append(metadata)
                total_mgmt_read_count += mgmt_read_count
                logger.debug(
                    f"BAM {os.path.basename(bam_path)}: {mgmt_read_count} MGMT reads"
                )
            else:
                logger.warning(
                    f"No MGMT reads found in BAM file: {os.path.basename(bam_path)}"
                )

        if not valid_bam_paths:
            logger.info(
                f"No MGMT reads found in any BAM files for sample {sample_id} - this is normal"
            )
            analysis_result["processing_steps"].append("no_mgmt_reads_found")
            analysis_result["status"] = "no_mgmt_reads"
            analysis_result["message"] = "No MGMT reads found in any BAM files"
            return analysis_result

        logger.info(
            f"Processing {len(valid_bam_paths)} valid BAM files out of {len(bam_paths)} total"
        )
        analysis_result["files_processed"] = len(valid_bam_paths)
        analysis_result["mgmt_read_count_from_preprocessing"] = total_mgmt_read_count
        analysis_result["processing_steps"].append("mgmt_reads_found")

        # Process each BAM file, accumulating MGMT reads
        mgmt_bam_output = os.path.join(sample_dir, "mgmt.bam")

        for i, (bam_path, metadata) in enumerate(
            zip(valid_bam_paths, valid_metadata_list)
        ):
            logger.info(
                f"Processing BAM file {i+1}/{len(valid_bam_paths)}: {os.path.basename(bam_path)}"
            )

            # Extract MGMT region using bedtools
            temp_bamfile = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)
            temp_bamfile.close()

            # Extract MGMT region to temporary file
            bedtools_cmd = [
                "bedtools",
                "intersect",
                "-a",
                bam_path,
                "-b",
                mgmt_bed_path,
            ]

            try:
                with open(temp_bamfile.name, "wb") as f:
                    result = subprocess.run(
                        bedtools_cmd,
                        stdout=f,
                        stderr=subprocess.PIPE,
                        timeout=300,  # 5 minute timeout
                    )

                if result.returncode != 0:
                    # Check if the output file was created and has content
                    if (
                        not os.path.exists(temp_bamfile.name)
                        or os.path.getsize(temp_bamfile.name) == 0
                    ):
                        logger.warning(
                            f"Bedtools failed to create valid BAM file for {os.path.basename(bam_path)}"
                        )
                        continue

                # Create index for the BAM file
                try:
                    pysam.index(temp_bamfile.name, f"{temp_bamfile.name}.bai")
                except Exception as e:
                    logger.warning(f"Failed to create BAM index: {e}")

                # Check if the extracted BAM has reads
                if (
                    pysam.AlignmentFile(temp_bamfile.name, "rb").count(until_eof=True)
                    > 0
                ):
                    if os.path.exists(mgmt_bam_output):
                        # Concatenate with existing mgmt.bam
                        temp_holder = tempfile.NamedTemporaryFile(
                            suffix=".bam", delete=False
                        )
                        temp_holder.close()

                        pysam.cat(
                            "-o", temp_holder.name, mgmt_bam_output, temp_bamfile.name
                        )
                        shutil.copy2(temp_holder.name, mgmt_bam_output)
                        try:
                            os.remove(f"{temp_holder.name}.bai")
                            os.remove(f"{temp_bamfile.name}.bai")
                        except FileNotFoundError:
                            pass
                    else:
                        # Copy to create new mgmt.bam
                        shutil.copy2(temp_bamfile.name, mgmt_bam_output)
                        try:
                            os.remove(f"{temp_bamfile.name}.bai")
                        except FileNotFoundError:
                            pass

                    logger.debug(
                        f"Accumulated MGMT reads from {os.path.basename(bam_path)}"
                    )
                else:
                    logger.warning(
                        f"No reads found in extracted MGMT region for {os.path.basename(bam_path)}"
                    )

            except subprocess.TimeoutExpired:
                logger.warning(
                    f"Bedtools extraction timed out for {os.path.basename(bam_path)}"
                )
            except FileNotFoundError:
                logger.error("bedtools not found in PATH")
                raise RuntimeError("bedtools not found in PATH")
            finally:
                # Clean up temporary file
                try:
                    if os.path.exists(temp_bamfile.name):
                        os.remove(temp_bamfile.name)
                except FileNotFoundError:
                    pass

        # Check if we have accumulated any MGMT reads
        if not os.path.exists(mgmt_bam_output) or os.path.getsize(mgmt_bam_output) == 0:
            analysis_result["processing_steps"].append("no_mgmt_region_reads")
            analysis_result["status"] = "no_mgmt_reads"
            analysis_result["message"] = "No MGMT reads accumulated from any BAM files"
            return analysis_result

        analysis_result["mgmt_bam_file"] = mgmt_bam_output
        analysis_result["processing_steps"].append("mgmt_region_extraction")

        # Generate visualization (if methylartist is available)
        try:
            # Sort and index the accumulated mgmt.bam file for methylartist
            sorted_mgmt_bam = os.path.join(sample_dir, "mgmt_sorted.bam")

            # Use safe sort and index function to prevent truncated file issues
            if safely_sort_and_index_bam(
                input_bam=mgmt_bam_output,
                output_bam=sorted_mgmt_bam,
                logger=logger,
                threads=4,
                verify_readable=True,
            ):
                logger.info(
                    f"Sorted and indexed accumulated MGMT BAM: {sorted_mgmt_bam}"
                )
                mgmt_bam_for_plot = sorted_mgmt_bam
            else:
                logger.warning(
                    f"Failed to sort/index accumulated MGMT BAM, using unsorted BAM"
                )
                mgmt_bam_for_plot = mgmt_bam_output

            plot_out = os.path.join(sample_dir, "final_mgmt.png")

            # Look for corresponding bed file for site annotations
            bed_file = os.path.join(sample_dir, "final_mgmt.bed")
            if not os.path.exists(bed_file):
                # Try alternative naming
                alt_bed = os.path.join(sample_dir, "final_mgmt_mgmt.bed")
                bed_file = alt_bed if os.path.exists(alt_bed) else None

            # Use the new safe methylartist wrapper
            methylartist_result = run_methylartist_safely(
                bam_file=mgmt_bam_for_plot,
                output_file=plot_out,
                interval="chr10:129466536-129467536",
                reference=reference,
                logger=logger,
                bed_file=bed_file if bed_file and os.path.exists(bed_file) else None,
            )

            if methylartist_result["success"]:
                analysis_result["processing_steps"].append("visualization")
                logger.info("Visualization complete")
                if methylartist_result["fallback_used"]:
                    logger.info("Visualization completed using fallback parameters")
            else:
                logger.warning(
                    f"Methylartist visualization failed: {methylartist_result['error_message']}"
                )
                # Don't treat this as a fatal error - analysis can continue without visualization

        except Exception as e:
            logger.warning(f"Error during visualization step: {e}")

        # Run matkit analysis on accumulated data (if available)
        try:
            logger.info("Running matkit on accumulated MGMT data")
            # Use the proper run_matkit function with robin package
            if run_matkit(work_dir, mgmt_bam_output, hv_path, sample_dir, "final"):
                analysis_result["processing_steps"].append("matkit_analysis")
                logger.info("Matkit analysis complete")
            else:
                logger.warning("Matkit analysis failed or robin package not available")

        except Exception as e:
            logger.error(f"Matkit analysis error: {e}")

        # Collect results and metadata
        analysis_result["tools_available"] = {
            "bedtools": shutil.which("bedtools") is not None,
            "robin_package": HVPATH is not None,
            "matkit": HVPATH is not None,  # matkit comes from robin package
            "r": shutil.which("Rscript") is not None,
            "methylartist": shutil.which("methylartist") is not None,
        }
        analysis_result["hv_path"] = hv_path
        analysis_result["hv_path_source"] = (
            "robin_module" if HVPATH and hv_path == HVPATH else "fallback"
        )

        analysis_result["processing_steps"].append("analysis_complete")
        logger.info(f"Multi-BAM MGMT analysis completed for {sample_id}")
        logger.info(
            f"Files successfully processed: {analysis_result['files_processed']}/{analysis_result['total_files']}"
        )
        logger.info(
            f"Total MGMT reads: {analysis_result['mgmt_read_count_from_preprocessing']}"
        )

        return analysis_result

    except Exception as e:
        logger.error(f"Error in multi-BAM MGMT analysis for {sample_id}: {e}")
        analysis_result["error_message"] = str(e)
        analysis_result["processing_steps"].append("analysis_failed")
        return analysis_result


def mgmt_handler(job, work_dir=None, reference=None):
    """
    Handler function for MGMT analysis jobs.
    This function processes BAM files for MGMT methylation analysis.

    Args:
        job: The workflow job containing file and metadata
        work_dir: Optional base directory for output (defaults to BAM file directory)
        reference: Optional path to reference genome for methylartist
    """
    # Get job-specific logger
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)

    # Check if this is a batched job
    batched_job = job.context.metadata.get("_batched_job")
    if batched_job:
        batch_size = batched_job.get_file_count()
        sample_id = batched_job.get_sample_id()
        batch_id = batched_job.batch_id
        logger.info(
            f"Processing MGMT batch: {batch_size} files for sample '{sample_id}' (batch_id: {batch_id})"
        )

        # Get all filepaths in the batch
        filepaths = batched_job.get_filepaths()

        # Log individual files in the batch
        for i, filepath in enumerate(filepaths):
            logger.info(
                f"  Batch file {i+1}/{batch_size}: {os.path.basename(filepath)}"
            )

        # Prepare metadata list for all BAM files in the batch
        metadata_list = []
        for i, bam_path in enumerate(filepaths):
            # Get metadata from preprocessing for this specific file
            file_metadata = batched_job.contexts[i].metadata.get("bam_metadata", {})

            # Get sample ID from preprocessing results for this specific file
            file_context = batched_job.contexts[i]
            file_sample_id = file_context.get_sample_id()

            # Use the sample ID from the file's context (which should have preprocessing results)
            if file_sample_id != "unknown":
                file_metadata["sample_id"] = file_sample_id
            else:
                file_metadata["sample_id"] = sample_id

            metadata_list.append(file_metadata)

        # Determine work directory for the batch
        if work_dir is None:
            # Default to first BAM file directory
            batch_work_dir = os.path.dirname(filepaths[0])
        else:
            # Use specified work directory, create if it doesn't exist
            os.makedirs(work_dir, exist_ok=True)
            batch_work_dir = work_dir
            logger.debug(f"Using specified work directory: {batch_work_dir}")

        # Process all BAM files in the batch using the new aggregated function
        logger.info(
            f"Processing {batch_size} BAM files as aggregated batch for sample '{sample_id}'"
        )
        batch_result = process_multiple_bams(
            bam_paths=filepaths,
            metadata_list=metadata_list,
            work_dir=batch_work_dir,
            logger=logger,
            reference=reference,
        )

        # Store batch results in job context (maintain compatibility with existing structure)
        job.context.add_metadata(
            "mgmt_analysis",
            {
                "batch_result": batch_result,  # Single aggregated result
                "batch_size": batch_size,
                "sample_id": sample_id,
                "batch_id": batch_id,
                "files_processed": batch_result.get("files_processed", batch_size),
                "total_files": batch_result.get("total_files", batch_size),
            },
        )

        logger.info(
            f"Completed MGMT batch processing: {batch_size} files for sample '{sample_id}'"
        )
        logger.info(
            f"Files successfully processed: {batch_result.get('files_processed', batch_size)}/{batch_result.get('total_files', batch_size)}"
        )

        # Check if this is an expected condition (no MGMT reads) vs actual error
        if batch_result.get("status") == "no_mgmt_reads":
            # This is an expected condition, not an error
            job.context.add_result(
                "mgmt_analysis",
                {
                    "status": "no_mgmt_reads",
                    "sample_id": sample_id,
                    "processing_steps": batch_result.get("processing_steps", []),
                    "message": batch_result.get(
                        "message", "No MGMT reads found in any BAM files"
                    ),
                    "mgmt_read_count_from_preprocessing": batch_result.get(
                        "mgmt_read_count_from_preprocessing", 0
                    ),
                    "files_processed": batch_result.get("files_processed", batch_size),
                    "total_files": batch_result.get("total_files", batch_size),
                },
            )
            logger.info(
                f"MGMT batch analysis completed - {batch_result.get('message', 'No MGMT reads found')}"
            )
        elif batch_result.get("error_message"):
            logger.error(
                f"Batch processing completed with errors: {batch_result['error_message']}"
            )
            job.context.add_error("mgmt_analysis", batch_result["error_message"])
        else:
            logger.info(
                "Batch processing completed successfully with aggregated MGMT analysis"
            )
            job.context.add_result(
                "mgmt_analysis",
                {
                    "status": "success",
                    "sample_id": sample_id,
                    "mgmt_bam_file": batch_result.get("mgmt_bam_file", ""),
                    "processing_steps": batch_result.get("processing_steps", []),
                    "tools_available": batch_result.get("tools_available", {}),
                    "mgmt_read_count_from_preprocessing": batch_result.get(
                        "mgmt_read_count_from_preprocessing", 0
                    ),
                    "files_processed": batch_result.get("files_processed", batch_size),
                    "total_files": batch_result.get("total_files", batch_size),
                },
            )

        return

    else:
        # Single file processing (backward compatibility)
        try:
            bam_path = job.context.filepath
            logger.info(f"Starting MGMT analysis for: {os.path.basename(bam_path)}")

            # Get metadata from preprocessing
            bam_metadata = job.context.metadata.get("bam_metadata", {})

            # Determine work directory
            if work_dir is None:
                # Default to BAM file directory
                work_dir = os.path.dirname(bam_path)
            else:
                # Use specified work directory, create if it doesn't exist
                os.makedirs(work_dir, exist_ok=True)
                logger.debug(f"Using specified work directory: {work_dir}")

            # Process the BAM file
            mgmt_result = process_bam_file(
                bam_path, bam_metadata, work_dir, reference=reference
            )

            # Store results in job context
            job.context.add_metadata("mgmt_analysis", mgmt_result.results)
            job.context.add_metadata(
                "mgmt_processing_steps", mgmt_result.processing_steps
            )

            if mgmt_result.error_message:
                # Check if this is an expected condition (no MGMT reads) vs actual error
                if mgmt_result.error_message == "No MGMT spanning reads found":
                    # This is an expected condition, not an error
                    job.context.add_result(
                        "mgmt_analysis",
                        {
                            "status": "no_mgmt_reads",
                            "sample_id": mgmt_result.sample_id,
                            "processing_steps": mgmt_result.processing_steps,
                            "message": mgmt_result.error_message,
                            "mgmt_read_count_from_preprocessing": mgmt_result.results.get(
                                "mgmt_read_count_from_preprocessing", 0
                            ),
                        },
                    )
                    logger.info(
                        f"MGMT analysis completed - {mgmt_result.error_message}"
                    )
                else:
                    # This is an actual error
                    job.context.add_error("mgmt_analysis", mgmt_result.error_message)
                    logger.error(f"MGMT analysis failed: {mgmt_result.error_message}")
            else:
                job.context.add_result(
                    "mgmt_analysis",
                    {
                        "status": "success",
                        "sample_id": mgmt_result.sample_id,
                        "analysis_time": mgmt_result.results.get("analysis_time", 0),
                        "mgmt_bam_file": mgmt_result.results.get("mgmt_bam_file", ""),
                        "processing_steps": mgmt_result.processing_steps,
                        "tools_available": mgmt_result.results.get(
                            "tools_available", {}
                        ),
                        "mgmt_read_count_from_preprocessing": mgmt_result.results.get(
                            "mgmt_read_count_from_preprocessing", 0
                        ),
                    },
                )
                logger.info(f"MGMT analysis complete for {os.path.basename(bam_path)}")
                logger.info(f"Sample ID: {mgmt_result.sample_id}")
                logger.info(
                    f"Analysis time: {mgmt_result.results.get('analysis_time', 0):.2f}s"
                )
                logger.debug(
                    f"Processing steps: {', '.join(mgmt_result.processing_steps)}"
                )
                logger.debug(
                    f"Output directory: {os.path.dirname(mgmt_result.results.get('mgmt_bam_file', ''))}"
                )

        except Exception as e:
            job.context.add_error("mgmt_analysis", str(e))
            logger.error(f"Error in MGMT analysis for {job.context.filepath}: {e}")
