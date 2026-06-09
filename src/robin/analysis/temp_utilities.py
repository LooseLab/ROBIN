#!/usr/bin/env python3
"""
Temporary utilities for robin BAM to parquet conversion.

This module provides utilities for merging modkit files and creating parquet files.
"""

import warnings
from typing import List, Optional
import gc
import os
import logging
import pickle
from datetime import datetime
import pandas as pd
import polars as pl
from contextlib import contextmanager
import tempfile
import json

# Suppress pkg_resources deprecation warnings from sorted_nearest
warnings.filterwarnings(
    "ignore", message="pkg_resources is deprecated", category=UserWarning
)
# Suppress matplotlib tight_layout warnings
warnings.filterwarnings(
    "ignore", message="The figure layout has changed to tight", category=UserWarning
)

try:
    import pyranges as pr
    from robin.analysis.utilities.mnp_flex import APIClient as MnpFlexClient #ToDo: Maintain to future integration.
except ImportError as e:
    logging.warning(f"Some dependencies not available: {e}")



# Simple cross-process file lock using POSIX flock when available (no-op on unsupported platforms)
try:
    import fcntl  # type: ignore
except Exception:
    fcntl = None  # type: ignore


@contextmanager
def _exclusive_file_lock(lock_path: str):
    """Provide an exclusive lock on a file path. Best-effort no-op when flock is unavailable."""
    lock_dir = os.path.dirname(lock_path) or "."
    os.makedirs(lock_dir, exist_ok=True)
    lock_file = open(lock_path, "a+")
    if fcntl is not None:
        try:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
        except Exception:
            # Continue without locking if flock fails
            pass
    try:
        yield
    finally:
        if fcntl is not None:
            try:
                fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)
            except Exception:
                pass
        lock_file.close()


def merge_modkit_files(
    new_files: List[str],
    existing_file: str,
    output_file: str,
    filter_bed_file: str,
    sample_id: str,
    output_dir: str,
    mnpflex_config: Optional[dict],
    num_bam_files_seen: int,
) -> None:
    """
    Merge modkit files with optimized column set and improved caching.

    This function uses only essential columns to reduce memory usage and processing time:
    - chrom, chromStart: Required for all classifiers
    - mod_code: Required for Sturgeon filtering
    - strand: Required for proper aggregation
    - valid_cov: Required for coverage calculation
    - percent_modified: Primary methylation data (required)
    - n_mod, n_canonical: Required for modification counts

    Args:
        new_files (List[str]): List of new modkit files to merge
        existing_file (str): Path to existing parquet file
        output_file (str): Path to output parquet file
        filter_bed_file (str): Path to BED file for filtering
        sample_id (str): Sample ID for organizing output files
        output_dir (str): Base output directory
        mnpflex_config (Optional[dict]): Configuration for MNP-FLEX integration
        num_bam_files_seen (int): Number of BAM files being processed
    """
    # Create sample-specific output directory
    sample_output_dir = os.path.join(output_dir, sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)

    # Define optimized schema with only essential columns
    essential_cols = [
        "chrom",
        "chromStart",
        "mod_code",
        "strand",
        "valid_cov",
        "percent_modified",
        "n_mod",
        "n_canonical",
    ]

    categorical_cols = ["chrom", "mod_code", "strand"]
    unsigned_int_cols = ["chromStart", "valid_cov", "n_mod", "n_canonical"]
    float_cols = ["percent_modified"]

    try:
        # Track cumulative BAM file count
        cumulative_bam_file_count = 0

        # Check if existing file has metadata about BAM file count
        if os.path.exists(existing_file):
            try:
                # Read existing metadata
                metadata_file = existing_file.replace(".parquet", "_metadata.json")
                if os.path.exists(metadata_file):
                    with open(metadata_file, "r") as f:
                        metadata = json.load(f)
                        cumulative_bam_file_count = metadata.get("bam_file_count", 0)
                        logging.info(
                            f"Found existing metadata with {cumulative_bam_file_count} BAM files"
                        )
            except Exception as e:
                logging.warning(f"Could not read existing metadata: {str(e)}")
                cumulative_bam_file_count = 0

        # Add the number of new BAM files being processed
        cumulative_bam_file_count += num_bam_files_seen

        logging.info(
            f"Total cumulative BAM files contributing to parquet: {cumulative_bam_file_count} (added {num_bam_files_seen} new files)"
        )

        
        # Cache or build PyRanges filter with improved caching
        # Use distinct cache for .txt (1-based converted) vs .gz (0-based) to avoid stale data
        cache_suffix = "_1based" if filter_bed_file.endswith(".txt") else ""
        cache_path = os.path.join(
            sample_output_dir,
            f"{os.path.basename(filter_bed_file)}{cache_suffix}.pgr_cache",
        )
        if os.path.exists(cache_path):
            with open(cache_path, "rb") as f:
                filter_ranges = pickle.load(f)
        else:
            comp = "gzip" if filter_bed_file.endswith(".gz") else None
            if filter_bed_file.endswith(".txt"):
                # parquet_filter.txt format: header "chr start end IlmnID" (spaces/tabs)
                # Illumina-derived files use 1-based coordinates; BED/PyRanges expect 0-based.
                bed_df = pd.read_csv(
                    filter_bed_file,
                    sep=r"\s+",
                    header=0,
                    dtype=str,
                )
                # Map to expected column names (handle chr/start/end)
                rename = {}
                for c in bed_df.columns:
                    if c.lower() == "chr":
                        rename[c] = "Chromosome"
                    elif c.lower() == "start":
                        rename[c] = "Start"
                    elif c.lower() == "end":
                        rename[c] = "End"
                bed_df = bed_df.rename(columns=rename)
                bed_df["Start"] = bed_df["Start"].astype(int) - 1  # 1-based -> 0-based
                bed_df["End"] = bed_df["End"].astype(int)  # 1-based end inclusive -> 0-based exclusive
            else:
                bed_df = pd.read_csv(
                    filter_bed_file,
                    sep="\t",
                    header=None,
                    names=["Chromosome", "Start", "End", "cg_label"],
                    compression=comp,
                    dtype={"Chromosome": str},
                )
            filter_ranges = pr.PyRanges(bed_df[["Chromosome", "Start", "End"]])
            with open(cache_path, "wb") as f:
                pickle.dump(filter_ranges, f)

        # BEDMethyl 18-column names for modkit output (tab-separated)
        full_cols = [
            "chrom",
            "chromStart",
            "chromEnd",
            "mod_code",
            "score_bed",
            "strand",
            "thickStart",
            "thickEnd",
            "color",
            "valid_cov",
            "percent_modified",
            "n_mod",
            "n_canonical",
            "n_othermod",
            "n_delete",
            "n_fail",
            "n_diff",
            "n_nocall",
        ]

        # Process new files: read with Polars, filter to CPG set, then combine
        new_frames: List[pl.DataFrame] = []
        for bed in new_files:
            try:
                if bed.endswith(".parquet"):
                    # Matkit wrote 8-column parquet; read directly (no CSV parse, types already correct)
                    pl_df = pl.read_parquet(bed, columns=essential_cols)
                    # Matkit may write chrom/mod_code/strand as binary; decode to Utf8.
                    # Column-wise decode via list comp avoids per-row map_elements overhead.
                    for c in categorical_cols:
                        if pl_df.schema.get(c) == pl.Binary:
                            raw_vals = pl_df[c].to_list()
                            decoded = [
                                (
                                    bytes(x).decode("utf-8", errors="replace")
                                    if x is not None
                                    else None
                                )
                                for x in raw_vals
                            ]
                            pl_df = pl_df.with_columns(
                                pl.Series(c, decoded, dtype=pl.Utf8).alias(c)
                            )
                else:
                    # Legacy BEDMethyl text: read CSV, select essential columns, cast
                    pl_df = pl.read_csv(
                        bed,
                        separator="\t",
                        has_header=False,
                        new_columns=full_cols,
                        infer_schema_length=0,
                    ).select(essential_cols)
                    pl_df = pl_df.with_columns(
                        [pl.col(c).cast(pl.UInt32, strict=False) for c in unsigned_int_cols]
                        + [pl.col(c).cast(pl.Float32, strict=False) for c in float_cols]
                    )

                # Skip empty files (e.g. parquet from a batch where no reads passed the QS filter)
                if pl_df.is_empty():
                    continue

                # Build PyRanges for intersection (only need chrom/start/end)
                pr_df = pl_df.rename({"chrom": "Chromosome", "chromStart": "Start"}).with_columns(
                    (pl.col("Start") + 1).alias("End")
                )
                gr = pr.PyRanges(pr_df.to_pandas()[["Chromosome", "Start", "End"]])
                inter = gr.intersect(filter_ranges).df

                # Join in Polars: keep only rows whose (chrom, chromStart) are in the intersection.
                # Cast join keys to match pl_df (str + UInt32); from_pandas can infer Categorical for chrom.
                inter_pl = pl.from_pandas(
                    inter[["Chromosome", "Start"]].rename(
                        columns={"Chromosome": "chrom", "Start": "chromStart"}
                    )
                ).with_columns(
                    pl.col("chrom").cast(pl.Utf8),
                    pl.col("chromStart").cast(pl.UInt32),
                )
                filt = pl_df.join(inter_pl, on=["chrom", "chromStart"], how="inner")
                new_frames.append(filt)
            except Exception as e:
                logging.error(f"Error processing file {bed}: {str(e)}")
                continue

        if not new_frames:
            logging.warning("No valid data to merge after filtering")
            return

        # Combine new data (all Polars)
        new_df = pl.concat(new_frames) if len(new_frames) > 1 else new_frames[0]

        # Critical section: read/merge/write parquet guarded by an exclusive lock
        lock_path = f"{output_file}.lock"
        with _exclusive_file_lock(lock_path):
            # Use a shared StringCache during the merge to avoid costly categorical re-encodings
            with pl.StringCache():
                # If no existing file, just save the new data
                if not os.path.exists(existing_file):
                    pl_df = new_df
                    # Atomic write: write to temp then replace
                    fd, tmp_path = tempfile.mkstemp(
                        prefix="lj_parquet_",
                        suffix=".parquet",
                        dir=os.path.dirname(output_file) or ".",
                    )
                    os.close(fd)
                    try:
                        pl_df.write_parquet(tmp_path)
                        os.replace(tmp_path, output_file)
                    finally:
                        try:
                            if os.path.exists(tmp_path):
                                os.remove(tmp_path)
                        except Exception:
                            pass

                    # Save metadata with cumulative BAM file count (atomic)
                    metadata = {
                        "bam_file_count": cumulative_bam_file_count,
                        "last_updated": datetime.now().isoformat(),
                        "sample_id": sample_id,
                        "files_added_in_this_update": num_bam_files_seen,
                        "column_format": "optimized",  # Mark as optimized format
                    }
                    metadata_file = output_file.replace(".parquet", "_metadata.json")
                    fdm, tmp_meta = tempfile.mkstemp(
                        prefix="lj_meta_",
                        suffix=".json",
                        dir=os.path.dirname(metadata_file) or ".",
                    )
                    os.close(fdm)
                    try:
                        with open(tmp_meta, "w") as f:
                            json.dump(metadata, f, indent=2)
                        os.replace(tmp_meta, metadata_file)
                    finally:
                        try:
                            if os.path.exists(tmp_meta):
                                os.remove(tmp_meta)
                        except Exception:
                            pass

                    logging.info(
                        f"Created new optimized parquet file with {cumulative_bam_file_count} cumulative BAM files"
                    )
                    return

                # Process existing data in chunks
                existing_df = pl.scan_parquet(existing_file)
                # New data is already Polars
                pl_new_df = new_df

                # Convert existing data to regular DataFrame for concatenation
                existing_df = existing_df.collect()

                # Check if existing file is in old format (18 columns) or new format (8 columns)
                is_old_format = (
                    len(existing_df.columns) > 10
                )  # More than 10 columns indicates old format

                if is_old_format:
                    # Convert old format to new format by selecting only essential columns
                    logging.info(
                        "Converting existing file from old format to optimized format"
                    )
                    existing_df = existing_df.select(essential_cols)

                # Ensure consistent data types and column order
                for c in categorical_cols:
                    if c in existing_df.columns and c in pl_new_df.columns:
                        existing_df = existing_df.with_columns(
                            pl.col(c).cast(pl.Categorical)
                        )
                        pl_new_df = pl_new_df.with_columns(
                            pl.col(c).cast(pl.Categorical)
                        )

                # Ensure columns are in the same order
                pl_new_df = pl_new_df.select(existing_df.columns)

                # Validate column names match
                if set(existing_df.columns) != set(pl_new_df.columns):
                    missing_cols = set(existing_df.columns) - set(pl_new_df.columns)
                    extra_cols = set(pl_new_df.columns) - set(existing_df.columns)
                    raise ValueError(
                        f"Column mismatch: missing {missing_cols}, extra {extra_cols}"
                    )

                # Combine existing and new data
                combined = pl.concat([existing_df, pl_new_df])

                # Define aggregation expressions for essential columns only
                exprs = [
                    pl.first("mod_code"),
                    pl.first("strand"),
                    pl.mean("percent_modified").alias("percent_modified"),
                    *[
                        pl.sum(c).alias(c)
                        for c in ["valid_cov", "n_mod", "n_canonical"]
                    ],
                ]

                # Perform groupby and aggregation
                grouped = combined.group_by(["chrom", "chromStart"]).agg(exprs)

                # Atomic write: write to temp then replace
                fd2, tmp_out = tempfile.mkstemp(
                    prefix="lj_parquet_",
                    suffix=".parquet",
                    dir=os.path.dirname(output_file) or ".",
                )
                os.close(fd2)
                try:
                    grouped.write_parquet(tmp_out)
                    os.replace(tmp_out, output_file)
                finally:
                    try:
                        if os.path.exists(tmp_out):
                            os.remove(tmp_out)
                    except Exception:
                        pass

                # Save metadata with updated cumulative BAM file count (atomic)
                metadata = {
                    "bam_file_count": cumulative_bam_file_count,
                    "last_updated": datetime.now().isoformat(),
                    "sample_id": sample_id,
                    "files_added_in_this_update": num_bam_files_seen,
                    "column_format": "optimized",  # Mark as optimized format
                }
                metadata_file = output_file.replace(".parquet", "_metadata.json")
                fd3, tmp_meta2 = tempfile.mkstemp(
                    prefix="lj_meta_",
                    suffix=".json",
                    dir=os.path.dirname(metadata_file) or ".",
                )
                os.close(fd3)
                try:
                    with open(tmp_meta2, "w") as f:
                        json.dump(metadata, f, indent=2)
                    os.replace(tmp_meta2, metadata_file)
                finally:
                    try:
                        if os.path.exists(tmp_meta2):
                            os.remove(tmp_meta2)
                    except Exception:
                        pass

                logging.info(
                    f"Updated optimized parquet file with {cumulative_bam_file_count} cumulative BAM files (added {num_bam_files_seen} in this update)"
                )

                logging.debug(
                    f"Merged with optimized Polars and cache saved to: {output_file}"
                )

    except Exception as e:
        logging.error(f"Error in merge_modkit_files: {str(e)}")
        raise
    finally:
        # Cleanup temporary files if needed
        if os.path.exists(cache_path) and not os.path.exists(filter_bed_file):
            try:
                os.remove(cache_path)
            except Exception as e:
                logging.error(f"Error removing cache file: {str(e)}")
        gc.collect()
