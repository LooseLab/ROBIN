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
from datetime import datetime
import numpy as np
import polars as pl
from contextlib import contextmanager
import tempfile
import json

from robin.analysis.utilities.modkit_merge import (
    aggregate_modkit_by_site,
)
from robin.analysis.utilities.panel_filter import get_panel_interval_index

# Suppress pkg_resources deprecation warnings from sorted_nearest
warnings.filterwarnings(
    "ignore", message="pkg_resources is deprecated", category=UserWarning
)
# Suppress matplotlib tight_layout warnings
warnings.filterwarnings(
    "ignore", message="The figure layout has changed to tight", category=UserWarning
)

try:
    from robin.analysis.utilities.mnp_flex import APIClient as MnpFlexClient #ToDo: Maintain to future integration.
except ImportError as e:
    logging.warning(f"Some dependencies not available: {e}")


# Simple cross-process file lock using POSIX flock when available (no-op on unsupported platforms)
try:
    import fcntl  # type: ignore
except Exception:
    fcntl = None  # type: ignore


@contextmanager
def exclusive_file_lock(lock_path: str):
    """Provide a blocking exclusive lock on a file path.

    The lock is advisory and therefore coordinates callers which use the same
    lock path. On POSIX, ``flock(LOCK_EX)`` waits until the current holder has
    finished. It is a best-effort no-op on platforms without ``fcntl``.
    """
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


# Backward-compatible private name for callers outside this module.
_exclusive_file_lock = exclusive_file_lock


def _filter_panel_sites(
    frame: pl.DataFrame,
    panel_index,
) -> pl.DataFrame:
    """Keep rows whose genomic positions are covered by the panel.

    This is the merge-stage equivalent of the early matkit filter.  It uses
    the cached interval arrays directly rather than converting the frame to
    Pandas and constructing a PyRanges object for every BAM output file.
    A boolean mask also means overlapping panel intervals cannot duplicate a
    site, matching ``panel_intersection_keys(...).drop_duplicates()``.
    """
    if frame.is_empty():
        return frame

    chrom_values = frame.get_column("chrom").to_numpy()
    positions = frame.get_column("chromStart").to_numpy()
    keep = np.zeros(frame.height, dtype=bool)

    for chrom in np.unique(chrom_values):
        row_indices = np.flatnonzero(chrom_values == chrom)
        intervals = panel_index.intervals_by_chrom.get(str(chrom))
        if intervals is None or len(row_indices) == 0:
            continue

        starts, _ends, prefix_max_ends = intervals
        row_positions = positions[row_indices]
        interval_indices = np.searchsorted(
            starts, row_positions, side="right"
        ) - 1
        valid = interval_indices >= 0
        if np.any(valid):
            valid_rows = row_indices[valid]
            valid_intervals = interval_indices[valid]
            keep[valid_rows] = (
                prefix_max_ends[valid_intervals] > row_positions[valid]
            )

    return frame.filter(pl.Series("panel_keep", keep))


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

        # Load the cached half-open interval index once.  The index is shared
        # in-process across merge calls and persisted on disk for later jobs.
        # It replaces the previous Pandas/PyRanges cache and intersection path.
        panel_index = get_panel_interval_index(
            filter_bed_file,
            cache_dir=sample_output_dir,
        )

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

                # Keep panel sites without converting through Pandas/PyRanges.
                # The interval mask is unique per input row, so overlapping
                # panel intervals cannot duplicate counts.
                filt = _filter_panel_sites(pl_df, panel_index)
                if filt.is_empty():
                    continue
                new_frames.append(filt)
            except Exception as e:
                logging.error(f"Error processing file {bed}: {str(e)}")
                continue

        if not new_frames:
            logging.warning("No valid data to merge after filtering")
            return

        # Combine new data (all Polars)
        new_df = pl.concat(new_frames) if len(new_frames) > 1 else new_frames[0]
        # A batch can contain multiple BAMs, and panel overlaps or legacy inputs can
        # also yield repeated rows. Collapse before the first write as well as on
        # subsequent updates.
        new_df = aggregate_modkit_by_site(new_df)

        # Critical section: read/merge/write parquet guarded by an exclusive lock
        lock_path = f"{output_file}.lock"
        with exclusive_file_lock(lock_path):
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

                # Aggregate by complete site identity and derive percent_modified
                # from the summed counts, rather than averaging per-BAM percentages.
                grouped = aggregate_modkit_by_site(combined)

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
        # Release temporary frame references promptly for large batches.
        gc.collect()
