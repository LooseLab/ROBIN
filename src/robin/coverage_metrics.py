"""Shared coverage calculations used by the GUI, TSV export and PDF reports."""

from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

import pandas as pd

from robin.reference_contigs import filter_dataframe_contigs


def filter_coverage_frames(
    coverage_df: pd.DataFrame,
    target_df: pd.DataFrame,
    scope: Optional[str] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Apply the configured reference-contig scope to coverage inputs.

    Files produced by different workflow versions do not all contain contig
    columns.  In that case the input is retained unchanged, matching the
    historical report behaviour while still allowing current files to be
    filtered consistently.
    """
    filtered_coverage = coverage_df
    if "#rname" in filtered_coverage.columns:
        filtered_coverage = filter_dataframe_contigs(
            filtered_coverage, "#rname", scope
        )

    filtered_target = target_df
    if "chrom" in filtered_target.columns:
        filtered_target = filter_dataframe_contigs(filtered_target, "chrom", scope)

    return filtered_coverage, filtered_target


def calculate_coverage_metrics(
    coverage_df: pd.DataFrame,
    target_df: pd.DataFrame,
    scope: Optional[str] = None,
) -> Dict[str, Any]:
    """Calculate global coverage, target coverage and enrichment.

    The returned data contains the filtered frames as well as the numeric
    metrics so report plots and summary values cannot silently use different
    populations of contigs.
    """
    filtered_coverage, filtered_target = filter_coverage_frames(
        coverage_df, target_df, scope
    )

    global_denominator = float(filtered_coverage["endpos"].sum())
    target_frame = filtered_target
    if "length" not in target_frame.columns:
        target_frame = target_frame.copy()
        target_frame["length"] = (
            target_frame["endpos"] - target_frame["startpos"] + 1
        ).astype(float)
    target_denominator = float(target_frame["length"].sum())

    global_coverage = (
        float(filtered_coverage["covbases"].sum()) / global_denominator
        if global_denominator > 0
        else None
    )
    target_coverage = (
        float(target_frame["bases"].sum()) / target_denominator
        if target_denominator > 0
        else None
    )
    enrichment = (
        target_coverage / global_coverage
        if global_coverage is not None
        and target_coverage is not None
        and global_coverage > 0
        else None
    )

    return {
        "global_coverage": global_coverage,
        "target_coverage": target_coverage,
        "enrichment": enrichment,
        "coverage_df": filtered_coverage,
        "target_df": target_frame,
    }
