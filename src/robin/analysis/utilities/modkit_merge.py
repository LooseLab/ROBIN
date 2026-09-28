"""
Reference-site merge helpers for cumulative modkit parquet (no PyRanges dependency).
"""

from __future__ import annotations

from typing import List

import pandas as pd
import polars as pl

# BEDMethyl/modkit site identity relative to the reference.
MERGE_SITE_KEYS: List[str] = ["chrom", "chromStart", "strand", "mod_code"]


def panel_intersection_keys(inter_df: pd.DataFrame) -> pd.DataFrame:
    """
    Unique (chrom, chromStart) hits from a PyRanges intersect result.

    Multiple panel intervals can overlap the same modkit position; dedupe before
    joining so per-site counts are not duplicated.
    """
    keys = inter_df.rename(columns={"Chromosome": "chrom", "Start": "chromStart"})[
        ["chrom", "chromStart"]
    ]
    return keys.drop_duplicates()


def aggregate_modkit_by_site(combined: pl.DataFrame) -> pl.DataFrame:
    """
    Collapse rows for the same reference site by summing counts and recomputing
    percent_modified from n_mod / valid_cov (not averaging stored percentages).
    """
    grouped = combined.group_by(MERGE_SITE_KEYS).agg(
        [
            pl.sum("valid_cov").alias("valid_cov"),
            pl.sum("n_mod").alias("n_mod"),
            pl.sum("n_canonical").alias("n_canonical"),
        ]
    )
    grouped = grouped.with_columns(
        pl.when(pl.col("valid_cov") > 0)
        .then(
            (pl.col("n_mod").cast(pl.Float64) / pl.col("valid_cov").cast(pl.Float64))
            * 100.0
        )
        .otherwise(0.0)
        .round(2)
        .cast(pl.Float32)
        .alias("percent_modified")
    )
    return grouped.select(
        [
            "chrom",
            "chromStart",
            "mod_code",
            "strand",
            "valid_cov",
            "percent_modified",
            "n_mod",
            "n_canonical",
        ]
    )
