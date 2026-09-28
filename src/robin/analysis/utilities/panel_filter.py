"""
CpG panel interval index for early methylation site filtering (parquet_filter.txt).
"""

from __future__ import annotations

import logging
import os
import pickle
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# In-process cache keyed by real path to the panel file.
_PANEL_INDEX_MEM_CACHE: Dict[str, "PanelIntervalIndex"] = {}

IntervalArrays = Tuple[np.ndarray, np.ndarray, np.ndarray]


@dataclass(frozen=True)
class PanelIntervalIndex:
    """
    Per-chromosome sorted half-open intervals [Start, End) for fast point lookup.

    Matches merge-time PyRanges intersect semantics on (chrom, 0-based position).
    """

    intervals_by_chrom: Dict[str, IntervalArrays]
    # Candidate CpG anchors derived from panel intervals. ``None`` means the
    # panel contains a broad interval for that chromosome and callers should
    # use their unrestricted/reference-derived path instead.
    cpg_candidates_by_chrom: Dict[str, Optional[np.ndarray]]

    def contains(self, chrom: str, pos: int) -> bool:
        """Return True if pos falls in any panel interval on chrom."""
        intervals = self.intervals_by_chrom.get(chrom)
        if intervals is None:
            return False
        starts, ends, prefix_max_ends = intervals
        n = len(starts)
        if n == 0:
            return False
        idx = int(np.searchsorted(starts, pos, side="right")) - 1
        # Any interval that can contain ``pos`` must start at or before it.
        # The prefix maximum lets us answer whether any such interval extends
        # past ``pos`` without scanning backwards through every preceding panel
        # interval for non-panel positions.
        return idx >= 0 and bool(prefix_max_ends[idx] > pos)

    def cpg_candidates(self, chrom: str, start: int, end: int) -> Optional[np.ndarray]:
        """Return candidate CpG anchors in ``[start, end)`` for this panel.

        Candidates are deliberately a superset: the caller validates that each
        candidate is an actual reference ``CG``. This avoids fetching whole
        chromosomes merely to discover CpGs that can never be emitted.
        """
        if chrom not in self.cpg_candidates_by_chrom:
            return np.empty(0, dtype=np.int32)
        candidates = self.cpg_candidates_by_chrom[chrom]
        if candidates is None:
            return None
        lo = int(np.searchsorted(candidates, start, side="left"))
        hi = int(np.searchsorted(candidates, end, side="left"))
        return candidates[lo:hi]


def panel_filter_cache_path(filter_bed_file: str, cache_dir: str) -> str:
    # Versioned suffix prevents reuse of caches created when .txt starts were
    # incorrectly decremented as though the bundled panel were 1-based.
    # v4 adds panel-derived CpG candidate anchors to the cached index.
    suffix = "_bed0_v4" if filter_bed_file.endswith(".txt") else ""
    return os.path.join(
        cache_dir,
        f"{os.path.basename(filter_bed_file)}{suffix}.interval_index.pkl",
    )


def load_panel_bed_dataframe(filter_bed_file: str) -> pd.DataFrame:
    """
    Load panel BED intervals as Chromosome / Start / End (0-based BED coordinates).

    parquet_filter.txt is a headered BED-like file using 0-based, half-open
    coordinates, consistent with the bundled classifier BED resources.
    """
    comp = "gzip" if filter_bed_file.endswith(".gz") else None
    if filter_bed_file.endswith(".txt"):
        bed_df = pd.read_csv(
            filter_bed_file,
            sep=r"\s+",
            header=0,
            dtype=str,
        )
        col_map = {"chr": "Chromosome", "start": "Start", "end": "End"}
        rename = {c: col_map[c.lower()] for c in bed_df.columns if c.lower() in col_map}
        bed_df = bed_df.rename(columns=rename)
        bed_df["Start"] = bed_df["Start"].astype(int)
        bed_df["End"] = bed_df["End"].astype(int)
    else:
        bed_df = pd.read_csv(
            filter_bed_file,
            sep="\t",
            header=None,
            names=["Chromosome", "Start", "End", "cg_label"],
            compression=comp,
            dtype={"Chromosome": str},
        )
    return bed_df[["Chromosome", "Start", "End"]]


def build_panel_interval_index(filter_bed_file: str) -> PanelIntervalIndex:
    bed_df = load_panel_bed_dataframe(filter_bed_file)
    intervals_by_chrom: Dict[str, IntervalArrays] = {}
    cpg_candidates_by_chrom: Dict[str, Optional[np.ndarray]] = {}
    for chrom, grp in bed_df.groupby("Chromosome", sort=False):
        starts = grp["Start"].to_numpy(dtype=np.int32, copy=True)
        ends = grp["End"].to_numpy(dtype=np.int32, copy=True)
        order = np.argsort(starts, kind="mergesort")
        starts = starts[order]
        ends = ends[order]
        intervals_by_chrom[str(chrom)] = (
            starts,
            ends,
            np.maximum.accumulate(ends),
        )
        lengths = ends - starts
        if len(lengths) and int(lengths.max()) <= 1000:
            # For an interval [start, end), a CpG anchor can be any output
            # position in the interval or one base before it (reverse-strand
            # cytosine). The production panel has one- and two-base intervals,
            # so this remains compact while preserving interval semantics.
            parts = [np.arange(s - 1, e, dtype=np.int32) for s, e in zip(starts, ends)]
            candidates = np.unique(np.concatenate(parts))
            cpg_candidates_by_chrom[str(chrom)] = candidates[candidates >= 0]
        else:
            cpg_candidates_by_chrom[str(chrom)] = None
    logger.info(
        "Built panel interval index: %s chromosomes, %s intervals from %s",
        len(intervals_by_chrom),
        len(bed_df),
        filter_bed_file,
    )
    return PanelIntervalIndex(
        intervals_by_chrom=intervals_by_chrom,
        cpg_candidates_by_chrom=cpg_candidates_by_chrom,
    )


def get_panel_interval_index(
    filter_bed_file: str,
    cache_dir: Optional[str] = None,
) -> PanelIntervalIndex:
    """Load or build a panel index, with optional on-disk pickle cache in cache_dir."""
    real_path = os.path.realpath(filter_bed_file)
    if real_path in _PANEL_INDEX_MEM_CACHE:
        return _PANEL_INDEX_MEM_CACHE[real_path]

    cache_path = (
        panel_filter_cache_path(filter_bed_file, cache_dir) if cache_dir else None
    )
    if cache_path and os.path.exists(cache_path):
        with open(cache_path, "rb") as f:
            index = pickle.load(f)
        logger.info("Loaded panel interval index from cache: %s", cache_path)
        _PANEL_INDEX_MEM_CACHE[real_path] = index
        return index

    index = build_panel_interval_index(filter_bed_file)
    if cache_path:
        os.makedirs(os.path.dirname(cache_path) or ".", exist_ok=True)
        with open(cache_path, "wb") as f:
            pickle.dump(index, f)
        logger.info("Wrote panel interval index cache: %s", cache_path)

    _PANEL_INDEX_MEM_CACHE[real_path] = index
    return index


def _load_panel_pyranges(filter_bed_file: str, cache_dir: Optional[str] = None):
    """Load or build cached PyRanges for panel intersect (same cache as merge path)."""
    import pyranges as pr

    cache_suffix = "_bed0_v2" if filter_bed_file.endswith(".txt") else ""
    cache_dir = cache_dir or "."
    cache_path = os.path.join(
        cache_dir,
        f"{os.path.basename(filter_bed_file)}{cache_suffix}.pgr_cache",
    )
    if os.path.exists(cache_path):
        with open(cache_path, "rb") as f:
            return pickle.load(f)
    bed_df = load_panel_bed_dataframe(filter_bed_file)
    filter_ranges = pr.PyRanges(bed_df)
    os.makedirs(cache_dir, exist_ok=True)
    with open(cache_path, "wb") as f:
        pickle.dump(filter_ranges, f)
    return filter_ranges


def filter_bedmethyl_to_panel(
    df: pd.DataFrame,
    filter_bed_file: str,
    cache_dir: Optional[str] = None,
) -> pd.DataFrame:
    """
    Keep only rows whose (chrom, chromStart) intersect panel intervals.

    Uses PyRanges intersect (same semantics as merge_modkit_files), not per-row Python loops.
    """
    if df.empty:
        return df
    import pyranges as pr
    from robin.analysis.utilities.modkit_merge import panel_intersection_keys

    filter_ranges = _load_panel_pyranges(filter_bed_file, cache_dir)
    pr_df = df.rename(columns={"chrom": "Chromosome", "chromStart": "Start"}).copy()
    pr_df["End"] = pr_df["Start"] + 1
    gr = pr.PyRanges(pr_df[["Chromosome", "Start", "End"]])
    inter = gr.intersect(filter_ranges).df
    keys = panel_intersection_keys(inter)
    return df.merge(keys, on=["chrom", "chromStart"], how="inner")
