"""Shared, display-independent CNV evidence for the GUI and PDF reports.

Coordinates are BED-style half-open intervals. CNV input must be the masked
log2(ploidy / expected) track, never the sample-minus-control difference track.
These review candidates do not assign absolute copy number or zygosity.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from robin.classification_config import get_cnv_thresholds
from robin.reference_contigs import is_canonical_contig

MIN_USABLE_FRACTION = 0.4
MIN_SEGMENT_BINS = 3
SEGMENT_WINDOW = 5
MIN_SEGMENT_BP = 250_000


def interval_statistics(values, bin_width, start, end):
    """Overlap-weighted mean and median of unique finite overlapping bins.

    The median weights each overlapping bin once. Usable fraction is measured in
    base pairs, including missing/out-of-track portions of the requested interval.
    """
    if bin_width <= 0 or start < 0 or end <= start:
        raise ValueError("Expected positive bin width and nonempty half-open interval")
    arr = np.asarray(values, dtype=float).ravel()
    first, stop = start // bin_width, (end + bin_width - 1) // bin_width
    indexes = np.arange(first, min(stop, arr.size), dtype=int)
    vals = arr[indexes]
    weights = np.maximum(
        0,
        np.minimum((indexes + 1) * bin_width, end)
        - np.maximum(indexes * bin_width, start),
    )
    finite = np.isfinite(vals) & (weights > 0)
    usable = float(weights[finite].sum()) / (end - start)
    return {
        "mean_log2": (
            float(np.average(vals[finite], weights=weights[finite]))
            if finite.any()
            else None
        ),
        "median_log2": float(np.median(vals[finite])) if finite.any() else None,
        "bin_count": int(finite.sum()),
        "total_bins": int(stop - first),
        "usable_fraction": usable,
    }


def _state(stats, chromosome, sex):
    if stats["usable_fraction"] < MIN_USABLE_FRACTION or stats["median_log2"] is None:
        return "NO_DATA"
    gain, loss = get_cnv_thresholds(chromosome, sex)
    return (
        "GAIN"
        if stats["median_log2"] >= gain
        else "LOSS" if stats["median_log2"] <= loss else "NORMAL"
    )


def _overlap_labels(frame, chrom, start, end, label):
    if frame.empty:
        return []
    hits = frame[
        (frame.chrom == chrom) & (frame.start_pos < end) & (frame.end_pos > start)
    ]
    return sorted(
        {
            part.strip()
            for value in hits[label].astype(str)
            for part in value.split(",")
            if part.strip()
        }
    )


def detect_regional_segments(
    track, bin_width, sex, cytobands, excluded, panel, *, min_segment_bp=MIN_SEGMENT_BP
):
    """Find sustained bin states independently of cytoband boundaries.

    A centred five-bin median suppresses isolated spikes. Missing or excluded bins
    split runs and are never bridged. Retain >=3 bins spanning >=250 kb with >=60% raw bins beyond
    the state threshold and an agreeing raw median. Endpoints are bin-resolution
    candidate boundaries, not validated breakpoints.
    """
    if bin_width <= 0:
        raise ValueError("bin_width must be positive")
    if bin_width > 10_000_000:
        return []
    events = []
    for chrom, data in track.items():
        if not is_canonical_contig(chrom):
            continue
        arr = np.asarray(data, dtype=float).copy()
        for _, row in excluded[excluded.chrom == chrom].iterrows():
            arr[
                max(0, int(row.start_pos) // bin_width) : (
                    int(row.end_pos) + bin_width - 1
                )
                // bin_width
            ] = np.nan
        gain, loss = get_cnv_thresholds(chrom, sex)
        valid = np.isfinite(arr)
        # Smooth separately within finite runs: no information crosses a gap.
        smoothed = np.full(arr.shape, np.nan)
        edges = np.flatnonzero(np.diff(np.r_[False, valid, False]))
        for lo, hi in zip(edges[::2], edges[1::2]):
            smoothed[lo:hi] = (
                pd.Series(arr[lo:hi])
                .rolling(SEGMENT_WINDOW, center=True, min_periods=1)
                .median()
                .to_numpy()
            )
        states = np.where(
            valid & (smoothed >= gain), 1, np.where(valid & (smoothed <= loss), -1, 0)
        )
        edges = np.r_[0, np.flatnonzero(np.diff(states)) + 1, len(states)]
        for lo, hi in zip(edges[:-1], edges[1:]):
            code = states[lo] if hi > lo else 0
            if not code or hi - lo < MIN_SEGMENT_BINS:
                continue
            start, end = int(lo * bin_width), int(hi * bin_width)
            # Clip a terminal bin to the reference chromosome end where known.
            bands = cytobands[cytobands.chrom == chrom]
            if not bands.empty:
                end = min(end, int(bands.end_pos.max()))
            if end <= start or end - start < min_segment_bp:
                continue
            stats = interval_statistics(arr, bin_width, start, end)
            state = "GAIN" if code == 1 else "LOSS"
            affected = float(
                np.mean(arr[lo:hi] >= gain if code == 1 else arr[lo:hi] <= loss)
            )
            if affected < 0.6 or _state(stats, chrom, sex) != state:
                continue
            band_hits = cytobands[
                (cytobands.chrom == chrom)
                & (cytobands.start_pos < end)
                & (cytobands.end_pos > start)
            ].sort_values("start_pos")
            names = list(dict.fromkeys(band_hits["name"].astype(str)))
            events.append(
                {
                    "chrom": chrom.removeprefix("chr"),
                    "chromosome": chrom,
                    "region": f"{chrom} "
                    + ("–".join([names[0], names[-1]]) if names else "segment"),
                    "start_pos": start,
                    "end_pos": end,
                    "start_mb": start / 1e6,
                    "end_mb": end / 1e6,
                    "length_mb": (end - start) / 1e6,
                    "mean_cnv": stats["mean_log2"],
                    "state": state,
                    "panel_genes": _overlap_labels(panel, chrom, start, end, "gene"),
                    "affected_fraction": affected,
                    **stats,
                }
            )
    return sorted(
        events,
        key=lambda e: (
            int(e["chrom"]) if e["chrom"].isdigit() else 100,
            e["chrom"],
            e["start_pos"],
        ),
    )


def load_panel_gene_intervals(panel):
    """Resolve panel symbols to gene bodies; never duplicate a combined target.

    Unresolved symbols remain explicit target intervals labelled TARGET_REGION.
    They must not be presented as separately resolved gene-body measurements.
    """
    from robin import resources

    symbols = {
        p.strip() for x in panel.get("gene", []) for p in str(x).split(",") if p.strip()
    }
    records = []
    path = Path(resources.__file__).parent / "all_genes3.bed"
    if path.exists():
        frame = pd.read_csv(
            path, sep="\t", header=None, names=["chrom", "start_pos", "end_pos", "gene"]
        )
        for _, row in frame.iterrows():
            for symbol in str(row.gene).split(","):
                if symbol.strip() in symbols:
                    records.append(
                        {
                            "chrom": row.chrom,
                            "start_pos": int(row.start_pos),
                            "end_pos": int(row.end_pos),
                            "gene": symbol.strip(),
                            "scope": "GENE_BODY",
                        }
                    )
    resolved = {r["gene"] for r in records}
    for _, row in panel.iterrows():
        missing = [
            x.strip() for x in str(row.gene).split(",") if x.strip() not in resolved
        ]
        if missing:
            records.append(
                {
                    "chrom": row.chrom,
                    "start_pos": int(row.start_pos),
                    "end_pos": int(row.end_pos),
                    "gene": ",".join(missing),
                    "scope": "TARGET_REGION",
                }
            )
    if not records:
        return pd.DataFrame(columns=["chrom", "start_pos", "end_pos", "gene", "scope"])
    return (
        pd.DataFrame(records)
        .groupby(["chrom", "gene", "scope"], as_index=False)
        .agg({"start_pos": "min", "end_pos": "max"})
    )


def target_depth_comparisons(coverage):
    """Equal-weight target-row baselines, excluding the tested interval.

    Ratios are supporting depth measurements, not normalized CNV or copy counts.
    """
    required = {"chrom", "startpos", "endpos", "name", "coverage"}
    if not required.issubset(coverage.columns):
        return []
    frame = coverage.copy().reset_index(drop=True)
    frame["coverage"] = pd.to_numeric(frame.coverage, errors="coerce")
    frame = frame[np.isfinite(frame.coverage) & (frame.coverage >= 0)]
    results = []
    for i, row in frame.iterrows():
        # Exclude all rows for the same interval, including aliases/duplicates.
        others = frame[
            ~(
                (frame.chrom == row.chrom)
                & (frame.startpos == row.startpos)
                & (frame.endpos == row.endpos)
            )
        ]
        record = {
            "chromosome": row.chrom,
            "start_pos": int(row.startpos),
            "end_pos": int(row.endpos),
            "target": str(row["name"]),
            "coverage": float(row.coverage),
        }
        for label, group in [
            ("panel", others),
            ("chromosome", others[others.chrom == row.chrom]),
        ]:
            mean = float(group.coverage.mean()) if len(group) else None
            median = float(group.coverage.median()) if len(group) else None
            record.update(
                {
                    f"{label}_n": len(group),
                    f"{label}_mean": mean,
                    f"{label}_median": median,
                    f"{label}_mean_ratio": (
                        float(row.coverage) / mean
                        if mean is not None and mean > 0
                        else None
                    ),
                    f"{label}_median_ratio": (
                        float(row.coverage) / median
                        if median is not None and median > 0
                        else None
                    ),
                }
            )
        results.append(record)
    return results


def build_cnv_evidence(
    track,
    bin_width,
    sex,
    cytobands,
    excluded,
    panel,
    coverage,
    gene_intervals=None,
    *,
    min_segment_bp=MIN_SEGMENT_BP,
):
    """Single evidence builder consumed by live GUI and report/CSV generation."""
    intervals = (
        load_panel_gene_intervals(panel) if gene_intervals is None else gene_intervals
    )
    genes = []
    for _, row in intervals.iterrows():
        if not is_canonical_contig(row.chrom):
            continue
        start, end = int(row.start_pos), int(row.end_pos)
        if end <= start or start < 0:
            continue
        values = np.asarray(track.get(row.chrom, []), dtype=float).copy()
        for _, mask in excluded[excluded.chrom == row.chrom].iterrows():
            values[
                max(0, int(mask.start_pos) // bin_width) : (
                    int(mask.end_pos) + bin_width - 1
                )
                // bin_width
            ] = np.nan
        stats = interval_statistics(values, bin_width, start, end)
        state = _state(stats, row.chrom, sex) if bin_width <= 10_000_000 else "NO_DATA"
        gain, loss = get_cnv_thresholds(row.chrom, sex)
        genes.append(
            {
                "chromosome": row.chrom,
                "gene": row.gene,
                "scope": row.get("scope", "GENE_BODY"),
                "start_pos": start,
                "end_pos": end,
                "state": state,
                "gain_threshold": gain,
                "loss_threshold": loss,
                "resolution": "SINGLE_BIN" if stats["total_bins"] == 1 else "BIN_LEVEL",
                **stats,
            }
        )
    return {
        "segments": detect_regional_segments(
            track,
            bin_width,
            sex,
            cytobands,
            excluded,
            panel,
            min_segment_bp=min_segment_bp,
        ),
        "genes": genes,
        "target_depth": target_depth_comparisons(coverage),
    }


EVIDENCE_NOTES = (
    "Candidate calls use masked log2(ploidy / expected). Gene state uses the median "
    "of overlapping bins and requires at least 40% usable interval. Gain/loss cutoffs "
    "are chromosome/sex specific (default +0.30 / -0.30). Regional candidates use a "
    "five-bin median, at least three consecutive bins spanning 250 kb and 60% raw-bin agreement. "
    "Boundaries are approximate at the analysis bin width. Target depth compares "
    "equal-weight intervals, excluding the tested interval, and retains target-specific "
    "biases. No absolute copy number, zygosity or independent confirmation is assigned."
)

GENE_COLUMNS = [
    ("gene", "Gene"),
    ("chromosome", "Chr"),
    ("scope", "Interval type"),
    ("start_pos", "Start bp"),
    ("end_pos", "End bp"),
    ("state", "State"),
    ("mean_log2", "Mean log2"),
    ("median_log2", "Median log2"),
    ("bins", "Usable bins"),
    ("usable_pct", "Usable %"),
    ("loss_threshold", "Loss cutoff"),
    ("gain_threshold", "Gain cutoff"),
    ("resolution", "Resolution"),
]
DEPTH_COLUMNS = [
    ("target", "Target"),
    ("chromosome", "Chr"),
    ("coverage", "Depth x"),
    ("panel_n", "Panel n"),
    ("panel_mean", "Panel mean x"),
    ("panel_median", "Panel median x"),
    ("panel_mean_ratio", "Ratio to panel mean"),
    ("panel_median_ratio", "Ratio to panel median"),
    ("chromosome_n", "Chr n"),
    ("chromosome_mean", "Chr mean x"),
    ("chromosome_median", "Chr median x"),
    ("chromosome_mean_ratio", "Ratio to chr mean"),
    ("chromosome_median_ratio", "Ratio to chr median"),
]


def evidence_table_rows(evidence):
    """Shared human-readable numeric fields for GUI, PDF and CSV."""

    def number(value):
        return (
            round(float(value), 3) if value is not None and np.isfinite(value) else None
        )

    genes = [
        {
            "gene": r["gene"],
            "chromosome": r["chromosome"],
            "scope": r["scope"],
            "start_pos": r["start_pos"],
            "end_pos": r["end_pos"],
            "state": r["state"],
            "mean_log2": number(r["mean_log2"]),
            "median_log2": number(r["median_log2"]),
            "bins": r["bin_count"],
            "usable_pct": round(100 * r["usable_fraction"], 1),
            "loss_threshold": r["loss_threshold"],
            "gain_threshold": r["gain_threshold"],
            "resolution": r["resolution"],
        }
        for r in evidence["genes"]
    ]
    depth = [
        {
            "target": r["target"],
            "chromosome": r["chromosome"],
            "start_pos": r["start_pos"],
            "end_pos": r["end_pos"],
            **{
                k: number(r[k])
                for k in [
                    "coverage",
                    "panel_mean",
                    "panel_median",
                    "panel_mean_ratio",
                    "panel_median_ratio",
                    "chromosome_mean",
                    "chromosome_median",
                    "chromosome_mean_ratio",
                    "chromosome_median_ratio",
                ]
            },
            "panel_n": r["panel_n"],
            "chromosome_n": r["chromosome_n"],
        }
        for r in evidence["target_depth"]
    ]
    return genes, depth
