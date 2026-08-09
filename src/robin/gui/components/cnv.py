from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple
from pathlib import Path

import asyncio
import json
import natsort
import numpy as np
from functools import lru_cache
import logging
import pickle
import time
import importlib.resources as importlib_resources
import pandas as pd

try:
    from nicegui import ui
except ImportError:  # pragma: no cover
    ui = None

from robin.gui.theme import (
    styled_table,
    register_theme_sync_callback,
    get_user_dark_mode,
    client_timer,
    stop_timer,
)
from robin.analysis.cnv_classification import (
    CNVEvent,
    detect_cnv_events,
    format_cnv_events_card_lines,
    format_cnv_events_section_summary,
)
from robin.analysis.cnv_analysis import (
    compute_cnv_log2_from_ploidy,
    downsample_cnv_for_plot,
    prepare_cnv_calling_track,
    resolve_cnv_calling_track,
    resolve_cnv_plot_bin_width,
)
from robin.analysis.cnv_regional import (
    SIGNIFICANT_CNV_STATES,
    analyze_cytoband_cnv,
    build_regional_cnv_events,
    format_regional_event_table_row,
    is_reportable_chromosome,
    load_panel_gene_bed,
)
from robin.analysis.itd_work import load_gene_target_coverage
from robin.classification_config import get_cnv_thresholds

# Same chromosome set as reporting (plotting.py): chr0–chr22, chrX, chrY only
CNV_PLOT_CONTIGS = frozenset(
    ["chr" + str(i) for i in range(0, 23)] + ["chrX", "chrY"]
)

_CNV_PLOT_BIN_KEY_DEFAULT = "Data default"
_CNV_PLOT_BIN_OPTIONS = {
    _CNV_PLOT_BIN_KEY_DEFAULT: "Data default",
    "500 kb": "500 kb",
    "1 Mb": "1 Mb",
    "2 Mb": "2 Mb",
    "5 Mb": "5 Mb",
    "10 Mb": "10 Mb",
}
_CNV_PLOT_BIN_KEY_TO_BP = {
    _CNV_PLOT_BIN_KEY_DEFAULT: None,
    "500 kb": 500_000,
    "1 Mb": 1_000_000,
    "2 Mb": 2_000_000,
    "5 Mb": 5_000_000,
    "10 Mb": 10_000_000,
}
_CNV_PLOT_BIN_KEYS_ORDERED = list(_CNV_PLOT_BIN_KEY_TO_BP.keys())

# ECharts JS: y values to at most 1 decimal place, trailing ".0" dropped, so
# axis ticks and zoom-slider handles never show full float precision.
_CNV_Y_VALUE_FORMATTER_JS = (
    "(value) => { const n = Number(value); "
    "return Number.isFinite(n) ? String(Number(n.toFixed(1))) : value; }"
)
_CONFIGURED_GENES_SERIES_NAME = "configured_genes_highlight"
_CONFIGURED_GENES_LABELS_SERIES_NAME = "configured_genes_labels"
_CNV_GENE_COVERAGE_FILTER_ALL = "all"
_CNV_GENE_COVERAGE_FILTER_OUTLIERS = "outliers"
_CNV_GENE_COVERAGE_FILTERS = (
    _CNV_GENE_COVERAGE_FILTER_ALL,
    _CNV_GENE_COVERAGE_FILTER_OUTLIERS,
)

# CNV scatter plot container height (Tailwind arbitrary values on the echart).
_CNV_ABS_HEIGHT_CLASS = "h-[22.5rem]"
_CNV_ABS_HEIGHT_CLASS_TALL = "h-[45rem]"
_CNV_GENE_GAIN_COLOR = "#DC2626"
_CNV_GENE_LOSS_COLOR = "#2563EB"
_CNV_GENE_OUTLIER_SD = 3.0
# Soft safety only — axis auto-scales to highlighted genes within this envelope.
_CNV_LOLLIPOP_LOG_Y_SOFT_CAP = 20.0
_CNV_LOLLIPOP_LINEAR_Y_SOFT_CAP_FACTOR = 20.0


@lru_cache(maxsize=16)
def _load_cnv_gene_locations(gene_names: tuple[str, ...]) -> tuple[Dict[str, Any], ...]:
    """Resolve configured gene symbols to GRCh38 intervals in the packaged gene BED."""
    requested = {name.casefold(): name for name in gene_names if name.strip()}
    if not requested:
        return ()

    intervals: Dict[tuple[str, str], Dict[str, Any]] = {}
    try:
        resources = importlib_resources.files("robin.resources")
        unresolved = set(requested)
        # all_genes3 supplies gene-body coordinates. unique_genes includes a few
        # panel aliases/non-coding genes absent from that reference.
        for resource_name in ("all_genes3.bed", "unique_genes.bed"):
            if not unresolved:
                break
            resource = resources / resource_name
            with resource.open("r", encoding="utf-8") as handle:
                for line in handle:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 4:
                        continue
                    chrom = fields[0].strip()
                    try:
                        start_pos, end_pos = int(fields[1]), int(fields[2])
                    except ValueError:
                        continue
                    for raw_symbol in fields[3].split(","):
                        symbol_key = raw_symbol.strip().casefold()
                        if symbol_key not in unresolved:
                            continue
                        requested_name = requested[symbol_key]
                        key = (symbol_key, chrom)
                        existing = intervals.get(key)
                        if existing is None:
                            intervals[key] = {
                                "gene": requested_name,
                                "chrom": chrom,
                                "start_pos": start_pos,
                                "end_pos": end_pos,
                            }
                        else:
                            existing["start_pos"] = min(
                                existing["start_pos"], start_pos
                            )
                            existing["end_pos"] = max(existing["end_pos"], end_pos)
            unresolved -= {key[0] for key in intervals}
    except Exception:
        logging.warning("Could not load CNV gene locations", exc_info=True)
        return ()

    order = {name.casefold(): index for index, name in enumerate(gene_names)}
    return tuple(
        sorted(
            intervals.values(),
            key=lambda row: (order.get(str(row["gene"]).casefold(), 10**9), row["chrom"]),
        )
    )


def _configured_genes_on_chrom(
    gene_locations: Sequence[Dict[str, Any]],
    selected: str,
) -> List[Dict[str, Any]]:
    """Return configured gene intervals for one chromosome (or all when selected is All)."""
    if selected == "All":
        return list(gene_locations)
    return [row for row in gene_locations if str(row["chrom"]) == selected]


def _find_configured_gene_interval(
    gene_locations: Sequence[Dict[str, Any]],
    *,
    selected: str,
    gene_name: str,
) -> Optional[Dict[str, Any]]:
    """Look up a configured gene interval on the current chromosome."""
    want = str(gene_name).casefold()
    for row in _configured_genes_on_chrom(gene_locations, selected):
        if str(row["gene"]).casefold() == want:
            return row
    return None


def _configured_gene_mark_series(
    gene_locations: Sequence[Dict[str, Any]],
    *,
    selected: str,
    chrom_offsets: Dict[str, float],
    dark: bool,
) -> Dict[str, Any]:
    """Build a dedicated ECharts overlay for configured CNV gene locations."""
    lines: List[Dict[str, Any]] = []
    areas: List[List[Dict[str, Any]]] = []
    # Invisible anchors so ECharts keeps the series (empty scatter markLines can vanish).
    anchors: List[List[float]] = []
    label_color = "#d8b4fe" if dark else "#6b21a8"
    line_color = "#a855f7" if dark else "#7e22ce"
    for row in gene_locations:
        chrom = str(row["chrom"])
        if selected != "All" and chrom != selected:
            continue
        if selected == "All" and chrom not in chrom_offsets:
            continue
        start_pos = float(row["start_pos"])
        end_pos = float(row["end_pos"])
        midpoint = (start_pos + end_pos) / 2.0
        offset = chrom_offsets.get(chrom, 0.0) if selected == "All" else 0.0
        x_pos = midpoint + offset
        gene = str(row["gene"])
        anchors.append([x_pos, 0.0])
        lines.append(
            {
                "name": gene,
                "xAxis": x_pos,
                "lineStyle": {
                    "color": line_color,
                    "width": 1 if selected == "All" else 1.5,
                    "opacity": 0.75,
                },
                "label": {
                    "show": True,
                    "formatter": gene,
                    "position": "insideEndTop",
                    "rotate": 90,
                    "fontSize": 9 if selected == "All" else 11,
                    "color": label_color,
                },
            }
        )
        # Chromosome view also gets a labelled span so names remain readable.
        if selected != "All":
            areas.append(
                [
                    {
                        "name": gene,
                        "xAxis": start_pos,
                        "itemStyle": {
                            "color": (
                                "rgba(168, 85, 247, 0.12)"
                                if dark
                                else "rgba(126, 34, 206, 0.10)"
                            )
                        },
                        "label": {
                            "show": True,
                            "position": "insideTop",
                            "formatter": gene,
                            "color": label_color,
                            "fontSize": 11,
                        },
                    },
                    {"xAxis": end_pos},
                ]
            )
    return {
        "type": "scatter",
        "name": _CONFIGURED_GENES_SERIES_NAME,
        "data": anchors,
        "symbolSize": 0,
        "silent": True,
        "zlevel": 2,
        "markLine": {
            "symbol": "none",
            "animation": False,
            "data": lines,
        },
        "markArea": {
            "silent": True,
            "data": areas,
        },
    }


def _upsert_configured_gene_series(
    chart: Any,
    gene_locations: Sequence[Dict[str, Any]],
    *,
    selected: str,
    chrom_offsets: Dict[str, float],
    dark: bool,
) -> None:
    """Ensure the configured-gene overlay is present and last among chart series."""
    series = chart.options.get("series")
    if not isinstance(series, list):
        return
    chart.options["series"] = [
        s
        for s in series
        if s.get("name")
        not in (
            _CONFIGURED_GENES_SERIES_NAME,
            _CONFIGURED_GENES_LABELS_SERIES_NAME,
        )
    ]
    if not gene_locations:
        return
    chart.options["series"].append(
        _configured_gene_mark_series(
            gene_locations,
            selected=selected,
            chrom_offsets=chrom_offsets,
            dark=dark,
        )
    )


def _is_configured_gene_cnv_outlier(
    cnv_val: float,
    mean_cnv: float,
    std_cnv: float,
) -> bool:
    """True when gene CNV differs from the local expected average by >3 SD."""
    if not np.isfinite(cnv_val) or not np.isfinite(mean_cnv):
        return False
    if not np.isfinite(std_cnv) or std_cnv < 1e-6:
        return abs(cnv_val - mean_cnv) > 0.5
    return abs(cnv_val - mean_cnv) > _CNV_GENE_OUTLIER_SD * std_cnv


def _configured_gene_region_cnv(
    values: np.ndarray,
    *,
    start_pos: float,
    end_pos: float,
    bin_width: int,
    use_max_abs: bool,
) -> Optional[float]:
    """Peak (or max-abs) CNV across bins overlapping a gene interval."""
    arr = np.asarray(values, dtype=float)
    if arr.size == 0 or bin_width <= 0:
        return None
    start_bin = max(0, int(start_pos // bin_width))
    end_bin = min(arr.size - 1, int(end_pos // bin_width))
    if end_bin < start_bin:
        return None
    region = arr[start_bin : end_bin + 1]
    finite = region[np.isfinite(region)]
    if finite.size == 0:
        return None
    if use_max_abs:
        return float(finite[np.nanargmax(np.abs(finite))])
    return float(np.nanmax(finite))


def _chrom_track_stats(
    values: np.ndarray,
    *,
    use_log: bool,
) -> Tuple[float, float]:
    """Mean/std for outlier + direction calls on one chromosome track."""
    vals = np.asarray(values, dtype=float)
    finite = vals[np.isfinite(vals)]
    default_mean = 0.0 if use_log else 2.0
    if finite.size == 0:
        return default_mean, 1.0
    return float(np.mean(finite)), float(np.std(finite))


def _configured_gene_cnv_direction(
    cnv_val: float,
    *,
    use_log: bool,
    mean_cnv: float,
) -> str:
    """Classify gene CNV as gain or loss relative to the expected average."""
    if use_log:
        return "gain" if cnv_val >= 0.0 else "loss"
    return "gain" if cnv_val >= mean_cnv else "loss"


def _normalise_coverage_to_cnv_axis(
    coverage: float,
    *,
    mean_cov: float,
    scale_mean_cnv: float,
    use_log: bool,
) -> Optional[float]:
    """Map target depth onto the CNV y-axis for visual comparison.

    Linear/ploidy: ``scale_mean_cnv * (cov / mean_cov)`` so mean coverage sits
    at the CNV mean. Log2: ``log2(cov / mean_cov)`` so mean coverage sits at 0.
    """
    if not np.isfinite(coverage) or coverage <= 0:
        return None
    if not np.isfinite(mean_cov) or mean_cov <= 0:
        return None
    ratio = float(coverage) / float(mean_cov)
    if use_log:
        return float(np.log2(ratio))
    if not np.isfinite(scale_mean_cnv):
        return None
    return float(scale_mean_cnv) * ratio


def _build_configured_gene_coverage_points(
    gene_locations: Sequence[Dict[str, Any]],
    *,
    selected: str,
    chrom_offsets: Dict[str, float],
    abs_plot_map: Dict[str, np.ndarray],
    bin_width: int,
    coverage_by_gene: Dict[str, float],
    filter_mode: str,
    use_log: bool,
    scale_mean_cnv: float,
) -> Tuple[List[Dict[str, Any]], Optional[float]]:
    """Join configured genes to target coverage for abs-chart lollipops."""
    if not gene_locations or not coverage_by_gene:
        return [], None

    cov_lookup = {str(k).casefold(): float(v) for k, v in coverage_by_gene.items()}
    cov_vals = [v for v in cov_lookup.values() if np.isfinite(v) and v > 0]
    mean_cov = float(np.mean(cov_vals)) if cov_vals else None
    if mean_cov is None:
        return [], None

    baseline_y = 0.0 if use_log else float(scale_mean_cnv)
    chrom_stats: Dict[str, Tuple[float, float]] = {}
    points: List[Dict[str, Any]] = []
    for row in _configured_genes_on_chrom(gene_locations, selected):
        gene = str(row["gene"])
        chrom = str(row["chrom"])
        if selected == "All" and chrom not in chrom_offsets and chrom not in abs_plot_map:
            continue
        coverage = cov_lookup.get(gene.casefold())
        if coverage is None or not np.isfinite(coverage):
            continue
        track = abs_plot_map.get(chrom)
        if track is None:
            continue
        if chrom not in chrom_stats:
            chrom_stats[chrom] = _chrom_track_stats(track, use_log=use_log)
        mean_cnv, std_cnv = chrom_stats[chrom]
        cnv_val = _configured_gene_region_cnv(
            track,
            start_pos=float(row["start_pos"]),
            end_pos=float(row["end_pos"]),
            bin_width=int(bin_width),
            use_max_abs=use_log,
        )
        if cnv_val is None:
            continue
        if filter_mode == _CNV_GENE_COVERAGE_FILTER_OUTLIERS and not (
            _is_configured_gene_cnv_outlier(cnv_val, mean_cnv, std_cnv)
        ):
            continue
        y_norm = _normalise_coverage_to_cnv_axis(
            float(coverage),
            mean_cov=mean_cov,
            scale_mean_cnv=float(scale_mean_cnv),
            use_log=use_log,
        )
        if y_norm is None:
            continue
        start_pos = float(row["start_pos"])
        end_pos = float(row["end_pos"])
        midpoint = (start_pos + end_pos) / 2.0
        offset = chrom_offsets.get(chrom, 0.0) if selected == "All" else 0.0
        points.append(
            {
                "gene": gene,
                "chrom": chrom,
                "x": midpoint + offset,
                "coverage": float(coverage),
                "y": float(y_norm),
                "baseline_y": float(baseline_y),
                "cnv_val": float(cnv_val),
                "direction": _configured_gene_cnv_direction(
                    float(cnv_val),
                    use_log=use_log,
                    mean_cnv=mean_cnv,
                ),
            }
        )
    return points, mean_cov


def _stagger_lollipop_label_distances(
    points: Sequence[Dict[str, Any]],
    *,
    base_distance: int = 10,
    step: int = 14,
    x_proximity_frac: float = 0.018,
) -> List[int]:
    """Push nearby gene labels apart vertically (pixel distance from the head)."""
    if not points:
        return []
    xs = [float(p["x"]) for p in points]
    x_span = max(xs) - min(xs) if len(xs) > 1 else 1.0
    proximity = max(x_span * x_proximity_frac, 1.0)
    order = sorted(range(len(points)), key=lambda idx: xs[idx])
    distances = [base_distance] * len(points)
    occupied: List[Tuple[float, int]] = []
    for idx in order:
        x_pos = xs[idx]
        lane = 0
        while any(
            abs(x_pos - other_x) < proximity and lane == other_lane
            for other_x, other_lane in occupied
        ):
            lane += 1
            if lane > 8:
                break
        occupied.append((x_pos, lane))
        distances[idx] = base_distance + lane * step
    return distances


def _lollipop_baseline_y_range(
    *,
    use_log: bool,
    scale_mean_cnv: float,
) -> Tuple[float, float]:
    """Minimum CNV-friendly window before expanding for highlighted genes."""
    if use_log:
        return (-2.0, 2.0)
    return (0.0, max(4.0, float(scale_mean_cnv) * 1.75))


def _lollipop_view_y_range(
    points: Sequence[Dict[str, Any]],
    *,
    use_log: bool,
    scale_mean_cnv: float,
) -> Tuple[float, float]:
    """Auto-scale the Y window so highlighted genes and nearby labels stay in range."""
    base_lo, base_hi = _lollipop_baseline_y_range(
        use_log=use_log,
        scale_mean_cnv=scale_mean_cnv,
    )
    gene_ys = [
        float(p["y"])
        for p in points
        if p.get("y") is not None and np.isfinite(float(p["y"]))
    ]
    if not gene_ys:
        span = max(base_hi - base_lo, 1.0)
        return float(base_lo - 0.1 * span), float(base_hi + 0.1 * span)

    g_lo = float(min(gene_ys))
    g_hi = float(max(gene_ys))
    lo = min(float(base_lo), g_lo)
    hi = max(float(base_hi), g_hi)
    span = max(hi - lo, 1.0)
    # Room for labels placed just above gains / below losses.
    label_pad = max(0.35, 0.18 * span)
    lo -= label_pad
    hi += label_pad
    if not use_log:
        lo = max(0.0, lo)
    return float(lo), float(hi)


def _soft_cap_lollipop_display_y(
    y_norm: float,
    *,
    use_log: bool,
    scale_mean_cnv: float,
) -> Tuple[float, bool]:
    """Only clip pathological extremes; normal gene values stay unscaled."""
    if use_log:
        lo, hi = -_CNV_LOLLIPOP_LOG_Y_SOFT_CAP, _CNV_LOLLIPOP_LOG_Y_SOFT_CAP
    else:
        lo = 0.0
        hi = max(12.0, float(scale_mean_cnv) * _CNV_LOLLIPOP_LINEAR_Y_SOFT_CAP_FACTOR)
    capped = bool(y_norm < lo or y_norm > hi)
    return float(np.clip(y_norm, lo, hi)), capped


def _layout_lollipop_label_placements(
    points: Sequence[Dict[str, Any]],
    *,
    y_lo: float,
    y_hi: float,
    x_proximity_frac: float = 0.022,
    n_lanes: int = 4,
) -> List[Dict[str, Any]]:
    """Place each gene badge beside its marker: above gains, below losses."""
    if not points:
        return []
    span = max(float(y_hi) - float(y_lo), 1e-6)
    base_offset = max(0.22, 0.08 * span)
    lane_step = max(0.16, 0.045 * span)
    lanes = n_lanes if len(points) > 8 else max(2, min(n_lanes, 3))

    xs = [float(p["x"]) for p in points]
    x_span = max(xs) - min(xs) if len(xs) > 1 else 1.0
    proximity = max(x_span * x_proximity_frac, 1.0)
    order = sorted(range(len(points)), key=lambda idx: xs[idx])

    placements: List[Optional[Dict[str, Any]]] = [None] * len(points)
    occupied: List[Tuple[float, str, int]] = []  # x, side, lane

    for idx in order:
        point = points[idx]
        x_pos = float(point["x"])
        y_head = float(point.get("y_disp", point.get("y", 0.0)))
        baseline = float(point.get("baseline_y", 0.0))
        direction = str(point.get("direction") or "")
        # Prefer the side away from the CNV baseline / toward the stem tip.
        if direction == "loss" or y_head < baseline:
            side = "below"
        else:
            side = "above"

        lane = 0
        while any(
            abs(x_pos - other_x) < proximity
            and side == other_side
            and lane == other_lane
            for other_x, other_side, other_lane in occupied
        ):
            lane += 1
            if lane >= lanes:
                lane = lanes - 1
                break
        occupied.append((x_pos, side, lane))

        offset = base_offset + lane * lane_step
        if side == "above":
            label_y = min(y_head + offset, float(y_hi) - 0.03 * span)
            label_position = "top"
        else:
            label_y = max(y_head - offset, float(y_lo) + 0.03 * span)
            label_position = "bottom"
        placements[idx] = {
            "y": float(label_y),
            "side": side,
            "position": label_position,
        }

    return [
        p if p is not None else {"y": 0.0, "side": "above", "position": "top"}
        for p in placements
    ]


def _configured_gene_coverage_lollipop_series(
    points: Sequence[Dict[str, Any]],
    *,
    use_log: bool,
    dark: bool,
    scale_mean_cnv: float,
    view_y_lo: float,
    view_y_hi: float,
    label_font_size: Optional[int] = None,
) -> List[Dict[str, Any]]:
    """Build head + callout-label series on the shared CNV axis."""
    if not points:
        return []

    label_bg = "rgba(15, 23, 42, 0.94)" if dark else "rgba(255, 255, 255, 0.97)"
    head_border = "#0f172a" if dark else "#ffffff"
    if label_font_size is None:
        label_font = 12 if len(points) <= 20 else 11
    else:
        label_font = int(label_font_size)

    prepared: List[Dict[str, Any]] = []
    for point in points:
        y_raw = float(point["y"])
        y_disp, capped = _soft_cap_lollipop_display_y(
            y_raw,
            use_log=use_log,
            scale_mean_cnv=scale_mean_cnv,
        )
        prepared.append({**point, "y_disp": y_disp, "capped": capped})

    placements = _layout_lollipop_label_placements(
        prepared,
        y_lo=view_y_lo,
        y_hi=view_y_hi,
    )

    heads: List[Dict[str, Any]] = []
    labels: List[Dict[str, Any]] = []
    stems: List[Any] = []
    leaders: List[Any] = []

    for point, place in zip(prepared, placements):
        x_pos = float(point["x"])
        y_disp = float(point["y_disp"])
        baseline = float(point.get("baseline_y", 0.0 if use_log else 2.0))
        label_y = float(place["y"])
        gene = str(point["gene"])
        cov = float(point["coverage"])
        capped = bool(point.get("capped"))
        color = (
            _CNV_GENE_GAIN_COLOR
            if point.get("direction") == "gain"
            else _CNV_GENE_LOSS_COLOR
        )
        heads.append(
            {
                "name": gene,
                "value": [x_pos, y_disp],
                "coverage": cov,
                "capped": capped,
                "itemStyle": {
                    "color": color,
                    "borderColor": head_border,
                    "borderWidth": 2,
                    "shadowBlur": 5,
                    "shadowColor": "rgba(0, 0, 0, 0.3)",
                },
                "label": {"show": False},
            }
        )
        labels.append(
            {
                "name": gene,
                "value": [x_pos, label_y],
                "coverage": cov,
                "capped": capped,
                "symbol": "roundRect",
                "symbolSize": [2, 2],
                "itemStyle": {
                    "color": "transparent",
                    "borderWidth": 0,
                },
                "label": {
                    "show": True,
                    "formatter": gene,
                    "position": place["position"],
                    "distance": 2,
                    "rotate": 0,
                    "fontSize": label_font,
                    "fontWeight": "bold",
                    "color": color,
                    "backgroundColor": label_bg,
                    "borderColor": color,
                    "borderWidth": 1.5,
                    "borderRadius": 4,
                    "padding": [1, 2],
                    "align": "center",
                    "verticalAlign": "middle",
                },
            }
        )
        stems.append(
            [
                {
                    "coord": [x_pos, baseline],
                    "lineStyle": {
                        "color": color,
                        "type": "solid",
                        "width": 2,
                        "opacity": 0.9,
                    },
                },
                {"coord": [x_pos, y_disp]},
            ]
        )
        # Short leader from the marker head to the gene badge.
        leaders.append(
            [
                {
                    "coord": [x_pos, y_disp],
                    "lineStyle": {
                        "color": color,
                        "type": "dashed",
                        "width": 1.25,
                        "opacity": 0.85,
                    },
                },
                {"coord": [x_pos, label_y]},
            ]
        )

    tooltip = {
        "trigger": "item",
        ":formatter": (
            "(params) => { const d = params.data || {}; "
            "const cov = Number(d.coverage); "
            "const v = params.value; "
            "const y = Array.isArray(v) ? Number(v[1]) : Number(v); "
            "const covTxt = Number.isFinite(cov) ? cov.toFixed(1) + 'x' : 'n/a'; "
            "const yTxt = Number.isFinite(y) ? y.toFixed(2) : ''; "
            "const cap = d.capped ? ' (capped)' : ''; "
            "return params.name + ': ' + covTxt + "
            "(yTxt ? ' (norm ' + yTxt + ')' : '') + cap; }"
        ),
    }

    head_series = {
        "type": "scatter",
        "name": _CONFIGURED_GENES_SERIES_NAME,
        "yAxisIndex": 0,
        "symbolSize": 13,
        "zlevel": 4,
        "z": 12,
        "clip": False,
        "animation": False,
        "animationDuration": 0,
        "progressive": 0,
        "data": heads,
        "markLine": {
            "symbol": "none",
            "animation": False,
            "animationDuration": 0,
            "silent": True,
            "z": 11,
            "data": stems + leaders,
        },
        "tooltip": tooltip,
    }
    label_series = {
        "type": "scatter",
        "name": _CONFIGURED_GENES_LABELS_SERIES_NAME,
        "yAxisIndex": 0,
        "symbolSize": 1,
        "zlevel": 5,
        "z": 14,
        "clip": False,
        "animation": False,
        "animationDuration": 0,
        "progressive": 0,
        "data": labels,
        "tooltip": tooltip,
        "silent": False,
    }
    return [head_series, label_series]


def _apply_cnv_abs_y_window(chart: Any, y_lo: float, y_hi: float) -> None:
    """Pin the abs-chart Y axis and slider to a marker-aware window."""
    try:
        dz_list = chart.options.get("dataZoom")
        if isinstance(dz_list, list) and len(dz_list) > 1 and isinstance(dz_list[1], dict):
            dz_list[1]["startValue"] = float(y_lo)
            dz_list[1]["endValue"] = float(y_hi)
            dz_list[1]["filterMode"] = "none"
            dz_list[1].pop("start", None)
            dz_list[1].pop("end", None)
    except Exception:
        pass
    try:
        y_axes = chart.options.get("yAxis")
        if isinstance(y_axes, list) and y_axes and isinstance(y_axes[0], dict):
            y_axes[0]["min"] = float(y_lo)
            y_axes[0]["max"] = float(y_hi)
            y_axes[0]["scale"] = False
    except Exception:
        pass


def _upsert_configured_gene_coverage_lollipops(
    chart: Any,
    points: Sequence[Dict[str, Any]],
    *,
    use_log: bool,
    dark: bool,
    scale_mean_cnv: float,
    label_font_size: Optional[int] = None,
) -> None:
    """Replace the abs-chart gene series with coverage lollipops (or remove it)."""
    series = chart.options.get("series")
    if not isinstance(series, list):
        return
    chart.options["series"] = [
        s
        for s in series
        if s.get("name")
        not in (
            _CONFIGURED_GENES_SERIES_NAME,
            _CONFIGURED_GENES_LABELS_SERIES_NAME,
        )
    ]
    if not points:
        # No markers: restore a CNV-friendly baseline window.
        base_lo, base_hi = _lollipop_baseline_y_range(
            use_log=use_log,
            scale_mean_cnv=scale_mean_cnv,
        )
        span = max(base_hi - base_lo, 1.0)
        _apply_cnv_abs_y_window(chart, base_lo, base_hi + 0.15 * span)
        return
    view_y_lo, view_y_hi = _lollipop_view_y_range(
        points,
        use_log=use_log,
        scale_mean_cnv=scale_mean_cnv,
    )
    _apply_cnv_abs_y_window(chart, view_y_lo, view_y_hi)
    chart.options["series"].extend(
        _configured_gene_coverage_lollipop_series(
            points,
            use_log=use_log,
            dark=dark,
            scale_mean_cnv=scale_mean_cnv,
            view_y_lo=view_y_lo,
            view_y_hi=view_y_hi,
            label_font_size=label_font_size,
        )
    )


def _set_cnv_abs_coverage_axis(chart: Any, *, show: bool = False) -> None:
    """Keep the unused right-hand axis hidden (lollipops share the CNV axis)."""
    y_axes = chart.options.get("yAxis")
    if not isinstance(y_axes, list) or len(y_axes) < 2:
        return
    axis = y_axes[1]
    if not isinstance(axis, dict):
        return
    axis["show"] = bool(show)
    axis.pop("max", None)
    grid = chart.options.get("grid")
    if isinstance(grid, dict):
        grid["right"] = "8%" if show else "5%"


def _cnv_plot_bin_key_from_ui(value: Any) -> str:
    """Resolve NiceGUI select value (key, index, or event payload) to a plot-bin option key."""
    if value is None:
        return _CNV_PLOT_BIN_KEY_DEFAULT
    if isinstance(value, dict):
        inner = value.get("value", value.get("label"))
        if inner is not None and inner is not value:
            return _cnv_plot_bin_key_from_ui(inner)
    if isinstance(value, int) and not isinstance(value, bool):
        if 0 <= value < len(_CNV_PLOT_BIN_KEYS_ORDERED):
            return _CNV_PLOT_BIN_KEYS_ORDERED[value]
    if isinstance(value, str):
        if value in _CNV_PLOT_BIN_KEY_TO_BP:
            return value
        for key, label in _CNV_PLOT_BIN_OPTIONS.items():
            if value == label:
                return key
    return _CNV_PLOT_BIN_KEY_DEFAULT


def _cnv_plot_bin_bp_from_ui(value: Any) -> Optional[int]:
    return _CNV_PLOT_BIN_KEY_TO_BP.get(_cnv_plot_bin_key_from_ui(value))


def _cnv_plot_bin_key_from_bp(bp: Optional[int]) -> str:
    if bp is None:
        return _CNV_PLOT_BIN_KEY_DEFAULT
    for key, width in _CNV_PLOT_BIN_KEY_TO_BP.items():
        if width == bp:
            return key
    return _CNV_PLOT_BIN_KEY_DEFAULT


def _cnv_gene_coverage_filter_from_ui(value: Any) -> str:
    """Normalize Coverage genes switch value to a filter mode key.

    Switch semantics: True / on = outliers only; False / off = all configured genes.
    """
    if isinstance(value, bool):
        return (
            _CNV_GENE_COVERAGE_FILTER_OUTLIERS
            if value
            else _CNV_GENE_COVERAGE_FILTER_ALL
        )
    if isinstance(value, dict):
        inner = value.get("value", value.get("label"))
        if inner is not None and inner is not value:
            return _cnv_gene_coverage_filter_from_ui(inner)
    vlow = str(value or "").strip().lower()
    if vlow in (
        _CNV_GENE_COVERAGE_FILTER_ALL,
        "all genes",
        "all",
        "false",
        "0",
        "off",
    ):
        return _CNV_GENE_COVERAGE_FILTER_ALL
    if vlow in (
        _CNV_GENE_COVERAGE_FILTER_OUTLIERS,
        "outliers only",
        "outliers",
        "≠ average",
        "!= average",
        "vs average",
        "not average",
        "true",
        "1",
        "on",
    ):
        return _CNV_GENE_COVERAGE_FILTER_OUTLIERS
    return _CNV_GENE_COVERAGE_FILTER_OUTLIERS


def _cnv_contig_ok(contig: str) -> bool:
    """True if contig should be included in CNV plots (matches report behaviour)."""
    return contig in CNV_PLOT_CONTIGS


def _unwrap_cnv_track_map(raw: Any) -> Optional[Dict[str, np.ndarray]]:
    if not isinstance(raw, dict):
        return None
    if "cnv" in raw:
        inner = raw["cnv"]
        return inner if isinstance(inner, dict) else None
    return raw


def _cnv_sex_estimate_label(xy_val: Any) -> str:
    try:
        s = str(xy_val).strip().upper()
        if s in ("MALE", "XY"):
            return "Male"
        if s in ("FEMALE", "XX"):
            return "Female"
    except Exception:
        pass
    return "Unknown"


def _recompute_cnv_log2_state(state: Dict[str, Any]) -> None:
    """log2(ploidy / expected copy number) from the same track as the ploidy plot."""
    sample = _unwrap_cnv_track_map(state.get("cnv"))
    if not sample:
        state.pop("cnv_log2", None)
        return
    state["cnv_log2"] = resolve_cnv_calling_track(
        sample,
        _cnv_sex_estimate_label(state.get("xy")),
    )


def _cnv_genome_x_extent_bp(
    cnv_map: Dict[str, Any],
    binw_analysis: int,
    selected: str = "All",
) -> int:
    """Full genomic span in bp for the CNV scatter x-axis (independent of plot bin)."""
    if selected == "All":
        return sum(
            len(cnv_map[c]) * int(binw_analysis)
            for c in natsort.natsorted(cnv_map.keys())
            if _cnv_contig_ok(c)
        )
    chr_cnv = cnv_map.get(selected)
    return len(chr_cnv) * int(binw_analysis) if chr_cnv is not None else 0


def _cnv_set_genome_x_axis(chart: Any, x_axis_max: int) -> None:
    """Pin the x-axis to the full chromosome span (not downsampled data extent)."""
    xa = chart.options.setdefault("xAxis", {})
    if not isinstance(xa, dict):
        return
    xa.pop("max", None)  # drop dataMax so merges cannot keep a stale auto scale
    xa["type"] = "value"
    xa["min"] = 0
    xa["max"] = int(x_axis_max)
    xa["scale"] = True


def _cnv_reset_genome_x_data_zoom(chart: Any) -> None:
    """Reset the horizontal dataZoom to the full pinned x-axis span."""
    try:
        dz_list = chart.options.get("dataZoom")
        if not isinstance(dz_list, list) or not dz_list:
            return
        dz = dz_list[0]
        if not isinstance(dz, dict):
            return
        dz.pop("startValue", None)
        dz.pop("endValue", None)
        dz["start"] = 0
        dz["end"] = 100
    except Exception:
        pass


def _cnv_echarts_option_to_json(obj: Any) -> Any:
    """Convert ECharts option fragments to JSON-serializable form."""
    if isinstance(obj, dict):
        return {k: _cnv_echarts_option_to_json(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_cnv_echarts_option_to_json(x) for x in obj]
    if isinstance(obj, (np.floating, np.integer)):
        v = float(obj) if isinstance(obj, np.floating) else int(obj)
        return None if (isinstance(obj, np.floating) and np.isnan(obj)) else v
    if isinstance(obj, float) and np.isnan(obj):
        return None
    return obj


def _cnv_extract_js_options(obj: Any, js_values: Dict[str, str]) -> Any:
    """Swap NiceGUI ``:key`` JS options for tokens, recording the JS source.

    The raw ``setOption`` push below never runs NiceGUI's dynamic-property
    conversion, so these have to be re-inserted as unquoted JS afterwards.
    """
    if isinstance(obj, dict):
        converted: Dict[Any, Any] = {}
        for key, value in obj.items():
            if isinstance(key, str) and key.startswith(":") and isinstance(value, str):
                token = f"__ROBIN_JS_{len(js_values)}__"
                js_values[token] = value
                converted[key[1:]] = token
            else:
                converted[key] = _cnv_extract_js_options(value, js_values)
        return converted
    if isinstance(obj, list):
        return [_cnv_extract_js_options(item, js_values) for item in obj]
    return obj


def _cnv_echart_push_update(chart: Any) -> None:
    """Replace the full ECharts option (NiceGUI merges by default and leaves stale scatter data)."""
    js_values: Dict[str, str] = {}
    options_clean = _cnv_echarts_option_to_json(
        _cnv_extract_js_options(chart.options, js_values)
    )
    opts_json = json.dumps(options_clean)
    # ``:setOption`` evaluates each argument as JS, so tokens become real functions.
    for token, js_source in js_values.items():
        opts_json = opts_json.replace(json.dumps(token), js_source)
    try:
        # Do not call chart.update() here: NiceGUI's update_chart uses setOption merge
        # unless the series count changes, which leaves stale per-chromosome scatter data
        # when only the plot bin width changes.
        chart.run_chart_method(":setOption", opts_json, '{"notMerge": true}')
    except Exception:
        logging.debug("CNV chart notMerge setOption failed", exc_info=True)


def _cnv_sample_relative_stats(
    cnv_map: Dict[str, np.ndarray],
    abs_plot_map: Optional[Dict[str, np.ndarray]],
    *,
    use_log: bool,
) -> Tuple[float, float]:
    """Sample-relative baseline from autosomes only, shared by all Up/Down coloring."""
    source_map = abs_plot_map if use_log and abs_plot_map else cnv_map
    default_mean = 0.0 if use_log else 2.0
    autosome_vals: List[float] = []
    for chrom, arr in source_map.items():
        if not chrom.startswith("chr") or not chrom[3:].isdigit():
            continue
        vals = np.asarray(arr, dtype=float)
        autosome_vals.extend(float(v) for v in vals if np.isfinite(v))
    if not autosome_vals:
        return default_mean, 1.0
    return float(np.mean(autosome_vals)), float(np.std(autosome_vals))


def _cnv_split_points_by_zscore(
    pts: List[List[float]],
    mean_val: float,
    std_val: float,
) -> Tuple[List[List[float]], List[List[float]], List[List[float]]]:
    """Partition points into high / low / normal by z-score."""
    high: List[List[float]] = []
    low: List[List[float]] = []
    norm: List[List[float]] = []
    for xi, vi in pts:
        z = (vi - mean_val) / std_val if std_val > 0 else 0.0
        (high if z > 0.5 else low if z < -0.5 else norm).append([xi, vi])
    return high, low, norm


def _build_cnv_track_scatter_series(
    track_map: Dict[str, np.ndarray],
    *,
    selected: str,
    binw_analysis: int,
    plot_bin_width: int,
    chrom_palette: List[str],
    filter_finite: bool = False,
) -> List[Dict[str, Any]]:
    """Build ECharts scatter series for a per-chromosome CNV track."""
    series: List[Dict[str, Any]] = []
    if selected == "All":
        offset_bp = 0
        dj = 0
        for contig, cnv in natsort.natsorted(track_map.items()):
            if not _cnv_contig_ok(contig):
                continue
            x_local, vals = downsample_cnv_for_plot(
                np.asarray(cnv), binw_analysis, int(plot_bin_width)
            )
            x_global = offset_bp + x_local
            if filter_finite:
                pts = [
                    [float(x), float(v)]
                    for x, v in zip(x_global.tolist(), vals.tolist())
                    if np.isfinite(v)
                ]
            else:
                pts = list(zip(x_global.tolist(), [float(v) for v in vals]))
            offset_bp += len(cnv) * binw_analysis
            series.append(
                {
                    "type": "scatter",
                    "name": contig,
                    "symbolSize": 3,
                    "itemStyle": {"color": chrom_palette[dj % len(chrom_palette)]},
                    "data": pts,
                }
            )
            dj += 1
    else:
        cnv = track_map.get(selected)
        if cnv is not None:
            x_local, vals = downsample_cnv_for_plot(
                np.asarray(cnv), binw_analysis, int(plot_bin_width)
            )
            if filter_finite:
                pts = [
                    [float(x), float(v)]
                    for x, v in zip(x_local.tolist(), vals.tolist())
                    if np.isfinite(v)
                ]
            else:
                pts = list(zip(x_local.tolist(), [float(v) for v in vals]))
            series.append(
                {
                    "type": "scatter",
                    "name": selected,
                    "symbolSize": 3,
                    "itemStyle": {"color": chrom_palette[0]},
                    "data": pts,
                }
            )
    return series


def _cnv_load_binary_payload(
    sample_dir: Path,
    *,
    cnv_dict_npy_changed: bool,
    cnv_npy_changed: bool,
    cnv3_npy_changed: bool,
    data_array_reload: bool,
    xy_pkl_changed: bool,
) -> Dict[str, Any]:
    """Load CNV numpy/pickle data from disk (no UI). Safe for ``asyncio.to_thread``."""
    out: Dict[str, Any] = {}
    cnv_dict_npy = sample_dir / "CNV_dict.npy"
    cnv_npy = sample_dir / "CNV.npy"
    cnv3_npy = sample_dir / "CNV3.npy"
    data_array_npy = sample_dir / "cnv_data_array.npy"
    xy_pkl = sample_dir / "XYestimate.pkl"

    if cnv_dict_npy_changed and cnv_dict_npy.exists():
        out["cnv_dict"] = np.load(cnv_dict_npy, allow_pickle=True).item()

    if xy_pkl_changed and xy_pkl.exists():
        try:
            with xy_pkl.open("rb") as f:
                out["xy"] = pickle.load(f)
        except Exception:
            pass

    if cnv_npy_changed and cnv_npy.exists():
        try:
            out["cnv"] = np.load(cnv_npy, allow_pickle=True).item()
        except Exception:
            out["cnv"] = None

    if cnv3_npy_changed and cnv3_npy.exists():
        try:
            out["cnv3"] = np.load(cnv3_npy, allow_pickle=True).item()
        except Exception:
            out["cnv3"] = None

    if data_array_reload and data_array_npy.exists():
        try:
            out["bp_array"] = np.load(data_array_npy, allow_pickle=True)
        except Exception:
            pass

    return out


def _is_dark_mode() -> bool:
    """Return normalized per-user dark mode."""
    return get_user_dark_mode(default=False)


def _cnv_chromosome_scatter_palette(dark: bool) -> List[str]:
    """Distinct scatter colours per chromosome.

    ECharts defaults include very dark greys that disappear on midnight backgrounds;
    dark mode uses lighter, saturated hues (design.md §5–§6).
    """
    if dark:
        return [
            "#34d399",
            "#38bdf8",
            "#fbbf24",
            "#fb7185",
            "#a78bfa",
            "#2dd4bf",
            "#f472b6",
            "#facc15",
            "#4ade80",
            "#60a5fa",
            "#f97316",
            "#e879f9",
            "#c084fc",
            "#22d3ee",
            "#fde047",
            "#93c5fd",
            "#f87171",
            "#bef264",
            "#5eead4",
            "#fcd34d",
            "#7dd3fc",
            "#fda4af",
            "#86efac",
            "#d8b4fe",
            "#eab308",
            "#67e8f9",
        ]
    return [
        "#059669",
        "#0284c7",
        "#b45309",
        "#dc2626",
        "#7c3aed",
        "#0d9488",
        "#db2777",
        "#ca8a04",
        "#16a34a",
        "#2563eb",
        "#ea580c",
        "#c026d3",
        "#9333ea",
        "#0891b2",
        "#ca8a04",
        "#3b82f6",
        "#ef4444",
        "#65a30d",
        "#14b8a6",
        "#eab308",
        "#0ea5e9",
        "#ec4899",
        "#4ade80",
        "#8b5cf6",
        "#ca8a04",
        "#06b6d4",
    ]


def _cnv_value_mode_colors(dark: bool) -> Tuple[str, str, str]:
    """High / low / in-range scatter colours (semantic green / rose / slate)."""
    if dark:
        return ("#34d399", "#fb7185", "#94a3b8")
    return ("#007AFF", "#FF3B30", "#8E8E93")


def _cnv_echart_palette(dark: bool) -> Dict[str, str]:
    """Axis, title, and tooltip colours for CNV scatter plots (design.md §5, §6)."""
    if dark:
        return {
            "text": "#e2e8f0",
            "muted": "#94a3b8",
            "axis_line": "#64748b",
            "split": "rgba(148, 163, 184, 0.52)",
            "tooltip_bg": "rgba(15, 23, 42, 0.96)",
            "tooltip_border": "#334155",
        }
    return {
        "text": "#0f172a",
        "muted": "#475569",
        "axis_line": "#94a3b8",
        "split": "rgba(71, 85, 105, 0.33)",
        "tooltip_bg": "rgba(255, 255, 255, 0.98)",
        "tooltip_border": "#e2e8f0",
    }


def _apply_cnv_echart_chrome(echart: Any, dark: bool) -> None:
    """Apply light/dark readable chrome without touching series data."""
    p = _cnv_echart_palette(dark)
    try:
        o = echart.options
        if not isinstance(o, dict):
            return
        o["backgroundColor"] = "transparent"
        o["textStyle"] = {"color": p["text"]}
        title = o.get("title")
        if isinstance(title, dict):
            title["textStyle"] = {"color": p["text"], "fontSize": 14}
        o["tooltip"] = {
            **(o.get("tooltip") or {}),
            "backgroundColor": p["tooltip_bg"],
            "borderColor": p["tooltip_border"],
            "textStyle": {"color": p["text"]},
        }
        xa = o.get("xAxis")
        if isinstance(xa, dict):
            xa["axisLine"] = {"lineStyle": {"color": p["axis_line"]}}
            xa["axisLabel"] = {**(xa.get("axisLabel") or {}), "color": p["muted"]}
            xa["splitLine"] = {"lineStyle": {"color": p["split"]}}
        ya_list = o.get("yAxis")
        if isinstance(ya_list, list):
            for ya in ya_list:
                if not isinstance(ya, dict):
                    continue
                ya["axisLine"] = {"lineStyle": {"color": p["axis_line"]}}
                ya["axisLabel"] = {
                    **(ya.get("axisLabel") or {}),
                    "color": p["muted"],
                    ":formatter": _CNV_Y_VALUE_FORMATTER_JS,
                }
                ya["nameTextStyle"] = {"color": p["muted"]}
                ya["splitLine"] = {"lineStyle": {"color": p["split"]}}
        leg = o.get("legend")
        if isinstance(leg, dict):
            leg["textStyle"] = {"color": p["muted"]}
        try:
            dz_list = o.get("dataZoom")
            if isinstance(dz_list, list):
                for dz in dz_list:
                    if not isinstance(dz, dict):
                        continue
                    if dz.get("yAxisIndex") is not None:
                        # Handle labels otherwise show the raw data min/max.
                        dz[":labelFormatter"] = _CNV_Y_VALUE_FORMATTER_JS
                    dz["borderColor"] = p["axis_line"]
                    dz["fillerColor"] = (
                        "rgba(51, 65, 85, 0.35)"
                        if dark
                        else "rgba(148, 163, 184, 0.2)"
                    )
                    dz["handleStyle"] = {
                        "color": p["text"],
                        "borderColor": p["axis_line"],
                    }
                    dz["moveHandleStyle"] = {"color": p["muted"]}
                    dz["emphasis"] = {
                        "handleStyle": {"borderColor": p["text"]},
                    }
                    dbb = dz.get("dataBackground") or {}
                    if isinstance(dbb, dict):
                        dbb["lineStyle"] = {
                            **(dbb.get("lineStyle") or {}),
                            "color": p["muted"],
                        }
                        dbb["areaStyle"] = {
                            **(dbb.get("areaStyle") or {}),
                            "color": (
                                "rgba(148, 163, 184, 0.12)"
                                if dark
                                else "rgba(71, 85, 105, 0.08)"
                            ),
                        }
                        dz["dataBackground"] = dbb
        except Exception:
            pass
    except Exception:
        pass


def add_cnv_section(launcher: Any, sample_dir: Path) -> None:
    """Build the CNV UI section and attach refresh timers.

    Uses `launcher._cnv_state` for per-sample cache/state.
    Expects CNV.npy, CNV3.npy, CNV_dict.npy, XYestimate.pkl, cnv_data_array.npy in sample folder.
    The top chart can toggle between absolute ploidy and log2(ploidy / expected copy number).

    Controls now trigger immediate refresh instead of waiting for timer updates.
    """
    configured_gene_names: tuple[str, ...] = ()
    try:
        from robin.workflow_config import get_cnv_genes, load_workflow_toml

        workflow_config = None
        workflow_toml_path = getattr(launcher, "workflow_toml_path", None)
        if workflow_toml_path:
            workflow_config = load_workflow_toml(Path(workflow_toml_path))
        configured_gene_names = get_cnv_genes(workflow_config)
    except Exception:
        logging.warning("Could not resolve [cnv].genes for GUI plots", exc_info=True)
    configured_gene_locations = _load_cnv_gene_locations(configured_gene_names)
    if configured_gene_names:
        found = {str(row["gene"]).casefold() for row in configured_gene_locations}
        missing = [name for name in configured_gene_names if name.casefold() not in found]
        if missing:
            logging.warning(
                "No packaged GRCh38 location found for configured CNV genes: %s",
                ", ".join(missing),
            )

    with ui.element("div").classes("w-full min-w-0").props("id=analysis-detail-cnv"):
        with ui.element("div").classes("classification-insight-shell w-full min-w-0"):
            ui.label("Copy number (CNV)").classes(
                "classification-insight-heading text-headline-small"
            )
            with ui.element("div").classes(
                "classification-insight-card w-full min-w-0"
            ):
                with ui.column().classes("w-full min-w-0 gap-2 p-2 md:p-3"):
                    with ui.row().classes("items-center gap-2 min-w-0"):
                        ui.icon("person").classes("classification-insight-icon")
                        ui.label("Genome-wide profile").classes(
                            "classification-insight-model flex-1 min-w-0"
                        )
                    cnv_status = ui.label("Status: Awaiting Data").classes(
                        "classification-insight-result w-full"
                    )
                    cnv_xy = ui.label("Genetic sex: --").classes(
                        "classification-insight-meta w-full"
                    )
                    cnv_whole_chr_summary = ui.label("Whole chromosome: --").classes(
                        "classification-insight-meta w-full"
                    )
                    cnv_arm_summary = ui.label("Arm-level: --").classes(
                        "classification-insight-meta w-full"
                    )
                    with ui.row().classes(
                        "w-full justify-end gap-4 flex-wrap items-baseline"
                    ):
                        cnv_bin = ui.label("Bin width: --").classes(
                            "classification-insight-meta"
                        )
                        cnv_var = ui.label("Variance: --").classes(
                            "classification-insight-meta"
                        )
            with ui.row().classes(
                "w-full gap-3 items-center mb-2 flex-wrap mt-2"
            ):
                ui.label("Chromosome").classes("classification-insight-meta")
                cnv_chrom_select = ui.select(options={"All": "All"}, value="All").style(
                    "width: 160px"
                )
                ui.label("Gene").classes("classification-insight-meta ml-2")
                cnv_gene_select = ui.select(options={"All": "All"}, value="All").style(
                    "width: 200px"
                )
                _cnv_ui_state = launcher._cnv_state.setdefault(str(sample_dir), {})
                from robin.gui.plotting_preferences import (
                    resolve_cnv_gui_color_mode,
                    resolve_cnv_gui_gene_coverage_filter,
                    resolve_cnv_gui_show_breakpoints,
                    resolve_cnv_gui_y_scale,
                )

                _plot_prefs = getattr(launcher, "plotting_preferences", None)
                if "gene_coverage_filter" not in _cnv_ui_state:
                    _cnv_ui_state["gene_coverage_filter"] = (
                        resolve_cnv_gui_gene_coverage_filter(_plot_prefs)
                    )
                if "color_mode" not in _cnv_ui_state:
                    _cnv_ui_state["color_mode"] = resolve_cnv_gui_color_mode(
                        _plot_prefs
                    )
                if "y_scale" not in _cnv_ui_state:
                    _cnv_ui_state["y_scale"] = resolve_cnv_gui_y_scale(_plot_prefs)
                if "show_bp" not in _cnv_ui_state:
                    _cnv_ui_state["show_bp"] = resolve_cnv_gui_show_breakpoints(
                        _plot_prefs
                    )
                if "double_height" not in _cnv_ui_state:
                    _cnv_ui_state["double_height"] = False

                ui.label("Coverage genes").classes("classification-insight-meta ml-2")
                _gene_cov_filter = _cnv_ui_state.get("gene_coverage_filter")
                if _gene_cov_filter not in _CNV_GENE_COVERAGE_FILTERS:
                    _gene_cov_filter = resolve_cnv_gui_gene_coverage_filter(_plot_prefs)
                _cnv_ui_state["gene_coverage_filter"] = _gene_cov_filter
                with ui.row().classes("items-center gap-1"):
                    ui.label("All").classes("classification-insight-meta")
                    cnv_gene_cov_filter = (
                        ui.switch(
                            value=_gene_cov_filter
                            == _CNV_GENE_COVERAGE_FILTER_OUTLIERS,
                        )
                        .props("dense")
                        .tooltip(
                            "Left: all configured genes · Right: outliers only"
                        )
                    )
                    cnv_gene_cov_filter.value = (
                        _gene_cov_filter == _CNV_GENE_COVERAGE_FILTER_OUTLIERS
                    )
                    ui.label("Outliers").classes("classification-insight-meta")
                ui.label("Color by").classes("classification-insight-meta ml-2")
                _color_mode = _cnv_ui_state.get("color_mode", "chromosome")
                if _color_mode not in ("chromosome", "value"):
                    _color_mode = resolve_cnv_gui_color_mode(_plot_prefs)
                _cnv_ui_state["color_mode"] = _color_mode
                with ui.row().classes("items-center gap-1"):
                    ui.label("Chromosome").classes("classification-insight-meta")
                    cnv_color = (
                        ui.switch(value=_color_mode == "value")
                        .props("dense")
                        .tooltip(
                            "Left: colour by chromosome · Right: gain/loss (up/down)"
                        )
                    )
                    cnv_color.value = _color_mode == "value"
                    ui.label("Up/Down").classes("classification-insight-meta")
                ui.label("Y-axis").classes("classification-insight-meta ml-2")
                default_y_scale = resolve_cnv_gui_y_scale(_plot_prefs)
                _y_scale = _cnv_ui_state.get("y_scale", default_y_scale)
                if _y_scale not in ("linear", "log"):
                    _y_scale = default_y_scale
                _cnv_ui_state["y_scale"] = _y_scale
                with ui.row().classes("items-center gap-1"):
                    ui.label("Linear").classes("classification-insight-meta")
                    cnv_scale = (
                        ui.switch(value=_y_scale == "log")
                        .props("dense")
                        .tooltip(
                            "Left: linear ploidy · Right: log2(ploidy / expected)"
                        )
                    )
                    cnv_scale.value = _y_scale == "log"
                    ui.label("Log2").classes("classification-insight-meta")
                ui.label("Plot bin").classes("classification-insight-meta ml-2")
                cnv_plot_bin = ui.select(
                    options=_CNV_PLOT_BIN_OPTIONS,
                    value=_CNV_PLOT_BIN_KEY_DEFAULT,
                ).style("width: 120px")
                cnv_bp_label = ui.label("Breakpoints").classes(
                    "classification-insight-meta ml-2"
                ).style("display: none")
                with ui.row().classes("items-center gap-1").style(
                    "display: none"
                ) as cnv_bp_row:
                    ui.label("Hide").classes("classification-insight-meta")
                    cnv_bp = (
                        ui.switch(value=bool(_cnv_ui_state.get("show_bp", True)))
                        .props("dense")
                        .tooltip(
                            "Show candidate breakpoint markers on single-chromosome view"
                        )
                    )
                    cnv_bp.value = bool(_cnv_ui_state.get("show_bp", True))
                    ui.label("Show").classes("classification-insight-meta")
                ui.label("Height").classes("classification-insight-meta ml-2")
                _double_height = bool(_cnv_ui_state.get("double_height", False))
                _cnv_ui_state["double_height"] = _double_height
                with ui.row().classes("items-center gap-1"):
                    ui.label("Normal").classes("classification-insight-meta")
                    cnv_height = (
                        ui.switch(value=_double_height)
                        .props("dense")
                        .tooltip(
                            "Left: standard scatter height · Right: double height"
                        )
                    )
                    cnv_height.value = _double_height
                    ui.label("Tall").classes("classification-insight-meta")
            with ui.element("div").classes("w-full target-coverage-panel__plot-wrap"):
                cnv_abs = ui.echart(
                    {
                        "backgroundColor": "transparent",
                        "title": {"text": "CNV scatter plot", "left": "center", "top": 10},
                        "grid": {
                            "left": "5%",
                            "right": "5%",
                            "bottom": "10%",
                            "top": "20%",
                            "containLabel": True,
                        },
                        "tooltip": {"trigger": "axis"},
                        "xAxis": {"type": "value", "min": 0},
                        "yAxis": [
                            {"type": "value", "name": "Ploidy"},
                            {
                                "type": "value",
                                "name": "Coverage (x)",
                                "position": "right",
                                "show": False,
                                "min": 0,
                            },
                        ],
                        "dataZoom": [
                            {"type": "slider", "xAxisIndex": [0]},
                            {
                                "type": "slider",
                                "yAxisIndex": [0],
                                "right": 20,
                                "startValue": 0,
                                "endValue": 6,
                            },
                        ],
                        "series": [
                            {"type": "scatter", "name": "CNV", "symbolSize": 3, "data": []},
                            {
                                "type": "scatter",
                                "name": "centromeres_highlight",
                                "data": [],
                                "symbolSize": 3,
                                "markArea": {
                                    "itemStyle": {"color": "rgba(135, 206, 250, 0.4)"},
                                    "data": [],
                                },
                            },
                            {
                                "type": "scatter",
                                "name": "cytobands_highlight",
                                "data": [],
                                "symbolSize": 3,
                                "markArea": {
                                    "itemStyle": {"color": "rgba(200, 200, 200, 0.4)"},
                                    "data": [],
                                },
                                "markLine": {"symbol": "none", "data": []},
                            },
                        ],
                    }
                ).classes(
                    f"w-full {_CNV_ABS_HEIGHT_CLASS_TALL if _double_height else _CNV_ABS_HEIGHT_CLASS} cnv-genome-abs-chart"
                )
            with ui.element("div").classes("w-full target-coverage-panel__plot-wrap mt-2"):
                cnv_diff = ui.echart(
                    {
                        "backgroundColor": "transparent",
                        "title": {"text": "Difference plot", "left": "center", "top": 10},
                        "grid": {
                            "left": "5%",
                            "right": "5%",
                            "bottom": "10%",
                            "top": "20%",
                            "containLabel": True,
                        },
                        "tooltip": {"trigger": "axis"},
                        "xAxis": {"type": "value", "min": 0},
                        "yAxis": [
                            {"type": "value", "name": "Relative"},
                            {
                                "type": "value",
                                "name": "Coverage (x)",
                                "position": "right",
                                "show": False,
                                "min": 0,
                            },
                        ],
                        "dataZoom": [
                            {"type": "slider", "xAxisIndex": [0]},
                            {
                                "type": "slider",
                                "yAxisIndex": [0],
                                "right": 20,
                                "startValue": -4,
                                "endValue": 4,
                            },
                        ],
                        "series": [
                            {
                                "type": "scatter",
                                "name": "CNV Δ",
                                "symbolSize": 3,
                                "data": [],
                            },
                            {
                                "type": "scatter",
                                "name": "centromeres_highlight",
                                "data": [],
                                "symbolSize": 3,
                                "markArea": {
                                    "itemStyle": {"color": "rgba(135, 206, 250, 0.4)"},
                                    "data": [],
                                },
                            },
                            {
                                "type": "scatter",
                                "name": "cytobands_highlight",
                                "data": [],
                                "symbolSize": 3,
                                "markArea": {
                                    "itemStyle": {"color": "rgba(200, 200, 200, 0.4)"},
                                    "data": [],
                                },
                                "markLine": {"symbol": "none", "data": []},
                            },
                        ],
                    }
                ).classes("w-full h-[22.5rem] cnv-genome-diff-chart")
            genome_charts = (cnv_abs, cnv_diff)

            ui.separator().classes("mgmt-detail-separator")
            regional_cnv_label = ui.label("Regional CNV events").classes(
                "target-coverage-panel__meta-label mt-2 mb-1"
            )
            regional_cnv_summary = ui.label("No regional CNV events detected").classes(
                "classification-insight-meta mb-2"
            )
            regional_cnv_columns = [
                {"name": "chrom", "label": "Chr", "field": "chrom", "sortable": True},
                {"name": "region", "label": "Region", "field": "region", "sortable": True},
                {
                    "name": "start_mb",
                    "label": "Start (Mb)",
                    "field": "start_mb",
                    "sortable": True,
                    "align": "right",
                },
                {
                    "name": "end_mb",
                    "label": "End (Mb)",
                    "field": "end_mb",
                    "sortable": True,
                    "align": "right",
                },
                {
                    "name": "length_mb",
                    "label": "Length (Mb)",
                    "field": "length_mb",
                    "sortable": True,
                    "align": "right",
                },
                {
                    "name": "mean_cnv",
                    "label": "Mean CNV",
                    "field": "mean_cnv",
                    "sortable": True,
                    "align": "right",
                },
                {
                    "name": "state",
                    "label": "State",
                    "field": "state",
                    "sortable": True,
                    "align": "center",
                },
                {"name": "panel_genes", "label": "Panel genes", "field": "panel_genes"},
            ]
            _, regional_cnv_table = styled_table(
                columns=regional_cnv_columns, rows=[], pagination=20, class_size="table-xs"
            )
            try:
                regional_cnv_table.props('multi-sort rows-per-page-options="[10,20,50,0]"')
            except Exception:
                pass

            ui.separator().classes("mgmt-detail-separator")
            ui.label("Arm / whole-chromosome CNV events").classes(
                "target-coverage-panel__meta-label mt-2 mb-1"
            )
            cnv_events_summary = ui.label("No CNV events detected").classes(
                "classification-insight-meta mb-2"
            )
        
        # CNV Events Table
        cnv_events_columns = [
            {"name": "chromosome", "label": "Chr", "field": "chromosome", "sortable": True},
            {"name": "event_type", "label": "Event Type", "field": "event_type", "sortable": True},
            {"name": "arm", "label": "Arm", "field": "arm", "sortable": True},
            {"name": "start_mb", "label": "Start (Mb)", "field": "start_mb", "sortable": True, "align": "right"},
            {"name": "end_mb", "label": "End (Mb)", "field": "end_mb", "sortable": True, "align": "right"},
            {"name": "length_mb", "label": "Length (Mb)", "field": "length_mb", "sortable": True, "align": "right"},
            {"name": "mean_cnv_str", "label": "Mean CNV", "field": "mean_cnv_str", "sortable": True, "align": "right"},
            {"name": "confidence", "label": "Confidence", "field": "confidence", "sortable": True, "align": "center"},
            {"name": "proportion_affected", "label": "% Affected", "field": "proportion_affected", "sortable": True, "align": "right"},
            {"name": "genes_str", "label": "Genes", "field": "genes_str"},
        ]
        _, cnv_events_table = styled_table(
            columns=cnv_events_columns, rows=[], pagination=20, class_size="table-xs"
        )
        try:
            cnv_events_table.props('multi-sort rows-per-page-options="[10,20,50,0]"')
        except Exception:
            pass

    # Adaptive thinning helpers
    MAX_POINTS_PER_CHART = 10000

    def _get_visible_range(chart, series_list):
        try:
            # Compute overall data span
            min_x = None
            max_x = None
            for s in series_list:
                data = s.get("data") or []
                if not data:
                    continue
                sx = data[0][0] if isinstance(data[0], (list, tuple)) else None
                ex = data[-1][0] if isinstance(data[-1], (list, tuple)) else None
                if sx is None or ex is None:
                    # Fallback compute min/max
                    for p in data:
                        try:
                            x = float(p[0])
                        except Exception:
                            continue
                        min_x = x if min_x is None else min(min_x, x)
                        max_x = x if max_x is None else max(max_x, x)
                else:
                    vmin = min(float(sx), float(ex))
                    vmax = max(float(sx), float(ex))
                    min_x = vmin if min_x is None else min(min_x, vmin)
                    max_x = vmax if max_x is None else max(max_x, vmax)
            dz = None
            try:
                dzo = chart.options.get("dataZoom")
                if isinstance(dzo, list) and dzo:
                    dz = dzo[0]
            except Exception:
                dz = None
            if not dz:
                return (min_x, max_x)
            # Prefer explicit values
            sv = dz.get("startValue") if isinstance(dz, dict) else None
            ev = dz.get("endValue") if isinstance(dz, dict) else None
            if sv is not None or ev is not None:
                left = float(sv) if sv is not None else min_x
                right = float(ev) if ev is not None else max_x
                return (left, right)
            # Fallback to percentage range
            sp = dz.get("start") if isinstance(dz, dict) else None
            ep = dz.get("end") if isinstance(dz, dict) else None
            if (
                (sp is not None or ep is not None)
                and min_x is not None
                and max_x is not None
            ):
                width = (
                    max_x - min_x if max_x is not None and min_x is not None else None
                )
                if width and width > 0:
                    left = min_x + (float(sp or 0) / 100.0) * width
                    right = min_x + (float(ep or 100) / 100.0) * width
                    return (left, right)
            return (min_x, max_x)
        except Exception:
            return (None, None)

    def _evenly_sample(seq, k):
        try:
            n = len(seq)
            if k <= 0 or n <= k:
                return list(seq)
            if k == 1:
                return [seq[n // 2]]
            # Choose k indices evenly across [0, n-1]
            return [seq[int(round(i * (n - 1) / (k - 1)))] for i in range(k)]
        except Exception:
            return list(seq)[:k]

    def _thin_chart_series(chart, max_points: int = MAX_POINTS_PER_CHART) -> None:
        try:
            series = chart.options.get("series", [])
            # Identify data series to thin (exclude overlays)
            data_idx = []
            data_series = []
            for idx, s in enumerate(series):
                name = s.get("name", "")
                if s.get("type") == "scatter" and name not in (
                    "centromeres_highlight",
                    "cytobands_highlight",
                    _CONFIGURED_GENES_SERIES_NAME,
                    _CONFIGURED_GENES_LABELS_SERIES_NAME,
                ):
                    data = s.get("data") or []
                    if isinstance(data, list) and data:
                        data_idx.append(idx)
                        data_series.append(s)
            if not data_series:
                return
            x_range = _get_visible_range(chart, data_series)
            # Gather visible counts and data within range
            vis_data = []
            total = 0
            left, right = x_range
            for s in data_series:
                pts = s.get("data") or []
                if left is not None and right is not None:
                    sub = [
                        p
                        for p in pts
                        if isinstance(p, (list, tuple)) and left <= float(p[0]) <= right
                    ]
                else:
                    sub = list(pts)
                vis_data.append(sub)
                total += len(sub)
            if total <= max_points:
                return
            # Allocate budgets proportional to visible counts with a small floor
            budgets = []
            
            for sub in vis_data:
                share = int(max(1, round((len(sub) / total) * max_points)))
                budgets.append(share)
            # Normalize budgets to exactly max_points
            adj = sum(budgets) - max_points
            i = 0
            while adj != 0 and budgets:
                if adj > 0 and budgets[i] > 1:
                    budgets[i] -= 1
                    adj -= 1
                elif adj < 0:
                    budgets[i] += 1
                    adj += 1
                i = (i + 1) % len(budgets)
            # Apply sampling and replace data (preserve points outside range sparsely)
            for (idx, s), sub, k in zip(zip(data_idx, data_series), vis_data, budgets):
                original = s.get("data") or []
                # Keep outside-range points sparsely so context remains when zoomed out/in
                if left is not None and right is not None:
                    outside = [
                        p
                        for p in original
                        if isinstance(p, (list, tuple))
                        and not (left <= float(p[0]) <= right)
                    ]
                    outside_keep = _evenly_sample(
                        outside, max(0, k // 10)
                    )  # at most 10% of budget
                else:
                    outside_keep = []
                inside_keep = _evenly_sample(sub, max(1, k - len(outside_keep)))
                new_data = inside_keep + outside_keep
                chart.options["series"][idx]["data"] = new_data
        except Exception:
            pass

    @lru_cache(maxsize=1)
    def _load_centromere_regions() -> Dict[str, List[Tuple[int, int, str]]]:
        """Load centromere/satellite regions from packaged resources.
        Returns mapping: chrom -> list of (start_bp, end_bp, name).
        """
        regions: Dict[str, List[Tuple[int, int, str]]] = {}
        try:
            res_path = (
                importlib_resources.files("robin.resources") / "cenSatRegions.bed"
            )
            with res_path.open("r") as fh:
                for line in fh:
                    parts = line.strip().split("\t")
                    if len(parts) < 4:
                        continue
                    chrom, start, end, name = (
                        parts[0],
                        int(parts[1]),
                        int(parts[2]),
                        parts[3],
                    )
                    regions.setdefault(chrom, []).append((start, end, name))
        except Exception:
            pass
        return regions

    @lru_cache(maxsize=1)
    def _load_cytobands_df() -> pd.DataFrame:
        try:
            res_path = importlib_resources.files("robin.resources") / "cytoBand.txt"
            df = pd.read_csv(
                res_path,
                sep="\t",
                header=None,
                names=["chrom", "start_pos", "end_pos", "name", "stain"],
            )
            return df
        except Exception:
            return pd.DataFrame(columns=["chrom", "start_pos", "end_pos", "name", "stain"])

    @lru_cache(maxsize=1)
    def _load_centromere_bed_df() -> pd.DataFrame:
        try:
            res_path = importlib_resources.files("robin.resources") / "cenSatRegions.bed"
            return pd.read_csv(
                res_path,
                sep="\t",
                header=None,
                names=["chrom", "start_pos", "end_pos", "name"],
                usecols=[0, 1, 2, 3],
            )
        except Exception:
            return pd.DataFrame(columns=["chrom", "start_pos", "end_pos", "name"])

    _EMPTY_GENE_BED = pd.DataFrame(columns=["chrom", "start_pos", "end_pos", "gene"])

    def _analyze_cytoband_cnv(
        cnv_data: Dict[str, np.ndarray],
        chromosome: str,
        bin_width: int,
        sex_estimate: str,
    ) -> pd.DataFrame:
        """Run shared regional cytoband analysis (same logic as PDF reports)."""
        try:
            return analyze_cytoband_cnv(
                cnv_data,
                chromosome,
                {"bin_width": int(bin_width)},
                _load_cytobands_df(),
                _load_centromere_bed_df(),
                _EMPTY_GENE_BED,
                sex_estimate,
            )
        except Exception:
            return pd.DataFrame()

    def _load_gene_bed(sample_dir: Path = None) -> pd.DataFrame:
        """Load gene BED file based on the analysis panel used for the sample"""
        try:
            # Determine which panel to use
            panel = ""  # No default fallback
            if sample_dir:
                try:
                    master_csv_path = sample_dir / "master.csv"
                    if master_csv_path.exists():
                        import pandas as pd
                        df = pd.read_csv(master_csv_path)
                        if not df.empty and "analysis_panel" in df.columns:
                            panel_val = df.iloc[0]["analysis_panel"]
                            if panel_val and str(panel_val).strip() != "":
                                panel = str(panel_val).strip()
                except Exception:
                    pass
            
            # Map panel to BED filename
            bed_filename = None
            if not panel:
                # No panel found - return empty DataFrame
                return pd.DataFrame(columns=["chrom", "start_pos", "end_pos", "gene"])
            elif panel == "rCNS2":
                bed_filename = "rCNS2_panel_name_uniq.bed"
            elif panel == "AML":
                bed_filename = "AML_panel_name_uniq.bed"
            else:
                # Check for custom panel
                bed_filename = f"{panel}_panel_name_uniq.bed"
            
            # Try to load the panel-specific BED file
            try:
                res_path = importlib_resources.files("robin.resources") / bed_filename
                if res_path.exists():
                    return pd.read_csv(
                        res_path,
                        sep="\t",
                        header=None,
                        names=["chrom", "start_pos", "end_pos", "gene"],
                    )
            except Exception:
                pass
            
            # Fallback to unique_genes.bed if panel-specific file not found
            try:
                res_path = importlib_resources.files("robin.resources") / "unique_genes.bed"
                if res_path.exists():
                    return pd.read_csv(
                        res_path,
                        sep="\t",
                        header=None,
                        names=["chrom", "start_pos", "end_pos", "gene"],
                    )
            except Exception:
                pass
                
        except Exception:
            pass
            
        return pd.DataFrame(columns=["chrom", "start_pos", "end_pos", "gene"])

    def _sex_label(xy_val: Any) -> str:
        try:
            s = str(xy_val).strip().upper()
            if s in ("MALE", "XY"):
                return "Male"
            if s in ("FEMALE", "XX"):
                return "Female"
        except Exception:
            pass
        return "Unknown"

    def _get_cytoband_cnv_summary(
        cnv_data: Dict[str, np.ndarray],
        chromosome: str,
        bin_width: int,
        sex_estimate: str,
    ) -> str:
        try:
            df = _analyze_cytoband_cnv(cnv_data, chromosome, bin_width, sex_estimate)
            if df.empty:
                return "No regional CNV events detected"
            significant = df[df["cnv_state"].isin(SIGNIFICANT_CNV_STATES)]
            gains = significant[significant["cnv_state"].isin({"GAIN", "HIGH_GAIN"})]
            losses = significant[significant["cnv_state"].isin({"LOSS", "DEEP_LOSS"})]
            parts: List[str] = []
            if not gains.empty:
                parts.append(
                    "Gains: "
                    + ", ".join(
                        f"{r['name']} ({r['mean_cnv']:.2f})"
                        for _, r in gains.iterrows()
                    )
                )
            if not losses.empty:
                parts.append(
                    "Losses: "
                    + ", ".join(
                        f"{r['name']} ({r['mean_cnv']:.2f})"
                        for _, r in losses.iterrows()
                    )
                )
            return "\n".join(parts) if parts else "No regional CNV events detected"
        except Exception:
            return "No CNV data available"

    def _compute_all_cytoband_df(
        cnv_data: Dict[str, np.ndarray], bin_width: int, sex_estimate: str
    ) -> pd.DataFrame:
        try:
            frames: List[pd.DataFrame] = []
            for chrom in natsort.natsorted(cnv_data.keys()):
                if not is_reportable_chromosome(chrom):
                    continue
                df = _analyze_cytoband_cnv(cnv_data, chrom, bin_width, sex_estimate)
                if not df.empty:
                    frames.append(df)
            if frames:
                out = pd.concat(frames, ignore_index=True)
                if not out.empty:
                    def _rank(label: Any) -> int:
                        try:
                            s = str(label)
                            if s.startswith("chr"):
                                s = s[3:]
                            mapping = {"X": 23, "Y": 24, "M": 25}
                            return int(s) if s.isdigit() else mapping.get(s, 1000)
                        except Exception:
                            return 1000

                    out["_chrom_rank"] = out["chrom"].map(_rank)
                    out = out.sort_values(["_chrom_rank", "start_pos"]).drop(
                        columns=["_chrom_rank"]
                    )
                return out
            return pd.DataFrame()
        except Exception:
            return pd.DataFrame()

    def _build_regional_rows(
        cytoband_df: pd.DataFrame, panel_genes_df: pd.DataFrame
    ) -> List[Dict[str, Any]]:
        try:
            events = build_regional_cnv_events(cytoband_df, panel_genes_df)
            return [format_regional_event_table_row(event) for event in events]
        except Exception:
            return []

    def _update_cnv_events_analysis(state: Dict[str, Any]) -> None:
        """Update CNV events analysis using centralized classification rules."""
        try:
            cnv_map = state.get("cnv")
            if isinstance(cnv_map, dict) and "cnv" in cnv_map:
                cnv_map = cnv_map["cnv"]

            if not cnv_map:
                cnv_events_table.rows = []
                cnv_events_summary.set_text("No CNV data available")
                cnv_whole_chr_summary.set_text("Whole chromosome: --")
                cnv_arm_summary.set_text("Arm-level: --")
                return

            binw = state.get("cnv_dict", {}).get("bin_width", 1000000)
            sex_lbl = _sex_label(state.get("xy"))
            data, calling_binw = prepare_cnv_calling_track(
                cnv_map,
                int(binw),
                _cnv_sex_estimate_label(state.get("xy")),
            )

            if data and calling_binw:
                # Load cytobands and genes
                cyto_df = _load_cytobands_df()
                gene_df = _load_gene_bed(sample_dir)
                
                # Detect CNV events using centralized rules
                events = detect_cnv_events(
                    cnv_data=data,
                    bin_width=int(calling_binw),
                    sex_estimate=sex_lbl,
                    cytobands_df=cyto_df,
                    gene_df=gene_df
                )
                
                # Update events table
                events_rows = []
                for event in events:
                    event_dict = event.to_dict()
                    # Format proportion as percentage
                    event_dict["proportion_affected"] = f"{event.proportion_affected:.1%}"
                    events_rows.append(event_dict)
                
                cnv_events_table.rows = events_rows
                try:
                    cnv_events_table.update()
                except Exception:
                    pass
                
                # Update summaries (insight card + events section)
                whole_text, arm_text = format_cnv_events_card_lines(events)
                cnv_whole_chr_summary.set_text(whole_text)
                cnv_arm_summary.set_text(arm_text)
                cnv_events_summary.set_text(format_cnv_events_section_summary(events))
            else:
                cnv_events_table.rows = []
                cnv_whole_chr_summary.set_text("Whole chromosome: --")
                cnv_arm_summary.set_text("Arm-level: --")
                cnv_events_summary.set_text("CNV data not available")
        except Exception as e:
            logging.error(f"Error updating CNV events analysis: {e}")
            cnv_events_table.rows = []
            cnv_whole_chr_summary.set_text("Whole chromosome: --")
            cnv_arm_summary.set_text("Arm-level: --")
            cnv_events_summary.set_text("Error analyzing CNV events")

    def _render_cnv_from_state(state: Dict[str, Any]) -> None:
        try:
            _recompute_cnv_log2_state(state)
            cnv_map = state.get("cnv")
            cnv3_map = state.get("cnv3")
            cnv_log2_map = _unwrap_cnv_track_map(state.get("cnv_log2"))
            if isinstance(cnv_map, dict) and "cnv" in cnv_map:
                cnv_map = cnv_map["cnv"]
            if isinstance(cnv3_map, dict) and "cnv" in cnv3_map:
                cnv3_map = cnv3_map["cnv"]
            if not cnv_map:
                return
            dark_ui = _is_dark_mode()
            chrom_palette = _cnv_chromosome_scatter_palette(dark_ui)
            col_high, col_low, col_norm = _cnv_value_mode_colors(dark_ui)
            chrom_divider = _cnv_echart_palette(dark_ui)["muted"]
            binw_analysis = state.get("cnv_dict", {}).get("bin_width", 1_000_000)
            plot_bin_width = resolve_cnv_plot_bin_width(
                int(binw_analysis),
                state.get("plot_bin_width"),
            )
            # Keep plot bin width dropdown in sync with state
            try:
                pb = state.get("plot_bin_width")
                want_key = _cnv_plot_bin_key_from_bp(pb)
                if getattr(cnv_plot_bin, "value", None) != want_key:
                    cnv_plot_bin.value = want_key
                    cnv_plot_bin.update()
            except Exception:
                pass
            try:
                want_filter = str(
                    state.get(
                        "gene_coverage_filter",
                        _CNV_GENE_COVERAGE_FILTER_OUTLIERS,
                    )
                )
                if want_filter not in _CNV_GENE_COVERAGE_FILTERS:
                    want_filter = _CNV_GENE_COVERAGE_FILTER_OUTLIERS
                # Switch is boolean: True = outliers, False = all.
                want_outliers = want_filter == _CNV_GENE_COVERAGE_FILTER_OUTLIERS
                current = getattr(cnv_gene_cov_filter, "value", None)
                if current is not want_outliers:
                    cnv_gene_cov_filter.value = want_outliers
                    cnv_gene_cov_filter.update()
            except Exception:
                pass
            selected = state.get("selected_chrom", "All")
            use_log = state.get("y_scale", "linear") == "log"
            abs_plot_map = (
                cnv_log2_map
                if use_log and cnv_log2_map
                else cnv_map
            )
            raw_color_mode = state.get("color_mode", "chromosome")
            # normalize color mode to expected keys
            lval = str(raw_color_mode).strip().lower()
            if lval in ("chromosome", "chromosomes"):
                color_mode = "chromosome"
            elif lval in (
                "value",
                "up/down",
                "updown",
                "up_down",
                "updown ",
                "up down",
            ):
                color_mode = "value"
            else:
                color_mode = "chromosome"

            # Diff chart never uses the right axis; abs coverage axis is set with lollipops.
            if len(cnv_diff.options.get("yAxis") or []) > 1:
                cnv_diff.options["yAxis"][1]["show"] = False

            cnv_abs.options["yAxis"][0]["type"] = "value"
            cnv_abs.options["yAxis"][0].pop("logBase", None)
            # Clear any previous pinned Y window; marker overlay re-applies auto-scale.
            cnv_abs.options["yAxis"][0].pop("min", None)
            cnv_abs.options["yAxis"][0].pop("max", None)
            if use_log:
                cnv_abs.options["yAxis"][0]["name"] = "Log2 ratio (ploidy / expected)"
                cnv_abs.options["title"]["text"] = "CNV scatter plot"
                cnv_abs.options["title"]["top"] = 4
                cnv_abs.options["title"]["subtext"] = (
                    "log2(ploidy / expected copy number); 0 = normal"
                )
                cnv_abs.options["grid"]["top"] = "26%"
            else:
                cnv_abs.options["yAxis"][0]["name"] = "Ploidy"
                cnv_abs.options["title"]["text"] = "CNV scatter plot"
                cnv_abs.options["title"]["top"] = 10
                cnv_abs.options["title"].pop("subtext", None)
                cnv_abs.options["grid"]["top"] = "20%"
            # Temporary baseline until gene-marker auto-scale runs at the end of render.
            try:
                _apply_cnv_abs_y_window(
                    cnv_abs,
                    -2.0 if use_log else 0.0,
                    2.0 if use_log else 6.0,
                )
            except Exception:
                pass
            # X-axis is always in genomic base pairs; use actual genome/chromosome length
            # so the scale does not change when plot bin width changes (dataMax would shrink
            # with fewer downsampled points).
            x_axis_max = _cnv_genome_x_extent_bp(cnv_map, int(binw_analysis), selected)
            _cnv_set_genome_x_axis(cnv_abs, x_axis_max)
            _cnv_set_genome_x_axis(cnv_diff, x_axis_max)
            # Clear any previous zoom constraints when viewing All
            if selected == "All":
                _cnv_reset_genome_x_data_zoom(cnv_abs)
                _cnv_reset_genome_x_data_zoom(cnv_diff)
            logging.debug(
                f"CNV render: selected={selected}, y_scale={state.get('y_scale')}, color_mode={state.get('color_mode')}"
            )
            sample_rel_mean, sample_rel_std = _cnv_sample_relative_stats(
                cnv_map,
                abs_plot_map if isinstance(abs_plot_map, dict) else None,
                use_log=use_log,
            )
            try:
                coverage_by_gene = load_gene_target_coverage(sample_dir)
            except Exception:
                logging.debug(
                    "Could not load gene target coverage for CNV lollipops",
                    exc_info=True,
                )
                coverage_by_gene = {}

            def _apply_abs_gene_coverage_overlay() -> None:
                """Coverage lollipops on the abs chart only (position markers if no coverage)."""
                from robin.gui.plotting_preferences import (
                    resolve_cnv_gene_label_font_size,
                )

                filter_mode = str(
                    state.get(
                        "gene_coverage_filter",
                        _CNV_GENE_COVERAGE_FILTER_OUTLIERS,
                    )
                )
                if filter_mode not in _CNV_GENE_COVERAGE_FILTERS:
                    filter_mode = _CNV_GENE_COVERAGE_FILTER_OUTLIERS
                plot_map = (
                    abs_plot_map if isinstance(abs_plot_map, dict) else cnv_map
                )
                points, _mean_cov = _build_configured_gene_coverage_points(
                    configured_gene_locations,
                    selected=selected,
                    chrom_offsets=chrom_offsets,
                    abs_plot_map=plot_map,
                    bin_width=int(binw_analysis),
                    coverage_by_gene=coverage_by_gene,
                    filter_mode=filter_mode,
                    use_log=use_log,
                    scale_mean_cnv=float(sample_rel_mean),
                )
                _set_cnv_abs_coverage_axis(cnv_abs, show=False)
                label_font_size = resolve_cnv_gene_label_font_size(
                    getattr(launcher, "plotting_preferences", None)
                )
                if points:
                    _upsert_configured_gene_coverage_lollipops(
                        cnv_abs,
                        points,
                        use_log=use_log,
                        dark=dark_ui,
                        scale_mean_cnv=float(sample_rel_mean),
                        label_font_size=label_font_size,
                    )
                    return
                _upsert_configured_gene_coverage_lollipops(
                    cnv_abs,
                    (),
                    use_log=use_log,
                    dark=dark_ui,
                    scale_mean_cnv=float(sample_rel_mean),
                )
                # Fall back to position markers when target coverage is unavailable.
                if configured_gene_locations and not coverage_by_gene:
                    _upsert_configured_gene_series(
                        cnv_abs,
                        configured_gene_locations,
                        selected=selected,
                        chrom_offsets=chrom_offsets,
                        dark=dark_ui,
                    )

            # Absolute plot
            series_abs = []
            # Prepare chromosome partitions for labels/areas when viewing All
            chrom_bounds = []  # list of (name, start_bp, end_bp)
            chrom_offsets: Dict[str, float] = {}
            if selected == "All":
                # X-axis is in genomic base pairs; chromosome bounds use actual lengths.
                offset_bp = 0
                for contig, cnv in natsort.natsorted(cnv_map.items()):
                    if not _cnv_contig_ok(contig):
                        continue
                    plot_cnv = abs_plot_map.get(contig, cnv)
                    x_local, vals = downsample_cnv_for_plot(
                        np.asarray(plot_cnv), binw_analysis, int(plot_bin_width)
                    )
                    x_global = offset_bp + x_local
                    if use_log:
                        pts = [
                            [float(x), float(v)]
                            for x, v in zip(x_global.tolist(), vals.tolist())
                            if np.isfinite(v)
                        ]
                    else:
                        pts = list(zip(x_global.tolist(), [float(v) for v in vals]))
                    start_bp = offset_bp
                    end_bp = offset_bp + len(cnv) * binw_analysis
                    chrom_offsets[contig] = start_bp
                    chrom_bounds.append((contig, start_bp, end_bp))
                    offset_bp = end_bp
                    if color_mode == "chromosome":
                        ci = len(
                            [s for s in series_abs if s.get("type") == "scatter"]
                        )
                        series_abs.append(
                            {
                                "type": "scatter",
                                "name": contig,
                                "symbolSize": 3,
                                "itemStyle": {
                                    "color": chrom_palette[
                                        ci % len(chrom_palette)
                                    ]
                                },
                                "data": pts,
                            }
                        )
                    else:
                        high, low, norm = _cnv_split_points_by_zscore(
                            pts,
                            sample_rel_mean,
                            sample_rel_std,
                        )
                        if high:
                            series_abs.append(
                                {
                                    "type": "scatter",
                                    "name": f"High {contig}",
                                    "symbolSize": 4,
                                    "itemStyle": {"color": col_high},
                                    "data": high,
                                }
                            )
                        if low:
                            series_abs.append(
                                {
                                    "type": "scatter",
                                    "name": f"Low {contig}",
                                    "symbolSize": 4,
                                    "itemStyle": {"color": col_low},
                                    "data": low,
                                }
                            )
                        if norm:
                            series_abs.append(
                                {
                                    "type": "scatter",
                                    "name": f"Normal {contig}",
                                    "symbolSize": 2,
                                    "itemStyle": {"color": col_norm},
                                    "data": norm,
                                }
                            )
            else:
                plot_cnv = abs_plot_map.get(selected)
                if plot_cnv is None:
                    plot_cnv = cnv_map.get(selected)
                if plot_cnv is not None:
                    x_local, vals = downsample_cnv_for_plot(
                        np.asarray(plot_cnv), binw_analysis, int(plot_bin_width)
                    )
                    if use_log:
                        pts = [
                            [float(x), float(v)]
                            for x, v in zip(x_local.tolist(), vals.tolist())
                            if np.isfinite(v)
                        ]
                    else:
                        pts = list(zip(x_local.tolist(), [float(v) for v in vals]))
                    if color_mode == "chromosome":
                        series_abs.append(
                            {
                                "type": "scatter",
                                "name": selected,
                                "symbolSize": 3,
                                "itemStyle": {"color": chrom_palette[0]},
                                "data": pts,
                            }
                        )
                    else:
                        high, low, norm = _cnv_split_points_by_zscore(
                            pts,
                            sample_rel_mean,
                            sample_rel_std,
                        )
                        if high:
                            series_abs.append(
                                {
                                    "type": "scatter",
                                    "name": f"High {selected}",
                                    "symbolSize": 4,
                                    "itemStyle": {"color": col_high},
                                    "data": high,
                                }
                            )
                        if low:
                            series_abs.append(
                                {
                                    "type": "scatter",
                                    "name": f"Low {selected}",
                                    "symbolSize": 4,
                                    "itemStyle": {"color": col_low},
                                    "data": low,
                                }
                            )
                        if norm:
                            series_abs.append(
                                {
                                    "type": "scatter",
                                    "name": f"Normal {selected}",
                                    "symbolSize": 2,
                                    "itemStyle": {"color": col_norm},
                                    "data": norm,
                                }
                            )
            # Preserve highlight series (centromeres, cytobands) and replace data series only
            keep = [
                s
                for s in cnv_abs.options["series"]
                if s.get("name") in ("centromeres_highlight", "cytobands_highlight")
            ]
            cnv_abs.options["series"] = series_abs + keep
            _apply_abs_gene_coverage_overlay()
            if use_log and series_abs:
                series_abs[0]["markLine"] = {
                    "symbol": "none",
                    "data": [
                        {
                            "yAxis": 0,
                            "lineStyle": {
                                "type": "dashed",
                                "color": "#888888",
                                "width": 1,
                            },
                        }
                    ],
                }
            # Build background chromosome areas and vertical labels when showing All
            try:
                if selected == "All" and chrom_bounds:
                    # Alternating shaded bands per chromosome for readability
                    areas_data = []
                    lines_data = []
                    for contig, start_bp, end_bp in chrom_bounds:
                        areas_data.append(
                            [{"xAxis": float(start_bp)}, {"xAxis": float(end_bp)}]
                        )
                        # Label at the center of chromosome region
                        center_bp = (start_bp + end_bp) / 2
                        lines_data.append(
                            {
                                "xAxis": float(center_bp),
                                "lineStyle": {
                                    "type": "dashed",
                                    "color": chrom_divider,
                                },
                                "label": {
                                    "show": True,
                                    "formatter": contig,
                                    "color": chrom_divider,
                                },
                            }
                        )

                    # Helper to set overlays by series name regardless of index
                    def _apply_overlays(
                        chart, band_areas, band_lines, centro_areas_list
                    ):
                        try:
                            idx_cyto = next(
                                (
                                    i
                                    for i, s in enumerate(chart.options["series"])
                                    if s.get("name") == "cytobands_highlight"
                                ),
                                None,
                            )
                            idx_centro = next(
                                (
                                    i
                                    for i, s in enumerate(chart.options["series"])
                                    if s.get("name") == "centromeres_highlight"
                                ),
                                None,
                            )
                            if idx_cyto is not None:
                                # Remove background shading per request; keep dashed vertical lines only
                                chart.options["series"][idx_cyto]["markArea"][
                                    "data"
                                ] = []
                                chart.options["series"][idx_cyto].setdefault(
                                    "markLine", {}
                                )
                                chart.options["series"][idx_cyto]["markLine"][
                                    "data"
                                ] = band_lines
                                # Disable animation for marker lines so they appear instantly
                                chart.options["series"][idx_cyto]["markLine"][
                                    "animation"
                                ] = False
                            # Do not show centromeres in All view
                            if idx_centro is not None:
                                chart.options["series"][idx_centro]["markArea"][
                                    "data"
                                ] = []
                        except Exception:
                            pass

                    _apply_overlays(cnv_abs, areas_data, lines_data, [])
                    _apply_overlays(cnv_diff, areas_data, lines_data, [])
                else:
                    # Clear overlays when focusing on a single chromosome
                    def _clear_overlays(chart):
                        try:
                            for s in chart.options["series"]:
                                if s.get("name") in (
                                    "cytobands_highlight",
                                    "centromeres_highlight",
                                ):
                                    if "markArea" in s and "data" in s["markArea"]:
                                        s["markArea"]["data"] = []
                                    if "markLine" in s and "data" in s["markLine"]:
                                        s["markLine"]["data"] = []
                        except Exception:
                            pass

                    _clear_overlays(cnv_abs)
                    _clear_overlays(cnv_diff)
                    # In single-chromosome view, draw centromeres, cytobands (by state), genes, and breakpoint candidates
                    if selected != "All":
                        try:
                            centro = _load_centromere_regions()
                            idx_centro_abs = next(
                                (
                                    i
                                    for i, s in enumerate(cnv_abs.options["series"])
                                    if s.get("name") == "centromeres_highlight"
                                ),
                                None,
                            )
                            idx_centro_diff = next(
                                (
                                    i
                                    for i, s in enumerate(cnv_diff.options["series"])
                                    if s.get("name") == "centromeres_highlight"
                                ),
                                None,
                            )
                            areas = []
                            for s, e, _n in centro.get(selected, []):
                                areas.append([{"xAxis": float(s)}, {"xAxis": float(e)}])
                            if idx_centro_abs is not None:
                                cnv_abs.options["series"][idx_centro_abs]["markArea"][
                                    "data"
                                ] = areas
                            if idx_centro_diff is not None:
                                cnv_diff.options["series"][idx_centro_diff]["markArea"][
                                    "data"
                                ] = areas
                        except Exception:
                            pass
                        # Cytobands colored by CNV state (from CNV3 values)
                        try:
                            idx_cyto_abs = next(
                                (
                                    i
                                    for i, s in enumerate(cnv_abs.options["series"])
                                    if s.get("name") == "cytobands_highlight"
                                ),
                                None,
                            )
                            if (
                                idx_cyto_abs is not None
                                and cnv3_map
                                and selected in cnv3_map
                            ):
                                cyto_df = _load_cytobands_df()
                                bands = cyto_df[cyto_df["chrom"] == selected]
                                vals = (
                                    np.array(cnv3_map[selected])
                                    if isinstance(
                                        cnv3_map[selected], (list, np.ndarray)
                                    )
                                    else np.array([])
                                )
                                band_areas = []
                                
                                # Get CNV events for this chromosome to highlight significant events
                                events = []
                                try:
                                    sex_lbl = _sex_label(state.get("xy"))
                                    gene_df = _load_gene_bed(sample_dir)
                                    calling_map, calling_binw = prepare_cnv_calling_track(
                                        cnv_map,
                                        int(binw_analysis),
                                        _cnv_sex_estimate_label(state.get("xy")),
                                    )
                                    call_vals = calling_map.get(selected)
                                    if call_vals is not None:
                                        events = detect_cnv_events(
                                            cnv_data={selected: np.asarray(call_vals)},
                                            bin_width=int(calling_binw),
                                            sex_estimate=sex_lbl,
                                            cytobands_df=cyto_df,
                                            gene_df=gene_df,
                                        )
                                except Exception:
                                    pass
                                
                                # Create event lookup for highlighting
                                event_regions = {}
                                for event in events:
                                    key = f"{event.start_pos}-{event.end_pos}"
                                    event_regions[key] = event
                                
                                for _, row in bands.iterrows():
                                    s_bp, e_bp = int(row["start_pos"]), int(row["end_pos"])
                                    s_bin = max(0, s_bp // binw_analysis)
                                    e_bin = min(len(vals) - 1, max(0, e_bp // binw_analysis))
                                    if len(vals) > 0 and e_bin >= s_bin:
                                        mean_val = float(
                                            np.mean(vals[s_bin : e_bin + 1])
                                        )
                                    else:
                                        mean_val = 0.0
                                    
                                    # Check if this region has a significant CNV event
                                    region_key = f"{s_bp}-{e_bp}"
                                    event = event_regions.get(region_key)
                                    
                                    fill_neutral = (
                                        "rgba(255, 255, 255, 0.07)"
                                        if dark_ui
                                        else "rgba(0, 0, 0, 0.03)"
                                    )
                                    if event:
                                        # Highlight significant events with stronger colors
                                        if event.event_type in ("GAIN", "WHOLE_CHR_GAIN"):
                                            color = "rgba(52, 199, 89, 0.3)"  # gains
                                        elif event.event_type in ("LOSS", "WHOLE_CHR_LOSS"):
                                            color = "rgba(255, 45, 85, 0.3)"  # losses
                                        else:
                                            color = fill_neutral
                                    else:
                                        # Standard cytoband coloring
                                        if mean_val > 0.5:
                                            color = "rgba(52, 199, 89, 0.12)"
                                        elif mean_val < -0.5:
                                            color = "rgba(255, 45, 85, 0.12)"
                                        else:
                                            color = fill_neutral

                                    band_areas.append(
                                        [
                                            {
                                                "name": str(row["name"]),
                                                "xAxis": float(s_bp),
                                                "itemStyle": {"color": color},
                                                "label": {
                                                    "show": True,
                                                    "position": "insideTop",
                                                    "color": (
                                                        "#cbd5e1"
                                                        if dark_ui
                                                        else "#555"
                                                    ),
                                                    "fontSize": 11,
                                                },
                                            },
                                            {"xAxis": float(e_bp)},
                                        ]
                                    )
                                cnv_abs.options["series"][idx_cyto_abs]["markArea"][
                                    "data"
                                ] = band_areas
                        except Exception:
                            pass
                        # Configured genes of interest (TOML [cnv].genes) + gene selector
                        try:
                            gene_opts = {"All": "All"}
                            chrom_genes = _configured_genes_on_chrom(
                                configured_gene_locations, selected
                            )
                            for row in chrom_genes:
                                gene_opts[str(row["gene"])] = str(row["gene"])

                            # Fall back to panel gene BED when no TOML genes are configured.
                            if len(gene_opts) == 1:
                                gene_df = _load_gene_bed(sample_dir)
                                gchr = gene_df[gene_df["chrom"] == selected]
                                for _, gr in gchr.iterrows():
                                    gene_opts[str(gr["gene"])] = str(gr["gene"])
                                if series_abs:
                                    main = series_abs[0]
                                    mark = []
                                    for _, gr in gchr.iterrows():
                                        mark.append(
                                            [
                                                {
                                                    "name": str(gr["gene"]),
                                                    "xAxis": float(gr["start_pos"]),
                                                    "label": {
                                                        "position": "insideTop",
                                                        "color": (
                                                            "#e2e8f0"
                                                            if dark_ui
                                                            else "#000"
                                                        ),
                                                        "fontSize": 11,
                                                    },
                                                },
                                                {"xAxis": float(gr["end_pos"])},
                                            ]
                                        )
                                    main.setdefault("markArea", {"data": []})
                                    main["markArea"]["data"] = (
                                        main["markArea"]["data"] or []
                                    ) + mark
                                    series_abs[0] = main

                            try:
                                cnv_gene_select.set_options(gene_opts)
                                current_gene = launcher._cnv_state.setdefault(
                                    str(sample_dir), {}
                                ).get("selected_gene", "All")
                                if current_gene not in gene_opts:
                                    launcher._cnv_state[str(sample_dir)][
                                        "selected_gene"
                                    ] = "All"
                                    cnv_gene_select.value = "All"
                            except Exception:
                                pass
                        except Exception:
                            pass

                        # Keep gene overlays on top after cytoband / breakpoint overlays.
                        _apply_abs_gene_coverage_overlay()
                        _upsert_configured_gene_series(
                            cnv_diff,
                            configured_gene_locations,
                            selected=selected,
                            chrom_offsets=chrom_offsets,
                            dark=dark_ui,
                        )
                        
                        # Breakpoint candidates as dashed vertical lines
                        try:
                            idx_cyto_abs = next(
                                (
                                    i
                                    for i, s in enumerate(cnv_abs.options["series"])
                                    if s.get("name") == "cytobands_highlight"
                                ),
                                None,
                            )
                            if (
                                idx_cyto_abs is not None
                                and state.get("bp_array") is not None
                            ):
                                arr = state["bp_array"]
                                pos = [
                                    int(r["end_pos"]) for r in arr if r["name"] == selected
                                ]
                                lines = [
                                    {
                                        "xAxis": float(p),
                                        "lineStyle": {
                                            "type": "dashed",
                                            "color": "#E0162B",
                                        },
                                    }
                                    for p in pos
                                ]
                                cnv_abs.options["series"][idx_cyto_abs].setdefault(
                                    "markLine", {"symbol": "none", "data": []}
                                )
                                cnv_abs.options["series"][idx_cyto_abs]["markLine"][
                                    "data"
                                ] = lines
                        except Exception:
                            pass
                        # Re-assert abs coverage lollipops after breakpoint markLines mutate series.
                        _apply_abs_gene_coverage_overlay()
                # debug label removed
            except Exception:
                pass
            # Adaptive thinning based on current zoom and cap total points
            _thin_chart_series(cnv_abs, MAX_POINTS_PER_CHART)

            # Apply gene zoom before updating chart
            try:
                sel_gene = launcher._cnv_state.setdefault(
                    str(sample_dir), {}
                ).get("selected_gene", "All")

                gene_interval = None
                if sel_gene and sel_gene != "All" and selected != "All":
                    gene_interval = _find_configured_gene_interval(
                        configured_gene_locations,
                        selected=selected,
                        gene_name=sel_gene,
                    )
                    if gene_interval is None:
                        gene_df = _load_gene_bed(sample_dir)
                        gchr = gene_df[gene_df["chrom"] == selected]
                        row = gchr[gchr["gene"] == sel_gene]
                        if not row.empty:
                            gene_interval = {
                                "start_pos": int(row.iloc[0]["start_pos"]),
                                "end_pos": int(row.iloc[0]["end_pos"]),
                            }

                if gene_interval is not None:
                    s_bp = int(gene_interval["start_pos"])
                    e_bp = int(gene_interval["end_pos"])
                    pad = 10 * binw_analysis
                    zoom_start = max(0, s_bp - pad)
                    zoom_end = e_bp + pad
                    try:
                        cnv_abs.options["dataZoom"][0].update(
                            {
                                "startValue": zoom_start,
                                "endValue": zoom_end,
                                "start": None,
                                "end": None,
                            }
                        )
                    except Exception:
                        pass
                else:
                    # Reset zoom when "All" is selected
                    try:
                        if isinstance(cnv_abs.options.get("dataZoom"), list) and cnv_abs.options["dataZoom"]:
                            dz = cnv_abs.options["dataZoom"][0]
                            dz.pop("startValue", None)
                            dz.pop("endValue", None)
                            dz.update({"start": 0, "end": 100})
                    except Exception:
                        pass
            except Exception:
                pass

            # Ensure abs coverage lollipops remain on top after thinning / overlay mutations.
            _apply_abs_gene_coverage_overlay()

            _apply_cnv_echart_chrome(cnv_abs, _is_dark_mode())
            _cnv_echart_push_update(cnv_abs)
            # Difference plot (linear CNV3)
            if cnv3_map:
                series_diff = _build_cnv_track_scatter_series(
                    cnv3_map,
                    selected=selected,
                    binw_analysis=int(binw_analysis),
                    plot_bin_width=int(plot_bin_width),
                    chrom_palette=chrom_palette,
                )
                try:
                    base_series = cnv_diff.options["series"]
                    keep = [
                        s
                        for s in base_series
                        if s.get("name")
                        in ("centromeres_highlight", "cytobands_highlight")
                    ]
                except Exception:
                    keep = []
                cnv_diff.options["series"] = series_diff + keep
                _upsert_configured_gene_series(
                    cnv_diff,
                    configured_gene_locations,
                    selected=selected,
                    chrom_offsets=chrom_offsets,
                    dark=dark_ui,
                )
                _thin_chart_series(cnv_diff, MAX_POINTS_PER_CHART)
                _upsert_configured_gene_series(
                    cnv_diff,
                    configured_gene_locations,
                    selected=selected,
                    chrom_offsets=chrom_offsets,
                    dark=dark_ui,
                )
                _apply_cnv_echart_chrome(cnv_diff, _is_dark_mode())
                _cnv_echart_push_update(cnv_diff)
            else:
                _apply_cnv_echart_chrome(cnv_diff, _is_dark_mode())
                _cnv_echart_push_update(cnv_diff)

            # Gene zoom on difference chart (single-chromosome view)
            try:
                sel_gene = launcher._cnv_state.setdefault(
                    str(sample_dir), {}
                ).get("selected_gene", "All")
                for rel_chart in (cnv_diff,):
                    gene_interval = None
                    if sel_gene and sel_gene != "All" and selected != "All":
                        gene_interval = _find_configured_gene_interval(
                            configured_gene_locations,
                            selected=selected,
                            gene_name=sel_gene,
                        )
                        if gene_interval is None:
                            gene_df = _load_gene_bed(sample_dir)
                            gchr = gene_df[gene_df["chrom"] == selected]
                            row = gchr[gchr["gene"] == sel_gene]
                            if not row.empty:
                                gene_interval = {
                                    "start_pos": int(row.iloc[0]["start_pos"]),
                                    "end_pos": int(row.iloc[0]["end_pos"]),
                                }
                    if gene_interval is not None:
                        s_bp = int(gene_interval["start_pos"])
                        e_bp = int(gene_interval["end_pos"])
                        pad = 10 * binw_analysis
                        zoom_start = max(0, s_bp - pad)
                        zoom_end = e_bp + pad
                        try:
                            rel_chart.options["dataZoom"][0].update(
                                {
                                    "startValue": zoom_start,
                                    "endValue": zoom_end,
                                    "start": None,
                                    "end": None,
                                }
                            )
                            _cnv_echart_push_update(rel_chart)
                        except Exception:
                            pass
                    else:
                        try:
                            if (
                                isinstance(rel_chart.options.get("dataZoom"), list)
                                and rel_chart.options["dataZoom"]
                            ):
                                dz = rel_chart.options["dataZoom"][0]
                                dz.pop("startValue", None)
                                dz.pop("endValue", None)
                                dz.update({"start": 0, "end": 100})
                        except Exception:
                            pass
            except Exception:
                pass

            # Regional CNV events table (same logic as PDF reports)
            try:
                selected = state.get("selected_chrom", "All")
                binw = state.get("cnv_dict", {}).get("bin_width", 1_000_000)
                sex_lbl = _sex_label(state.get("xy"))
                if isinstance(cnv3_map, dict):
                    data = cnv3_map
                    source = "cnv3"
                else:
                    data = cnv_map if isinstance(cnv_map, dict) else None
                    source = "cnv"
                if data and binw:
                    cache_key = (
                        f"{source}:{state.get(source+'_m')}:{int(binw)}:{sex_lbl}"
                    )
                    if state.get("cyto_cache_key") != cache_key:
                        panel_name, panel_genes_df = load_panel_gene_bed(str(sample_dir))
                        df_all = _compute_all_cytoband_df(data, int(binw), sex_lbl)
                        state["cyto_df_all"] = df_all
                        state["panel_name"] = panel_name
                        state["panel_genes_df"] = panel_genes_df
                        state["cyto_cache_key"] = cache_key
                    df_all = state.get("cyto_df_all")
                    panel_genes_df = state.get("panel_genes_df")
                    if not isinstance(panel_genes_df, pd.DataFrame):
                        _, panel_genes_df = load_panel_gene_bed(str(sample_dir))
                    panel_name = state.get("panel_name")
                    regional_title = "Regional CNV events"
                    if panel_name:
                        regional_title += f" ({panel_name} panel genes)"
                    regional_cnv_label.set_text(regional_title)

                    if isinstance(df_all, pd.DataFrame) and not df_all.empty:
                        if selected and selected != "All":
                            df_show = df_all[df_all["chrom"] == selected]
                        else:
                            df_show = df_all
                        regional_rows = _build_regional_rows(df_show, panel_genes_df)
                        regional_cnv_table.rows = regional_rows
                        try:
                            regional_cnv_table.update()
                        except Exception:
                            pass
                        panel_gene_count = sum(
                            1 for row in regional_rows if row.get("panel_genes") != "—"
                        )
                        if selected and selected != "All":
                            regional_cnv_summary.set_text(
                                _get_cytoband_cnv_summary(
                                    data, selected, int(binw), sex_lbl
                                )
                            )
                        elif regional_rows:
                            summary = (
                                f"Detected {len(regional_rows)} regional CNV events"
                            )
                            if panel_gene_count:
                                summary += f"; {panel_gene_count} with panel genes"
                            regional_cnv_summary.set_text(summary)
                        else:
                            regional_cnv_summary.set_text(
                                "No regional CNV events detected"
                            )
                    else:
                        regional_cnv_table.rows = []
                        try:
                            regional_cnv_table.update()
                        except Exception:
                            pass
                        regional_cnv_summary.set_text(
                            "No regional CNV events detected"
                        )
                else:
                    regional_cnv_summary.set_text("CNV data not available")
            except Exception:
                pass
        except Exception:
            pass

    def _prepare_cnv_refresh() -> Optional[Dict[str, Any]]:
        """Main-thread debounce, UI sync, mtime checks. Returns None if no work needed."""
        key = str(sample_dir)
        state = launcher._cnv_state.get(key, {})

        current_time = time.time()
        last_refresh = state.get("_last_refresh", 0)
        force_color_refresh = state.get("_force_color_refresh", False)
        force_gene_refresh = state.get("_force_gene_refresh", False)
        force_chrom_refresh = state.get("_force_chrom_refresh", False)
        force_gene_cov_filter_refresh = state.get(
            "_force_gene_cov_filter_refresh", False
        )
        force_ui_refresh = state.get("_force_ui_refresh", False)
        force_refresh = (
            force_color_refresh
            or force_gene_refresh
            or force_chrom_refresh
            or force_gene_cov_filter_refresh
            or force_ui_refresh
        )
        # Never debounce away an explicit UI control click.
        if current_time - last_refresh < 0.1 and not force_refresh:
            return None
        state["_last_refresh"] = current_time

        ui_changed = False
        try:
            ui_sel = getattr(cnv_chrom_select, "value", None)
            if ui_sel and ui_sel != state.get("selected_chrom"):
                state["selected_chrom"] = ui_sel
                ui_changed = True
            ui_scale = getattr(cnv_scale, "value", None)
            if ui_scale is not None:
                # Switch: True = log2 ratio, False = linear ploidy.
                if isinstance(ui_scale, bool):
                    want_scale = "log" if ui_scale else "linear"
                else:
                    scale_key = str(ui_scale).strip().lower()
                    if scale_key in (
                        "log",
                        "log2",
                        "log2 ratio",
                        "log2 ratio (ploidy / expected)",
                        "true",
                        "1",
                    ):
                        want_scale = "log"
                    elif scale_key in ("linear", "ploidy", "false", "0"):
                        want_scale = "linear"
                    else:
                        want_scale = (
                            "log" if "log" in scale_key else state.get("y_scale", "linear")
                        )
                if want_scale != state.get("y_scale"):
                    state["y_scale"] = want_scale
                    ui_changed = True
            ui_bp = getattr(cnv_bp, "value", None)
            if ui_bp is not None:
                # Switch: True = show breakpoints, False = hide.
                if isinstance(ui_bp, bool):
                    desired = ui_bp
                else:
                    desired = str(ui_bp).strip().lower() in (
                        "show",
                        "true",
                        "1",
                        "on",
                    )
                if desired != state.get("show_bp", True):
                    state["show_bp"] = desired
                    ui_changed = True
            ui_color = getattr(cnv_color, "value", None)
            current_color_mode = state.get("color_mode", "chromosome")
            if ui_color is not None:
                # Switch: True = up/down, False = chromosome.
                if isinstance(ui_color, bool):
                    want_color = "value" if ui_color else "chromosome"
                else:
                    vlow = str(ui_color).strip().lower()
                    if vlow in ("value", "up/down", "updown", "true", "1"):
                        want_color = "value"
                    else:
                        want_color = "chromosome"
                if want_color != current_color_mode:
                    state["color_mode"] = want_color
                    ui_changed = True
            ui_plot_bin = getattr(cnv_plot_bin, "value", None)
            if ui_plot_bin is not None:
                want_bin = _cnv_plot_bin_bp_from_ui(ui_plot_bin)
                if want_bin != state.get("plot_bin_width"):
                    state["plot_bin_width"] = want_bin
                    ui_changed = True
            ui_gene_cov = getattr(cnv_gene_cov_filter, "value", None)
            if ui_gene_cov is not None:
                # Coerce legacy string values that may still be on the widget.
                if isinstance(ui_gene_cov, bool):
                    want_filter = (
                        _CNV_GENE_COVERAGE_FILTER_OUTLIERS
                        if ui_gene_cov
                        else _CNV_GENE_COVERAGE_FILTER_ALL
                    )
                else:
                    want_filter = _cnv_gene_coverage_filter_from_ui(ui_gene_cov)
                if want_filter != state.get(
                    "gene_coverage_filter", _CNV_GENE_COVERAGE_FILTER_OUTLIERS
                ):
                    state["gene_coverage_filter"] = want_filter
                    ui_changed = True
                # Keep the widget strictly boolean after any legacy string value.
                want_outliers = want_filter == _CNV_GENE_COVERAGE_FILTER_OUTLIERS
                if ui_gene_cov is not want_outliers:
                    try:
                        cnv_gene_cov_filter.value = want_outliers
                    except Exception:
                        pass
        except Exception:
            pass

        is_fresh_visit = "last_visit_time" not in state
        if is_fresh_visit:
            state["last_visit_time"] = time.time()

        cnv_npy = sample_dir / "CNV.npy"
        cnv3_npy = sample_dir / "CNV3.npy"
        cnv_dict_npy = sample_dir / "CNV_dict.npy"
        data_array_npy = sample_dir / "cnv_data_array.npy"
        xy_pkl = sample_dir / "XYestimate.pkl"

        cnv_npy_mtime = cnv_npy.stat().st_mtime if cnv_npy.exists() else 0
        cnv3_npy_mtime = cnv3_npy.stat().st_mtime if cnv3_npy.exists() else 0
        cnv_dict_npy_mtime = cnv_dict_npy.stat().st_mtime if cnv_dict_npy.exists() else 0
        data_array_npy_mtime = data_array_npy.stat().st_mtime if data_array_npy.exists() else 0
        xy_pkl_mtime = xy_pkl.stat().st_mtime if xy_pkl.exists() else 0

        prev_cnv_npy_mtime = state.get("cnv_m", 0)
        prev_cnv3_npy_mtime = state.get("cnv3_m", 0)
        prev_cnv_dict_npy_mtime = state.get("dict_m", 0)
        prev_data_array_npy_mtime = state.get("bp_array_mtime", 0)
        prev_xy_pkl_mtime = state.get("xy_m", 0)

        cnv_npy_changed = prev_cnv_npy_mtime != cnv_npy_mtime
        cnv3_npy_changed = prev_cnv3_npy_mtime != cnv3_npy_mtime
        cnv_dict_npy_changed = prev_cnv_dict_npy_mtime != cnv_dict_npy_mtime
        data_array_npy_changed = prev_data_array_npy_mtime != data_array_npy_mtime
        xy_pkl_changed = prev_xy_pkl_mtime != xy_pkl_mtime

        files_changed = (
            cnv_npy_changed
            or cnv3_npy_changed
            or cnv_dict_npy_changed
            or data_array_npy_changed
            or xy_pkl_changed
        )

        needs_update = (
            is_fresh_visit
            or files_changed
            or ui_changed
            or force_refresh
            or not state.get("_rendered_once")
        )

        if not needs_update:
            logging.debug("[CNV] ⏭ Skipping CNV update - no changes detected")
            launcher._cnv_state[key] = state
            return None

        reasons = []
        if is_fresh_visit:
            reasons.append("fresh_visit")
        if cnv_npy_changed:
            reasons.append("CNV.npy")
        if cnv3_npy_changed:
            reasons.append("CNV3.npy")
        if cnv_dict_npy_changed:
            reasons.append("CNV_dict.npy")
        if data_array_npy_changed:
            reasons.append("cnv_data_array.npy")
        if xy_pkl_changed:
            reasons.append("XYestimate.pkl")
        if ui_changed:
            reasons.append("ui_changed")
        if force_color_refresh:
            reasons.append("force_color_refresh")
        if force_gene_refresh:
            reasons.append("force_gene_refresh")
        if force_chrom_refresh:
            reasons.append("force_chrom_refresh")
        if force_gene_cov_filter_refresh:
            reasons.append("force_gene_cov_filter_refresh")
        if force_ui_refresh:
            reasons.append("force_ui_refresh")
        if not state.get("_rendered_once"):
            reasons.append("first_render")
        logging.debug(f"[CNV] Update needed. Reasons: {', '.join(reasons)}")

        data_array_reload = data_array_npy_changed or is_fresh_visit
        need_load = (
            cnv_dict_npy_changed
            or cnv_npy_changed
            or cnv3_npy_changed
            or data_array_reload
            or xy_pkl_changed
        )

        return {
            "state": state,
            "key": key,
            "ui_changed": ui_changed,
            "is_fresh_visit": is_fresh_visit,
            "cnv_npy": cnv_npy,
            "cnv3_npy": cnv3_npy,
            "cnv_dict_npy": cnv_dict_npy,
            "data_array_npy": data_array_npy,
            "xy_pkl": xy_pkl,
            "cnv_npy_mtime": cnv_npy_mtime,
            "cnv3_npy_mtime": cnv3_npy_mtime,
            "cnv_dict_npy_mtime": cnv_dict_npy_mtime,
            "data_array_npy_mtime": data_array_npy_mtime,
            "xy_pkl_mtime": xy_pkl_mtime,
            "cnv_dict_npy_changed": cnv_dict_npy_changed,
            "cnv_npy_changed": cnv_npy_changed,
            "cnv3_npy_changed": cnv3_npy_changed,
            "data_array_npy_changed": data_array_npy_changed,
            "xy_pkl_changed": xy_pkl_changed,
            "data_array_reload": data_array_reload,
            "need_load": need_load,
        }

    def _apply_breakpoint_marklines(
        chart,
        selected: str,
        state: Dict[str, Any],
        breakpoint_lines: List[int],
        *,
        reference_at_zero: bool = False,
    ) -> None:
        """Attach breakpoint x-lines (and optional y=0 reference) to the main data series."""
        try:
            current_series = [
                s
                for s in chart.options["series"]
                if not s.get("name", "").startswith("Breakpoint")
            ]
            if selected != "All" and state.get("show_bp", True) and breakpoint_lines:
                if current_series:
                    mark_line_data: List[Dict[str, Any]] = []
                    if reference_at_zero:
                        mark_line_data.append(
                            {
                                "yAxis": 0,
                                "lineStyle": {
                                    "type": "dashed",
                                    "color": "#888888",
                                    "width": 1,
                                },
                            }
                        )
                    for bp_pos in breakpoint_lines:
                        mark_line_data.append(
                            {
                                "xAxis": bp_pos,
                                "lineStyle": {
                                    "type": "dashed",
                                    "color": "#ff6b6b",
                                    "width": 3,
                                },
                            }
                        )
                    current_series[0]["markLine"] = {
                        "data": mark_line_data,
                        "symbol": "none",
                        "lineStyle": {"type": "dashed", "color": "#ff6b6b", "width": 3},
                    }
            elif reference_at_zero and current_series:
                current_series[0]["markLine"] = {
                    "symbol": "none",
                    "data": [
                        {
                            "yAxis": 0,
                            "lineStyle": {
                                "type": "dashed",
                                "color": "#888888",
                                "width": 1,
                            },
                        }
                    ],
                }
            else:
                if current_series:
                    current_series[0].pop("markLine", None)
            chart.options["series"] = current_series
            _upsert_configured_gene_series(
                chart,
                configured_gene_locations,
                selected=selected,
                chrom_offsets={},
                dark=_is_dark_mode(),
            )
            _apply_cnv_echart_chrome(chart, _is_dark_mode())
            _cnv_echart_push_update(chart)
        except Exception:
            pass

    def _apply_cnv_refresh_after_load(
        plan: Dict[str, Any], payload: Dict[str, Any]
    ) -> None:
        """Merge disk payload into state and update charts (main thread only)."""
        p = plan
        state = p["state"]
        key = p["key"]
        ui_changed = p["ui_changed"]
        is_fresh_visit = p["is_fresh_visit"]
        cnv_dict_npy = p["cnv_dict_npy"]
        cnv_npy = p["cnv_npy"]
        cnv3_npy = p["cnv3_npy"]
        data_array_npy = p["data_array_npy"]
        xy_pkl = p["xy_pkl"]
        cnv_npy_mtime = p["cnv_npy_mtime"]
        cnv3_npy_mtime = p["cnv3_npy_mtime"]
        cnv_dict_npy_mtime = p["cnv_dict_npy_mtime"]
        data_array_npy_mtime = p["data_array_npy_mtime"]
        xy_pkl_mtime = p["xy_pkl_mtime"]
        cnv_dict_npy_changed = p["cnv_dict_npy_changed"]
        cnv_npy_changed = p["cnv_npy_changed"]
        cnv3_npy_changed = p["cnv3_npy_changed"]
        data_array_npy_changed = p["data_array_npy_changed"]
        data_array_reload = p["data_array_reload"]
        xy_pkl_changed = p["xy_pkl_changed"]

        changed = ("cnv" in payload or "cnv3" in payload or "xy" in payload)

        if cnv_dict_npy.exists():
            m = cnv_dict_npy_mtime
            if cnv_dict_npy_changed and "cnv_dict" in payload:
                state["cnv_dict"] = payload["cnv_dict"]
                cnv_bin.set_text(
                    f"Bin width: {state['cnv_dict'].get('bin_width', '--'):,}"
                )
                cnv_var.set_text(
                    f"Variance: {state['cnv_dict'].get('variance','--'):.3f}"
                    if isinstance(state["cnv_dict"].get("variance"), (int, float))
                    else "Variance: --"
                )
                state["dict_m"] = m
        if xy_pkl.exists():
            m = xy_pkl_mtime
            if xy_pkl_changed:
                try:
                    if "xy" in payload:
                        xy = payload["xy"]
                        cnv_xy.set_text(f"Genetic sex: {xy}")
                        state["xy"] = xy
                except Exception:
                    pass
                state["xy_m"] = m

        if "cnv" in payload:
            state["cnv"] = payload["cnv"]
            state["cnv_m"] = cnv_npy_mtime
        if "cnv3" in payload:
            state["cnv3"] = payload["cnv3"]
            state["cnv3_m"] = cnv3_npy_mtime
        _recompute_cnv_log2_state(state)

        if state.get("cnv"):
            if changed:
                cnv_status.set_text("Status: CNV data loaded")
            if changed or not state.get("_chrom_opts_set"):
                try:
                    cnv_map = state["cnv"].get("cnv", state["cnv"])
                    chrom_opts = {"All": "All"}
                    for contig in natsort.natsorted(cnv_map.keys()):
                        if _cnv_contig_ok(contig):
                            chrom_opts[contig] = contig
                    cnv_chrom_select.set_options(chrom_opts)
                    chrom_opts = {"All": "All"}
                    for contig in natsort.natsorted(cnv_map.keys()):
                        if _cnv_contig_ok(contig):
                            chrom_opts[contig] = contig
                    cnv_chrom_select.set_options(chrom_opts)
                    sel = launcher._cnv_state.setdefault(str(sample_dir), {}).get(
                        "selected_chrom", "All"
                    )
                    if sel not in chrom_opts:
                        sel = "All"
                    cnv_chrom_select.value = sel
                    try:
                        cnv_chrom_select.update()
                    except Exception:
                        pass
                    state["_chrom_opts_set"] = True
                except Exception:
                    pass
            force_color_refresh = state.get("_force_color_refresh", False)
            force_gene_refresh = state.get("_force_gene_refresh", False)
            force_chrom_refresh = state.get("_force_chrom_refresh", False)
            force_gene_cov_filter_refresh = state.get(
                "_force_gene_cov_filter_refresh", False
            )
            force_ui_refresh = state.get("_force_ui_refresh", False)
            if (
                changed
                or ui_changed
                or not state.get("_rendered_once")
                or force_color_refresh
                or force_gene_refresh
                or force_chrom_refresh
                or force_gene_cov_filter_refresh
                or force_ui_refresh
            ):
                _render_cnv_from_state(state)
                state["_rendered_once"] = True
                if force_color_refresh:
                    state["_force_color_refresh"] = False
                if force_gene_refresh:
                    state["_force_gene_refresh"] = False
                if force_chrom_refresh:
                    state["_force_chrom_refresh"] = False
                if force_gene_cov_filter_refresh:
                    state["_force_gene_cov_filter_refresh"] = False
                if force_ui_refresh:
                    state["_force_ui_refresh"] = False

            _update_cnv_events_analysis(state)

        if data_array_npy.exists() and (data_array_npy_changed or is_fresh_visit):
            try:
                arr = payload.get("bp_array")
                if arr is None:
                    arr = np.load(data_array_npy, allow_pickle=True)
                if hasattr(arr, "dtype") and "name" in arr.dtype.names:
                    state["bp_array"] = arr
                    selected = launcher._cnv_state.setdefault(
                        str(sample_dir), {}
                    ).get("selected_chrom", "All")
                    breakpoint_lines = []
                    for r in arr:
                        if selected == "All" or r["name"] == selected:
                            start_pos = int(r["start"])
                            end_pos = int(r["end"])
                            breakpoint_lines.append((start_pos + end_pos) // 2)
                    _apply_breakpoint_marklines(
                        cnv_diff, selected, state, breakpoint_lines
                    )
                state["bp_array_mtime"] = data_array_npy_mtime
            except Exception:
                pass
        elif data_array_npy.exists() and not data_array_npy_changed:
            if state.get("bp_array") is not None:
                try:
                    selected = launcher._cnv_state.setdefault(
                        str(sample_dir), {}
                    ).get("selected_chrom", "All")
                    arr = state["bp_array"]
                    breakpoint_lines = []
                    for r in arr:
                        if selected == "All" or r["name"] == selected:
                            start_pos = int(r["start"])
                            end_pos = int(r["end"])
                            breakpoint_lines.append((start_pos + end_pos) // 2)
                    _apply_breakpoint_marklines(
                        cnv_diff, selected, state, breakpoint_lines
                    )
                except Exception:
                    pass

        state["cnv_m"] = cnv_npy_mtime
        state["cnv3_m"] = cnv3_npy_mtime
        state["dict_m"] = cnv_dict_npy_mtime
        state["xy_m"] = xy_pkl_mtime
        state["last_visit_time"] = state.get("last_visit_time", time.time())
        state["cnv_plot_theme_dark"] = _is_dark_mode()

        launcher._cnv_state[key] = state
        _update_breakpoints_visibility()

    async def _refresh_cnv_async() -> None:
        try:
            plan = _prepare_cnv_refresh()
            if plan is None:
                return
            if plan["need_load"]:
                payload = await asyncio.to_thread(
                    _cnv_load_binary_payload,
                    sample_dir,
                    cnv_dict_npy_changed=plan["cnv_dict_npy_changed"],
                    cnv_npy_changed=plan["cnv_npy_changed"],
                    cnv3_npy_changed=plan["cnv3_npy_changed"],
                    data_array_reload=plan["data_array_reload"],
                    xy_pkl_changed=plan["xy_pkl_changed"],
                )
            else:
                payload = {}
            _apply_cnv_refresh_after_load(plan, payload)
        except Exception:
            pass

    def _refresh_cnv_sync(sample_dir: Path, launcher: Any) -> None:
        """Synchronous CNV refresh (no event loop)."""
        try:
            plan = _prepare_cnv_refresh()
            if plan is None:
                return
            if plan["need_load"]:
                payload = _cnv_load_binary_payload(
                    sample_dir,
                    cnv_dict_npy_changed=plan["cnv_dict_npy_changed"],
                    cnv_npy_changed=plan["cnv_npy_changed"],
                    cnv3_npy_changed=plan["cnv3_npy_changed"],
                    data_array_reload=plan["data_array_reload"],
                    xy_pkl_changed=plan["xy_pkl_changed"],
                )
            else:
                payload = {}
            _apply_cnv_refresh_after_load(plan, payload)
        except Exception:
            pass

    def _refresh_cnv() -> None:
        """Refresh CNV data; offload binary loads when a loop is running."""
        try:
            if not sample_dir or not sample_dir.exists():
                logging.warning(f"[CNV] Sample directory not found: {sample_dir}")
                return
            try:
                asyncio.get_running_loop()
            except RuntimeError:
                _refresh_cnv_sync(sample_dir, launcher)
                return
            asyncio.create_task(_refresh_cnv_async())
        except Exception as e:
            logging.exception(f"[CNV] Refresh failed: {e}")

    def _update_breakpoints_visibility() -> None:
        """Show/hide breakpoints controls based on chromosome selection."""
        try:
            key = str(sample_dir)
            state = launcher._cnv_state.get(key, {})
            selected = state.get("selected_chrom", "All")
            
            # Show breakpoints controls only when viewing individual chromosomes
            should_show = selected != "All"
            
            try:
                display_value = "block" if should_show else "none"
                cnv_bp_label.style(f"display: {display_value}")
                cnv_bp_row.style(f"display: {display_value}")
            except Exception:
                pass
        except Exception:
            pass

    # Bind control events
    try:

        def _val(ev, default=None):
            # Prefer boolean switch values before legacy toggle label parsing.
            if hasattr(ev, "value") and isinstance(ev.value, bool):
                return ev.value
            args = getattr(ev, "args", None)
            if isinstance(args, bool):
                return args
            if isinstance(args, (list, tuple)) and args and isinstance(args[0], bool):
                return args[0]
            # Handle select/toggle events like [index, {'value': X, 'label': 'Y'}]
            if isinstance(args, list) and len(args) >= 2 and isinstance(args[1], dict):
                return args[1].get("label", default)
            # Handle direct value objects like {'value': 2, 'label': 'GNB1'}
            if hasattr(ev, "value") and isinstance(ev.value, dict):
                return ev.value.get("label", default)
            # Handle args that are directly a dictionary with label
            if isinstance(args, dict) and "label" in args:
                return args.get("label", default)
            # Fallback to standard value extraction
            if hasattr(ev, "value"):
                return ev.value
            return args if args is not None else default

        def _switch_bool(ev, *, default: bool = False) -> bool:
            """Read a ui.switch boolean, falling back to the widget value if needed."""
            raw = _val(ev, None)
            if isinstance(raw, bool):
                return raw
            if raw is None:
                return default
            return str(raw).strip().lower() in ("true", "1", "on", "show")

        def _force_redraw_with_marker_autoscale() -> None:
            """Every CNV control click should re-render and re-fit Y to markers."""
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            st["_force_ui_refresh"] = True
            ui.timer(0.05, _refresh_cnv, once=True)

        def _on_chrom(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            st["selected_chrom"] = _val(ev, "All") or "All"
            st["_force_chrom_refresh"] = True  # Force refresh for chromosome selection
            logging.debug(f"CNV select changed -> {st['selected_chrom']}")
            # Update breakpoints visibility
            _update_breakpoints_visibility()
            # reset x zoom when switching scope
            try:
                for chart in genome_charts:
                    if (
                        isinstance(chart.options.get("dataZoom"), list)
                        and chart.options["dataZoom"]
                    ):
                        chart.options["dataZoom"][0].pop("startValue", None)
                        chart.options["dataZoom"][0].pop("endValue", None)
            except Exception:
                pass
            _force_redraw_with_marker_autoscale()

        def _on_scale(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            st["y_scale"] = "log" if _switch_bool(ev, default=False) else "linear"
            # Prefer live widget value if event parsing failed.
            try:
                if isinstance(getattr(cnv_scale, "value", None), bool):
                    st["y_scale"] = "log" if cnv_scale.value else "linear"
            except Exception:
                pass
            _force_redraw_with_marker_autoscale()

        def _on_bp(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            st["show_bp"] = _switch_bool(ev, default=True)
            try:
                if isinstance(getattr(cnv_bp, "value", None), bool):
                    st["show_bp"] = bool(cnv_bp.value)
            except Exception:
                pass
            _force_redraw_with_marker_autoscale()

        def _apply_cnv_abs_height(tall: bool) -> None:
            """Toggle scatter plot container height and ask ECharts to resize."""
            try:
                if tall:
                    cnv_abs.classes(
                        remove=_CNV_ABS_HEIGHT_CLASS, add=_CNV_ABS_HEIGHT_CLASS_TALL
                    )
                else:
                    cnv_abs.classes(
                        remove=_CNV_ABS_HEIGHT_CLASS_TALL, add=_CNV_ABS_HEIGHT_CLASS
                    )
            except Exception:
                logging.debug("CNV scatter height class update failed", exc_info=True)
            try:
                cnv_abs.run_chart_method("resize")
            except Exception:
                logging.debug("CNV scatter resize after height change failed", exc_info=True)

        def _on_height(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            tall = _switch_bool(ev, default=False)
            try:
                if isinstance(getattr(cnv_height, "value", None), bool):
                    tall = bool(cnv_height.value)
            except Exception:
                pass
            st["double_height"] = tall
            _apply_cnv_abs_height(tall)

        def _on_color(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            st["color_mode"] = "value" if _switch_bool(ev, default=False) else "chromosome"
            try:
                if isinstance(getattr(cnv_color, "value", None), bool):
                    st["color_mode"] = "value" if cnv_color.value else "chromosome"
            except Exception:
                pass
            # Force a refresh by setting a flag that bypasses the state sync logic
            st["_force_color_refresh"] = True
            _force_redraw_with_marker_autoscale()

        def _on_plot_bin(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            v = getattr(ev, "args", None) if hasattr(ev, "args") else getattr(ev, "value", None)
            if v is None and hasattr(ev, "value"):
                v = ev.value
            st["plot_bin_width"] = _cnv_plot_bin_bp_from_ui(v)
            st["_force_chrom_refresh"] = True  # force re-render with new bin width
            _force_redraw_with_marker_autoscale()

        # Bind both native change and model-value updates for robustness
        cnv_chrom_select.on("change", _on_chrom)
        cnv_chrom_select.on("update:model-value", _on_chrom)

        # Gene selection zoom
        def _on_gene(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            selected_gene = _val(ev, "All") or "All"
            st["selected_gene"] = selected_gene
            st["_force_gene_refresh"] = True  # Force refresh for gene selection
            _force_redraw_with_marker_autoscale()

        cnv_gene_select.on("change", _on_gene)
        cnv_gene_select.on("update:model-value", _on_gene)

        def _on_gene_cov_filter(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            # Prefer the live widget boolean; event args can be ambiguous.
            enabled = True
            try:
                raw = getattr(cnv_gene_cov_filter, "value", None)
                if isinstance(raw, bool):
                    enabled = raw
                else:
                    enabled = _switch_bool(ev, default=True)
            except Exception:
                enabled = _switch_bool(ev, default=True)
            st["gene_coverage_filter"] = (
                _CNV_GENE_COVERAGE_FILTER_OUTLIERS
                if enabled
                else _CNV_GENE_COVERAGE_FILTER_ALL
            )
            try:
                if cnv_gene_cov_filter.value is not enabled:
                    cnv_gene_cov_filter.value = enabled
            except Exception:
                pass
            st["_force_gene_cov_filter_refresh"] = True
            _force_redraw_with_marker_autoscale()

        cnv_gene_cov_filter.on("change", _on_gene_cov_filter)
        cnv_gene_cov_filter.on("update:model-value", _on_gene_cov_filter)
        cnv_scale.on("change", _on_scale)
        cnv_scale.on("update:model-value", _on_scale)
        cnv_plot_bin.on("change", _on_plot_bin)
        cnv_plot_bin.on("update:model-value", _on_plot_bin)
        cnv_bp.on("change", _on_bp)
        cnv_bp.on("update:model-value", _on_bp)
        cnv_height.on("change", _on_height)
        cnv_height.on("update:model-value", _on_height)
        cnv_color.on("change", _on_color)
        cnv_color.on("update:model-value", _on_color)
    except Exception:
        pass

    def _sync_cnv_echarts_theme_if_needed() -> None:
        """Re-apply axis/tooltip/title colours when the user toggles light/dark mode."""
        try:
            dark = _is_dark_mode()
        except Exception:
            dark = False
        key = str(sample_dir)
        st = launcher._cnv_state.get(key, {})
        if st.get("cnv_plot_theme_dark") == dark:
            return
        _apply_cnv_echart_chrome(cnv_abs, dark)
        _apply_cnv_echart_chrome(cnv_diff, dark)
        try:
            _cnv_echart_push_update(cnv_abs)
            _cnv_echart_push_update(cnv_diff)
        except Exception:
            pass
        launcher._cnv_state.setdefault(key, {})["cnv_plot_theme_dark"] = dark

    # Start the refresh timer (every 30 seconds)
    refresh_timer = client_timer(30.0, _refresh_cnv, active=True, immediate=False)
    client_timer(0.5, _refresh_cnv, once=True)
    unregister_cnv_theme_sync = register_theme_sync_callback(
        _sync_cnv_echarts_theme_if_needed,
        element=cnv_abs,
        interval_s=0.5,
        immediate=True,
    )
    try:
        ui.context.client.on_disconnect(
            lambda: (stop_timer(refresh_timer), unregister_cnv_theme_sync())
        )
    except Exception:
        pass
