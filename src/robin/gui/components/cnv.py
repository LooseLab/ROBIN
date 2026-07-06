from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple
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


def _cnv_echart_push_update(chart: Any) -> None:
    """Replace the full ECharts option (NiceGUI merges by default and leaves stale scatter data)."""
    options_clean = _cnv_echarts_option_to_json(chart.options)
    opts_json = json.dumps(options_clean)
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
                ya["axisLabel"] = {**(ya.get("axisLabel") or {}), "color": p["muted"]}
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
                ui.label("Color by").classes("classification-insight-meta ml-2")
                cnv_color = ui.toggle(
                    options={"chromosome": "Chromosome", "value": "Up/Down"},
                    value="chromosome",
                ).classes("mt-1")
                ui.label("Y-axis").classes("classification-insight-meta ml-2")
                cnv_scale = ui.toggle(
                    options={"linear": "Linear", "log": "Log2 ratio"},
                    value="linear",
                ).classes("mt-1")
                ui.label("Plot bin").classes("classification-insight-meta ml-2")
                cnv_plot_bin = ui.select(
                    options=_CNV_PLOT_BIN_OPTIONS,
                    value=_CNV_PLOT_BIN_KEY_DEFAULT,
                ).style("width: 120px")
                cnv_bp_label = ui.label("Breakpoints").classes(
                    "classification-insight-meta ml-2"
                ).style("display: none")
                cnv_bp = ui.toggle(
                    options={"hide": "Hide", "show": "Show"}, value="show"
                ).classes("mt-1").style("display: none")
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
                                "name": "Breakpoint density",
                                "position": "right",
                            },
                        ],
                        "dataZoom": [
                            {"type": "slider", "xAxisIndex": [0]},
                            {
                                "type": "slider",
                                "yAxisIndex": [0, 1],
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
                ).classes("w-full h-72 cnv-genome-abs-chart")
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
                                "name": "Breakpoint density",
                                "position": "right",
                            },
                        ],
                        "dataZoom": [
                            {"type": "slider", "xAxisIndex": [0]},
                            {
                                "type": "slider",
                                "yAxisIndex": [0, 1],
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
                ).classes("w-full h-72 cnv-genome-diff-chart")
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
            
            # Show/hide breakpoint density y-axis based on selection
            should_show_breakpoint_density = selected != "All"
            for chart in genome_charts:
                if len(chart.options["yAxis"]) > 1:
                    # Index 1 is the "Breakpoint density" axis
                    chart.options["yAxis"][1]["show"] = should_show_breakpoint_density
            
            cnv_abs.options["yAxis"][0]["type"] = "value"
            cnv_abs.options["yAxis"][0].pop("logBase", None)
            if use_log:
                cnv_abs.options["yAxis"][0]["name"] = "Log2 ratio (ploidy / expected)"
                cnv_abs.options["title"]["text"] = "CNV scatter plot"
                cnv_abs.options["title"]["top"] = 4
                cnv_abs.options["title"]["subtext"] = (
                    "log2(ploidy / expected copy number); 0 = normal"
                )
                cnv_abs.options["grid"]["top"] = "26%"
                try:
                    y_dz = cnv_abs.options["dataZoom"][1]
                    y_dz["startValue"] = -2
                    y_dz["endValue"] = 2
                except Exception:
                    pass
            else:
                cnv_abs.options["yAxis"][0]["name"] = "Ploidy"
                cnv_abs.options["title"]["text"] = "CNV scatter plot"
                cnv_abs.options["title"]["top"] = 10
                cnv_abs.options["title"].pop("subtext", None)
                cnv_abs.options["grid"]["top"] = "20%"
                try:
                    y_dz = cnv_abs.options["dataZoom"][1]
                    y_dz["startValue"] = 0
                    y_dz["endValue"] = 6
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
                        # Genes of interest and gene selector options
                        try:
                            gene_df = _load_gene_bed(sample_dir)
                            gchr = gene_df[gene_df["chrom"] == selected]
                            # limit to reduce clutter; still add many labels
                            gene_opts = {"All": "All"}
                            for _, gr in gchr.iterrows():
                                gene_opts[str(gr["gene"])] = str(gr["gene"])
                            try:
                                cnv_gene_select.set_options(gene_opts)
                            except Exception as e:
                                pass
                            
                            # annotate genes on main series
                            if series_abs:
                                main = series_abs[0]
                                # attach gene regions via markArea on main series after replacement
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
                                                }
                                            },
                                            {"xAxis": float(gr["end_pos"])},
                                        ]
                                    )
                                main.setdefault("markArea", {"data": []})
                                main["markArea"]["data"] = (
                                    main["markArea"]["data"] or []
                                ) + mark
                                # ensure series_abs[0] updated
                                series_abs[0] = main
                            
                        except Exception as e:
                            pass
                        
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
                
                if sel_gene and sel_gene != "All":
                    # Load gene data for zoom if not already loaded
                    gene_df = _load_gene_bed(sample_dir)
                    gchr = gene_df[gene_df["chrom"] == selected]
                    
                    row = gchr[gchr["gene"] == sel_gene]
                    if not row.empty:
                        s_bp = int(row.iloc[0]["start_pos"])
                        e_bp = int(row.iloc[0]["end_pos"])
                        # Use 10x analysis bin width for padding (in bp)
                        pad = 10 * binw_analysis
                        zoom_start = max(0, s_bp - pad)
                        zoom_end = e_bp + pad
                        try:
                            cnv_abs.options["dataZoom"][0].update(
                                {
                                    "startValue": zoom_start,
                                    "endValue": zoom_end,
                                    "start": None,  # Remove percentage-based zoom
                                    "end": None,    # Remove percentage-based zoom
                                }
                            )
                        except Exception as e:
                            pass
                else:
                    # Reset zoom when "All" is selected
                    try:
                        if isinstance(cnv_abs.options.get("dataZoom"), list) and cnv_abs.options["dataZoom"]:
                            dz = cnv_abs.options["dataZoom"][0]
                            dz.pop("startValue", None)
                            dz.pop("endValue", None)
                            dz.update({"start": 0, "end": 100})
                    except Exception as e:
                        pass
            except Exception as e:
                pass
            
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
                _thin_chart_series(cnv_diff, MAX_POINTS_PER_CHART)
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
                    if sel_gene and sel_gene != "All" and selected != "All":
                        gene_df = _load_gene_bed(sample_dir)
                        gchr = gene_df[gene_df["chrom"] == selected]
                        row = gchr[gchr["gene"] == sel_gene]
                        if not row.empty:
                            s_bp = int(row.iloc[0]["start_pos"])
                            e_bp = int(row.iloc[0]["end_pos"])
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
        if current_time - last_refresh < 0.1:
            return None
        state["_last_refresh"] = current_time

        ui_changed = False
        try:
            ui_sel = getattr(cnv_chrom_select, "value", None)
            if ui_sel and ui_sel != state.get("selected_chrom"):
                state["selected_chrom"] = ui_sel
                ui_changed = True
            ui_scale = getattr(cnv_scale, "value", None)
            if ui_scale and ui_scale != state.get("y_scale"):
                state["y_scale"] = ui_scale
                ui_changed = True
            ui_bp = getattr(cnv_bp, "value", None)
            if ui_bp is not None:
                desired = ui_bp == "show"
                if desired != state.get("show_bp", True):
                    state["show_bp"] = desired
                    ui_changed = True
            ui_color = getattr(cnv_color, "value", None)
            current_color_mode = state.get("color_mode", "chromosome")
            if ui_color and ui_color != current_color_mode:
                state["color_mode"] = ui_color
                ui_changed = True
            ui_plot_bin = getattr(cnv_plot_bin, "value", None)
            if ui_plot_bin is not None:
                want_bin = _cnv_plot_bin_bp_from_ui(ui_plot_bin)
                if want_bin != state.get("plot_bin_width"):
                    state["plot_bin_width"] = want_bin
                    ui_changed = True
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

        force_color_refresh = state.get("_force_color_refresh", False)
        force_gene_refresh = state.get("_force_gene_refresh", False)
        force_chrom_refresh = state.get("_force_chrom_refresh", False)

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
            or force_color_refresh
            or force_gene_refresh
            or force_chrom_refresh
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
            if (
                changed
                or ui_changed
                or not state.get("_rendered_once")
                or force_color_refresh
                or force_gene_refresh
                or force_chrom_refresh
            ):
                _render_cnv_from_state(state)
                state["_rendered_once"] = True
                if force_color_refresh:
                    state["_force_color_refresh"] = False
                if force_gene_refresh:
                    state["_force_gene_refresh"] = False
                if force_chrom_refresh:
                    state["_force_chrom_refresh"] = False

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
                cnv_bp.style(f"display: {display_value}")
            except Exception:
                pass
        except Exception:
            pass

    # Bind control events
    try:

        def _val(ev, default=None):
            # Handle toggle events which have args like [index, {'value': X, 'label': 'Y'}]
            if hasattr(ev, "args") and ev.args and isinstance(ev.args, list) and len(ev.args) >= 2:
                if isinstance(ev.args[1], dict):
                    # Extract the label from the toggle event structure
                    return ev.args[1].get("label", default)
            # Handle direct value objects like {'value': 2, 'label': 'GNB1'}
            if hasattr(ev, "value") and isinstance(ev.value, dict):
                return ev.value.get("label", default)
            # Handle args that are directly a dictionary with label
            if hasattr(ev, "args") and isinstance(ev.args, dict) and "label" in ev.args:
                return ev.args.get("label", default)
            # Fallback to standard value extraction
            return (
                getattr(ev, "value", None)
                if hasattr(ev, "value")
                else (getattr(ev, "args", None) or default)
            )

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
            # Trigger immediate refresh to update all UI elements
            ui.timer(0.1, _refresh_cnv, once=True)

        def _on_scale(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            st["y_scale"] = _val(ev, "linear") or "linear"
            # Trigger immediate refresh to update all UI elements
            ui.timer(0.1, _refresh_cnv, once=True)

        def _on_bp(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            show_bp = _val(ev, "show") == "show"
            st["show_bp"] = show_bp
            # Trigger immediate refresh to update all UI elements
            ui.timer(0.1, _refresh_cnv, once=True)

        def _on_color(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            val = _val(ev, "chromosome") or "chromosome"
            # Accept either keys or labels from the toggle
            vlow = str(val).strip().lower()
            if vlow in ("chromosome", "chromosomes"):
                st["color_mode"] = "chromosome"
            elif vlow in ("value", "up/down", "updown", "up_down", "up down"):
                st["color_mode"] = "value"
            else:
                st["color_mode"] = "chromosome"
            # Force a refresh by setting a flag that bypasses the state sync logic
            st["_force_color_refresh"] = True
            
            # Trigger immediate refresh to update all UI elements
            ui.timer(0.1, _refresh_cnv, once=True)

        def _on_plot_bin(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            v = getattr(ev, "args", None) if hasattr(ev, "args") else getattr(ev, "value", None)
            if v is None and hasattr(ev, "value"):
                v = ev.value
            st["plot_bin_width"] = _cnv_plot_bin_bp_from_ui(v)
            st["_force_chrom_refresh"] = True  # force re-render with new bin width
            ui.timer(0.1, _refresh_cnv, once=True)

        # Bind both native change and model-value updates for robustness
        cnv_chrom_select.on("change", _on_chrom)
        cnv_chrom_select.on("update:model-value", _on_chrom)

        # Gene selection zoom
        def _on_gene(ev):
            st = launcher._cnv_state.setdefault(str(sample_dir), {})
            selected_gene = _val(ev, "All") or "All"
            st["selected_gene"] = selected_gene
            st["_force_gene_refresh"] = True  # Force refresh for gene selection
            # Trigger immediate refresh to update all UI elements
            ui.timer(0.1, _refresh_cnv, once=True)

        cnv_gene_select.on("change", _on_gene)
        cnv_gene_select.on("update:model-value", _on_gene)
        cnv_scale.on("change", _on_scale)
        cnv_scale.on("update:model-value", _on_scale)
        cnv_plot_bin.on("change", _on_plot_bin)
        cnv_plot_bin.on("update:model-value", _on_plot_bin)
        cnv_bp.on("change", _on_bp)
        cnv_bp.on("update:model-value", _on_bp)
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
    refresh_timer = ui.timer(30.0, _refresh_cnv, active=True, immediate=False)
    ui.timer(0.5, _refresh_cnv, once=True)
    unregister_cnv_theme_sync = register_theme_sync_callback(
        _sync_cnv_echarts_theme_if_needed,
        element=cnv_abs,
        interval_s=0.5,
        immediate=True,
    )
    try:
        ui.context.client.on_disconnect(
            lambda: (refresh_timer.deactivate(), unregister_cnv_theme_sync())
        )
    except Exception:
        pass
