from __future__ import annotations

import asyncio
from pathlib import Path
from typing import Any, Dict, List, Optional
import csv
import logging
import zlib

# Classification charts — design.md §8.8 (ECharts; brand / diagnostic palette)
_BRAND_GREEN = "#10b981"
_BRAND_GREEN_NEON = "#22c55e"
# Distinct hues for tumour / class rows (not brand green — reserved for top/leader).
# Light UI: deeper 600–700 shades so bars/lines read on white/light cards.
_CLASS_TUMOUR_PALETTE_LIGHT = [
    "#0369a1",  # sky-700
    "#6d28d9",  # violet-700
    "#c2410c",  # orange-700
    "#be185d",  # pink-700
    "#0f766e",  # teal-700
    "#a16207",  # yellow-700
    "#4338ca",  # indigo-700
    "#7e22ce",  # purple-700
    "#be123c",  # rose-700
    "#0e7490",  # cyan-700
]
# Dark UI: brighter mid tones so series pop on charcoal surfaces.
_CLASS_TUMOUR_PALETTE_DARK = [
    "#38bdf8",  # sky-400
    "#a78bfa",  # violet-400
    "#fb923c",  # orange-400
    "#f472b6",  # pink-400
    "#2dd4bf",  # teal-400
    "#eab308",  # yellow-500
    "#818cf8",  # indigo-400
    "#c084fc",  # purple-400
    "#fb7185",  # rose-400
    "#22d3ee",  # cyan-400
]
_BAR_RADIUS = 8
_BAR_WIDTH_PX = 22
_TOOLTIP_DARK_BG = "rgba(15, 23, 42, 0.96)"
_TOOLTIP_DARK_BORDER = "#475569"
_TOOLTIP_DARK_TEXT = "#f8fafc"
# Medium / high confidence band fills (markArea) — keep low alpha so plots stay readable.
_ZONE_GREEN = "rgba(16, 185, 129, 0.07)"
_ZONE_SLATE = "rgba(100, 116, 139, 0.06)"
_ZONE_GREEN_DARK = "rgba(34, 197, 94, 0.1)"
_ZONE_SLATE_DARK = "rgba(148, 163, 184, 0.08)"

try:
    from nicegui import ui
except ImportError:  # pragma: no cover
    ui = None


def _is_dark_mode() -> bool:
    """Quasar body--dark is driven by app storage (see theme.frame)."""
    try:
        from nicegui import app

        return bool(app.storage.user.get("dark_mode"))
    except Exception:
        return False


def _echart_surface_palette() -> Dict[str, str]:
    """Axis, title, legend, and guide-line colors for ECharts on light vs dark UI."""
    if _is_dark_mode():
        return {
            "title": "#f1f5f9",
            "axis": "#cbd5e1",
            "axis_muted": "#94a3b8",
            "axis_line": "#64748b",
            "split_line": "rgba(148, 163, 184, 0.38)",
            "legend": "#cbd5e1",
            "bar_label": "#e2e8f0",
            "mark_line": "#94a3b8",
            "mark_label": "#cbd5e1",
            "axis_pointer": "#94a3b8",
        }
    # Light surfaces: higher contrast than mid-slate on white (WCAG-ish body text tones).
    return {
        "title": "#0f172a",
        "axis": "#1e293b",
        "axis_muted": "#475569",
        "axis_line": "#94a3b8",
        "split_line": "rgba(71, 85, 105, 0.33)",
        "legend": "#1e293b",
        "bar_label": "#0f172a",
        "mark_line": "#475569",
        "mark_label": "#334155",
        "axis_pointer": "#475569",
    }


def _apply_palette_to_bar_chart(bar: Any, palette: Dict[str, str]) -> None:
    """Mutate bar ECharts options for readable labels on current theme."""
    o = bar.options
    al = palette.get("axis_line", palette["axis"])
    o.setdefault("title", {}).setdefault("textStyle", {})["color"] = palette["title"]
    o.setdefault("xAxis", {}).setdefault("axisLabel", {})["color"] = palette["axis"]
    o.setdefault("xAxis", {}).setdefault("axisLine", {}).setdefault("lineStyle", {})[
        "color"
    ] = al
    o.setdefault("xAxis", {}).setdefault("axisTick", {}).setdefault("lineStyle", {})[
        "color"
    ] = al
    o.setdefault("yAxis", {}).setdefault("axisLabel", {})["color"] = palette["axis"]
    o.setdefault("yAxis", {}).setdefault("axisLine", {}).setdefault("lineStyle", {})[
        "color"
    ] = al
    o.setdefault("yAxis", {}).setdefault("axisTick", {}).setdefault("lineStyle", {})[
        "color"
    ] = al
    if o.get("series"):
        s0 = o["series"][0]
        s0.setdefault("label", {})["color"] = palette["bar_label"]
        ml = s0.get("markLine")
        if ml and isinstance(ml, dict):
            ml.setdefault("lineStyle", {})["color"] = palette["mark_line"]
            ml.setdefault("label", {})["color"] = palette["mark_label"]


def _apply_palette_to_ts_chart(ts: Any, palette: Dict[str, str]) -> None:
    """Mutate time-series ECharts options for readable labels on current theme."""
    o = ts.options
    al = palette.get("axis_line", palette["axis"])
    o.setdefault("title", {}).setdefault("textStyle", {})["color"] = palette["title"]
    o.setdefault("legend", {}).setdefault("textStyle", {})["color"] = palette["legend"]
    xa = o.get("xAxis")
    if isinstance(xa, dict):
        xa.setdefault("axisLabel", {})["color"] = palette["axis"]
        xa.setdefault("axisLine", {}).setdefault("lineStyle", {})["color"] = al
        xa.setdefault("axisTick", {}).setdefault("lineStyle", {})["color"] = al
    ya = o.setdefault("yAxis", {})
    ya.setdefault("axisLabel", {})["color"] = palette["axis"]
    ya.setdefault("axisLine", {}).setdefault("lineStyle", {})["color"] = al
    ya.setdefault("axisTick", {}).setdefault("lineStyle", {})["color"] = al
    ya.setdefault("splitLine", {}).setdefault("lineStyle", {})["color"] = palette[
        "split_line"
    ]
    ya.setdefault("splitLine", {}).setdefault("lineStyle", {})["type"] = "dashed"
    tip = o.setdefault("tooltip", {})
    tip.setdefault("axisPointer", {}).setdefault("lineStyle", {})["color"] = palette[
        "axis_pointer"
    ]


def _active_tumour_palette() -> List[str]:
    """Categorical bar/line colours: deep on light UI, brighter on dark UI."""
    return _CLASS_TUMOUR_PALETTE_DARK if _is_dark_mode() else _CLASS_TUMOUR_PALETTE_LIGHT


def _tumour_category_color(index: int) -> str:
    """Distinct colour by index (fallback)."""
    pal = _active_tumour_palette()
    return pal[index % len(pal)]


def _color_for_class_name(name: str) -> str:
    """Stable hue per class label (bars + lines) using adler32 — same name → same colour."""
    pal = _active_tumour_palette()
    raw = zlib.adler32(name.encode("utf-8", errors="replace")) & 0xFFFFFFFF
    return pal[raw % len(pal)]

try:
    from robin.classification_config import (
        CLASSIFIER_CONFIDENCE_THRESHOLDS,
        DEFAULT_CONFIDENCE_THRESHOLDS,
    )
except ImportError:  # pragma: no cover
    CLASSIFIER_CONFIDENCE_THRESHOLDS = {}
    DEFAULT_CONFIDENCE_THRESHOLDS = {"high": 80.0, "medium": 50.0, "low": 0.0}


def _confidence_thresholds_for_classifier(classifier_key: str) -> Dict[str, float]:
    """Per-classifier medium / high bands (percent), from classification_config."""
    raw = CLASSIFIER_CONFIDENCE_THRESHOLDS.get(
        classifier_key, DEFAULT_CONFIDENCE_THRESHOLDS
    )
    return {
        "medium": float(raw.get("medium", 50.0)),
        "high": float(raw.get("high", 80.0)),
    }


def _echart_media_responsive_bar() -> List[Dict[str, Any]]:
    """Narrow viewports: more vertical space, readable category labels."""
    return [
        {
            "query": {"maxWidth": 640},
            "option": {
                "title": {"textStyle": {"fontSize": 13}},
                "grid": {
                    "left": "2%",
                    "right": "10%",
                    "top": "16%",
                    "bottom": "10%",
                    "containLabel": True,
                },
                "yAxis": {
                    "axisLabel": {
                        "fontSize": 11,
                        "lineHeight": 14,
                        "width": 100,
                        "overflow": "truncate",
                    }
                },
                "xAxis": {"axisLabel": {"fontSize": 10}},
            },
        },
        {
            "query": {"maxWidth": 400},
            "option": {
                "grid": {"bottom": "12%", "top": "18%"},
                "yAxis": {
                    "axisLabel": {
                        "fontSize": 10,
                        "width": 88,
                        "overflow": "truncate",
                    }
                },
            },
        },
    ]


def _echart_media_responsive_ts(palette: Dict[str, str]) -> List[Dict[str, Any]]:
    """Narrow viewports: legend below chart (full plot width); avoid side legend squeeze."""
    leg_color = palette.get("legend", "#475569")
    return [
        {
            "query": {"maxWidth": 640},
            "option": {
                "title": {"textStyle": {"fontSize": 13}},
                "legend": {
                    "type": "scroll",
                    "orient": "horizontal",
                    "left": "center",
                    "right": "auto",
                    "top": "auto",
                    "bottom": 4,
                    "width": "92%",
                    "itemGap": 8,
                    "itemWidth": 10,
                    "itemHeight": 10,
                    "textStyle": {"fontSize": 10, "color": leg_color},
                },
                "grid": {
                    "left": "10%",
                    "right": "6%",
                    "top": "20%",
                    "bottom": "30%",
                    "containLabel": True,
                },
            },
        },
        {
            "query": {"maxWidth": 400},
            "option": {
                "legend": {
                    "bottom": 2,
                    "textStyle": {"fontSize": 9, "color": leg_color},
                },
                "grid": {"bottom": "34%", "left": "12%"},
            },
        },
    ]


def _bar_chart_mark_line(
    medium: float, high: float, palette: Dict[str, str]
) -> Dict[str, Any]:
    """Vertical reference lines on the horizontal bar chart (value axis)."""
    return {
        "silent": True,
        "symbol": "none",
        "lineStyle": {
            "color": palette["mark_line"],
            "width": 1,
            "type": "dashed",
        },
        "label": {"color": palette["mark_label"], "fontSize": 11},
        "data": [
            {
                "xAxis": medium,
                "label": {
                    "show": True,
                    "formatter": f"Medium ({medium:.0f}%)",
                    "position": "insideEndTop",
                },
            },
            {
                "xAxis": high,
                "label": {
                    "show": True,
                    "formatter": f"High ({high:.0f}%)",
                    "position": "insideEndTop",
                },
            },
        ],
    }


# Import section visibility helpers
try:
    from robin.gui.config import get_enabled_classification_steps, CLASSIFICATION_STEPS
except ImportError:
    get_enabled_classification_steps = lambda steps: {"sturgeon", "nanodx", "random_forest", "pannanodx"}
    CLASSIFICATION_STEPS = {
        "sturgeon": "Sturgeon",
        "nanodx": "NanoDX",
        "random_forest": "Random Forest",
        "pannanodx": "PanNanoDX",
    }


def add_classification_section(sample_dir: Path, launcher: Any = None) -> None:
    """Build the Classification section (Sturgeon, NanoDX, PanNanoDX, RF)."""
    # Get workflow steps from launcher if available
    workflow_steps = launcher.workflow_steps if launcher and hasattr(launcher, 'workflow_steps') else None
    enabled_classification_steps = get_enabled_classification_steps(workflow_steps)
    
    # Map workflow step names to tool display names
    tool_to_step_map = {
        "Sturgeon": "sturgeon",
        "NanoDX": "nanodx",
        "PanNanoDX": "pannanodx",
        "Random Forest": "random_forest",
    }
    
    with ui.element("div").classes("classification-insight-shell w-full min-w-0").props(
        "id=classification-section"
    ):
        ui.label("Classification").classes(
            "classification-insight-heading text-headline-small"
        )
        tool_to_file = {
            "Sturgeon": {"file": "sturgeon_scores.csv", "mode": "fraction"},
            "NanoDX": {"file": "NanoDX_scores.csv", "mode": "fraction"},
            "PanNanoDX": {"file": "PanNanoDX_scores.csv", "mode": "fraction"},
            "Random Forest": {"file": "random_forest_scores.csv", "mode": "percent"},
        }
        charts: Dict[str, Dict[str, Any]] = {}
        for tool_name, cfg in tool_to_file.items():
            # Check if this tool should be shown
            tool_step = tool_to_step_map.get(tool_name)
            if workflow_steps and tool_step and tool_step not in enabled_classification_steps:
                continue  # Skip this tool if not enabled
            exp = (
                ui.expansion(tool_name, icon="analytics")
                .classes("w-full")
                .props(f'id=classification-detail-{tool_step}')
            )
            with exp:
                summary_labels = None
                tool_icon = {
                    "Sturgeon": "psychology",
                    "NanoDX": "biotech",
                    "PanNanoDX": "science",
                    "Random Forest": "forest",
                }[tool_name]
                if tool_name == "Sturgeon":
                    with ui.element("div").classes(
                        "classification-insight-card w-full min-w-0 mb-2"
                    ):
                        with ui.row().classes(
                            "w-full items-start gap-2 p-2 md:p-3 flex-wrap"
                        ):
                            ui.icon(tool_icon).classes("classification-insight-icon")
                            with ui.column().classes("gap-1 flex-1 min-w-0"):
                                st_class = ui.label(
                                    "Sturgeon classification: Unknown"
                                ).classes("classification-insight-result text-sm")
                                st_conf = ui.label("Confidence: --%").classes(
                                    "classification-insight-meta"
                                )
                                st_probes = ui.label("Probes: --").classes(
                                    "classification-insight-meta"
                                )
                        summary_labels = {
                            "class": st_class,
                            "conf": st_conf,
                            "probes": st_probes,
                        }
                elif tool_name in ("NanoDX", "PanNanoDX"):
                    with ui.element("div").classes(
                        "classification-insight-card w-full min-w-0 mb-2"
                    ):
                        with ui.row().classes(
                            "w-full items-start gap-2 p-2 md:p-3 flex-wrap"
                        ):
                            ui.icon(tool_icon).classes("classification-insight-icon")
                            with ui.column().classes("gap-1 flex-1 min-w-0"):
                                ndx_class = ui.label(
                                    f"{tool_name} classification: Unknown"
                                ).classes("classification-insight-result text-sm")
                                ndx_conf = ui.label("Confidence: --%").classes(
                                    "classification-insight-meta"
                                )
                                ndx_feats = ui.label("Probes: --").classes(
                                    "classification-insight-meta"
                                )
                        summary_labels = {
                            "class": ndx_class,
                            "conf": ndx_conf,
                            "probes": ndx_feats,
                        }
                elif tool_name == "Random Forest":
                    with ui.element("div").classes(
                        "classification-insight-card w-full min-w-0 mb-2"
                    ):
                        with ui.row().classes(
                            "w-full items-start gap-2 p-2 md:p-3 flex-wrap"
                        ):
                            ui.icon(tool_icon).classes("classification-insight-icon")
                            with ui.column().classes("gap-1 flex-1 min-w-0"):
                                rf_class = ui.label(
                                    "Forest classification: Unknown"
                                ).classes("classification-insight-result text-sm")
                                rf_conf = ui.label("Confidence: --%").classes(
                                    "classification-insight-meta"
                                )
                                rf_feats = ui.label("Features: --").classes(
                                    "classification-insight-meta"
                                )
                        summary_labels = {
                            "class": rf_class,
                            "conf": rf_conf,
                            "probes": rf_feats,
                        }
                ui.label(f"{tool_name} current classification").classes(
                    "classification-insight-meta"
                )
                _thr = _confidence_thresholds_for_classifier(tool_to_step_map[tool_name])
                pal = _echart_surface_palette()
                bar = ui.echart(
                    {
                        "backgroundColor": "transparent",
                        "title": {
                            "text": f"{tool_name} (Top classes)",
                            "left": "center",
                            "top": 8,
                            "textStyle": {
                                "fontSize": 15,
                                "color": pal["title"],
                                "fontWeight": 600,
                            },
                        },
                        "tooltip": {
                            "trigger": "axis",
                            "axisPointer": {"type": "shadow"},
                            "formatter": "{b}\nConfidence: {c}%",
                            "backgroundColor": _TOOLTIP_DARK_BG,
                            "borderColor": _TOOLTIP_DARK_BORDER,
                            "borderWidth": 1,
                            "textStyle": {"color": _TOOLTIP_DARK_TEXT},
                        },
                        "grid": {
                            "left": "3%",
                            "right": "8%",
                            "bottom": "5%",
                            "top": "20%",
                            "containLabel": True,
                        },
                        "xAxis": {
                            "type": "value",
                            "min": 0,
                            "max": 100,
                            "interval": 20,
                            "axisLabel": {
                                "formatter": "{value}%",
                                "color": pal["axis"],
                            },
                        },
                        "yAxis": {
                            "type": "category",
                            "inverse": True,
                            "data": [],
                            "axisLabel": {
                                "fontSize": 11,
                                "color": pal["axis"],
                            },
                        },
                        "series": [
                            {
                                "type": "bar",
                                "data": [],
                                "barWidth": _BAR_WIDTH_PX,
                                "barMaxWidth": _BAR_WIDTH_PX + 10,
                                "itemStyle": {"borderRadius": _BAR_RADIUS},
                                "label": {
                                    "show": True,
                                    "position": "right",
                                    "formatter": "{c}%",
                                    "color": pal["bar_label"],
                                },
                                "markLine": _bar_chart_mark_line(
                                    _thr["medium"], _thr["high"], pal
                                ),
                            }
                        ],
                        "media": _echart_media_responsive_bar(),
                    }
                ).classes(
                    "w-full min-h-[300px] h-[340px] sm:min-h-[320px] sm:h-80"
                )
                ui.label(f"{tool_name} confidence over time").classes(
                    "classification-insight-meta mt-2"
                )
                ts = ui.echart(
                    {
                        "backgroundColor": "transparent",
                        "title": {
                            "text": f"{tool_name} (confidence over time)",
                            "left": "center",
                            "top": 6,
                            "textStyle": {
                                "fontSize": 15,
                                "color": pal["title"],
                                "fontWeight": 600,
                            },
                        },
                        "tooltip": {
                            "trigger": "axis",
                            "backgroundColor": _TOOLTIP_DARK_BG,
                            "borderColor": _TOOLTIP_DARK_BORDER,
                            "borderWidth": 1,
                            "textStyle": {"color": _TOOLTIP_DARK_TEXT},
                            "axisPointer": {
                                "type": "line",
                                "lineStyle": {"color": pal["axis_pointer"]},
                            },
                        },
                        "legend": {
                            "type": "scroll",
                            "orient": "vertical",
                            "right": 8,
                            "top": 52,
                            "width": 140,
                            "height": 110,
                            "itemGap": 6,
                            "itemWidth": 10,
                            "itemHeight": 10,
                            "textStyle": {
                                "fontSize": 11,
                                "color": pal["legend"],
                            },
                        },
                        "grid": {
                            "left": "5%",
                            "right": "20%",
                            "bottom": "12%",
                            "top": "22%",
                            "containLabel": True,
                        },
                        "xAxis": {
                            "type": "time",
                            "splitNumber": 4,
                            "axisLabel": {
                                "hideOverlap": True,
                                "color": pal["axis"],
                            },
                        },
                        "yAxis": {
                            "type": "value",
                            "min": 0,
                            "max": 100,
                            "axisLabel": {
                                "formatter": "{value}%",
                                "color": pal["axis"],
                            },
                            "splitLine": {
                                "show": True,
                                "lineStyle": {
                                    "type": "dashed",
                                    "color": pal["split_line"],
                                },
                            },
                        },
                        "series": [],
                        "media": _echart_media_responsive_ts(pal),
                    }
                ).classes(
                    "w-full min-h-[320px] h-[380px] sm:min-h-[280px] sm:h-72"
                )
                charts[tool_name] = {
                    "bar": bar,
                    "ts": ts,
                    "file": cfg["file"],
                    "mode": cfg["mode"],
                    "summary": summary_labels,
                    "expansion": exp,
                    "last_mtime": None,
                }

    _last_dark_sig: list = [None]

    def _sync_classification_theme(force: bool = False) -> None:
        """Re-apply axis/title colors when user toggles light/dark (storage-backed).

        ``force=True`` reapplies after client/storage settle so first paint matches theme.
        """
        try:
            from nicegui import app

            cur = bool(app.storage.user.get("dark_mode"))
        except Exception:
            cur = False
        if not force and _last_dark_sig[0] == cur:
            return
        _last_dark_sig[0] = cur
        for tname, ch in charts.items():
            try:
                fn = ch.get("file")
                fp = (sample_dir / fn) if (sample_dir and fn) else None
                if fp and fp.exists() and fn:
                    # Rebuild series so tumour categorical colours match light/dark palette.
                    _update_charts_from_file(tname, fn)
                else:
                    pal = _echart_surface_palette()
                    _apply_palette_to_bar_chart(ch["bar"], pal)
                    _apply_palette_to_ts_chart(ch["ts"], pal)
                    ck = tool_to_step_map.get(tname, "sturgeon")
                    thr = _confidence_thresholds_for_classifier(ck)
                    ch["bar"].options["series"][0]["markLine"] = _bar_chart_mark_line(
                        thr["medium"], thr["high"], pal
                    )
                    ch["bar"].update()
                    ch["ts"].update()
            except Exception:
                pass

    from robin.gui.theme import register_theme_sync_callback

    unregister_classification_theme_sync = register_theme_sync_callback(
        _sync_classification_theme,
        element=next(iter(charts.values()))["bar"] if charts else None,
        interval_s=1.0,
        immediate=True,
    )

    def _read_scores_csv(csv_path: Path, mode: str):
        try:
            with csv_path.open("r", newline="") as fh:
                reader = csv.DictReader(fh)
                rows = list(reader)
                if not rows:
                    return None
                numeric_keys = []
                for k in reader.fieldnames or []:
                    if k.lower() in {"timestamp", "number_probes"}:
                        continue
                    try:
                        float(rows[-1].get(k, ""))
                        numeric_keys.append(k)
                    except Exception:
                        continue
                if not numeric_keys:
                    return None
                last_scores = {}
                for k in numeric_keys:
                    try:
                        v = float(rows[-1].get(k, 0))
                        if mode == "fraction":
                            v *= 100.0
                        last_scores[k] = round(v, 2)
                    except Exception:
                        pass
                time_key = (
                    "timestamp" if "timestamp" in (reader.fieldnames or []) else None
                )
                x_labels = []
                series_map: Dict[str, List[List[Any]]] = {k: [] for k in numeric_keys}
                for idx, r in enumerate(rows):
                    if time_key:
                        raw_x = r.get(time_key)
                        try:
                            x_val = float(raw_x)
                        except Exception:
                            x_val = raw_x
                    else:
                        x_val = idx + 1
                    x_labels.append(x_val)
                    for k in numeric_keys:
                        try:
                            vv = float(r.get(k, 0))
                            if mode == "fraction":
                                vv *= 100.0
                            series_map[k].append([x_val, round(vv, 2)])
                        except Exception:
                            series_map[k].append([x_val, 0])
                if time_key:

                    def _sort_key(p):
                        try:
                            return float(p[0])
                        except Exception:
                            return str(p[0])

                    for k in series_map:
                        series_map[k] = sorted(series_map[k], key=_sort_key)
                number_probes: Optional[int] = None
                try:
                    lr = rows[-1]
                    for nk in ("number_probes", "Number_probes"):
                        raw = lr.get(nk)
                        if raw is not None and str(raw).strip() != "":
                            number_probes = int(float(raw))
                            break
                except Exception:
                    number_probes = None
                return {
                    "last": last_scores,
                    "x": x_labels,
                    "series": series_map,
                    "has_time": bool(time_key),
                    "number_probes": number_probes,
                }
        except Exception:
            return None

    def _bar_values_to_rich_data(
        labels: List[str], values: List[float]
    ) -> List[Dict[str, Any]]:
        """Top-scoring bar (last index) = brand green + glow; others = distinct hues."""
        out: List[Dict[str, Any]] = []
        n = len(values)
        if len(labels) != n:
            labels = [""] * n
        for i, v in enumerate(values):
            is_top = n > 0 and i == n - 1
            if is_top:
                bar_color = _BRAND_GREEN
            else:
                lab = labels[i] if i < len(labels) else ""
                bar_color = _color_for_class_name(lab) if lab else _tumour_category_color(i)
            out.append(
                {
                    "value": v,
                    "itemStyle": {
                        "color": bar_color,
                        "borderRadius": _BAR_RADIUS,
                        **(
                            {
                                "shadowBlur": 8,
                                "shadowColor": "rgba(16,185,129,0.45)",
                            }
                            if is_top
                            else {}
                        ),
                    },
                }
            )
        return out

    def _update_charts_from_file(
        tool_name: str,
        file_name: str,
        preloaded: Optional[Dict[str, Any]] = None,
    ):
        try:
            file_path = sample_dir / file_name if sample_dir else None
            if not file_path or not file_path.exists():
                return
            mode = charts[tool_name].get("mode", "percent")
            data = preloaded if preloaded is not None else _read_scores_csv(file_path, mode)
            if not data:
                return
            bar = charts[tool_name]["bar"]
            ts = charts[tool_name]["ts"]
            classifier_key = tool_to_step_map.get(tool_name, "sturgeon")
            thr = _confidence_thresholds_for_classifier(classifier_key)
            medium_t = thr["medium"]
            high_t = thr["high"]
            pal = _echart_surface_palette()
            _z_slate = _ZONE_SLATE_DARK if _is_dark_mode() else _ZONE_SLATE
            _z_green = _ZONE_GREEN_DARK if _is_dark_mode() else _ZONE_GREEN

            last = data["last"]
            top = sorted(last.items(), key=lambda kv: kv[1], reverse=True)[:10]
            labels = [k for k, _ in top][::-1]
            values = [round(v, 2) for _, v in top][::-1]
            bar.options["yAxis"]["data"] = labels
            bar.options["series"][0]["data"] = _bar_values_to_rich_data(labels, values)
            bar.options["series"][0]["markLine"] = _bar_chart_mark_line(
                medium_t, high_t, pal
            )
            _apply_palette_to_bar_chart(bar, pal)
            bar.update()

            has_time = bool(data.get("has_time"))
            keys_top5 = [k for k, _ in top[:5]]
            leader = keys_top5[0] if keys_top5 else None
            secondaries = keys_top5[1:] if len(keys_top5) > 1 else []

            series_out: List[Dict[str, Any]] = [
                {
                    "name": "zone_medium_high",
                    "type": "line",
                    "data": [],
                    "silent": True,
                    "tooltip": {"show": False},
                    "symbol": "none",
                    "lineStyle": {"width": 0, "opacity": 0},
                    "markArea": {
                        "silent": True,
                        "itemStyle": {"color": _z_slate, "borderWidth": 0},
                        "data": [[{"yAxis": medium_t}, {"yAxis": high_t}]],
                    },
                    "z": 1,
                },
                {
                    "name": "zone_high_100",
                    "type": "line",
                    "data": [],
                    "silent": True,
                    "tooltip": {"show": False},
                    "symbol": "none",
                    "lineStyle": {"width": 0, "opacity": 0},
                    "markArea": {
                        "silent": True,
                        "itemStyle": {"color": _z_green, "borderWidth": 0},
                        "data": [[{"yAxis": high_t}, {"yAxis": 100}]],
                    },
                    "z": 1,
                },
                {
                    "name": "thr_lines",
                    "type": "line",
                    "data": [],
                    "silent": True,
                    "tooltip": {"show": False},
                    "symbol": "none",
                    "lineStyle": {"width": 0, "opacity": 0},
                    "markLine": {
                        "silent": True,
                        "symbol": "none",
                        "lineStyle": {
                            "color": pal["mark_line"],
                            "width": 1,
                            "type": "dashed",
                        },
                        "label": {"color": pal["mark_label"], "fontSize": 10},
                        "data": [
                            {
                                "yAxis": medium_t,
                                "label": {
                                    "formatter": f"Medium ({medium_t:.0f}%)",
                                    "show": True,
                                    "position": "end",
                                },
                            },
                            {
                                "yAxis": high_t,
                                "label": {
                                    "formatter": f"High ({high_t:.0f}%)",
                                    "show": True,
                                    "position": "end",
                                },
                            },
                        ],
                    },
                    "z": 2,
                },
            ]

            def _line_pts(k: str) -> Any:
                if has_time:
                    return data["series"][k]
                return [y for _, y in data["series"][k]]

            for k in secondaries:
                lc = _color_for_class_name(k)
                series_out.append(
                    {
                        "name": k,
                        "type": "line",
                        # Legend swatch + line must match (ECharts uses series.color for legend).
                        "color": lc,
                        "smooth": True,
                        "symbol": "none",
                        "animation": False,
                        "data": _line_pts(k),
                        "lineStyle": {
                            "width": 1.5,
                            "type": "dashed",
                            "color": lc,
                        },
                        "itemStyle": {"color": lc},
                        "z": 3,
                    }
                )

            if leader:
                series_out.append(
                    {
                        "name": leader,
                        "type": "line",
                        "color": _BRAND_GREEN_NEON,
                        "smooth": True,
                        "symbol": "none",
                        "animation": False,
                        "data": _line_pts(leader),
                        "lineStyle": {
                            "width": 3.5,
                            "color": _BRAND_GREEN_NEON,
                            "shadowBlur": 14,
                            "shadowColor": "rgba(34,197,94,0.55)",
                        },
                        "itemStyle": {"color": _BRAND_GREEN_NEON},
                        "z": 10,
                    }
                )

            ts.options["series"] = series_out
            ts.options["legend"]["data"] = keys_top5
            ts.options["tooltip"] = {
                "trigger": "axis",
                "backgroundColor": _TOOLTIP_DARK_BG,
                "borderColor": _TOOLTIP_DARK_BORDER,
                "borderWidth": 1,
                "textStyle": {"color": _TOOLTIP_DARK_TEXT},
                "axisPointer": {
                    "type": "line",
                    "lineStyle": {"color": pal["axis_pointer"]},
                },
            }

            if has_time:
                x_min = None
                x_max = None
                for k in keys_top5:
                    for p in data["series"][k]:
                        try:
                            xv = float(p[0])
                            x_min = xv if x_min is None else min(x_min, xv)
                            x_max = xv if x_max is None else max(x_max, xv)
                        except Exception:
                            pass
                xa = {
                    "type": "time",
                    "splitNumber": 4,
                    "axisLabel": {"hideOverlap": True, "color": pal["axis"]},
                }
                if x_min is not None and x_max is not None and x_max > x_min:
                    xa["min"] = x_min
                    xa["max"] = x_max
                ts.options["xAxis"] = xa
            else:
                ts.options["xAxis"] = {
                    "type": "category",
                    "data": data["x"],
                    "axisLabel": {"color": pal["axis"]},
                }

            _apply_palette_to_ts_chart(ts, pal)
            ts.update()
            # Summaries
            if charts[tool_name].get("summary"):
                labels_map = charts[tool_name]["summary"]
                if labels_map and top:
                    best_label, best_value = top[0]
                    if tool_name == "Random Forest":
                        labels_map["class"].set_text(
                            f"Forest classification: {best_label}"
                        )
                    else:
                        labels_map["class"].set_text(
                            f"{tool_name} classification: {best_label}"
                        )
                    labels_map["conf"].set_text(f"Confidence: {best_value:.2f}%")
                    try:
                        npv = data.get("number_probes")
                        if npv is not None:
                            label = (
                                "Features"
                                if tool_name == "Random Forest"
                                else "Probes"
                            )
                            labels_map["probes"].set_text(f"{label}: {int(npv)}")
                    except Exception:
                        pass
                try:
                    exp = charts[tool_name].get("expansion")
                    if exp and top:
                        best_label, best_value = top[0]
                        exp.props(
                            f'label="{tool_name} — {best_label} ({round(best_value, 2)}%)"'
                        )
                except Exception:
                    pass
        except Exception:
            pass

    def _refresh_classification_sync_impl() -> None:
        """Refresh classification data (synchronous; for contexts without asyncio loop)."""
        try:
            if not sample_dir or not sample_dir.exists():
                logging.warning(f"[Classification] Sample directory not found: {sample_dir}")
                return

            import time

            file_mtimes = {}
            files_changed = {}
            any_changes = False

            for tool_name, cfg in tool_to_file.items():
                file_path = sample_dir / cfg["file"] if sample_dir else None
                if file_path and file_path.exists():
                    mtime = file_path.stat().st_mtime
                    file_mtimes[tool_name] = mtime
                    prev_mtime = charts[tool_name].get("last_mtime", 0)
                    if prev_mtime is None or mtime > prev_mtime:
                        files_changed[tool_name] = True
                        any_changes = True
                    else:
                        files_changed[tool_name] = False
                else:
                    file_mtimes[tool_name] = 0
                    files_changed[tool_name] = False

            is_fresh_visit = any(
                charts[tool_name].get("last_mtime") is None
                for tool_name in tool_to_file.keys()
            )

            if not any_changes and not is_fresh_visit:
                logging.debug(
                    "[Classification] ⏭ Skipping classification update - no file changes detected"
                )
                return

            if any_changes or is_fresh_visit:
                reasons = []
                if is_fresh_visit:
                    reasons.append("fresh_visit")
                for tool_name, changed in files_changed.items():
                    if changed:
                        reasons.append(f"{tool_to_file[tool_name]['file']}")
                if reasons:
                    logging.debug(
                        f"[Classification] Update needed. Reasons: {', '.join(reasons)}"
                    )

            for tool_name, cfg in tool_to_file.items():
                file_path = sample_dir / cfg["file"] if sample_dir else None
                if file_path and files_changed.get(tool_name, False):
                    _check_and_update_file(tool_name, cfg["file"], file_path, charts)
                    if tool_name in file_mtimes:
                        charts[tool_name]["last_mtime"] = file_mtimes[tool_name]
        except Exception as e:
            logging.exception(f"[Classification] Refresh failed: {e}")

    async def _refresh_classification_async() -> None:
        """Parse score CSVs off the event loop, then update ECharts on the main thread."""
        try:
            if not sample_dir or not sample_dir.exists():
                logging.warning(f"[Classification] Sample directory not found: {sample_dir}")
                return

            import time

            file_mtimes = {}
            files_changed = {}
            any_changes = False

            for tool_name, cfg in tool_to_file.items():
                file_path = sample_dir / cfg["file"] if sample_dir else None
                if file_path and file_path.exists():
                    mtime = file_path.stat().st_mtime
                    file_mtimes[tool_name] = mtime
                    prev_mtime = charts[tool_name].get("last_mtime", 0)
                    if prev_mtime is None or mtime > prev_mtime:
                        files_changed[tool_name] = True
                        any_changes = True
                    else:
                        files_changed[tool_name] = False
                else:
                    file_mtimes[tool_name] = 0
                    files_changed[tool_name] = False

            is_fresh_visit = any(
                charts[tool_name].get("last_mtime") is None
                for tool_name in tool_to_file.keys()
            )

            if not any_changes and not is_fresh_visit:
                logging.debug(
                    "[Classification] ⏭ Skipping classification update - no file changes detected"
                )
                return

            if any_changes or is_fresh_visit:
                reasons = []
                if is_fresh_visit:
                    reasons.append("fresh_visit")
                for tool_name, changed in files_changed.items():
                    if changed:
                        reasons.append(f"{tool_to_file[tool_name]['file']}")
                if reasons:
                    logging.debug(
                        f"[Classification] Update needed. Reasons: {', '.join(reasons)}"
                    )

            for tool_name, cfg in tool_to_file.items():
                file_path = sample_dir / cfg["file"] if sample_dir else None
                if not file_path or not files_changed.get(tool_name, False):
                    continue
                mode = charts[tool_name].get("mode", "percent")
                data = await asyncio.to_thread(_read_scores_csv, file_path, mode)
                if data:
                    _update_charts_from_file(
                        tool_name, cfg["file"], preloaded=data
                    )
                if tool_name in file_mtimes:
                    charts[tool_name]["last_mtime"] = file_mtimes[tool_name]
        except Exception as e:
            logging.exception(f"[Classification] Refresh failed: {e}")

    def _refresh_classification() -> None:
        """Refresh classification data."""
        try:
            asyncio.get_running_loop()
        except RuntimeError:
            _refresh_classification_sync_impl()
            return
        asyncio.create_task(_refresh_classification_async())

    def _check_and_update_file(tool_name: str, filename: str, file_path: Path, charts: Dict[str, Any]) -> None:
        """Check file and update charts if needed.
        
        Note: This function is now called only when the file has changed,
        so we can skip the mtime check here (it's done upstream).
        """
        try:
            if file_path.exists():
                _update_charts_from_file(tool_name, filename)
        except Exception:
            pass

    # Start the refresh timer (every 30 seconds)
    refresh_timer = ui.timer(30.0, _refresh_classification, active=True, immediate=False)
    # Initial refresh after the page is rendered
    ui.timer(0.5, _refresh_classification, once=True)
    try:
        ui.context.client.on_disconnect(
            lambda: (refresh_timer.deactivate(), unregister_classification_theme_sync())
        )
    except Exception:
        pass
