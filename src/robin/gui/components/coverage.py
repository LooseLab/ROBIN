from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import json
import time
import os

import natsort
import numpy as np
import pandas as pd
import logging

try:
    from nicegui import ui, app
except ImportError:  # pragma: no cover
    ui = None
    app = None

from robin.gui.theme import styled_table, register_theme_sync_callback

# Shared paged table renderer to avoid materializing full row lists for large DataFrames.
def _render_paged_df_table(
    df: pd.DataFrame,
    columns: List[Dict[str, Any]],
    *,
    pagination_default: int = 100,
    class_size: str = "table-xs",
    search_placeholder: str = "Search...",
) -> Any:
    """Quasar server-side pagination over a bounded ``source_df`` slice (Python is source of truth)."""
    from robin.gui.theme import (
        clamp_qtable_server_pagination,
        styled_server_paged_table,
        wire_qtable_server_pagination_handlers,
    )

    preview_mode = len(df) > 50_000
    source_df = df.iloc[:5_000].copy() if preview_mode else df
    n_rows = len(source_df)
    page_state: Dict[str, Any] = {
        "filtered_positions": list(range(n_rows)),
    }

    total0 = len(page_state["filtered_positions"])
    pd_default = int(pagination_default)
    init_pagination = clamp_qtable_server_pagination(
        {
            "sortBy": None,
            "descending": False,
            "page": 1,
            "rowsPerPage": pd_default,
            "rowsNumber": total0,
        },
        rows_number=total0,
        rows_per_page_default=pd_default,
    )
    _, table = styled_server_paged_table(
        columns=columns,
        rows=[],
        pagination=init_pagination,
        row_key="__row_idx",
        class_size=class_size,
        rows_per_page_options=[25, 50, 100, 250],
    )

    count_label = ui.label("").classes("text-xs text-gray-500")
    if preview_mode:
        ui.label(
            f"Preview mode: showing first {len(source_df):,} rows of {len(df):,}."
        ).classes("text-xs text-amber-600")

    def _fill_from_pagination(pag: Dict[str, Any]) -> None:
        total = len(page_state["filtered_positions"])
        pag = clamp_qtable_server_pagination(
            pag,
            rows_number=total,
            rows_per_page_default=pd_default,
        )
        rpp = int(pag["rowsPerPage"])
        page = int(pag["page"])
        start = (page - 1) * rpp
        end = start + rpp
        positions = page_state["filtered_positions"]
        slice_pos = positions[start:end]
        page_rows = (
            source_df.iloc[slice_pos].reset_index(drop=True).to_dict("records")
            if slice_pos
            else []
        )
        for i, pos in enumerate(slice_pos):
            page_rows[i]["__row_idx"] = int(pos)
        table.rows = page_rows
        table.pagination = pag
        count_label.text = (
            f"{total} rows match search ({len(df)} total in source)"
            if total
            else f"0 rows match ({len(df)} total in source)"
        )
        table.update()

    wire_qtable_server_pagination_handlers(table, _fill_from_pagination)

    def _apply_search(term: str) -> None:
        txt = str(term or "").strip().lower()
        if not txt:
            page_state["filtered_positions"] = list(range(n_rows))
        else:
            filtered: List[int] = []
            for pos in range(n_rows):
                row = source_df.iloc[pos]
                for col in source_df.columns:
                    value = row.get(col)
                    if value is not None and txt in str(value).lower():
                        filtered.append(pos)
                        break
            page_state["filtered_positions"] = filtered
        pag = clamp_qtable_server_pagination(
            dict(table.pagination),
            rows_number=len(page_state["filtered_positions"]),
            rows_per_page_default=pd_default,
        )
        pag["page"] = 1
        _fill_from_pagination(pag)

    with table.add_slot("top-right"):
        search_input = ui.input(placeholder=search_placeholder).props(
            "type=search dense clearable"
        )
        search_input.on("update:model-value", lambda e: _apply_search(getattr(e, "value", "")))

    for col in table.columns:
        col["sortable"] = False

    _fill_from_pagination(init_pagination)
    try:
        def _cleanup_table_state() -> None:
            page_state["filtered_positions"] = []
            table.rows = []

        ui.context.client.on_disconnect(_cleanup_table_state)
    except Exception:
        pass
    return table

# --- Coverage charts (design.md §9.5): on/off colours, mean line, gene palette ---
_COV_GENE_LINE_PALETTE = [
    "#06b6d4",  # cyan
    "#c026d3",  # magenta / fuchsia
    "#ca8a04",  # gold / amber
    "#2563eb",  # blue
    "#db2777",  # rose
]
_COV_MAX_DIM_GENES = 25


def _cov_dim_muted_line() -> str:
    """Dim in-range gene traces: readable on light and dark plot surfaces."""
    if _cov_ui_dark():
        return "#94a3b8"
    return "#64748b"


def _cov_band_fill_rgba() -> str:
    """Mean ±2σ band fill behind target coverage over time."""
    if _cov_ui_dark():
        return "rgba(148,163,184,0.18)"
    return "rgba(148,163,184,0.12)"


def _cov_tooltip_option() -> Dict[str, Any]:
    """Tooltip chrome for coverage ECharts (JSON-serialisable)."""
    p = _cov_chrome_palette()
    return {
        "backgroundColor": p["tooltip_bg"],
        "borderColor": p["tooltip_border"],
        "borderWidth": 1,
        "textStyle": {"color": p["tooltip_text"]},
    }


def _cov_boxplot_item_style() -> Dict[str, Any]:
    """Box + whisker stroke/fill for chromosome overview (readable in light/dark)."""
    if _cov_ui_dark():
        return {
            "color": "rgba(148, 163, 184, 0.22)",
            "borderColor": "#cbd5e1",
            "borderWidth": 1.2,
        }
    return {
        "color": "rgba(71, 85, 105, 0.12)",
        "borderColor": "#334155",
        "borderWidth": 1.2,
    }


def _cov_scatter_outlier_style() -> Dict[str, Any]:
    """Outlier points on boxplot: stroke matches surface."""
    if _cov_ui_dark():
        return {"color": "#38bdf8", "borderColor": "#0f172a", "borderWidth": 1}
    return {"color": "#0284c7", "borderColor": "#ffffff", "borderWidth": 1}


def _cov_ui_dark() -> bool:
    try:
        return bool(app.storage.user.get("dark_mode"))
    except Exception:
        return False


def _cov_on_off_colors() -> Tuple[str, str]:
    """On-target (brand green family) vs off-target (muted) for chromosome bars."""
    if _cov_ui_dark():
        return ("#22c55e", "#475569")
    # Light surfaces: off-target must stay visible on white / --md-surface (not slate-300).
    return ("#10b981", "#64748b")


def _cov_mean_line_color() -> str:
    """Thick dashed mean line: high contrast on plot surface."""
    return "#f8fafc" if _cov_ui_dark() else "#1e293b"


def _cov_chrome_palette() -> Dict[str, str]:
    """Titles / axes / legend / tooltip for coverage ECharts (matches classification charts)."""
    if _cov_ui_dark():
        return {
            "title": "#f1f5f9",
            "axis": "#cbd5e1",
            "split": "rgba(148, 163, 184, 0.35)",
            "legend": "#cbd5e1",
            "tooltip_bg": "rgba(15, 23, 42, 0.96)",
            "tooltip_border": "#475569",
            "tooltip_text": "#f8fafc",
        }
    return {
        "title": "#0f172a",
        "axis": "#1e293b",
        "split": "rgba(51, 65, 85, 0.45)",
        "legend": "#1e293b",
        "tooltip_bg": "rgba(255, 255, 255, 0.96)",
        "tooltip_border": "#e2e8f0",
        "tooltip_text": "#0f172a",
    }


def _cov_area_fill_primary() -> str:
    """Semi-transparent fill under cumulative coverage line."""
    if _cov_ui_dark():
        return "rgba(34, 197, 94, 0.22)"
    return "rgba(16, 185, 129, 0.18)"


def _cov_line_primary() -> str:
    if _cov_ui_dark():
        return "#22c55e"
    return "#10b981"


def _short_chr_label(name: str) -> str:
    """Minimalist chromosome labels: chr1 -> 1, chrX -> X."""
    t = str(name).strip()
    if t.lower().startswith("chr"):
        t = t[3:]
    return t


def _cov_gene_line_color(idx: int) -> str:
    return _COV_GENE_LINE_PALETTE[idx % len(_COV_GENE_LINE_PALETTE)]


def _cov_legend_label(text: str, max_len: int = 44) -> str:
    """Short label for ECharts legend (JSON-safe options cannot use JS formatters)."""
    t = (text or "").replace("\n", " ").strip()
    if len(t) <= max_len:
        return t
    return t[: max_len - 1] + "…"


def _apply_coverage_target_cov_chrome(ec: Any) -> None:
    """Re-apply axis/title/legend/tooltip when light/dark toggles."""
    try:
        p = _cov_chrome_palette()
        tt = _cov_tooltip_option()
        o = ec.options
        o["textStyle"] = {"color": p["axis"]}
        o["title"]["textStyle"] = {"fontSize": 15, "color": p["title"]}
        o["legend"]["textStyle"] = {"color": p["legend"]}
        o["tooltip"] = {
            "trigger": "axis",
            "axisPointer": {"type": "shadow"},
            **tt,
        }
        o["xAxis"]["axisLabel"]["color"] = p["axis"]
        o["xAxis"]["axisLine"]["lineStyle"]["color"] = p["axis"]
        o["yAxis"]["nameTextStyle"]["color"] = p["axis"]
        o["yAxis"]["axisLabel"]["color"] = p["axis"]
        o["yAxis"]["splitLine"]["lineStyle"]["color"] = p["split"]
    except Exception:
        pass


def _apply_coverage_time_series_chrome(ec: Any) -> None:
    """Cumulative coverage over time — chrome + line/area from theme."""
    try:
        p = _cov_chrome_palette()
        tt = _cov_tooltip_option()
        ln = _cov_line_primary()
        fill = _cov_area_fill_primary()
        o = ec.options
        o["textStyle"] = {"color": p["axis"]}
        o["title"]["textStyle"] = {"fontSize": 15, "color": p["title"]}
        o["title"]["subtextStyle"] = {"fontSize": 11, "color": p["axis"]}
        o["tooltip"] = {
            "trigger": "axis",
            "axisPointer": {
                "type": "cross",
                "lineStyle": {"color": p["axis"], "width": 1},
                "crossStyle": {"color": p["axis"]},
            },
            **tt,
        }
        o["xAxis"]["axisLabel"] = {"color": p["axis"]}
        o["xAxis"]["axisLine"] = {"lineStyle": {"color": p["axis"]}}
        o["xAxis"]["splitLine"] = {"lineStyle": {"color": p["split"]}}
        o["yAxis"]["nameTextStyle"] = {"color": p["axis"]}
        o["yAxis"]["axisLabel"] = {"color": p["axis"]}
        o["yAxis"]["splitLine"] = {
            "lineStyle": {"color": p["split"], "type": "dashed"},
        }
        if o.get("series") and len(o["series"]) > 0:
            o["series"][0]["lineStyle"] = {"width": 2.5, "color": ln}
            o["series"][0]["itemStyle"] = {"color": ln}
            o["series"][0]["areaStyle"] = {"color": fill}
    except Exception:
        pass


def _apply_coverage_boxplot_chrome(bp: Any) -> None:
    """Target coverage boxplot + scatter outliers — readable in light/dark."""
    try:
        p = _cov_chrome_palette()
        tt = _cov_tooltip_option()
        box_st = _cov_boxplot_item_style()
        sc_st = _cov_scatter_outlier_style()
        o = bp.options
        o["textStyle"] = {"color": p["axis"]}
        o.setdefault("title", {})
        o["title"]["textStyle"] = {"fontSize": 16, "color": p["title"]}
        o["title"].setdefault("subtextStyle", {})
        o["title"]["subtextStyle"] = {"fontSize": 12, "color": p["axis"]}
        o["tooltip"] = {"trigger": "item", **tt}
        o.setdefault("legend", {})
        o["legend"]["textStyle"] = {"color": p["legend"]}
        o.setdefault("xAxis", {})
        o["xAxis"]["name"] = o["xAxis"].get("name", "Chromosome")
        # Center the axis title — default "end" puts the label at the right and it clips.
        o["xAxis"]["nameLocation"] = "middle"
        o["xAxis"]["nameGap"] = int(o["xAxis"].get("nameGap") or 40)
        o["xAxis"]["nameTextStyle"] = {"color": p["axis"]}
        o["xAxis"]["axisLabel"] = {
            **(o["xAxis"].get("axisLabel") or {}),
            "color": p["axis"],
        }
        o["xAxis"]["axisLine"] = {"lineStyle": {"color": p["axis"]}}
        o.setdefault("yAxis", {})
        o["yAxis"]["name"] = o["yAxis"].get("name", "Coverage (x)")
        o["yAxis"]["nameTextStyle"] = {"color": p["axis"]}
        o["yAxis"]["axisLabel"] = {"color": p["axis"]}
        o["yAxis"]["splitLine"] = {
            "show": True,
            "lineStyle": {"color": p["split"], "type": "dashed"},
        }
        for s in o.get("series") or []:
            t = s.get("type")
            if t == "boxplot":
                s["itemStyle"] = box_st
            elif t == "scatter":
                nm = str(s.get("name") or "")
                if nm in ("outliers", "global outliers"):
                    s["itemStyle"] = {**sc_st, **(s.get("itemStyle") or {})}
                    if s.get("label", {}).get("show"):
                        lab = dict(s.get("label") or {})
                        lab["color"] = p["axis"]
                        s["label"] = lab
                elif nm == "Gene Coverage":
                    ex = s.get("itemStyle") or {}
                    s["itemStyle"] = {
                        "color": ex.get("color", "#3b82f6"),
                        "opacity": ex.get("opacity", 0.7),
                    }
                    em = s.get("emphasis") or {}
                    eis = dict(em.get("itemStyle") or {})
                    eis["borderColor"] = p["title"]
                    eis["borderWidth"] = 2
                    s["emphasis"] = {**em, "itemStyle": eis}
                    if s.get("label", {}).get("show"):
                        lab = dict(s.get("label") or {})
                        lab["color"] = p["axis"]
                        s["label"] = lab
    except Exception:
        pass


def _apply_target_coverage_time_analysis_chrome(ec: Any) -> None:
    """Outlier / mean / band chart — sync chrome and theme-dependent series colours."""
    try:
        p = _cov_chrome_palette()
        tt = _cov_tooltip_option()
        mean_c = _cov_mean_line_color()
        band_fill = _cov_band_fill_rgba()
        dim_c = _cov_dim_muted_line()
        o = ec.options
        o["textStyle"] = {"color": p["axis"]}
        o.setdefault("title", {})
        o["title"]["textStyle"] = {"fontSize": 15, "color": p["title"]}
        o["title"]["subtextStyle"] = {"fontSize": 11, "color": p["axis"]}
        o["tooltip"] = {
            "trigger": "axis",
            "axisPointer": {
                "type": "cross",
                "lineStyle": {"color": p["axis"], "width": 1},
                "crossStyle": {"color": p["axis"]},
            },
            "confine": True,
            "appendToBody": True,
            **tt,
        }
        o.setdefault("legend", {})
        o["legend"]["pageIconColor"] = p["title"]
        o["legend"]["pageTextStyle"] = {"color": p["axis"]}
        o["legend"]["textStyle"] = {"fontSize": 11, "color": p["legend"]}
        o.setdefault("xAxis", {})
        o["xAxis"]["axisLabel"] = {"color": p["axis"]}
        o["xAxis"]["axisLine"] = {"lineStyle": {"color": p["axis"]}}
        o["xAxis"]["splitLine"] = {"lineStyle": {"color": p["split"]}}
        o.setdefault("yAxis", {})
        o["yAxis"]["axisLabel"] = {"color": p["axis"]}
        o["yAxis"]["axisLine"] = {"lineStyle": {"color": p["axis"]}}
        o["yAxis"]["splitLine"] = {
            "lineStyle": {"color": p["split"], "type": "dashed"},
        }
        o["yAxis"]["nameTextStyle"] = {"color": p["axis"]}
        for s in o.get("series") or []:
            name = str(s.get("name") or "")
            if name == "_band":
                ma = s.get("markArea") or {}
                ist = ma.get("itemStyle") or {}
                ist["color"] = band_fill
                ma["itemStyle"] = ist
                s["markArea"] = ma
            elif name == "Mean (population)":
                s["lineStyle"] = {
                    **(s.get("lineStyle") or {}),
                    "width": 4,
                    "type": "dashed",
                    "color": mean_c,
                }
                s["itemStyle"] = {**(s.get("itemStyle") or {}), "color": mean_c}
            elif "in-range" in name:
                s["itemStyle"] = {
                    **(s.get("itemStyle") or {}),
                    "color": dim_c,
                    "opacity": 0.42,
                }
                s["lineStyle"] = {
                    **(s.get("lineStyle") or {}),
                    "width": 1.5,
                    "type": "dashed",
                    "color": dim_c,
                    "opacity": 0.42,
                }
    except Exception:
        pass


def _apply_lga_gene_bar_chrome(ec: Any) -> None:
    """Top genes by coverage bar chart (LGA overview)."""
    try:
        p = _cov_chrome_palette()
        tt = _cov_tooltip_option()
        bar_fill = "#60a5fa" if _cov_ui_dark() else "#3b82f6"
        shadow = "rgba(0,0,0,0.45)" if _cov_ui_dark() else "rgba(0,0,0,0.3)"
        o = ec.options
        o["textStyle"] = {"color": p["axis"]}
        o["title"]["textStyle"] = {
            "fontSize": 16,
            "color": p["title"],
            "fontWeight": "bold",
        }
        o["title"]["subtextStyle"] = {"fontSize": 12, "color": p["axis"]}
        o["tooltip"] = {
            "trigger": "axis",
            **tt,
        }
        o["xAxis"]["axisLabel"] = {
            **(o["xAxis"].get("axisLabel") or {}),
            "color": p["axis"],
        }
        o["xAxis"]["nameTextStyle"] = {"color": p["axis"]}
        o["xAxis"]["axisLine"] = {"lineStyle": {"color": p["axis"]}}
        o["yAxis"]["nameTextStyle"] = {"color": p["axis"]}
        o["yAxis"]["axisLabel"] = {
            **(o["yAxis"].get("axisLabel") or {}),
            "color": p["axis"],
        }
        o["yAxis"]["splitLine"] = {
            "show": True,
            "lineStyle": {"color": p["split"], "type": "dashed"},
        }
        if o.get("dataZoom"):
            for dz in o["dataZoom"]:
                if isinstance(dz, dict):
                    dz["textStyle"] = {"color": p["axis"]}
                    dz["borderColor"] = p["tooltip_border"]
                    dz["handleStyle"] = {"color": p["axis"]}
        if o.get("series") and len(o["series"]) > 0:
            o["series"][0]["itemStyle"] = {"color": bar_fill}
            o["series"][0]["emphasis"] = {
                "itemStyle": {
                    "shadowBlur": 10,
                    "shadowColor": shadow,
                }
            }
    except Exception:
        pass


def add_igv_viewer(launcher: Any, sample_dir: Path) -> None:
    """Add the IGV viewer section to a page.
    
    Args:
        launcher: The GUI launcher instance for state management
        sample_dir: Path to the sample directory
    """
    # IGV viewer section - only show if target.bam exists
    target_bam = sample_dir / "target.bam"
    if not (target_bam.exists() and target_bam.is_file()):
        with ui.element("div").classes(
            "classification-insight-shell w-full min-w-0"
        ).props("id=sample-details-igv"):
            ui.label("IGV browser").classes(
                "classification-insight-heading text-headline-small"
            )
            with ui.element("div").classes("classification-insight-card w-full min-w-0"):
                with ui.column().classes("w-full min-w-0 gap-2 p-2 md:p-3"):
                    with ui.row().classes("items-center gap-2 min-w-0"):
                        ui.icon("biotech").classes("classification-insight-icon")
                        ui.label("target.bam required").classes(
                            "classification-insight-model flex-1 min-w-0"
                        )
                    ui.label(
                        "IGV requires target.bam. Run target analysis for this sample first."
                    ).classes("classification-insight-foot")
        return

    with ui.element("div").classes(
        "classification-insight-shell w-full min-w-0"
    ).props("id=sample-details-igv"):
        ui.label("IGV browser").classes(
            "classification-insight-heading text-headline-small"
        )
        ui.label("Genome: hg38 · Interactive alignment view").classes(
            "classification-insight-meta mb-2"
        )
        igv_div = ui.element("div").classes(
            "w-full h-[900px] border border-slate-200/90 rounded-lg overflow-hidden "
            "dark:border-slate-700/80"
        )
        igv_div._props["id"] = "igv-container"
        igv_status = ui.label("Initializing IGV viewer…").classes(
            "classification-insight-meta"
        )

        igv_lib_status = ui.label("IGV library: checking…").classes(
            "classification-insight-foot"
        )
        
        # Initialize IGV browser immediately on page load
        def _initialize_igv_browser():
            """Initialize the IGV browser immediately."""
            js_init = """
                (function() {
                    try {
                        console.log('[IGV] Starting browser initialization...');
                        const el = document.getElementById('igv-container');
                        if (!el) {
                            console.error('[IGV] Container element not found');
                            return 'missing-element';
                        }

                        if (window.lj_igv && window.lj_igv_browser_ready) {
                            console.log('[IGV] Browser already exists');
                            return 'already-ready';
                        }

                        if (window.lj_igv_initializing) {
                            console.log('[IGV] Initialization already in progress');
                            return 'initializing';
                        }

                        // Clear any existing content
                        el.innerHTML = '';

                        window.lj_igv_initializing = true;
                        window.lj_igv_browser_ready = false;

                        console.log('[IGV] Creating browser with hg38 genome...');
                        const options = {
                            genome: 'hg38',
                            showRuler: true,
                            showNavigation: true,
                            showCenterGuide: true,
                            defaultLocus: 'chr1:1-500000'
                        };

                        igv.createBrowser(el, options)
                            .then(function(browser) {
                                console.log('[IGV] Browser created successfully');
                                window.lj_igv = browser;
                                window.lj_igv_browser_ready = true;
                                window.lj_igv_initializing = false;

                                const statusEl = document.querySelector('[data-igv-status]');
                                if (statusEl) {
                                    statusEl.textContent = 'IGV browser ready';
                                }

                                console.log('[IGV] Browser ready, target BED will be loaded with BAM');
                            })
                            .catch(function(error) {
                                console.error('[IGV] Browser creation failed:', error);
                                window.lj_igv_browser_ready = false;
                                window.lj_igv_initializing = false;

                                const statusEl = document.querySelector('[data-igv-status]');
                                if (statusEl) {
                                    statusEl.textContent = 'IGV browser creation failed: ' + (error && error.message ? error.message : error);
                                }
                            });

                        return 'started';
                    } catch (error) {
                        console.error('[IGV] Initialization error:', error);
                        window.lj_igv_browser_ready = false;
                        window.lj_igv_initializing = false;
                        return 'error';
                    }
                })();
            """
            try:
                init_result = ui.run_javascript(js_init, timeout=5.0)
                if init_result == "missing-element":
                    igv_status.set_text("IGV container element not found.")
                elif init_result == "error":
                    igv_status.set_text("IGV initialization failed (see console).")
                elif init_result == "started":
                    igv_status.set_text("IGV browser initialization started...")
                elif init_result == "initializing":
                    igv_status.set_text("IGV browser is still initializing...")
                elif init_result == "already-ready":
                    igv_status.set_text("IGV browser already ready.")
            except Exception as e:
                print(f"Error initializing IGV browser: {e}")
                igv_status.set_text(f"Error initializing IGV: {e}")
        
        # Add data attribute for status updates
        igv_status._props["data-igv-status"] = True
        
        # Initialize browser after a short delay to ensure DOM is ready
        ui.timer(0.5, _initialize_igv_browser, once=True)

        # IGV state management - prevent unnecessary redraws
        def _get_igv_state():
            try:
                key = str(sample_dir)
                if hasattr(launcher, "_coverage_state"):
                    if key not in launcher._coverage_state:
                        launcher._coverage_state[key] = {}
                    return launcher._coverage_state[key]
                else:
                    launcher._coverage_state = {}
                    launcher._coverage_state[key] = {}
                    return launcher._coverage_state[key]
            except Exception:
                # Return empty state as fallback
                return {}

        def _is_igv_ready():
            """Check if IGV browser is properly initialized and ready."""
            try:
                state = _get_igv_state()
                return (
                    state.get("igv_initialized", False)
                    and state.get("igv_browser_ready", False)
                    and state.get("igv_loaded_bam") is not None
                )
            except Exception:
                return False

        def _set_igv_ready(bam_url: str):
            """Mark IGV as ready with the current BAM."""
            state = _get_igv_state()
            state["igv_initialized"] = True
            state["igv_browser_ready"] = True
            state["igv_loaded_bam"] = bam_url

        def _clear_igv_state():
            """Reset IGV state to force reinitialization."""
            state = _get_igv_state()
            state["igv_initialized"] = False
            state["igv_browser_ready"] = False
            state["igv_loaded_bam"] = None
            state["igv_loading"] = False

        # Check for existing IGV BAM files and load them
        def _check_existing_igv_bam():
            try:
                # Wait a bit for browser to be ready
                js_check_ready = """
                    return window.lj_igv && window.lj_igv_browser_ready === true && 
                           typeof window.lj_igv.loadTrack === 'function';
                """
                
                # Poll for browser readiness
                max_attempts = 10
                for attempt in range(max_attempts):
                    try:
                        result = ui.run_javascript(js_check_ready, timeout=2.0)
                        if result is True:
                            break
                    except Exception:
                        pass
                    time.sleep(0.5)
                
                # Determine candidate BAMs (prefer standardized IGV path)
                igv_dir = sample_dir / "igv"
                clair_dir = sample_dir / "clair3"
                candidates = [
                    igv_dir / "igv_ready.bam",
                    clair_dir / "sorted_targets_exceeding_rerun.bam",
                    clair_dir / "sorted_targets_exceeding.bam",
                    sample_dir / "target.bam",
                ]

                bam_path = None
                for p in candidates:
                    bai = Path(f"{p}.bai")
                    if p.exists() and bai.exists():
                        bam_path = p
                        break

                if bam_path is not None:
                    igv_status.set_text(f"Found BAM: {bam_path.name}, loading...")
                    # Load the BAM file
                    _load_bam_track_simple(bam_path)
                else:
                    igv_status.set_text(
                        "IGV browser ready. No BAM files found yet."
                    )

            except Exception as e:
                igv_status.set_text(f"Error checking for BAM files: {e}")
        
        # Simple function to load BAM track
        def _load_bam_track_simple(bam_path: Path):
            """Simple function to load a BAM track into IGV."""
            try:
                # Mount directory
                if bam_path.parent == sample_dir / "igv":
                    mount = f"/samples/{sample_dir.name}/igv"
                elif bam_path.parent == sample_dir / "clair3":
                    mount = f"/samples/{sample_dir.name}/clair3"
                else:
                    mount = f"/samples/{sample_dir.name}"

                try:
                    if app is not None:
                        app.add_static_files(mount, str(bam_path.parent))
                except Exception:
                    pass

                bam_url = f"{mount}/{bam_path.name}"
                bai_url = f"{mount}/{bam_path.name}.bai"
                
                js_load = f"""
                    (function loadBamTrack() {{
                        try {{
                            if (!window.lj_igv || !window.lj_igv_browser_ready) {{
                                console.log('[IGV] Browser not ready yet, retrying in 500ms...');
                                setTimeout(loadBamTrack, 500);
                                return;
                            }}
                            
                            // Verify browser has loadTrack method
                            if (typeof window.lj_igv.loadTrack !== 'function') {{
                                console.error('[IGV] Browser does not have loadTrack method');
                                const statusEl = document.querySelector('[data-igv-status]');
                                if (statusEl) {{
                                    statusEl.textContent = 'Error: Browser not fully initialized';
                                }}
                                return;
                            }}
                            
                            console.log('[IGV] Loading BAM track:', '{bam_path.name}');
                            console.log('[IGV] BAM URL:', '{bam_url}');
                            console.log('[IGV] BAI URL:', '{bai_url}');
                            
                            const trackConfig = {{
                                name: '{sample_dir.name}',
                                url: '{bam_url}',
                                indexURL: '{bai_url}',
                                format: 'bam',
                                type: 'alignment',
                                order: Number.MAX_VALUE,
                                visibilityWindow: 500000,
                                height: 600,
                                autoScale: true,
                                colorBy: 'tag',
                                tag: 'SA'
                            }};
                            
                            console.log('[IGV] Track config:', trackConfig);
                            
                            window.lj_igv.loadTrack(trackConfig)
                                .then(function(trackView) {{
                                    console.log('[IGV] Track loaded successfully:', trackView);
                                    
                                    // Navigate to a region to make data visible (IGV needs to be zoomed in)
                                    // Navigate to chr1:1-500000 to match visibility window
                                    try {{
                                        window.lj_igv.search('chr1:1-500000');
                                        console.log('[IGV] Navigated to chr1:1-500000');
                                    }} catch (navError) {{
                                        console.log('[IGV] Could not navigate, user can zoom manually');
                                    }}
                                    
                                    const statusEl = document.querySelector('[data-igv-status]');
                                    if (statusEl) {{
                                        statusEl.textContent = 'BAM loaded: {bam_path.name} - Zoom in to see alignments';
                                    }}
                                }})
                                .catch(function(error) {{
                                    console.error('[IGV] Track load failed:', error);
                                    console.error('[IGV] Error details:', error.message, error.stack);
                                    
                                    const statusEl = document.querySelector('[data-igv-status]');
                                    if (statusEl) {{
                                        statusEl.textContent = 'Failed to load BAM: ' + (error.message || String(error));
                                    }}
                                }});
                        }} catch (error) {{
                            console.error('[IGV] Error in loadBamTrack function:', error);
                            const statusEl = document.querySelector('[data-igv-status]');
                            if (statusEl) {{
                                statusEl.textContent = 'Error: ' + error.message;
                            }}
                        }}
                    }})();
                """
                ui.run_javascript(js_load, timeout=60.0)
                
            except Exception as e:
                igv_status.set_text(f"Error loading BAM: {e}")
                print(f"Error in _load_bam_track_simple: {e}")

        # Function to check if IGV library is available
        def _check_igv_library():
            """Check if the IGV JavaScript library is available."""
            try:
                js_check = """
                    try {
                        if (typeof igv === 'undefined') {
                            console.error('IGV library not loaded');
                            return false;
                        }
                
                        // Check if IGV is fully initialized
                        if (typeof igv.createBrowser !== 'function') {
                            console.error('IGV library not fully initialized');
                            return false;
                        }
                
                        console.log('IGV library is available and ready');
                        return true;
                    } catch (e) {
                        console.error('Error checking IGV library:', e);
                        return false;
                    }
                """

                try:
                    result = ui.run_javascript(js_check, timeout=10.0)
                    # Update the status indicator
                    if result:
                        igv_lib_status.set_text("IGV library: ✓ Loaded and ready")
                        igv_lib_status.classes("text-xs text-green-600")
                    else:
                        igv_lib_status.set_text("IGV library: ✗ Not ready")
                        igv_lib_status.classes("text-xs text-red-600")
                    return result
                except Exception as e:
                    igv_lib_status.set_text("IGV library: ✗ Error checking")
                    igv_lib_status.classes("text-xs text-red-600")
                    print(f"Error checking IGV library: {e}")
                    return False
            except Exception:
                return False

        # Function to retry IGV creation
        def _retry_igv_creation(bam_path: Path):
            """Retry IGV browser creation after a failure."""
            try:
                igv_status.set_text("Retrying IGV browser creation...")
                _load_igv_bam(bam_path)
            except Exception as e:
                igv_status.set_text(f"Retry failed: {e}")
                print(f"IGV retry error: {e}")

        # Function to wait for element to be ready
        def _wait_for_element_ready():
            """Wait for the IGV div element to be properly rendered."""
            js_wait = """
                return new Promise((resolve, reject) => {
                    const el = document.getElementById('igv-container');
                    if (!el) {
                        reject(new Error('Element not found'));
                        return;
                    }
            
                    if (el.offsetWidth > 0 && el.offsetHeight > 0) {
                        resolve(true);
                        return;
                    }
            
                    // Wait for element to be ready
                    const checkReady = () => {
                        if (el.offsetWidth > 0 && el.offsetHeight > 0) {
                            resolve(true);
                        } else {
                            setTimeout(checkReady, 100);
                        }
                    };
                    checkReady();
                });
            """
            try:
                return ui.run_javascript(js_wait, timeout=10.0)
            except Exception:
                return False

        # Function to load BAM into IGV
        def _load_igv_bam(bam_path: Path):
            try:
                # Get state first
                state = _get_igv_state()

                # Prevent multiple simultaneous IGV loading attempts
                if state.get("igv_loading", False):
                    igv_status.set_text("IGV is already being loaded, please wait...")
                    return

                # First check if IGV library is available
                if not _check_igv_library():
                    igv_status.set_text(
                        "IGV library not yet loaded. Please wait and try again."
                    )
                    return

                # Check if we already have this BAM loaded
                bam_url = f"/samples/{sample_dir.name}/{bam_path.name}"
                if state.get("igv_loaded_bam") == bam_url and _is_igv_ready():
                    igv_status.set_text(f"BAM {bam_path.name} already loaded in IGV.")
                    return

                # Mark that we're loading IGV
                state["igv_loading"] = True

                # Mount directory once
                if bam_path.parent == sample_dir / "igv":
                    mount = f"/samples/{sample_dir.name}/igv"
                elif bam_path.parent == sample_dir / "clair3":
                    mount = f"/samples/{sample_dir.name}/clair3"
                else:
                    mount = f"/samples/{sample_dir.name}"

                try:
                    if app is not None:
                        app.add_static_files(mount, str(bam_path.parent))
                except Exception:
                    # Already mounted or in tests without app
                    pass

                # Define URLs
                bam_url = f"{mount}/{bam_path.name}"
                bai_url = f"{mount}/{bam_path.name}.bai"

                # Check if we need to create a new browser or just add tracks
                if not state.get("igv_initialized") or not state.get(
                    "igv_browser_ready"
                ):
                    # Wait for element to be ready before creating IGV browser
                    if not _wait_for_element_ready():
                        igv_status.set_text("Waiting for IGV element to be ready...")
                        ui.timer(0.5, lambda: _load_igv_bam(bam_path), once=True)
                        return

                    # Initialize JavaScript flag before creating browser
                    js_init_flag = """
                        window.lj_igv_browser_ready = false;
                        console.log('Initialized IGV browser ready flag to false');
                    """
                    ui.run_javascript(js_init_flag, timeout=5.0)

                    # Check if browser already exists in JavaScript
                    js_check_existing = """
                        if (window.lj_igv && window.lj_igv_browser_ready === true) {
                            console.log('IGV browser already exists, skipping creation');
                            return true;
                        }
                        const el = document.getElementById('igv-container');
                        if (el && el.children.length > 0) {
                            console.log('IGV container already has children, browser may exist');
                            return true;
                        }
                        return false;
                    """
                    try:
                        existing = ui.run_javascript(js_check_existing, timeout=5.0)
                        if existing:
                            # Browser might already exist, check Python state after a delay
                            time.sleep(0.5)
                            if state.get("igv_initialized") and state.get("igv_browser_ready"):
                                # Browser is ready, just add the track
                                pass  # Will fall through to track loading
                            else:
                                # Wait a bit more and check again
                                ui.timer(1.0, lambda: _load_igv_bam(bam_path), once=True)
                                return
                    except Exception:
                        pass

                    # Create new IGV browser
                    js_create = """
                        try {
                            console.log('[DEBUG] Creating IGV browser...');
                            const el = document.getElementById('igv-container');
                            if (!el) {
                                throw new Error('Target element not found');
                            }
                            
                            // Check again if browser already exists
                            if (window.lj_igv && window.lj_igv_browser_ready === true) {
                                console.log('[DEBUG] Browser already exists, skipping creation');
                                return;
                            }
                            
                            const options = { genome: 'hg38' };
                            console.log('[DEBUG] IGV options:', options);
                            
                            window.lj_igv_browser_ready = false;
                            igv.createBrowser(el, options).then(b => {
                                window.lj_igv = b; 
                                window.lj_igv_browser_ready = true;
                                console.log('[DEBUG] IGV browser created successfully');
                            }).catch(error => {
                                console.error('[DEBUG] IGV browser creation failed:', error);
                                window.lj_igv_browser_ready = false;
                            });
                        } catch (error) {
                            console.error('[DEBUG] Error in IGV browser creation:', error);
                            window.lj_igv_browser_ready = false;
                        }
                    """

                    try:
                        ui.run_javascript(js_create, timeout=30.0)
                        igv_status.set_text("Creating IGV browser...")
                        
                        # Wait a bit for the async browser creation to start
                        time.sleep(0.5)
                        
                        # Poll for browser readiness
                        max_poll_attempts = 20
                        poll_interval = 1.0
                        browser_ready = False
                        
                        for attempt in range(max_poll_attempts):
                            js_check_ready = """
                                return window.lj_igv && window.lj_igv_browser_ready === true && 
                                       typeof window.lj_igv.loadTrack === 'function';
                            """
                            try:
                                result = ui.run_javascript(js_check_ready, timeout=5.0)
                                if result is True:
                                    browser_ready = True
                                    break
                            except Exception:
                                pass
                            
                            time.sleep(poll_interval)
                        
                        if browser_ready:
                            _set_igv_ready(bam_url)
                            igv_status.set_text(f"IGV browser ready with {bam_path.name}")
                            state["igv_loading"] = False
                            
                            # Now load the track
                            js_add_track = f"""
                                if (window.lj_igv && window.lj_igv_browser_ready) {{
                                    console.log('[DEBUG] Adding track to newly created browser...');
                                    const track = {{ name: '{sample_dir.name}', url: '{bam_url}', indexURL: '{bai_url}', format: 'bam', type: 'alignment', height: 600, autoScale: true, colorBy: 'tag', tag: 'SA' }};
                                    window.lj_igv.loadTrack(track).then(() => {{
                                        console.log('[DEBUG] Track loaded successfully');
                                    }}).catch(error => {{
                                        console.error('[DEBUG] Error loading track:', error);
                                    }});
                                }}
                            """
                            ui.run_javascript(js_add_track, timeout=30.0)
                        else:
                            igv_status.set_text("IGV browser creation timed out")
                            state["igv_loading"] = False
                            ui.timer(2.0, lambda: _retry_igv_creation(bam_path), once=True)
                            _clear_igv_state()
                            
                    except Exception as e:
                        igv_status.set_text(f"Failed to create IGV browser: {e}")
                        print(f"IGV browser creation error: {e}")
                        state["igv_loading"] = False
                        ui.timer(2.0, lambda: _retry_igv_creation(bam_path), once=True)
                        _clear_igv_state()
                else:
                    # Browser exists, just add/update the track
                    js_add_track = f"""
                        if (window.lj_igv && window.lj_igv_browser_ready) {{
                            console.log('Adding track to existing IGV browser...');
                    
                            const track = {{ name: '{sample_dir.name}', url: '{bam_url}', indexURL: '{bai_url}', format: 'bam', type: 'alignment', height: 600, autoScale: true }};
                            console.log('Track to load:', track);
                    
                            try {{
                            window.lj_igv.loadTrack(track);
                            console.log('Track loaded successfully');
                            }} catch (error) {{
                            console.error('Error loading track:', error);
                            }}
                        }}
                    """

                    try:
                        ui.run_javascript(js_add_track, timeout=30.0)
                        _set_igv_ready(bam_url)
                        igv_status.set_text(f"Track updated in IGV: {bam_path.name}")
                        # Clear loading flag on success
                        state["igv_loading"] = False
                    except Exception as e:
                        igv_status.set_text(f"Failed to update track: {e}")
                        # Clear loading flag on failure
                        state["igv_loading"] = False

            except Exception as e:
                igv_status.set_text(f"Error loading IGV: {e}")
                _clear_igv_state()
                # Clear loading flag on error
                if "state" in locals():
                    state["igv_loading"] = False

        # Function to refresh IGV BAM check (useful after regeneration)
        def _refresh_igv_check():
            try:
                # Only refresh if IGV isn't ready or if we need to check for new files
                if _is_igv_ready():
                    # Check if the current BAM still exists
                    current_bam = _get_igv_state().get("igv_loaded_bam")
                    if current_bam:
                        # Extract filename from URL and check if it still exists
                        bam_name = current_bam.split("/")[-1]
                        igv_dir = sample_dir / "igv"
                        clair_dir = sample_dir / "clair3"
                        candidates = [
                            igv_dir / bam_name,
                            clair_dir / bam_name,
                            sample_dir / bam_name,
                        ]

                        if any(
                            p.exists() and Path(f"{p}.bai").exists() for p in candidates
                        ):
                            igv_status.set_text("IGV is ready and BAM file is current.")
                            return

                # If we get here, we need to refresh
                _clear_igv_state()
                _check_existing_igv_bam()
            except Exception as e:
                print(f"Error refreshing IGV check: {e}")

        # Function to clear data tracks from IGV (keeps reference tracks)
        def _clear_igv_tracks():
            # Initialize state variable
            state = None
            try:
                # First check if IGV is actually ready
                if not _is_igv_ready():
                    igv_status.set_text("IGV is not ready - cannot clear tracks.")
                    return

                # Simple JavaScript to clear tracks
                js_clear = """
                    try {
                        console.log('Attempting to clear IGV tracks...');
                
                        if (window.lj_igv && typeof window.lj_igv.getTracks === 'function') {
                            const tracks = window.lj_igv.getTracks();
                            console.log('Found tracks:', tracks.length);
                    
                            if (tracks.length > 0) {
                            tracks.forEach((track, index) => {
                                try {
                                    if (typeof window.lj_igv.removeTrack === 'function') {
                                        window.lj_igv.removeTrack(track);
                                        console.log('Removed track:', index);
                                    }
                                } catch (e) {
                                    console.error('Error removing track:', e);
                                }
                            });
                            console.log('Track clearing completed');
                            }
                        } else {
                            console.log('No IGV browser or getTracks method available');
                        }
                    } catch (e) {
                        console.error('Error in track clearing:', e);
                    }
                """

                ui.run_javascript(js_clear, timeout=15.0)

                # Update status and state
                igv_status.set_text(
                    "IGV data tracks cleared. Reference tracks preserved."
                )

                # Don't clear the IGV state - just mark that data tracks are cleared
                state = _get_igv_state()
                if state:
                    state["data_tracks_cleared"] = True

            except Exception as e:
                igv_status.set_text(f"Error clearing tracks: {e}")
                print(f"Error in _clear_igv_tracks: {e}")

        # Function to wait for BAM file to stabilize (not being written to)
        def _wait_for_bam_ready(bam_path: Path, max_wait_time: int = 30) -> bool:
            """
            Wait for BAM file to be ready (not actively being written to).
            Returns True if file is ready, False if timeout.
            """
            import time
    
            start_time = time.time()
            last_size = -1
            stable_count = 0
            required_stable_checks = 2  # Need 2 consecutive checks with same size
    
            while time.time() - start_time < max_wait_time:
                try:
                    if not bam_path.exists():
                        time.sleep(0.5)
                        continue
                
                    current_size = bam_path.stat().st_size
                
                    # Check if file size has changed
                    if current_size == last_size:
                        stable_count += 1
                        if stable_count >= required_stable_checks:
                            # File size is stable for consecutive checks
                            return True
                    else:
                        stable_count = 0  # Reset counter if size changed
                
                    last_size = current_size
                    time.sleep(0.5)
                
                except (OSError, IOError):
                    # File might be temporarily unavailable
                    time.sleep(0.5)
                    continue
        
            # Timeout reached
            return False

        # Function to reload BAM data into existing IGV browser
        def _reload_bam_track():
            # Initialize state variable
            state = None
            try:
                # First check if IGV is actually ready
                if not _is_igv_ready():
                    # Try to find and load a BAM file to initialize IGV
                    igv_dir = sample_dir / "igv"
                    clair_dir = sample_dir / "clair3"
                    candidates = [
                        igv_dir / "igv_ready.bam",
                        clair_dir / "sorted_targets_exceeding_rerun.bam",
                        clair_dir / "sorted_targets_exceeding.bam",
                        sample_dir / "target.bam",
                    ]
                    
                    bam_path = None
                    for p in candidates:
                        bai = Path(f"{p}.bai")
                        if p.exists() and bai.exists():
                            bam_path = p
                            break
                    
                    if bam_path:
                        igv_status.set_text("IGV not ready, initializing...")
                        _load_igv_bam(bam_path)
                        # Wait a bit for initialization, then try again
                        ui.timer(3.0, _reload_bam_track, once=True)
                        return
                    else:
                        igv_status.set_text("IGV is not ready - cannot reload BAM.")
                        return

                # Find the current BAM file
                state = _get_igv_state()
                current_bam_url = state.get("igv_loaded_bam")

                if not current_bam_url:
                    igv_status.set_text("No BAM file currently loaded - cannot reload.")
                    return

                # Extract BAM name for display
                bam_name = current_bam_url.split("/")[-1]
        
                # Determine the actual BAM file path
                bam_path = None
                if current_bam_url.startswith("/samples/"):
                    # Extract sample directory and BAM name
                    parts = current_bam_url.split("/")
                    if len(parts) >= 4:
                        sample_dir_name = parts[2]
                        bam_file_name = parts[3]
                        bam_path = sample_dir / bam_file_name
        
                # Wait for BAM file to be ready if it exists
                if bam_path and bam_path.exists():
                    igv_status.set_text("Checking if BAM file is ready...")
                    if not _wait_for_bam_ready(bam_path):
                        igv_status.set_text("BAM file is still being updated. Please wait and try again.")
                        ui.notify("BAM file is still being updated. Please wait and try again.", type="warning")
                        return
                
                    # Also check if BAI file is ready
                    bai_path = bam_path.with_suffix(bam_path.suffix + ".bai")
                    if bai_path.exists():
                        if not _wait_for_bam_ready(bai_path):
                            igv_status.set_text("BAM index is still being updated. Please wait and try again.")
                            ui.notify("BAM index is still being updated. Please wait and try again.", type="warning")
                            return

                # JavaScript to reload the BAM track
                js_reload = f"""
                    try {{
                                                if (window.lj_igv) {{
                    // Try to get tracks list
                    let tracks = [];
                    try {{
                    // Try multiple methods to get the tracks list
                    if (typeof window.lj_igv.getTracks === 'function') {{
                        tracks = window.lj_igv.getTracks();
                    }} else if (window.lj_igv.trackList && Array.isArray(window.lj_igv.trackList)) {{
                        tracks = window.lj_igv.trackList;
                    }} else if (window.lj_igv.tracks && Array.isArray(window.lj_igv.tracks)) {{
                        tracks = window.lj_igv.tracks;
                    }} else if (window.lj_igv.roster && Array.isArray(window.lj_igv.roster.trackViews)) {{
                        tracks = window.lj_igv.roster.trackViews;
                    }} else {{
                        // Try to remove by name directly
                        if (typeof window.lj_igv.removeTrackByName === 'function') {{
                            try {{
                                window.lj_igv.removeTrackByName('{bam_name}');
                            }} catch (e) {{
                                // Ignore errors
                            }}
                        }}
                    }}
                
                    // Remove alignment tracks
                    for (let track of tracks) {{
                        try {{
                            // Only remove BAM/alignment tracks, not reference tracks
                            if (track && (track.format === 'bam' || track.format === 'cram' || track.type === 'alignment')) {{
                                // Try multiple removal methods
                                if (typeof window.lj_igv.removeTrack === 'function') {{
                                    window.lj_igv.removeTrack(track);
                                }} else if (typeof window.lj_igv.removeTrackByName === 'function') {{
                                    window.lj_igv.removeTrackByName(track.name);
                                }}
                            }}
                        }} catch (e) {{
                            // Ignore errors
                        }}
                    }}
                    }} catch (e) {{
                    // Ignore errors
                    }}
            
                    // Create the track configuration
                    const trackConfig = {{
                    name: '{bam_name}',
                    url: '{current_bam_url}',
                    indexURL: '{current_bam_url}.bai',
                    format: 'bam',
                    visibilityWindow: 1000000,
                    type: 'alignment',
                    height: 600,
                    autoScale: true,
                    colorBy: 'tag',
                    tag: 'SA'
                    }};
            
                    // Add the track to the browser
                    if (typeof window.lj_igv.loadTrack === 'function') {{
                    window.lj_igv.loadTrack(trackConfig);
                    }} else {{
                    // Try alternative method
                    if (typeof window.lj_igv.addTrack === 'function') {{
                        window.lj_igv.addTrack(trackConfig);
                    }}
                    }}
            
                    // Force a redraw
                    if (typeof window.lj_igv.redraw === 'function') {{
                    window.lj_igv.redraw();
                    }} else if (typeof window.lj_igv.update === 'function') {{
                    window.lj_igv.update();
                    }}
                }}
            }} catch (e) {{
                // Ignore errors
            }}
                """

                ui.run_javascript(js_reload, timeout=20.0)

                # Update status
                igv_status.set_text(f"Loading BAM track: {bam_name}")
        
                # Change button text to "Reload BAM" after first load
                reload_bam_button.set_text("Reload BAM")

                # Mark that data tracks are no longer cleared
                if state:
                    state["data_tracks_cleared"] = False
                
                # Auto-load the BED file after BAM is loaded
                _load_target_bed()

            except Exception as e:
                igv_status.set_text(f"Error reloading BAM: {e}")
                print(f"Error in _reload_bam_track: {e}")

        # Function to load target BED file into IGV
        def _load_target_bed():
            """Load the target BED file as a track in IGV"""
            try:
                # First check if IGV is actually ready
                if not _is_igv_ready():
                    igv_status.set_text("IGV is not ready - cannot load BED file.")
                    return

                # Get the target panel information
                state = _get_igv_state()
        
                # Check if BED file has already been loaded
                if state.get("bed_file_loaded", False):
                    return  # BED file already loaded, don't load again
        
                # Try to read the panel from master.csv
                target_panel = None
                bed_file_path = None
                try:
                    master_csv_path = sample_dir / "master.csv"
                    if master_csv_path.exists():
                        import pandas as pd
                        df = pd.read_csv(master_csv_path)
                        if not df.empty and "analysis_panel" in df.columns:
                            target_panel = str(df.iloc[0]["analysis_panel"]).strip()
                    
                            # Map panel to BED filename
                            bed_file_mapping = {
                            "rCNS2": "rCNS2_panel_name_uniq.bed",
                            "AML": "AML_panel_name_uniq.bed",
                            "Sarcoma": "Sarcoma_panel_name_uniq.bed"
                            }
                    
                            bed_filename = bed_file_mapping.get(target_panel, f"{target_panel}_panel_name_uniq.bed")
                    
                            # Try to find the BED file in robin resources
                            try:
                                from robin import resources
                                bed_file_path = os.path.join(
                                    os.path.dirname(os.path.abspath(resources.__file__)),
                                    bed_filename
                                )
                                if not os.path.exists(bed_file_path):
                                    bed_file_path = None
                            except Exception:
                                pass
                    
                            # Fallback paths
                            if not bed_file_path:
                                possible_paths = [
                                bed_filename,
                                f"data/{bed_filename}",
                                f"/usr/local/share/{bed_filename}",
                            ]
                            for path in possible_paths:
                                if os.path.exists(path):
                                    bed_file_path = path
                                    break
                except Exception as e:
                    print(f"Error reading panel information: {e}")

                if not bed_file_path or not os.path.exists(bed_file_path):
                    igv_status.set_text(f"Could not find BED file for panel: {target_panel}")
                    ui.notify(f"Target BED file not found for panel: {target_panel}", type="warning")
                    return

                # Mount the BED file's directory
                bed_dir = os.path.dirname(bed_file_path)
                bed_name = os.path.basename(bed_file_path)
                bed_url = f"/robin_resources/{bed_name}"
        
                try:
                    if app is not None:
                        app.add_static_files("/robin_resources", bed_dir)
                except Exception:
                    # Already mounted or in tests without app
                    pass

                # JavaScript to load the BED file as a track
                js_load_bed = f"""
                    try {{
                        if (window.lj_igv) {{
                            // Create the track configuration for the BED file
                            const trackConfig = {{
                            name: '{bed_name}',
                            url: '{bed_url}',
                            format: 'bed',
                            indexed: false,
                            displayMode: 'COLLAPSED',
                            height: 40,
                            color: '#0072B2'
                            }};
                    
                            // Load the track into IGV
                            if (typeof window.lj_igv.loadTrack === 'function') {{
                            window.lj_igv.loadTrack(trackConfig);
                            }} else {{
                            // Try alternative method
                            if (typeof window.lj_igv.addTrack === 'function') {{
                                window.lj_igv.addTrack(trackConfig);
                            }}
                            }}
                        }}
                    }} catch (e) {{
                        // Ignore errors
                    }}
                """

                ui.run_javascript(js_load_bed, timeout=20.0)
        
                # Mark that BED file has been loaded
                state["bed_file_loaded"] = True
        
                # Update status
                igv_status.set_text(f"Loaded BED file: {bed_name}")
                ui.notify(f"Target BED file loaded: {bed_name}", type="positive")

            except Exception as e:
                igv_status.set_text(f"Error loading BED: {e}")
                print(f"Error in _load_target_bed: {e}")

        # Function to check console for errors
        def _check_console_errors():
            """Check the browser console for any JavaScript errors."""
            js_check_console = """
            try {
                const errors = [];

                // Check for basic IGV functionality
                if (typeof igv === 'undefined') {
                    errors.push('IGV library not loaded');
                } else if (typeof igv.createBrowser !== 'function') {
                    errors.push('IGV createBrowser method not available');
                } else {
                    errors.push(`IGV version: ${igv.version || 'unknown'}`);
                }
        
                return errors;
            } catch (e) {
                return ['Error checking console: ' + e.message];
            }
            """

            try:
                errors = ui.run_javascript(js_check_console, timeout=10.0)
                if errors and len(errors) > 0:
                    error_text = "Console errors found:\n" + "\n".join(errors)
                    ui.notify(error_text, type="warning")
                    print(error_text)
                else:
                    ui.notify("No console errors found", type="positive")
            except Exception as e:
                ui.notify(f"Error checking console: {e}", type="negative")

        # Function to create an empty IGV browser
        def _create_empty_igv():
            """Create an empty IGV browser with no tracks."""
            try:
                if not _check_igv_library():
                    ui.notify("IGV library not available", type="warning")
                    return

                # First, let's try to find the element and wait for it to be ready
                js_find_element = """
                    console.log('Looking for IGV div element...');
                    console.log('Expected ID: igv-container');
            
                    // List all divs on the page
                    const allDivs = document.querySelectorAll('div');
                    console.log('Total divs found:', allDivs.length);
            
                    allDivs.forEach((div, i) => {
                        console.log(`Div ${i}:`, {
                            id: div.id,
                            className: div.className,
                            tagName: div.tagName,
                            offsetWidth: div.offsetWidth,
                            offsetHeight: div.offsetHeight
                        });
                    });
            
                    // Try to find our specific element
                    let el = document.getElementById('igv-container');
                    if (!el) {
                        console.log('Element not found by ID, trying by class...');
                        const candidates = document.querySelectorAll('.w-full.h-\\[600px\\].border');
                        console.log('Found candidates by class:', candidates.length);
                        if (candidates.length > 0) {
                            el = candidates[0];
                            console.log('Using first candidate by class');
                        }
                    }
            
                    if (el) {
                        console.log('Element found:', el);
                        console.log('Element dimensions:', el.offsetWidth, 'x', el.offsetHeight);
                        return true;
                    } else {
                        console.log('No suitable element found');
                        return false;
                    }
                """

                # Check if element exists first
                element_found = ui.run_javascript(js_find_element, timeout=10.0)
                if not element_found:
                    ui.notify(
                        "IGV div element not found. Please wait for the page to load completely.",
                        type="warning",
                    )
                    return

                # Now create the IGV browser
                js_empty = """
                    try {
                        console.log('Creating empty IGV browser...');
                
                        // Find the element again
                        let el = document.getElementById('igv-container');
                        if (!el) {
                            el = document.querySelector('.w-full.h-\\[600px\\].border');
                        }
                
                        if (!el) {
                            throw new Error('Element still not found');
                        }
                
                        const options = { genome: 'hg38' };
                        console.log('Creating empty IGV browser with options:', options);
                
                        igv.createBrowser(el, options).then(b => {
                            window.lj_igv = b; 
                            window.lj_igv_browser_ready = true;
                            console.log('Empty IGV browser created successfully');
                        }).catch(error => {
                            console.error('Empty IGV browser creation failed:', error);
                            window.lj_igv_browser_ready = false;
                        });
                    } catch (error) {
                        console.error('Error creating empty IGV:', error);
                    }
                """

                ui.run_javascript(js_empty, timeout=30.0)
                ui.notify("Empty IGV browser creation initiated", type="info")

            except Exception as e:
                ui.notify(f"Error creating empty IGV: {e}", type="negative")
                print(f"Error in empty IGV: {e}")

        # Function to debug IGV state
        def _debug_igv_state():
            try:
                js_debug = """
                console.log('=== IGV Debug Info ===');
                console.log('window.lj_igv:', window.lj_igv);
                console.log('window.lj_igv_browser_ready:', window.lj_igv_browser_ready);
        
                if (window.lj_igv) {
                    console.log('IGV browser exists');
                    console.log('Available methods:', Object.getOwnPropertyNames(window.lj_igv));
                } else {
                    console.log('No IGV browser found');
                }
        
                console.log('Python state check requested');
                """

                ui.run_javascript(js_debug, timeout=10.0)

                # Also show Python state
                state = _get_igv_state()
                debug_info = f"Python IGV State: initialized={state.get('igv_initialized')}, browser_ready={state.get('igv_browser_ready')}, bam={state.get('igv_loaded_bam')}, tracks_cleared={state.get('tracks_cleared')}"
                igv_status.set_text(debug_info)

            except Exception as e:
                igv_status.set_text(f"Debug error: {e}")
                print(f"Error in debug: {e}")

        # IGV control buttons
        def _simple_reload_bam():
            """Simple reload function."""
            _check_existing_igv_bam()
        
        reload_bam_button = ui.button("Load BAM", on_click=_simple_reload_bam)

        # Add IGV library status check timer once after a delay
        ui.timer(3.0, _check_igv_library, once=True)
        
        # Check for existing BAM files after browser is initialized (give it time to create)
        ui.timer(3.0, _check_existing_igv_bam, once=True)

        # BAM generation buttons
        async def _trigger_build_sorted_bam(force_regenerate: bool = False) -> None:
            try:
                # Debug: Check what's in the launcher
                debug_info = f"launcher type: {type(launcher)}, workflow_runner: {getattr(launcher, 'workflow_runner', 'None')}"

                if force_regenerate:
                    ui.notify(
                        "Force regenerate mode: will recreate IGV BAM even if it exists",
                        type="info",
                    )

                # Check if target.bam exists (required for IGV BAM generation)
                target_bam = sample_dir / "target.bam"
                if not target_bam.exists():
                    ui.notify(
                        "No target.bam file found. Please run target analysis first.",
                        type="warning",
                    )
                    return

                # Check if IGV BAM already exists (only for normal generation, not force regenerate)
                igv_bam = sample_dir / "igv" / "igv_ready.bam"
                if (
                    not force_regenerate
                    and igv_bam.exists()
                    and (sample_dir / "igv" / "igv_ready.bam.bai").exists()
                ):
                    ui.notify("IGV BAM already exists and is ready.", type="positive")
                    return

                if force_regenerate:
                    ui.notify(
                        "Proceeding with force regenerate - existing files will be overwritten",
                        type="info",
                    )

                # Check if we have a workflow runner
                if (
                    not hasattr(launcher, "workflow_runner")
                    or not launcher.workflow_runner
                ):
                    ui.notify(
                        f"No workflow runner available. Debug: {debug_info}. GUI is running but workflow integration is not available. This may be a configuration issue.",
                        type="warning",
                    )
                    return

                # Submit the IGV BAM generation job to the workflow system
                try:
                    runner = launcher.workflow_runner
                    sample_id = sample_dir.name

                    # Check if it's a Ray workflow or Simple workflow
                    if hasattr(runner, "submit_sample_job"):
                        # Simple workflow
                        success = runner.submit_sample_job(
                            str(sample_dir), "igv_bam", sample_id, force_regenerate
                        )
                    elif hasattr(runner, "manager") and hasattr(
                        runner.manager, "submit_sample_job"
                    ):
                        # Ray workflow - need to get the coordinator
                        # This is more complex for Ray, so we'll provide a fallback
                        ui.notify(
                            "Ray workflow detected. IGV BAM generation should happen automatically "
                            "after target analysis completes. If you need to regenerate it, "
                            "please restart the workflow or check the logs.",
                            type="info",
                        )
                        return
                    else:
                        ui.notify(
                            "Unknown workflow type. Cannot submit IGV BAM job.",
                            type="warning",
                        )
                        return

                    if success:
                        action = "regenerated" if force_regenerate else "generated"
                        ui.notify(
                            f"IGV BAM {action} job submitted to workflow queue!",
                            type="positive",
                        )
                        # Refresh the IGV check after successful job submission
                        _refresh_igv_check()
                    else:
                        ui.notify(
                            "Failed to submit IGV BAM generation job. Check logs for details.",
                            type="warning",
                        )

                except Exception as e:
                    ui.notify(f"Error submitting IGV BAM job: {e}", type="negative")

            except Exception as e:
                try:
                    ui.notify(f"Error checking IGV BAM status: {e}", type="negative")
                except Exception:
                    # Client may have disconnected, ignore UI errors
                    pass


def add_coverage_section(launcher: Any, sample_dir: Path) -> None:
    """Build the Coverage UI section and attach refresh timers.

    This mirrors the existing inline implementation but lives in a reusable module.
    """
    # Check for development environment variable to show/hide testing features
    is_development_mode = os.environ.get("ROBIN_DEV_MODE", "").lower() in ("1", "true", "yes", "on")
    
    with ui.card().classes("w-full").props("id=analysis-detail-coverage"):
        ui.label("Coverage").classes("text-lg font-semibold mb-2")
        # Summary row (quality + metrics)
        with ui.row().classes("w-full items-center justify-between mb-2"):
            with ui.column().classes("gap-1"):
                ui.label("Coverage Analysis").classes("text-sm font-medium")
                with ui.row().classes("items-center gap-2"):
                    cov_quality_label = ui.label("Quality: --").classes(
                        "text-gray-600 font-medium"
                    )
                    cov_quality_badge = ui.label("--x").classes(
                        "px-2 py-1 rounded bg-gray-100 text-gray-600"
                    )
            with ui.column().classes("gap-1 items-end"):
                cov_global_lbl = ui.label("Global Estimated Coverage: --x").classes(
                    "text-sm text-gray-600"
                )
                cov_target_lbl = ui.label("Targets Estimated Coverage: --x").classes(
                    "text-sm text-gray-600"
                )
                cov_enrich_lbl = ui.label("Estimated enrichment: --x").classes(
                    "text-sm text-gray-600"
                )
        with ui.card().classes("w-full"):
            # Per Chromosome Target Coverage (grouped bar — design.md §9.5)
            _cp0 = _cov_chrome_palette()
            echart_target_cov = ui.echart(
                    {
                        "backgroundColor": "transparent",
                        "textStyle": {"color": _cp0["axis"]},
                        "title": {
                            "text": "Per Chromosome Target Coverage",
                            "left": "center",
                            "top": 10,
                            "textStyle": {"fontSize": 15, "color": _cp0["title"]},
                        },
                        "legend": {
                            "data": ["On Target", "Off Target"],
                            "left": 10,
                            "top": "center",
                            "orient": "vertical",
                            "itemGap": 10,
                            "textStyle": {"color": _cp0["legend"]},
                        },
                        "tooltip": {
                            "trigger": "axis",
                            "axisPointer": {"type": "shadow"},
                            **_cov_tooltip_option(),
                        },
                        "grid": {
                            "left": "12%",
                            "right": "8%",
                            "bottom": "12%",
                            "top": "20%",
                            "containLabel": True,
                        },
                        "xAxis": {
                            "type": "category",
                            "data": [],
                            "axisLabel": {
                                "rotate": 0,
                                "interval": 0,
                                "fontSize": 11,
                                "color": _cp0["axis"],
                            },
                            "axisLine": {"lineStyle": {"color": _cp0["axis"]}},
                        },
                        "yAxis": {
                            "type": "value",
                            "name": "Coverage (×)",
                            "nameTextStyle": {"color": _cp0["axis"]},
                            "axisLabel": {"color": _cp0["axis"]},
                            "splitLine": {
                                "lineStyle": {"color": _cp0["split"], "type": "dashed"},
                            },
                        },
                        "series": [],
                    }
                ).classes("w-full h-64")

        with ui.card().classes("w-full"):
            ui.label("Coverage Over Time").classes("text-lg font-semibold mb-2")
            _cp1 = _cov_chrome_palette()
            _ln = _cov_line_primary()
            _fill = _cov_area_fill_primary()
            echart_time = ui.echart(
                {
                    "backgroundColor": "transparent",
                    "textStyle": {"color": _cp1["axis"]},
                    "title": {
                        "text": "Cumulative coverage depth",
                        "subtext": "Estimated coverage (×) over time",
                        "left": "center",
                        "top": 8,
                        "textStyle": {"fontSize": 15, "color": _cp1["title"]},
                        "subtextStyle": {"fontSize": 11, "color": _cp1["axis"]},
                    },
                    "tooltip": {
                        "trigger": "axis",
                        "axisPointer": {
                            "type": "cross",
                            "lineStyle": {"color": _cp1["axis"], "width": 1},
                            "crossStyle": {"color": _cp1["axis"]},
                        },
                        **_cov_tooltip_option(),
                    },
                    "grid": {
                        "left": "8%",
                        "right": "6%",
                        "bottom": "12%",
                        "top": "22%",
                        "containLabel": True,
                    },
                    "xAxis": {
                        "type": "time",
                        "axisLabel": {"color": _cp1["axis"]},
                        "axisLine": {"lineStyle": {"color": _cp1["axis"]}},
                        "splitLine": {"lineStyle": {"color": _cp1["split"]}},
                    },
                    "yAxis": {
                        "type": "value",
                        "name": "Estimated coverage (×)",
                        "nameTextStyle": {"color": _cp1["axis"]},
                        "axisLabel": {"color": _cp1["axis"]},
                        "splitLine": {
                            "lineStyle": {"color": _cp1["split"], "type": "dashed"},
                        },
                    },
                    "series": [
                        {
                            "name": "Coverage",
                            "type": "line",
                            "smooth": True,
                            "showSymbol": False,
                            "data": [],
                            "lineStyle": {"width": 2.5, "color": _ln},
                            "itemStyle": {"color": _ln},
                            "areaStyle": {"color": _fill},
                        }
                    ],
                }
            ).classes("w-full h-72")

        # Target Coverage Over Time Analysis
        with ui.card().classes("w-full mt-4"):
            ui.label("Target Coverage Over Time Analysis").classes("text-lg font-semibold mb-2")
            ui.label(
                "Analyze mean target coverage over time and identify significant outliers. "
                "Outliers are detected using standard deviation method: values more than 2 SD above or below the mean for each gene."
            ).classes("text-sm text-gray-600 mb-4")
            
            def _get_coverage_state():
                key = str(sample_dir)
                if hasattr(launcher, "_coverage_state"):
                    if key not in launcher._coverage_state:
                        launcher._coverage_state[key] = {}
                    return launcher._coverage_state[key]
                launcher._coverage_state = {}
                launcher._coverage_state[key] = {}
                return launcher._coverage_state[key]
            
            coverage_state = _get_coverage_state()
            stored_outlier_limit = coverage_state.get("target_cov_outlier_limit", 10)
            try:
                stored_outlier_limit = int(stored_outlier_limit)
            except (TypeError, ValueError):
                stored_outlier_limit = 10

            # Chart container
            target_coverage_time_chart_state = {"container": None, "chart": None}
            outlier_limit_state = {"value": max(1, stored_outlier_limit)}
            
            def _plot_target_coverage_over_time():
                """Load target_coverage_time.csv and plot mean coverage with outlier detection"""
                try:
                    time_coverage_file = sample_dir / "target_coverage_time.csv"
                    if not time_coverage_file.exists():
                        if target_coverage_time_chart_state["container"]:
                            with target_coverage_time_chart_state["container"]:
                                ui.label("No target_coverage_time.csv file found.").classes("text-gray-600")
                        return
                    
                    # Load data
                    df = pd.read_csv(time_coverage_file)
                    if df.empty:
                        if target_coverage_time_chart_state["container"]:
                            with target_coverage_time_chart_state["container"]:
                                ui.label("No data available in target_coverage_time.csv").classes("text-gray-600")
                        return
                    if not {'chrom', 'startpos', 'endpos', 'name'}.issubset(df.columns):
                        if target_coverage_time_chart_state["container"]:
                            with target_coverage_time_chart_state["container"]:
                                ui.label("target_coverage_time.csv missing chrom/startpos/endpos/name columns").classes("text-gray-600")
                        return
                    
                    # Convert timestamp to datetime (milliseconds to datetime)
                    df['datetime'] = pd.to_datetime(df['timestamp'], unit='ms')
                    
                    # Create unique target identifier (chrom, startpos, endpos) — vectorized
                    df["target_key"] = list(
                        zip(
                            df["chrom"].astype(str),
                            df["startpos"].astype(np.int64),
                            df["endpos"].astype(np.int64),
                        )
                    )
                    # Display labels (same format for all rows; avoids slow df.apply)
                    _start_mb = (df["startpos"].astype(np.float64) / 1e6).round(2)
                    df["target_label"] = (
                        df["name"].astype(str)
                        + " ("
                        + df["chrom"].astype(str)
                        + " "
                        + _start_mb.astype(str)
                        + " Mb)"
                    )
                    
                    # Calculate mean coverage per timepoint
                    mean_coverage = df.groupby('timestamp')['coverage'].mean().reset_index()
                    mean_coverage['datetime'] = pd.to_datetime(mean_coverage['timestamp'], unit='ms')
                    
                    # Detect outliers using standard deviation method (Z-score)
                    def detect_outliers_sd(series, num_sd=2.0):
                        """Detect outliers using standard deviation method (Z-score)
                        
                        Args:
                            series: Pandas Series of values
                            num_sd: Number of standard deviations from mean (default: 2.0)
                        
                        Returns:
                            Boolean Series indicating which values are outliers
                        """
                        mean = series.mean()
                        std = series.std()
                        
                        # Handle case where std is 0 (all values are the same)
                        if std == 0:
                            return pd.Series([False] * len(series), index=series.index)
                        
                        lower_bound = mean - num_sd * std
                        upper_bound = mean + num_sd * std
                        return (series < lower_bound) | (series > upper_bound)
                    
                    # Detect outliers: compare each gene's coverage to the distribution of ALL genes
                    # This detects genes that are outliers compared to other genes, not just timepoints
                    outliers = []
                    
                    # Calculate global statistics across all genes at each timepoint
                    # This allows us to detect genes that are outliers relative to the population
                    for timestamp in df['timestamp'].unique():
                        timepoint_data = df[df['timestamp'] == timestamp].copy()
                        if len(timepoint_data) < 3:  # Need at least 3 genes to calculate SD
                            continue
                        
                        # Calculate mean and SD across all genes at this timepoint
                        global_mean = timepoint_data['coverage'].mean()
                        global_std = timepoint_data['coverage'].std()
                        
                        if global_std == 0:  # All genes have same coverage
                            continue
                        
                        # Detect outliers: genes > 2 SD from global mean at this timepoint
                        lower_bound = global_mean - 2.0 * global_std
                        upper_bound = global_mean + 2.0 * global_std
                        
                        outlier_mask = (timepoint_data['coverage'] < lower_bound) | (timepoint_data['coverage'] > upper_bound)
                        outlier_points = timepoint_data[outlier_mask]
                        
                        for _, row in outlier_points.iterrows():
                            tk = row['target_key']
                            tk = (str(tk[0]), int(tk[1]), int(tk[2])) if isinstance(tk, (list, tuple)) else (str(row['chrom']), int(row['startpos']), int(row['endpos']))
                            outliers.append({
                                'target_key': tk,
                                'target_label': row['target_label'],
                                'gene': row['name'],
                                'timestamp': row['timestamp'],
                                'datetime': row['datetime'],
                                'coverage': row['coverage'],
                                'reads': row['reads'],
                                'reads_per_length': row['reads_per_length'],
                                'type': 'high' if row['coverage'] > global_mean else 'low',
                                'global_mean': global_mean,
                                'global_std': global_std
                            })
                    
                    outliers_df = pd.DataFrame(outliers) if outliers else pd.DataFrame()
                    
                    # Identify unique outlier targets (regions that are outliers at any timepoint)
                    outlier_targets = []
                    if not outliers_df.empty:
                        # Preserve order: (target_key, target_label) for consistent display
                        seen = set()
                        for _, row in outliers_df.iterrows():
                            key = row['target_key']
                            if isinstance(key, (list, tuple)):
                                key = tuple(key)
                            if key not in seen:
                                seen.add(key)
                                outlier_targets.append((key, row['target_label']))
                    
                    try:
                        outlier_limit = int(outlier_limit_state["value"])
                    except (TypeError, ValueError):
                        outlier_limit = 10
                    outlier_limit = max(1, outlier_limit)
                    
                    # Prepare data for ECharts (design.md §9.5 — mean, ±2 SD band, dim in-range, outliers)
                    outlier_key_set = set()
                    for k, _ in outlier_targets:
                        kk = tuple(k) if isinstance(k, (list, tuple)) else k
                        outlier_key_set.add(kk)

                    mean_coverage_data = [
                        [int(ts), float(cov)]
                        for ts, cov in zip(
                            mean_coverage["timestamp"], mean_coverage["coverage"]
                        )
                    ]

                    mc_vals = mean_coverage["coverage"].astype(float)
                    m_mu = float(mc_vals.mean())
                    m_sd = float(mc_vals.std()) if len(mc_vals) > 1 else 0.0
                    band_low = max(0.0, m_mu - 2.0 * m_sd)
                    band_high = max(band_low + 1e-9, m_mu + 2.0 * m_sd)
                    band_fill = _cov_band_fill_rgba()

                    band_series = {
                        "name": "_band",
                        "type": "line",
                        "data": [],
                        "silent": True,
                        "showInLegend": False,
                        "z": 0,
                        "markArea": {
                            "silent": True,
                            "itemStyle": {"color": band_fill},
                            "data": [[{"yAxis": band_low}, {"yAxis": band_high}]],
                        },
                    }

                    mean_c = _cov_mean_line_color()
                    mean_series = {
                        "name": "Mean (population)",
                        "type": "line",
                        "smooth": True,
                        "data": mean_coverage_data,
                        "symbol": "none",
                        "z": 80,
                        "lineStyle": {"width": 4, "type": "dashed", "color": mean_c},
                        "itemStyle": {"color": mean_c},
                    }

                    outlier_gene_series: List[Dict[str, Any]] = []
                    if outlier_targets:
                        target_max_cov = df.groupby("target_key")["coverage"].max()

                        def _sort_key(item: Any) -> float:
                            key, _ = item
                            k = tuple(key) if isinstance(key, (list, tuple)) else key
                            return float(target_max_cov.get(k, 0))

                        outlier_targets_sorted = sorted(
                            outlier_targets,
                            key=_sort_key,
                            reverse=True,
                        )
                        targets_to_plot = outlier_targets_sorted[:outlier_limit]

                        for idx, (target_key, target_label) in enumerate(targets_to_plot):
                            key = (
                                tuple(target_key)
                                if isinstance(target_key, (list, tuple))
                                else target_key
                            )
                            mask = (
                                (df["chrom"].astype(str) == str(key[0]))
                                & (df["startpos"] == key[1])
                                & (df["endpos"] == key[2])
                            )
                            target_data = df[mask].copy()
                            if len(target_data) < 2:
                                continue
                            target_data = target_data.sort_values("timestamp")
                            series_data = [
                                [int(ts), float(cov)]
                                for ts, cov in zip(
                                    target_data["timestamp"], target_data["coverage"]
                                )
                            ]
                            out_sub = outliers_df[
                                (
                                    outliers_df["target_key"].apply(
                                        lambda x: (
                                            tuple(x)
                                            if isinstance(x, (list, tuple))
                                            else x
                                        )
                                    )
                                    == key
                                )
                            ]
                            is_high = (
                                out_sub["type"].value_counts().get("high", 0)
                                > out_sub["type"].value_counts().get("low", 0)
                            )
                            gc = _cov_gene_line_color(idx)
                            outlier_gene_series.append(
                                {
                                    "name": _cov_legend_label(target_label),
                                    "type": "line",
                                    "smooth": True,
                                    "data": series_data,
                                    "symbol": "none",
                                    "z": 12,
                                    "itemStyle": {"color": gc},
                                    "lineStyle": {
                                        "width": 2.5,
                                        "type": "dashed" if not is_high else "solid",
                                        "color": gc,
                                    },
                                    "label": {"show": False},
                                    "labelLayout": {"hideOverlap": True},
                                    "emphasis": {
                                        "focus": "series",
                                        "lineStyle": {"width": 4},
                                    },
                                    "blur": {"lineStyle": {"opacity": 0.15}},
                                }
                            )

                    dim_gene_series: List[Dict[str, Any]] = []
                    for tkey, group in df.groupby("target_key", sort=False):
                        if len(dim_gene_series) >= _COV_MAX_DIM_GENES:
                            break
                        kk = tuple(tkey) if isinstance(tkey, (list, tuple)) else tkey
                        if kk in outlier_key_set:
                            continue
                        row0 = group.sort_values("timestamp").iloc[0]
                        target_label = row0["target_label"]
                        key = kk
                        target_data = group.copy()
                        if len(target_data) < 2:
                            continue
                        target_data = target_data.sort_values("timestamp")
                        series_data = [
                            [int(ts), float(cov)]
                            for ts, cov in zip(
                                target_data["timestamp"], target_data["coverage"]
                            )
                        ]
                        dim_gene_series.append(
                            {
                                "name": _cov_legend_label(
                                    f"{target_label} · in-range"
                                ),
                                "type": "line",
                                "smooth": True,
                                "data": series_data,
                                "symbol": "none",
                                "z": 4,
                                "itemStyle": {
                                    "color": _cov_dim_muted_line(),
                                    "opacity": 0.42,
                                },
                                "lineStyle": {
                                    "width": 1.5,
                                    "type": "dashed",
                                    "color": _cov_dim_muted_line(),
                                    "opacity": 0.42,
                                },
                                "label": {"show": False},
                                "emphasis": {
                                    "focus": "series",
                                    "lineStyle": {"width": 3, "opacity": 1},
                                },
                                "blur": {"lineStyle": {"opacity": 0.08}},
                            }
                        )

                    # Create/update chart
                    if target_coverage_time_chart_state["container"] is None:
                        target_coverage_time_chart_state["container"] = ui.column().classes("w-full")
                    
                    with target_coverage_time_chart_state["container"]:
                        # Clear existing content
                        target_coverage_time_chart_state["container"].clear()
                        
                        # Summary statistics
                        with ui.row().classes("w-full mb-4 gap-3"):
                            ui.label(f"Total timepoints: {len(mean_coverage)}").classes("text-sm")
                            ui.label(f"Total targets: {len(df['target_key'].unique())}").classes("text-sm")
                            ui.label(f"Outliers detected: {len(outliers_df)}").classes("text-sm")
                        
                        # Chart
                        _cp_t = _cov_chrome_palette()
                        chart_config = {
                            "backgroundColor": "transparent",
                            "textStyle": {"color": _cp_t["axis"]},
                            "title": {
                                "text": "Individual target coverage (outlier analysis)",
                                "subtext": "Mean ±2σ band; in-range targets dimmed; outliers in colour",
                                "left": "center",
                                "top": 6,
                                "textStyle": {
                                    "fontSize": 15,
                                    "color": _cp_t["title"],
                                },
                                "subtextStyle": {
                                    "fontSize": 11,
                                    "color": _cp_t["axis"],
                                },
                            },
                            # No JS formatter here — options are JSON-serialied; functions
                            # become strings and break the tooltip (see classification charts).
                            "tooltip": {
                                "trigger": "axis",
                                "axisPointer": {
                                    "type": "cross",
                                    "lineStyle": {"color": _cp_t["axis"], "width": 1},
                                    "crossStyle": {"color": _cp_t["axis"]},
                                },
                                "confine": True,
                                "appendToBody": True,
                                **_cov_tooltip_option(),
                            },
                            "legendHoverLink": True,
                            # Legend below plot (horizontal scroll) — avoids overlapping Y-axis / grid.
                            "legend": {
                                "type": "scroll",
                                "orient": "horizontal",
                                "left": "center",
                                "bottom": 2,
                                "itemGap": 10,
                                "itemWidth": 10,
                                "itemHeight": 10,
                                "width": "96%",
                                "pageIconColor": _cp_t["title"],
                                "pageTextStyle": {"color": _cp_t["axis"]},
                                "textStyle": {
                                    "fontSize": 11,
                                    "color": _cp_t["legend"],
                                },
                            },
                            "grid": {
                                "left": "10%",
                                "right": "6%",
                                "top": "18%",
                                "bottom": "32%",
                                "containLabel": True,
                            },
                            "xAxis": {
                                "type": "time",
                                "name": "Time",
                                "nameGap": 30,
                                "axisLabel": {"color": _cp_t["axis"]},
                                "axisLine": {"lineStyle": {"color": _cp_t["axis"]}},
                                "splitLine": {"lineStyle": {"color": _cp_t["split"]}},
                            },
                            "yAxis": {
                                "type": "value",
                                "name": "Coverage (×)",
                                "nameGap": 50,
                                "axisLabel": {"color": _cp_t["axis"]},
                                "axisLine": {"lineStyle": {"color": _cp_t["axis"]}},
                                "splitLine": {
                                    "lineStyle": {"color": _cp_t["split"], "type": "dashed"},
                                },
                            },
                            "series": [band_series]
                            + dim_gene_series
                            + outlier_gene_series
                            + [mean_series],
                        }
                        
                        target_coverage_time_chart_state["chart"] = ui.echart(
                            chart_config
                        ).classes("w-full min-h-[420px] h-[480px]")
                        _apply_target_coverage_time_analysis_chrome(
                            target_coverage_time_chart_state["chart"]
                        )
                        
                        # Show summary of outlier targets
                        if outlier_targets:
                            outlier_count = len(outlier_targets)
                            ui.label(f"Showing profiles for {min(outlier_count, outlier_limit)} outlier targets (out of {outlier_count} total)").classes("text-sm text-gray-600 mt-2")
                        else:
                            ui.label("No significant outliers detected.").classes("text-gray-600 mt-2")
                
                except Exception as e:
                    logging.error(f"Error plotting target coverage over time: {e}")
                    import traceback
                    logging.error(traceback.format_exc())
                    if target_coverage_time_chart_state["container"]:
                        with target_coverage_time_chart_state["container"]:
                            ui.label(f"Error: {str(e)}").classes("text-red-600")
            
            # Defer plot work so the sample page shell and websocket can finish first
            ui.timer(0.05, _plot_target_coverage_over_time, once=True)
            
            def _set_outlier_limit(e) -> None:
                try:
                    outlier_limit_state["value"] = int(e.value)
                except (TypeError, ValueError):
                    outlier_limit_state["value"] = 10
                outlier_limit_state["value"] = max(1, outlier_limit_state["value"])
                coverage_state["target_cov_outlier_limit"] = outlier_limit_state["value"]
                _plot_target_coverage_over_time()
            
            # Refresh button + outlier limit control
            with ui.row().classes("w-full mt-4 items-center gap-3"):
                ui.label("Outliers to show").classes("text-sm text-gray-600")
                ui.number(
                    value=outlier_limit_state["value"],
                    min=1,
                    max=50,
                    step=1,
                    format="%.0f",
                    on_change=_set_outlier_limit,
                ).props("dense").classes("w-24")
                ui.button(
                    "Refresh Analysis",
                    on_click=_plot_target_coverage_over_time
                ).classes("ml-2")

        with ui.column().classes("w-full target-coverage-panel").props(
            "id=analysis-detail-target-coverage"
        ):
            # Add target panel legend
            def _get_target_panel_info():
                """Get the target panel information from master.csv with fallback mechanisms"""
                try:
                    # Read from master.csv
                    master_csv_path = sample_dir / "master.csv"
                    if master_csv_path.exists():
                        import pandas as pd
                        df = pd.read_csv(master_csv_path)
                        if not df.empty and "analysis_panel" in df.columns:
                            panel = df.iloc[0]["analysis_panel"]
                            if panel and str(panel).strip() != "" and str(panel).strip().lower() != "nan":
                                return str(panel).strip()
                    
                    # Fallback 1: Try to detect panel from BED files in the sample directory
                    bed_files = list(sample_dir.glob("*.bed"))
                    if bed_files:
                        # Look for panel-specific BED files
                        for bed_file in bed_files:
                            bed_name = bed_file.stem.lower()
                            if "rcns2" in bed_name or "rCNS2" in bed_name:
                                return "rCNS2"
                            elif "aml" in bed_name:
                                return "AML"
                            elif "sarcoma" in bed_name:
                                return "Sarcoma"
                    
                    # Fallback 2: Try to detect from target analysis output files
                    target_files = list(sample_dir.glob("*target*.csv")) + list(sample_dir.glob("*coverage*.csv"))
                    if target_files:
                        # This is a heuristic - if we have target analysis files, 
                        # we can assume it's likely a known panel
                        return "Unknown Panel"
                    
                    return ""  # No panel found
                except Exception as e:
                    _log_notify(f"Exception in _get_target_panel_info: {e}", level="error", notify=False)
                    return ""  # No default fallback
            
            target_panel = _get_target_panel_info()
            
            # Panel legend with color coding
            panel_colors = {
                "rCNS2": ("bg-blue-100", "text-blue-800", "rCNS2 Panel"),
                "AML": ("bg-green-100", "text-green-800", "AML Panel"), 
                "Sarcoma": ("bg-orange-100", "text-orange-800", "Sarcoma Panel"),
                "Unknown Panel": ("bg-yellow-100", "text-yellow-800", "Unknown Panel")
            }
            
            if not target_panel:
                # No panel found - show warning
                panel_color_classes, panel_text_classes, panel_display_name = ("bg-red-100", "text-red-800", "Panel Not Found")
            else:
                panel_color_classes, panel_text_classes, panel_display_name = panel_colors.get(
                    target_panel, ("bg-gray-100", "text-gray-800", f"{target_panel} Panel")
                )
            
            # Store panel information in state for use by other functions
            key = str(sample_dir)
            if key not in launcher._coverage_state:
                launcher._coverage_state[key] = {}
            launcher._coverage_state[key]["panel_display_name"] = panel_display_name

            with ui.row().classes(
                "w-full items-center justify-between mb-2 gap-2 flex-wrap"
            ):
                ui.label("Target Coverage").classes("text-lg font-semibold")
                target_coverage_back_button = ui.button(
                    "← Back to overview",
                    on_click=lambda: _show_target_coverage_overview(),
                ).props("flat dense no-caps outline").classes("hidden")

            with ui.row().classes("w-full items-center gap-2 mb-2"):
                ui.label("Panel:").classes("text-sm font-medium text-gray-600")
                ui.label(panel_display_name).classes(f"px-2 py-1 rounded text-sm font-medium {panel_color_classes} {panel_text_classes}")
                ui.label("•").classes("text-gray-400")
                ui.label("Target regions defined by gene panel").classes("text-xs text-gray-500")
            
            # Add detailed panel information in an expansion
            with ui.expansion().classes("w-full mb-2").props("icon=info dense"):
                ui.label("Panel Details").classes("text-sm font-medium mb-2")
                
                # Get the BED file name for the current panel
                # Use the same logic as target_analysis._find_target_bed to determine BED filename
                bed_filename = None
                if target_panel == "rCNS2":
                    bed_filename = "rCNS2_panel_name_uniq.bed"
                elif target_panel == "AML":
                    bed_filename = "AML_panel_name_uniq.bed"
                elif target_panel and target_panel.strip():
                    # For custom panels, construct the expected filename
                    bed_filename = f"{target_panel}_panel_name_uniq.bed"
                else:
                    bed_filename = "Unknown"
                
                # Verify the BED file exists in robin.resources
                if bed_filename and bed_filename != "Unknown":
                    try:
                        from robin import resources
                        resources_dir = os.path.dirname(os.path.abspath(resources.__file__))
                        bed_path = os.path.join(resources_dir, bed_filename)
                        if not os.path.exists(bed_path):
                            # File doesn't exist, but we'll still show the expected filename
                            # This helps users understand what file should be present
                            pass
                    except Exception:
                        # If we can't verify, still show the expected filename
                        pass
                
                with ui.column().classes("gap-1 text-sm"):
                    with ui.row().classes("items-center gap-2"):
                        ui.label("Panel:").classes("font-medium w-20")
                        ui.label(target_panel).classes("font-mono")
                    
                    with ui.row().classes("items-center gap-2"):
                        ui.label("BED File:").classes("font-medium w-20")
                        ui.label(bed_filename).classes("font-mono text-xs")
                    
                    with ui.row().classes("items-center gap-2"):
                        ui.label("Description:").classes("font-medium w-20")
                        panel_descriptions = {
                            "rCNS2": "Central Nervous System genes (244 regions)",
                            "AML": "Acute Myeloid Leukemia genes (1,181 regions)", 
                            "Sarcoma": "Sarcoma-specific gene panel",
                            "Unknown Panel": "Panel type could not be determined"
                        }
                        ui.label(panel_descriptions.get(target_panel, "Custom gene panel")).classes("text-xs")
            
            # Define helper functions before chart creation
            def _show_chromosome_scatter(chromosome: str) -> None:
                """Show scatter plot for individual genes on a specific chromosome"""
                try:
                    # Get the bed coverage data
                    bed_cov = sample_dir / "bed_coverage_main.csv"
                    if not bed_cov.exists():
                        ui.notify("No coverage data available", type="warning")
                        return
                    
                    # Read the data
                    df = pd.read_csv(bed_cov)
                    if "coverage" not in df.columns:
                        df["length"] = (df["endpos"] - df["startpos"] + 1).astype(float)
                        df["coverage"] = df["bases"] / df["length"]
                    
                    # Filter for the selected chromosome
                    chrom_data = df[df["chrom"].astype(str) == chromosome].copy()
                    
                    if chrom_data.empty:
                        ui.notify(f"No data found for chromosome {chromosome}", type="warning")
                        return
                    
                    # Sort by position for better visualization
                    chrom_data = chrom_data.sort_values("startpos")
                    
                    # Build unique display labels for genes that appear multiple times (same name, different regions)
                    name_counts = chrom_data["name"].value_counts()
                    scatter_data = []
                    for _, row in chrom_data.iterrows():
                        name = row["name"]
                        if name_counts.get(name, 0) > 1:
                            # Append position (in Mb) to distinguish duplicate gene regions
                            start_mb = row["startpos"] / 1e6
                            end_mb = row["endpos"] / 1e6
                            if abs(end_mb - start_mb) < 0.01:
                                label = f"{name} ({start_mb:.2f} Mb)"
                            else:
                                label = f"{name} ({start_mb:.2f}-{end_mb:.2f} Mb)"
                        else:
                            label = name
                        scatter_data.append([label, row["coverage"], row["startpos"], row["endpos"]])
                    
                    # Update chart to show scatter plot
                    target_boxplot.options["title"]["text"] = f"Gene Coverage - {chromosome}"
                    target_boxplot.options["title"]["subtext"] = f"{len(scatter_data)} targets"
                    
                    # Update x-axis to show gene names (with position when duplicated)
                    gene_names = [d[0] for d in scatter_data]
                    target_boxplot.options["xAxis"]["data"] = gene_names
                    target_boxplot.options["xAxis"]["axisLabel"]["rotate"] = 45
                    target_boxplot.options["xAxis"]["axisLabel"]["interval"] = 0
                    target_boxplot.options["xAxis"]["name"] = "Gene / region"
                    target_boxplot.options["xAxis"]["nameLocation"] = "middle"
                    target_boxplot.options["xAxis"]["nameGap"] = 44
                    
                    # Update series to show scatter plot
                    _cp_sc = _cov_chrome_palette()
                    target_boxplot.options["series"] = [
                        {
                            "name": "Gene Coverage",
                            "type": "scatter",
                            "id": "coverage_data",  # Add ID for universal transition
                            "data": scatter_data,
                            "symbolSize": 8,
                            "universalTransition": True,  # Enable universal transition
                            "animationDurationUpdate": 1000,  # Set transition duration
                            "itemStyle": {
                                "color": "#3b82f6",
                                "opacity": 0.7
                            },
                            "emphasis": {
                                "itemStyle": {
                                    "color": "#1d4ed8",
                                    "opacity": 1,
                                    "borderColor": _cp_sc["title"],
                                    "borderWidth": 2
                                }
                            },
                            "label": {
                                "show": True,
                                "position": "top",
                                "formatter": "{@[1]:.2f}x",
                                "fontSize": 10,
                                "color": _cp_sc["axis"],
                            },
                            "tooltip": {
                                ":formatter": "function(params) { const data = params.data; return 'Gene: ' + data[0] + '<br/>Coverage: ' + Number(data[1]).toFixed(2) + 'x<br/>Position: ' + data[2].toLocaleString() + '-' + data[3].toLocaleString(); }"
                            }
                        }
                    ]
                    
                    # Hide legend for scatter view
                    target_boxplot.options["legend"]["show"] = False
                    _apply_coverage_boxplot_chrome(target_boxplot)
                    target_boxplot.update()
                    
                    # Show back button (remove Tailwind hidden — classes("", replace="") is unreliable)
                    try:
                        target_coverage_back_button.classes(remove="hidden")
                    except Exception:
                        target_coverage_back_button.set_visibility(True)
                    
                    _log_notify(f"Showing gene coverage for {chromosome}", level="info", notify=False)
                    
                except Exception as e:
                    _log_notify(f"Failed to show chromosome scatter: {e}", level="error", notify=True)
            
            def _show_target_coverage_overview() -> None:
                """Return to the original box plot overview"""
                try:
                    # Reload the original box plot data by calling _update_boxplot
                    bed_cov = sample_dir / "bed_coverage_main.csv"
                    if bed_cov.exists():
                        df = pd.read_csv(bed_cov)
                        _update_boxplot(df, panel_display_name)
                    else:
                        ui.notify("No coverage data available", type="warning")
                        return
                    
                    # Hide back button
                    try:
                        target_coverage_back_button.classes("hidden")
                    except Exception:
                        target_coverage_back_button.set_visibility(False)
                    
                    _log_notify("Returned to target coverage overview", level="info", notify=False)
                    
                except Exception as e:
                    _log_notify(f"Failed to show overview: {e}", level="error", notify=True)
            
            # Define click handler function before chart creation
            def handle_boxplot_click(params):
                """Handle clicks on the ECharts box plot and show chromosome scatter"""
                try:
                    if params.series_name == 'box plot' and params.data:
                        chromosome = params.name
                        ui.notify(f"Showing coverage for chromosome: {chromosome}", type="info")
                        _log_notify(f"User clicked on chromosome: {chromosome}", level="info", notify=False)
                        _show_chromosome_scatter(chromosome)
                except Exception as e:
                    _log_notify(f"Error handling chart click: {e}", level="warning", notify=False)
            
            _cp_bp = _cov_chrome_palette()
            target_boxplot = ui.echart(
                {
                    "backgroundColor": "transparent",
                    "textStyle": {"color": _cp_bp["axis"]},
                    "title": {
                        "text": f"Target Coverage ({panel_display_name})",
                        "left": "center",
                        "top": 10,
                        "textStyle": {"fontSize": 16, "color": _cp_bp["title"]},
                    },
                    "grid": {
                        "left": "15%",
                        "right": "8%",
                        "bottom": "18%",
                        "top": "20%",
                        "containLabel": True,
                    },
                    "dataset": [
                        {
                            "id": "raw",
                            "dimensions": [
                                "chrom",
                                "min",
                                "Q1",
                                "median",
                                "Q3",
                                "max",
                                "chrom_index",
                            ],
                            "source": [],
                        },
                        {
                            "id": "rawdata",
                            "dimensions": ["chrom", "coverage", "name"],
                            "source": [],
                        },
                        {
                            "id": "outliers",
                            "dimensions": ["chrom", "coverage", "name"],
                            "source": [],
                        },
                        {
                            "id": "globaloutliers",
                            "dimensions": ["chrom", "coverage", "name"],
                            "source": [],
                        },
                    ],
                    "xAxis": {
                        "type": "category",
                        "name": "Chromosome",
                        "nameLocation": "middle",
                        "nameGap": 40,
                        "axisLabel": {"interval": 0, "rotate": 45},
                    },
                    "yAxis": {"type": "value", "name": "Coverage (x)"},
                    "legend": {
                        "left": 10,
                        "top": "center",
                        "orient": "vertical",
                        "itemGap": 10,
                        "selected": {
                            "box plot": True,
                            "outliers": True,
                            "global outliers": True,
                            "raw data": False,
                        },
                    },
                    "series": [
                        {
                            "id": "target_box_series",
                            "name": "box plot",
                            "type": "boxplot",
                            "itemStyle": _cov_boxplot_item_style(),
                            "datasetId": "raw",
                            "encode": {
                                "x": "chrom",
                                "y": ["min", "Q1", "median", "Q3", "max"],
                                "itemName": ["chrom"],
                                "tooltip": ["min", "Q1", "median", "Q3", "max"],
                            },
                        },
                        {
                            "id": "target_outliers_series",
                            "name": "outliers",
                            "type": "scatter",
                            "datasetId": "outliers",
                            "symbolSize": 6,
                            "itemStyle": _cov_scatter_outlier_style(),
                            "label": {
                                "show": True,
                                "position": "right",
                                "formatter": "{@name}",
                            },
                            "encode": {
                                "x": "chrom",
                                "y": "coverage",
                                "label": ["name"],
                                "tooltip": ["name", "coverage"],
                            },
                        },
                        {
                            "id": "target_global_outliers_series",
                            "name": "global outliers",
                            "type": "scatter",
                            "datasetId": "globaloutliers",
                            "symbolSize": 6,
                            "itemStyle": _cov_scatter_outlier_style(),
                            "label": {
                                "show": True,
                                "position": "right",
                                "formatter": "{@name}",
                            },
                            "encode": {
                                "x": "chrom",
                                "y": "coverage",
                                "label": ["name"],
                                "tooltip": ["name", "coverage"],
                            },
                        },
                    ],
                },
                on_point_click=handle_boxplot_click
            ).classes("w-full h-80 target-coverage-boxplot")
            _apply_coverage_boxplot_chrome(target_boxplot)
        with ui.card().classes("w-full"):
            with ui.row().classes("items-center gap-3"):
                target_search = ui.input("Search targets…").props(
                    "borderless dense clearable"
                )
            _, target_cov_table = styled_table(
                columns=[
                    {
                        "name": "chrom",
                        "label": "Chrom",
                        "field": "chrom",
                        "sortable": True,
                    },
                    {
                        "name": "startpos",
                        "label": "Start",
                        "field": "startpos",
                        "sortable": True,
                    },
                    {
                        "name": "endpos",
                        "label": "End",
                        "field": "endpos",
                        "sortable": True,
                    },
                    {
                        "name": "name",
                        "label": "Name",
                        "field": "name",
                        "sortable": True,
                    },
                    {
                        "name": "coverage",
                        "label": "Coverage (x)",
                        "field": "coverage",
                        "sortable": True,
                    },
                    {
                        "name": "outlier_status",
                        "label": "Outlier (mean ± 2 SD)",
                        "field": "outlier_status",
                        "sortable": True,
                    },
                ],
                rows=[],
                pagination=20,
                class_size="table-xs",
            )
            try:
                target_cov_table.props(
                    'multi-sort rows-per-page-options="[10,20,50,0]"'
                )
                target_search.bind_value(target_cov_table, "filter")
            except Exception:
                pass
            target_cov_table.add_slot(
                "body-cell-coverage",
                """
    <q-td key="coverage" :props="props">
    <q-badge :color="props.value >= 30 ? 'green' : props.value >= 20 ? 'blue' : props.value >= 10 ? 'orange' : 'red'">
        {{ Number(props.value).toFixed(2) }}
    </q-badge>
    </q-td>
    """,
            )
            target_cov_table.add_slot(
                "body-cell-outlier_status",
                """
    <q-td key="outlier_status" :props="props">
    <span v-if="props.value">{{ props.value }}</span>
    <span v-else class="text-gray-400">—</span>
    </q-td>
    """,
            )

            # IGV viewer has been moved to the details page
            # Use add_igv_viewer() function to add it to other pages
            if False:  # Disabled - IGV viewer moved to details page
                with ui.card().classes("w-full"):
                    ui.label("IGV").classes("text-lg font-semibold mb-2")
                    igv_div = ui.element("div").classes("w-full h-[600px] border")
                    igv_div._props["id"] = "igv-container"
                    igv_status = ui.label("Checking for IGV-ready BAM files...").classes(
                        "text-sm text-gray-600"
                    )

                    # Add IGV library status indicator
                    igv_lib_status = ui.label("IGV library: Checking...").classes(
                        "text-xs text-gray-500"
                    )

                    # IGV state management - prevent unnecessary redraws
                    def _get_igv_state():
                        try:
                            key = str(sample_dir)
                            if hasattr(launcher, "_coverage_state"):
                                if key not in launcher._coverage_state:
                                    launcher._coverage_state[key] = {}
                                return launcher._coverage_state[key]
                            else:
                                launcher._coverage_state = {}
                                launcher._coverage_state[key] = {}
                                return launcher._coverage_state[key]
                        except Exception:
                            # Return empty state as fallback
                            return {}

                def _is_igv_ready():
                    """Check if IGV browser is properly initialized and ready."""
                    try:
                        state = _get_igv_state()
                        return (
                            state.get("igv_initialized", False)
                            and state.get("igv_browser_ready", False)
                            and state.get("igv_loaded_bam") is not None
                        )
                    except Exception:
                        return False

                def _set_igv_ready(bam_url: str):
                    """Mark IGV as ready with the current BAM."""
                    state = _get_igv_state()
                    state["igv_initialized"] = True
                    state["igv_browser_ready"] = True
                    state["igv_loaded_bam"] = bam_url

                def _clear_igv_state():
                    """Reset IGV state to force reinitialization."""
                    state = _get_igv_state()
                    state["igv_initialized"] = False
                    state["igv_browser_ready"] = False
                    state["igv_loaded_bam"] = None
                    state["igv_loading"] = False

                # Immediate check for existing IGV BAM files
                def _check_existing_igv_bam():
                    try:
                        # Only proceed if IGV isn't already ready
                        if _is_igv_ready():
                            igv_status.set_text("IGV is already loaded and ready.")
                            return

                        # Determine candidate BAMs (prefer standardized IGV path)
                        igv_dir = sample_dir / "igv"
                        clair_dir = sample_dir / "clair3"
                        candidates = [
                            igv_dir / "igv_ready.bam",
                            clair_dir / "sorted_targets_exceeding_rerun.bam",
                            clair_dir / "sorted_targets_exceeding.bam",
                            sample_dir / "target.bam",
                        ]

                        bam_path = None
                        for p in candidates:
                            bai = Path(f"{p}.bai")
                            if p.exists() and bai.exists():
                                bam_path = p
                                break

                        if bam_path is not None:
                            igv_status.set_text(f"Found existing BAM: {bam_path.name}")
                            # Trigger immediate IGV loading
                            ui.timer(0.1, lambda: _load_igv_bam(bam_path), once=True)
                        else:
                            igv_status.set_text(
                                "No IGV-ready BAM found yet (need coordinate-sorted BAM with .bai)."
                            )

                    except Exception as e:
                        igv_status.set_text(f"Error checking for BAM files: {e}")

                # Function to check if IGV library is available
                def _check_igv_library():
                    """Check if the IGV JavaScript library is available."""
                    try:
                        js_check = """
                            try {
                                if (typeof igv === 'undefined') {
                                    console.error('IGV library not loaded');
                                    return false;
                                }
                        
                                // Check if IGV is fully initialized
                                if (typeof igv.createBrowser !== 'function') {
                                    console.error('IGV library not fully initialized');
                                    return false;
                                }
                        
                                console.log('IGV library is available and ready');
                                return true;
                            } catch (e) {
                                console.error('Error checking IGV library:', e);
                                return false;
                            }
                        """

                        try:
                            result = ui.run_javascript(js_check, timeout=10.0)
                            # Update the status indicator
                            if result:
                                igv_lib_status.set_text("IGV library: ✓ Loaded and ready")
                                igv_lib_status.classes("text-xs text-green-600")
                            else:
                                igv_lib_status.set_text("IGV library: ✗ Not ready")
                                igv_lib_status.classes("text-xs text-red-600")
                            return result
                        except Exception as e:
                            igv_lib_status.set_text("IGV library: ✗ Error checking")
                            igv_lib_status.classes("text-xs text-red-600")
                            print(f"Error checking IGV library: {e}")
                            return False
                    except Exception:
                        return False

                # Function to retry IGV creation
                def _retry_igv_creation(bam_path: Path):
                    """Retry IGV browser creation after a failure."""
                    try:
                        igv_status.set_text("Retrying IGV browser creation...")
                        _load_igv_bam(bam_path)
                    except Exception as e:
                        igv_status.set_text(f"Retry failed: {e}")
                        print(f"IGV retry error: {e}")

                # Function to wait for element to be ready
                def _wait_for_element_ready():
                    """Wait for the IGV div element to be properly rendered."""
                    js_wait = """
                    return new Promise((resolve, reject) => {
                        const el = document.getElementById('igv-container');
                        if (!el) {
                            reject(new Error('Element not found'));
                            return;
                        }
                
                        if (el.offsetWidth > 0 && el.offsetHeight > 0) {
                            resolve(true);
                            return;
                        }
                
                        // Wait for element to be ready
                        const checkReady = () => {
                            if (el.offsetWidth > 0 && el.offsetHeight > 0) {
                                resolve(true);
                            } else {
                                setTimeout(checkReady, 100);
                            }
                        };
                        checkReady();
                    });
                    """
                    try:
                        return ui.run_javascript(js_wait, timeout=10.0)
                    except Exception:
                        return False

                # Function to load BAM into IGV
                def _load_igv_bam(bam_path: Path):
                    try:
                        # Get state first
                        state = _get_igv_state()

                        # Prevent multiple simultaneous IGV loading attempts
                        if state.get("igv_loading", False):
                            igv_status.set_text("IGV is already being loaded, please wait...")
                            return

                        # First check if IGV library is available
                        if not _check_igv_library():
                            igv_status.set_text(
                                "IGV library not yet loaded. Please wait and try again."
                            )
                            return

                        # Check if we already have this BAM loaded
                        bam_url = f"/samples/{sample_dir.name}/{bam_path.name}"
                        if state.get("igv_loaded_bam") == bam_url and _is_igv_ready():
                            igv_status.set_text(f"BAM {bam_path.name} already loaded in IGV.")
                            return

                        # Mark that we're loading IGV
                        state["igv_loading"] = True

                        # Mount directory once
                        if bam_path.parent == sample_dir / "igv":
                            mount = f"/samples/{sample_dir.name}/igv"
                        elif bam_path.parent == sample_dir / "clair3":
                            mount = f"/samples/{sample_dir.name}/clair3"
                        else:
                            mount = f"/samples/{sample_dir.name}"

                        try:
                            if app is not None:
                                app.add_static_files(mount, str(bam_path.parent))
                        except Exception:
                            # Already mounted or in tests without app
                            pass

                        # Define URLs
                        bam_url = f"{mount}/{bam_path.name}"
                        bai_url = f"{mount}/{bam_path.name}.bai"

                        # Check if we need to create a new browser or just add tracks
                        if not state.get("igv_initialized") or not state.get(
                            "igv_browser_ready"
                        ):
                            # Wait for element to be ready before creating IGV browser
                            if not _wait_for_element_ready():
                                igv_status.set_text("Waiting for IGV element to be ready...")
                                ui.timer(0.5, lambda: _load_igv_bam(bam_path), once=True)
                                return

                            # Create new IGV browser
                            js_create = """
                                try {
                                    console.log('Creating minimal IGV browser...');
                                    const el = document.getElementById('igv-container');
                                    if (!el) {
                                    throw new Error('Target element not found');
                                    }
                            
                                    const options = { genome: 'hg38' };
                                    console.log('IGV options:', options);
                            
                                    igv.createBrowser(el, options).then(b => {
                                    window.lj_igv = b; 
                                    window.lj_igv_browser_ready = true;
                                    console.log('IGV browser created successfully');
                                    }).catch(error => {
                                    console.error('IGV browser creation failed:', error);
                                    window.lj_igv_browser_ready = false;
                                    });
                                } catch (error) {
                                    console.error('Error in IGV browser creation:', error);
                                }
                            """

                            try:
                                ui.run_javascript(js_create, timeout=30.0)
                                _set_igv_ready(bam_url)
                                igv_status.set_text(f"IGV browser created with {bam_path.name}")
                                # Clear loading flag on success
                                state["igv_loading"] = False
                            except Exception as e:
                                igv_status.set_text(f"Failed to create IGV browser: {e}")
                                print(f"IGV browser creation error: {e}")
                                # Clear loading flag on failure
                                state["igv_loading"] = False

                                # Try to retry after a delay
                                ui.timer(2.0, lambda: _retry_igv_creation(bam_path), once=True)
                                _clear_igv_state()
                        else:
                            # Browser exists, just add/update the track
                            js_add_track = f"""
                                if (window.lj_igv && window.lj_igv_browser_ready) {{
                                    console.log('Adding track to existing IGV browser...');
                            
                                    const track = {{ name: '{sample_dir.name}', url: '{bam_url}', indexURL: '{bai_url}', format: 'bam', type: 'alignment', height: 600, autoScale: true, colorBy: 'tag', tag: 'SA' }};
                                    console.log('Track to load:', track);
                            
                                    try {{
                                    window.lj_igv.loadTrack(track);
                                    console.log('Track loaded successfully');
                                    }} catch (error) {{
                                    console.error('Error loading track:', error);
                                    }}
                                }}
                            """

                            try:
                                ui.run_javascript(js_add_track, timeout=30.0)
                                _set_igv_ready(bam_url)
                                igv_status.set_text(f"Track updated in IGV: {bam_path.name}")
                                # Clear loading flag on success
                                state["igv_loading"] = False
                            except Exception as e:
                                igv_status.set_text(f"Failed to update track: {e}")
                                # Clear loading flag on failure
                                state["igv_loading"] = False

                    except Exception as e:
                        igv_status.set_text(f"Error loading IGV: {e}")
                        _clear_igv_state()
                        # Clear loading flag on error
                        if "state" in locals():
                            state["igv_loading"] = False

                    # Check for existing files once after a delay to ensure IGV library is loaded
                    ui.timer(2.0, _check_existing_igv_bam, once=True)

                # Function to refresh IGV BAM check (useful after regeneration)
                def _refresh_igv_check():
                    try:
                        # Only refresh if IGV isn't ready or if we need to check for new files
                        if _is_igv_ready():
                            # Check if the current BAM still exists
                            current_bam = _get_igv_state().get("igv_loaded_bam")
                            if current_bam:
                                # Extract filename from URL and check if it still exists
                                bam_name = current_bam.split("/")[-1]
                                igv_dir = sample_dir / "igv"
                                clair_dir = sample_dir / "clair3"
                                candidates = [
                                    igv_dir / bam_name,
                                    clair_dir / bam_name,
                                    sample_dir / bam_name,
                                ]

                                if any(
                                    p.exists() and Path(f"{p}.bai").exists() for p in candidates
                                ):
                                    igv_status.set_text("IGV is ready and BAM file is current.")
                                    return

                        # If we get here, we need to refresh
                        _clear_igv_state()
                        _check_existing_igv_bam()
                    except Exception as e:
                        print(f"Error refreshing IGV check: {e}")

                # Function to clear data tracks from IGV (keeps reference tracks)
                def _clear_igv_tracks():
                    # Initialize state variable
                    state = None
                    try:
                        # First check if IGV is actually ready
                        if not _is_igv_ready():
                            igv_status.set_text("IGV is not ready - cannot clear tracks.")
                            return

                        # Simple JavaScript to clear tracks
                        js_clear = """
                            try {
                                console.log('Attempting to clear IGV tracks...');
                        
                                if (window.lj_igv && typeof window.lj_igv.getTracks === 'function') {
                                    const tracks = window.lj_igv.getTracks();
                                    console.log('Found tracks:', tracks.length);
                            
                                    if (tracks.length > 0) {
                                    tracks.forEach((track, index) => {
                                        try {
                                            if (typeof window.lj_igv.removeTrack === 'function') {
                                                window.lj_igv.removeTrack(track);
                                                console.log('Removed track:', index);
                                            }
                                        } catch (e) {
                                            console.error('Error removing track:', e);
                                        }
                                    });
                                    console.log('Track clearing completed');
                                    }
                                } else {
                                    console.log('No IGV browser or getTracks method available');
                                }
                            } catch (e) {
                                console.error('Error in track clearing:', e);
                            }
                        """

                        ui.run_javascript(js_clear, timeout=15.0)

                        # Update status and state
                        igv_status.set_text(
                            "IGV data tracks cleared. Reference tracks preserved."
                        )

                        # Don't clear the IGV state - just mark that data tracks are cleared
                        state = _get_igv_state()
                        if state:
                            state["data_tracks_cleared"] = True

                    except Exception as e:
                        igv_status.set_text(f"Error clearing tracks: {e}")
                        print(f"Error in _clear_igv_tracks: {e}")

                # Function to wait for BAM file to stabilize (not being written to)
                def _wait_for_bam_ready(bam_path: Path, max_wait_time: int = 30) -> bool:
                    """
                    Wait for BAM file to be ready (not actively being written to).
                    Returns True if file is ready, False if timeout.
                    """
                    import time
            
                    start_time = time.time()
                    last_size = -1
                    stable_count = 0
                    required_stable_checks = 2  # Need 2 consecutive checks with same size
            
                    while time.time() - start_time < max_wait_time:
                        try:
                            if not bam_path.exists():
                                time.sleep(0.5)
                                continue
                    
                            current_size = bam_path.stat().st_size
                    
                            # Check if file size has changed
                            if current_size == last_size:
                                stable_count += 1
                                if stable_count >= required_stable_checks:
                                    # File size is stable for consecutive checks
                                    return True
                            else:
                                stable_count = 0  # Reset counter if size changed
                    
                            last_size = current_size
                            time.sleep(0.5)
                    
                        except (OSError, IOError):
                            # File might be temporarily unavailable
                            time.sleep(0.5)
                            continue
            
                    # Timeout reached
                    return False
        
                # Function to reload BAM data into existing IGV browser
                def _reload_bam_track():
                    # Initialize state variable
                    state = None
                    try:
                        # First check if IGV is actually ready
                        if not _is_igv_ready():
                            igv_status.set_text("IGV is not ready - cannot reload BAM.")
                            return

                        # Find the current BAM file
                        state = _get_igv_state()
                        current_bam_url = state.get("igv_loaded_bam")

                        if not current_bam_url:
                            igv_status.set_text("No BAM file currently loaded - cannot reload.")
                            return

                        # Extract BAM name for display
                        bam_name = current_bam_url.split("/")[-1]
                
                        # Determine the actual BAM file path
                        bam_path = None
                        if current_bam_url.startswith("/samples/"):
                            # Extract sample directory and BAM name
                            parts = current_bam_url.split("/")
                            if len(parts) >= 4:
                                sample_dir_name = parts[2]
                                bam_file_name = parts[3]
                                bam_path = sample_dir / bam_file_name
                
                        # Wait for BAM file to be ready if it exists
                        if bam_path and bam_path.exists():
                            igv_status.set_text("Checking if BAM file is ready...")
                            if not _wait_for_bam_ready(bam_path):
                                igv_status.set_text("BAM file is still being updated. Please wait and try again.")
                                ui.notify("BAM file is still being updated. Please wait and try again.", type="warning")
                                return
                    
                            # Also check if BAI file is ready
                            bai_path = bam_path.with_suffix(bam_path.suffix + ".bai")
                            if bai_path.exists():
                                if not _wait_for_bam_ready(bai_path):
                                    igv_status.set_text("BAM index is still being updated. Please wait and try again.")
                                    ui.notify("BAM index is still being updated. Please wait and try again.", type="warning")
                                    return

                        # JavaScript to reload the BAM track
                        js_reload = f"""
                            try {{
                                                        if (window.lj_igv) {{
                                    // Try to get tracks list
                                    let tracks = [];
                                    try {{
                                    // Try multiple methods to get the tracks list
                                    if (typeof window.lj_igv.getTracks === 'function') {{
                                        tracks = window.lj_igv.getTracks();
                                    }} else if (window.lj_igv.trackList && Array.isArray(window.lj_igv.trackList)) {{
                                        tracks = window.lj_igv.trackList;
                                    }} else if (window.lj_igv.tracks && Array.isArray(window.lj_igv.tracks)) {{
                                        tracks = window.lj_igv.tracks;
                                    }} else if (window.lj_igv.roster && Array.isArray(window.lj_igv.roster.trackViews)) {{
                                        tracks = window.lj_igv.roster.trackViews;
                                    }} else {{
                                        // Try to remove by name directly
                                        if (typeof window.lj_igv.removeTrackByName === 'function') {{
                                            try {{
                                                window.lj_igv.removeTrackByName('{bam_name}');
                                            }} catch (e) {{
                                                // Ignore errors
                                            }}
                                        }}
                                    }}
                                
                                    // Remove alignment tracks
                                    for (let track of tracks) {{
                                        try {{
                                            // Only remove BAM/alignment tracks, not reference tracks
                                            if (track && (track.format === 'bam' || track.format === 'cram' || track.type === 'alignment')) {{
                                                // Try multiple removal methods
                                                if (typeof window.lj_igv.removeTrack === 'function') {{
                                                    window.lj_igv.removeTrack(track);
                                                }} else if (typeof window.lj_igv.removeTrackByName === 'function') {{
                                                    window.lj_igv.removeTrackByName(track.name);
                                                }}
                                            }}
                                        }} catch (e) {{
                                            // Ignore errors
                                        }}
                                    }}
                                    }} catch (e) {{
                                    // Ignore errors
                                    }}
                            
                                    // Create the track configuration
                                    const trackConfig = {{
                                    name: '{bam_name}',
                                    url: '{current_bam_url}',
                                    indexURL: '{current_bam_url}.bai',
                                    format: 'bam',
                                    visibilityWindow: 1000000,
                                    type: 'alignment',
                                    height: 600,
                                    autoScale: true
                                    }};
                            
                                    // Add the track to the browser
                                    if (typeof window.lj_igv.loadTrack === 'function') {{
                                    window.lj_igv.loadTrack(trackConfig);
                                    }} else {{
                                    // Try alternative method
                                    if (typeof window.lj_igv.addTrack === 'function') {{
                                        window.lj_igv.addTrack(trackConfig);
                                    }}
                                    }}
                            
                                    // Force a redraw
                                    if (typeof window.lj_igv.redraw === 'function') {{
                                    window.lj_igv.redraw();
                                    }} else if (typeof window.lj_igv.update === 'function') {{
                                    window.lj_igv.update();
                                    }}
                                }}
                            }} catch (e) {{
                                // Ignore errors
                            }}
                        """

                        ui.run_javascript(js_reload, timeout=20.0)

                        # Update status
                        igv_status.set_text(f"Loading BAM track: {bam_name}")
                
                        # Change button text to "Reload BAM" after first load
                        reload_bam_button.set_text("Reload BAM")

                        # Mark that data tracks are no longer cleared
                        if state:
                            state["data_tracks_cleared"] = False
                    
                        # Auto-load the BED file after BAM is loaded
                        _load_target_bed()

                    except Exception as e:
                        igv_status.set_text(f"Error reloading BAM: {e}")
                        print(f"Error in _reload_bam_track: {e}")

                # Function to load target BED file into IGV
                def _load_target_bed():
                    """Load the target BED file as a track in IGV"""
                    try:
                        # First check if IGV is actually ready
                        if not _is_igv_ready():
                            igv_status.set_text("IGV is not ready - cannot load BED file.")
                            return

                        # Get the target panel information
                        state = _get_igv_state()
                
                        # Check if BED file has already been loaded
                        if state.get("bed_file_loaded", False):
                            return  # BED file already loaded, don't load again
                
                        # Try to read the panel from master.csv
                        target_panel = None
                        bed_file_path = None
                        try:
                            master_csv_path = sample_dir / "master.csv"
                            if master_csv_path.exists():
                                import pandas as pd
                                df = pd.read_csv(master_csv_path)
                                if not df.empty and "analysis_panel" in df.columns:
                                    target_panel = str(df.iloc[0]["analysis_panel"]).strip()
                            
                                    # Map panel to BED filename
                                    bed_file_mapping = {
                                    "rCNS2": "rCNS2_panel_name_uniq.bed",
                                    "AML": "AML_panel_name_uniq.bed",
                                    "Sarcoma": "Sarcoma_panel_name_uniq.bed"
                                    }
                            
                                    bed_filename = bed_file_mapping.get(target_panel, f"{target_panel}_panel_name_uniq.bed")
                            
                                    # Try to find the BED file in robin resources
                                    try:
                                        from robin import resources
                                        bed_file_path = os.path.join(
                                            os.path.dirname(os.path.abspath(resources.__file__)),
                                            bed_filename
                                        )
                                        if not os.path.exists(bed_file_path):
                                            bed_file_path = None
                                    except Exception:
                                        pass
                            
                                    # Fallback paths
                                    if not bed_file_path:
                                        possible_paths = [
                                        bed_filename,
                                        f"data/{bed_filename}",
                                        f"/usr/local/share/{bed_filename}",
                                    ]
                                    for path in possible_paths:
                                        if os.path.exists(path):
                                            bed_file_path = path
                                            break
                        except Exception as e:
                            print(f"Error reading panel information: {e}")

                        if not bed_file_path or not os.path.exists(bed_file_path):
                            igv_status.set_text(f"Could not find BED file for panel: {target_panel}")
                            ui.notify(f"Target BED file not found for panel: {target_panel}", type="warning")
                            return

                        # Mount the BED file's directory
                        bed_dir = os.path.dirname(bed_file_path)
                        bed_name = os.path.basename(bed_file_path)
                        bed_url = f"/robin_resources/{bed_name}"
                
                        try:
                            if app is not None:
                                app.add_static_files("/robin_resources", bed_dir)
                        except Exception:
                            # Already mounted or in tests without app
                            pass

                        # JavaScript to load the BED file as a track
                        js_load_bed = f"""
                            try {{
                                if (window.lj_igv) {{
                                    // Create the track configuration for the BED file
                                    const trackConfig = {{
                                    name: '{bed_name}',
                                    url: '{bed_url}',
                                    format: 'bed',
                                    indexed: false,
                                    displayMode: 'COLLAPSED',
                                    height: 40,
                                    color: '#0072B2'
                                    }};
                            
                                    // Load the track into IGV
                                    if (typeof window.lj_igv.loadTrack === 'function') {{
                                    window.lj_igv.loadTrack(trackConfig);
                                    }} else {{
                                    // Try alternative method
                                    if (typeof window.lj_igv.addTrack === 'function') {{
                                        window.lj_igv.addTrack(trackConfig);
                                    }}
                                    }}
                                }}
                            }} catch (e) {{
                                // Ignore errors
                            }}
                        """

                        ui.run_javascript(js_load_bed, timeout=20.0)
                
                        # Mark that BED file has been loaded
                        state["bed_file_loaded"] = True
                
                        # Update status
                        igv_status.set_text(f"Loaded BED file: {bed_name}")
                        ui.notify(f"Target BED file loaded: {bed_name}", type="positive")

                    except Exception as e:
                        igv_status.set_text(f"Error loading BED: {e}")
                        print(f"Error in _load_target_bed: {e}")

                # Function to check console for errors
                def _check_console_errors():
                    """Check the browser console for any JavaScript errors."""
                    js_check_console = """
                    try {
                        const errors = [];

                        // Check for basic IGV functionality
                        if (typeof igv === 'undefined') {
                            errors.push('IGV library not loaded');
                        } else if (typeof igv.createBrowser !== 'function') {
                            errors.push('IGV createBrowser method not available');
                        } else {
                            errors.push(`IGV version: ${igv.version || 'unknown'}`);
                        }
                
                        return errors;
                    } catch (e) {
                        return ['Error checking console: ' + e.message];
                    }
                    """

                    try:
                        errors = ui.run_javascript(js_check_console, timeout=10.0)
                        if errors and len(errors) > 0:
                            error_text = "Console errors found:\n" + "\n".join(errors)
                            ui.notify(error_text, type="warning")
                            print(error_text)
                        else:
                            ui.notify("No console errors found", type="positive")
                    except Exception as e:
                        ui.notify(f"Error checking console: {e}", type="negative")

                # Function to create an empty IGV browser
                def _create_empty_igv():
                    """Create an empty IGV browser with no tracks."""
                    try:
                        if not _check_igv_library():
                            ui.notify("IGV library not available", type="warning")
                            return

                        # First, let's try to find the element and wait for it to be ready
                        js_find_element = """
                            console.log('Looking for IGV div element...');
                            console.log('Expected ID: igv-container');
                    
                            // List all divs on the page
                            const allDivs = document.querySelectorAll('div');
                            console.log('Total divs found:', allDivs.length);
                    
                            allDivs.forEach((div, i) => {
                                console.log(`Div ${i}:`, {
                                    id: div.id,
                                    className: div.className,
                                    tagName: div.tagName,
                                    offsetWidth: div.offsetWidth,
                                    offsetHeight: div.offsetHeight
                                });
                            });
                    
                            // Try to find our specific element
                            let el = document.getElementById('igv-container');
                            if (!el) {
                                console.log('Element not found by ID, trying by class...');
                                const candidates = document.querySelectorAll('.w-full.h-\\[600px\\].border');
                                console.log('Found candidates by class:', candidates.length);
                                if (candidates.length > 0) {
                                    el = candidates[0];
                                    console.log('Using first candidate by class');
                                }
                            }
                    
                            if (el) {
                                console.log('Element found:', el);
                                console.log('Element dimensions:', el.offsetWidth, 'x', el.offsetHeight);
                                return true;
                            } else {
                                console.log('No suitable element found');
                                return false;
                            }
                        """

                        # Check if element exists first
                        element_found = ui.run_javascript(js_find_element, timeout=10.0)
                        if not element_found:
                            ui.notify(
                                "IGV div element not found. Please wait for the page to load completely.",
                                type="warning",
                            )
                            return

                        # Now create the IGV browser
                        js_empty = """
                            try {
                                console.log('Creating empty IGV browser...');
                        
                                // Find the element again
                                let el = document.getElementById('igv-container');
                                if (!el) {
                                    el = document.querySelector('.w-full.h-\\[600px\\].border');
                                }
                        
                                if (!el) {
                                    throw new Error('Element still not found');
                                }
                        
                                const options = { genome: 'hg38' };
                                console.log('Creating empty IGV browser with options:', options);
                        
                                igv.createBrowser(el, options).then(b => {
                                    window.lj_igv = b; 
                                    window.lj_igv_browser_ready = true;
                                    console.log('Empty IGV browser created successfully');
                                }).catch(error => {
                                    console.error('Empty IGV browser creation failed:', error);
                                    window.lj_igv_browser_ready = false;
                                });
                            } catch (error) {
                                console.error('Error creating empty IGV:', error);
                            }
                        """

                        ui.run_javascript(js_empty, timeout=30.0)
                        ui.notify("Empty IGV browser creation initiated", type="info")

                    except Exception as e:
                        ui.notify(f"Error creating empty IGV: {e}", type="negative")
                        print(f"Error in empty IGV: {e}")

                # Function to debug IGV state
                def _debug_igv_state():
                    try:
                        js_debug = """
                        console.log('=== IGV Debug Info ===');
                        console.log('window.lj_igv:', window.lj_igv);
                        console.log('window.lj_igv_browser_ready:', window.lj_igv_browser_ready);
                
                        if (window.lj_igv) {
                            console.log('IGV browser exists');
                            console.log('Available methods:', Object.getOwnPropertyNames(window.lj_igv));
                        } else {
                            console.log('No IGV browser found');
                        }
                
                        console.log('Python state check requested');
                        """

                        ui.run_javascript(js_debug, timeout=10.0)

                        # Also show Python state
                        state = _get_igv_state()
                        debug_info = f"Python IGV State: initialized={state.get('igv_initialized')}, browser_ready={state.get('igv_browser_ready')}, bam={state.get('igv_loaded_bam')}, tracks_cleared={state.get('tracks_cleared')}"
                        igv_status.set_text(debug_info)

                    except Exception as e:
                        igv_status.set_text(f"Debug error: {e}")
                        print(f"Error in debug: {e}")

                # IGV control buttons
                reload_bam_button = ui.button("Load BAM", on_click=_reload_bam_track)

                # Add IGV library status check timer once after a delay
                ui.timer(3.0, _check_igv_library, once=True)

                # BAM generation buttons
                async def _trigger_build_sorted_bam(force_regenerate: bool = False) -> None:
                    try:
                        # Debug: Check what's in the launcher
                        debug_info = f"launcher type: {type(launcher)}, workflow_runner: {getattr(launcher, 'workflow_runner', 'None')}"

                        if force_regenerate:
                            ui.notify(
                                "Force regenerate mode: will recreate IGV BAM even if it exists",
                                type="info",
                            )

                        # Check if target.bam exists (required for IGV BAM generation)
                        target_bam = sample_dir / "target.bam"
                        if not target_bam.exists():
                            ui.notify(
                                "No target.bam file found. Please run target analysis first.",
                                type="warning",
                            )
                            return

                        # Check if IGV BAM already exists (only for normal generation, not force regenerate)
                        igv_bam = sample_dir / "igv" / "igv_ready.bam"
                        if (
                            not force_regenerate
                            and igv_bam.exists()
                            and (sample_dir / "igv" / "igv_ready.bam.bai").exists()
                        ):
                            ui.notify("IGV BAM already exists and is ready.", type="positive")
                            return

                        if force_regenerate:
                            ui.notify(
                                "Proceeding with force regenerate - existing files will be overwritten",
                                type="info",
                            )

                        # Check if we have a workflow runner
                        if (
                            not hasattr(launcher, "workflow_runner")
                            or not launcher.workflow_runner
                        ):
                            ui.notify(
                                f"No workflow runner available. Debug: {debug_info}. GUI is running but workflow integration is not available. This may be a configuration issue.",
                                type="warning",
                            )
                            return

                        # Submit the IGV BAM generation job to the workflow system
                        try:
                            runner = launcher.workflow_runner
                            sample_id = sample_dir.name

                            # Check if it's a Ray workflow or Simple workflow
                            if hasattr(runner, "submit_sample_job"):
                                # Simple workflow
                                success = runner.submit_sample_job(
                                    str(sample_dir), "igv_bam", sample_id, force_regenerate
                                )
                            elif hasattr(runner, "manager") and hasattr(
                                runner.manager, "submit_sample_job"
                            ):
                                # Ray workflow - need to get the coordinator
                                # This is more complex for Ray, so we'll provide a fallback
                                ui.notify(
                                    "Ray workflow detected. IGV BAM generation should happen automatically "
                                    "after target analysis completes. If you need to regenerate it, "
                                    "please restart the workflow or check the logs.",
                                    type="info",
                                )
                                return
                            else:
                                ui.notify(
                                    "Unknown workflow type. Cannot submit IGV BAM job.",
                                    type="warning",
                                )
                                return

                            if success:
                                action = "regenerated" if force_regenerate else "generated"
                                ui.notify(
                                    f"IGV BAM {action} job submitted to workflow queue!",
                                    type="positive",
                                )
                                # Refresh the IGV check after successful job submission
                                _refresh_igv_check()
                            else:
                                ui.notify(
                                    "Failed to submit IGV BAM generation job. Check logs for details.",
                                    type="warning",
                                )

                        except Exception as e:
                            ui.notify(f"Error submitting IGV BAM job: {e}", type="negative")

                    except Exception as e:
                        try:
                            ui.notify(f"Error checking IGV BAM status: {e}", type="negative")
                        except Exception:
                            # Client may have disconnected, ignore UI errors
                            pass



    # SNP Analysis section (only shown in development mode)
    if is_development_mode:
        with ui.card().classes("w-full"):
            ui.label("SNP Analysis").classes("text-lg font-semibold mb-2")

            # SNP Analysis controls and status
            with ui.row().classes("w-full items-center justify-between mb-4"):
                with ui.column().classes("gap-2"):
                    ui.label("Variant Calling").classes("text-sm font-medium")
                    snp_status_label = ui.label("Ready to run SNP analysis").classes(
                        "text-sm text-gray-600"
                    )

                    # Add file availability status
                    snp_files_status = ui.label("Checking required files...").classes(
                        "text-xs text-gray-500"
                    )

                with ui.column().classes("gap-2 items-end"):
                    # SNP Analysis button
                    snp_analysis_button = ui.button(
                        "Run SNP Analysis",
                        icon="play_arrow",
                        on_click=lambda: _trigger_snp_analysis(),
                    ).props("color=primary data-snp-analysis-button")

                    # Force regenerate checkbox
                    force_regenerate_checkbox = ui.checkbox(
                        "Force regenerate", value=False
                    ).props("dense")

            # SNP Analysis requirements check
            with ui.expansion().classes("w-full").props("icon=info"):
                ui.label("Requirements").classes("text-sm font-medium mb-2")
                with ui.column().classes("gap-1 text-sm"):
                    ui.label(
                        "• target.bam - Target regions BAM file (generated by target analysis)"
                    )
                    ui.label(
                        "• targets_exceeding_threshold.bed - BED file defining regions (auto-generated by target analysis)"
                    )
                    ui.label(
                        "• Reference genome (hg38 recommended) - provided via CLI --reference option"
                    )
                    ui.label("• Docker running with hkubal/clairs-to:latest image")
                    ui.label("• snpEff and SnpSift installed")

                # Add helpful note about file generation
                with ui.row().classes("w-full mt-2 p-2 bg-blue-50 rounded"):
                    ui.icon("info").classes("text-blue-600")
                    ui.label(
                        "Note: Both target.bam and targets_exceeding_threshold.bed are automatically generated when you run target analysis. You don't need to create them manually."
                    ).classes("text-sm text-blue-800")

            # Lazy-load SNP/INDEL UI when Results expansion opens (avoids full VCF parse on page load)
            snp_results_lazy_state: Dict[str, Any] = {"sig": None}

            def _snp_results_file_signature() -> Tuple[Optional[int], Optional[int]]:
                clair_dir = sample_dir / "clair3"
                snp_vcf = clair_dir / "snpsift_output.vcf"
                indel_vcf = clair_dir / "snpsift_indel_output.vcf"

                def _mt(p: Path) -> Optional[int]:
                    try:
                        return p.stat().st_mtime_ns if p.exists() else None
                    except OSError:
                        return None

                return (_mt(snp_vcf), _mt(indel_vcf))

            snp_results_expansion = ui.expansion("Results", icon="assessment").classes(
                "w-full"
            )
            with snp_results_expansion:
                snp_results_status = (
                    ui.label("Expand this section to load variant results.")
                    .classes("text-sm text-gray-600")
                    .props("data-snp-results-status")
                )
                snp_results_container = ui.column().classes("w-full")

            def _check_snp_results():
                try:
                    clair_dir = sample_dir / "clair3"
                    snp_vcf = clair_dir / "snpsift_output.vcf"
                    indel_vcf = clair_dir / "snpsift_indel_output.vcf"
                    snp_csv = clair_dir / "snpsift_output.vcf.csv"
                    indel_csv = clair_dir / "snpsift_indel_output.vcf.csv"

                    if snp_vcf.exists() and indel_vcf.exists():
                        snp_results_status.set_text("SNP analysis completed successfully!")

                        # Clear previous results
                        snp_results_container.clear()

                        # Display results summary
                        with snp_results_container:
                            with ui.row().classes("w-full gap-3"):
                                with ui.card().classes("flex-1"):
                                    ui.label("SNPs").classes("text-sm font-medium")
                                    if snp_csv.exists():
                                        try:
                                            snp_df = pd.read_csv(snp_csv)
                                            ui.label(f"Total SNPs: {len(snp_df)}").classes(
                                                "text-xs text-gray-600"
                                            )
                                        except Exception:
                                            ui.label("SNP data available").classes(
                                                "text-xs text-green-600"
                                            )
                                    else:
                                        ui.label("SNP data available").classes(
                                            "text-xs text-green-600"
                                        )

                                with ui.card().classes("flex-1"):
                                    ui.label("INDELs").classes("text-sm font-medium")
                                    if indel_csv.exists():
                                        try:
                                            indel_df = pd.read_csv(indel_csv)
                                            ui.label(
                                                f"Total INDELs: {len(indel_df)}"
                                            ).classes("text-xs text-gray-600")
                                        except Exception:
                                            ui.label("INDEL data available").classes(
                                                "text-xs text-green-600"
                                            )
                                    else:
                                        ui.label("INDEL data available").classes(
                                            "text-xs text-green-600"
                                        )

                        # Update button state
                        snp_analysis_button.set_text("Rerun SNP Analysis")
                        snp_analysis_button.props("color=secondary")

                        # Add detailed results viewer
                        with ui.expansion().classes("w-full").props("icon=table_chart"):
                            ui.label("Detailed Results").classes("text-sm font-medium mb-2")

                            # Tabs for SNPs and INDELs
                            with ui.tabs().classes("w-full"):# as tabs:
                                with ui.tab("SNPs", icon="dna"):
                                    _display_variant_table(snp_csv, "SNP", clair_dir)

                                with ui.tab("INDELs", icon="straighten"):
                                    _display_variant_table(indel_csv, "INDEL", clair_dir)

                    else:
                        snp_results_status.set_text("No SNP analysis results found")

                except Exception as e:
                    snp_results_status.set_text(f"Error checking results: {e}")

            def _on_snp_results_expansion_change(e) -> None:
                if not e.value:
                    return
                sig = _snp_results_file_signature()
                if sig == snp_results_lazy_state["sig"]:
                    return
                _check_snp_results()
                snp_results_lazy_state["sig"] = sig

            snp_results_expansion.on_value_change(_on_snp_results_expansion_change)

            # Helper function to detect gzipped files
            def _is_gzipped(file_path):
                """Check if a file is gzipped by reading first few bytes"""
                try:
                    with open(file_path, "rb") as f:
                        return f.read(2) == b"\x1f\x8b"
                except Exception as e:
                    logging.debug(f"   Coverage: <access denied>: {e}")
                    return False

            # Function to display variant table with filtering and details
            def _display_variant_table(csv_file, variant_type, clair_dir):
                """Display variant data in a comprehensive table with filtering"""
                try:
                    # Try to read the VCF file directly instead of the malformed CSV
                    vcf_file = None
                    if variant_type == "SNP":
                        vcf_file = clair_dir / "snpsift_output.vcf"
                    else:
                        vcf_file = clair_dir / "snpsift_indel_output.vcf"

                    if not vcf_file.exists():
                        ui.label(f"No {variant_type} VCF data available").classes(
                            "text-sm text-gray-500"
                        )
                        return

                    # Read the VCF data directly
                    df = None
                    try:
                        # Check if file is compressed
                        if str(vcf_file).endswith(".gz") or _is_gzipped(vcf_file):
                            import gzip

                            with gzip.open(vcf_file, "rt") as f:
                                lines = [
                                    line.strip() for line in f if not line.startswith("#")
                                ]
                        else:
                            with open(vcf_file, "r") as f:
                                lines = [
                                    line.strip() for line in f if not line.startswith("#")
                                ]

                        if not lines:
                            ui.label(f"No {variant_type} variants found").classes(
                                "text-sm text-gray-500"
                            )
                            return

                        # Parse VCF data properly
                        data = []
                        for line in lines:
                            fields = line.split("\t")
                            if len(fields) >= 8:  # Minimum VCF format
                                data.append(
                                    fields[:8]
                                )  # CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO

                        # Create DataFrame with proper column names
                        df = pd.DataFrame(
                            data,
                            columns=[
                                "CHROM",
                                "POS",
                                "ID",
                                "REF",
                                "ALT",
                                "QUAL",
                                "FILTER",
                                "INFO",
                            ],
                        )

                    except Exception as e:
                        ui.label(f"Error reading {variant_type} VCF: {e}").classes(
                            "text-sm text-red-500"
                        )
                        return

                    if df.empty:
                        ui.label(f"No {variant_type} variants found").classes(
                            "text-sm text-gray-500"
                        )
                        return

                    # Display summary statistics
                    with ui.row().classes("w-full gap-3 mb-4"):
                        with ui.card().classes("flex-1"):
                            ui.label(f"Total {variant_type}s").classes(
                                "text-sm font-medium"
                            )
                            ui.label(f"{len(df)}").classes(
                                "text-2xl font-bold text-blue-600"
                            )

                        with ui.card().classes("flex-1"):
                            ui.label("Passing Filters").classes("text-sm font-medium")
                            # Check if we have filter data in the FILTER column
                            if "FILTER" in df.columns:
                                try:
                                    # Count PASS variants
                                    passing = sum(
                                        1 for val in df["FILTER"] if str(val) == "PASS"
                                    )
                                    ui.label(f"{passing}").classes(
                                        "text-2xl font-bold text-green-600"
                                    )
                                except Exception as e:
                                    logging.debug(f"   Coverage: <access denied>: {e}")
                                    ui.label("N/A").classes(
                                        "text-2xl font-bold text-gray-400"
                                    )
                            else:
                                ui.label("N/A").classes("text-2xl font-bold text-gray-400")

                        with ui.card().classes("flex-1"):
                            ui.label("Avg Quality").classes("text-sm font-medium")
                            # Check if we have quality data in the QUAL column
                            if "QUAL" in df.columns:
                                try:
                                    # Calculate average quality from QUAL column
                                    qual_values = pd.to_numeric(df["QUAL"], errors="coerce")
                                    avg_qual = qual_values.mean()
                                    if pd.notna(avg_qual):
                                        ui.label(f"{avg_qual:.1f}").classes(
                                            "text-2xl font-bold text-purple-600"
                                        )
                                    else:
                                        ui.label("N/A").classes(
                                            "text-2xl font-bold text-gray-400"
                                        )
                                except Exception as e:
                                    logging.debug(f"   Coverage: <access denied>: {e}")
                                    ui.label("N/A").classes(
                                        "text-2xl font-bold text-gray-400"
                                    )
                            else:
                                ui.label("N/A").classes("text-2xl font-bold text-gray-400")

                    # Add filtering controls
                    with ui.row().classes("w-full gap-3 mb-4"):
                        # Quality filter
                        quality_threshold = ui.input(
                            value="0", label="Min Quality"
                        ).classes("w-32")

                        # Filter status
                        filter_status = ui.select(
                            options=[
                                "All",
                                "PASS",
                                "NonSomatic",
                                "LowQual",
                                "LowAltBQ",
                                "LowAltMQ",
                            ],
                            value="All",
                            label="Filter Status",
                        ).classes("w-40")

                        # Apply filters button
                        filter_button = ui.button(
                            "Apply Filters", icon="filter_list"
                        ).classes("w-32")

                    # Create filtered dataframe
                    filtered_df = df.copy()

                    def apply_filters():
                        nonlocal filtered_df
                        filtered_df = df.copy()

                        # Apply quality filter
                        try:
                            quality_val = float(quality_threshold.value)
                            if quality_val > 0 and "QUAL" in filtered_df.columns:
                                # Filter by quality in QUAL column
                                qual_values = pd.to_numeric(
                                    filtered_df["QUAL"], errors="coerce"
                                )
                                quality_mask = qual_values >= quality_val
                                filtered_df = filtered_df[quality_mask]
                        except (ValueError, TypeError):
                            pass  # Invalid quality value, skip filtering

                        # Apply filter status
                        if filter_status.value != "All" and "FILTER" in filtered_df.columns:
                            # Filter by status in FILTER column
                            filter_mask = filtered_df["FILTER"] == filter_status.value
                            filtered_df = filtered_df[filter_mask]

                        # Update table
                        if variant_table:
                            variant_table.clear()
                            new_table = _create_variant_table(filtered_df, variant_type)
                            if new_table:
                                variant_table.replace(new_table)

                        # Update summary
                        ui.notify(
                            f"Showing {len(filtered_df)} of {len(df)} {variant_type}s"
                        )

                    filter_button.on_click(apply_filters)

                    # Create conventional table using ui.table.from_pandas
                    def _create_variant_table(data_df, v_type):
                        """Create a conventional table from pandas DataFrame"""
                        if data_df.empty:
                            ui.label(f"No {v_type}s match the current filters").classes(
                                "text-sm text-gray-500"
                            )
                            return None

                        # Prepare DataFrame for display
                        display_df = data_df.copy()

                        # The DataFrame already has proper VCF column names
                        # Select relevant columns for display
                        display_columns = [
                            "CHROM",
                            "POS",
                            "REF",
                            "ALT",
                            "QUAL",
                            "FILTER",
                            "INFO",
                        ]
                        if all(col in display_df.columns for col in display_columns):
                            display_df = display_df[display_columns]
                        else:
                            ui.label(f"Invalid VCF format for {v_type}s").classes(
                                "text-sm text-red-500"
                            )
                            return None

                        # Clean up data
                        display_df = display_df.replace({np.nan: None})

                        # Add quality color coding
                        def color_quality(qual):
                            try:
                                qual_val = float(qual)
                                if qual_val >= 30:
                                    return f"<span style='color: #059669; font-weight: bold;'>{qual_val:.1f}</span>"
                                elif qual_val >= 20:
                                    return f"<span style='color: #d97706; font-weight: bold;'>{qual_val:.1f}</span>"
                                else:
                                    return f"<span style='color: #dc2626; font-weight: bold;'>{qual_val:.1f}</span>"
                            except Exception as e:
                                logging.debug(f"   Coverage: <access denied>: {e}")
                                return str(qual)

                        # Apply quality formatting
                        display_df["QUAL"] = display_df["QUAL"].apply(color_quality)

                        # Add filter color coding
                        def color_filter(filter_val):
                            if filter_val == "PASS":
                                return f"<span style='color: #059669; font-weight: bold;'>{filter_val}</span>"
                            else:
                                return f"<span style='color: #dc2626; font-weight: bold;'>{filter_val}</span>"

                        display_df["FILTER"] = display_df["FILTER"].apply(color_filter)

                        # Truncate INFO field for display
                        def truncate_info(info_val):
                            if info_val and len(str(info_val)) > 50:
                                return str(info_val)[:50] + "..."
                            return str(info_val)

                        display_df["INFO"] = display_df["INFO"].apply(truncate_info)

                        # Create columns definition from DataFrame
                        columns = []
                        for col in display_df.columns:
                            columns.append({
                                "name": col,
                                "label": col,
                                "field": col,
                                "sortable": True
                            })
                        
                        # Create paged table (slice rows per page, avoid full row materialization).
                        variant_table = _render_paged_df_table(
                            display_df,
                            columns,
                            pagination_default=100,
                            class_size="table-xs",
                            search_placeholder=f"Search {v_type}s...",
                        )

                        # Make all columns sortable
                        for col in variant_table.columns:
                            col["sortable"] = True

                        # Add table controls
                        with variant_table.add_slot("top-left"):

                            def toggle_fs():
                                variant_table.toggle_fullscreen()
                                fs_button.props(
                                    "icon=fullscreen_exit"
                                    if variant_table.is_fullscreen
                                    else "icon=fullscreen"
                                )

                            fs_button = ui.button(
                                "Toggle fullscreen",
                                icon="fullscreen",
                                on_click=toggle_fs,
                            ).props("flat")

                            # Column visibility toggle
                            with ui.button(icon="menu").props("flat"):
                                with ui.menu(), ui.column().classes("gap-0 p-2"):
                                    for column in variant_table.columns:
                                        ui.switch(
                                            column["label"],
                                            value=True,
                                            on_change=lambda e, column=column: toggle_column(
                                                column, e.value
                                            ),
                                        )

                        # Column toggle function
                        def toggle_column(column: dict, visible: bool) -> None:
                            column["classes"] = "" if visible else "hidden"
                            column["headerClasses"] = "" if visible else "hidden"
                            variant_table.update()

                        return variant_table

                    # Create and display the table
                    variant_table = _create_variant_table(filtered_df, variant_type)

                    # Add export functionality
                    with ui.row().classes("w-full mt-4"):
                        ui.button(
                            f"Export {variant_type}s to CSV",
                            icon="download",
                            on_click=lambda: _export_variants(filtered_df, variant_type),
                        ).classes("w-full")

                        # Add row count display
                        ui.label(
                            f"Showing {len(filtered_df)} of {len(df)} {variant_type}s"
                        ).classes("text-sm text-gray-600 ml-auto")

                except Exception as e:
                    ui.label(f"Error displaying {variant_type} table: {e}").classes(
                        "text-sm text-red-500"
                    )

            # Function to show variant details
            def _show_variant_details(variant_row, variant_type, clair_dir):
                """Show detailed information for a specific variant"""
                try:
                    with ui.dialog() as dialog, ui.card().classes(
                        "robin-dialog-surface p-4 md:p-5 w-full max-w-2xl max-h-[85vh] overflow-auto"
                    ):
                        ui.label(f"{variant_type} details").classes(
                            "classification-insight-heading text-headline-small mb-2"
                        )

                        with ui.column().classes("gap-2"):
                            for col, val in variant_row.items():
                                if pd.notna(val) and str(val).strip():
                                    with ui.row().classes("w-full items-start gap-2"):
                                        ui.label(f"{col}:").classes(
                                            "text-xs font-medium min-w-[8rem] shrink-0 "
                                            "target-coverage-panel__meta-label"
                                        )
                                        ui.label(str(val)).classes(
                                            "text-xs flex-1 break-all whitespace-pre-wrap"
                                        )

                        ui.button("Close", on_click=dialog.close).props(
                            "color=primary no-caps outline"
                        ).classes("w-full mt-4")

                except Exception as e:
                    ui.notify(f"Error showing variant details: {e}", type="error")

            # Function to export variants
            def _export_variants(data_df, variant_type):
                """Export filtered variants to CSV"""
                try:
                    import tempfile

                    # Create temporary file
                    with tempfile.NamedTemporaryFile(
                        mode="w", suffix=".csv", delete=False
                    ) as f:
                        data_df.to_csv(f.name, index=False)
                        temp_path = f.name

                    # Download the file
                    ui.download(temp_path, filename=f"{variant_type.lower()}_variants.csv")

                    # Clean up
                    os.unlink(temp_path)

                    ui.notify(f"{variant_type} variants exported successfully!")

                except Exception as e:
                    ui.notify(f"Export failed: {e}", type="error")

            # Function to check file availability
            def _check_snp_file_availability():
                try:
                    target_bam = sample_dir / "target.bam"
                    targets_bed = sample_dir / "targets_exceeding_threshold.bed"

                    if not target_bam.exists():
                        snp_files_status.set_text(
                            " Missing target.bam - run target analysis first"
                        )
                        snp_files_status.classes(replace="text-xs text-red-500")
                        snp_analysis_button.disable()
                        return False

                    if not targets_bed.exists():
                        snp_files_status.set_text(
                            " Missing targets_exceeding_threshold.bed - run target analysis first"
                        )
                        snp_files_status.classes(replace="text-xs text-red-500")
                        snp_analysis_button.disable()
                        return False

                    # Check if BED file has content
                    try:
                        with open(targets_bed, "r") as f:
                            bed_content = f.read().strip()

                        if not bed_content:
                            snp_files_status.set_text(
                                "BED file is empty - no regions exceed threshold"
                            )
                            snp_files_status.classes(replace="text-xs text-yellow-500")
                            snp_analysis_button.disable()
                            return False
                        else:
                            snp_files_status.set_text("All required files available")
                            snp_files_status.classes(replace="text-xs text-green-500")
                            snp_analysis_button.enable()
                            return True

                    except Exception as e:
                        snp_files_status.set_text(f" Error reading BED file: {e}")
                        snp_files_status.classes(replace="text-xs text-red-500")
                        snp_analysis_button.disable()
                        return False

                except Exception as e:
                    snp_files_status.set_text(f" Error checking files: {e}")
                    snp_files_status.classes(replace="text-xs text-red-500")
                    snp_analysis_button.disable()
                    return False

            # Check file availability initially
            _check_snp_file_availability()

            # Function to trigger SNP analysis
            def _trigger_snp_analysis():
                try:
                    # Check prerequisites
                    target_bam = sample_dir / "target.bam"
                    targets_bed = sample_dir / "targets_exceeding_threshold.bed"

                    if not target_bam.exists():
                        ui.notify(
                            "Missing target.bam file. Please run target analysis first.",
                            type="warning",
                        )
                        return

                    if not targets_bed.exists():
                        ui.notify(
                            "Missing targets_exceeding_threshold.bed file. Please run target analysis first.",
                            type="warning",
                        )
                        return

                    # Check if BED file has content (not empty)
                    try:
                        with open(targets_bed, "r") as f:
                            bed_content = f.read().strip()

                        if not bed_content:
                            ui.notify(
                                "targets_exceeding_threshold.bed file is empty. No regions exceed coverage threshold.",
                                type="warning",
                            )
                            return

                    except Exception as e:
                        ui.notify(
                            f"Error reading targets_exceeding_threshold.bed file: {e}",
                            type="error",
                        )
                        return

                    # Check for reference genome
                    # Remove the unnecessary fallback reference check since we get it from CLI
                    # The reference genome is passed via CLI and stored in workflow_runner.reference

                    # Update UI state
                    snp_results_lazy_state["sig"] = None
                    snp_analysis_button.disable()
                    snp_status_label.set_text("Starting SNP analysis...")
                    snp_status_label.classes(replace="text-sm text-blue-600")

                    # Create and submit SNP analysis job
                    try:
                        # Import required modules
                        from robin.analysis.target_analysis import snp_analysis_handler
                        from robin.workflow_simple import Job, WorkflowContext

                        # Create job metadata with reference genome
                        # Get reference genome from workflow runner (set by CLI)
                        reference_genome = None

                        # Try to get reference from workflow runner
                        print("=== REFERENCE GENOME DEBUGGING ===")
                        print(f"launcher type: {type(launcher)}")
                        print(
                            f"launcher dir: {[attr for attr in dir(launcher) if not attr.startswith('_')]}"
                        )
                        print(
                            f"hasattr(launcher, 'workflow_runner'): {hasattr(launcher, 'workflow_runner')}"
                        )

                        if hasattr(launcher, "workflow_runner"):
                            print(f"launcher.workflow_runner: {launcher.workflow_runner}")
                            if launcher.workflow_runner:
                                print(
                                    f"workflow_runner type: {type(launcher.workflow_runner)}"
                                )
                                print(
                                    f"workflow_runner dir: {[attr for attr in dir(launcher.workflow_runner) if not attr.startswith('_')]}"
                                )
                                print(
                                    f"hasattr(workflow_runner, 'reference'): {hasattr(launcher.workflow_runner, 'reference')}"
                                )

                                if hasattr(launcher.workflow_runner, "reference"):
                                    reference_genome = launcher.workflow_runner.reference
                                    print(f"workflow_runner.reference: {reference_genome}")
                                    if reference_genome:
                                        print(
                                            f"SUCCESS: Using reference genome from workflow runner: {reference_genome}"
                                        )
                                    else:
                                        print(
                                            " FAILED: No reference genome in workflow runner"
                                        )
                                else:
                                    print(" Workflow runner has no reference attribute")
                            else:
                                print(" Workflow runner is None")
                        else:
                            print(" No workflow runner attribute found")

                        print("=== END REFERENCE GENOME DEBUGGING ===")

                        # Fallback to environment variable
                        if not reference_genome:
                            env_reference = os.environ.get("robin_REFERENCE")
                            if env_reference and os.path.exists(env_reference):
                                reference_genome = env_reference
                                print(
                                    f"Using reference genome from environment: {reference_genome}"
                                )

                        if not reference_genome:
                            ui.notify(
                                "No reference genome found. SNP calling may fail. Please ensure --reference is provided via CLI.",
                                type="warning",
                            )
                        else:
                            ui.notify(
                                f"Reference genome found: {os.path.basename(reference_genome)}",
                                type="positive",
                            )

                        # Use monitored_directory as work_dir, or fall back to sample_dir if not available
                        work_dir = (
                            launcher.monitored_directory
                            if launcher.monitored_directory
                            else str(sample_dir)
                        )

                        metadata = {
                            "work_dir": work_dir,
                            "threads": 4,
                            "force_regenerate": force_regenerate_checkbox.value,
                            "reference": reference_genome,
                        }

                        print("=== JOB METADATA DEBUGGING ===")
                        print(f"Created metadata: {metadata}")
                        print(f"reference_genome value: {reference_genome}")
                        print("=== END JOB METADATA DEBUGGING ===")

                        # Create workflow context
                        context = WorkflowContext(
                            filepath=str(sample_dir), metadata=metadata
                        )

                        # Add sample_id method
                        def get_sample_id():
                            return sample_dir.name

                        context.get_sample_id = get_sample_id
                        context.add_result = lambda key, value: None  # Mock for now
                        context.add_error = lambda key, value: None  # Mock for now

                        # Create job
                        job = Job(
                            job_id=hash(f"snp_analysis_{sample_dir.name}")
                            % 1000000,  # Simple hash-based ID
                            job_type="snp_analysis",
                            context=context,
                            origin="fast",
                            workflow=["fast:snp_analysis"],
                        )

                        # Submit SNP analysis as a workflow job (non-blocking)
                        try:
                            # Get the workflow runner from the launcher
                            if (
                                hasattr(launcher, "workflow_runner")
                                and launcher.workflow_runner is not None
                            ):
                                workflow_runner = launcher.workflow_runner

                                # Submit the job through the workflow system
                                if hasattr(workflow_runner, "submit_snp_analysis_job"):
                                    # Use the dedicated SNP analysis method
                                    success = workflow_runner.submit_snp_analysis_job(
                                        sample_dir=work_dir,
                                        sample_id=sample_id,
                                        reference=(
                                            str(reference_genome)
                                            if reference_genome
                                            else None
                                        ),
                                        threads=4,
                                        force_regenerate=False,
                                    )

                                    if success:
                                        ui.notify(
                                            "SNP analysis job submitted to workflow queue",
                                            type="info",
                                        )
                                        snp_status_label.set_text(
                                            "Job submitted to workflow queue"
                                        )
                                        snp_status_label.classes(
                                            replace="text-sm text-blue-600"
                                        )
                                    else:
                                        raise Exception(
                                            "Failed to submit SNP analysis job to workflow"
                                        )

                                elif hasattr(workflow_runner, "submit_sample_job"):
                                    # Fallback to generic job submission
                                    success = workflow_runner.submit_sample_job(
                                        sample_dir=work_dir,
                                        job_type="snp_analysis",
                                        sample_id=sample_id,
                                    )

                                    if success:
                                        ui.notify(
                                            "SNP analysis job submitted to workflow queue",
                                            type="info",
                                        )
                                        snp_status_label.set_text(
                                            "Job submitted to workflow queue"
                                        )
                                        snp_status_label.classes(
                                            replace="text-sm text-blue-600"
                                        )
                                    else:
                                        raise Exception(
                                            "Failed to submit SNP analysis job to workflow"
                                        )
                                else:
                                    raise Exception(
                                        "Workflow runner does not support SNP analysis job submission"
                                    )
                            else:
                                raise Exception("No workflow runner available")

                        except Exception as e:
                            # Fallback to direct handler call if workflow submission fails
                            ui.notify(
                                f"Workflow submission failed, falling back to direct execution: {e}",
                                type="warning",
                            )

                            def run_snp_analysis_fallback():
                                try:
                                    # Use the work_dir variable we calculated earlier
                                    snp_analysis_handler(job, work_dir=work_dir)

                                    # Update UI on completion using the main thread
                                    print("SNP analysis completed successfully!")

                                except Exception as e:
                                    # Log error instead of trying to update UI from background thread
                                    ui.run_javascript(
                                        f"""
                                        // Update status
                                        var statusElement = document.querySelector('{snp_status_label.id}');
                                        if (statusElement) {{
                                            statusElement.textContent = "SNP analysis failed: " + "{str(e)}";
                                            statusElement.className = "text-sm text-red-600";
                                        }}
                                        
                                        // Re-enable button
                                        var buttonElement = document.querySelector('{snp_analysis_button.id}');
                                        if (buttonElement) {{
                                            buttonElement.disabled = false;
                                        }}
                                    """
                                    )

                            # Run in background thread as fallback
                            import threading

                            thread = threading.Thread(
                                target=run_snp_analysis_fallback, daemon=True
                            )
                            thread.start()

                        ui.notify("SNP analysis started in background", type="info")

                    except Exception as e:
                        ui.notify(f"Failed to start SNP analysis: {e}", type="error")
                        snp_analysis_button.enable()
                        snp_status_label.set_text("Failed to start SNP analysis")
                        snp_status_label.classes(replace="text-sm text-red-600")

                except Exception as e:
                    ui.notify(f"Error triggering SNP analysis: {e}", type="error")
                    snp_analysis_button.enable()
                    snp_status_label.set_text("Error triggering SNP analysis")
                    snp_status_label.classes(replace="text-sm text-red-600")

    # Lightweight Gene Analysis section (only shown in development mode)
    if is_development_mode:
        with ui.card().classes("w-full"):
            ui.label("Lightweight Variant Analysis").classes(
                "text-lg font-semibold mb-2"
            )

            # Lightweight Variant Analysis controls and status
            with ui.row().classes("w-full items-center justify-between mb-4"):
                with ui.column().classes("gap-2"):
                    ui.label("Variant Analysis").classes("text-sm font-medium")
                    lga_status_label = ui.label(
                        "Ready to run lightweight variant analysis"
                    ).classes("text-sm text-gray-600")

                    # Add file availability status
                    lga_files_status = ui.label("Checking required files...").classes(
                        "text-xs text-gray-500"
                    )

                with ui.column().classes("gap-2 items-end"):
                    # Lightweight Variant Analysis button
                    lga_analysis_button = ui.button(
                        "Run Variant Analysis",
                        icon="dna",
                        on_click=lambda: _trigger_lightweight_gene_analysis(),
                    ).props("color=primary data-lga-analysis-button")

                    # Force regenerate checkbox
                    lga_force_regenerate_checkbox = ui.checkbox(
                        "Force regenerate", value=False
                    ).props("dense")

            # Lightweight Variant Analysis requirements check
            with ui.expansion().classes("w-full").props("icon=info"):
                ui.label("Requirements").classes("text-sm font-medium mb-2")
                with ui.column().classes("gap-1 text-sm"):
                    ui.label(
                        "• target.bam - Target regions BAM file (generated by target analysis)"
                    )
                    ui.label(
                        "• targets_exceeding_threshold.bed - BED file defining regions (auto-generated by target analysis)"
                    )
                    ui.label(
                        "• ClinVar VCF file - Pathogenic variants database (optional, will use default)"
                    )
                    ui.label("• Reference genome (optional) - for enhanced analysis")

                # Add helpful note about file generation
                with ui.row().classes("w-full mt-2 p-2 bg-blue-50 rounded"):
                    ui.icon("info").classes("text-blue-600")
                    ui.label(
                        "Note: This analysis automatically finds ClinVar pathogenic variants that intersect with your BED regions and analyzes coverage at those sites. It's much faster than full variant calling."
                    ).classes("text-sm text-blue-800")

            # Lightweight Variant Analysis results display
            with ui.expansion().classes("w-full").props("icon=assessment"):
                ui.label("Results").classes("text-sm font-medium mb-2")

                # Results status
                lga_results_status = (
                    ui.label("No lightweight variant analysis results yet")
                    .classes("text-sm text-gray-600")
                    .props("data-lga-results-status")
                )

                # Results container
                lga_results_container = ui.column().classes("w-full")

                # Function to check and display lightweight gene analysis results
                def _check_lga_results():
                    try:
                        # Check for JSON results file in the correct location
                        # The JSON file is written to {sample_dir}/{sample_id}/lightweight_gene_analysis_results.json
                        sample_id = sample_dir.name
                        json_report = (
                            sample_dir
                            / sample_id
                            / "lightweight_gene_analysis_results.json"
                        )

                        # Also check the old location as fallback
                        json_report_fallback = (
                            sample_dir / "lightweight_gene_analysis_results.json"
                        )

                        if json_report.exists():
                            lga_results_status.set_text(
                                "Lightweight variant analysis completed successfully!"
                            )
                            _display_lga_results(json_report)
                        elif json_report_fallback.exists():
                            lga_results_status.set_text(
                                "Lightweight variant analysis completed successfully! (legacy location)"
                            )
                            _display_lga_results(json_report_fallback)
                        else:
                            lga_results_status.set_text(
                                "No lightweight variant analysis results found"
                            )

                    except Exception as e:
                        lga_results_status.set_text(f"Error checking results: {e}")
                        import traceback

                        traceback.print_exc()

                def _display_lga_results(json_report):
                    """Display the results interface when JSON file is found"""
                    # Clear previous results
                    lga_results_container.clear()

                    # Display results summary
                    with lga_results_container:
                        with ui.row().classes("w-full gap-3"):
                            with ui.card().classes("flex-1"):
                                ui.label("JSON Results").classes("text-sm font-medium")
                                ui.label("Available").classes("text-xs text-green-600")

                            with ui.card().classes("flex-1"):
                                ui.label("Analysis Complete").classes("text-sm font-medium")
                                ui.label("Ready to view").classes("text-xs text-green-600")

                        # Add view results button
                        with ui.row().classes("w-full mt-4"):
                            view_button = ui.button(
                                "View Analysis Results",
                                icon="analytics",
                                on_click=lambda: _view_lga_json_results_inline(
                                    str(json_report)
                                ),
                            ).classes("flex-1")

                            # Add debug info
                            print(f" GUI: Created view button with path: {json_report}")
                            print(f" GUI: Button object: {view_button}")

                    # Update button state
                    lga_analysis_button.set_text("Rerun Variant Analysis")
                    lga_analysis_button.props("color=secondary")

                # Function to view JSON analysis results inline
                def _view_lga_json_results_inline(json_path):
                    try:
                        import json

                        print(f" GUI: Function called with path: {json_path}")
                        print(f" GUI: Path type: {type(json_path)}")

                        # Convert to Path object if it's a string
                        if isinstance(json_path, str):
                            json_path = Path(json_path)

                        print(f" GUI: Attempting to load JSON from: {json_path}")
                        print(f" GUI: File exists: {json_path.exists()}")
                        print(
                            f" GUI: File size: {json_path.stat().st_size if json_path.exists() else 'N/A'} bytes"
                        )

                        # Read and parse JSON
                        with open(json_path, "r") as f:
                            data = json.load(f)

                        print(
                            f" GUI: Successfully loaded JSON with {len(data.get('genes', {}))} genes"
                        )

                        # Clear previous results and display inline
                        lga_results_container.clear()

                        with lga_results_container:
                            # Header
                            ui.label("Lightweight Variant Analysis Results").classes(
                                "text-2xl font-bold mb-6"
                            )

                            # Metadata summary
                            with ui.card().classes("w-full mb-6"):
                                ui.label("Analysis Summary").classes(
                                    "text-xl font-semibold mb-4"
                                )
                                metadata = data.get("metadata", {})

                                # Linear layout - each metric on its own row
                                with ui.column().classes("w-full gap-3"):
                                    with ui.row().classes(
                                        "w-full items-center justify-between p-4 bg-blue-50 rounded"
                                    ):
                                        ui.label("Total Genes Analyzed").classes(
                                            "text-lg font-medium text-gray-700"
                                        )
                                        ui.label(
                                            f"{metadata.get('total_genes_analyzed', 0)}"
                                        ).classes("text-3xl font-bold text-blue-600")

                                    with ui.row().classes(
                                        "w-full items-center justify-between p-4 bg-purple-50 rounded"
                                    ):
                                        ui.label("Total Pathogenic Variants Found").classes(
                                            "text-lg font-medium text-gray-700"
                                        )
                                        ui.label(
                                            f"{metadata.get('total_variants_found', 0):,}"
                                        ).classes("text-3xl font-bold text-purple-600")

                                    with ui.row().classes(
                                        "w-full items-center justify-between p-4 bg-green-50 rounded"
                                    ):
                                        ui.label("Genes with Good Coverage (≥10x)").classes(
                                            "text-lg font-medium text-gray-700"
                                        )
                                        ui.label(
                                            f"{metadata.get('genes_with_good_coverage_variants', 0)}"
                                        ).classes("text-3xl font-bold text-green-600")

                                    with ui.row().classes(
                                        "w-full items-center justify-between p-4 bg-orange-50 rounded"
                                    ):
                                        ui.label("Genes with Low Coverage (<10x)").classes(
                                            "text-lg font-medium text-gray-700"
                                        )
                                        low_coverage = metadata.get(
                                            "total_genes_analyzed", 0
                                        ) - metadata.get(
                                            "genes_with_good_coverage_variants", 0
                                        )
                                        ui.label(f"{low_coverage}").classes(
                                            "text-3xl font-bold text-orange-600"
                                        )

                            # Gene analysis sections - linear layout, no tabs
                            ui.label("Gene Analysis Results").classes(
                                "text-2xl font-bold mt-8 mb-6"
                            )

                            # Display each section sequentially
                            _display_gene_overview(data)
                            _display_high_coverage_genes(data)
                            _display_low_coverage_genes(data)
                            _display_gene_details(data)

                            # Add hide results button
                            with ui.row().classes("w-full mt-6 justify-center"):
                                ui.button(
                                    "Hide Results",
                                    icon="visibility_off",
                                    on_click=lambda: _hide_lga_results(),
                                ).props("color=secondary")

                        print(" GUI: Results displayed inline successfully")

                    except Exception as e:
                        ui.notify(f"Error viewing results: {e}", type="error")
                        import traceback

                        print(f"Error in _view_lga_json_results_inline: {e}")
                        traceback.print_exc()

                # Function to hide results and restore original view
                def _hide_lga_results():
                    """Hide the inline results and restore the original simple view"""
                    try:
                        print(" GUI: Hiding results...")

                        # Clear the results container
                        lga_results_container.clear()

                        # Restore the original simple view
                        with lga_results_container:
                            with ui.row().classes("w-full gap-3"):
                                with ui.card().classes("flex-1"):
                                    ui.label("JSON Results").classes("text-sm font-medium")
                                    ui.label("Available").classes("text-xs text-green-600")

                                with ui.card().classes("flex-1"):
                                    ui.label("Analysis Complete").classes(
                                        "text-sm font-medium"
                                    )
                                    ui.label("Ready to view").classes(
                                        "text-xs text-green-600"
                                    )

                            # Add view results button
                            with ui.row().classes("w-full mt-4"):
                                ui.button(
                                    "View Analysis Results",
                                    icon="analytics",
                                    on_click=lambda: _view_lga_json_results_inline(
                                        str(
                                            Path(sample_dir)
                                            / sample_dir.name
                                            / "lightweight_gene_analysis_results.json"
                                        )
                                    ),
                                ).classes("flex-1")

                        print(" GUI: Results hidden successfully")

                    except Exception as e:
                        print(f" GUI: Error hiding results: {e}")
                        import traceback

                        traceback.print_exc()

                # Helper functions for displaying JSON results
                def _display_gene_overview(data):
                    """Display gene overview with statistics and charts"""
                    try:
                        metadata = data.get("metadata", {})
                        genes = data.get("genes", {})

                        # Coverage distribution chart
                        coverage_data = []
                        for gene_name, gene_data in genes.items():
                            coverage_stats = gene_data.get("coverage_statistics", {})
                            mean_cov = coverage_stats.get("mean_coverage", 0)
                            coverage_data.append(
                                {
                                    "gene": gene_name,
                                    "mean_coverage": mean_cov,
                                    "total_variants": gene_data.get("summary", {}).get(
                                        "total_variants", 0
                                    ),
                                    "high_coverage_variants": gene_data.get(
                                        "summary", {}
                                    ).get("high_coverage_variants", 0),
                                }
                            )

                        # Sort by mean coverage
                        coverage_data.sort(key=lambda x: x["mean_coverage"], reverse=True)

                        # Create coverage distribution chart
                        with ui.card().classes("w-full mb-4"):
                            ui.label("Coverage Distribution by Gene").classes(
                                "text-lg font-semibold mb-2"
                            )

                            # Top 20 genes by coverage
                            top_genes = coverage_data[:20]
                            gene_names = [g["gene"] for g in top_genes]
                            coverage_values = [g["mean_coverage"] for g in top_genes]

                            # Create enhanced bar chart using echart
                            _lga_p = _cov_chrome_palette()
                            coverage_chart = ui.echart(
                                {
                                    "backgroundColor": "transparent",
                                    "textStyle": {"color": _lga_p["axis"]},
                                    "title": {
                                        "text": "Top 20 Genes by Mean Coverage",
                                        "left": "center",
                                        "top": 10,
                                        "textStyle": {
                                            "fontSize": 16,
                                            "color": _lga_p["title"],
                                            "fontWeight": "bold",
                                        },
                                        "subtext": "Higher coverage = More reliable variant detection",
                                        "subtextStyle": {
                                            "fontSize": 12,
                                            "color": _lga_p["axis"],
                                        },
                                    },
                                    "tooltip": {
                                        "trigger": "axis",
                                        "formatter": "Gene: {b}<br/>Coverage: {c}x",
                                        **_cov_tooltip_option(),
                                    },
                                    "grid": {
                                        "left": "12%",
                                        "right": "8%",
                                        "bottom": "20%",
                                        "top": "30%",
                                        "containLabel": True,
                                    },
                                    "xAxis": {
                                        "type": "category",
                                        "data": gene_names,
                                        "axisLabel": {
                                            "rotate": 45,
                                            "fontSize": 11,
                                            "fontWeight": "bold",
                                            "color": _lga_p["axis"],
                                        },
                                        "name": "Gene Names",
                                        "nameLocation": "middle",
                                        "nameGap": 50,
                                        "nameTextStyle": {"color": _lga_p["axis"]},
                                        "axisLine": {"lineStyle": {"color": _lga_p["axis"]}},
                                    },
                                    "yAxis": {
                                        "type": "value",
                                        "name": "Mean Coverage (x)",
                                        "nameLocation": "middle",
                                        "nameGap": 40,
                                        "nameTextStyle": {"color": _lga_p["axis"]},
                                        "axisLabel": {
                                            "fontSize": 12,
                                            "fontWeight": "bold",
                                            "color": _lga_p["axis"],
                                        },
                                        "splitLine": {
                                            "show": True,
                                            "lineStyle": {
                                                "color": _lga_p["split"],
                                                "type": "dashed",
                                            },
                                        },
                                    },
                                    "series": [
                                        {
                                            "type": "bar",
                                            "data": coverage_values,
                                            "itemStyle": {
                                                "color": (
                                                    "#60a5fa"
                                                    if _cov_ui_dark()
                                                    else "#3b82f6"
                                                )
                                            },
                                            "barWidth": "60%",
                                            "emphasis": {
                                                "itemStyle": {
                                                    "shadowBlur": 10,
                                                    "shadowColor": (
                                                        "rgba(0,0,0,0.45)"
                                                        if _cov_ui_dark()
                                                        else "rgba(0,0,0,0.3)"
                                                    ),
                                                }
                                            },
                                        }
                                    ],
                                    "dataZoom": [
                                        {
                                            "type": "slider",
                                            "show": True,
                                            "xAxisIndex": [0],
                                            "start": 0,
                                            "end": 100,
                                            "textStyle": {"color": _lga_p["axis"]},
                                            "borderColor": _lga_p["tooltip_border"],
                                            "handleStyle": {"color": _lga_p["axis"]},
                                        }
                                    ],
                                }
                            ).classes("w-full h-80")
                            try:
                                launcher._coverage_state.setdefault(
                                    str(sample_dir), {}
                                )["lga_gene_coverage_chart"] = coverage_chart
                            except Exception:
                                pass
                            _apply_lga_gene_bar_chrome(coverage_chart)

                            # Add coverage threshold indicators
                            with ui.column().classes("w-full mt-3 gap-2"):
                                ui.label("Coverage Thresholds:").classes(
                                    "text-sm font-medium text-gray-600"
                                )

                                with ui.row().classes(
                                    "items-center gap-3 p-2 bg-green-50 rounded"
                                ):
                                    ui.element("div").classes("w-4 h-4 bg-green-500 rounded")
                                    ui.label(
                                        "≥10x (Good Coverage) - Reliable variant detection"
                                    ).classes("text-sm text-gray-700")

                                with ui.row().classes(
                                    "items-center gap-3 p-2 bg-orange-50 rounded"
                                ):
                                    ui.element("div").classes("w-4 h-4 bg-orange-500 rounded")
                                    ui.label(
                                        "5-10x (Moderate Coverage) - Limited reliability"
                                    ).classes("text-sm text-gray-700")

                                with ui.row().classes(
                                    "items-center gap-3 p-2 bg-red-50 rounded"
                                ):
                                    ui.element("div").classes("w-4 h-4 bg-red-500 rounded")
                                    ui.label(
                                        "<5x (Low Coverage) - Poor reliability"
                                    ).classes(
                                        "text-sm text-gray-700"
                                    )

                        # Summary statistics
                        with ui.card().classes("w-full mb-4"):
                            ui.label("Summary Statistics").classes(
                                "text-lg font-semibold mb-2"
                            )

                            # Calculate additional statistics
                            total_coverage = sum(g["mean_coverage"] for g in coverage_data)
                            avg_coverage = (
                                total_coverage / len(coverage_data) if coverage_data else 0
                            )
                            high_coverage_count = sum(
                                1 for g in coverage_data if g["mean_coverage"] >= 10
                            )
                            low_coverage_count = sum(
                                1 for g in coverage_data if g["mean_coverage"] < 10
                            )

                            # Linear layout for statistics
                            with ui.column().classes("w-full gap-3"):
                                with ui.row().classes(
                                    "w-full items-center justify-between p-3 bg-blue-50 rounded"
                                ):
                                    ui.label("Average Coverage Across All Genes").classes(
                                        "text-base font-medium text-gray-700"
                                    )
                                    ui.label(f"{avg_coverage:.1f}x").classes(
                                        "text-2xl font-bold text-blue-600"
                                    )

                                with ui.row().classes(
                                    "w-full items-center justify-between p-3 bg-green-50 rounded"
                                ):
                                    ui.label("Genes with High Coverage (≥10x)").classes(
                                        "text-base font-medium text-gray-700"
                                    )
                                    ui.label(f"{high_coverage_count}").classes(
                                        "text-2xl font-bold text-green-600"
                                    )

                                with ui.row().classes(
                                    "w-full items-center justify-between p-3 bg-orange-50 rounded"
                                ):
                                    ui.label("Genes with Low Coverage (<10x)").classes(
                                        "text-base font-medium text-gray-700"
                                    )
                                    ui.label(f"{low_coverage_count}").classes(
                                        "text-2xl font-bold text-orange-600"
                                    )

                                with ui.row().classes(
                                    "w-full items-center justify-between p-3 bg-purple-50 rounded"
                                ):
                                    ui.label("Total Pathogenic Variants").classes(
                                        "text-base font-medium text-gray-700"
                                    )
                                    ui.label(
                                        f"{metadata.get('total_variants_found', 0):,}"
                                    ).classes("text-2xl font-bold text-purple-600")

                        # Gene list with quick stats
                        with ui.card().classes("w-full"):
                            ui.label("Gene List with Coverage Status").classes(
                                "text-lg font-semibold mb-2"
                            )

                            # Create searchable table
                            gene_table_data = []
                            for gene_name, gene_data in genes.items():
                                summary = gene_data.get("summary", {})
                                coverage_stats = gene_data.get("coverage_statistics", {})

                                gene_table_data.append(
                                    {
                                        "gene": gene_name,
                                        "mean_coverage": f"{coverage_stats.get('mean_coverage', 0):.1f}x",
                                        "total_variants": summary.get("total_variants", 0),
                                        "high_coverage_variants": summary.get(
                                            "high_coverage_variants", 0
                                        ),
                                        "low_coverage_variants": summary.get(
                                            "low_coverage_variants", 0
                                        ),
                                        "status": summary.get(
                                            "pathogenic_status", "unknown"
                                        )
                                        .replace("_", " ")
                                        .title(),
                                    }
                                )

                            # Sort by mean coverage
                            gene_table_data.sort(
                                key=lambda x: float(x["mean_coverage"].replace("x", "")),
                                reverse=True,
                            )

                            # Create table using styled_table for consistency
                            from robin.gui.theme import styled_table
                            
                            # Create columns definition
                            columns = []
                            for col in gene_table_data[0].keys():
                                columns.append({
                                    "name": col,
                                    "label": col.replace("_", " ").title(),
                                    "field": col,
                                    "sortable": True
                                })
                            
                            # Create styled table
                            table_container, gene_table = styled_table(
                                columns=columns,
                                rows=gene_table_data,
                                pagination=25,
                                class_size="table-xs"
                            )

                            # Add search functionality
                            with gene_table.add_slot("top-right"):
                                with ui.input(placeholder="Search genes...").props(
                                    "type=search"
                                ).bind_value(gene_table, "filter").add_slot("append"):
                                    ui.icon("search")

                            # Make columns sortable
                            for col in gene_table.columns:
                                col["sortable"] = True

                    except Exception as e:
                        ui.label(f"Error displaying gene overview: {e}").classes(
                            "text-sm text-red-500"
                        )
                        import traceback

                        traceback.print_exc()

                def _display_high_coverage_genes(data):
                    """Display genes with high coverage variants"""
                    try:
                        genes = data.get("genes", {})

                        # Filter for genes with high coverage variants
                        high_coverage_genes = {}
                        for gene_name, gene_data in genes.items():
                            summary = gene_data.get("summary", {})
                            if summary.get("high_coverage_variants", 0) > 0:
                                high_coverage_genes[gene_name] = gene_data

                        if not high_coverage_genes:
                            ui.label("No genes with high coverage variants found.").classes(
                                "text-sm text-gray-500"
                            )
                            return

                        ui.label(
                            f"Found {len(high_coverage_genes)} genes with high coverage variants"
                        ).classes("text-lg font-semibold mb-4")

                        # Display each gene with its high coverage variants
                        for gene_name, gene_data in high_coverage_genes.items():
                            with ui.expansion().classes("w-full mb-2").props(
                                f"icon=dna label={gene_name}"
                            ):
                                summary = gene_data.get("summary", {})
                                coverage_stats = gene_data.get("coverage_statistics", {})
                                variants = gene_data.get("variants", [])

                                # Gene summary
                                with ui.row().classes("w-full gap-3 mb-2"):
                                    with ui.card().classes("flex-1"):
                                        ui.label("Mean Coverage").classes(
                                            "text-sm text-gray-600"
                                        )
                                        ui.label(
                                            f"{coverage_stats.get('mean_coverage', 0):.1f}x"
                                        ).classes("text-lg font-bold")

                                    with ui.card().classes("flex-1"):
                                        ui.label("High Coverage Variants").classes(
                                            "text-sm text-gray-600"
                                        )
                                        ui.label(
                                            f"{summary.get('high_coverage_variants', 0)}"
                                        ).classes("text-lg font-bold text-green-600")

                                    with ui.card().classes("flex-1"):
                                        ui.label("Total Variants").classes(
                                            "text-sm text-gray-600"
                                        )
                                        ui.label(
                                            f"{summary.get('total_variants', 0)}"
                                        ).classes("text-lg font-bold")

                                # High coverage variants table - only show variants where alternate allele is detected
                                high_cov_variants = [
                                    v
                                    for v in variants
                                    if v.get("coverage_at_variant", 0) >= 10
                                    and v.get("evidence_analysis", {}).get(
                                        "alternate_support", 0
                                    )
                                    > 0
                                ]

                                if high_cov_variants:
                                    ui.label("High Coverage Variants (≥10x):").classes(
                                        "text-sm font-medium mb-2"
                                    )

                                    # Create table for high coverage variants
                                    variant_data = []
                                    for variant in high_cov_variants[
                                        :20
                                    ]:  # Limit to first 20
                                        evidence = variant.get("evidence_analysis", {})
                                        # Create genomic locus in format chr1:119,915,631
                                        chromosome = variant.get("chromosome", "N")
                                        position = variant.get("position", 0)
                                        genomic_locus = (
                                            f"{chromosome}:{position:,}"
                                            if chromosome != "N" and position
                                            else "N/A"
                                        )

                                        variant_data.append(
                                            {
                                                "genomic_locus": genomic_locus,
                                                "position": f"{position:,}",
                                                "reference": variant.get("reference", "N"),
                                                "alternate": variant.get("alternate", "N"),
                                                "variant_type": evidence.get(
                                                    "variant_type", "unknown"
                                                ).upper(),
                                                "clinical_significance": variant.get(
                                                    "clinical_significance", "unknown"
                                                )
                                                .replace("_", " ")
                                                .title(),
                                                "disease": variant.get(
                                                    "disease_name", "unknown"
                                                )[:50]
                                                + (
                                                    "..."
                                                    if len(variant.get("disease_name", ""))
                                                    > 50
                                                    else ""
                                                ),
                                                "coverage": f"{variant.get('coverage_at_variant', 0)}x",
                                                "alt_support": evidence.get(
                                                    "alternate_support", 0
                                                ),
                                                "vaf": (
                                                    f"{evidence.get('variant_allele_frequency', 0):.1%}"
                                                    if evidence.get(
                                                        "variant_allele_frequency"
                                                    )
                                                    is not None
                                                    else "N/A"
                                                ),
                                                "evidence_level": evidence.get(
                                                    "evidence_level", "unknown"
                                                )
                                                .replace("_", " ")
                                                .title(),
                                            }
                                        )

                                    if variant_data:
                                        # Create table with proper columns using pandas DataFrame
                                        df = pd.DataFrame(variant_data)
                                        # Create styled table for consistency
                                        from robin.gui.theme import styled_table
                                        
                                        # Create columns definition
                                        columns = []
                                        for col in df.columns:
                                            columns.append({
                                                "name": col,
                                                "label": col.replace("_", " ").title(),
                                                "field": col,
                                                "sortable": True
                                            })
                                        
                                        variant_table = _render_paged_df_table(
                                            df,
                                            columns,
                                            pagination_default=100,
                                            class_size="table-xs",
                                            search_placeholder="Search variants...",
                                        )

                                        # Make columns sortable
                                        for col in variant_table.columns:
                                            col["sortable"] = True

                                        if len(high_cov_variants) > 20:
                                            ui.label(
                                                f"Showing first 20 of {len(high_cov_variants)} high coverage variants"
                                            ).classes("text-xs text-gray-500 mt-2")

                    except Exception as e:
                        ui.label(f"Error displaying high coverage genes: {e}").classes(
                            "text-sm text-red-500"
                        )
                        import traceback

                        traceback.print_exc()

                def _display_low_coverage_genes(data):
                    """Display genes with only low coverage variants"""
                    try:
                        genes = data.get("genes", {})

                        # Filter for genes with only low coverage variants
                        low_coverage_genes = {}
                        for gene_name, gene_data in genes.items():
                            summary = gene_data.get("summary", {})
                            if (
                                summary.get("high_coverage_variants", 0) == 0
                                and summary.get("total_variants", 0) > 0
                            ):
                                low_coverage_genes[gene_name] = gene_data

                        if not low_coverage_genes:
                            ui.label(
                                "No genes with only low coverage variants found."
                            ).classes("text-sm text-gray-500")
                            return

                        ui.label(
                            f"Found {len(low_coverage_genes)} genes with only low coverage variants"
                        ).classes("text-lg font-semibold mb-4")

                        # Create summary table
                        low_cov_data = []
                        for gene_name, gene_data in low_coverage_genes.items():
                            summary = gene_data.get("summary", {})
                            coverage_stats = gene_data.get("coverage_statistics", {})

                            low_cov_data.append(
                                {
                                    "gene": gene_name,
                                    "mean_coverage": f"{coverage_stats.get('mean_coverage', 0):.1f}x",
                                    "total_variants": summary.get("total_variants", 0),
                                    "low_coverage_variants": summary.get(
                                        "low_coverage_variants", 0
                                    ),
                                    "status": summary.get("pathogenic_status", "unknown")
                                    .replace("_", " ")
                                    .title(),
                                }
                            )

                        # Sort by mean coverage (lowest first)
                        low_cov_data.sort(
                            key=lambda x: float(x["mean_coverage"].replace("x", ""))
                        )

                        # Create table using styled_table for consistency
                        from robin.gui.theme import styled_table
                        
                        # Create columns definition
                        columns = []
                        for col in low_cov_data[0].keys():
                            columns.append({
                                "name": col,
                                "label": col.replace("_", " ").title(),
                                "field": col,
                                "sortable": True
                            })
                        
                        # Create styled table
                        table_container, low_cov_table = styled_table(
                            columns=columns,
                            rows=low_cov_data,
                            pagination=25,
                            class_size="table-xs"
                        )

                        # Add search functionality
                        with low_cov_table.add_slot("top-right"):
                            with ui.input(placeholder="Search genes...").props(
                                "type=search"
                            ).bind_value(low_cov_table, "filter").add_slot("append"):
                                ui.icon("search")

                        # Make columns sortable
                        for col in low_cov_table.columns:
                            col["sortable"] = True

                    except Exception as e:
                        ui.label(f"Error displaying low coverage genes: {e}").classes(
                            "text-sm text-red-500"
                        )
                        import traceback

                        traceback.print_exc()

                def _display_gene_details(data):
                    """Display detailed gene information with variant analysis"""
                    try:
                        genes = data.get("genes", {})

                        # Gene selector
                        with ui.row().classes("w-full items-center gap-3 mb-4"):
                            ui.label("Select Gene:").classes("text-sm font-medium")
                            gene_selector = ui.select(
                                options=list(genes.keys()),
                                value=list(genes.keys())[0] if genes else None,
                                label="Gene",
                            ).classes("w-64")

                        # Gene details display
                        gene_details_container = ui.column().classes("w-full")

                        def update_gene_details():
                            gene_details_container.clear()

                            selected_gene = gene_selector.value
                            if not selected_gene or selected_gene not in genes:
                                return

                            gene_data = genes[selected_gene]
                            summary = gene_data.get("summary", {})
                            coverage_stats = gene_data.get("coverage_statistics", {})
                            variants = gene_data.get("variants", [])

                            with gene_details_container:
                                # Gene header
                                with ui.card().classes("w-full mb-4"):
                                    ui.label(
                                        f"{selected_gene} - Gene Analysis Details"
                                    ).classes("text-xl font-semibold mb-2")

                                    # Linear layout for gene stats
                                    with ui.column().classes("w-full gap-3"):
                                        with ui.row().classes(
                                            "w-full items-center justify-between p-3 bg-blue-50 rounded"
                                        ):
                                            ui.label("Mean Coverage").classes(
                                                "text-base font-medium text-gray-700"
                                            )
                                            ui.label(
                                                f"{coverage_stats.get('mean_coverage', 0):.1f}x"
                                            ).classes("text-2xl font-bold text-blue-600")

                                        with ui.row().classes(
                                            "w-full items-center justify-between p-3 bg-purple-50 rounded"
                                        ):
                                            ui.label("Total Variants").classes(
                                                "text-base font-medium text-gray-700"
                                            )
                                            ui.label(
                                                f"{summary.get('total_variants', 0)}"
                                            ).classes("text-2xl font-bold text-purple-600")

                                        with ui.row().classes(
                                            "w-full items-center justify-between p-3 bg-green-50 rounded"
                                        ):
                                            ui.label(
                                                "High Coverage Variants (≥10x)"
                                            ).classes("text-base font-medium text-gray-700")
                                            ui.label(
                                                f"{summary.get('high_coverage_variants', 0)}"
                                            ).classes("text-2xl font-bold text-green-600")

                                        with ui.row().classes(
                                            "w-full items-center justify-between p-3 bg-orange-50 rounded"
                                        ):
                                            ui.label(
                                                "Low Coverage Variants (<10x)"
                                            ).classes("text-base font-medium text-gray-700")
                                            ui.label(
                                                f"{summary.get('low_coverage_variants', 0)}"
                                            ).classes("text-2xl font-bold text-orange-600")

                                # Variants table
                                if variants:
                                    ui.label("Variant Details").classes(
                                        "text-lg font-semibold mb-2"
                                    )

                                    # Create comprehensive variants table - only show variants where alternate allele is detected
                                    variant_data = []
                                    for variant in variants:
                                        evidence = variant.get("evidence_analysis", {})
                                        # Only include variants where alternate allele is actually detected
                                        if evidence.get("alternate_support", 0) > 0:
                                            # Create genomic locus in format chr1:119,915,631
                                            chromosome = variant.get("chromosome", "N")
                                            position = variant.get("position", 0)
                                            genomic_locus = (
                                                f"{chromosome}:{position:,}"
                                                if chromosome != "N" and position
                                                else "N/A"
                                            )

                                            variant_data.append(
                                                {
                                                    "genomic_locus": genomic_locus,
                                                    "position": f"{position:,}",
                                                    "chromosome": chromosome,
                                                    "reference": variant.get(
                                                        "reference", "N"
                                                    ),
                                                    "alternate": variant.get(
                                                        "alternate", "N"
                                                    ),
                                                    "variant_type": evidence.get(
                                                        "variant_type", "unknown"
                                                    ).upper(),
                                                    "clinical_significance": variant.get(
                                                        "clinical_significance", "unknown"
                                                    )
                                                    .replace("_", " ")
                                                    .title(),
                                                    "disease": variant.get(
                                                        "disease_name", "unknown"
                                                    )[:40]
                                                    + (
                                                        "..."
                                                        if len(
                                                            variant.get("disease_name", "")
                                                        )
                                                        > 40
                                                        else ""
                                                    ),
                                                    "coverage": f"{variant.get('coverage_at_variant', 0)}x",
                                                    "coverage_status": (
                                                        "High"
                                                        if variant.get(
                                                            "coverage_at_variant", 0
                                                        )
                                                        >= 10
                                                        else "Low"
                                                    ),
                                                    "alt_support": evidence.get(
                                                        "alternate_support", 0
                                                    ),
                                                    "evidence_level": evidence.get(
                                                        "evidence_level", "unknown"
                                                    )
                                                    .replace("_", " ")
                                                    .title(),
                                                    "zygosity": evidence.get(
                                                        "zygosity", "unknown"
                                                    )
                                                    .replace("_", " ")
                                                    .title(),
                                                    "vaf": (
                                                        f"{evidence.get('variant_allele_frequency', 0):.1%}"
                                                        if evidence.get(
                                                            "variant_allele_frequency"
                                                        )
                                                        is not None
                                                        else "N/A"
                                                    ),
                                                }
                                            )

                                    # Create table with proper columns using pandas DataFrame
                                    df = pd.DataFrame(variant_data)
                                    # Create styled table for consistency
                                    from robin.gui.theme import styled_table
                                    
                                    # Create columns definition
                                    columns = []
                                    for col in df.columns:
                                        columns.append({
                                            "name": col,
                                            "label": col.replace("_", " ").title(),
                                            "field": col,
                                            "sortable": True
                                        })
                                    
                                    detailed_variant_table = _render_paged_df_table(
                                        df,
                                        columns,
                                        pagination_default=100,
                                        class_size="table-xs",
                                        search_placeholder="Search variants...",
                                    )

                                    # Make columns sortable
                                    for col in detailed_variant_table.columns:
                                        col["sortable"] = True

                                    # Add export functionality
                                    with ui.row().classes("w-full mt-4"):
                                        ui.button(
                                            "Export Gene Variants to CSV",
                                            icon="download",
                                            on_click=lambda: _export_gene_variants(
                                                variant_data, selected_gene
                                            ),
                                        ).classes("w-full")
                                else:
                                    ui.label("No variants found for this gene.").classes(
                                        "text-sm text-gray-500"
                                    )

                        # Initial display
                        update_gene_details()

                        # Update when gene selection changes
                        gene_selector.on("change", update_gene_details)

                    except Exception as e:
                        ui.label(f"Error displaying gene details: {e}").classes(
                            "text-sm text-red-500"
                        )
                        import traceback

                        traceback.print_exc()

                # Helper function to export gene variants
                def _export_gene_variants(variant_data, gene_name):
                    """Export gene variants to CSV"""
                    try:
                        import tempfile

                        # Create temporary file
                        with tempfile.NamedTemporaryFile(
                            mode="w", suffix=".csv", delete=False
                        ) as f:
                            df = pd.DataFrame(variant_data)
                            df.to_csv(f.name, index=False)
                            temp_path = f.name

                        # Download the file
                        ui.download(temp_path, filename=f"{gene_name}_variants.csv")

                        # Clean up
                        os.unlink(temp_path)

                        ui.notify(f"{gene_name} variants exported successfully!")

                    except Exception as e:
                        ui.notify(f"Export failed: {e}", type="error")

                # Check for existing results
                _check_lga_results()

        # Function to check file availability for lightweight variant analysis
        def _check_lga_file_availability():
            try:
                target_bam = sample_dir / "target.bam"
                targets_bed = sample_dir / "targets_exceeding_threshold.bed"

                if not target_bam.exists():
                    lga_files_status.set_text(
                        " Missing target.bam - run target analysis first"
                    )
                    lga_files_status.classes(replace="text-xs text-red-500")
                    lga_analysis_button.disable()
                    return False

                if not targets_bed.exists():
                    lga_files_status.set_text(
                        " Missing targets_exceeding_threshold.bed - run target analysis first"
                    )
                    lga_files_status.classes(replace="text-xs text-red-500")
                    lga_analysis_button.disable()
                    return False

                # Check if BED file has content
                try:
                    with open(targets_bed, "r") as f:
                        bed_content = f.read().strip()

                    if not bed_content:
                        lga_files_status.set_text(
                            "BED file is empty - no regions exceed threshold"
                        )
                        lga_files_status.classes(replace="text-xs text-yellow-500")
                        lga_analysis_button.disable()
                        return False
                    else:
                        lga_files_status.set_text("All required files available")
                        lga_files_status.classes(replace="text-xs text-green-500")
                        lga_analysis_button.enable()
                        return True

                except Exception as e:
                    lga_files_status.set_text(f" Error reading BED file: {e}")
                    lga_files_status.classes(replace="text-xs text-red-500")
                    lga_analysis_button.disable()
                    return False

            except Exception as e:
                lga_files_status.set_text(f" Error checking files: {e}")
                lga_files_status.classes(replace="text-xs text-red-500")
                lga_analysis_button.disable()
                return False

        # Check file availability initially
        _check_lga_file_availability()

        # Function to trigger lightweight variant analysis
        def _trigger_lightweight_gene_analysis():
            try:
                # Check prerequisites
                target_bam = sample_dir / "target.bam"
                targets_bed = sample_dir / "targets_exceeding_threshold.bed"

                if not target_bam.exists():
                    ui.notify(
                        "Missing target.bam file. Please run target analysis first.",
                        type="warning",
                    )
                    return

                if not targets_bed.exists():
                    ui.notify(
                        "Missing targets_exceeding_threshold.bed file. Please run target analysis first.",
                        type="warning",
                    )
                    return

                # Check if BED file has content (not empty)
                try:
                    with open(targets_bed, "r") as f:
                        bed_content = f.read().strip()

                    if not bed_content:
                        ui.notify(
                            "targets_exceeding_threshold.bed file is empty. No regions exceed coverage threshold.",
                            type="warning",
                        )
                        return

                except Exception as e:
                    ui.notify(
                        f"Error reading targets_exceeding_threshold.bed file: {e}",
                        type="error",
                    )
                    return

                # Update UI state
                lga_analysis_button.disable()
                lga_status_label.set_text("Starting lightweight variant analysis...")
                lga_status_label.classes(replace="text-sm text-blue-600")

                print(" GUI: Starting lightweight variant analysis...")

                # Create and submit lightweight variant analysis job
                try:
                    # Import required modules
                    print(" GUI: Importing required modules...")
                    from robin.analysis.lightweight_gene_analysis import (
                        lightweight_gene_analysis_handler,
                    )
                    from robin.workflow_simple import Job, WorkflowContext

                    # Get reference genome from workflow runner (set by CLI)
                    reference_genome = None

                    if (
                        hasattr(launcher, "workflow_runner")
                        and launcher.workflow_runner
                    ):
                        if hasattr(launcher.workflow_runner, "reference"):
                            reference_genome = launcher.workflow_runner.reference

                    # Fallback to environment variable
                    if not reference_genome:
                        env_reference = os.environ.get("robin_REFERENCE")
                        if env_reference and os.path.exists(env_reference):
                            reference_genome = env_reference

                    print(f" GUI: Reference genome: {reference_genome}")

                    # Use monitored_directory as work_dir, or fall back to sample_dir if not available
                    work_dir = (
                        launcher.monitored_directory
                        if launcher.monitored_directory
                        else str(sample_dir)
                    )
                    print(f" GUI: Work directory: {work_dir}")

                    metadata = {
                        "work_dir": work_dir,
                        "bed_path": str(targets_bed),
                        "reference": reference_genome,
                        "force_regenerate": lga_force_regenerate_checkbox.value,
                    }

                    print(f" GUI: Job metadata: {metadata}")

                    # Create workflow context
                    print(" GUI: Creating workflow context...")
                    context = WorkflowContext(
                        filepath=str(sample_dir), metadata=metadata
                    )

                    # Add sample_id method
                    def get_sample_id():
                        return sample_dir.name

                    context.get_sample_id = get_sample_id
                    context.add_result = lambda key, value: None  # Mock for now
                    context.add_error = lambda key, value: None  # Mock for now

                    # Create job
                    print(" GUI: Creating job object...")
                    job = Job(
                        job_id=hash(f"lightweight_gene_analysis_{sample_dir.name}")
                        % 1000000,  # Simple hash-based ID
                        job_type="lightweight_gene_analysis",
                        context=context,
                        origin="fast",
                        workflow=["fast:lightweight_gene_analysis"],
                    )

                    print(f" GUI: Job created with ID: {job.job_id}")

                    # For now, let's force direct execution to see our debug prints
                    # TODO: Re-enable workflow submission once we understand how it works
                    print(" Forcing direct execution to see debug output...")

                    # Use a shared variable to communicate between threads
                    lga_status = {"status": "running", "message": "", "error": None}

                    def run_lga_analysis_direct():
                        try:
                            print(
                                "Starting lightweight gene analysis in background thread..."
                            )
                            # Use the work_dir variable we calculated earlier
                            lightweight_gene_analysis_handler(job, work_dir=work_dir)

                            # Update UI on completion using the main thread
                            print(
                                "Lightweight variant analysis completed successfully!"
                            )
                            lga_status["status"] = "completed"
                            lga_status["message"] = (
                                "Lightweight variant analysis completed successfully!"
                            )

                        except Exception as e:
                            print(f" Error in background thread: {e}")
                            lga_status["status"] = "error"
                            lga_status["error"] = str(e)

                    # Run in background thread
                    import threading

                    thread = threading.Thread(
                        target=run_lga_analysis_direct, daemon=True
                    )
                    thread.start()

                    # Check for completion status periodically
                    def check_lga_status():
                        if lga_status["status"] == "completed":
                            try:
                                lga_status_label.set_text(lga_status["message"])
                                lga_status_label.classes(
                                    replace="text-sm text-green-600"
                                )
                                lga_analysis_button.enable()
                                lga_analysis_button.set_text("Rerun Variant Analysis")
                                lga_analysis_button.props("color=secondary")
                            except Exception as e:
                                print(f" Error updating UI: {e}")
                            return False  # Stop checking
                        elif lga_status["status"] == "error":
                            try:
                                lga_status_label.set_text(
                                    f"Lightweight variant analysis failed: {lga_status['error']}"
                                )
                                lga_status_label.classes(replace="text-sm text-red-600")
                                lga_analysis_button.enable()
                            except Exception as e:
                                print(f" Error updating UI: {e}")
                            return False  # Stop checking
                        return True  # Keep checking

                    # Start checking status every 0.5 seconds
                    ui.timer(0.5, check_lga_status, active=True)

                    ui.notify(
                        "Lightweight variant analysis started in background (direct execution)",
                        type="info",
                    )

                    # Comment out the old workflow submission code for now
                    # TODO: Re-enable workflow submission once we understand how it works
                    pass

                except Exception as e:
                    ui.notify(
                        f"Failed to start lightweight variant analysis: {e}",
                        type="error",
                    )
                    lga_analysis_button.enable()
                    lga_status_label.set_text(
                        "Failed to start lightweight variant analysis"
                    )
                    lga_status_label.classes(replace="text-sm text-red-600")

            except Exception as e:
                ui.notify(
                    f"Error triggering lightweight variant analysis: {e}", type="error"
                )
                lga_analysis_button.enable()
                lga_status_label.set_text(
                    "Error triggering lightweight variant analysis"
                )
                lga_status_label.classes(replace="text-sm text-red-600")

    # Helpers
    def _log_notify(message: str, level: str = "warning", notify: bool = False) -> None:
        try:
            # Python logging
            log_func = getattr(logging, level, logging.warning)
            log_func(f"[coverage] {message}")
            # Append to GUI log buffer if available
            try:
                launcher._log_buffer.append(f"[coverage] {message}\n")  # type: ignore[attr-defined]
                if hasattr(launcher, "log_area") and launcher.log_area is not None:  # type: ignore[attr-defined]
                    launcher.log_area.set_value("".join(launcher._log_buffer))  # type: ignore[attr-defined]
            except Exception:
                pass
            # User notification (optional)
            if notify and ui is not None:
                n_type = "warning" if level in {"warning", "debug", "info"} else "error"
                try:
                    ui.notify(message, type=n_type)
                except Exception:
                    pass
        except Exception:
            # Never let logging itself break the UI
            pass

    def _update_target_cov(cov_df: pd.DataFrame, bed_df: pd.DataFrame) -> None:
        try:

            def chr_key(label: str) -> int:
                s = str(label)
                if s.startswith("chr"):
                    s = s[3:]
                mapping = {"X": 23, "Y": 24}
                try:
                    return int(s)
                except Exception:
                    return mapping.get(s, 1000)

            bed_df = bed_df.copy()
            bed_df["length"] = (bed_df["endpos"] - bed_df["startpos"] + 1).astype(float)
            grouped = (
                bed_df.groupby("chrom")
                .agg({"bases": "sum", "length": "sum"})
                .reset_index()
            )
            grouped["meandepth"] = grouped["bases"] / grouped["length"]
            pattern = r"^chr([0-9]+|X|Y)$"
            name_col = "#rname" if "#rname" in cov_df.columns else "rname"
            temp_df = cov_df[cov_df[name_col].astype(str).str.match(pattern)]
            temp_df = temp_df[temp_df[name_col] != "chrM"]
            names = sorted(temp_df[name_col].astype(str).unique(), key=chr_key)
            
            if not names:
                _log_notify(
                    "No valid chromosome data found for target coverage chart",
                    level="warning",
                    notify=False,
                )
                return
            
            names_short = [_short_chr_label(n) for n in names]
            echart_target_cov.options["xAxis"]["data"] = names_short
            grouped = grouped.set_index("chrom").reindex(names).reset_index()
            
            # Get off-target data - group by chromosome and take mean if multiple entries
            temp_df_grouped = temp_df.groupby(name_col)["meandepth"].mean()
            off_target_data = [
                float(temp_df_grouped.get(chrom, 0.0)) if pd.notna(temp_df_grouped.get(chrom, 0.0)) else 0.0
                for chrom in names
            ]
            
            # Get on-target data
            on_target_data = [
                float(v) if pd.notna(v) else 0.0 
                for v in grouped["meandepth"].fillna(0).tolist()
            ]
            
            _on, _off = _cov_on_off_colors()
            echart_target_cov.options["series"] = [
                {
                    "type": "bar",
                    "name": "On Target",
                    "barWidth": "32%",
                    "barGap": "15%",
                    "data": on_target_data,
                    "itemStyle": {"color": _on},
                },
                {
                    "type": "bar",
                    "name": "Off Target",
                    "barWidth": "32%",
                    "data": off_target_data,
                    "itemStyle": {"color": _off},
                },
            ]
            _apply_coverage_target_cov_chrome(echart_target_cov)
            echart_target_cov.update()
        except Exception as e:
            _log_notify(
                f"Target vs Off-target coverage update failed: {e}",
                level="warning",
                notify=False,
            )

    def _echarts_option_to_json(obj: Any) -> Any:
        """Convert ECharts options to JSON-serializable form (numpy types, NaN)."""
        if isinstance(obj, dict):
            return {k: _echarts_option_to_json(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [_echarts_option_to_json(x) for x in obj]
        if isinstance(obj, (np.floating, np.integer)):
            v = float(obj) if isinstance(obj, np.floating) else int(obj)
            return None if (isinstance(obj, np.floating) and np.isnan(obj)) else v
        if isinstance(obj, float) and np.isnan(obj):
            return None
        return obj

    def _update_boxplot(bed_df: pd.DataFrame, panel_display_name: str = "Target Coverage") -> None:
        try:
            df = bed_df.copy()
            if "coverage" not in df.columns:
                df["length"] = (df["endpos"] - df["startpos"] + 1).astype(float)
                df["coverage"] = df["bases"] / df["length"]
            # Normalize chrom to canonical form so "chr1" and "1" map to one category per chromosome
            chrom_str = df["chrom"].astype(str).str.strip()
            df["chrom"] = chrom_str.where(chrom_str.str.startswith("chr"), "chr" + chrom_str)
            chroms = natsort.natsorted(df["chrom"].unique())
            chrom_lookup = {chrom: idx for idx, chrom in enumerate(chroms)}
            df["chrom_index"] = df["chrom"].map(chrom_lookup)
            agg = (
                df.groupby("chrom")
                .agg(
                    min=("coverage", "min"),
                    Q1=("coverage", lambda x: float(np.percentile(x, 25))),
                    median=("coverage", "median"),
                    Q3=("coverage", lambda x: float(np.percentile(x, 75))),
                    max=("coverage", "max"),
                    chrom_index=("chrom_index", "first"),
                )
                .reset_index()
            )
            agg["chrom"] = pd.Categorical(agg["chrom"], categories=chroms, ordered=True)
            agg = agg.sort_values("chrom").reset_index(drop=True)
            result = [
                ["chrom", "min", "Q1", "median", "Q3", "max", "chrom_index"]
            ] + agg.values.tolist()

            def iqr_bounds(sub):
                q1 = np.percentile(sub["coverage"], 25)
                q3 = np.percentile(sub["coverage"], 75)
                iqr = q3 - q1
                return q1 - 1.5 * iqr, q3 + 1.5 * iqr

            out_rows = []
            for c in chroms:
                sub = df[df["chrom"].astype(str) == c]
                lb, ub = iqr_bounds(sub)
                out = sub[(sub["coverage"] < lb) | (sub["coverage"] > ub)]
                out_rows += out[["chrom", "coverage", "name"]].values.tolist()
            lb_g, ub_g = iqr_bounds(df)
            glob = df[(df["coverage"] < lb_g) | (df["coverage"] > ub_g)]
            glob_rows = glob[["chrom", "coverage", "name"]].values.tolist()
            target_boxplot.options["dataset"][0]["source"] = result
            target_boxplot.options["dataset"][1]["source"] = df[
                ["chrom", "coverage", "name"]
            ].values.tolist()
            target_boxplot.options["dataset"][2]["source"] = out_rows
            target_boxplot.options["dataset"][3]["source"] = glob_rows
            target_boxplot.options["xAxis"]["data"] = chroms
            target_boxplot.options["xAxis"]["name"] = "Chromosome"
            target_boxplot.options["xAxis"]["nameLocation"] = "middle"
            target_boxplot.options["xAxis"]["nameGap"] = 40

            # Restore original series configuration (explicit ids so ECharts replaceMerge
            # replaces series correctly on update, avoiding duplicate boxes per chromosome)
            target_boxplot.options["series"] = [
                {
                    "id": "target_box_series",
                    "name": "box plot",
                    "type": "boxplot",
                    "itemStyle": _cov_boxplot_item_style(),
                    "datasetId": "raw",
                    "universalTransition": True,
                    "animationDurationUpdate": 1000,
                    "encode": {
                        "x": "chrom",
                        "y": ["min", "Q1", "median", "Q3", "max"],
                        "itemName": ["chrom"],
                        "tooltip": ["min", "Q1", "median", "Q3", "max"],
                    },
                },
                {
                    "id": "target_outliers_series",
                    "name": "outliers",
                    "type": "scatter",
                    "datasetId": "outliers",
                    "symbolSize": 6,
                    "itemStyle": _cov_scatter_outlier_style(),
                    "label": {
                        "show": True,
                        "position": "right",
                        "formatter": "{@name}",
                    },
                    "encode": {
                        "x": "chrom",
                        "y": "coverage",
                        "label": ["name"],
                        "tooltip": ["name", "coverage"],
                    },
                },
                {
                    "id": "target_global_outliers_series",
                    "name": "global outliers",
                    "type": "scatter",
                    "datasetId": "globaloutliers",
                    "symbolSize": 6,
                    "itemStyle": _cov_scatter_outlier_style(),
                    "label": {
                        "show": True,
                        "position": "right",
                        "formatter": "{@name}",
                    },
                    "encode": {
                        "x": "chrom",
                        "y": "coverage",
                        "label": ["name"],
                        "tooltip": ["name", "coverage"],
                    },
                },
            ]
            
            # Restore original title and legend
            target_boxplot.options["title"]["text"] = f"Target Coverage ({panel_display_name})"
            target_boxplot.options["title"]["subtext"] = ""
            target_boxplot.options["legend"]["show"] = True
            _apply_coverage_boxplot_chrome(target_boxplot)
            
            target_boxplot.update()

            # Force ECharts to replace series/dataset instead of merging, so we don't get
            # duplicate boxes per chromosome on first load (merge can append series).
            try:
                options_clean = _echarts_option_to_json(target_boxplot.options)
                options_json = json.dumps(options_clean)
                options_escaped = json.dumps(options_json)
                ui.run_javascript(
                    f"""
                    (function() {{
                      var el = document.querySelector('.target-coverage-boxplot');
                      if (!el) return;
                      var chart = echarts.getInstanceByDom(el);
                      if (!chart) {{
                        var child = el.querySelector('div');
                        if (child) chart = echarts.getInstanceByDom(child);
                      }}
                      if (!chart) return;
                      try {{
                        var options = JSON.parse({options_escaped});
                        chart.setOption(options, {{ replaceMerge: ['series', 'dataset'] }});
                      }} catch (e) {{ console.warn('Target coverage replaceMerge:', e); }}
                    }})();
                    """
                )
            except Exception as js_err:
                logging.debug("Target coverage replaceMerge JS skip: %s", js_err)
        except Exception as e:
            _log_notify(
                f"Target boxplot update failed: {e}", level="warning", notify=False
            )

    def _update_time(npy_path: Path) -> None:
        try:
            arr = np.load(npy_path)
            raw = arr.tolist()
            # Ensure [[ts_ms, cov], ...] for time axis
            echart_time.options["series"][0]["data"] = raw
            _apply_coverage_time_series_chrome(echart_time)
            echart_time.update()
        except Exception as e:
            _log_notify(
                f"Coverage-over-time update failed: {e}", level="warning", notify=False
            )

    def _update_target_table(df: pd.DataFrame) -> None:
        try:
            dfx = df.copy()
            if "coverage" not in dfx.columns:
                dfx["length"] = (dfx["endpos"] - dfx["startpos"] + 1).astype(float)
                dfx["coverage"] = dfx["bases"] / dfx["length"]
            dfx["coverage"] = dfx["coverage"].astype(float).round(2)

            # Compute outlier status using 2 SD method (matches Target Coverage Over Time plot)
            outlier_status = []
            if len(dfx) >= 3:
                mean_cov = dfx["coverage"].mean()
                std_cov = dfx["coverage"].std()
                if std_cov > 0:
                    lower = mean_cov - 2.0 * std_cov
                    upper = mean_cov + 2.0 * std_cov
                    for cov in dfx["coverage"]:
                        if cov > upper:
                            outlier_status.append("High")
                        elif cov < lower:
                            outlier_status.append("Low")
                        else:
                            outlier_status.append("")
                else:
                    outlier_status = [""] * len(dfx)
            else:
                outlier_status = [""] * len(dfx)

            dfx["outlier_status"] = outlier_status
            rows = dfx[
                ["chrom", "startpos", "endpos", "name", "coverage", "outlier_status"]
            ].to_dict(orient="records")
            target_cov_table.rows = rows
            target_cov_table.update()
        except Exception as e:
            _log_notify(
                f"Target table update failed: {e}", level="warning", notify=False
            )

    async def _refresh_coverage_async() -> None:
        """Refresh coverage data asynchronously."""
        try:
            # Check directory existence
            if not sample_dir or not sample_dir.exists():
                _log_notify(
                    f"Sample directory not found: {sample_dir}",
                    level="warning",
                    notify=True,
                )
                return
            
            # Run coverage refresh synchronously (already in background)
            _refresh_coverage_sync(sample_dir, launcher)
                
        except Exception as e:
            _log_notify(
                f"Unexpected coverage refresh error: {e}", level="error", notify=True
            )

    def _refresh_coverage_sync(sample_dir: Path, launcher: Any) -> None:
        """Synchronous coverage refresh - runs in background thread."""
        try:
            key = str(sample_dir)
            state = launcher._coverage_state.get(key, {})
            
            # Check if this is a fresh page visit (no last_visit_time means new page)
            # Force updates on fresh page visits regardless of mtime
            is_fresh_visit = "last_visit_time" not in state
            if is_fresh_visit:
                state["last_visit_time"] = time.time()
            
            cov_main = sample_dir / "coverage_main.csv"
            bed_cov = sample_dir / "bed_coverage_main.csv"
            target_cov = sample_dir / "target_coverage.csv"
            cov_time = sample_dir / "coverage_time_chart.npy"
            
            # Check modification times for all files
            cov_main_mtime = cov_main.stat().st_mtime if cov_main.exists() else 0
            bed_cov_mtime = bed_cov.stat().st_mtime if bed_cov.exists() else 0
            target_cov_mtime = target_cov.stat().st_mtime if target_cov.exists() else 0
            cov_time_mtime = cov_time.stat().st_mtime if cov_time.exists() else 0
            
            # Get previous state values (use sentinel values if not present)
            prev_cov_main_mtime = state.get("cov_main_mtime", 0)
            prev_bed_cov_mtime = state.get("bed_cov_mtime", 0)
            prev_target_cov_mtime = state.get("target_cov_mtime", 0)
            prev_cov_time_mtime = state.get("cov_time_mtime", 0)
            
            # Check if any file has changed
            cov_main_changed = prev_cov_main_mtime != cov_main_mtime
            bed_cov_changed = prev_bed_cov_mtime != bed_cov_mtime
            target_cov_changed = prev_target_cov_mtime != target_cov_mtime
            cov_time_changed = prev_cov_time_mtime != cov_time_mtime
            
            # Also check if data is missing from state (needs initial load)
            data_missing = (
                (cov_main.exists() and "cov_df" not in state) or
                (bed_cov.exists() and "bed_df" not in state)
            )
            
            # Determine if any update is needed
            needs_update = (
                is_fresh_visit
                or cov_main_changed
                or bed_cov_changed
                or target_cov_changed
                or cov_time_changed
                or data_missing
            )
            
            # Early exit if nothing has changed (except for periodic checks that run less frequently)
            if not needs_update:
                # Still check for SNP/LGA/IGV status updates (these have their own throttling)
                # But skip all file reading and UI updates
                logging.debug(f"[Coverage] ⏭ Skipping coverage update - no file changes detected")
                # Continue to SNP/LGA/IGV checks below (they have their own throttling)
                # But we'll skip the file processing section
            else:
                # Log what changed
                reasons = []
                if is_fresh_visit:
                    reasons.append("fresh_visit")
                if cov_main_changed:
                    reasons.append("coverage_main.csv")
                if bed_cov_changed:
                    reasons.append("bed_coverage_main.csv")
                if target_cov_changed:
                    reasons.append("target_coverage.csv")
                if cov_time_changed:
                    reasons.append("coverage_time_chart.npy")
                if data_missing:
                    reasons.append("data_missing")
                logging.debug(f"[Coverage] Update needed. Reasons: {', '.join(reasons)}")
            # Only process files if update is needed
            if needs_update:
                if cov_main.exists():
                    m = cov_main_mtime
                    # Force update on fresh page visit or if mtime changed or data not in state
                    if is_fresh_visit or cov_main_changed or "cov_df" not in state:
                        try:
                            cov_df = pd.read_csv(cov_main)
                            # Backfill meandepth if missing
                            if "meandepth" not in cov_df.columns and {
                                "covbases",
                                "endpos",
                            }.issubset(set(cov_df.columns)):
                                cov_df = cov_df.copy()
                                with np.errstate(divide="ignore", invalid="ignore"):
                                    cov_df["meandepth"] = (
                                        cov_df["covbases"]
                                        / cov_df["endpos"].replace(0, np.nan)
                                    ).fillna(0)
                            state["cov_df"] = cov_df
                            state["cov_main_mtime"] = m  # only set on success
                        except Exception as e:
                            _log_notify(
                                f"Failed to read coverage_main.csv: {e}",
                                level="error",
                                notify=True,
                            )
                else:
                    _log_notify("coverage_main.csv not found", level="debug", notify=False)
                if bed_cov.exists():
                    m = bed_cov_mtime
                    # Force update on fresh page visit or if mtime changed or data not in state
                    if is_fresh_visit or bed_cov_changed or "bed_df" not in state:
                        try:
                            bed_df = pd.read_csv(bed_cov)
                            state["bed_df"] = bed_df
                            panel_display_name = state.get("panel_display_name", "Target Coverage")
                            _update_boxplot(bed_df, panel_display_name)
                            if state.get("cov_df") is not None:
                                _update_target_cov(state["cov_df"], bed_df)
                            state["bed_cov_mtime"] = m  # only set on success
                        except Exception as e:
                            _log_notify(
                                f"Failed to read bed_coverage_main.csv: {e}",
                                level="error",
                                notify=True,
                            )
                else:
                    _log_notify(
                        "bed_coverage_main.csv not found", level="debug", notify=False
                    )
                if target_cov.exists():
                    m = target_cov_mtime
                    # Force update on fresh page visit or if mtime changed
                    if is_fresh_visit or target_cov_changed:
                        try:
                            tdf = pd.read_csv(target_cov)
                            _update_target_table(tdf)
                            state["target_cov_mtime"] = m  # only set on success
                        except Exception as e:
                            _log_notify(
                                f"Failed to read target_coverage.csv: {e}",
                                level="error",
                                notify=True,
                            )
                else:
                    _log_notify(
                        "target_coverage.csv not found", level="debug", notify=False
                    )
                if cov_time.exists():
                    m = cov_time_mtime
                    # Force update on fresh page visit or if mtime changed
                    if is_fresh_visit or cov_time_changed:
                        try:
                            _update_time(cov_time)
                            state["cov_time_mtime"] = m
                        except Exception as e:
                            _log_notify(
                                f"Failed to load coverage_time_chart.npy: {e}",
                                level="warning",
                                notify=False,
                            )
                # Summary - only recalculate if data changed
                try:
                    global_cov = None
                    target_cov_v = None
                    enrich_v = None
                    if state.get("cov_df") is not None:
                        cdf = state["cov_df"]
                        if (
                            "covbases" in cdf.columns
                            and "endpos" in cdf.columns
                            and cdf["endpos"].sum() > 0
                        ):
                            global_cov = float(cdf["covbases"].sum()) / float(
                                cdf["endpos"].sum()
                            )
                    if state.get("bed_df") is not None:
                        bdf = state["bed_df"].copy()
                        if "length" not in bdf.columns:
                            bdf["length"] = (bdf["endpos"] - bdf["startpos"] + 1).astype(
                                float
                            )
                        if bdf["length"].sum() > 0:
                            target_cov_v = float(bdf["bases"].sum()) / float(
                                bdf["length"].sum()
                            )
                    if (
                        global_cov is not None
                        and target_cov_v is not None
                        and global_cov > 0
                    ):
                        enrich_v = target_cov_v / global_cov
                    if target_cov_v is not None:
                        if target_cov_v >= 30:
                            q_text, q_cls, q_bg = (
                                "Sufficient",
                                "text-emerald-700",
                                "bg-emerald-100",
                            )
                        elif target_cov_v >= 15:
                            q_text, q_cls, q_bg = (
                                "Moderate",
                                "text-amber-700",
                                "bg-amber-100",
                            )
                        else:
                            q_text, q_cls, q_bg = (
                                "Low",
                                "text-rose-700",
                                "bg-rose-100",
                            )
                        cov_quality_label.set_text(f"Quality: {q_text}")
                        try:
                            cov_quality_label.classes(
                                replace=f"text-sm font-medium {q_cls}"
                            )
                            cov_quality_badge.classes(
                                replace=f"px-2 py-1 rounded {q_bg} {q_cls}"
                            )
                        except Exception:
                            pass
                        cov_quality_badge.set_text(f"{target_cov_v:.2f}x")
                    if global_cov is not None:
                        cov_global_lbl.set_text(
                            f"Global Estimated Coverage: {global_cov:.2f}x"
                        )
                    if target_cov_v is not None:
                        cov_target_lbl.set_text(
                            f"Targets Estimated Coverage: {target_cov_v:.2f}x"
                        )
                    if enrich_v is not None:
                        cov_enrich_lbl.set_text(f"Estimated enrichment: {enrich_v:.2f}x")
                except Exception as e:
                    _log_notify(
                        f"Coverage summary update failed: {e}",
                        level="warning",
                        notify=False,
                    )
            # Update state with current modification times
            if needs_update:
                state["cov_main_mtime"] = cov_main_mtime
                state["bed_cov_mtime"] = bed_cov_mtime
                state["target_cov_mtime"] = target_cov_mtime
                state["cov_time_mtime"] = cov_time_mtime
                state["last_visit_time"] = state.get("last_visit_time", time.time())
            
            launcher._coverage_state[key] = state

            # SNP Analysis results check - only update if necessary
            try:
                if (
                    state.get("last_snp_check", 0) < time.time() - 30
                ):  # Check every 30 seconds
                    # Check for SNP analysis results
                    clair_dir = sample_dir / "clair3"
                    snp_vcf = clair_dir / "snpsift_output.vcf"
                    indel_vcf = clair_dir / "snpsift_indel_output.vcf"

                    if snp_vcf.exists() and indel_vcf.exists():
                        # Update SNP results status if it exists
                        try:
                            # Find the SNP results status label and update it
                            snp_results_elements = document.querySelectorAll(
                                "[data-snp-results-status]"
                            )
                            if snp_results_elements.length > 0:
                                for element in snp_results_elements:
                                    element.textContent = (
                                        "SNP analysis completed successfully!"
                                    )
                                    element.className = "text-sm text-green-600"
                        except Exception:
                            pass

                        # Update button state if it exists
                        try:
                            snp_button_elements = document.querySelectorAll(
                                "[data-snp-analysis-button]"
                            )
                            if snp_button_elements.length > 0:
                                for element in snp_button_elements:
                                    element.textContent = "Rerun SNP Analysis"
                                    element.disabled = false
                                    element.className = "q-btn q-btn--standard q-btn--rectangle q-btn--secondary"
                        except Exception:
                            pass

                    # Also check file availability for SNP analysis
                    try:
                        target_bam = sample_dir / "target.bam"
                        targets_bed = sample_dir / "targets_exceeding_threshold.bed"

                        if target_bam.exists() and targets_bed.exists():
                            # Check if BED file has content
                            try:
                                with open(targets_bed, "r") as f:
                                    bed_content = f.read().strip()

                                if bed_content:
                                    # Files are available, enable button if not already enabled
                                    try:
                                        snp_button_elements = document.querySelectorAll(
                                            "[data-snp-analysis-button]"
                                        )
                                        if snp_button_elements.length > 0:
                                            for element in snp_button_elements:
                                                if element.disabled:
                                                    element.disabled = false
                                                    element.className = "q-btn q-btn--standard q-btn--rectangle q-btn--primary"
                                    except Exception:
                                        pass
                            except Exception:
                                pass
                    except Exception:
                        pass

                    state["last_snp_check"] = time.time()
            except Exception:
                # Never let SNP logic break coverage refresh
                pass

            # Lightweight Variant Analysis results check - only update if necessary
            try:
                if (
                    state.get("last_lga_check", 0) < time.time() - 30
                ):  # Check every 30 seconds
                    # Also check file availability for lightweight variant analysis
                    try:
                        target_bam = sample_dir / "target.bam"
                        targets_bed = sample_dir / "targets_exceeding_threshold.bed"

                    
                    except Exception as e:
                        logging.debug(f"   LGA: <access denied>: {e}")
                        pass

                    state["last_lga_check"] = time.time()
            except Exception:
                # Never let LGA logic break coverage refresh
                pass

            # IGV management - only update status if already ready, don't auto-load
            try:
                if _is_igv_ready():
                    # IGV is ready, just update status occasionally
                    if (
                        state.get("last_igv_status_update", 0) < time.time() - 60
                    ):  # Update status every minute
                        try:
                            igv_status.set_text("IGV is ready and loaded.")
                            state["last_igv_status_update"] = time.time()
                        except Exception:
                            pass
            except Exception:
                # Never let IGV logic break coverage refresh
                pass

            # After file/UI updates, re-apply chart chrome so light mode stays correct on first load
            # (theme timers can run before series/options are set; refresh runs at 0.5s).
            try:
                _sync_coverage_echarts_theme(True)
            except Exception:
                pass
        except Exception as e:
            _log_notify(
                f"Unexpected coverage refresh error: {e}", level="error", notify=True
            )

    _last_cov_theme_sig: list = [None]

    def _sync_coverage_echarts_theme(force: bool = False) -> None:
        """Re-apply coverage chart chrome when the user toggles light/dark mode.

        On first paint, ``app.storage.user['dark_mode']`` and ``body--dark`` can lag
        behind chart options built at layout time; ``force=True`` reapplies the current
        palette so plots match the visible theme without toggling.
        """
        try:
            cur = bool(app.storage.user.get("dark_mode"))
        except Exception:
            cur = False
        if not force and _last_cov_theme_sig[0] == cur:
            return
        _last_cov_theme_sig[0] = cur
        try:
            _apply_coverage_target_cov_chrome(echart_target_cov)
            echart_target_cov.update()
        except Exception:
            pass
        try:
            _apply_coverage_time_series_chrome(echart_time)
            echart_time.update()
        except Exception:
            pass
        try:
            _apply_coverage_boxplot_chrome(target_boxplot)
            target_boxplot.update()
        except Exception:
            pass
        try:
            tc = target_coverage_time_chart_state.get("chart")
            if tc is not None:
                _apply_target_coverage_time_analysis_chrome(tc)
                tc.update()
        except Exception:
            pass
        try:
            g = launcher._coverage_state.get(str(sample_dir), {}).get(
                "lga_gene_coverage_chart"
            )
            if g is not None:
                _apply_lga_gene_bar_chrome(g)
                g.update()
        except Exception:
            pass

    # Trigger an immediate refresh on page load, then continue with periodic refreshes
    refresh_timer = ui.timer(
        30.0, _refresh_coverage_async, active=True, immediate=False
    )  # Periodic refresh every 30 seconds
    ui.timer(0.5, _refresh_coverage_async, once=True)
    # Theme sync is centralized via a single global timer in theme.py.
    unregister_cov_theme_sync = register_theme_sync_callback(
        _sync_coverage_echarts_theme,
        element=echart_target_cov,
        interval_s=1.0,
        immediate=True,
    )
    try:
        ui.context.client.on_disconnect(
            lambda: (refresh_timer.deactivate(), unregister_cov_theme_sync())
        )
    except Exception:
        pass
