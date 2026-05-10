from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional
import asyncio
import threading
import json
import csv
import time
from datetime import datetime
import hashlib
import logging
from collections import OrderedDict

try:
    from nicegui import ui
except ImportError:  # pragma: no cover
    ui = None

from robin.classification_config import get_confidence_ui_tier
from robin.gui.config import (
    get_confidence_level as get_classifier_confidence_level,
    is_section_enabled,
    get_enabled_classification_steps,
    CLASSIFICATION_STEPS,
)


_SUMMARY_CACHE: "OrderedDict[str, Dict[str, Any]]" = OrderedDict()
_SUMMARY_CACHE_LOCK = threading.Lock()
_SUMMARY_CACHE_MAX_SAMPLES = 64
_SUMMARY_CACHE_MAX_AGE_SECONDS = 15 * 60


def _evict_summary_cache_locked(now_ts: Optional[float] = None) -> None:
    now_ts = now_ts or time.time()
    stale_keys = [
        key
        for key, payload in _SUMMARY_CACHE.items()
        if (now_ts - float(payload.get("_cached_at", 0.0))) > _SUMMARY_CACHE_MAX_AGE_SECONDS
    ]
    for key in stale_keys:
        _SUMMARY_CACHE.pop(key, None)
    while len(_SUMMARY_CACHE) > _SUMMARY_CACHE_MAX_SAMPLES:
        _SUMMARY_CACHE.popitem(last=False)


def _get_summary_cache(sample_dir: Path) -> Dict[str, Any]:
    key = str(sample_dir)
    with _SUMMARY_CACHE_LOCK:
        payload = _SUMMARY_CACHE.get(key)
        if not payload:
            return {}
        now_ts = time.time()
        if (now_ts - float(payload.get("_cached_at", 0.0))) > _SUMMARY_CACHE_MAX_AGE_SECONDS:
            _SUMMARY_CACHE.pop(key, None)
            return {}
        _SUMMARY_CACHE.move_to_end(key)
        return {k: v for k, v in payload.items() if k != "_cached_at"}


def _set_summary_cache(sample_dir: Path, data: Dict[str, Any]) -> None:
    key = str(sample_dir)
    with _SUMMARY_CACHE_LOCK:
        payload = dict(data)
        payload["_cached_at"] = time.time()
        _SUMMARY_CACHE[key] = payload
        _SUMMARY_CACHE.move_to_end(key)
        _evict_summary_cache_locked()


def _clear_summary_cache(sample_dir: Path) -> None:
    key = str(sample_dir)
    with _SUMMARY_CACHE_LOCK:
        _SUMMARY_CACHE.pop(key, None)


def _run_summary_cell(
    icon: str,
    label: str,
    value: str,
    *,
    col_class: str = "",
    hint: str = "",
) -> None:
    """One Run Summary grid cell; colors from styles.css (§7, .body--dark — not Tailwind dark:)."""
    cell = ("run-summary-cell min-w-0 " + col_class).strip()
    with ui.element("div").classes(cell):
        with ui.element("div").classes("run-summary-cell-inner"):
            ui.icon(icon).classes("run-summary-icon")
            with ui.element("div").classes("run-summary-stack"):
                ui.label(label).classes("run-summary-label")
                vl = ui.label(value).classes("run-summary-value")
                if hint:
                    safe = hint.replace('"', "'")
                    vl.props(f'title="{safe}"')


@ui.refreshable
def _run_info_section(sample_dir: Path, sample_id: str):
    """Run Summary block (design.md §7): responsive grid, light/dark tokens."""
    cache = _get_summary_cache(sample_dir)
    run_info = cache.get(
        "run_info",
        {
            "run_time": "Loading...",
            "model": "Loading...",
            "device": "Loading...",
            "flow_cell": "Loading...",
            "panel": "Loading...",
            "bam_passed": "Not available",
            "bam_failed": "Not available",
            "total_bases": "Not available",
            "mapped_reads": "Not available",
            "unmapped_reads": "Not available",
            "bam_batches": "Not available",
        },
    )

    rt = run_info.get("run_time", "Not available")
    model = run_info.get("model", "Not available")
    device = run_info.get("device", "Not available")
    flow = run_info.get("flow_cell", "Not available")
    panel = run_info.get("panel", "Not available")

    # (icon, label, value, col_span); optional tooltip text via hints map
    hints = {
        "Run time": "run_info_run_time (from master.csv; shown in local time with seconds)",
        "Device": "run_info_device",
        "Flow cell": "run_info_flow_cell",
        "Analysis panel": "analysis_panel",
        "Basecall model": "run_info_model",
        "Sample ID": "Unique sample identifier",
        "BAM passed": "counter_bam_passed",
        "BAM failed": "counter_bam_failed",
        "Total bases": "counter_bases_count",
        "Mapped reads": "counter_mapped_count",
        "Unmapped reads": "counter_unmapped_count",
        "BAM batches": "bam_tracking_counter / bam_tracking_total_files",
    }
    cells: List[tuple] = [
        ("schedule", "Run time", rt, ""),
        ("smartphone", "Device", device, ""),
        ("tag", "Flow cell", flow, ""),
        ("science", "Analysis panel", panel, ""),
        ("settings", "Basecall model", model, "run-summary-cell--basecall"),
        ("pin", "Sample ID", sample_id, "run-summary-cell--sample"),
    ]
    for icn, lab, key in (
        ("check_circle", "BAM passed", "bam_passed"),
        ("cancel", "BAM failed", "bam_failed"),
        ("analytics", "Total bases", "total_bases"),
        ("map", "Mapped reads", "mapped_reads"),
        ("link_off", "Unmapped reads", "unmapped_reads"),
        ("layers", "BAM batches", "bam_batches"),
    ):
        val = run_info.get(key, "Not available")
        if val != "Not available":
            cells.append((icn, lab, val, ""))

    # §7.2 / §7.3 — shell is flex column + full-width grid (see styles.css; avoid ui.column shrink)
    with ui.element("div").classes("run-summary-shell w-full min-w-0"):
        ui.label("Run summary").classes("run-summary-title text-headline-small")
        with ui.element("div").classes("run-summary-grid"):
            for icn, lab, val, span in cells:
                ht = hints.get(lab, "")
                _run_summary_cell(
                    icn, lab, val, col_class=span, hint=ht
                )


@ui.refreshable
def _classification_section(sample_dir: Path, launcher: Any = None):
    """Create refreshable classification section."""
    cache = _get_summary_cache(sample_dir)
    classification_data = cache.get(
        "classification_data",
        {
            "sturgeon": {
                "classification": "Loading...",
                "confidence": 0.0,
                "confidence_level": "Loading...",
                "features": 0,
            },
            "nanodx": {
                "classification": "Loading...",
                "confidence": 0.0,
                "confidence_level": "Loading...",
                "features": 0,
            },
            "pannanodx": {
                "classification": "Loading...",
                "confidence": 0.0,
                "confidence_level": "Loading...",
                "features": 0,
            },
            "random_forest": {
                "classification": "Loading...",
                "confidence": 0.0,
                "confidence_level": "Loading...",
                "features": 0,
            },
        },
    )
    
    # Get workflow steps from launcher if available
    workflow_steps = launcher.workflow_steps if launcher and hasattr(launcher, 'workflow_steps') else None
    enabled_classification_steps = get_enabled_classification_steps(workflow_steps)
    
    # Only show section if at least one classification step is enabled
    if workflow_steps and not enabled_classification_steps:
        return
    
    with ui.element("div").classes("classification-insight-shell w-full min-w-0"):
        ui.label("Classification details").classes(
            "classification-insight-heading text-headline-small"
        )
        with ui.element("div").classes("classification-insight-grid"):
            # Sturgeon
            if not workflow_steps or "sturgeon" in enabled_classification_steps:
                sturgeon_data = classification_data.get("sturgeon", {})
                _create_classification_dashboard_card_with_data(
                    "Sturgeon",
                    sturgeon_data.get("classification", "Not available"),
                    sturgeon_data.get("confidence", 0.0),
                    sturgeon_data.get("confidence_level", "Not available"),
                    sturgeon_data.get("features", 0),
                    "psychology",
                    "CNS tumor classification",
                    "sturgeon",
                )

            # NanoDX
            if not workflow_steps or "nanodx" in enabled_classification_steps:
                nanodx_data = classification_data.get("nanodx", {})
                _create_classification_dashboard_card_with_data(
                    "NanoDX",
                    nanodx_data.get("classification", "Not available"),
                    nanodx_data.get("confidence", 0.0),
                    nanodx_data.get("confidence_level", "Not available"),
                    nanodx_data.get("features", 0),
                    "biotech",
                    "CNS tumor classification",
                    "nanodx",
                )

            # PanNanoDX
            if not workflow_steps or "pannanodx" in enabled_classification_steps:
                pannanodx_data = classification_data.get("pannanodx", {})
                _create_classification_dashboard_card_with_data(
                    "PanNanoDX",
                    pannanodx_data.get("classification", "Not available"),
                    pannanodx_data.get("confidence", 0.0),
                    pannanodx_data.get("confidence_level", "Not available"),
                    pannanodx_data.get("features", 0),
                    "science",
                    "Pan-cancer tumor classification",
                    "pannanodx",
                )

            # Random Forest
            if not workflow_steps or "random_forest" in enabled_classification_steps:
                rf_data = classification_data.get("random_forest", {})
                _create_classification_dashboard_card_with_data(
                    "Random Forest",
                    rf_data.get("classification", "Not available"),
                    rf_data.get("confidence", 0.0),
                    rf_data.get("confidence_level", "Not available"),
                    rf_data.get("features", 0),
                    "forest",
                    "CNS Machine learning classification",
                    "random_forest",
                )


@ui.refreshable
def _analysis_section(sample_dir: Path, launcher: Any = None):
    """Create refreshable analysis section."""
    # Get workflow steps from launcher if available
    workflow_steps = launcher.workflow_steps if launcher and hasattr(launcher, 'workflow_steps') else None
    cache = _get_summary_cache(sample_dir)
    analysis_data = cache.get("analysis_data", {})
    coverage_data = analysis_data.get("coverage", {}) if not workflow_steps or is_section_enabled("target", workflow_steps) else {}
    cnv_data = analysis_data.get("cnv", {}) if not workflow_steps or is_section_enabled("cnv", workflow_steps) else {}
    mgmt_data = analysis_data.get("mgmt", {}) if not workflow_steps or is_section_enabled("mgmt", workflow_steps) else {}
    fusion_data = analysis_data.get("fusion", {}) if not workflow_steps or is_section_enabled("fusion", workflow_steps) else {}
    
    # Check if any analysis section should be shown
    should_show_target = not workflow_steps or is_section_enabled("target", workflow_steps)
    should_show_cnv = not workflow_steps or is_section_enabled("cnv", workflow_steps)
    should_show_mgmt = not workflow_steps or is_section_enabled("mgmt", workflow_steps)
    should_show_fusion = not workflow_steps or is_section_enabled("fusion", workflow_steps)
    
    if workflow_steps and not any([should_show_target, should_show_cnv, should_show_mgmt, should_show_fusion]):
        return

    with ui.element("div").classes("classification-insight-shell w-full min-w-0"):
        ui.label("Analysis details").classes(
            "classification-insight-heading text-headline-small"
        )
        with ui.element("div").classes("classification-insight-grid"):
            # Coverage Analysis (target)
            if should_show_target:
                _create_coverage_dashboard_card_with_data(
                    coverage_data.get("quality", "Not available"),
                    coverage_data.get("target_coverage", "Not available"),
                    coverage_data.get("global_coverage", "Not available"),
                    coverage_data.get("enrichment", "Not available"),
                    anchor_key="coverage",
                )

            # CNV Analysis
            if should_show_cnv:
                _create_cnv_dashboard_card_with_data(
                    cnv_data.get("genetic_sex", "Not available"),
                    cnv_data.get("bin_width", "Not available"),
                    cnv_data.get("variance", "Not available"),
                    cnv_data.get("gained", 0),
                    cnv_data.get("lost", 0),
                    anchor_key="cnv",
                )

            # MGMT Analysis
            if should_show_mgmt:
                _create_mgmt_dashboard_card_with_data(
                    mgmt_data.get("status", "Not available"),
                    mgmt_data.get("methylation_percent", "Not available"),
                    mgmt_data.get("average_methylation", "Not available"),
                    mgmt_data.get("prediction_score", "Not available"),
                    mgmt_data.get("cpg_sites", "Not available"),
                    anchor_key="mgmt",
                )

            # Fusion Analysis
            if should_show_fusion:
                panel = _get_analysis_panel(sample_dir)
                _create_fusion_dashboard_card_with_data(
                    panel,
                    fusion_data.get("target_fusions", 0),
                    fusion_data.get("genome_fusions", 0),
                    fusion_data.get("target_pairs", 0),
                    fusion_data.get("target_groups", 0),
                    fusion_data.get("genome_pairs", 0),
                    fusion_data.get("genome_groups", 0),
                    anchor_key="fusion",
                )


def add_summary_section(sample_dir: Path, sample_id: str, launcher: Any = None) -> None:
    """Build the Summary section at the top of the sample detail page using ui.refreshable."""

    # Create refreshable sections synchronously for immediate rendering
    _run_info_section(sample_dir, sample_id)
    _classification_section(sample_dir, launcher)
    _analysis_section(sample_dir, launcher)

    async def _refresh_summary_cache_async() -> None:
        try:
            data = await asyncio.to_thread(
                _refresh_summary_cache_sync, sample_dir, sample_id, launcher
            )
            _set_summary_cache(sample_dir, data)
            _run_info_section.refresh()
            _classification_section.refresh()
            _analysis_section.refresh()
        except Exception as e:
            logging.debug(f"[Summary] Async refresh failed: {e}")

    def _schedule_refresh() -> None:
        try:
            asyncio.create_task(_refresh_summary_cache_async())
        except RuntimeError:
            # If no running loop, fall back to sync refresh (still updates cache)
            data = _refresh_summary_cache_sync(sample_dir, sample_id, launcher)
            _set_summary_cache(sample_dir, data)
            _run_info_section.refresh()
            _classification_section.refresh()
            _analysis_section.refresh()

    # Initial async load shortly after page render
    ui.timer(0.2, _refresh_summary_cache_async, once=True)
    # Periodic refresh timer (every 30 seconds)
    refresh_timer = ui.timer(
        30.0, _refresh_summary_cache_async, active=True, immediate=False
    )
    try:
        def _on_disconnect_cleanup() -> None:
            refresh_timer.deactivate()
            _clear_summary_cache(sample_dir)

        ui.context.client.on_disconnect(_on_disconnect_cleanup)
    except Exception:
        pass


def _refresh_summary_cache_sync(
    sample_dir: Path, sample_id: str, launcher: Any = None
) -> Dict[str, Any]:
    """Compute summary data in a background thread."""
    workflow_steps = launcher.workflow_steps if launcher and hasattr(launcher, "workflow_steps") else None
    data: Dict[str, Any] = {
        "run_info": _extract_run_information(sample_dir, sample_id),
        "classification_data": _extract_classification_data(sample_dir),
        "analysis_data": {},
    }

    # Analysis data (only compute enabled sections)
    analysis_data: Dict[str, Any] = {}
    if not workflow_steps or is_section_enabled("target", workflow_steps):
        analysis_data["coverage"] = _extract_coverage_data(sample_dir)
    if not workflow_steps or is_section_enabled("cnv", workflow_steps):
        analysis_data["cnv"] = _extract_cnv_data(sample_dir)
    if not workflow_steps or is_section_enabled("mgmt", workflow_steps):
        analysis_data["mgmt"] = _extract_mgmt_data(sample_dir)
    if not workflow_steps or is_section_enabled("fusion", workflow_steps):
        analysis_data["fusion"] = _extract_fusion_data(sample_dir)

    data["analysis_data"] = analysis_data
    return data


def _scroll_to_classification_detail(anchor_key: str) -> None:
    """Scroll to the detailed classification panel for this tool (design.md §8)."""
    safe = "".join(c for c in anchor_key if c.isalnum() or c in "_-")
    if not safe:
        return
    try:
        ui.run_javascript(
            f'document.getElementById("classification-detail-{safe}")'
            f'?.scrollIntoView({{behavior:"smooth",block:"start"}})'
        )
    except Exception:
        pass


def _scroll_to_analysis_detail(anchor_key: str) -> None:
    """Scroll to the detailed analysis section (coverage, CNV, MGMT, fusion). design.md §9."""
    safe = "".join(c for c in anchor_key if c.isalnum() or c in "_-")
    if not safe:
        return
    try:
        ui.run_javascript(
            f'document.getElementById("analysis-detail-{safe}")'
            f'?.scrollIntoView({{behavior:"smooth",block:"start"}})'
        )
    except Exception:
        pass


def _create_classification_dashboard_card_with_data(
    title: str,
    classification: str,
    confidence: float,
    confidence_level: str,
    features: int,
    icon: str,
    description: str,
    anchor_key: str,
) -> None:
    """Classification Insight card — design.md §8; click scrolls to detail expansion."""
    pct = max(0.0, min(100.0, float(confidence)))
    tier = get_confidence_ui_tier(anchor_key, pct)

    def _on_card_click() -> None:
        _scroll_to_classification_detail(anchor_key)

    with (
        ui.element("div")
        .classes(
            "classification-insight-card w-full min-w-0 "
            "cursor-pointer transition-all duration-200 select-none"
        )
        .on("click", _on_card_click)
    ):
        with ui.column().classes("w-full min-w-0 gap-2 p-2 md:p-3"):
            with ui.row().classes("items-center gap-2 min-w-0"):
                ui.icon(icon).classes("classification-insight-icon")
                ui.label(title).classes("classification-insight-model flex-1 min-w-0")
            _cls_title = classification.replace('"', "'")
            ui.label(classification).classes(
                "classification-insight-result w-full"
            ).props(f'title="{_cls_title}"')
            with ui.element("div").classes("classification-insight-track w-full"):
                ui.element("div").classes(
                    f"classification-insight-bar classification-insight-bar--{tier}"
                ).style(f"width: {pct}%")
            with ui.row().classes(
                "w-full justify-between items-baseline gap-2 flex-wrap"
            ):
                ui.label(f"Confidence: {confidence:.1f}%").classes(
                    "classification-insight-meta"
                )
                ui.label(f"Features: {features:,}").classes(
                    "classification-insight-meta"
                )
            ui.label(confidence_level).classes(
                f"classification-insight-level classification-insight-level--{tier}"
            )
            ui.label(description).classes("classification-insight-foot")


def _create_classification_dashboard_card(title: str, classification: str, icon: str, description: str) -> Dict[str, Any]:
    """Create a compact classification dashboard card with detailed information. Returns labels for updating."""
    with ui.card().classes("flex-1 bg-white rounded-lg shadow-sm border border-gray-200 p-4"):
        with ui.row().classes("flex items-center justify-between mb-2"):
            # Title
            ui.label(title).classes("text-sm font-medium text-gray-600")
            
            # Icon in circular background
            with ui.row().classes("w-7 h-7 bg-blue-100 rounded-full flex items-center justify-center"):
                ui.icon(icon).classes("w-3.5 h-3.5 text-blue-600")
        
        # Main classification result with confidence badge
        with ui.row().classes("flex items-center justify-between mb-1"):
            classification_label = ui.label(classification).classes("text-xl font-bold text-gray-900")
            # Confidence level badge
            confidence_badge = ui.label("Loading...").classes("px-2 py-1 text-xs font-medium rounded-full bg-gray-100 text-gray-600")
        
        # Compact details in a single row
        with ui.row().classes("flex items-center justify-between mb-1"):
            confidence_label = ui.label("Confidence: Loading...").classes("text-xs font-medium text-gray-700")
            features_label = ui.label("Features: Loading...").classes("text-xs text-gray-500")
        
        # Description
        ui.label(description).classes("text-xs text-gray-500")
        
        return {
            "classification": classification_label,
            "confidence": confidence_label,
            "confidence_badge": confidence_badge,
            "features": features_label
        }


def _create_classification_card(
    title: str,
    classification: str,
    confidence: str,
    confidence_level: str,
    features: str = "",
    features_label: str = "",
    model: str = "",
) -> Dict[str, Any]:
    """Create a classification summary card. Returns labels for updating."""
    labels = {}
    with ui.card().classes("flex-1 elevation-4 rounded-xl bg-gradient-to-br from-blue-50 to-indigo-50 border-l-4 border-blue-500"):
        ui.label(f"{title}").classes("font-bold text-blue-800 mb-2")
        ui.separator().classes().style("border: 1px solid var(--md-primary)")
        labels["classification"] = ui.label(f"Class: {classification}").classes("font-bold text-medium text-blue-600")
        confidence_color = _get_confidence_color(confidence_level)
        labels["confidence"] = ui.label(f"Confidence: {confidence}%").classes(f"text-sm text-{confidence_color}-600")
        labels["confidence_level"] = ui.label(confidence_level).classes(f"text-sm text-{confidence_color}-600")
        if model:
            ui.label(f"Model: {model}").classes("text-sm text-blue-600")
        if features and features_label:
            labels["features"] = ui.label(f"Features: {features_label}").classes("text-sm text-blue-600")

    return labels


def _parse_target_coverage_x(target_coverage: str) -> float | None:
    try:
        return float(str(target_coverage).replace("x", "").strip())
    except (ValueError, AttributeError, TypeError):
        return None


def _coverage_depth_tier(coverage_num: float) -> str:
    if coverage_num >= 30:
        return "high"
    if coverage_num >= 20:
        return "medium"
    return "low"


def _create_coverage_dashboard_card_with_data(
    quality: str,
    target_coverage: str,
    global_coverage: str,
    enrichment: str,
    anchor_key: str = "coverage",
) -> None:
    """Coverage insight card — same shell as classification (design.md §9)."""
    cov_num = _parse_target_coverage_x(target_coverage)
    if cov_num is not None:
        tier = _coverage_depth_tier(cov_num)
    else:
        tier = "low"

    def _on_click() -> None:
        _scroll_to_analysis_detail(anchor_key)

    _q = quality.replace('"', "'")
    with (
        ui.element("div")
        .classes(
            "classification-insight-card w-full min-w-0 "
            "cursor-pointer transition-all duration-200 select-none"
        )
        .on("click", _on_click)
    ):
        with ui.column().classes("w-full min-w-0 gap-2 p-2 md:p-3"):
            with ui.row().classes("items-center gap-2 min-w-0"):
                ui.icon("analytics").classes("classification-insight-icon")
                ui.label("Coverage analysis").classes(
                    "classification-insight-model flex-1 min-w-0"
                )
            ui.label(quality).classes("classification-insight-result w-full").props(
                f'title="{_q}"'
            )
            with ui.row().classes(
                "w-full justify-between items-baseline gap-2 flex-wrap"
            ):
                ui.label(f"Target: {target_coverage}").classes(
                    f"classification-insight-level classification-insight-level--{tier}"
                )
            with ui.column().classes("w-full gap-1"):
                ui.label(f"Global: {global_coverage}").classes(
                    "classification-insight-meta"
                )
                ui.label(f"Targets: {target_coverage}").classes(
                    "classification-insight-meta"
                )
                ui.label(f"Enrichment: {enrichment}").classes(
                    "classification-insight-meta"
                )
            with ui.row().classes("gap-1 flex-wrap"):
                ui.label("≥30x").classes(
                    "analysis-insight-pill analysis-insight-pill--emerald"
                )
                ui.label("≥20x").classes("analysis-insight-pill analysis-insight-pill--sky")
                ui.label("≥10x").classes(
                    "analysis-insight-pill analysis-insight-pill--amber"
                )
                ui.label("<10x").classes("analysis-insight-pill analysis-insight-pill--rose")
            ui.label("Target panel coverage quality and depth").classes(
                "classification-insight-foot"
            )


def _create_coverage_dashboard_card() -> Dict[str, Any]:
    """Create a compact coverage analysis dashboard card. Returns labels for updating."""
    with ui.card().classes("flex-1 bg-white rounded-lg shadow-sm border border-gray-200 p-4"):
        with ui.row().classes("flex items-center justify-between mb-2"):
            # Title
            ui.label("Coverage Analysis").classes("text-sm font-medium text-gray-600")
            
            # Icon in circular background
            with ui.row().classes("w-7 h-7 bg-blue-100 rounded-full flex items-center justify-center"):
                ui.icon("analytics").classes("w-3.5 h-3.5 text-blue-600")
        
        # Main quality result with badge
        with ui.row().classes("flex items-center justify-between mb-1"):
            quality_label = ui.label("Loading...").classes("text-xl font-bold text-gray-900")
            # Coverage badge
            coverage_badge = ui.label("--x").classes("px-2 py-1 text-xs font-medium rounded-full bg-gray-100 text-gray-600")
        
        # Coverage details in compact layout
        with ui.column().classes("mb-2"):
            global_coverage_label = ui.label("Global: Loading...").classes("text-xs text-gray-600")
            target_coverage_label = ui.label("Targets: Loading...").classes("text-xs text-gray-600")
            enrichment_label = ui.label("Enrichment: Loading...").classes("text-xs text-gray-600")
        
        # Coverage thresholds as small badges
        with ui.row().classes("gap-1 flex-wrap"):
            ui.label("≥30x").classes("px-1 py-0.5 text-xs bg-green-100 text-green-800 rounded")
            ui.label("≥20x").classes("px-1 py-0.5 text-xs bg-blue-100 text-blue-800 rounded")
            ui.label("≥10x").classes("px-1 py-0.5 text-xs bg-yellow-100 text-yellow-800 rounded")
            ui.label("<10x").classes("px-1 py-0.5 text-xs bg-red-100 text-red-800 rounded")
        
        return {
            "quality": quality_label,
            "coverage_badge": coverage_badge,
            "global_coverage": global_coverage_label,
            "target_coverage": target_coverage_label,
            "enrichment": enrichment_label
        }


def _create_cnv_dashboard_card_with_data(
    genetic_sex: str,
    bin_width: str,
    variance: str,
    gained: int,
    lost: int,
    anchor_key: str = "cnv",
) -> None:
    """CNV insight card — design.md §9."""
    total = int(gained) + int(lost)
    total_denom = str(total) if total > 0 else "—"

    def _on_click() -> None:
        _scroll_to_analysis_detail(anchor_key)

    _gs = genetic_sex.replace('"', "'")
    with (
        ui.element("div")
        .classes(
            "classification-insight-card w-full min-w-0 "
            "cursor-pointer transition-all duration-200 select-none"
        )
        .on("click", _on_click)
    ):
        with ui.column().classes("w-full min-w-0 gap-2 p-2 md:p-3"):
            with ui.row().classes("items-center gap-2 min-w-0"):
                ui.icon("person").classes("classification-insight-icon")
                ui.label("Copy number (CNV)").classes(
                    "classification-insight-model flex-1 min-w-0"
                )
            ui.label(genetic_sex).classes("classification-insight-result w-full").props(
                f'title="{_gs}"'
            )
            ui.label(f"Gained vs total: {gained} / {total_denom}").classes(
                "classification-insight-meta w-full"
            )
            with ui.column().classes("w-full gap-1"):
                ui.label(f"Bin width: {bin_width}").classes("classification-insight-meta")
                ui.label(f"Variance: {variance}").classes("classification-insight-meta")
            with ui.row().classes("gap-1 flex-wrap"):
                ui.label(f"Gained: {gained}").classes(
                    "analysis-insight-pill analysis-insight-pill--emerald"
                )
                ui.label(f"Lost: {lost}").classes(
                    "analysis-insight-pill analysis-insight-pill--rose"
                )
            ui.label(
                "Copy number across the genome with breakpoint detection"
            ).classes("classification-insight-foot")


def _create_cnv_dashboard_card() -> Dict[str, Any]:
    """Create a compact CNV analysis dashboard card. Returns labels for updating."""
    with ui.card().classes("flex-1 bg-white rounded-lg shadow-sm border border-gray-200 p-4"):
        with ui.row().classes("flex items-center justify-between mb-2"):
            # Title
            ui.label("Copy Number Analysis").classes("text-sm font-medium text-gray-600")
            
            # Icon in circular background
            with ui.row().classes("w-7 h-7 bg-purple-100 rounded-full flex items-center justify-center"):
                ui.icon("person").classes("w-3.5 h-3.5 text-purple-600")
        
        # Main genetic sex result
        genetic_sex_label = ui.label("Loading...").classes("text-xl font-bold text-gray-900 mb-1")
        
        # Analysis details in compact layout
        with ui.column().classes("mb-2"):
            bin_width_label = ui.label("Bin Width: Loading...").classes("text-xs text-gray-600")
            variance_label = ui.label("Variance: Loading...").classes("text-xs text-gray-600")
        
        # CNV counts as badges
        with ui.row().classes("gap-2 mb-1"):
            gained_badge = ui.label("Gained: --").classes("px-2 py-1 text-xs font-medium rounded-full bg-green-100 text-green-800")
            lost_badge = ui.label("Lost: --").classes("px-2 py-1 text-xs font-medium rounded-full bg-red-100 text-red-800")
        
        # Description
        ui.label("Copy number analysis across genome with breakpoint detection").classes("text-xs text-gray-500")
        
        return {
            "genetic_sex": genetic_sex_label,
            "bin_width": bin_width_label,
            "variance": variance_label,
            "gained": gained_badge,
            "lost": lost_badge
        }


def _create_mgmt_dashboard_card_with_data(
    status: str,
    methylation_percent: str,
    average_methylation: str,
    prediction_score: str,
    cpg_sites: str,
    anchor_key: str = "mgmt",
) -> None:
    """MGMT insight card — design.md §9."""
    try:
        meth_num = float(str(methylation_percent).replace("%", "").strip())
    except (ValueError, AttributeError, TypeError):
        meth_num = None

    if meth_num is not None:
        if meth_num > 10:
            tier = "high"
        elif meth_num > 5:
            tier = "medium"
        else:
            tier = "low"
    else:
        tier = "low"

    def _on_click() -> None:
        _scroll_to_analysis_detail(anchor_key)

    _st = status.replace('"', "'")
    with (
        ui.element("div")
        .classes(
            "classification-insight-card w-full min-w-0 "
            "cursor-pointer transition-all duration-200 select-none"
        )
        .on("click", _on_click)
    ):
        with ui.column().classes("w-full min-w-0 gap-2 p-2 md:p-3"):
            with ui.row().classes("items-center gap-2 min-w-0"):
                ui.icon("science").classes("classification-insight-icon")
                ui.label("MGMT methylation").classes(
                    "classification-insight-model flex-1 min-w-0"
                )
            ui.label(status).classes("classification-insight-result w-full").props(
                f'title="{_st}"'
            )
            with ui.row().classes(
                "w-full justify-between items-baseline gap-2 flex-wrap"
            ):
                ui.label(f"Methylation: {methylation_percent}").classes(
                    f"classification-insight-level classification-insight-level--{tier}"
                )
                ui.label("Promoter %").classes("classification-insight-meta")
            with ui.column().classes("w-full gap-1"):
                ui.label(f"Average: {average_methylation}").classes(
                    "classification-insight-meta"
                )
                ui.label(f"Score: {prediction_score}").classes(
                    "classification-insight-meta"
                )
            ui.label(
                f"Status from methylation at {cpg_sites} CpG sites"
            ).classes("classification-insight-foot")


def _create_mgmt_dashboard_card() -> Dict[str, Any]:
    """Create a compact MGMT analysis dashboard card. Returns labels for updating."""
    with ui.card().classes("flex-1 bg-white rounded-lg shadow-sm border border-gray-200 p-4"):
        with ui.row().classes("flex items-center justify-between mb-2"):
            # Title
            ui.label("MGMT Analysis").classes("text-sm font-medium text-gray-600")
            
            # Icon in circular background
            with ui.row().classes("w-7 h-7 bg-orange-100 rounded-full flex items-center justify-center"):
                ui.icon("science").classes("w-3.5 h-3.5 text-orange-600")
        
        # Main status result with badge
        with ui.row().classes("flex items-center justify-between mb-1"):
            status_label = ui.label("Loading...").classes("text-xl font-bold text-gray-900")
            # Methylation badge
            methylation_badge = ui.label("--%").classes("px-2 py-1 text-xs font-medium rounded-full bg-gray-100 text-gray-600")
        
        # Analysis details in compact layout
        with ui.column().classes("mb-2"):
            average_methylation_label = ui.label("Average: Loading...").classes("text-xs text-gray-600")
            prediction_score_label = ui.label("Score: Loading...").classes("text-xs text-gray-600")
        
        # Description
        cpg_sites_label = ui.label("MGMT status determined from methylation analysis of -- CpG sites").classes("text-xs text-gray-500")
        
        return {
            "status": status_label,
            "methylation_badge": methylation_badge,
            "average_methylation": average_methylation_label,
            "prediction_score": prediction_score_label,
            "cpg_sites": cpg_sites_label
        }


def _create_fusion_dashboard_card_with_data(
    panel: str,
    target_fusions: int,
    genome_fusions: int,
    target_pairs: int = 0,
    target_groups: int = 0,
    genome_pairs: int = 0,
    genome_groups: int = 0,
    anchor_key: str = "fusion",
) -> None:
    """Fusion insight card — design.md §9 (no single % metric; no bar)."""

    def _on_click() -> None:
        _scroll_to_analysis_detail(anchor_key)

    if target_pairs > 0 or target_groups > 0:
        target_badge = f"{target_pairs} pairs, {target_groups} groups"
    else:
        target_badge = f"{target_fusions} target fusions"

    if genome_pairs > 0 or genome_groups > 0:
        genome_line = f"{genome_pairs} pairs, {genome_groups} groups"
    else:
        genome_line = f"{genome_fusions} genome-wide fusions"

    _panel = panel.replace('"', "'")
    main_line = f"Panel: {panel}"
    with (
        ui.element("div")
        .classes(
            "classification-insight-card w-full min-w-0 "
            "cursor-pointer transition-all duration-200 select-none"
        )
        .on("click", _on_click)
    ):
        with ui.column().classes("w-full min-w-0 gap-2 p-2 md:p-3"):
            with ui.row().classes("items-center gap-2 min-w-0"):
                ui.icon("merge").classes("classification-insight-icon")
                ui.label("Fusion analysis").classes(
                    "classification-insight-model flex-1 min-w-0"
                )
            ui.label(main_line).classes("classification-insight-result w-full").props(
                f'title="{_panel}"'
            )
            with ui.row().classes(
                "w-full justify-between items-start gap-2 flex-wrap"
            ):
                ui.label(target_badge).classes(
                    "analysis-insight-pill analysis-insight-pill--sky"
                )
            with ui.column().classes("w-full gap-1"):
                ui.label(genome_line).classes("classification-insight-meta")
            ui.label(
                "Candidates from reads with supplementary alignments"
            ).classes("classification-insight-foot")


def _create_fusion_dashboard_card() -> Dict[str, Any]:
    """Create a compact fusion analysis dashboard card. Returns labels for updating."""
    with ui.card().classes("flex-1 bg-white rounded-lg shadow-sm border border-gray-200 p-4"):
        with ui.row().classes("flex items-center justify-between mb-2"):
            # Title
            ui.label("Fusion Analysis").classes("text-sm font-medium text-gray-600")
            
            # Icon in circular background
            with ui.row().classes("w-7 h-7 bg-green-100 rounded-full flex items-center justify-center"):
                ui.icon("merge").classes("w-3.5 h-3.5 text-green-600")
        
        # Panel info and main result
        with ui.row().classes("flex items-center justify-between mb-1"):
            panel_label = ui.label("Panel: --").classes("text-sm font-medium text-gray-700")
            target_fusions_badge = ui.label("-- target fusions").classes("px-2 py-1 text-xs font-medium rounded-full bg-blue-100 text-blue-800")
        
        # Analysis details in compact layout
        with ui.column().classes("mb-2"):
            genome_fusions_label = ui.label("-- genome wide fusions").classes("text-xs text-gray-600")
        
        # Description
        ui.label("Fusion candidates identified from reads with supplementary alignments").classes("text-xs text-gray-500")
        
        return {
            "panel": panel_label,
            "target_fusions": target_fusions_badge,
            "genome_fusions": genome_fusions_label
        }






def _fmt_master_csv_integer(v: Optional[str]) -> Optional[str]:
    """Format integer counters from master.csv without float precision loss on large values."""
    if v is None or not str(v).strip():
        return None
    s = str(v).replace(",", "").strip()
    try:
        if "e" in s.lower():
            n = int(float(s))
        elif "." in s:
            n = int(float(s))
        else:
            n = int(s)
        return f"{n:,}"
    except (ValueError, OverflowError, TypeError):
        out = str(v).strip()
        return out or None


def _get_analysis_panel(sample_dir: Path) -> str:
    """Get the analysis panel from master.csv"""
    try:
        master_csv_path = sample_dir / "master.csv"
        if master_csv_path.exists():
            with master_csv_path.open("r", newline="") as fh:
                reader = csv.DictReader(fh)
                row = next(reader, None)
            if row:
                panel = row.get("analysis_panel", "")
                if panel and str(panel).strip() != "":
                    return str(panel).strip()
        # No fallback to rCNS2 - return empty string if not found
        return ""
    except Exception:
        # No fallback to rCNS2 - return empty string if error
        return ""


def _apply_master_csv_to_run_info(sample_dir: Path, run_info: Dict[str, str]) -> None:
    """Merge first row of master.csv into run_info.

    Maps typical columns to Run summary tiles:

    - ``run_info_run_time`` → run time (local clock, with seconds)
    - ``run_info_device`` / ``devices`` → Device
    - ``run_info_model`` / ``basecall_models`` → Basecall model
    - ``run_info_flow_cell`` / ``flowcell_ids`` → Flow cell
    - ``analysis_panel`` → Analysis panel
    - ``counter_bam_passed`` / ``counter_bam_failed`` / ``counter_bases_count`` /
      ``counter_mapped_count`` / ``counter_unmapped_count`` → BAM / bases tiles
    - ``bam_tracking_counter`` / ``bam_tracking_total_files`` → BAM batches (current / total)
    """
    master_path = sample_dir / "master.csv"
    if not master_path.exists():
        return

    def get_ci(r: Dict[str, str], key: str) -> Optional[str]:
        for k, v in r.items():
            if k and k.strip().lower() == key.lower():
                return v
        return None

    try:
        with master_path.open("r", newline="") as f:
            reader = csv.DictReader(f)
            row = next(reader, None)
        if not row:
            return

        model_val = get_ci(row, "run_info_model") or get_ci(row, "basecall_models")
        if model_val and str(model_val).strip():
            run_info["model"] = str(model_val).strip()

        device_val = get_ci(row, "run_info_device") or get_ci(row, "devices")
        if device_val and str(device_val).strip():
            run_info["device"] = str(device_val).strip()

        flowcell_val = get_ci(row, "run_info_flow_cell") or get_ci(
            row, "flowcell_ids"
        )
        if flowcell_val and str(flowcell_val).strip():
            run_info["flow_cell"] = str(flowcell_val).strip()

        panel_val = get_ci(row, "analysis_panel")
        if panel_val and str(panel_val).strip():
            run_info["panel"] = str(panel_val).strip()

        run_time_val = get_ci(row, "run_info_run_time")
        if run_time_val and str(run_time_val).strip():
            try:
                raw_ts = str(run_time_val).strip().replace("Z", "+00:00")
                dt = datetime.fromisoformat(raw_ts)
                if dt.tzinfo is not None:
                    dt = dt.astimezone()
                run_info["run_time"] = dt.strftime("%Y-%m-%d %H:%M:%S")
            except Exception:
                pass

        for csv_key, ri_key in (
            ("counter_bam_passed", "bam_passed"),
            ("counter_bam_failed", "bam_failed"),
            ("counter_bases_count", "total_bases"),
            ("counter_mapped_count", "mapped_reads"),
            ("counter_unmapped_count", "unmapped_reads"),
        ):
            raw = get_ci(row, csv_key)
            if raw is not None and str(raw).strip() != "":
                fv = _fmt_master_csv_integer(raw)
                if fv:
                    run_info[ri_key] = fv

        btc = get_ci(row, "bam_tracking_counter")
        btt = get_ci(row, "bam_tracking_total_files")
        parts = []
        if btc is not None and str(btc).strip() != "":
            parts.append(_fmt_master_csv_integer(btc) or str(btc).strip())
        if btt is not None and str(btt).strip() != "":
            parts.append(_fmt_master_csv_integer(btt) or str(btt).strip())
        if len(parts) == 2:
            run_info["bam_batches"] = f"{parts[0]} / {parts[1]}"
        elif len(parts) == 1:
            run_info["bam_batches"] = parts[0]
    except Exception:
        pass


def _extract_run_information(sample_dir: Path, sample_id: str) -> Dict[str, str]:
    """Extract run information from sample directory."""
    run_info = {
        "run_time": "Not available",
        "model": "Not available",
        "device": "Not available",
        "flow_cell": "Not available",
        "panel": "Not available",
        "bam_passed": "Not available",
        "bam_failed": "Not available",
        "total_bases": "Not available",
        "mapped_reads": "Not available",
        "unmapped_reads": "Not available",
        "bam_batches": "Not available",
    }

    try:
        # Try to get run time from file modification times
        if sample_dir.exists():
            latest_time = 0
            for file_path in sample_dir.rglob("*"):
                if file_path.is_file():
                    mtime = file_path.stat().st_mtime
                    if mtime > latest_time:
                        latest_time = mtime
            if latest_time > 0:
                run_info["run_time"] = datetime.fromtimestamp(latest_time).strftime(
                    "%Y-%m-%d %H:%M"
                )

        # Try to extract from various metadata files (search recursively)
        metadata_files = [
            "sample_metadata.json",
            "run_info.json",
            "sequencing_summary.txt",
        ]

        for metadata_file in metadata_files:
            candidate_paths: List[Path] = []
            direct_path = sample_dir / metadata_file
            if direct_path.exists():
                candidate_paths = [direct_path]
            else:
                candidate_paths = list(sample_dir.rglob(metadata_file))

            for file_path in candidate_paths:
                if metadata_file.endswith(".json"):
                    try:
                        with open(file_path, "r") as f:
                            data = json.load(f)
                            if "model" in data and data["model"]:
                                run_info["model"] = str(data["model"]).strip()
                            if "device" in data and data["device"]:
                                run_info["device"] = str(data["device"]).strip()
                            if "flow_cell" in data and data["flow_cell"]:
                                run_info["flow_cell"] = str(data["flow_cell"]).strip()
                    except Exception:
                        pass
                # Stop early if we have all fields
                if all(
                    run_info[k] != "Not available"
                    for k in ["run_time", "model", "device", "flow_cell"]
                ):
                    break
            # Stop early if complete
            if all(
                run_info[k] != "Not available"
                for k in ["run_time", "model", "device", "flow_cell"]
            ):
                break
    except Exception:
        pass

    _apply_master_csv_to_run_info(sample_dir, run_info)

    return run_info


def _extract_coverage_data(sample_dir: Path) -> Dict[str, Any]:
    """Extract coverage analysis data using the same logic as the coverage component."""
    coverage_data = {
        "quality": "Not available",
        "global_coverage": "Not available",
        "target_coverage": "Not available",
        "enrichment": "Not available",
    }

    try:
        # Use the same file reading logic as the coverage component
        cov_main = sample_dir / "coverage_main.csv"
        bed_cov = sample_dir / "bed_coverage_main.csv"

        global_cov = None
        target_cov_v = None
        enrich_v = None

        # Read coverage_main.csv for global coverage (same as coverage component)
        if cov_main.exists():
            try:
                import pandas as pd

                cov_df = pd.read_csv(
                    cov_main,
                    usecols=lambda c: c in {"covbases", "endpos"},
                    dtype={"covbases": "float64", "endpos": "float64"},
                )

                # Calculate global coverage using the same formula as coverage component
                if (
                    "covbases" in cov_df.columns
                    and "endpos" in cov_df.columns
                    and cov_df["endpos"].sum() > 0
                ):
                    global_cov = float(cov_df["covbases"].sum()) / float(
                        cov_df["endpos"].sum()
                    )
                    coverage_data["global_coverage"] = f"{global_cov:.2f}x"
            except Exception as e:
                logging.debug(f"   MGMT: <access denied>: {e}")
                pass

        # Read bed_coverage_main.csv for target coverage (same as coverage component)
        if bed_cov.exists():
            try:
                import pandas as pd

                bed_df = pd.read_csv(
                    bed_cov,
                    usecols=lambda c: c in {"bases", "length", "endpos", "startpos"},
                    dtype={
                        "bases": "float64",
                        "length": "float64",
                        "endpos": "float64",
                        "startpos": "float64",
                    },
                )

                # Calculate target coverage using the same formula as coverage component
                if "bases" in bed_df.columns:
                    # Calculate length if not present (same as coverage component)
                    if "length" not in bed_df.columns:
                        bed_df["length"] = (
                            bed_df["endpos"] - bed_df["startpos"] + 1
                        ).astype(float)

                    if bed_df["length"].sum() > 0:
                        target_cov_v = float(bed_df["bases"].sum()) / float(
                            bed_df["length"].sum()
                        )
                        coverage_data["target_coverage"] = f"{target_cov_v:.2f}x"

                        # Determine quality based on target coverage (same logic as coverage component)
                        if target_cov_v >= 30:
                            coverage_data["quality"] = "Excellent"
                        elif target_cov_v >= 20:
                            coverage_data["quality"] = "Good"
                        elif target_cov_v >= 10:
                            coverage_data["quality"] = "Moderate"
                        else:
                            coverage_data["quality"] = "Insufficient"
            except Exception as e:
                logging.debug(f"   MGMT: <access denied>: {e}")
                pass

        # Calculate enrichment if we have both values (same logic as coverage component)
        if global_cov is not None and target_cov_v is not None and global_cov > 0:
            enrich_v = target_cov_v / global_cov
            coverage_data["enrichment"] = f"{enrich_v:.2f}x"

    except Exception as e:
        logging.debug(f"   MGMT: <access denied>: {e}")
        pass

    return coverage_data


def _extract_cnv_data(sample_dir: Path) -> Dict[str, Any]:
    """Extract CNV analysis data using the same methods as the CNV component."""
    cnv_data = {
        "genetic_sex": "Not available",
        "bin_width": "Not available",
        "variance": "Not available",
        "gained": 0,
        "lost": 0,
    }

    try:
        # Load CNV_dict.npy for bin width and variance (same as CNV component)
        cnv_dict_npy = sample_dir / "CNV_dict.npy"
        if cnv_dict_npy.exists():
            try:
                import numpy as np

                cnv_dict = np.load(cnv_dict_npy, allow_pickle=True).item()
                if isinstance(cnv_dict, dict):
                    if "bin_width" in cnv_dict:
                        bin_width = cnv_dict["bin_width"]
                        if isinstance(bin_width, (int, float)):
                            cnv_data["bin_width"] = f"{bin_width:,}"
                    if "variance" in cnv_dict:
                        variance = cnv_dict["variance"]
                        if isinstance(variance, (int, float)):
                            cnv_data["variance"] = f"{variance:.3f}"
            except Exception as e:
                logging.debug(f"   MGMT: <access denied>: {e}")
                pass

        # Load XYestimate.pkl for genetic sex (same as CNV component)
        xy_pkl = sample_dir / "XYestimate.pkl"
        if xy_pkl.exists():
            try:
                import pickle

                with xy_pkl.open("rb") as f:
                    xy = pickle.load(f)
                if xy:
                    cnv_data["genetic_sex"] = str(xy)
            except Exception as e:
                logging.debug(f"   MGMT: <access denied>: {e}")
                pass

        # Look for CNV result files for gained/lost counts
        cnv_files = [
            "cnv_analysis_counter.txt",
            "cnv_results.csv",
            "cnv_summary.txt",
            "cnv_analysis_results.pkl",
        ]

        for cnv_file in cnv_files:
            file_path = sample_dir / cnv_file
            if file_path.exists():
                if cnv_file.endswith(".txt"):
                    try:
                        with open(file_path, "r") as f:
                            content = f.read().strip()
                            # Parse counter file for basic CNV info
                            if content.isdigit():
                                cnv_data["gained"] = int(
                                    content
                                )  # Assume this is gained regions
                                cnv_data["lost"] = 0  # Default value
                    except Exception as e:
                        logging.debug(f"   MGMT: <access denied>: {e}")
                        pass
                elif cnv_file.endswith(".csv"):
                    try:
                        with open(file_path, "r") as f:
                            reader = csv.DictReader(f)
                            for row in reader:
                                if "genetic_sex" in row:
                                    cnv_data["genetic_sex"] = row["genetic_sex"]
                                if "bin_width" in row:
                                    cnv_data["bin_width"] = row["bin_width"]
                                if "variance" in row:
                                    cnv_data["variance"] = row["variance"]
                                if "gained" in row:
                                    cnv_data["gained"] = int(row["gained"])
                                if "lost" in row:
                                    cnv_data["lost"] = int(row["lost"])
                    except Exception as e:
                        logging.debug(f"   MGMT: <access denied>: {e}")
                        pass
                break

    except Exception as e:
        logging.debug(f"   MGMT: <access denied>: {e}")
        pass

    return cnv_data


def _extract_mgmt_data(sample_dir: Path) -> Dict[str, Any]:
    """Extract MGMT analysis data using the same logic as the MGMT component."""
    mgmt_data = {
        "status": "Not available",
        "methylation_percent": "Not available",
        "average_methylation": "Not available",
        "prediction_score": "Not available",
        "cpg_sites": "Not available",
    }

    try:
        # Use the same file discovery logic as the MGMT component
        # First check for final_mgmt.csv files (highest priority)
        final_files = list(sample_dir.glob("final_mgmt.csv"))
        if final_files:
            latest_csv = final_files[0]
        else:
            # Fallback to numeric-prefixed files
            csv_files = list(sample_dir.glob("*_mgmt.csv"))
            if not csv_files:
                return mgmt_data

            # Find the latest MGMT CSV file (same logic as MGMT component)
            def _count_from_name(p: Path) -> int:
                try:
                    return int(p.name.split("_")[0])
                except Exception as e:
                    logging.debug(f"   MGMT: <access denied>: {e}")
                    return -1

            latest_csv = max(csv_files, key=_count_from_name)

        # Read the latest MGMT CSV file
        try:
            import pandas as pd

            required_cols = {"status", "pred", "average"}
            df = pd.read_csv(
                latest_csv,
                usecols=lambda c: c in required_cols,
                nrows=1,
                dtype=str,
            )

            # Extract data using the same column names as the MGMT component
            if "status" in df.columns:
                status = str(df.get("status", pd.Series(["Unknown"])).iloc[0])
                mgmt_data["status"] = status
            elif "pred" in df.columns:
                # Convert prediction score to status (same logic as MGMT component)
                try:
                    pred_val = float(df.get("pred", pd.Series([0.0])).iloc[0])
                    if pred_val > 0.1:
                        mgmt_data["status"] = "Methylated"
                    else:
                        mgmt_data["status"] = "Unmethylated"
                except Exception as e:
                    logging.debug(f"   MGMT: <access denied>: {e}")
                    pass

            if "average" in df.columns:
                try:
                    avg_val = float(df.get("average", pd.Series([0.0])).iloc[0])
                    mgmt_data["average_methylation"] = f"{avg_val:.1f}%"
                    mgmt_data["methylation_percent"] = f"{avg_val:.1f}%"
                except Exception as e:
                    logging.debug(f"   MGMT: <access denied>: {e}")
                    pass

            if "pred" in df.columns:
                try:
                    pred_val = float(df.get("pred", pd.Series([0.0])).iloc[0])
                    mgmt_data["prediction_score"] = f"{pred_val:.1f}%"
                except Exception as e:
                    logging.debug(f"   MGMT: <access denied>: {e}")
                    pass

            # Count CPG sites from the BED file if available
            current_count = _count_from_name(latest_csv)
            bed_path = sample_dir / f"{current_count}_mgmt.bed"
            if not bed_path.exists():
                alt_bed = sample_dir / f"{current_count}_mgmt_mgmt.bed"
                bed_path = alt_bed if alt_bed.exists() else bed_path

            if bed_path.exists():
                try:
                    # Count lines in BED file to estimate CPG sites
                    with open(bed_path, "r") as f:
                        line_count = sum(1 for line in f if line.strip())
                    mgmt_data["cpg_sites"] = line_count
                except Exception as e:
                    logging.debug(f"   MGMT: <access denied>: {e}")
                    pass

        except Exception as e:
            logging.debug(f"   MGMT: <access denied>: {e}")
            pass

    except Exception as e:
        logging.debug(f"   MGMT: <access denied>: {e}")
        pass

    return mgmt_data


def _extract_classification_data(sample_dir: Path) -> Dict[str, Any]:
    """Extract classification data from various score files."""
    classification_data = {
        "sturgeon": {
            "classification": "Not available",
            "confidence": 0.0,
            "confidence_level": "Not available",
            "features": 0,
        },
        "nanodx": {
            "classification": "Not available",
            "confidence": 0.0,
            "confidence_level": "Not available",
            "features": 0,
        },
        "pannanodx": {
            "classification": "Not available",
            "confidence": 0.0,
            "confidence_level": "Not available",
            "features": 0,
        },
        "random_forest": {
            "classification": "Not available",
            "confidence": 0.0,
            "confidence_level": "Not available",
            "features": 0,
        },
    }

    try:
        # Extract Sturgeon data
        sturgeon_file = sample_dir / "sturgeon_scores.csv"
        if sturgeon_file.exists():
            try:
                with open(sturgeon_file, "r") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        # Get the latest row (last one)
                        features = int(row.get("number_probes", 0))
                        # Find the highest scoring classification
                        max_score = 0.0
                        best_class = "Unknown"
                        for col, value in row.items():
                            if col not in ["timestamp", "number_probes"]:
                                try:
                                    score = float(value)
                                    if score > max_score:
                                        max_score = score
                                        best_class = col
                                except Exception as e:
                                    logging.debug(f"   Sturgeon: <access denied>: {e}")
                                    pass

                        classification_data["sturgeon"] = {
                            "classification": best_class,
                            "confidence": max_score * 100,
                            "confidence_level": _get_confidence_level(max_score * 100, "sturgeon"),
                            "features": features,
                        }
            except Exception as e:
                logging.debug(f"   Sturgeon: <access denied>: {e}")
                pass

        # Extract NanoDX data
        nanodx_file = sample_dir / "NanoDX_scores.csv"
        if nanodx_file.exists():
            try:
                with open(nanodx_file, "r") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        # Get the latest row (last one)
                        features = int(row.get("number_probes", 0))
                        # Find the highest scoring classification
                        max_score = 0.0
                        best_class = "Unknown"
                        for col, value in row.items():
                            if col not in ["timestamp", "number_probes"]:
                                try:
                                    score = float(value)
                                    if score > max_score:
                                        max_score = score
                                        best_class = col
                                except Exception as e:
                                    logging.debug(f"   NanoDX: <access denied>: {e}")
                                    pass

                        classification_data["nanodx"] = {
                            "classification": best_class,
                            "confidence": max_score * 100,
                            "confidence_level": _get_confidence_level(max_score * 100, "nanodx"),
                            "features": features,
                        }
            except Exception as e:
                logging.debug(f"   NanoDX: <access denied>: {e}")
                pass

        # Extract PanNanoDX data
        pannanodx_file = sample_dir / "PanNanoDX_scores.csv"
        if pannanodx_file.exists():
            try:
                with open(pannanodx_file, "r") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        # Get the latest row (last one)
                        features = int(row.get("number_probes", 0))
                        # Find the highest scoring classification
                        max_score = 0.0
                        best_class = "Unknown"
                        for col, value in row.items():
                            if col not in ["timestamp", "number_probes"]:
                                try:
                                    score = float(value)
                                    if score > max_score:
                                        max_score = score
                                        best_class = col
                                except Exception as e:
                                    logging.debug(f"   PanNanoDX: <access denied>: {e}")
                                    pass

                        classification_data["pannanodx"] = {
                            "classification": best_class,
                            "confidence": max_score * 100,
                            "confidence_level": _get_confidence_level(max_score * 100, "pannanodx"),
                            "features": features,
                        }
            except Exception as e:
                logging.debug(f"   PanNanoDX: <access denied>: {e}")
                pass

        # Extract Random Forest data
        rf_file = sample_dir / "random_forest_scores.csv"
        if rf_file.exists():
            try:
                with open(rf_file, "r") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        # Get the latest row (last one)
                        features = int(row.get("number_probes", 0))
                        # Find the highest scoring classification
                        max_score = 0.0
                        best_class = "Unknown"
                        for col, value in row.items():
                            if col not in ["timestamp", "number_probes"]:
                                try:
                                    score = float(value)
                                    if score > max_score:
                                        max_score = score
                                        best_class = col
                                except Exception as e:
                                    logging.debug(f"   Random Forest: <access denied>: {e}")
                                    pass

                        # Some random forest outputs are already in percent (0-100),
                        # while others are fractional (0-1). Normalize to percent.
                        confidence_percent = (
                            max_score * 100 if max_score <= 1 else max_score
                        )
                        classification_data["random_forest"] = {
                            "classification": best_class,
                            "confidence": confidence_percent,
                            "confidence_level": _get_confidence_level(
                                confidence_percent, "random_forest"
                            ),
                            "features": features,
                        }
            except Exception as e:
                logging.debug(f"   Random Forest: <access denied>: {e}")
                pass

    except Exception as e:
        logging.debug(f"   Random Forest: <access denied>: {e}")
        pass

    return classification_data


def _extract_fusion_data(sample_dir: Path) -> Dict[str, Any]:
    """Extract fusion analysis data from generated summary files."""
    fusion_data = {
        "target_fusions": 0, 
        "genome_fusions": 0,
        "target_pairs": 0,
        "target_groups": 0,
        "genome_pairs": 0,
        "genome_groups": 0
    }
    
    logging.info(f"[Summary] _extract_fusion_data() called with sample_dir: {sample_dir}")

    try:
        # Debug: List all fusion-related files
        fusion_files = list(sample_dir.glob("*fusion*"))
        logging.info(f"[Summary] Found fusion files: {[f.name for f in fusion_files]}")
        
        # Debug: Check if genome-wide processed file exists
        genome_file = sample_dir / "fusion_candidates_all_processed.pkl"
        logging.info(f"[Summary] Genome-wide processed file exists: {genome_file.exists()}")
        if genome_file.exists():
            logging.info(f"[Summary] Genome-wide processed file size: {genome_file.stat().st_size} bytes")
        # First try to read from the new fusion_summary.csv file
        summary_file = sample_dir / "fusion_summary.csv"
        if summary_file.exists():
            try:
                with open(summary_file, "r") as f:
                    content = f.read()
                    
                with open(summary_file, "r") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        fusion_data["target_fusions"] = int(row.get("target_fusions", 0))
                        fusion_data["genome_fusions"] = int(row.get("genome_fusions", 0))
                        break
                
                # Check if the summary file has incorrect genome-wide count (0)
                if fusion_data["genome_fusions"] == 0:
                    # Try to regenerate the summary file with correct data
                    try:
                        from robin.gui.components.fusion import _generate_summary_files_from_pickle
                        if _generate_summary_files_from_pickle(sample_dir, force_regenerate=True):
                            # Re-read the regenerated summary file
                            with open(summary_file, "r") as f:
                                reader = csv.DictReader(f)
                                for row in reader:
                                    fusion_data["target_fusions"] = int(row.get("target_fusions", 0))
                                    fusion_data["genome_fusions"] = int(row.get("genome_fusions", 0))
                                    break
                            return fusion_data
                    except Exception as e:
                        pass  # Don't return early, continue to pickle file loading
                else:
                    logging.info(f"[Summary] Fusion data extracted from summary file - target: {fusion_data['target_fusions']}, genome: {fusion_data['genome_fusions']}")
                    logging.info(f"[Summary] Summary file path: {summary_file}")
                    return fusion_data
            except Exception as e:
                logging.debug(f"   Fusion: Failed to read summary file: {e}")
        
        # If summary files don't exist, try to generate them from pickle files
        try:
            from robin.gui.components.fusion import _generate_summary_files_from_pickle
            if _generate_summary_files_from_pickle(sample_dir, force_regenerate=False):
                # Try reading the newly generated summary file
                if summary_file.exists():
                    with open(summary_file, "r") as f:
                        reader = csv.DictReader(f)
                        for row in reader:
                            fusion_data["target_fusions"] = int(row.get("target_fusions", 0))
                            fusion_data["genome_fusions"] = int(row.get("genome_fusions", 0))
                            break
                    logging.info(f"[Summary] Fusion data extracted from generated summary file - target: {fusion_data['target_fusions']}, genome: {fusion_data['genome_fusions']}")
                    return fusion_data
        except Exception as e:
            logging.debug(f"   Fusion: Failed to generate summary files from pickle: {e}")
        
        # If still no data, try to load directly from pickle files and count gene_pairs
        try:
            from robin.gui.components.fusion import _load_processed_pickle, _count_unique_fusion_pairs, _count_unique_fusion_groups
            target_file = sample_dir / "fusion_candidates_master_processed.pkl"
            genome_file = sample_dir / "fusion_candidates_all_processed.pkl"
            
            
            target_data = _load_processed_pickle(target_file)
            genome_data = _load_processed_pickle(genome_file)
            
            
            if target_data and isinstance(target_data, dict):
                # Use filtered counts to match what's displayed in fusion section
                fusion_data["target_fusions"] = _count_unique_fusion_pairs(target_data)
                fusion_data["target_pairs"] = _count_unique_fusion_pairs(target_data)
                fusion_data["target_groups"] = _count_unique_fusion_groups(target_data)
            
            if genome_data and isinstance(genome_data, dict):
                # Use filtered counts to match what's displayed in fusion section
                fusion_data["genome_fusions"] = _count_unique_fusion_pairs(genome_data)
                fusion_data["genome_pairs"] = _count_unique_fusion_pairs(genome_data)
                fusion_data["genome_groups"] = _count_unique_fusion_groups(genome_data)
            
            logging.info(f"[Summary] Fusion data loaded directly from pickle files - target: {fusion_data['target_fusions']} fusions, {fusion_data['target_pairs']} pairs, {fusion_data['target_groups']} groups; genome: {fusion_data['genome_fusions']} fusions, {fusion_data['genome_pairs']} pairs, {fusion_data['genome_groups']} groups")
            return fusion_data
        except Exception as e:
            logging.debug(f"   Fusion: Failed to load directly from pickle files: {e}")
        
        # Fallback to individual fusion_results.csv file
        fusion_results_file = sample_dir / "fusion_results.csv"
        if fusion_results_file.exists():
            try:
                with open(fusion_results_file, "r") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        fusion_data["target_fusions"] = int(row.get("target_fusions", 0))
                        fusion_data["genome_fusions"] = int(row.get("genome_fusions", 0))
                        break
                logging.debug(f"[Summary] Fusion data extracted from results file - target: {fusion_data['target_fusions']}, genome: {fusion_data['genome_fusions']}")
                return fusion_data
            except Exception as e:
                logging.debug(f"   Fusion: Failed to read results file: {e}")
        
        # Final fallback to legacy sv_count.txt file
        sv_count_file = sample_dir / "sv_count.txt"
        if sv_count_file.exists():
            try:
                with open(sv_count_file, "r") as f:
                    content = f.read().strip()
                    if content.isdigit():
                        # Legacy behavior - assume this is genome-wide count
                        fusion_data["genome_fusions"] = int(content)
                        fusion_data["target_fusions"] = 0
                logging.debug(f"[Summary] Fusion data extracted from legacy sv_count file - target: {fusion_data['target_fusions']}, genome: {fusion_data['genome_fusions']}")
            except Exception as e:
                logging.debug(f"   Fusion: Failed to read legacy sv_count file: {e}")

    except Exception as e:
        logging.debug(f"   Fusion: <access denied>: {e}")

    return fusion_data


def _get_confidence_color(confidence_level: str) -> str:
    """Get the appropriate color class for confidence level."""
    if "High" in confidence_level:
        return "text-green-600"
    elif "Medium" in confidence_level:
        return "text-blue-600"
    elif "Low" in confidence_level:
        return "text-orange-600"
    else:
        return "text-gray-600"


def _get_quality_color(quality: str) -> str:
    """Get the appropriate color class for quality status."""
    if quality == "Excellent":
        return "text-green-600"
    elif quality == "Good":
        return "text-blue-600"
    elif quality == "Moderate":
        return "text-orange-600"
    elif quality == "Insufficient":
        return "text-red-600"
    else:
        return "text-gray-600"


def _get_coverage_badge_color(coverage: str) -> str:
    """Get the appropriate badge color for coverage value."""
    try:
        cov_val = float(coverage.replace("x", ""))
        if cov_val >= 30:
            return "positive"
        elif cov_val >= 20:
            return "primary"
        elif cov_val >= 10:
            return "warning"
        else:
            return "negative"
    except Exception as e:
        logging.debug(f"   Coverage: <access denied>: {e}")
        return "grey"


def _get_methylation_color(status: str) -> str:
    """Get the appropriate color class for methylation status."""
    if "Methylated" in status:
        return "text-green-600"
    elif "Unmethylated" in status:
        return "text-orange-600"
    else:
        return "text-gray-600"


def _get_confidence_level(confidence: float, classifier: str = "default") -> str:
    """Get the confidence level based on classifier-specific thresholds."""
    return get_classifier_confidence_level(classifier, confidence)


def _get_confidence_badge_classes(confidence: float) -> str:
    """Get CSS classes for confidence level badge based on confidence score."""
    if confidence >= 80:
        return "px-2 py-1 text-xs font-medium rounded-full bg-green-100 text-green-800"
    elif confidence >= 50:
        return "px-2 py-1 text-xs font-medium rounded-full bg-yellow-100 text-yellow-800"
    else:
        return "px-2 py-1 text-xs font-medium rounded-full bg-red-100 text-red-800"


def _get_methylation_badge_color(methylation: str) -> str:
    """Get the appropriate badge color for methylation value."""
    try:
        meth_val = float(methylation.replace("%", ""))
        if meth_val > 10:
            return "positive"
        elif meth_val > 5:
            return "warning"
        else:
            return "negative"
    except Exception as e:
        logging.debug(f"   Methylation: <access denied>: {e}")
        return "grey"


