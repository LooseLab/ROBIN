from __future__ import annotations

import logging
import os
import threading
from typing import Any, Callable, Dict, List, Optional
from pathlib import Path
from robin.analysis.snp_processing import VariantDisplayStore, VariantTableFilters
from robin.gui.client_notify import run_javascript_when_connected, schedule_after_page_sent
from robin.utils.clinvar_manager import compare_sample_clinvar_to_installed

try:
    from nicegui import ui
except ImportError:  # pragma: no cover
    ui = None

logger = logging.getLogger(__name__)

# Default SNP/INDEL table and detail views for ClinVar-annotated variants.
VARIANT_TABLE_FIELDS = [
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "Gene_Name",
    "HGVS.p",
    "Annotation",
    "Annotation_Impact",
    "CLNSIG",
    "ONC",
    "SCI",
    "ONCDN",
    "SCIDN",
    "FILTER",
    "QUAL",
    "GT",
    "is_clinvar_significant",
    "details",
    "action",
]
VARIANT_WIDE_FIELDS = {
    "Gene_Name",
    "Annotation",
    "HGVS.p",
    "CLNSIG",
    "ONC",
    "SCI",
    "ONCDN",
    "SCIDN",
}
VARIANT_DETAIL_FIELDS = [
    "CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "Gene_Name",
    "HGVS.c",
    "HGVS.p",
    "Annotation",
    "Annotation_Impact",
    "is_clinvar_significant",
    "CLNSIG",
    "ONC",
    "SCI",
    "ONCDN",
    "SCIDN",
    "CLNDN",
    "is_pathogenic",
    "is_vus",
    "FILTER",
    "QUAL",
    "GT",
]
VARIANT_COLUMN_LABELS = {
    "is_clinvar_significant": "ClinVar significant",
    "is_pathogenic": "Germline pathogenic",
    "is_vus": "VUS",
    "ONCDN": "Oncogenic disease",
    "SCIDN": "Somatic disease",
    "CLNDN": "Germline disease",
    "CLNSIG": "CLNSIG (germline)",
    "ONC": "ONC (oncogenic)",
    "SCI": "SCI (somatic tier)",
}


def _apply_variant_column_labels(columns: List[Dict[str, Any]]) -> None:
    for col in columns:
        field = col.get("field")
        if field in VARIANT_COLUMN_LABELS:
            col["label"] = VARIANT_COLUMN_LABELS[field]


def navigate_igv_to_snp(chrom: str, pos: int, flank: int = 100) -> None:
    """
    Navigate IGV browser to a specific SNP location.
    
    Args:
        chrom: Chromosome name (e.g., "chr1")
        pos: Position on the chromosome
        flank: Number of bases to include on each side (default: 100 bp)
    """
    try:
        # Ensure chromosome name has 'chr' prefix if needed
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        
        # Calculate window around the SNP
        start = max(1, pos - flank)
        end = pos + flank
        
        region = f"{chrom}:{start}-{end}"
        
        # Escape region string for JavaScript
        escaped_region = region.replace('"', '\\"').replace("'", "\\'")
        
        js_navigate = f"""
            (function() {{
                try {{
                    if (window.lj_igv && window.lj_igv_browser_ready) {{
                        console.log('[IGV] Navigating to SNP region: {escaped_region}');
                        window.lj_igv.search('{escaped_region}');
                    }} else {{
                        console.warn('[IGV] Browser not ready yet, will navigate when ready');
                        setTimeout(function() {{
                            if (window.lj_igv && window.lj_igv_browser_ready) {{
                                window.lj_igv.search('{escaped_region}');
                            }}
                        }}, 1000);
                    }}
                }} catch (error) {{
                    console.error('[IGV] Error navigating to SNP: ' + error);
                }}
            }})();
        """
        
        run_javascript_when_connected(js_navigate)
        
    except Exception as e:
        logger.error(f"Error navigating IGV to SNP {chrom}:{pos}: {e}")


def _resolve_reference_genome(launcher: Any) -> Optional[str]:
    reference_genome = None
    workflow_runner = getattr(launcher, "workflow_runner", None)
    if workflow_runner is not None:
        reference_genome = getattr(workflow_runner, "reference", None)
    if not reference_genome:
        env_reference = os.environ.get("robin_REFERENCE")
        if env_reference and os.path.exists(env_reference):
            reference_genome = env_reference
    return reference_genome


def _submit_snp_workflow_job(
    launcher: Any,
    sample_dir: Path,
    *,
    annotation_only: bool = False,
    force_regenerate: bool = False,
) -> bool:
    """Queue SNP analysis or annotation-only rerun through the workflow runner."""
    from robin.analysis.target_analysis import snp_analysis_handler
    from robin.workflow_simple import Job, WorkflowContext

    work_dir = str(sample_dir.parent)
    sample_id = sample_dir.name
    reference_genome = _resolve_reference_genome(launcher)

    metadata: Dict[str, Any] = {
        "work_dir": work_dir,
        "threads": 4,
        "force_regenerate": force_regenerate,
        "annotation_only": annotation_only,
        "reference": reference_genome,
    }

    try:
        master_csv = sample_dir / "master.csv"
        if master_csv.exists():
            import csv

            with master_csv.open("r", newline="", encoding="utf-8") as fh:
                reader = csv.DictReader(fh)
                first_row = next(reader, None)
                if first_row:
                    panel = str(first_row.get("analysis_panel", "") or "").strip()
                    if panel:
                        metadata["target_panel"] = panel
    except Exception:
        pass

    context = WorkflowContext(filepath=str(sample_dir), metadata=metadata)
    context.get_sample_id = lambda: sample_id
    context.add_result = lambda key, value: None
    context.add_error = lambda key, value: None

    job = Job(
        job_id=hash(f"snp_analysis_{sample_id}_{annotation_only}") % 1000000,
        job_type="snp_analysis",
        context=context,
        origin="manual",
        workflow=["slow:snp_analysis"],
    )

    workflow_runner = getattr(launcher, "workflow_runner", None)
    if workflow_runner is not None and hasattr(workflow_runner, "submit_snp_analysis_job"):
        try:
            success = workflow_runner.submit_snp_analysis_job(
                sample_dir=str(sample_dir),
                sample_id=sample_id,
                reference=str(reference_genome) if reference_genome else None,
                threads=4,
                force_regenerate=force_regenerate,
                annotation_only=annotation_only,
            )
            if success:
                return True
        except Exception as exc:
            logger.warning("Workflow SNP submission failed: %s", exc)

    def _run_handler() -> None:
        try:
            snp_analysis_handler(job, work_dir=work_dir)
        except Exception as exc:
            logger.error("Direct SNP handler failed: %s", exc)

    threading.Thread(target=_run_handler, daemon=True).start()
    return True


def _add_clinvar_annotation_controls(
    launcher: Any,
    sample_dir: Path,
    *,
    compact: bool = False,
    on_status_change: Optional[Callable[[str, str], None]] = None,
    clinvar_status: Optional[Dict[str, Any]] = None,
) -> None:
    """Render ClinVar provenance labels and optional re-annotation trigger."""
    if ui is None:
        return

    status = clinvar_status if clinvar_status is not None else compare_sample_clinvar_to_installed(sample_dir)
    installed_label = status.get("installed_label", "ClinVar (version unknown)")
    sample_label = status.get("sample_label", "Annotation release not recorded")
    is_stale = bool(status.get("is_stale"))
    can_reannotate = bool(status.get("can_reannotate"))

    def _set_status(text: str, tone: str = "meta") -> None:
        if on_status_change is not None:
            on_status_change(text, tone)

    with ui.row().classes(
        "w-full gap-2 mb-2 flex-wrap items-center"
        if not compact
        else "w-full gap-2 mb-2 flex-wrap items-center"
    ):
        ui.label(f"Annotated with: {sample_label}").classes("classification-insight-meta")
        ui.label(f"Installed: {installed_label}").classes(
            "classification-insight-level classification-insight-level--low w-auto"
            if is_stale
            else "classification-insight-meta"
        )
        if is_stale:
            ui.label("Newer ClinVar is installed").classes(
                "classification-insight-level classification-insight-level--low w-auto"
            )

        if not compact:

            def _start_reannotation() -> None:
                if not can_reannotate:
                    ui.notify(
                        "Existing Clair3 outputs were not found. Run full SNP analysis first.",
                        type="warning",
                    )
                    return
                submitted = _submit_snp_workflow_job(
                    launcher,
                    sample_dir,
                    annotation_only=True,
                )
                if submitted:
                    ui.notify(
                        "ClinVar re-annotation started (snpEff/SnpSift). "
                        "Refresh this page when the job completes.",
                        type="info",
                    )
                    _set_status("ClinVar re-annotation running…", "running")
                else:
                    ui.notify("Could not start ClinVar re-annotation.", type="negative")
                    _set_status("Failed to start ClinVar re-annotation", "error")

            reannotate_button = ui.button(
                "Re-annotate with current ClinVar",
                on_click=_start_reannotation,
            ).props("dense no-caps outline color=secondary")
            if not can_reannotate:
                reannotate_button.disable()
                reannotate_button.props('title="Requires clair3/output_done.vcf.gz"')


VARIANT_TABLE_PAGE_MAX = 250
_UI_ONLY_FIELDS = {"action", "details", "__row_id", "__row_idx"}


def _to_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip().replace(",", "")
    if text == "":
        return None
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def _compact_variant_row(
    row: Dict[str, Any],
    visible_fields: List[str],
    *,
    wide_fields: set[str],
    max_field_length: int = 80,
) -> Dict[str, Any]:
    compact: Dict[str, Any] = {}
    for field_name in visible_fields:
        value = row.get(field_name, "")
        text = "" if value is None else str(value)
        if field_name in wide_fields and len(text) > max_field_length:
            text = f"{text[:max_field_length - 1]}..."
        compact[field_name] = text
    compact["details"] = " "
    compact["action"] = " "
    compact["__row_idx"] = int(row.get("__row_idx", -1))
    return compact


def _mount_paged_variant_table(
    store: VariantDisplayStore,
    *,
    details_title: str,
    event_prefix: str,
    navigate_region: Callable[[str], None],
) -> None:
    """Render one paged SNP/INDEL table backed by ``VariantDisplayStore``."""
    from robin.gui.theme import (
        clamp_qtable_server_pagination,
        styled_server_paged_table,
        wire_qtable_server_pagination_handlers,
    )

    columns = list(store.columns or [])
    column_lookup = {
        col.get("field"): col
        for col in columns
        if isinstance(col, dict) and col.get("field")
    }
    visible_fields = [
        field_name for field_name in VARIANT_TABLE_FIELDS if field_name in column_lookup
    ]
    if not visible_fields:
        visible_fields = [
            col.get("field")
            for col in columns[:12]
            if isinstance(col, dict) and col.get("field") not in _UI_ONLY_FIELDS
        ]

    display_columns = [column_lookup[field_name].copy() for field_name in visible_fields if field_name in column_lookup]
    if "details" not in {col.get("field") for col in display_columns}:
        display_columns.append(
            {"name": "details", "label": "Details", "field": "details", "sortable": False}
        )
    if "action" not in {col.get("field") for col in display_columns}:
        display_columns.append(
            {"name": "action", "label": "View in IGV", "field": "action", "sortable": False}
        )
    _apply_variant_column_labels(display_columns)

    total_variants = store.total_variants
    summary = store.summary
    pathogenic_count = int(summary.get("pathogenic_variants") or 0)
    significant_count = int(
        summary.get("clinvar_significant_variants", pathogenic_count) or 0
    )

    with ui.row().classes("w-full gap-3 mb-3 flex-wrap items-baseline"):
        ui.label(f"Total variants: {total_variants:,}").classes("classification-insight-meta")
        if total_variants > 5_000:
            ui.label(
                "Large dataset: paging and filters query Parquet on the server "
                "(only the current page is sent to the browser)."
            ).classes("classification-insight-meta w-full")
        if significant_count > 0:
            ui.label(f"ClinVar significant variants: {significant_count:,}").classes(
                "classification-insight-level classification-insight-level--low w-auto"
            )
        elif pathogenic_count > 0:
            ui.label(f"Pathogenic variants: {pathogenic_count:,}").classes(
                "classification-insight-level classification-insight-level--low w-auto"
            )

    has_dp = store.has_column("DP")
    with ui.row().classes("w-full gap-2 mb-2 flex-wrap items-end"):
        pass_only = ui.checkbox("PASS only").props("dense")
        significant_only = ui.checkbox("ClinVar significant only").props("dense")
        min_qual = ui.number("Min QUAL", value=None).props("dense outlined clearable").classes("w-32")
        min_dp = (
            ui.number("Min DP", value=None).props("dense outlined clearable").classes("w-32")
            if has_dp
            else None
        )
        search = ui.input("Search (gene/variant)").props(
            "dense outlined clearable debounce=400"
        ).classes("w-64")
        reset_button = ui.button("Reset").props("dense no-caps")

    page_state: Dict[str, Any] = {
        "filters": VariantTableFilters(search_fields=tuple(visible_fields)),
    }

    def _current_filters() -> VariantTableFilters:
        return VariantTableFilters(
            pass_only=bool(getattr(pass_only, "value", False)),
            significant_only=bool(getattr(significant_only, "value", False)),
            min_qual=_to_float(getattr(min_qual, "value", None)),
            min_dp=_to_float(getattr(min_dp, "value", None)) if has_dp else None,
            search_text=str(getattr(search, "value", "") or ""),
            search_fields=tuple(visible_fields),
        )

    if page_state["filters"].is_active:
        _rows0, total_filtered0 = store.page(page_state["filters"], offset=0, limit=0)
    else:
        total_filtered0 = total_variants
    init_pagination = clamp_qtable_server_pagination(
        {
            "sortBy": None,
            "descending": False,
            "page": 1,
            "rowsPerPage": 100,
            "rowsNumber": total_filtered0,
        },
        rows_number=total_filtered0,
        rows_per_page_default=100,
        rows_per_page_max=VARIANT_TABLE_PAGE_MAX,
    )
    _, table = styled_server_paged_table(
        columns=display_columns,
        rows=[],
        pagination=init_pagination,
        row_key="__row_idx",
        class_size="table-xs",
    )
    filtered_count_label = ui.label(
        f"{total_filtered0:,} variants match filters (of {total_variants:,} total)"
    ).classes("classification-insight-meta")

    def _fill_from_pagination(pag: Dict[str, Any]) -> None:
        filters = page_state["filters"]
        pag = clamp_qtable_server_pagination(
            pag,
            rows_number=int(pag.get("rowsNumber") or total_variants),
            rows_per_page_default=100,
            rows_per_page_max=VARIANT_TABLE_PAGE_MAX,
        )
        rpp = int(pag["rowsPerPage"])
        page = int(pag["page"])
        start = (page - 1) * rpp
        rows, total_filtered = store.page(filters, offset=start, limit=rpp)
        pag = clamp_qtable_server_pagination(
            pag,
            rows_number=total_filtered,
            rows_per_page_default=100,
            rows_per_page_max=VARIANT_TABLE_PAGE_MAX,
        )
        if int(pag["page"]) != page or int(pag["rowsPerPage"]) != rpp:
            rpp = int(pag["rowsPerPage"])
            page = int(pag["page"])
            start = (page - 1) * rpp
            rows, total_filtered = store.page(filters, offset=start, limit=rpp)
        table.rows = [
            _compact_variant_row(
                row,
                visible_fields,
                wide_fields=VARIANT_WIDE_FIELDS,
            )
            for row in rows
        ]
        pag["rowsNumber"] = total_filtered
        table.pagination = pag
        filtered_count_label.text = (
            f"{total_filtered:,} variants match filters (of {total_variants:,} total)"
        )
        table.update()

    wire_qtable_server_pagination_handlers(table, _fill_from_pagination)

    def _apply_filters() -> None:
        page_state["filters"] = _current_filters()
        pag = clamp_qtable_server_pagination(
            dict(table.pagination),
            rows_number=0,
            rows_per_page_default=100,
            rows_per_page_max=VARIANT_TABLE_PAGE_MAX,
        )
        pag["page"] = 1
        _fill_from_pagination(pag)

    pass_only.on("update:model-value", lambda _e: _apply_filters())
    significant_only.on("update:model-value", lambda _e: _apply_filters())
    min_qual.on("update:model-value", lambda _e: _apply_filters())
    if has_dp and min_dp is not None:
        min_dp.on("update:model-value", lambda _e: _apply_filters())
    search.on("update:model-value", lambda _e: _apply_filters())
    reset_button.on_click(
        lambda: (
            setattr(pass_only, "value", False),
            setattr(significant_only, "value", False),
            setattr(min_qual, "value", None),
            has_dp and setattr(min_dp, "value", None),
            setattr(search, "value", ""),
            _apply_filters(),
        )
    )

    with ui.dialog() as details_dialog, ui.card().classes(
        "robin-dialog-surface w-[95vw] max-w-6xl max-h-[85vh] overflow-auto "
        "p-4 md:p-5"
    ):
        ui.label(details_title).classes("classification-insight-heading text-headline-small")
        ui.separator().classes("mgmt-detail-separator")
        details_container = ui.column().classes("w-full gap-2")
        with ui.row().classes("w-full justify-end pt-2"):
            ui.button("Close", on_click=details_dialog.close).props(
                "color=primary no-caps outline"
            )

    def show_variant_details(row_idx: int) -> None:
        row_data = store.fetch_row(int(row_idx))
        if not row_data:
            ui.notify("Variant details not found.", type="warning")
            return
        details_container.clear()
        ordered_fields = list(VARIANT_DETAIL_FIELDS) + sorted(
            [
                key
                for key in row_data.keys()
                if key not in VARIANT_DETAIL_FIELDS and key not in _UI_ONLY_FIELDS
            ]
        )
        with details_container:
            for field_name in ordered_fields:
                value = row_data.get(field_name, "")
                if value is None or str(value) == "":
                    continue
                label = VARIANT_COLUMN_LABELS.get(field_name, field_name.replace("_", " "))
                with ui.row().classes("w-full items-start gap-2"):
                    ui.label(f"{label}:").classes("text-xs font-semibold min-w-[180px]")
                    ui.label(str(value)).classes("text-xs whitespace-pre-wrap break-all flex-1")
        details_dialog.open()

    table.add_slot(
        "body-cell-details",
        f"""
<q-td key="details" :props="props">
  <q-btn
    icon="description"
    size="sm"
    dense
    flat
    color="secondary"
    @click="$parent.$emit('{event_prefix}-show-details', props.row.__row_idx)"
    title="Show full details"
  />
</q-td>
""",
    )
    table.add_slot(
        "body-cell-action",
        f"""
<q-td key="action" :props="props">
  <q-btn
    icon="visibility"
    size="sm"
    dense
    flat
    color="primary"
    @click="$parent.$emit('{event_prefix}-view-igv', props.row.__row_idx)"
    title="View in IGV"
  />
</q-td>
""",
    )

    def on_view_igv(e: Any) -> None:
        try:
            row_idx = getattr(e, "args", None)
            if row_idx is None:
                return
            row = store.fetch_row(int(row_idx))
            if not row:
                return
            chrom = str(row.get("CHROM", "")).strip()
            pos_text = str(row.get("POS", "")).replace(",", "").strip()
            if chrom and pos_text:
                navigate_region(f"{chrom}:{pos_text}")
        except Exception as exc:
            logger.debug("Error handling %s IGV view: %s", event_prefix, exc)

    def on_show_details(e: Any) -> None:
        try:
            row_idx = getattr(e, "args", None)
            if row_idx is not None:
                show_variant_details(int(row_idx))
        except Exception as exc:
            logger.debug("Error handling %s details view: %s", event_prefix, exc)

    table.on(f"{event_prefix}-view-igv", on_view_igv)
    table.on(f"{event_prefix}-show-details", on_show_details)
    for col in table.columns:
        col["sortable"] = False
    _fill_from_pagination(init_pagination)
    try:
        def _cleanup() -> None:
            table.rows = []
        ui.context.client.on_disconnect(_cleanup)
    except Exception:
        pass


def add_snp_section(launcher: Any, sample_dir: Path) -> None:
    """
    Add SNP analysis section to the sample details page.

    Args:
        launcher: The GUI launcher instance
        sample_dir: Path to the sample directory
    """
    if not sample_dir or not sample_dir.exists():
        return

    clair3_dir = sample_dir / "clair3"
    if not clair3_dir.exists():
        return

    snp_store = VariantDisplayStore.open_snp(clair3_dir)
    if snp_store is None:
        with ui.element("div").classes("classification-insight-shell w-full min-w-0"):
            ui.label("SNP analysis").classes(
                "classification-insight-heading text-headline-small"
            )
            with ui.element("div").classes("classification-insight-card w-full min-w-0"):
                with ui.column().classes("w-full min-w-0 gap-2 p-2 md:p-3"):
                    ui.label(
                        "Precomputed SNP display data was not found. "
                        "Ensure SNP analysis has completed."
                    ).classes("classification-insight-meta")
        return

    snp_regions_map = snp_store.regions_map

    def navigate_to_snp_region(snp_key: str) -> None:
        if snp_key in snp_regions_map:
            snp_data = snp_regions_map[snp_key]
            navigate_igv_to_snp(snp_data["chrom"], snp_data["pos"])
            return
        try:
            chrom, pos_str = snp_key.split(":", 1)
            navigate_igv_to_snp(chrom, int(str(pos_str).replace(",", "")))
        except (TypeError, ValueError):
            logger.debug("Invalid SNP key for IGV navigation: %s", snp_key)

    clinvar_status = compare_sample_clinvar_to_installed(sample_dir)

    with ui.element("div").classes("classification-insight-shell w-full min-w-0"):
        ui.label("SNP analysis").classes(
            "classification-insight-heading text-headline-small"
        )

        snp_clinvar_status_label: Dict[str, Any] = {"element": None}

        def _update_snp_clinvar_status(text: str, tone: str = "meta") -> None:
            label = snp_clinvar_status_label.get("element")
            if label is None:
                return
            label.text = text
            classes = "classification-insight-meta"
            if tone == "running":
                classes = "text-sm text-blue-600"
            elif tone == "error":
                classes = "classification-insight-level classification-insight-level--low w-full"
            label.classes(replace=classes)

        with ui.element("div").classes("classification-insight-card w-full min-w-0"):
            with ui.column().classes("w-full min-w-0 gap-2 p-2 md:p-3"):
                _add_clinvar_annotation_controls(
                    launcher,
                    sample_dir,
                    on_status_change=_update_snp_clinvar_status,
                    clinvar_status=clinvar_status,
                )
                snp_clinvar_status_label["element"] = ui.label("").classes(
                    "classification-insight-meta"
                )

        _mount_paged_variant_table(
            snp_store,
            details_title="Variant details",
            event_prefix="snp",
            navigate_region=navigate_to_snp_region,
        )

    ui.separator().classes("mgmt-detail-separator")

    with ui.element("div").classes("classification-insight-shell w-full min-w-0"):
        ui.label("INDEL analysis").classes(
            "classification-insight-heading text-headline-small"
        )
        _add_clinvar_annotation_controls(
            launcher,
            sample_dir,
            compact=True,
            clinvar_status=clinvar_status,
        )

        indel_store = VariantDisplayStore.open_indel(clair3_dir, vcf_fallback=False)
        if indel_store is None:
            indel_host = ui.column().classes("w-full")
            with indel_host:
                ui.label("Loading INDEL variants…").classes("classification-insight-meta")

            def _mount_indel_fallback() -> None:
                indel_host.clear()
                fallback_store = VariantDisplayStore.open_indel(clair3_dir)
                if fallback_store is None:
                    with indel_host:
                        ui.label(
                            "INDEL VCF was not found. Run SNP analysis to generate INDEL output."
                        ).classes("classification-insight-meta")
                    return
                if fallback_store.total_variants == 0:
                    with indel_host:
                        ui.label("No INDEL variants were found.").classes(
                            "classification-insight-meta"
                        )
                    return

                def navigate_to_indel_region(indel_key: str) -> None:
                    try:
                        chrom, pos_str = indel_key.split(":", 1)
                        navigate_igv_to_snp(chrom, int(str(pos_str).replace(",", "")))
                    except (TypeError, ValueError):
                        logger.debug("Invalid INDEL key for IGV navigation: %s", indel_key)

                with indel_host:
                    _mount_paged_variant_table(
                        fallback_store,
                        details_title="INDEL details",
                        event_prefix="indel",
                        navigate_region=navigate_to_indel_region,
                    )

            schedule_after_page_sent(_mount_indel_fallback)
            return
        if indel_store.total_variants == 0:
            ui.label("No INDEL variants were found.").classes("classification-insight-meta")
            return

        def navigate_to_indel_region(indel_key: str) -> None:
            try:
                chrom, pos_str = indel_key.split(":", 1)
                navigate_igv_to_snp(chrom, int(str(pos_str).replace(",", "")))
            except (TypeError, ValueError):
                logger.debug("Invalid INDEL key for IGV navigation: %s", indel_key)

        _mount_paged_variant_table(
            indel_store,
            details_title="INDEL details",
            event_prefix="indel",
            navigate_region=navigate_to_indel_region,
        )
