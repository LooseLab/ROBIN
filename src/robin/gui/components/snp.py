from __future__ import annotations

from typing import Any, Dict, List, Optional
from pathlib import Path
import logging
import json
from robin.analysis.snp_processing import parse_vcf

try:
    from nicegui import ui
except ImportError:  # pragma: no cover
    ui = None

logger = logging.getLogger(__name__)

_TABLE_PREVIEW_THRESHOLD = 50_000
_TABLE_PREVIEW_LIMIT = 5_000


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
        
        ui.run_javascript(js_navigate, timeout=5.0)
        
    except Exception as e:
        logger.error(f"Error navigating IGV to SNP {chrom}:{pos}: {e}")


def add_snp_section(launcher: Any, sample_dir: Path) -> None:
    """
    Add SNP analysis section to the sample details page.
    
    Args:
        launcher: The GUI launcher instance
        sample_dir: Path to the sample directory
    """
    if not sample_dir or not sample_dir.exists():
        return
    
    # Look for VCF files in clair3 directory
    clair3_dir = sample_dir / "clair3"
    if not clair3_dir.exists():
        return
    
    # Look for snpsift_output.vcf (preferred) or other VCF files
    display_file = clair3_dir / "snpsift_output_display.json"

    if not display_file.exists():
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

    try:
        with display_file.open("r", encoding="utf-8") as f_in:
            snp_display = json.load(f_in)
    except Exception as exc:
        logger.error(f"Failed to load SNP display data: {exc}")
        with ui.element("div").classes("classification-insight-shell w-full min-w-0"):
            ui.label("SNP analysis").classes(
                "classification-insight-heading text-headline-small"
            )
            with ui.element("div").classes("classification-insight-card w-full min-w-0"):
                with ui.column().classes("w-full min-w-0 gap-2 p-2 md:p-3"):
                    ui.label(
                        "Could not load SNP variant data. Check logs for details."
                    ).classes("classification-insight-level classification-insight-level--low w-full")
        return

    columns: List[Dict[str, Any]] = snp_display.get("columns", [])
    rows_all: List[Dict[str, Any]] = snp_display.get("rows_all", [])
    summary: Dict[str, Any] = snp_display.get("summary", {})
    snp_regions_map: Dict[str, Dict[str, int]] = snp_display.get("snp_regions_map", {})

    total_variants = summary.get("total_variants", len(rows_all))
    pathogenic_count = summary.get("pathogenic_variants")
    if pathogenic_count is None:
        pathogenic_count = sum(
            1
            for row in rows_all
            if str(row.get("is_pathogenic", "")).strip().upper()
            in {"YES", "TRUE", "1", "PATHOGENIC"}
        )

    js_snp_regions_json = json.dumps(snp_regions_map)

    def navigate_to_snp_region(snp_key: str) -> None:
        if snp_key in snp_regions_map:
            snp_data = snp_regions_map[snp_key]
            navigate_igv_to_snp(snp_data["chrom"], snp_data["pos"])
            return
        # Fallback for non-pathogenic rows that are not in snp_regions_map.
        try:
            chrom, pos_str = snp_key.split(":", 1)
            navigate_igv_to_snp(chrom, int(str(pos_str).replace(",", "")))
        except (TypeError, ValueError):
            logger.debug(f"Invalid SNP key for IGV navigation: {snp_key}")

    js_init_snp_regions = f"""
        (function() {{
            window.snpRegionsMap = window.snpRegionsMap || {{}};
            Object.assign(window.snpRegionsMap, {js_snp_regions_json});
            console.log('[SNP] Initialized SNP regions map with', Object.keys(window.snpRegionsMap).length, 'pathogenic SNPs');
        }})();
    """
    ui.run_javascript(js_init_snp_regions, timeout=5.0)

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

    def _is_truthy(value: Any) -> bool:
        if isinstance(value, bool):
            return value
        return str(value).strip().upper() in {"YES", "TRUE", "1", "PATHOGENIC"}

    # Keep the SNP table readable by showing a concise default column set.
    preferred_display_fields = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "Gene_Name",
        "HGVS.p",
        "Annotation",
        "Annotation_Impact",
        "CLNSIG",
        "FILTER",
        "QUAL",
        "GT",
        "is_pathogenic",
        "details",
        "action",
    ]
    wide_fields = {"Gene_Name", "Annotation", "HGVS.p", "CLNSIG"}
    max_field_length = 80

    column_lookup = {
        col.get("field"): col for col in columns if isinstance(col, dict) and col.get("field")
    }

    visible_fields = [
        field for field in preferred_display_fields if field in column_lookup
    ]
    if not visible_fields:
        visible_fields = [
            col.get("field") for col in columns[:12] if isinstance(col, dict) and col.get("field")
        ]

    display_columns = [column_lookup[field].copy() for field in visible_fields]
    if "details" not in {col.get("field") for col in display_columns}:
        display_columns.append(
            {
                "name": "details",
                "label": "Details",
                "field": "details",
                "sortable": False,
            }
        )
    if "action" not in {col.get("field") for col in display_columns}:
        display_columns.append(
            {
                "name": "action",
                "label": "View in IGV",
                "field": "action",
                "sortable": False,
            }
        )

    def _compact_row(row: Dict[str, Any]) -> Dict[str, Any]:
        compact: Dict[str, Any] = {}
        for field in visible_fields:
            value = row.get(field, "")
            text = "" if value is None else str(value)
            if field in wide_fields and len(text) > max_field_length:
                text = f"{text[:max_field_length - 1]}..."
            compact[field] = text
        compact["details"] = " "
        compact["action"] = " "
        return compact

    snp_preview_mode = len(rows_all) > _TABLE_PREVIEW_THRESHOLD
    snp_rows_source = rows_all[:_TABLE_PREVIEW_LIMIT] if snp_preview_mode else rows_all

    page_state: Dict[str, Any] = {
        "filtered_indices": list(range(len(snp_rows_source))),
    }

    with ui.element("div").classes("classification-insight-shell w-full min-w-0"):
        ui.label("SNP analysis").classes(
            "classification-insight-heading text-headline-small"
        )

        with ui.row().classes("w-full gap-3 mb-3 flex-wrap items-baseline"):
            shown_total = len(snp_rows_source)
            ui.label(f"Total variants: {total_variants}").classes(
                "classification-insight-meta"
            )
            if snp_preview_mode:
                ui.label(
                    f"Preview mode: showing first {shown_total:,} rows (dataset too large). Use filters/search to narrow."
                ).classes("classification-insight-level classification-insight-level--low w-full")
            if pathogenic_count > 0:
                ui.label(f"Pathogenic variants: {pathogenic_count}").classes(
                    "classification-insight-level classification-insight-level--low w-auto"
                )

        from robin.gui.theme import (
            clamp_qtable_server_pagination,
            styled_server_paged_table,
            wire_qtable_server_pagination_handlers,
        )

        snp_has_dp = "DP" in column_lookup

        with ui.row().classes("w-full gap-2 mb-2 flex-wrap items-end"):
            snp_pass_only = ui.checkbox("PASS only").props("dense")
            snp_pathogenic_only = ui.checkbox("Pathogenic only").props("dense")
            snp_min_qual = ui.number("Min QUAL", value=None).props(
                "dense outlined clearable"
            ).classes("w-32")
            if snp_has_dp:
                snp_min_dp = ui.number("Min DP", value=None).props(
                    "dense outlined clearable"
                ).classes("w-32")
            else:
                snp_min_dp = None
            snp_search = ui.input("Search (gene/variant)").props(
                "dense outlined clearable"
            ).classes("w-64")
            snp_reset_button = ui.button("Reset").props("dense no-caps")

        _snp_total_filtered = len(page_state["filtered_indices"])
        snp_init_pagination = clamp_qtable_server_pagination(
            {
                "sortBy": None,
                "descending": False,
                "page": 1,
                "rowsPerPage": 100,
                "rowsNumber": _snp_total_filtered,
            },
            rows_number=_snp_total_filtered,
            rows_per_page_default=100,
        )
        table_container, snp_table = styled_server_paged_table(
            columns=display_columns,
            rows=[],
            pagination=snp_init_pagination,
            row_key="__row_idx",
            class_size="table-xs",
        )
        snp_filtered_count_label = ui.label(
            f"{_snp_total_filtered} variants match filters ({len(snp_rows_source)} loaded)"
        ).classes("classification-insight-meta")

        def _fill_snp_from_pagination(pag: Dict[str, Any]) -> None:
            total_filtered = len(page_state["filtered_indices"])
            pag = clamp_qtable_server_pagination(
                pag,
                rows_number=total_filtered,
                rows_per_page_default=100,
            )
            rpp = int(pag["rowsPerPage"])
            page = int(pag["page"])
            start = (page - 1) * rpp
            end = start + rpp
            filtered_indices = page_state["filtered_indices"]
            rows_out: List[Dict[str, Any]] = []
            for idx in filtered_indices[start:end]:
                rows_out.append(_compact_row(snp_rows_source[idx]))
                rows_out[-1]["__row_idx"] = idx
            snp_table.rows = rows_out
            snp_table.pagination = pag
            snp_filtered_count_label.text = (
                f"{total_filtered} variants match filters ({len(snp_rows_source)} loaded)"
                if total_filtered
                else f"0 variants match filters ({len(snp_rows_source)} loaded)"
            )
            snp_table.update()

        wire_qtable_server_pagination_handlers(snp_table, _fill_snp_from_pagination)

        def _apply_snp_filters() -> None:
            pass_only = bool(getattr(snp_pass_only, "value", False))
            pathogenic_only = bool(getattr(snp_pathogenic_only, "value", False))
            min_qual = _to_float(getattr(snp_min_qual, "value", None))
            min_dp = _to_float(getattr(snp_min_dp, "value", None)) if snp_has_dp else None
            search_text = str(getattr(snp_search, "value", "") or "").strip().lower()

            filtered_indices: List[int] = []
            for idx, full_row in enumerate(snp_rows_source):

                if pass_only and str(full_row.get("FILTER", "")).strip().upper() != "PASS":
                    continue
                if pathogenic_only and not _is_truthy(full_row.get("is_pathogenic", "")):
                    continue

                qual = _to_float(full_row.get("QUAL"))
                if min_qual is not None and (qual is None or qual < min_qual):
                    continue

                if snp_has_dp and min_dp is not None:
                    dp = _to_float(full_row.get("DP"))
                    if dp is None or dp < min_dp:
                        continue

                if search_text:
                    haystack_parts = []
                    for field in visible_fields:
                        value = full_row.get(field, "")
                        if value is not None:
                            haystack_parts.append(str(value))
                    haystack = " ".join(haystack_parts).lower()
                    if search_text not in haystack:
                        continue

                filtered_indices.append(idx)

            page_state["filtered_indices"] = filtered_indices
            pag = clamp_qtable_server_pagination(
                dict(snp_table.pagination),
                rows_number=len(filtered_indices),
                rows_per_page_default=100,
            )
            pag["page"] = 1
            _fill_snp_from_pagination(pag)

        snp_pass_only.on("update:model-value", lambda _e: _apply_snp_filters())
        snp_pathogenic_only.on("update:model-value", lambda _e: _apply_snp_filters())
        snp_min_qual.on("update:model-value", lambda _e: _apply_snp_filters())
        if snp_has_dp and snp_min_dp is not None:
            snp_min_dp.on("update:model-value", lambda _e: _apply_snp_filters())
        snp_search.on("update:model-value", lambda _e: _apply_snp_filters())
        snp_reset_button.on_click(
            lambda: (
                setattr(snp_pass_only, "value", False),
                setattr(snp_pathogenic_only, "value", False),
                setattr(snp_min_qual, "value", None),
                snp_has_dp and setattr(snp_min_dp, "value", None),
                setattr(snp_search, "value", ""),
                _apply_snp_filters(),
            )
        )

        if any(col.get("field") in {"action", "details"} for col in display_columns):
            try:
                with ui.dialog() as details_dialog, ui.card().classes(
                    "robin-dialog-surface w-[95vw] max-w-6xl max-h-[85vh] overflow-auto "
                    "p-4 md:p-5"
                ):
                    ui.label("Variant details").classes(
                        "classification-insight-heading text-headline-small"
                    )
                    ui.separator().classes("mgmt-detail-separator")
                    details_container = ui.column().classes("w-full gap-2")
                    with ui.row().classes("w-full justify-end pt-2"):
                        ui.button("Close", on_click=details_dialog.close).props(
                            "color=primary no-caps outline"
                        )

                def show_variant_details(row_idx: int) -> None:
                    if row_idx < 0 or row_idx >= len(snp_rows_source):
                        ui.notify("Variant details not found.", type="warning")
                        return
                    row_data = snp_rows_source[row_idx]
                    if not row_data:
                        ui.notify("Variant details not found.", type="warning")
                        return

                    details_container.clear()
                    detail_fields = [
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
                        "CLNSIG",
                        "FILTER",
                        "QUAL",
                        "GT",
                    ]
                    ui_only_fields = {"action", "details", "__row_id"}
                    ordered_fields = detail_fields + sorted(
                        [
                            k
                            for k in row_data.keys()
                            if k not in detail_fields and k not in ui_only_fields
                        ]
                    )
                    with details_container:
                        for field in ordered_fields:
                            value = row_data.get(field, "")
                            if value is None or str(value) == "":
                                continue
                            with ui.row().classes("w-full items-start gap-2"):
                                ui.label(f"{field}:").classes(
                                    "text-xs font-semibold min-w-[180px]"
                                )
                                ui.label(str(value)).classes(
                                    "text-xs whitespace-pre-wrap break-all flex-1"
                                )
                    details_dialog.open()

                snp_table.add_slot(
                    "body-cell-details",
                    """
<q-td key="details" :props="props">
  <q-btn
    icon="description"
    size="sm"
    dense
    flat
    color="secondary"
    @click="$parent.$emit('snp-show-details', props.row.__row_idx)"
    title="Show full details"
  />
</q-td>
""",
                )

                snp_table.add_slot(
                    "body-cell-action",
                    """
<q-td key="action" :props="props">
  <q-btn 
    icon="visibility" 
    size="sm" 
    dense 
    flat 
    color="primary"
    @click="$parent.$emit('snp-view-igv', props.row.__row_idx)"
    title="View in IGV"
  />
</q-td>
""",
                )

                def on_snp_view_igv(e):
                    try:
                        row_idx = getattr(e, "args", None)
                        if row_idx is None:
                            return
                        row_idx = int(row_idx)
                        if row_idx < 0 or row_idx >= len(snp_rows_source):
                            return

                        row = snp_rows_source[row_idx]
                        chrom = str(row.get("CHROM", "")).strip()
                        pos_text = str(row.get("POS", "")).replace(",", "").strip()
                        if chrom and pos_text:
                            snp_key = f"{chrom}:{pos_text}"
                            logging.debug(f"[SNP] Button clicked for: {snp_key}")
                            navigate_to_snp_region(snp_key)
                    except Exception as ex:
                        logger.debug(f"Error handling SNP IGV view: {ex}")

                def on_snp_show_details(e):
                    try:
                        row_idx = getattr(e, "args", None)
                        if row_idx is not None:
                            show_variant_details(int(row_idx))
                    except Exception as ex:
                        logger.debug(f"Error handling SNP details view: {ex}")

                snp_table.on("snp-view-igv", on_snp_view_igv)
                snp_table.on("snp-show-details", on_snp_show_details)
            except Exception as ex:
                logger.warning(f"Could not add action button slot: {ex}")

        for col in snp_table.columns:
            col["sortable"] = False
        _fill_snp_from_pagination(snp_init_pagination)
        try:
            def _cleanup_snp_page() -> None:
                page_state["filtered_indices"] = []
                snp_table.rows = []
                snp_rows_source.clear()
            ui.context.client.on_disconnect(_cleanup_snp_page)
        except Exception:
            pass

    ui.separator().classes("mgmt-detail-separator")

    # INDEL table for More Details page.
    indel_vcf = clair3_dir / "snpsift_indel_output.vcf"
    with ui.element("div").classes("classification-insight-shell w-full min-w-0"):
        ui.label("INDEL analysis").classes(
            "classification-insight-heading text-headline-small"
        )

        if not indel_vcf.exists():
            ui.label(
                "INDEL VCF was not found. Run SNP analysis to generate INDEL output."
            ).classes("classification-insight-meta")
            return

        indel_df = parse_vcf(indel_vcf)
        if indel_df is None:
            ui.label(
                "Could not parse INDEL VCF data. Check logs for details."
            ).classes("classification-insight-level classification-insight-level--low w-full")
            return

        if indel_df.empty:
            ui.label("No INDEL variants were found.").classes(
                "classification-insight-meta"
            )
            return

        indel_column_lookup = {
            str(field): {
                "name": str(field),
                "label": str(field).replace("_", " "),
                "field": str(field),
                "sortable": True,
            }
            for field in indel_df.columns
        }
        indel_visible_fields = [
            field for field in preferred_display_fields if field in indel_column_lookup
        ]
        if not indel_visible_fields:
            indel_visible_fields = list(indel_column_lookup.keys())[:12]

        indel_display_columns = [
            indel_column_lookup[field].copy() for field in indel_visible_fields
        ]
        if "details" not in {col.get("field") for col in indel_display_columns}:
            indel_display_columns.append(
                {
                    "name": "details",
                    "label": "Details",
                    "field": "details",
                    "sortable": False,
                }
            )
        if "action" not in {col.get("field") for col in indel_display_columns}:
            indel_display_columns.append(
                {
                    "name": "action",
                    "label": "View in IGV",
                    "field": "action",
                    "sortable": False,
                }
            )

        def _indel_cell_to_text(value: Any) -> str:
            if value is None:
                return ""
            if isinstance(value, bool):
                return "Yes" if value else "No"
            return str(value)

        def _indel_row_text_map(idx: int) -> Dict[str, str]:
            row_data: Dict[str, str] = {}
            row = indel_df.iloc[idx]
            for field in indel_df.columns:
                row_data[str(field)] = _indel_cell_to_text(row.get(field))
            # Normalize boolean display consistency for filtering and details.
            row_data["is_pathogenic"] = (
                "Yes"
                if _is_truthy(row_data.get("is_pathogenic", ""))
                else "No"
            )
            return row_data

        def _compact_indel_row(idx: int) -> Dict[str, Any]:
            row_map = _indel_row_text_map(idx)
            compact: Dict[str, Any] = {}
            for field in indel_visible_fields:
                value = row_map.get(field, "")
                text = "" if value is None else str(value)
                if field in wide_fields and len(text) > max_field_length:
                    text = f"{text[:max_field_length - 1]}..."
                compact[field] = text
            compact["details"] = " "
            compact["action"] = " "
            compact["__row_idx"] = idx
            return compact

        total_indel_rows_all = int(len(indel_df))
        indel_preview_mode = total_indel_rows_all > _TABLE_PREVIEW_THRESHOLD
        total_indel_rows = (
            _TABLE_PREVIEW_LIMIT if indel_preview_mode else total_indel_rows_all
        )
        pathogenic_indel_count = 0
        if "is_pathogenic" in indel_df.columns:
            try:
                pathogenic_indel_count = int(
                    indel_df["is_pathogenic"]
                    .astype(str)
                    .str.strip()
                    .str.upper()
                    .isin({"YES", "TRUE", "1", "PATHOGENIC"})
                    .sum()
                )
            except Exception:
                pathogenic_indel_count = 0

        with ui.row().classes("w-full gap-3 mb-3 flex-wrap items-baseline"):
            ui.label(f"Total variants: {total_indel_rows_all}").classes(
                "classification-insight-meta"
            )
            if indel_preview_mode:
                ui.label(
                    f"Preview mode: showing first {total_indel_rows:,} rows (dataset too large). Use filters/search to narrow."
                ).classes("classification-insight-level classification-insight-level--low w-full")
            if pathogenic_indel_count:
                ui.label(
                    f"Pathogenic variants: {pathogenic_indel_count}"
                ).classes(
                    "classification-insight-level classification-insight-level--low w-auto"
                )

        from robin.gui.theme import (
            clamp_qtable_server_pagination,
            styled_server_paged_table,
            wire_qtable_server_pagination_handlers,
        )

        indel_has_dp = "DP" in indel_column_lookup

        with ui.row().classes("w-full gap-2 mb-2 flex-wrap items-end"):
            indel_pass_only = ui.checkbox("PASS only").props("dense")
            indel_pathogenic_only = ui.checkbox("Pathogenic only").props("dense")
            indel_min_qual = ui.number("Min QUAL", value=None).props(
                "dense outlined clearable"
            ).classes("w-32")
            if indel_has_dp:
                indel_min_dp = ui.number("Min DP", value=None).props(
                    "dense outlined clearable"
                ).classes("w-32")
            else:
                indel_min_dp = None
            indel_search = ui.input("Search (gene/variant)").props(
                "dense outlined clearable"
            ).classes("w-64")
            indel_reset_button = ui.button("Reset").props("dense no-caps")

        indel_page_state: Dict[str, Any] = {
            "filtered_indices": list(range(total_indel_rows)),
        }
        _indel_tf0 = len(indel_page_state["filtered_indices"])
        indel_init_pagination = clamp_qtable_server_pagination(
            {
                "sortBy": None,
                "descending": False,
                "page": 1,
                "rowsPerPage": 100,
                "rowsNumber": _indel_tf0,
            },
            rows_number=_indel_tf0,
            rows_per_page_default=100,
        )
        _, indel_table = styled_server_paged_table(
            columns=indel_display_columns,
            rows=[],
            pagination=indel_init_pagination,
            row_key="__row_idx",
            class_size="table-xs",
        )
        indel_filtered_count_label = ui.label(
            f"{_indel_tf0} variants match filters ({total_indel_rows} loaded)"
        ).classes("classification-insight-meta")

        def _fill_indel_from_pagination(pag: Dict[str, Any]) -> None:
            total_filtered = len(indel_page_state["filtered_indices"])
            pag = clamp_qtable_server_pagination(
                pag,
                rows_number=total_filtered,
                rows_per_page_default=100,
            )
            rpp = int(pag["rowsPerPage"])
            page = int(pag["page"])
            start = (page - 1) * rpp
            end = start + rpp
            filtered_indices = indel_page_state["filtered_indices"]
            indel_table.rows = [
                _compact_indel_row(idx) for idx in filtered_indices[start:end]
            ]
            indel_table.pagination = pag
            indel_filtered_count_label.text = (
                f"{total_filtered} variants match filters ({total_indel_rows} loaded)"
                if total_filtered
                else f"0 variants match filters ({total_indel_rows} loaded)"
            )
            indel_table.update()

        wire_qtable_server_pagination_handlers(indel_table, _fill_indel_from_pagination)

        def _apply_indel_filters() -> None:
            pass_only = bool(getattr(indel_pass_only, "value", False))
            pathogenic_only = bool(getattr(indel_pathogenic_only, "value", False))
            min_qual = _to_float(getattr(indel_min_qual, "value", None))
            min_dp = (
                _to_float(getattr(indel_min_dp, "value", None)) if indel_has_dp else None
            )
            search_text = str(getattr(indel_search, "value", "") or "").strip().lower()

            filtered_indices: List[int] = []
            for idx in range(total_indel_rows):
                full_row = _indel_row_text_map(idx)

                if pass_only and str(full_row.get("FILTER", "")).strip().upper() != "PASS":
                    continue
                if pathogenic_only and not _is_truthy(full_row.get("is_pathogenic", "")):
                    continue

                qual = _to_float(full_row.get("QUAL"))
                if min_qual is not None and (qual is None or qual < min_qual):
                    continue

                if indel_has_dp and min_dp is not None:
                    dp = _to_float(full_row.get("DP"))
                    if dp is None or dp < min_dp:
                        continue

                if search_text:
                    haystack_parts = []
                    for field in indel_visible_fields:
                        value = full_row.get(field, "")
                        if value is not None:
                            haystack_parts.append(str(value))
                    haystack = " ".join(haystack_parts).lower()
                    if search_text not in haystack:
                        continue

                filtered_indices.append(idx)

            indel_page_state["filtered_indices"] = filtered_indices
            pag = clamp_qtable_server_pagination(
                dict(indel_table.pagination),
                rows_number=len(filtered_indices),
                rows_per_page_default=100,
            )
            pag["page"] = 1
            _fill_indel_from_pagination(pag)

        indel_pass_only.on("update:model-value", lambda _e: _apply_indel_filters())
        indel_pathogenic_only.on("update:model-value", lambda _e: _apply_indel_filters())
        indel_min_qual.on("update:model-value", lambda _e: _apply_indel_filters())
        if indel_has_dp and indel_min_dp is not None:
            indel_min_dp.on("update:model-value", lambda _e: _apply_indel_filters())
        indel_search.on("update:model-value", lambda _e: _apply_indel_filters())
        indel_reset_button.on_click(
            lambda: (
                setattr(indel_pass_only, "value", False),
                setattr(indel_pathogenic_only, "value", False),
                setattr(indel_min_qual, "value", None),
                indel_has_dp and setattr(indel_min_dp, "value", None),
                setattr(indel_search, "value", ""),
                _apply_indel_filters(),
            )
        )

        with ui.dialog() as indel_details_dialog, ui.card().classes(
            "robin-dialog-surface w-[95vw] max-w-6xl max-h-[85vh] overflow-auto "
            "p-4 md:p-5"
        ):
            ui.label("INDEL details").classes(
                "classification-insight-heading text-headline-small"
            )
            ui.separator().classes("mgmt-detail-separator")
            indel_details_container = ui.column().classes("w-full gap-2")
            with ui.row().classes("w-full justify-end pt-2"):
                ui.button("Close", on_click=indel_details_dialog.close).props(
                    "color=primary no-caps outline"
                )

        def show_indel_details(row_idx: int) -> None:
            if row_idx < 0 or row_idx >= total_indel_rows:
                ui.notify("INDEL details not found.", type="warning")
                return
            row_data = _indel_row_text_map(row_idx)
            indel_details_container.clear()
            detail_fields = [
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
                "CLNSIG",
                "FILTER",
                "QUAL",
                "GT",
            ]
            ui_only_fields = {"action", "details", "__row_id", "__row_idx"}
            ordered_fields = detail_fields + sorted(
                [
                    k
                    for k in row_data.keys()
                    if k not in detail_fields and k not in ui_only_fields
                ]
            )
            with indel_details_container:
                for field in ordered_fields:
                    value = row_data.get(field, "")
                    if value in (None, ""):
                        continue
                    with ui.row().classes("w-full items-start gap-2"):
                        ui.label(f"{field}:").classes("text-xs font-semibold min-w-[180px]")
                        ui.label(str(value)).classes("text-xs whitespace-pre-wrap break-all flex-1")
            indel_details_dialog.open()

        indel_table.add_slot(
            "body-cell-details",
            """
<q-td key="details" :props="props">
  <q-btn
    icon="description"
    size="sm"
    dense
    flat
    color="secondary"
    @click="$parent.$emit('indel-show-details', props.row.__row_idx)"
    title="Show full details"
  />
</q-td>
""",
        )

        indel_table.add_slot(
            "body-cell-action",
            """
<q-td key="action" :props="props">
  <q-btn
    icon="visibility"
    size="sm"
    dense
    flat
    color="primary"
    @click="$parent.$emit('indel-view-igv', props.row.__row_idx)"
    title="View in IGV"
  />
</q-td>
""",
        )

        def navigate_to_indel_region(indel_key: str) -> None:
            # Direct `chrom:pos` navigation for all rows.
            try:
                chrom, pos_str = indel_key.split(":", 1)
                navigate_igv_to_snp(chrom, int(str(pos_str).replace(",", "")))
            except (TypeError, ValueError):
                logger.debug(f"Invalid INDEL key for IGV navigation: {indel_key}")

        def on_indel_view_igv(e):
            try:
                row_idx = getattr(e, "args", None)
                if row_idx is None:
                    return
                row_idx = int(row_idx)
                if row_idx < 0 or row_idx >= total_indel_rows:
                    return
                row = _indel_row_text_map(row_idx)
                chrom = str(row.get("CHROM", "")).strip()
                pos_text = str(row.get("POS", "")).replace(",", "").strip()
                if not chrom or not pos_text:
                    return
                indel_key = f"{chrom}:{pos_text}"
                navigate_to_indel_region(indel_key)
            except Exception as ex:
                logger.debug(f"Error handling INDEL IGV view: {ex}")

        def on_indel_show_details(e):
            try:
                row_idx = getattr(e, "args", None)
                if row_idx is not None:
                    show_indel_details(int(row_idx))
            except Exception as ex:
                logger.debug(f"Error handling INDEL details view: {ex}")

        indel_table.on("indel-view-igv", on_indel_view_igv)
        indel_table.on("indel-show-details", on_indel_show_details)

        for col in indel_table.columns:
            col["sortable"] = False
        _fill_indel_from_pagination(indel_init_pagination)
        try:
            def _cleanup_indel_page() -> None:
                indel_page_state["filtered_indices"] = []
                indel_table.rows = []
            ui.context.client.on_disconnect(_cleanup_indel_page)
        except Exception:
            pass
