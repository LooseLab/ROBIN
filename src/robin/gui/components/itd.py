"""ITD / insertion hotspot results pane."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
from nicegui import ui

from robin.analysis.itd_work import (
    ITD_EVENT_COLUMN_LABELS,
    ITD_EVENT_DISPLAY_COLUMNS,
    normalize_itd_events_df,
)
from robin.gui.components.snp import navigate_igv_to_snp
from robin.gui.theme import styled_table, client_timer

logger = logging.getLogger(__name__)

_SUMMARY_PAGE_SIZE = 25
_EVENTS_PAGE_SIZE = 25

_SUMMARY_DISPLAY_COLUMNS = [
    "gene",
    "label",
    "chrom",
    "start",
    "end",
    "n_events",
    "max_support",
    "max_vaf",
]

_SUMMARY_COLUMN_LABELS = {
    "gene": "Gene",
    "label": "Label",
    "chrom": "Chrom",
    "start": "Start",
    "end": "End",
    "n_events": "Events",
    "max_support": "Max support",
    "max_vaf": "Max VAF",
}


def _load_events(sample_dir: Path) -> Optional[pd.DataFrame]:
    path = sample_dir / "itd_events.csv"
    if not path.is_file():
        return None
    try:
        return normalize_itd_events_df(pd.read_csv(path))
    except Exception as exc:
        logger.debug("Could not read %s: %s", path, exc)
        return None


def _load_summary(sample_dir: Path) -> Optional[pd.DataFrame]:
    path = sample_dir / "itd_summary.csv"
    if not path.is_file():
        return None
    try:
        df = pd.read_csv(path)
        cols = [c for c in _SUMMARY_DISPLAY_COLUMNS if c in df.columns]
        return df[cols] if cols else df
    except Exception as exc:
        logger.debug("Could not read %s: %s", path, exc)
        return None


def _columns_from_df(
    df: pd.DataFrame,
    *,
    labels: Optional[Dict[str, str]] = None,
    preferred: Optional[List[str]] = None,
    extra_columns: Optional[List[Dict[str, Any]]] = None,
) -> List[Dict[str, Any]]:
    label_map = labels or {}
    cols = [c for c in (preferred or []) if c in df.columns]
    if not cols:
        cols = list(df.columns)
    columns = [
        {
            "name": str(col),
            "label": label_map.get(str(col), str(col).replace("_", " ").title()),
            "field": str(col),
            "sortable": True,
            "align": "left",
        }
        for col in cols
    ]
    if extra_columns:
        columns.extend(extra_columns)
    return columns


def _format_event_rows(df: pd.DataFrame) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for idx, row in df.iterrows():
        record = row.to_dict()
        # Stable unique key across genes / positions for Quasar row_key.
        record["__row_key"] = (
            f"{record.get('gene', '')}|{record.get('chrom', '')}|"
            f"{record.get('position', '')}|{record.get('length', '')}|{idx}"
        )
        for key in ("vaf", "coverage"):
            if key in record and pd.notna(record[key]):
                try:
                    record[key] = round(float(record[key]), 4)
                except (TypeError, ValueError):
                    pass
        record["action"] = " "
        rows.append(record)
    return rows


def _add_table_search(table: Any, placeholder: str) -> None:
    try:
        with table.add_slot("top-right"):
            with ui.input(placeholder=placeholder).props(
                "type=search dense clearable"
            ).bind_value(table, "filter").add_slot("append"):
                ui.icon("search")
    except Exception:
        pass


def _igv_flank_for_event(length: Any, default: int = 150) -> int:
    """Pad IGV window so longer ITDs remain visible around the anchor."""
    try:
        event_len = int(float(length))
    except (TypeError, ValueError):
        return default
    return max(default, event_len + 50)


def navigate_igv_to_itd(chrom: str, position: int, length: Any = None) -> None:
    """Centre the embedded IGV browser on an ITD / insertion locus."""
    navigate_igv_to_snp(chrom, position, flank=_igv_flank_for_event(length))


def _wire_events_igv(table: Any, rows: List[Dict[str, Any]]) -> None:
    """Add View-in-IGV action column handlers for called ITD events."""
    row_by_key = {str(row.get("__row_key")): row for row in rows}

    table.add_slot(
        "body-cell-action",
        """
<q-td key="action" :props="props">
  <q-btn
    icon="visibility"
    size="sm"
    dense
    flat
    color="primary"
    @click="$parent.$emit('itd-view-igv', props.row.__row_key)"
    title="View in IGV"
  />
</q-td>
""",
    )

    def on_itd_view_igv(e: Any) -> None:
        try:
            row_key = getattr(e, "args", None)
            if row_key is None:
                return
            row = row_by_key.get(str(row_key))
            if not row:
                return
            chrom = str(row.get("chrom", "")).strip()
            pos_text = str(row.get("position", "")).replace(",", "").strip()
            if not chrom or not pos_text:
                return
            navigate_igv_to_itd(chrom, int(pos_text), row.get("length"))
        except Exception as ex:
            logger.debug("Error handling ITD IGV view: %s", ex)

    table.on("itd-view-igv", on_itd_view_igv)

    def on_row_click(e: Any) -> None:
        try:
            args = getattr(e, "args", None)
            row = None
            if isinstance(args, dict):
                row = args
            elif isinstance(args, (list, tuple)) and args:
                candidate = args[0]
                if isinstance(candidate, dict):
                    row = candidate
            if not row:
                return
            row_key = str(row.get("__row_key", ""))
            mapped = row_by_key.get(row_key, row)
            chrom = str(mapped.get("chrom", "")).strip()
            pos_text = str(mapped.get("position", "")).replace(",", "").strip()
            if chrom and pos_text:
                navigate_igv_to_itd(chrom, int(pos_text), mapped.get("length"))
        except Exception as ex:
            logger.debug("Error handling ITD row click: %s", ex)

    try:
        table.on("rowClick", on_row_click)
    except Exception:
        pass


def add_itd_section(launcher: Any, sample_dir: Path) -> None:
    """Render the ITDs pane for a sample directory (sample page or More details)."""
    sample_dir = Path(sample_dir)
    events_available = (sample_dir / "itd_events.csv").is_file() or (
        sample_dir / "itd_summary.csv"
    ).is_file()
    if not events_available:
        with ui.element("div").classes("classification-insight-shell w-full min-w-0"):
            ui.label("ITDs / insertions").classes(
                "classification-insight-heading text-headline-small"
            )
            ui.label(
                "No ITD results yet. Run ITD analysis for this sample to populate events."
            ).classes("classification-insight-meta")
        return

    with ui.element("div").classes("classification-insight-shell w-full min-w-0"):
        ui.label("ITDs / insertions").classes(
            "classification-insight-heading text-headline-small"
        )
        ui.label(
            "CIGAR insertion calls in curated hotspots and/or panel gene "
            "intervals (see workflow [itd] region_mode). Restricted to the "
            "active target panel. Coverage is mean gene depth when available "
            "(else local hotspot depth) and is the VAF denominator. "
            "Use View in IGV (or click a called-event row) to inspect the locus."
        ).classes("classification-insight-meta")

        status = ui.label("Loading…").classes("classification-insight-meta")
        with ui.row().classes("w-full items-center q-gutter-sm q-mt-sm"):
            show_empty = ui.checkbox("Show genes with no events").props("dense")
            show_empty.value = False
        table_host = ui.column().classes("w-full")

        def refresh() -> None:
            table_host.clear()
            summary = _load_summary(sample_dir)
            events = _load_events(sample_dir)

            if summary is None and events is None:
                status.set_text("No ITD results yet (waiting for analysis).")
                return

            n_events = 0 if events is None else len(events)
            n_summary = 0 if summary is None else len(summary)
            n_called_genes = (
                int((summary["n_events"] > 0).sum())
                if summary is not None
                and not summary.empty
                and "n_events" in summary.columns
                else n_events
            )
            qc_note = ""
            if (sample_dir / "itd_event_qc.csv").is_file():
                qc_note = " · read QC available (itd_event_qc.csv)"
            status.set_text(
                f"{n_events} called event(s) across {n_called_genes} gene(s)"
                + (f" ({n_summary} scanned)" if n_summary else "")
                + qc_note
            )

            with table_host:
                if summary is not None and not summary.empty:
                    ui.label("Gene summary").classes("text-subtitle2 q-mt-sm")
                    view = summary
                    if (
                        not show_empty.value
                        and "n_events" in summary.columns
                    ):
                        view = summary[summary["n_events"] > 0]
                    if view.empty:
                        ui.label(
                            "No genes with called events "
                            "(enable “Show genes with no events” to list all scanned genes)."
                        ).classes("text-body2 text-grey-7")
                    else:
                        rows = view.to_dict("records")
                        _, summary_table = styled_table(
                            columns=_columns_from_df(
                                view,
                                labels=_SUMMARY_COLUMN_LABELS,
                                preferred=_SUMMARY_DISPLAY_COLUMNS,
                            ),
                            rows=rows,
                            pagination=_SUMMARY_PAGE_SIZE,
                            class_size="table-xs",
                            row_key="gene",
                        )
                        _add_table_search(summary_table, "Search genes…")

                if events is not None and not events.empty:
                    ui.label("Called events").classes("text-subtitle2 q-mt-md")
                    # Highest support / VAF first for panel-scale tables.
                    sort_cols = [
                        c for c in ("support", "vaf", "gene") if c in events.columns
                    ]
                    ordered = (
                        events.sort_values(sort_cols, ascending=[False] * len(sort_cols))
                        if sort_cols
                        else events
                    )
                    event_rows = _format_event_rows(ordered)
                    _, events_table = styled_table(
                        columns=_columns_from_df(
                            ordered,
                            labels=ITD_EVENT_COLUMN_LABELS,
                            preferred=ITD_EVENT_DISPLAY_COLUMNS,
                            extra_columns=[
                                {
                                    "name": "action",
                                    "label": "View in IGV",
                                    "field": "action",
                                    "sortable": False,
                                    "align": "center",
                                }
                            ],
                        ),
                        rows=event_rows,
                        pagination=_EVENTS_PAGE_SIZE,
                        class_size="table-xs",
                        row_key="__row_key",
                    )
                    _add_table_search(events_table, "Search events…")
                    _wire_events_igv(events_table, event_rows)
                elif n_events == 0:
                    ui.label(
                        "No insertions passed length / support / VAF filters."
                    ).classes("text-body2 text-grey-7")

        show_empty.on("update:model-value", lambda _e: refresh())
        refresh()
        client_timer(30.0, refresh)
