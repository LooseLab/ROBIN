"""ITD / insertion hotspot results pane."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
from nicegui import ui

from robin.gui.theme import styled_table, client_timer

logger = logging.getLogger(__name__)

_SUMMARY_PAGE_SIZE = 25
_EVENTS_PAGE_SIZE = 25


def _load_events(sample_dir: Path) -> Optional[pd.DataFrame]:
    path = sample_dir / "itd_events.csv"
    if not path.is_file():
        return None
    try:
        return pd.read_csv(path)
    except Exception as exc:
        logger.debug("Could not read %s: %s", path, exc)
        return None


def _load_summary(sample_dir: Path) -> Optional[pd.DataFrame]:
    path = sample_dir / "itd_summary.csv"
    if not path.is_file():
        return None
    try:
        return pd.read_csv(path)
    except Exception as exc:
        logger.debug("Could not read %s: %s", path, exc)
        return None


def _columns_from_df(df: pd.DataFrame) -> List[Dict[str, Any]]:
    return [
        {
            "name": str(col),
            "label": str(col),
            "field": str(col),
            "sortable": True,
            "align": "left",
        }
        for col in df.columns
    ]


def _format_event_rows(df: pd.DataFrame) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for idx, row in df.iterrows():
        record = row.to_dict()
        # Stable unique key across genes / positions for Quasar row_key.
        record["__row_key"] = (
            f"{record.get('gene', '')}|{record.get('chrom', '')}|"
            f"{record.get('position', '')}|{record.get('length', '')}|{idx}"
        )
        for key in ("vaf", "gene_coverage", "coverage"):
            if key in record and pd.notna(record[key]):
                try:
                    record[key] = round(float(record[key]), 4)
                except (TypeError, ValueError):
                    pass
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


def add_itd_section(launcher: Any, sample_dir: Path) -> None:
    """Render the ITDs pane for a sample directory."""
    sample_dir = Path(sample_dir)
    with ui.card().classes("w-full"):
        ui.label("ITDs / insertions").classes("text-h6")
        ui.label(
            "CIGAR insertion calls in curated hotspots and/or panel gene "
            "intervals (see workflow [itd] region_mode). Restricted to the "
            "active target panel."
        ).classes("text-caption text-grey-7")

        status = ui.label("Loading…").classes("text-body2")
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
                            columns=_columns_from_df(view),
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
                    _, events_table = styled_table(
                        columns=_columns_from_df(ordered),
                        rows=_format_event_rows(ordered),
                        pagination=_EVENTS_PAGE_SIZE,
                        class_size="table-xs",
                        row_key="__row_key",
                    )
                    _add_table_search(events_table, "Search events…")
                elif n_events == 0:
                    ui.label(
                        "No insertions passed length / support / VAF filters."
                    ).classes("text-body2 text-grey-7")

        show_empty.on("update:model-value", lambda _e: refresh())
        refresh()
        client_timer(30.0, refresh)
