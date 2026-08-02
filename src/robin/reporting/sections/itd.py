"""
ITD / insertion hotspot section for ROBIN PDF reports.
"""

from __future__ import annotations

import logging
import os
from typing import Optional

import pandas as pd
from reportlab.platypus import Paragraph, Spacer

from robin.analysis.itd_work import (
    ITD_EVENT_COLUMN_LABELS,
    ITD_EVENT_DISPLAY_COLUMNS,
    normalize_itd_events_df,
)
from .base import ReportSection

logger = logging.getLogger(__name__)


class ItdSection(ReportSection):
    """Section containing ITD / insertion hotspot calls when available."""

    def __init__(self, report):
        super().__init__(report)
        self.events: Optional[pd.DataFrame] = None
        self.summary: Optional[pd.DataFrame] = None
        self._analysis_available = False
        self._load_data()

    def _load_data(self) -> None:
        output = self.report.output
        events_path = os.path.join(output, "itd_events.csv")
        summary_path = os.path.join(output, "itd_summary.csv")

        try:
            if os.path.isfile(events_path):
                self._analysis_available = True
                self.events = normalize_itd_events_df(pd.read_csv(events_path))
            if os.path.isfile(summary_path):
                self._analysis_available = True
                summary = pd.read_csv(summary_path)
                preferred = [
                    "gene",
                    "label",
                    "chrom",
                    "start",
                    "end",
                    "n_events",
                    "max_support",
                    "max_vaf",
                ]
                cols = [c for c in preferred if c in summary.columns]
                self.summary = summary[cols] if cols else summary
        except Exception as exc:
            logger.error("Failed to load ITD report data: %s", exc, exc_info=True)
            self.events = None
            self.summary = None

    def _events_table_data(self, events: pd.DataFrame) -> list:
        headers = [
            ITD_EVENT_COLUMN_LABELS.get(col, col.replace("_", " ").title())
            for col in ITD_EVENT_DISPLAY_COLUMNS
            if col in events.columns
        ]
        fields = [col for col in ITD_EVENT_DISPLAY_COLUMNS if col in events.columns]
        rows = [headers]
        ordered = events
        sort_cols = [c for c in ("support", "vaf", "gene") if c in events.columns]
        if sort_cols:
            ordered = events.sort_values(sort_cols, ascending=[False] * len(sort_cols))

        for _, row in ordered.iterrows():
            cells = []
            for field in fields:
                value = row.get(field, "")
                if pd.isna(value):
                    cells.append("")
                elif field in {"vaf", "coverage"}:
                    try:
                        cells.append(f"{float(value):.4f}")
                    except (TypeError, ValueError):
                        cells.append(str(value))
                else:
                    cells.append(str(value))
            rows.append(cells)
        return rows

    def add_content(self):
        """Add the ITD / insertion section when analysis outputs are present."""
        logger.debug("Starting ITD section content generation")

        if not self._analysis_available:
            # Analysis not run for this sample — omit the section entirely.
            return

        self.elements.append(
            Paragraph("ITDs / Insertions", self.styles.styles["Heading1"])
        )
        self.elements.append(Spacer(1, 6))

        events = self.events if self.events is not None else pd.DataFrame()
        summary = self.summary if self.summary is not None else pd.DataFrame()
        n_events = len(events)
        n_called_genes = (
            int((summary["n_events"] > 0).sum())
            if not summary.empty and "n_events" in summary.columns
            else int(events["gene"].nunique())
            if n_events and "gene" in events.columns
            else 0
        )

        if n_events:
            summary_text = (
                f"Found {n_events} called ITD / insertion event(s) across "
                f"{n_called_genes} gene(s)."
            )
        else:
            summary_text = (
                "ITD / insertion analysis completed; no events passed "
                "length / support / VAF filters."
            )

        self.elements.append(Paragraph(summary_text, self.styles.styles["Normal"]))
        self.elements.append(Spacer(1, 4))

        self.summary_elements.append(
            Paragraph("ITDs / Insertions", self.styles.styles["Heading3"])
        )
        self.summary_elements.append(
            Paragraph(summary_text, self.styles.styles["Normal"])
        )

        note = (
            "Coverage is mean gene depth from target coverage when available "
            "(otherwise local hotspot depth) and is used as the VAF denominator."
        )
        self.elements.append(Paragraph(note, self.styles.styles["Normal"]))
        self.elements.append(Spacer(1, 6))

        if n_events:
            self.elements.append(
                Paragraph("Called events", self.styles.styles["Heading2"])
            )
            table = self.create_table(
                self._events_table_data(events),
                compact=True,
                font_size=7,
            )
            self.elements.append(table)
            self.elements.append(Spacer(1, 6))
            self.export_frames["itd_events"] = events.copy()
        else:
            self.elements.append(
                Paragraph(
                    "No insertions passed calling filters.",
                    self.styles.styles["Normal"],
                )
            )
            self.export_frames["itd_events"] = pd.DataFrame(
                columns=ITD_EVENT_DISPLAY_COLUMNS
            )

        if not summary.empty:
            called = (
                summary[summary["n_events"] > 0]
                if "n_events" in summary.columns
                else summary
            )
            self.export_frames["itd_summary"] = called.copy() if not called.empty else summary.head(0)
