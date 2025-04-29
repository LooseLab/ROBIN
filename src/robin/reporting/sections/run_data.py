"""
run_data.py

This module contains the run data summary section of the report.
"""

import logging
from datetime import datetime
from reportlab.platypus import Paragraph, Spacer, Table, TableStyle, PageBreak
from reportlab.lib.units import inch
from reportlab.lib.styles import ParagraphStyle
from .base import ReportSection

logger = logging.getLogger(__name__)


class RunDataSection(ReportSection):
    """Section containing the run data summary."""

    def _format_number(self, n):
        """Format number with commas for readability."""
        try:
            return "{:,}".format(int(float(n)))
        except (ValueError, TypeError):
            return str(n)

    def _format_timestamp(self, timestamp_str):
        """Convert timestamp string to readable format."""
        clean_ts = str(timestamp_str).strip("[]'\"")
        try:
            dt = datetime.strptime(clean_ts.split("+")[0], "%Y-%m-%dT%H:%M:%S.%f")
            return dt.strftime("%B %d, %Y %H:%M:%S")
        except Exception as e:
            logger.debug(f"Error formatting timestamp {timestamp_str}: {e}")
            return str(timestamp_str)

    def _convert_to_space_separated_string(self, value):
        """Convert various data types to a space-separated string."""
        if isinstance(value, (list, tuple)):
            return " ".join(map(str, value))
        elif isinstance(value, dict):
            return " ".join(f"{k}:{v}" for k, v in value.items())
        else:
            return str(value)

    def _create_info_table(self, data, title):
        """Create a formatted table for information display."""
        # Create header style
        header_style = ParagraphStyle(
            "HeaderStyle",
            parent=self.styles.styles["Normal"],
            fontName="Helvetica-Bold",
            fontSize=8,
            textColor=self.styles.COLORS["primary"],
            leading=10,
            spaceBefore=0,
            spaceAfter=0,
        )

        # Create cell style
        cell_style = ParagraphStyle(
            "CellStyle",
            parent=self.styles.styles["Normal"],
            fontSize=8,
            leading=10,
            spaceBefore=0,
            spaceAfter=0,
        )

        table_data = [[Paragraph(title, header_style)]]

        for key, value in data:
            table_data.append(
                [Paragraph(f"{key}:", cell_style), Paragraph(str(value), cell_style)]
            )

        table = Table(table_data, colWidths=[2 * inch, 4 * inch])
        table.setStyle(
            TableStyle(
                [
                    # Inherit modern table style
                    *self.MODERN_TABLE_STYLE._cmds,
                    # Preserve specific alignments
                    ("ALIGN", (0, 0), (0, -1), "LEFT"),  # Labels left-aligned
                    ("ALIGN", (1, 0), (1, -1), "LEFT"),  # Values left-aligned
                ]
            )
        )
        return table

    def add_content(self):
        """Add the run data content to the report."""
        logger.debug("Starting run data section content generation")

        # Add page break before detailed section
        self.elements.append(PageBreak())  # Add section header
        self.elements.append(
            Paragraph("Run Data Summary", self.styles.styles["Heading1"])
        )
        self.elements.append(Spacer(1, 12))

        # Get master data
        masterdf = self.report.masterdf
        if masterdf is None:
            logger.warning("No master data available")
            self.elements.append(
                Paragraph("No run data available", self.styles.styles["Normal"])
            )
            return

        try:
            # Debug logging
            logger.debug(f"Master data columns: {masterdf.columns.tolist()}")
            logger.debug(f"Master data shape: {masterdf.shape}")

            # Get the first row of data (since master.csv has all data in one row)
            master_data = masterdf.iloc[0]

            # Sample Information
            sample_info = [
                ("Sample ID", self.report.sample_id),
                ("Run Start", self._format_timestamp(master_data["run_time"])),
                ("Sequencing Device", master_data["devices"]),
                ("Flowcell ID", master_data["flowcell_ids"]),
                ("Basecalling Model", master_data["basecall_models"]),
            ]
            self.elements.append(
                self._create_info_table(sample_info, "Sample Information")
            )
            self.elements.append(Spacer(1, 12))

            # Run Statistics
            stats_info = [
                (
                    "Total Files",
                    self._format_number(
                        master_data["counter_bam_passed"]
                        + master_data["counter_bam_failed"]
                    ),
                ),
                (
                    "Passed Files",
                    self._format_number(master_data["counter_bam_passed"]),
                ),
                (
                    "Failed Files",
                    self._format_number(master_data["counter_bam_failed"]),
                ),
                (
                    "Total Bases",
                    self._format_number(master_data["counter_bases_count"]),
                ),
                (
                    "Mapped Reads",
                    self._format_number(master_data["counter_mapped_reads_num"]),
                ),
                (
                    "Unmapped Reads",
                    self._format_number(master_data["counter_unmapped_reads_num"]),
                ),
            ]
            self.elements.append(self._create_info_table(stats_info, "Run Statistics"))
            self.elements.append(Spacer(1, 12))

        except Exception as e:
            logger.error(f"Error generating run data section: {str(e)}")
            self.elements.append(
                Paragraph(
                    f"Error generating run data: {str(e)}",
                    self.styles.styles["Normal"],
                )
            )
