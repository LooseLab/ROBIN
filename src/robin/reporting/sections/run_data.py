"""
run_data.py

This module contains the run data summary section of the report.
"""

import logging
from datetime import datetime
from reportlab.platypus import Paragraph, Spacer, PageBreak
from .base import ReportSection
from robin.reporting.components import metadata_table, section_heading, warning_box
from robin.analysis.bam_preprocessor import (
    _get_modbase_model_warning,
    _get_modbase_model_warning_level,
    _is_unresolved_modbase_model,
)

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

    @staticmethod
    def _master_value(master_data, key):
        """Return a displayable master.csv value, including for older files."""
        if key not in master_data.index:
            return "Missing"
        value = master_data.get(key)
        if value is None:
            return "Missing"
        try:
            if value != value:
                return "Missing"
        except Exception:
            pass
        text = str(value).strip()
        return text if text else "Missing"

    def _create_info_table(self, data, title):
        """Create a formatted metadata table using the shared component."""
        return metadata_table(data, self.styles, title=title, width=self.report.doc.width)

    def add_content(self):
        """Add the run data content to the report."""
        logger.debug("Starting run data section content generation")

        # Add page break before detailed section
        self.elements.append(PageBreak())  # Add section header
        self.elements.extend(section_heading("Run Data Summary", self.styles))
        self.elements.append(Spacer(1, 6))

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
            modbase_model = self._master_value(master_data, "modbase_models")
            modbase_model_display = (
                "Unknown"
                if _is_unresolved_modbase_model(modbase_model)
                else modbase_model
            )
            modbase_warning = _get_modbase_model_warning(
                None if modbase_model == "Missing" else modbase_model
            )

            if modbase_warning:
                warning_level = _get_modbase_model_warning_level(
                    None if modbase_model == "Missing" else modbase_model
                )
                label = (
                    "Methylation model note"
                    if warning_level == "info"
                    else "Methylation model warning"
                )
                warning_text = f"<b>{label}:</b> {modbase_warning}"
                style_name = "Normal" if warning_level == "info" else "Warning"
                self.summary_elements.extend(warning_box(warning_text, self.styles, markup=True))
                self.elements.extend(warning_box(warning_text, self.styles, markup=True))

            # Sample Information
            sample_info = [
                ("Sample ID", self.report.sample_id),
                ("Run Start", self._format_timestamp(master_data["run_time"])),
                ("Sequencing Device", master_data["devices"]),
                ("Flowcell ID", master_data["flowcell_ids"]),
                (
                    "Basecalling Model",
                    self._master_value(master_data, "basecall_models"),
                ),
                ("Modbase Model", modbase_model_display),
            ]
            self.elements.append(
                self._create_info_table(sample_info, "Sample Information")
            )
            self.elements.append(Spacer(1, 6))

            # Sample identifiers (if decrypted and provided)
            if getattr(self.report, "sample_identifiers", None):
                si = self.report.sample_identifiers
                sample_info = [
                    ("Sample ID", si.get("sample_id", "") or self.report.sample_id or "—"),
                    ("Test ID", si.get("test_id", "") or "—"),
                    ("First name", si.get("first_name", "") or "—"),
                    ("Last name", si.get("last_name", "") or "—"),
                    ("Date of birth", si.get("dob", "") or "—"),
                ]
                if si.get("nhs_number"):
                    sample_info.append(("Hospital Number", si["nhs_number"]))
                if si.get("notes"):
                    sample_info.append(("Notes", si["notes"]))
                self.elements.append(
                    self._create_info_table(sample_info, "Sample Identifiers")
                )
                self.elements.append(Spacer(1, 6))

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
            self.elements.append(Spacer(1, 6))

        except Exception as e:
            logger.error(f"Error generating run data section: {str(e)}")
            self.elements.append(
                Paragraph(
                    f"Error generating run data: {str(e)}",
                    self.styles.styles["Normal"],
                )
            )
