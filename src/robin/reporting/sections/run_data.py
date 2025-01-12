"""
run_data.py

This module contains the run data summary section of the report.
"""

import os
import logging
import pickle
from datetime import datetime
from reportlab.platypus import Paragraph, Spacer, Table, TableStyle
from reportlab.lib import colors
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
            'HeaderStyle',
            parent=self.styles.styles["Normal"],
            fontName='Helvetica-Bold',
            fontSize=8,
            textColor=self.styles.COLORS["primary"],
            leading=10,
            spaceBefore=0,
            spaceAfter=0
        )
        
        # Create cell style
        cell_style = ParagraphStyle(
            'CellStyle',
            parent=self.styles.styles["Normal"],
            fontSize=8,
            leading=10,
            spaceBefore=0,
            spaceAfter=0
        )
        
        table_data = [[Paragraph(title, header_style)]]
        
        for key, value in data:
            table_data.append([
                Paragraph(f"{key}:", cell_style),
                Paragraph(str(value), cell_style)
            ])

        table = Table(table_data, colWidths=[2*inch, 4*inch])
        table.setStyle(TableStyle([
            # Inherit modern table style
            *self.MODERN_TABLE_STYLE._cmds,
            
            # Preserve specific alignments
            ('ALIGN', (0, 0), (0, -1), 'LEFT'),    # Labels left-aligned
            ('ALIGN', (1, 0), (1, -1), 'LEFT'),    # Values left-aligned
        ]))
        return table

    def add_content(self):
        """Add the run data content to the report."""
        logger.debug("Starting run data section content generation")
        
        # Add section header
        self.elements.append(Paragraph("Run Data Summary", self.styles.styles["Heading1"]))
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
            # Get sample-specific data
            sample_id = self.report.sample_id
            masterdf_dict = eval(masterdf[masterdf.index == "samples"][1]["samples"])[sample_id]

            # Sample Information
            sample_info = [
                ("Sample ID", sample_id),
                ("Run Start", self._format_timestamp(masterdf_dict['run_time'])),
                ("Target Panel", ' '.join(masterdf.loc[(masterdf.index == 'target_panel')][1].values))
            ]
            self.elements.append(self._create_info_table(sample_info, "Sample Information"))
            self.elements.append(Spacer(1, 12))

            # Device Details
            device_info = [
                ("Sequencing Device", self._convert_to_space_separated_string(masterdf_dict['devices'])),
                ("Flowcell ID", self._convert_to_space_separated_string(masterdf_dict['flowcell_ids'])),
                ("Basecalling Model", self._convert_to_space_separated_string(masterdf_dict['basecall_models']))
            ]
            self.elements.append(self._create_info_table(device_info, "Device Details"))
            self.elements.append(Spacer(1, 12))

            # File Locations
            file_info = [
                ("Run Directory", ' '.join(masterdf.loc[(masterdf.index == 'watchfolder')][1].values)),
                ("Output Directory", ' '.join(masterdf.loc[(masterdf.index == 'output')][1].values)),
                ("Reference", ' '.join(masterdf.loc[(masterdf.index == 'reference')][1].values))
            ]
            self.elements.append(self._create_info_table(file_info, "File Locations"))
            self.elements.append(Spacer(1, 12))

            # File Sources
            file_sources = [
                ("Master Data", os.path.join(self.report.output, "master.csv")),
            ]
            self.elements.append(self._create_info_table(file_sources, "File Sources"))
            self.elements.append(Spacer(1, 12))

            # Sequencing Statistics
            if 'file_counters' in masterdf_dict:
                counters = masterdf_dict['file_counters']
                seq_stats = [
                    ("BAM Files", f"{self._format_number(counters.get('bam_passed', 0))} passed, "
                               f"{self._format_number(counters.get('bam_failed', 0))} failed"),
                    ("Mapped Reads", f"{self._format_number(counters.get('mapped_count', 0))} total "
                                  f"({self._format_number(counters.get('pass_mapped_count', 0))} passed, "
                                  f"{self._format_number(counters.get('fail_mapped_count', 0))} failed)"),
                    ("Unmapped Reads", f"{self._format_number(counters.get('unmapped_count', 0))} total "
                                    f"({self._format_number(counters.get('pass_unmapped_count', 0))} passed, "
                                    f"{self._format_number(counters.get('fail_unmapped_count', 0))} failed)"),
                    ("Total Bases", f"{self._format_number(counters.get('bases_count', 0))} "
                                 f"({self._format_number(counters.get('pass_bases_count', 0))} passed, "
                                 f"{self._format_number(counters.get('fail_bases_count', 0))} failed)")
                ]
                self.elements.append(self._create_info_table(seq_stats, "Sequencing Statistics"))
                self.elements.append(Spacer(1, 12))

            # Add summary to summary section
            summary_text = (
                f"Sample {sample_id} was sequenced using {self._convert_to_space_separated_string(masterdf_dict['devices'])} "
                f"with flowcell {self._convert_to_space_separated_string(masterdf_dict['flowcell_ids'])}. "
                f"Run started {self._format_timestamp(masterdf_dict['run_time'])}."
            )
            self.summary_elements.append(Paragraph("Run Information", self.styles.styles["Heading3"]))
            self.summary_elements.append(Paragraph(summary_text, self.styles.styles["Normal"]))
            self.summary_elements.append(Spacer(1, 12))

        except Exception as e:
            logger.error(f"Error generating run data section: {str(e)}")
            logger.debug("Exception details:", exc_info=True)
            self.elements.append(
                Paragraph("Error generating run data section", self.styles.styles["Normal"])
            )
