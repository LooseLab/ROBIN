"""
report.py

This module contains the main report class that coordinates the generation of the PDF report.
"""

import os
import logging
import pandas as pd
from reportlab.platypus import SimpleDocTemplate
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import inch
from .styling.styles import ReportStyles
from robin import fonts

logger = logging.getLogger(__name__)


class RobinReport:
    """Main class for generating ROBIN PDF reports."""

    def __init__(self, filename, output):
        """Initialize the report generator.

        Args:
            filename: Output PDF filename
            output: Directory containing analysis output files
        """
        self.filename = filename
        self.output = output
        self.sample_id = os.path.basename(os.path.normpath(output))

        # Handle filename with None prefix
        if filename.startswith("None"):
            final_folder = os.path.basename(os.path.normpath(output))
            self.filename = filename.replace("None", final_folder, 1)
            self.sample_id = final_folder

        # Initialize styling
        self.fonts_dir = os.path.join(os.path.dirname(os.path.abspath(fonts.__file__)))
        self.styles = ReportStyles(self.fonts_dir)

        # Initialize document elements
        self.elements_summary = []
        self.elements = []
        self.end_of_report_elements = []

        # Load master data
        self.masterdf = self._load_master_data()
        self.centreID = self._get_centre_id()

        # Create document
        self.doc = self._create_document()

        # Initialize sections
        self.sections = []
        self._initialize_sections()

    def _load_master_data(self):
        """Load the master data file if it exists."""
        master_path = os.path.join(self.output, "master.csv")
        if os.path.exists(master_path):
            return pd.read_csv(master_path, index_col=0, header=None)
        return None

    def _get_centre_id(self):
        """Get the centre ID from master data."""
        try:
            return (
                self.masterdf.loc["centreID"][1] if self.masterdf is not None else None
            )
        except KeyError:
            return None

    def _create_document(self):
        """Create the PDF document with proper margins and settings."""
        return SimpleDocTemplate(
            self.filename,
            pagesize=A4,
            rightMargin=0.75 * inch,
            leftMargin=0.75 * inch,
            topMargin=1.35 * inch,
            bottomMargin=0.75 * inch,
        )

    def _initialize_sections(self):
        """Initialize all report sections."""
        # Import sections here to avoid circular imports
        from .sections.classification import ClassificationSection
        from .sections.cnv import CNVSection
        from .sections.fusion import FusionSection
        from .sections.coverage import CoverageSection
        from .sections.mgmt import MGMTSection
        from .sections.run_data import RunDataSection
        from .sections.disclaimer import DisclaimerSection

        # Add sections in order
        self.sections = [
            ClassificationSection(self),
            CNVSection(self),
            FusionSection(self),
            CoverageSection(self),
            MGMTSection(self),
            RunDataSection(self),
            DisclaimerSection(self),
        ]

    def generate_report(self):
        """Generate the complete PDF report.

        Returns:
            Path to the generated PDF file
        """
        try:
            logger.info("Starting report generation")

            # Process each section
            for section in self.sections:
                try:
                    section.add_content()
                    summary_elements, main_elements = section.get_elements()
                    self.elements_summary.extend(summary_elements)
                    self.elements.extend(main_elements)
                except Exception as e:
                    logger.error(
                        f"Error processing section {section.__class__.__name__}: {e}",
                        exc_info=True,
                    )

            # Combine all elements
            logger.info("Combining elements for final PDF")
            final_elements = (
                self.elements_summary + self.elements + self.end_of_report_elements
            )

            # Build the PDF
            from .header_footer import header_footer_canvas_factory

            self.doc.multiBuild(
                final_elements,
                canvasmaker=header_footer_canvas_factory(
                    self.sample_id, self.centreID, self.styles, self.fonts_dir
                ),
            )

            logger.info(f"PDF created: {self.filename}")
            return self.filename
        except Exception as e:
            logger.error(f"Error generating report: {e}", exc_info=True)
            raise


def create_pdf(filename, output):
    """Create a PDF report from ROBIN analysis results.

    Args:
        filename: Output PDF filename
        output: Directory containing analysis output files

    Returns:
        Path to the generated PDF file
    """
    report = RobinReport(filename, output)
    return report.generate_report()
