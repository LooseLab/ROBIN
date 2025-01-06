"""
disclaimer.py

This module contains the disclaimer section of the report.
"""

from reportlab.platypus import PageBreak, Paragraph
from ..sections.base import ReportSection

class DisclaimerSection(ReportSection):
    """Section containing the research use disclaimer."""

    def add_content(self):
        """Add the disclaimer content to the report."""
        self.elements.append(PageBreak())
        self.add_section_header("Disclaimer")

        disclaimer_text = (
            "This report and the data contained within it are intended for research use only and "
            "should not be used for direct diagnostic purposes. The methylation-based classification "
            "and other analyses provided here may be considered by neuropathologists as supplementary "
            "information in the context of comprehensive diagnostic assessment, which should include "
            "clinical history, radiological findings, and complete histopathological and molecular evaluation. "
            "The final interpretation and diagnosis should always be made by qualified healthcare professionals "
            "based on all available information."
        )

        self.elements.append(
            Paragraph(disclaimer_text, self.styles.styles["Normal"])
        ) 