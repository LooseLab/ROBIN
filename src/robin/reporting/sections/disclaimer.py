"""
disclaimer.py

This module contains the disclaimer section of the report.
"""

from reportlab.platypus import PageBreak, Paragraph, Spacer
from ..sections.base import ReportSection
from .disclaimer_text import EXTENDED_DISCLAIMER_TEXT


class DisclaimerSection(ReportSection):
    """Section containing the research use disclaimer."""

    def add_content(self):
        """Add the disclaimer content to the report."""
        self.elements.append(PageBreak())
        self.add_section_header("Disclaimer")

        # Split the text into paragraphs and create separate Paragraph objects
        paragraphs = EXTENDED_DISCLAIMER_TEXT.split('\n\n')
        for p in paragraphs:
            self.elements.append(Paragraph(p, self.styles.styles["Normal"]))
            self.elements.append(Spacer(1, 12))  # Add some space between paragraphs
