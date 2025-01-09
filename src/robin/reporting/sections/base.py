"""
base.py

This module contains the base class for report sections.
"""

from abc import ABC, abstractmethod
from reportlab.platypus import Paragraph, Spacer, PageBreak, Image, Table
import logging

logger = logging.getLogger(__name__)


class ReportSection(ABC):
    """Base class for all report sections."""

    def __init__(self, report):
        """Initialize the report section.

        Args:
            report: The parent report object containing styles and common functionality
        """
        self.report = report
        self.styles = report.styles
        self.elements = []
        self.summary_elements = []

    @abstractmethod
    def add_content(self):
        """Add the section's content to the report.

        This method should be implemented by each section to add its specific content.
        It should populate both self.elements and self.summary_elements as needed.
        """
        pass

    def add_section_header(self, title, level=1):
        """Add a section header with consistent styling.

        Args:
            title: The header text
            level: The heading level (1 or 2)
        """
        self.elements.append(Spacer(1, 6))
        self.elements.append(Paragraph(title, self.styles.styles[f"Heading{level}"]))
        self.elements.append(Spacer(1, 4))

    def add_section_divider(self):
        """Add a visual divider between sections."""
        from reportlab.platypus import HRFlowable

        self.elements.append(
            HRFlowable(
                width="100%",
                thickness=1,
                color=self.styles.COLORS["border"],
                spaceBefore=12,
                spaceAfter=12,
            )
        )

    def add_summary_card(self, content):
        """Add a summary card with consistent styling.

        Args:
            content: The card content (can include HTML-like formatting)
        """
        self.summary_elements.append(
            Paragraph(content, self.styles.styles["SummaryCard"])
        )
        self.summary_elements.append(Spacer(1, 8))

    def create_table(self, data, style=None, col_widths=None):
        """Create a table with consistent styling.

        Args:
            data: List of lists containing table data
            style: Optional TableStyle object (defaults to MODERN_TABLE_STYLE)
            col_widths: Optional list of column widths

        Returns:
            Table object with applied styling
        """
        if style is None:
            style = self.styles.MODERN_TABLE_STYLE

        table = Table(data, colWidths=col_widths)
        table.setStyle(style)
        return table

    def add_figure(self, img, caption=None, width=None, height=None):
        """Add a figure with optional caption.

        Args:
            img: Image data or path
            caption: Optional caption text
            width: Optional width override
            height: Optional height override
        """
        self.elements.append(Spacer(1, 8))

        if width is not None or height is not None:
            self.elements.append(Image(img, width=width, height=height))
        else:
            self.elements.append(Image(img))

        if caption:
            self.elements.append(Spacer(1, 2))
            self.elements.append(Paragraph(caption, self.styles.styles["Caption"]))
        self.elements.append(Spacer(1, 8))

    def log_error(self, message, exception=None):
        """Log an error with consistent formatting.

        Args:
            message: The error message
            exception: Optional exception object
        """
        if exception:
            logger.error(f"{message}: {str(exception)}", exc_info=True)
        else:
            logger.error(message)

    def get_elements(self):
        """Get all elements for this section.

        Returns:
            Tuple of (summary_elements, main_elements)
        """
        return self.summary_elements, self.elements
