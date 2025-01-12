"""
base.py

This module contains the base class for report sections.
"""

from abc import ABC, abstractmethod
import logging
from reportlab.platypus import Paragraph, Spacer, Table, TableStyle, Image
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.lib.styles import ParagraphStyle

logger = logging.getLogger(__name__)

class ReportSection(ABC):
    """Base class for report sections."""

    def __init__(self, report):
        """Initialize the section with report reference."""
        self.report = report
        self.styles = report.styles
        self.elements = []
        self.summary_elements = []
        
        # Define standard table styles
        self.MODERN_TABLE_STYLE = TableStyle([
            # Header styling
            ('BACKGROUND', (0, 0), (-1, 0), self.styles.COLORS["background"]),
            ('TEXTCOLOR', (0, 0), (-1, 0), self.styles.COLORS["primary"]),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 8),  # Consistent font size
            
            # Cell styling
            ('BACKGROUND', (0, 1), (-1, -1), colors.white),
            ('TEXTCOLOR', (0, 1), (-1, -1), self.styles.COLORS["text"]),
            ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
            
            # Grid styling
            ('GRID', (0, 0), (-1, -1), 0.25, self.styles.COLORS["border"]),
            ('LINEBELOW', (0, 0), (-1, 0), 0.5, self.styles.COLORS["primary"]),  # Thicker header bottom line
            
            # Spacing
            ('TOPPADDING', (0, 0), (-1, -1), 4),      # Reduced padding
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4),    # Reduced padding
            ('LEFTPADDING', (0, 0), (-1, -1), 4),      # Reduced padding
            ('RIGHTPADDING', (0, 0), (-1, -1), 4),     # Reduced padding
            
            # Alignment
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),      # Top alignment for better text wrapping
            
            # Alternating row colors for better readability
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, self.styles.COLORS["background"]]),
        ])
        
        # Header style for tables
        self.TABLE_HEADER_STYLE = ParagraphStyle(
            'TableHeader',
            parent=self.styles.styles["Normal"],
            fontName='Helvetica-Bold',
            fontSize=8,
            textColor=self.styles.COLORS["primary"],
            leading=10,
            spaceBefore=0,
            spaceAfter=0
        )
        
        # Cell style for tables
        self.TABLE_CELL_STYLE = ParagraphStyle(
            'TableCell',
            parent=self.styles.styles["Normal"],
            fontSize=8,
            leading=10,
            spaceBefore=0,
            spaceAfter=0,
            wordWrap='LTR'  # Enable left-to-right word wrapping
        )

    @abstractmethod
    def add_content(self):
        """Add the section content to the report."""
        pass

    def add_section_header(self, title, level=1):
        """Add a section header with consistent styling."""
        self.elements.append(Spacer(1, 6))
        self.elements.append(
            Paragraph(title, self.styles.styles[f"Heading{level}"])
        )
        self.elements.append(Spacer(1, 4))

    def add_summary_card(self, content):
        """Add a summary card with consistent styling."""
        self.summary_elements.append(
            Paragraph(content, self.styles.styles["SummaryCard"])
        )
        self.summary_elements.append(Spacer(1, 8))

    def create_table(self, data, col_widths=None, repeat_rows=1, auto_col_width=True):
        """Create a table with consistent styling and automatic column width calculation.
        
        Args:
            data: List of lists containing table data
            col_widths: Optional list of column widths
            repeat_rows: Number of header rows to repeat on new pages
            auto_col_width: Whether to automatically calculate column widths
            
        Returns:
            Table object with applied styling
        """
        # Convert all data to Paragraphs with proper styling
        formatted_data = []
        for i, row in enumerate(data):
            formatted_row = []
            for cell in row:
                if isinstance(cell, (int, float)):
                    cell = str(cell)
                style = self.TABLE_HEADER_STYLE if i == 0 else self.TABLE_CELL_STYLE
                formatted_row.append(Paragraph(str(cell), style))
            formatted_data.append(formatted_row)

        # Calculate column widths if not provided and auto_col_width is True
        if col_widths is None and auto_col_width:
            page_width = self.report.doc.width
            num_cols = len(data[0])
            
            # Calculate minimum widths based on content
            min_widths = [0] * num_cols
            for row in formatted_data:
                for i, cell in enumerate(row):
                    min_widths[i] = max(min_widths[i], len(str(cell)) * 5)  # Approximate character width
            
            # Ensure minimum width and adjust to page width
            min_widths = [max(0.4 * inch, w) for w in min_widths]
            total_width = sum(min_widths)
            
            if total_width > page_width:
                # Scale down if total width exceeds page width
                scale = page_width / total_width
                col_widths = [w * scale for w in min_widths]
            else:
                # Distribute remaining space proportionally
                extra_space = page_width - total_width
                col_widths = [w + (extra_space / num_cols) for w in min_widths]

        # Create table with calculated or provided column widths
        table = Table(formatted_data, colWidths=col_widths, repeatRows=repeat_rows)
        table.setStyle(self.MODERN_TABLE_STYLE)
        return table

    def add_figure(self, img, caption=None, width=None, height=None):
        """Add a figure with optional caption."""
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
