"""
header_footer.py

This module contains the class for adding headers and footers to the PDF report.
"""

import os
from datetime import datetime

from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas

from robin.__init__ import __version__
from robin.gui import images

VERSION = __version__


def header_footer_canvas_factory(
    sample_id,
    centreID,
    styles,
    fonts_dir,
    generated_by=None,
    generated_at=None,
    robin_commit=None,
):
    """Factory function to create a header/footer canvas class."""

    class HeaderFooterCanvas(canvas.Canvas):
        def __init__(self, *args, **kwargs):
            canvas.Canvas.__init__(self, *args, **kwargs)
            self.pages = []
            self._sync_page_size()

        def _sync_page_size(self) -> None:
            """Use the active page size (portrait or landscape)."""
            pagesize = getattr(self, "_pagesize", None) or A4
            self.width, self.height = pagesize

        def showPage(self):
            self.pages.append(dict(self.__dict__))
            self._startPage()

        def save(self):
            # page_count = len(self.pages)
            for page in self.pages:
                self.__dict__.update(page)
                self._sync_page_size()
                self._draw_header_footer()
                canvas.Canvas.showPage(self)
            canvas.Canvas.save(self)

        def _draw_header_footer(self):
            self.saveState()

            # Header
            # Draw header divider line with ROBIN green theme (M3 primary color)
            self.setStrokeColor(
                colors.HexColor("#4F9153")
            )  # ROBIN theme green (M3 primary)
            self.setLineWidth(4)
            # Calculate margins for 90% width
            margin = self.width * 0.05  # 5% margin on each side
            self.line(
                margin,
                self.height - 1.0 * inch,
                self.width - margin,
                self.height - 1.0 * inch,
            )  # Header divider line

            # Add ROBIN logo
            logo_path = os.path.join(
                os.path.dirname(os.path.abspath(images.__file__)),
                "ROBIN_logo_small.png",
            )
            if os.path.exists(logo_path):
                img = PILImage.open(logo_path)
                if img.mode != "RGBA":
                    img = img.convert("RGBA")
                img_width, img_height = img.size
                aspect = img_height / float(img_width)
                self.drawImage(
                    logo_path,
                    self.width - 1.5 * inch,
                    self.height - 0.95 * inch,
                    width=0.85 * inch,
                    height=0.85 * inch * aspect,
                    mask="auto",
                )

            # Add title and Research Use Only warning with M3 colors
            self.setFont("FiraSans", 16)
            self.setFillColor(colors.HexColor("#1C1B1F"))  # M3 on_surface color
            self.drawString(0.5 * inch, self.height - 0.4 * inch, "ROBIN Report")
            self.setFont("FiraSans", 12)
            self.setFillColor(colors.HexColor("#F44336"))  # M3 error color for warning
            self.drawString(0.5 * inch, self.height - 0.65 * inch, "RESEARCH USE ONLY")

            # Add sample ID and centre ID with M3 colors
            self.setFont("FiraSans", 10)
            self.setFillColor(colors.HexColor("#49454F"))  # M3 on_surface_variant color

            # Calculate text widths
            sample_id_text = f"Sample ID: {sample_id}"
            centre_id_text = f"Centre ID: {centreID}" if centreID else ""

            # Get text widths in points
            sample_id_width = self.stringWidth(sample_id_text, "FiraSans", 10)

            # Draw Sample ID
            self.drawString(0.5 * inch, self.height - 0.9 * inch, sample_id_text)

            # Draw Centre ID only if it exists and ensure it doesn't overlap
            if centreID:
                # Calculate minimum spacing between IDs (0.25 inch)
                min_spacing = 0.25 * inch
                # Calculate x position for Centre ID
                centre_id_x = max(
                    0.5 * inch
                    + sample_id_width
                    + min_spacing,  # Right after Sample ID with minimum spacing
                    2.75 * inch,  # Original position if there's enough space
                )
                # Only draw if it fits on the page (leaving 1.75 inch from right edge for logo)
                if (
                    centre_id_x + self.stringWidth(centre_id_text, "FiraSans", 10)
                    < self.width - 1.75 * inch
                ):
                    self.drawString(
                        centre_id_x, self.height - 0.9 * inch, centre_id_text
                    )

            # Footer divider line
            self.setStrokeColor(
                colors.HexColor("#4F9153")
            )  # ROBIN theme green (M3 primary)
            self.setLineWidth(2)
            self.line(
                margin, 0.4 * inch, self.width - margin, 0.4 * inch
            )  # Footer divider line

            # Footer text with M3 colors
            self.setFont("FiraSans", 8)
            self.setFillColor(colors.HexColor("#49454F"))  # M3 on_surface_variant color
            timestamp = generated_at or datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            page_text = f"Page {self._pageNumber} | Generated: {timestamp}"
            if generated_by:
                page_text += f" by {generated_by}"
            page_text += f" | Version: {VERSION}"
            if robin_commit:
                page_text += f" | Commit: {robin_commit[:12]}"
            self.drawString(0.5 * inch, 0.15 * inch, page_text)

            warning_text = "RESEARCH USE ONLY"
            self.setFillColor(colors.HexColor("#F44336"))  # M3 error color for warning
            self.drawRightString(self.width - 0.5 * inch, 0.15 * inch, warning_text)

            self.restoreState()

    return HeaderFooterCanvas
