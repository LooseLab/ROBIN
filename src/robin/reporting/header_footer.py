"""
header_footer.py

This module contains the class for adding headers and footers to the PDF report.
"""

from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import inch
from datetime import datetime
import os
from robin import images
from reportlab.lib import colors
from PIL import Image as PILImage
import io

from robin.__about__ import __version__

VERSION = __version__


def header_footer_canvas_factory(sample_id, centreID, styles, fonts_dir):
    """Factory function to create a header/footer canvas class."""

    class HeaderFooterCanvas(canvas.Canvas):
        def __init__(self, *args, **kwargs):
            canvas.Canvas.__init__(self, *args, **kwargs)
            self.pages = []
            self.width, self.height = A4

        def showPage(self):
            self.pages.append(dict(self.__dict__))
            self._startPage()

        def save(self):
            page_count = len(self.pages)
            for page in self.pages:
                self.__dict__.update(page)
                self._draw_header_footer()
                canvas.Canvas.showPage(self)
            canvas.Canvas.save(self)

        def _draw_header_footer(self):
            self.saveState()

            # Header
            # Draw header divider line with ROBIN green theme
            self.setStrokeColor(colors.HexColor("#4F9153"))  # ROBIN theme green
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

            # Add title and Research Use Only warning with black text
            self.setFont("FiraSans", 16)
            self.setFillColor(colors.black)
            self.drawString(0.5 * inch, self.height - 0.4 * inch, "ROBIN Report")
            self.setFont("FiraSans", 12)
            self.setFillColor(colors.HexColor("#FF0000"))  # Keep warning in red
            self.drawString(0.5 * inch, self.height - 0.65 * inch, "RESEARCH USE ONLY")

            # Add sample ID and centre ID with black text
            self.setFont("FiraSans", 10)
            self.setFillColor(colors.black)
            self.drawString(
                0.5 * inch, self.height - 0.9 * inch, f"Sample ID: {sample_id}"
            )
            if centreID:
                self.drawString(
                    2.75 * inch, self.height - 0.9 * inch, f"Centre ID: {centreID}"
                )

            # Footer divider line
            self.setStrokeColor(colors.HexColor("#4F9153"))  # ROBIN theme green
            self.setLineWidth(2)
            self.line(
                margin, 0.4 * inch, self.width - margin, 0.4 * inch
            )  # Footer divider line

            # Footer text with black text
            self.setFont("FiraSans", 8)
            self.setFillColor(colors.black)
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            page_text = (
                f"Page {self._pageNumber} | Generated: {timestamp} | Version: {VERSION}"
            )
            self.drawString(0.5 * inch, 0.15 * inch, page_text)

            warning_text = "RESEARCH USE ONLY"
            self.setFillColor(colors.HexColor("#FF0000"))  # Keep warning in red
            self.drawRightString(self.width - 0.5 * inch, 0.15 * inch, warning_text)

            self.restoreState()

    return HeaderFooterCanvas
