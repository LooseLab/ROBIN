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
    """Creates a canvas with headers and footers."""
    def header_footer_canvas(canvas, doc):
        canvas.saveState()
        width, height = A4

        # Header
        # Draw header bar with ROBIN green theme
        canvas.setFillColor(colors.HexColor('#4F9153'))  # ROBIN theme green
        canvas.rect(0, height - 1.2*inch, width, 1.2*inch, fill=1, stroke=0)

        # Add ROBIN logo
        logo_path = os.path.join(os.path.dirname(os.path.abspath(images.__file__)), "ROBIN_logo_small.png")
        if os.path.exists(logo_path):
            # Open and convert logo to RGBA if it isn't already
            img = PILImage.open(logo_path)
            if img.mode != 'RGBA':
                img = img.convert('RGBA')
            
            # Calculate logo dimensions while maintaining aspect ratio
            img_width, img_height = img.size
            aspect = img_height / float(img_width)
            # Position logo in top right with reduced size
            canvas.drawImage(
                logo_path,
                width - 1.75*inch,
                height - 1.1*inch,
                width=1.0*inch,
                height=1.0*inch*aspect,
                mask='auto'
            )

        # Add title and Research Use Only warning
        canvas.setFont("FiraSans-Bold", 16)
        canvas.setFillColor(colors.white)
        canvas.drawString(0.75*inch, height - 0.5*inch, "ROBIN Report")
        canvas.setFont("FiraSans-Bold", 12)
        canvas.setFillColor(colors.HexColor('#FF0000'))  # Red color for warning
        canvas.drawString(0.75*inch, height - 0.8*inch, "RESEARCH USE ONLY")

        # Add sample ID and centre ID
        canvas.setFont("FiraSans-Bold", 10)
        canvas.setFillColor(colors.white)
        canvas.drawString(0.75*inch, height - 1.1*inch, f"Sample ID: {sample_id}")
        if centreID:
            canvas.drawString(3*inch, height - 1.1*inch, f"Centre ID: {centreID}")

        # Footer
        # Draw footer bar with ROBIN green theme
        canvas.setFillColor(colors.HexColor('#4F9153'))  # ROBIN theme green
        canvas.rect(0, 0, width, 0.5*inch, fill=1, stroke=0)
        
        # Add page numbers, timestamp, and version on the left
        canvas.setFont("FiraSans", 8)
        canvas.setFillColor(colors.white)
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        page_text = f"Page {doc.page} | Generated: {timestamp} | Version: {VERSION}"
        canvas.drawString(0.5*inch, 0.2*inch, page_text)
        
        # Add Research Use Only warning on the right
        warning_text = "RESEARCH USE ONLY"
        canvas.drawRightString(width - 0.5*inch, 0.2*inch, warning_text)

        canvas.restoreState()

    class HeaderFooterCanvas(canvas.Canvas):
        def __init__(self, *args, **kwargs):
            canvas.Canvas.__init__(self, *args, **kwargs)
            self.pages = []

        def showPage(self):
            self.pages.append(dict(self.__dict__))
            self._startPage()

        def save(self):
            page_count = len(self.pages)
            for page in self.pages:
                self.__dict__.update(page)
                self._drawPageContent()
                canvas.Canvas.showPage(self)
            canvas.Canvas.save(self)

        def _drawPageContent(self):
            header_footer_canvas(self, self._doctemplate)

    return HeaderFooterCanvas
