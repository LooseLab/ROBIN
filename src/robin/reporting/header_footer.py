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


class HeaderFooterCanvas(canvas.Canvas):
    """
    A custom canvas class for adding headers and footers to the PDF report.
    """

    def __init__(self, sample_id, centreID, styles, fonts_dir, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.pages = []
        self.sample_id = sample_id  # Store the sample_id
        self.centreID = centreID  # Store the centreID
        self.styles = styles  # Store the styles
        self.fonts_dir = fonts_dir  # Store the fonts directory

    def showPage(self):
        self.pages.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        page_count = len(self.pages)
        for page in self.pages:
            self.__dict__.update(page)
            self.draw_canvas(page_count)
            canvas.Canvas.showPage(self)
        canvas.Canvas.save(self)

    def draw_canvas(self, page_count):
        width, height = A4
        
        # Add header background - using ROBIN green theme
        self.setFillColor(colors.HexColor('#4F9153'))  # Match the ROBIN theme green
        self.rect(0, height-1.0*inch, width, 1.0*inch, fill=True, stroke=0)
        
        # Calculate center point of header (header spans 1 inch)
        header_center = height - 0.5*inch
        
        # Add ROBIN logo - centered vertically
        logo_path = os.path.join(os.path.dirname(images.__file__), "ROBIN_logo_small.png")  # Use the same logo as theme
        if os.path.exists(logo_path):
            # Open and convert logo to RGBA if it isn't already
            img = PILImage.open(logo_path)
            if img.mode != 'RGBA':
                img = img.convert('RGBA')
                
            # Create a new image with theme-colored background
            bg = PILImage.new('RGBA', img.size, (79, 145, 83, 255))  # #4F9153 in RGB
            composite = PILImage.alpha_composite(bg, img)
            
            # Save to temporary file
            temp_path = os.path.join(os.path.dirname(logo_path), 'temp_logo.png')
            composite.save(temp_path, format='PNG')
            
            # Draw the processed logo
            logo_offset = 25.0
            self.drawImage(temp_path, width-1.4*inch, header_center - logo_offset, 
                         width=0.9*inch, height=0.7*inch, preserveAspectRatio=True)
            
            # Clean up temporary file
            os.remove(temp_path)
        
        # Title text - white text on green background with red warning
        header1 = Paragraph(
            '<font name="FiraSans-Bold" color="#FFFFFF" size="14">ROBIN Reports</font>'
            '<font name="FiraSans-Bold" color="#FF0000" size="12"> RESEARCH USE ONLY</font>', 
            self.styles["Bold"]
        )
        w, h = header1.wrap(width - 3*inch, inch)
        header1.drawOn(self, 0.5*inch, header_center + 1.5*12)  # Shift up by 1.5 lines (12 points per line)
        
        # Metadata section - white text on green background
        metadata = [
            f"Sample ID: {self.sample_id}",
            f"Centre ID: {self.centreID}",
            datetime.now().strftime("Report Generated: %Y-%m-%d %H:%M:%S")
        ]
        y_position = header_center - 0.1*inch
        for line in metadata:
            self.setFont("FiraSans", 9)
            self.setFillColor(colors.white)  # White text on green background
            self.drawString(0.5*inch, y_position, line)
            y_position -= 12

        # Footer - using ROBIN green theme
        self.setFillColor(colors.HexColor('#4F9153'))
        self.rect(0, 0, width, 0.35*inch, fill=True, stroke=0)
        
        # Footer text - white on green
        page = f"Sample: {self.sample_id} | Centre: {self.centreID} | ROBIN Version: v{VERSION} | Page {self._pageNumber} of {page_count}"
        self.setFont("FiraSans", 8)
        self.setFillColor(colors.white)
        self.drawString(width/2 - len(page)*2.5, 0.15*inch, page)


def header_footer_canvas_factory(sample_id, centreID, styles, fonts_dir):
    """
    A factory function for creating HeaderFooterCanvas objects.

    Args:
        sample_id (str): The sample ID.
        styles: The styles object.
        fonts_dir: The fonts directory.

    Returns:
        function: A function for creating HeaderFooterCanvas objects.
    """

    def create_canvas(*args, **kwargs):
        return HeaderFooterCanvas(sample_id, centreID, styles, fonts_dir, *args, **kwargs)

    return create_canvas
