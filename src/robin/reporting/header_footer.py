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

class HeaderFooterCanvas(canvas.Canvas):
    """
    A custom canvas class for adding headers and footers to the PDF report.
    """

    def __init__(self, sample_id, styles, fonts_dir, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.pages = []
        self.sample_id = sample_id  # Store the sample_id
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

        # Add first line of the header in bold
        header1 = Paragraph("R.O.B.I.N Reports... RESEARCH USE ONLY", self.styles["Bold"])
        w, h = header1.wrap(width - 2 * inch, inch)
        header1.drawOn(self, inch, height - h - inch + 36)

        # Add second line of the header
        header2 = Paragraph(f"Sample ID: {self.sample_id}", self.styles["Normal"])
        w, h = header2.wrap(width - 2 * inch, inch)
        header2.drawOn(self, inch, height - h - inch + 24)

        # Add third line of the header with date and time
        report_date = datetime.now().strftime(
            "R.O.B.I.N. Report Generated: %Y-%m-%d %H:%M:%S"
        )
        header3 = Paragraph(report_date, self.styles["Smaller"])
        w, h = header3.wrap(width - 2 * inch, inch)
        header3.drawOn(self, inch, height - h - inch + 12)

        # Add logo to the top right corner of the header
        logo_path = os.path.join(
            os.path.dirname(os.path.abspath(images.__file__)), "ROBIN_logo_small.png"
        )
        # logo_path = "src/robin/images/Robin_logo_small.png"  # Replace with the path to your logo
        max_logo_size = 50  # Maximum width and height in pixels
        self.drawImage(
            logo_path,
            width - max_logo_size - inch,
            height - max_logo_size - inch + 36,
            width=max_logo_size,
            height=max_logo_size,
            preserveAspectRatio=True,
            mask="auto",
        )

        # Add footer
        page = f"SampleID: {self.sample_id} - Page {self._pageNumber} of {page_count}"
        x = 190
        self.saveState()
        self.setStrokeColorRGB(0, 0, 0)
        self.setLineWidth(0.5)
        self.line(66, 78, A4[0] - 66, 78)
        self.setFont("FiraSans", 7)
        self.drawString(A4[0] - x, 65, page)
        self.restoreState()


def header_footer_canvas_factory(sample_id, styles, fonts_dir):
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
        return HeaderFooterCanvas(sample_id, styles, fonts_dir, *args, **kwargs)

    return create_canvas
