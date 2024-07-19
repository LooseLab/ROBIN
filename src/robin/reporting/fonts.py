"""
fonts.py

This module handles font registration for the PDF report generation.
"""

import os
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase import pdfmetrics
from robin import fonts


def register_fonts(fonts_dir):
    """
    Registers the custom fonts for use in the PDF.

    Args:
        fonts_dir (str): The directory where the font files are located.
    """
    pdfmetrics.registerFont(
        TTFont(
            "FiraSans",
            os.path.join(
                os.path.dirname(os.path.abspath(fonts.__file__)),
                "fira-sans-v16-latin-regular.ttf",
            ),
        )
    )
    pdfmetrics.registerFont(
        TTFont(
            "FiraSans-Bold",
            os.path.join(
                os.path.dirname(os.path.abspath(fonts.__file__)),
                "fira-sans-v16-latin-700.ttf",
            ),
        )
    )
