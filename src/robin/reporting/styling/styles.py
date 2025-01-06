"""
styles.py

This module contains styling-related code for the ROBIN report generation.
"""

from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.colors import Color, HexColor, white
from reportlab.lib.enums import TA_LEFT, TA_RIGHT, TA_CENTER, TA_JUSTIFY
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import TableStyle
import os
import logging

logger = logging.getLogger(__name__)

class ReportStyles:
    """Encapsulates all styling-related functionality for the report."""
    
    def __init__(self, fonts_dir):
        """Initialize styling with font directory path."""
        self.fonts_dir = fonts_dir
        self.use_default_font = True
        self.register_fonts()
        
        # Initialize styles
        self.styles = getSampleStyleSheet()
        self._setup_colors()
        self._setup_styles()
        self._setup_table_style()

    def register_fonts(self):
        """Register custom fonts for the report."""
        font_files = {
            "FiraSans": "fira-sans-v16-latin-regular.ttf",
            "FiraSans-Medium": "fira-sans-v16-latin-500.ttf",
            "FiraSans-Bold": "fira-sans-v16-latin-700.ttf",
            "FiraMono": "fira-mono-v14-latin-regular.ttf",
        }

        # Default to Helvetica if Fira Sans is not available
        self.use_default_font = True

        for font_name, font_file in font_files.items():
            font_path = os.path.join(self.fonts_dir, font_file)
            if os.path.exists(font_path):
                try:
                    pdfmetrics.registerFont(TTFont(font_name, font_path))
                    if font_name == "FiraSans":
                        self.use_default_font = False
                except Exception as e:
                    logger.warning(f"Could not register font {font_name}: {e}")
            else:
                logger.warning(f"Font file not found: {font_path}")

        if self.use_default_font:
            logger.info("Using Helvetica as fallback font")

    def _setup_colors(self):
        """Set up color palette for consistent branding."""
        self.COLORS = {
            "primary": HexColor("#4F9153"),      # ROBIN theme green
            "secondary": HexColor("#367AB3"),    # Blue accent
            "text": HexColor("#37474f"),
            "success": HexColor("#2e7d32"),
            "warning": HexColor("#f57c00"),
            "error": HexColor("#c62828"),
            "muted": HexColor("#546e7a"),
            "border": HexColor("#e9ecef"),
            "background": HexColor("#f8f9fa"),
        }

    def _setup_styles(self):
        """Set up document styles with modern design."""
        # Get a fresh stylesheet and override any italic styles with regular font
        self.styles = getSampleStyleSheet()
        for style_name in self.styles.byName:
            if "Italic" in style_name:
                self.styles[style_name].fontName = (
                    "FiraSans" if not self.use_default_font else "Helvetica"
                )

        # Base fonts to use
        if self.use_default_font:
            base_font = "Helvetica"
            bold_font = "Helvetica-Bold"
            medium_font = "Helvetica-Bold"  # Helvetica doesn't have medium weight
            mono_font = "Courier"
        else:
            base_font = "FiraSans"
            bold_font = "FiraSans-Bold"
            medium_font = "FiraSans-Medium"
            mono_font = "FiraMono"

        # Update base styles
        for style_name in ["Title", "Heading1", "Heading2", "Normal"]:
            if style_name == "Title":
                self.styles[style_name].fontName = bold_font
                self.styles[style_name].fontSize = 18
                self.styles[style_name].spaceAfter = 14
                self.styles[style_name].spaceBefore = 7
                self.styles[style_name].textColor = self.COLORS["primary"]
                self.styles[style_name].alignment = 1
            elif style_name == "Heading1":
                self.styles[style_name].fontName = bold_font
                self.styles[style_name].fontSize = 14
                self.styles[style_name].spaceAfter = 7
                self.styles[style_name].spaceBefore = 4
                self.styles[style_name].textColor = self.COLORS["primary"]
            elif style_name == "Heading2":
                self.styles[style_name].fontName = medium_font
                self.styles[style_name].fontSize = 12
                self.styles[style_name].spaceAfter = 6
                self.styles[style_name].spaceBefore = 4
                self.styles[style_name].textColor = self.COLORS["secondary"]
            else:  # Normal
                self.styles[style_name].fontName = base_font
                self.styles[style_name].fontSize = 9
                self.styles[style_name].leading = 14
                self.styles[style_name].spaceBefore = 2
                self.styles[style_name].spaceAfter = 2
                self.styles[style_name].textColor = self.COLORS["text"]
                # Add support for HTML-like tags
                self.styles[style_name].allowWidows = 0
                self.styles[style_name].allowOrphans = 0
                self.styles[style_name].wordWrap = 'CJK'

        # Add custom styles
        custom_styles = {
            "Smaller": {
                "fontSize": 7,
                "leading": 11,
                "textColor": self.COLORS["text"],
                "fontName": base_font,
            },
            "Bold": {
                "fontSize": 10,
                "leading": 12,
                "textColor": self.COLORS["primary"],
                "fontName": bold_font,
            },
            "Emphasis": {
                "fontSize": 10,
                "leading": 12,
                "textColor": self.COLORS["secondary"],
                "fontName": medium_font,
            },
            "SummaryCard": {
                "fontSize": 11,
                "leading": 14,
                "textColor": self.COLORS["primary"],
                "backColor": self.COLORS["background"],
                "borderColor": self.COLORS["border"],
                "borderWidth": 1,
                "borderPadding": 8,
                "spaceBefore": 8,
                "spaceAfter": 8,
                "bulletIndent": 0,
                "leftIndent": 8,
                "rightIndent": 8,
                "fontName": medium_font,
                "allowWidows": 0,
                "allowOrphans": 0,
                "wordWrap": "CJK",
                "alignment": 0,  # Left alignment
            },
            "Metric": {
                "fontSize": 12,
                "leading": 16,
                "textColor": self.COLORS["secondary"],
                "alignment": 1,
                "spaceBefore": 4,
                "spaceAfter": 4,
                "fontName": bold_font,
            },
            "Caption": {
                "fontSize": 9,
                "leading": 11,
                "textColor": self.COLORS["muted"],
                "alignment": 1,
                "spaceBefore": 4,
                "spaceAfter": 12,
                "fontName": base_font,
            },
            "MonoText": {
                "fontSize": 9,
                "leading": 12,
                "fontName": mono_font,
                "textColor": self.COLORS["text"],
                "backColor": self.COLORS["background"],
            },
        }

        for style_name, style_props in custom_styles.items():
            style_props["parent"] = self.styles["Normal"]
            # Remove existing style if it exists
            if style_name in self.styles:
                del self.styles[style_name]
            self.styles.add(ParagraphStyle(name=style_name, **style_props))

    def _setup_table_style(self):
        """Set up the default table style."""
        base_font = "FiraSans" if not self.use_default_font else "Helvetica"
        bold_font = "FiraSans-Bold" if not self.use_default_font else "Helvetica-Bold"

        self.MODERN_TABLE_STYLE = TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), self.COLORS["background"]),
            ("TEXTCOLOR", (0, 0), (-1, 0), self.COLORS["primary"]),
            ("FONTNAME", (0, 0), (-1, 0), bold_font),  # Header uses bold
            ("FONTNAME", (0, 1), (-1, -1), base_font),  # Body uses regular
            ("FONTSIZE", (0, 0), (-1, 0), 10),
            ("TOPPADDING", (0, 0), (-1, 0), 12),
            ("BOTTOMPADDING", (0, 0), (-1, 0), 12),
            ("LINEBELOW", (0, 0), (-1, 0), 1, self.COLORS["border"]),
            ("BACKGROUND", (0, 1), (-1, -1), white),
            ("TEXTCOLOR", (0, 1), (-1, -1), self.COLORS["text"]),
            ("FONTSIZE", (0, 1), (-1, -1), 9),
            ("TOPPADDING", (0, 1), (-1, -1), 8),
            ("BOTTOMPADDING", (0, 1), (-1, -1), 8),
            ("LEFTPADDING", (0, 0), (-1, -1), 16),
            ("RIGHTPADDING", (0, 0), (-1, -1), 16),
            ("GRID", (0, 0), (-1, -1), 0.5, self.COLORS["border"]),
            ("ALIGN", (0, 0), (-1, -1), "LEFT"),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [white, self.COLORS["background"]]),
        ]) 