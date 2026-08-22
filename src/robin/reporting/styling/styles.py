"""
styles.py

This module contains styling-related code for the ROBIN report generation.
"""

import logging
import os

from reportlab.lib.colors import HexColor
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import TableStyle

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
        """Set up Material Design 3 color palette for consistent branding."""
        self.COLORS = {
            # Primary colors (ROBIN green theme)
            "primary": HexColor("#4F9153"),  # ROBIN theme green (unchanged)
            "primary_container": HexColor("#D1E7DD"),  # Light green container
            "on_primary": HexColor("#FFFFFF"),  # White text on primary
            "on_primary_container": HexColor(
                "#0D3B1E"
            ),  # Dark text on primary container
            # Secondary colors (M3 purple)
            "secondary": HexColor("#6750A4"),  # M3 purple
            "secondary_container": HexColor("#EADDFF"),  # Light purple container
            "on_secondary": HexColor("#FFFFFF"),  # White text on secondary
            "on_secondary_container": HexColor(
                "#21005D"
            ),  # Dark text on secondary container
            # Surface colors
            "surface": HexColor("#FEFBFF"),  # M3 surface
            "surface_variant": HexColor("#E7E0EC"),  # M3 surface variant
            "on_surface": HexColor("#1C1B1F"),  # M3 on surface
            "on_surface_variant": HexColor("#49454F"),  # M3 on surface variant
            # Outline colors
            "outline": HexColor("#79747E"),  # M3 outline
            "outline_variant": HexColor("#CAC4D0"),  # M3 outline variant
            # Semantic colors
            "success": HexColor("#4CAF50"),  # M3 success green
            "warning": HexColor("#FF9800"),  # M3 warning orange
            "error": HexColor("#F44336"),  # M3 error red
            "info": HexColor("#2196F3"),  # M3 info blue
            # Legacy colors (for backward compatibility)
            "text": HexColor("#1C1B1F"),  # M3 on surface
            "muted": HexColor("#49454F"),  # M3 on surface variant
            "border": HexColor("#CAC4D0"),  # M3 outline variant
            "background": HexColor("#FEFBFF"),  # M3 surface
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

        # Update base styles with enhanced M3 typography
        for style_name in ["Title", "Heading1", "Heading2", "Normal"]:
            if style_name == "Title":
                self.styles[style_name].fontName = bold_font
                self.styles[style_name].fontSize = 20  # Enhanced from 18
                self.styles[style_name].leading = 28  # M3 leading ratio
                self.styles[style_name].spaceAfter = 16
                self.styles[style_name].spaceBefore = 8
                self.styles[style_name].textColor = self.COLORS["primary"]
                self.styles[style_name].alignment = 1
                self.styles[style_name].allowWidows = 0
                self.styles[style_name].allowOrphans = 0
            elif style_name == "Heading1":
                self.styles[style_name].fontName = bold_font
                self.styles[style_name].fontSize = 16
                self.styles[style_name].leading = 24
                self.styles[style_name].spaceAfter = 6
                self.styles[style_name].spaceBefore = 4
                self.styles[style_name].textColor = self.COLORS["primary"]
                self.styles[style_name].allowWidows = 0
                self.styles[style_name].allowOrphans = 0
            elif style_name == "Heading2":
                self.styles[style_name].fontName = medium_font
                self.styles[style_name].fontSize = 14
                self.styles[style_name].leading = 20
                self.styles[style_name].spaceAfter = 6
                self.styles[style_name].spaceBefore = 4
                self.styles[style_name].textColor = self.COLORS["secondary"]
                self.styles[style_name].allowWidows = 0
                self.styles[style_name].allowOrphans = 0
            else:  # Normal - body text throughout report
                self.styles[style_name].fontName = base_font
                self.styles[style_name].fontSize = 10
                self.styles[style_name].leading = 14
                self.styles[style_name].spaceBefore = 3
                self.styles[style_name].spaceAfter = 3
                self.styles[style_name].textColor = self.COLORS["on_surface"]
                # Enhanced typography settings
                self.styles[style_name].allowWidows = 0
                self.styles[style_name].allowOrphans = 0
                self.styles[style_name].wordWrap = "CJK"
                self.styles[style_name].hyphenationLang = "en_US"

        # Add custom styles with enhanced Material Design 3 typography
        custom_styles = {
            "Smaller": {
                "fontSize": 8,  # Enhanced from 7
                "leading": 12,  # M3 leading ratio
                "textColor": self.COLORS["on_surface_variant"],
                "fontName": base_font,
                "spaceAfter": 2,
                "spaceBefore": 1,
            },
            "Bold": {
                "fontSize": 12,  # Enhanced from 10
                "leading": 16,  # M3 leading ratio
                "textColor": self.COLORS["primary"],
                "fontName": bold_font,
                "spaceAfter": 4,
                "spaceBefore": 2,
            },
            "Emphasis": {
                "fontSize": 12,  # Enhanced from 10
                "leading": 16,  # M3 leading ratio
                "textColor": self.COLORS["secondary"],
                "fontName": medium_font,
                "spaceAfter": 4,
                "spaceBefore": 2,
            },
            "SummaryCard": {
                "fontSize": 10,  # Body text
                "leading": 14,
                "textColor": self.COLORS["on_surface"],
                "backColor": self.COLORS["surface"],
                "borderColor": self.COLORS["outline_variant"],
                "borderWidth": 1,
                "borderPadding": 12,
                "spaceBefore": 12,
                "spaceAfter": 12,
                "bulletIndent": 0,
                "leftIndent": 12,
                "rightIndent": 12,
                "fontName": base_font,
                "allowWidows": 0,
                "allowOrphans": 0,
                "wordWrap": "CJK",
                "alignment": 0,  # Left alignment
            },
            "Metric": {
                "fontSize": 14,  # Enhanced from 12
                "leading": 20,  # M3 leading ratio
                "textColor": self.COLORS["secondary"],
                "alignment": 1,
                "spaceBefore": 6,  # Enhanced from 4
                "spaceAfter": 6,  # Enhanced from 4
                "fontName": bold_font,
            },
            "Caption": {
                "fontSize": 10,
                "leading": 14,
                "textColor": self.COLORS["on_surface_variant"],
                "alignment": 1,
                "spaceBefore": 2,
                "spaceAfter": 4,
                "fontName": base_font,
            },
            "MonoText": {
                "fontSize": 10,  # Enhanced from 9
                "leading": 14,  # M3 leading ratio
                "fontName": mono_font,
                "textColor": self.COLORS["on_surface"],
                "backColor": self.COLORS["surface_variant"],
                "spaceBefore": 3,
                "spaceAfter": 3,
            },
            # Enhanced M3 typography scale
            "DisplayLarge": {
                "fontSize": 24,  # Enhanced from 22
                "leading": 32,  # M3 leading ratio
                "textColor": self.COLORS["on_surface"],
                "fontName": bold_font,
                "spaceAfter": 20,  # Enhanced from 16
                "spaceBefore": 10,  # Enhanced from 8
                "allowWidows": 0,
                "allowOrphans": 0,
            },
            "DisplayMedium": {
                "fontSize": 16,  # Reduced for Summary heading
                "leading": 20,  # Adjusted to match fontSize
                "textColor": self.COLORS["on_surface"],
                "fontName": bold_font,
                "spaceAfter": 12,  # Adjusted spacing
                "spaceBefore": 6,  # Adjusted spacing
                "allowWidows": 0,
                "allowOrphans": 0,
            },
            "HeadlineLarge": {
                "fontSize": 20,
                "leading": 28,
                "textColor": self.COLORS["on_surface"],
                "fontName": bold_font,
                "spaceAfter": 8,
                "spaceBefore": 4,
                "allowWidows": 0,
                "allowOrphans": 0,
            },
            "HeadlineMedium": {
                "fontSize": 18,  # Enhanced from 16
                "leading": 24,  # M3 leading ratio
                "textColor": self.COLORS["on_surface"],
                "fontName": bold_font,
                "spaceAfter": 12,  # Enhanced from 10
                "spaceBefore": 6,  # Enhanced from 5
                "allowWidows": 0,
                "allowOrphans": 0,
            },
            "HeadlineSmall": {
                "fontSize": 16,  # Enhanced from 14
                "leading": 20,  # M3 leading ratio
                "textColor": self.COLORS["on_surface"],
                "fontName": bold_font,
                "spaceAfter": 10,  # Enhanced from 8
                "spaceBefore": 5,  # Enhanced from 4
                "allowWidows": 0,
                "allowOrphans": 0,
            },
            "LabelLarge": {
                "fontSize": 10,
                "leading": 14,
                "textColor": self.COLORS["on_surface"],
                "fontName": medium_font,
                "spaceAfter": 4,
                "spaceBefore": 2,
            },
            "LabelMedium": {
                "fontSize": 9,
                "leading": 12,
                "textColor": self.COLORS["on_surface_variant"],
                "fontName": medium_font,
                "spaceAfter": 3,
                "spaceBefore": 1,
            },
            "LabelSmall": {
                "fontSize": 8,
                "leading": 10,
                "textColor": self.COLORS["on_surface_variant"],
                "fontName": medium_font,
                "spaceAfter": 2,
                "spaceBefore": 1,
            },
            "Alert": {
                "fontSize": 10,
                "leading": 14,
                "textColor": self.COLORS["on_surface"],
                "backColor": self.COLORS["surface_variant"],
                "borderColor": self.COLORS["outline"],
                "borderWidth": 1,
                "borderPadding": 6,
                "spaceBefore": 6,
                "spaceAfter": 6,
                "leftIndent": 6,
                "rightIndent": 6,
                "fontName": base_font,
            },
            "Success": {
                "fontSize": 10,
                "leading": 14,
                "textColor": self.COLORS["success"],
                "backColor": self.COLORS["surface"],
                "borderColor": self.COLORS["success"],
                "borderWidth": 1,
                "borderPadding": 6,
                "spaceBefore": 6,
                "spaceAfter": 6,
                "leftIndent": 6,
                "rightIndent": 6,
                "fontName": base_font,
            },
            "Warning": {
                "fontSize": 10,
                "leading": 14,
                "textColor": self.COLORS["warning"],
                "backColor": self.COLORS["surface"],
                "borderColor": self.COLORS["warning"],
                "borderWidth": 1,
                "borderPadding": 6,
                "spaceBefore": 6,
                "spaceAfter": 6,
                "leftIndent": 6,
                "rightIndent": 6,
                "fontName": base_font,
            },
            "Error": {
                "fontSize": 10,
                "leading": 14,
                "textColor": self.COLORS["error"],
                "backColor": self.COLORS["surface"],
                "borderColor": self.COLORS["error"],
                "borderWidth": 1,
                "borderPadding": 6,
                "spaceBefore": 6,
                "spaceAfter": 6,
                "leftIndent": 6,
                "rightIndent": 6,
                "fontName": base_font,
            },
            # Advanced M3 component styles
            "Card": {
                "fontSize": 12,
                "leading": 18,
                "textColor": self.COLORS["on_surface"],
                "backColor": self.COLORS["surface"],
                "borderColor": self.COLORS["outline_variant"],
                "borderWidth": 1,
                "borderPadding": 16,
                "spaceBefore": 16,
                "spaceAfter": 16,
                "leftIndent": 16,
                "rightIndent": 16,
                "fontName": base_font,
                "allowWidows": 0,
                "allowOrphans": 0,
                "wordWrap": "CJK",
            },
            "ElevatedCard": {
                "fontSize": 12,
                "leading": 18,
                "textColor": self.COLORS["on_surface"],
                "backColor": self.COLORS["surface"],
                "borderColor": self.COLORS["outline"],
                "borderWidth": 2,
                "borderPadding": 20,
                "spaceBefore": 20,
                "spaceAfter": 20,
                "leftIndent": 20,
                "rightIndent": 20,
                "fontName": base_font,
                "allowWidows": 0,
                "allowOrphans": 0,
                "wordWrap": "CJK",
            },
            "InfoCard": {
                # Body text (10pt) for compact header box
                "fontSize": 10,
                "leading": 14,
                "textColor": self.COLORS["on_surface"],
                "backColor": self.COLORS["primary_container"],
                "borderColor": self.COLORS["primary"],
                "borderWidth": 1,
                "borderPadding": 6,
                "spaceBefore": 6,
                "spaceAfter": 6,
                "leftIndent": 6,
                "rightIndent": 6,
                "fontName": base_font,
                "allowWidows": 0,
                "allowOrphans": 0,
                "wordWrap": "CJK",
            },
            "CodeBlock": {
                "fontSize": 10,
                "leading": 14,
                "textColor": self.COLORS["on_surface"],
                "backColor": self.COLORS["surface_variant"],
                "borderColor": self.COLORS["outline_variant"],
                "borderWidth": 1,
                "borderPadding": 12,
                "spaceBefore": 8,
                "spaceAfter": 8,
                "leftIndent": 12,
                "rightIndent": 12,
                "fontName": mono_font,
                "allowWidows": 0,
                "allowOrphans": 0,
                "wordWrap": "CJK",
            },
            "Quote": {
                "fontSize": 11,
                "leading": 16,
                "textColor": self.COLORS["on_surface_variant"],
                "backColor": self.COLORS["surface"],
                "borderColor": self.COLORS["secondary"],
                "borderWidth": 3,
                "borderPadding": 12,
                "spaceBefore": 12,
                "spaceAfter": 12,
                "leftIndent": 16,
                "rightIndent": 12,
                "fontName": base_font,
                "allowWidows": 0,
                "allowOrphans": 0,
                "wordWrap": "CJK",
                "alignment": 0,
            },
            "Callout": {
                "fontSize": 12,
                "leading": 18,
                "textColor": self.COLORS["on_surface"],
                "backColor": self.COLORS["secondary_container"],
                "borderColor": self.COLORS["secondary"],
                "borderWidth": 1,
                "borderPadding": 16,
                "spaceBefore": 16,
                "spaceAfter": 16,
                "leftIndent": 16,
                "rightIndent": 16,
                "fontName": medium_font,
                "allowWidows": 0,
                "allowOrphans": 0,
                "wordWrap": "CJK",
            },
        }

        for style_name, style_props in custom_styles.items():
            style_props["parent"] = self.styles["Normal"]
            # Remove existing style if it exists
            if style_name in self.styles:
                del self.styles[style_name]
            self.styles.add(ParagraphStyle(name=style_name, **style_props))

    def _setup_table_style(self):
        """Set up Material Design 3 table styles."""
        base_font = "FiraSans" if not self.use_default_font else "Helvetica"
        bold_font = "FiraSans-Bold" if not self.use_default_font else "Helvetica-Bold"

        # Enhanced main table style with advanced M3 design
        self.MODERN_TABLE_STYLE = TableStyle(
            [
                # Header styling with enhanced visual hierarchy
                ("BACKGROUND", (0, 0), (-1, 0), self.COLORS["primary"]),
                ("TEXTCOLOR", (0, 0), (-1, 0), self.COLORS["on_primary"]),
                ("FONTNAME", (0, 0), (-1, 0), bold_font),
                ("FONTSIZE", (0, 0), (-1, 0), 11),
                ("TOPPADDING", (0, 0), (-1, 0), 8),
                ("BOTTOMPADDING", (0, 0), (-1, 0), 8),
                ("LINEBELOW", (0, 0), (-1, 0), 3, self.COLORS["primary"]),
                # Body styling with alternating row colors
                ("BACKGROUND", (0, 1), (-1, -1), self.COLORS["surface"]),
                ("TEXTCOLOR", (0, 1), (-1, -1), self.COLORS["on_surface"]),
                ("FONTNAME", (0, 1), (-1, -1), base_font),
                ("FONTSIZE", (0, 1), (-1, -1), 10),
                ("TOPPADDING", (0, 1), (-1, -1), 6),
                ("BOTTOMPADDING", (0, 1), (-1, -1), 6),
                # Padding and alignment
                ("LEFTPADDING", (0, 0), (-1, -1), 12),
                ("RIGHTPADDING", (0, 0), (-1, -1), 12),
                ("ALIGN", (0, 0), (-1, -1), "LEFT"),
                # Subtle alternating row backgrounds for better readability
                (
                    "ROWBACKGROUNDS",
                    (0, 1),
                    (-1, -1),
                    [self.COLORS["surface"], self.COLORS["surface_variant"]],
                ),
                # Enhanced grid with better contrast
                ("GRID", (0, 0), (-1, -1), 0.75, self.COLORS["outline"]),
            ]
        )

        # Compact table style for dense data
        self.COMPACT_TABLE_STYLE = TableStyle(
            [
                # Header styling with secondary color theme
                ("BACKGROUND", (0, 0), (-1, 0), self.COLORS["secondary"]),
                ("TEXTCOLOR", (0, 0), (-1, 0), self.COLORS["on_secondary"]),
                ("FONTNAME", (0, 0), (-1, 0), bold_font),
                ("FONTSIZE", (0, 0), (-1, 0), 10),
                ("TOPPADDING", (0, 0), (-1, 0), 6),
                ("BOTTOMPADDING", (0, 0), (-1, 0), 6),
                ("LINEBELOW", (0, 0), (-1, 0), 2, self.COLORS["secondary"]),
                # Body styling with subtle alternating rows
                ("BACKGROUND", (0, 1), (-1, -1), self.COLORS["surface"]),
                ("TEXTCOLOR", (0, 1), (-1, -1), self.COLORS["on_surface"]),
                ("FONTNAME", (0, 1), (-1, -1), base_font),
                ("FONTSIZE", (0, 1), (-1, -1), 9),
                ("TOPPADDING", (0, 1), (-1, -1), 5),
                ("BOTTOMPADDING", (0, 1), (-1, -1), 5),
                # Padding and alignment
                ("LEFTPADDING", (0, 0), (-1, -1), 10),
                ("RIGHTPADDING", (0, 0), (-1, -1), 10),
                ("ALIGN", (0, 0), (-1, -1), "LEFT"),
                # Subtle alternating row backgrounds
                (
                    "ROWBACKGROUNDS",
                    (0, 1),
                    (-1, -1),
                    [self.COLORS["surface"], self.COLORS["surface_variant"]],
                ),
                # Enhanced grid
                ("GRID", (0, 0), (-1, -1), 0.5, self.COLORS["outline_variant"]),
            ]
        )

        # Summary table style for key metrics
        self.SUMMARY_TABLE_STYLE = TableStyle(
            [
                # Header styling
                ("BACKGROUND", (0, 0), (-1, 0), self.COLORS["primary"]),
                ("TEXTCOLOR", (0, 0), (-1, 0), self.COLORS["on_primary"]),
                ("FONTNAME", (0, 0), (-1, 0), bold_font),
                ("FONTSIZE", (0, 0), (-1, 0), 11),
                ("TOPPADDING", (0, 0), (-1, 0), 8),
                ("BOTTOMPADDING", (0, 0), (-1, 0), 8),
                # Body styling
                ("BACKGROUND", (0, 1), (-1, -1), self.COLORS["surface"]),
                ("TEXTCOLOR", (0, 1), (-1, -1), self.COLORS["on_surface"]),
                ("FONTNAME", (0, 1), (-1, -1), base_font),
                ("FONTSIZE", (0, 1), (-1, -1), 10),
                ("TOPPADDING", (0, 1), (-1, -1), 6),
                ("BOTTOMPADDING", (0, 1), (-1, -1), 6),
                # Padding and alignment
                ("LEFTPADDING", (0, 0), (-1, -1), 12),
                ("RIGHTPADDING", (0, 0), (-1, -1), 12),
                ("ALIGN", (0, 0), (-1, -1), "LEFT"),
                # Border
                ("BOX", (0, 0), (-1, -1), 1, self.COLORS["outline"]),
            ]
        )
