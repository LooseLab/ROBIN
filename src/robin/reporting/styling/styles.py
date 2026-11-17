"""
styles.py

This module contains styling-related code for the ROBIN report generation.
"""

from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.platypus import TableStyle
import logging
from reportlab.lib.colors import HexColor

from robin.reporting.fonts import register_report_fonts
from robin.reporting.theme import DEFAULT_REPORT_THEME

logger = logging.getLogger(__name__)


class ReportStyles:
    """Encapsulates all styling-related functionality for the report."""

    def __init__(self, fonts_dir, theme=None, fonts=None):
        """Initialize styling with font directory path."""
        self.fonts_dir = fonts_dir
        self.theme = theme or DEFAULT_REPORT_THEME
        self.fonts = fonts or register_report_fonts(fonts_dir)
        self.use_default_font = self.fonts.regular == "Helvetica"

        # Initialize styles
        self.styles = getSampleStyleSheet()
        self._setup_colors()
        self._setup_styles()
        self._setup_table_style()

    def register_fonts(self):
        """Register fonts and refresh the resolved font names."""
        self.fonts = register_report_fonts(self.fonts_dir)
        self.use_default_font = self.fonts.regular == "Helvetica"

    def _setup_colors(self):
        """Set up Material Design 3 color palette for consistent branding."""
        theme = self.theme
        self.COLORS = {
            # Primary colors (ROBIN green theme)
            "primary": HexColor(theme.brand),
            "primary_container": HexColor(theme.surface_muted),
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
            "surface": HexColor(theme.surface),
            "surface_variant": HexColor(theme.surface_muted),
            "surface_muted": HexColor(theme.surface_muted),
            "on_surface": HexColor(theme.text),
            "on_surface_variant": HexColor(theme.muted),
            # Outline colors
            "outline": HexColor(theme.muted),
            "outline_variant": HexColor(theme.border),
            # Semantic colors
            "success": HexColor(theme.positive),
            "warning": HexColor(theme.warning),
            "error": HexColor(theme.negative),
            "info": HexColor(theme.information),
            # Legacy colors (for backward compatibility)
            "text": HexColor(theme.text),
            "muted": HexColor(theme.muted),
            "border": HexColor(theme.border),
            "background": HexColor(theme.surface),
        }

    def _setup_styles(self):
        """Set up document styles with modern design."""
        # Get a fresh stylesheet and override any italic styles with regular font
        self.styles = getSampleStyleSheet()
        for style_name in self.styles.byName:
            if "Italic" in style_name:
                self.styles[style_name].fontName = (
                    self.fonts.regular
                )

        # Base fonts to use
        base_font = self.fonts.regular
        bold_font = self.fonts.bold
        medium_font = self.fonts.medium
        mono_font = self.fonts.mono

        # Update base styles with enhanced M3 typography
        for style_name in ["Title", "Heading1", "Heading2", "Heading3", "Normal"]:
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
                self.styles[style_name].leading = 20
                self.styles[style_name].spaceAfter = 4
                self.styles[style_name].spaceBefore = 2
                self.styles[style_name].textColor = self.COLORS["primary"]
                self.styles[style_name].allowWidows = 0
                self.styles[style_name].allowOrphans = 0
            elif style_name == "Heading2":
                self.styles[style_name].fontName = medium_font
                self.styles[style_name].fontSize = 13
                self.styles[style_name].leading = 17
                self.styles[style_name].spaceAfter = 4
                self.styles[style_name].spaceBefore = 3
                self.styles[style_name].textColor = self.COLORS["on_surface"]
                self.styles[style_name].allowWidows = 0
                self.styles[style_name].allowOrphans = 0
            elif style_name == "Heading3":
                self.styles[style_name].fontName = medium_font
                self.styles[style_name].fontSize = 11
                self.styles[style_name].leading = 14
                self.styles[style_name].spaceAfter = 3
                self.styles[style_name].spaceBefore = 2
                self.styles[style_name].textColor = self.COLORS["on_surface"]
                self.styles[style_name].keepWithNext = True
            else:  # Normal - body text throughout report
                self.styles[style_name].fontName = base_font
                self.styles[style_name].fontSize = 10
                self.styles[style_name].leading = 14
                self.styles[style_name].spaceBefore = 1
                self.styles[style_name].spaceAfter = 2
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
                "fontName": medium_font,
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

        self._setup_report_role_styles(
            base_font=base_font, medium_font=medium_font, bold_font=bold_font, mono_font=mono_font
        )

    def _setup_report_role_styles(self, *, base_font, medium_font, bold_font, mono_font):
        """Add the small, authoritative vocabulary used by new components."""
        t = self.theme
        c = self.COLORS
        roles = {
            "ReportTitle": dict(fontName=bold_font, fontSize=20, leading=24, textColor=c["text"], spaceAfter=8, keepWithNext=True),
            "SectionHeading": dict(fontName=bold_font, fontSize=t.heading_size, leading=18, textColor=c["primary"], spaceBefore=8, spaceAfter=4, keepWithNext=True),
            "SubsectionHeading": dict(fontName=medium_font, fontSize=t.subheading_size, leading=15, textColor=c["text"], spaceBefore=5, spaceAfter=3, keepWithNext=True),
            "Body": dict(fontName=base_font, fontSize=t.body_size, leading=t.body_leading, textColor=c["text"], spaceBefore=1, spaceAfter=3),
            "Small": dict(fontName=base_font, fontSize=t.small_size, leading=10, textColor=c["muted"], spaceBefore=1, spaceAfter=2),
            "Caption": dict(fontName=base_font, fontSize=t.small_size, leading=10, textColor=c["muted"], spaceBefore=2, spaceAfter=4, alignment=1),
            "TableHeader": dict(fontName=medium_font, fontSize=t.table_size, leading=10, textColor=c["text"], spaceBefore=0, spaceAfter=0),
            "TableCell": dict(fontName=base_font, fontSize=t.table_size, leading=10, textColor=c["text"], spaceBefore=0, spaceAfter=0, wordWrap="LTR"),
            "TableCellCompact": dict(fontName=base_font, fontSize=8, leading=9, textColor=c["text"], spaceBefore=0, spaceAfter=0, wordWrap="LTR"),
            "Technical": dict(fontName=mono_font, fontSize=t.small_size, leading=10, textColor=c["text"], spaceBefore=1, spaceAfter=2, wordWrap="LTR"),
            "Finding": dict(fontName=medium_font, fontSize=11, leading=14, textColor=c["text"], spaceBefore=2, spaceAfter=3, keepWithNext=True),
            "Warning": dict(fontName=base_font, fontSize=t.body_size, leading=t.body_leading, textColor=c["warning"], backColor=c["surface_muted"], borderColor=c["warning"], borderWidth=0.75, borderPadding=5, spaceBefore=3, spaceAfter=5),
        }
        for name, props in roles.items():
            if name in self.styles:
                self.styles.byName.pop(name, None)
            self.styles.add(ParagraphStyle(name=name, parent=self.styles["Normal"], **props))

    def _setup_table_style(self):
        """Set up Material Design 3 table styles."""
        base_font = self.fonts.regular
        bold_font = self.fonts.bold

        # Enhanced main table style with advanced M3 design
        self.MODERN_TABLE_STYLE = TableStyle(
            [
                # Header styling with enhanced visual hierarchy
                ("BACKGROUND", (0, 0), (-1, 0), self.COLORS["primary"]),
                ("TEXTCOLOR", (0, 0), (-1, 0), self.COLORS["on_primary"]),
                ("FONTNAME", (0, 0), (-1, 0), bold_font),
                ("FONTSIZE", (0, 0), (-1, 0), 11),
                ("TOPPADDING", (0, 0), (-1, 0), 4),
                ("BOTTOMPADDING", (0, 0), (-1, 0), 4),
                ("LINEBELOW", (0, 0), (-1, 0), 3, self.COLORS["primary"]),
                # Body styling with alternating row colors
                ("BACKGROUND", (0, 1), (-1, -1), self.COLORS["surface"]),
                ("TEXTCOLOR", (0, 1), (-1, -1), self.COLORS["on_surface"]),
                ("FONTNAME", (0, 1), (-1, -1), base_font),
                ("FONTSIZE", (0, 1), (-1, -1), 10),
                ("TOPPADDING", (0, 1), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 1), (-1, -1), 3),
                # Padding and alignment
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
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
                ("TOPPADDING", (0, 0), (-1, 0), 3),
                ("BOTTOMPADDING", (0, 0), (-1, 0), 3),
                ("LINEBELOW", (0, 0), (-1, 0), 2, self.COLORS["secondary"]),
                # Body styling with subtle alternating rows
                ("BACKGROUND", (0, 1), (-1, -1), self.COLORS["surface"]),
                ("TEXTCOLOR", (0, 1), (-1, -1), self.COLORS["on_surface"]),
                ("FONTNAME", (0, 1), (-1, -1), base_font),
                ("FONTSIZE", (0, 1), (-1, -1), 9),
                ("TOPPADDING", (0, 1), (-1, -1), 2),
                ("BOTTOMPADDING", (0, 1), (-1, -1), 2),
                # Padding and alignment
                ("LEFTPADDING", (0, 0), (-1, -1), 5),
                ("RIGHTPADDING", (0, 0), (-1, -1), 5),
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
                ("TOPPADDING", (0, 0), (-1, 0), 4),
                ("BOTTOMPADDING", (0, 0), (-1, 0), 4),
                # Body styling
                ("BACKGROUND", (0, 1), (-1, -1), self.COLORS["surface"]),
                ("TEXTCOLOR", (0, 1), (-1, -1), self.COLORS["on_surface"]),
                ("FONTNAME", (0, 1), (-1, -1), base_font),
                ("FONTSIZE", (0, 1), (-1, -1), 10),
                ("TOPPADDING", (0, 1), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 1), (-1, -1), 3),
                # Padding and alignment
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("ALIGN", (0, 0), (-1, -1), "LEFT"),
                # Border
                ("BOX", (0, 0), (-1, -1), 1, self.COLORS["outline"]),
            ]
        )
