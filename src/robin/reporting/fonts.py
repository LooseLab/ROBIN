"""Font registration and fallback resolution for PDF reports."""

from dataclasses import dataclass
import logging
import os

from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class ReportFonts:
    """The actual font names available to every report component."""

    regular: str
    medium: str
    bold: str
    mono: str


def register_report_fonts(fonts_dir: str) -> ReportFonts:
    """Register bundled fonts and return usable names, with safe fallbacks."""

    candidates = {
        "regular": ("FiraSans", "fira-sans-v16-latin-regular.ttf", "Helvetica"),
        "medium": ("FiraSans-Medium", "fira-sans-v16-latin-500.ttf", "Helvetica-Bold"),
        "bold": ("FiraSans-Bold", "fira-sans-v16-latin-700.ttf", "Helvetica-Bold"),
        "mono": ("FiraMono", "fira-mono-v14-latin-regular.ttf", "Courier"),
    }
    resolved = {}
    for role, (font_name, filename, fallback) in candidates.items():
        if font_name in pdfmetrics.getRegisteredFontNames():
            resolved[role] = font_name
            continue
        path = os.path.join(fonts_dir, filename)
        try:
            if os.path.exists(path):
                pdfmetrics.registerFont(TTFont(font_name, path))
                resolved[role] = font_name
                continue
        except Exception as exc:  # pragma: no cover - depends on local font files
            logger.warning("Could not register %s: %s", font_name, exc)
        resolved[role] = fallback

    return ReportFonts(**resolved)


# Backwards-compatible entry point used by older callers.
def register_fonts(fonts_dir):
    return register_report_fonts(fonts_dir)
