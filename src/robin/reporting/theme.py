"""Authoritative visual tokens for ROBIN PDF reports.

The PDF theme is deliberately independent from the browser theme.  It uses
the same semantic meanings where useful, but the values are chosen for print
contrast and predictable ReportLab rendering.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class ReportTheme:
    """Colours, typography and spacing shared by PDF report components."""

    brand: str = "#15803D"
    text: str = "#0F172A"
    muted: str = "#475569"
    border: str = "#E2E8F0"
    surface: str = "#FFFFFF"
    surface_muted: str = "#F8FAFC"
    positive: str = "#15803D"
    negative: str = "#DC2626"
    warning: str = "#B45309"
    information: str = "#2563EB"
    font_regular: str = "FiraSans"
    font_medium: str = "FiraSans-Medium"
    font_bold: str = "FiraSans-Bold"
    font_mono: str = "FiraMono"
    body_size: float = 10
    body_leading: float = 14
    small_size: float = 8.5
    table_size: float = 9
    heading_size: float = 16
    subheading_size: float = 12


DEFAULT_REPORT_THEME = ReportTheme()
