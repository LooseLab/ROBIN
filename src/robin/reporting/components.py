"""Reusable ReportLab components for ROBIN reports.

Components own the mechanics of spacing, escaping and table formatting so
sections can focus on preparing scientific results.
"""

from collections.abc import Iterable, Sequence
from html import escape

from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.platypus import Image, KeepTogether, Paragraph, Spacer, Table, TableStyle


def _styles(report_styles):
    return report_styles.styles if hasattr(report_styles, "styles") else report_styles


def _style(report_styles, name, fallback="Normal"):
    styles = _styles(report_styles)
    try:
        return styles[name]
    except KeyError:
        return styles[fallback]


def _color(report_styles, name, fallback):
    return getattr(report_styles, "COLORS", {}).get(name, fallback)


def _paragraph(value, style, *, markup=False):
    text = "" if value is None else str(value)
    return Paragraph(text if markup else escape(text), style)


def section_heading(title, report_styles):
    """Return a major section heading with a small, deliberate gap."""
    return [Paragraph(escape(str(title)), _style(report_styles, "SectionHeading", "Heading1")), Spacer(1, 4)]


def finding_block(summary, report_styles):
    """Return a short finding paragraph kept with its immediate explanation."""
    return [KeepTogether([Paragraph(escape(str(summary)), _style(report_styles, "Finding"))]), Spacer(1, 2)]


def metadata_table(data: Sequence[tuple[str, object]], report_styles, *, title=None, width=None):
    """Build a wrapped metadata table with an optional spanning title row."""
    rows = []
    header_style = _style(report_styles, "TableHeader")
    cell_style = _style(report_styles, "TableCell")
    subsection_style = _style(report_styles, "SubsectionHeading", "Heading2")
    if title:
        rows.append([_paragraph(title, subsection_style)])
    rows.extend(
        [[_paragraph(label, header_style), _paragraph(value, cell_style)] for label, value in data]
    )
    table = Table(rows, colWidths=[1.85 * inch, (width or 6.15 * inch) - 1.85 * inch], repeatRows=1 if title else 0)
    commands = [
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 7),
        ("RIGHTPADDING", (0, 0), (-1, -1), 7),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("GRID", (0, 0), (-1, -1), 0.35, _color(report_styles, "border", colors.lightgrey)),
        ("ROWBACKGROUNDS", (0, 1 if title else 0), (-1, -1), [_color(report_styles, "surface", colors.white), _color(report_styles, "surface_variant", colors.whitesmoke)]),
    ]
    if title:
        commands.extend([
            ("SPAN", (0, 0), (-1, 0)),
            ("BACKGROUND", (0, 0), (-1, 0), _color(report_styles, "surface_muted", colors.whitesmoke)),
            ("LINEBELOW", (0, 0), (-1, 0), 0.8, _color(report_styles, "primary", colors.black)),
        ])
    table.setStyle(TableStyle(commands))
    return table


def metric_grid(metrics: Sequence[tuple[str, object]], report_styles, *, width=None, columns=2):
    """Build a compact label/value grid for headline metrics."""
    rows = []
    cells = [[_paragraph(label, _style(report_styles, "Small")), _paragraph(value, _style(report_styles, "Finding"))] for label, value in metrics]
    for offset in range(0, len(cells), columns):
        row = []
        for label, value in cells[offset:offset + columns]:
            row.extend([label, value])
        while len(row) < columns * 2:
            row.extend(["", ""])
        rows.append(row)
    table = Table(rows, colWidths=[(width or 8 * inch) * 0.18, (width or 8 * inch) * 0.32] * columns)
    table.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LINEBELOW", (0, 0), (-1, -1), 0.35, _color(report_styles, "border", colors.lightgrey)),
    ]))
    return table


def results_table(columns: Sequence[str], rows: Iterable[Sequence[object]], report_styles, *, column_proportions=None, width=None, compact=False, numeric_columns=()):
    """Build a splittable results table using explicit column proportions."""
    styles = _styles(report_styles)
    cell_style = _style(report_styles, "TableCellCompact" if compact else "TableCell")
    header_style = _style(report_styles, "TableHeader")
    data = [[_paragraph(column, header_style) for column in columns]]
    data.extend([[_paragraph(cell, cell_style) for cell in row] for row in rows])
    available = width or 8 * inch
    if column_proportions is None:
        column_proportions = tuple(1 / len(columns) for _ in columns)
    if len(column_proportions) != len(columns):
        raise ValueError("column_proportions must have one value per column")
    total = sum(column_proportions)
    if total <= 0:
        raise ValueError("column_proportions must sum to a positive value")
    widths = [available * proportion / total for proportion in column_proportions]
    table = Table(data, colWidths=widths, repeatRows=1, splitByRow=1, hAlign="LEFT")
    commands = [
        ("BACKGROUND", (0, 0), (-1, 0), _color(report_styles, "surface_muted", colors.whitesmoke)),
        ("TEXTCOLOR", (0, 0), (-1, 0), _color(report_styles, "text", colors.black)),
        ("LINEBELOW", (0, 0), (-1, 0), 0.8, _color(report_styles, "primary", colors.black)),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [_color(report_styles, "surface", colors.white), _color(report_styles, "surface_variant", colors.whitesmoke)]),
        ("GRID", (0, 0), (-1, -1), 0.35, _color(report_styles, "border", colors.lightgrey)),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("ALIGN", (0, 0), (-1, -1), "LEFT"),
        ("ALIGN", (0, 0), (-1, -1), "LEFT"),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 4 if compact else 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4 if compact else 5),
    ]
    for index in numeric_columns:
        commands.append(("ALIGN", (index, 1), (index, -1), "RIGHT"))
    table.setStyle(TableStyle(commands))
    return table


def figure_with_caption(figure, caption, report_styles, *, width=None, height=None):
    """Return a figure and a consistently styled caption."""
    styles = _styles(report_styles)
    image = figure if hasattr(figure, "wrap") else Image(figure, width=width, height=height)
    return [image, Paragraph(escape(str(caption)), _style(report_styles, "Caption")), Spacer(1, 4)]


def warning_box(message, report_styles, *, markup=False):
    """Return a print-readable warning box."""
    return [Paragraph(str(message) if markup else escape(str(message)), _style(report_styles, "Warning")), Spacer(1, 2)]


def unavailable_result(label, reason, report_styles):
    """Make an explicit unavailable state instead of an ambiguous empty block."""
    styles = _styles(report_styles)
    text = f"{label}: {reason}"
    return [Paragraph(escape(text), _style(report_styles, "Small")), Spacer(1, 4)]
