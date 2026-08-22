"""ReportLab helpers for MNP-Flex hierarchical classifier summaries."""

from __future__ import annotations

from typing import Any, Dict, List, Optional
from xml.sax.saxutils import escape

from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import Paragraph, Spacer, Table, TableStyle

from robin.analysis.mnpflex_hierarchy import (
    format_mnpflex_hierarchy_score,
    iter_mnpflex_hierarchy_nodes,
    mnpflex_hierarchy_has_content,
)


def append_mnpflex_top_path_paragraph(
    elements: List[Any],
    hierarchy: List[Dict[str, Any]],
    *,
    styles: Any,
    format_score,
) -> None:
    """Append the best-scoring hierarchical path as a paragraph."""
    from robin.analysis.mnpflex_hierarchy import best_mnpflex_hierarchy_path

    best_score, best_path = best_mnpflex_hierarchy_path(hierarchy)
    if not best_path:
        return
    elements.append(
        Paragraph(
            f"<b>Top path</b>: {escape(' > '.join(best_path))} "
            f"({escape(format_score(best_score))})",
            styles.styles["Normal"],
        )
    )
    elements.append(Spacer(1, 4))


def append_mnpflex_classifier_prediction(
    elements: List[Any],
    hierarchy: List[Dict[str, Any]],
    *,
    styles: Any,
    page_width: float,
    include_heading: bool = True,
    heading_style: str = "Heading3",
) -> None:
    """Append a nested classifier prediction tree matching the Epignostix GUI."""
    if not mnpflex_hierarchy_has_content(hierarchy):
        return

    if include_heading:
        elements.append(
            Paragraph("Classifier prediction", styles.styles[heading_style])
        )
        elements.append(Spacer(1, 4))

    label_width = max(page_width - 0.9 * inch, 3.5 * inch)
    score_width = 0.85 * inch
    indent_per_depth = 12

    label_style = ParagraphStyle(
        "MNPFlexHierarchyLabel",
        parent=styles.styles["Normal"],
        fontSize=9,
        leading=12,
        spaceBefore=0,
        spaceAfter=0,
    )
    label_bold_style = ParagraphStyle(
        "MNPFlexHierarchyLabelBold",
        parent=label_style,
        fontName="Helvetica-Bold",
    )
    description_style = ParagraphStyle(
        "MNPFlexHierarchyDescription",
        parent=styles.styles["Smaller"],
        fontSize=8,
        leading=11,
        textColor=colors.HexColor("#64748B"),
        spaceBefore=0,
        spaceAfter=0,
    )
    score_style = ParagraphStyle(
        "MNPFlexHierarchyScore",
        parent=styles.styles["Normal"],
        fontSize=9,
        leading=12,
        alignment=2,
        spaceBefore=0,
        spaceAfter=0,
    )

    tree_rows: List[List[Any]] = []
    tree_styles: List[tuple] = []
    row_idx = 0

    for depth, _path, node in iter_mnpflex_hierarchy_nodes(hierarchy):
        group = escape((node.get("group") or "Unknown").strip() or "Unknown")
        score = escape(format_mnpflex_hierarchy_score(node.get("score")))
        left_pad = 6 + depth * indent_per_depth
        text_style = label_bold_style if depth == 0 else label_style

        tree_rows.append(
            [
                Paragraph(group, text_style),
                Paragraph(f"<b>{score}</b>", score_style),
            ]
        )
        tree_styles.extend(
            [
                ("LEFTPADDING", (0, row_idx), (0, row_idx), left_pad),
                ("RIGHTPADDING", (1, row_idx), (1, row_idx), 6),
                ("ALIGN", (1, row_idx), (1, row_idx), "RIGHT"),
                ("VALIGN", (0, row_idx), (-1, row_idx), "TOP"),
                ("TOPPADDING", (0, row_idx), (-1, row_idx), 3),
                ("BOTTOMPADDING", (0, row_idx), (-1, row_idx), 3),
            ]
        )
        row_idx += 1

        description = (node.get("description") or "").strip()
        if description:
            tree_rows.append([Paragraph(escape(description), description_style), ""])
            tree_styles.extend(
                [
                    ("SPAN", (0, row_idx), (1, row_idx)),
                    ("LEFTPADDING", (0, row_idx), (0, row_idx), left_pad + 8),
                    ("TOPPADDING", (0, row_idx), (-1, row_idx), 0),
                    ("BOTTOMPADDING", (0, row_idx), (-1, row_idx), 6),
                ]
            )
            row_idx += 1

    tree_table = Table(
        tree_rows,
        colWidths=[label_width, score_width],
        hAlign="LEFT",
    )
    tree_table.setStyle(
        TableStyle(
            [
                ("BOX", (0, 0), (-1, -1), 0.5, colors.HexColor("#E2E8F0")),
                ("BACKGROUND", (0, 0), (-1, -1), colors.white),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ]
            + tree_styles
        )
    )
    elements.append(tree_table)
    elements.append(Spacer(1, 6))
