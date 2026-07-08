from __future__ import annotations

from functools import lru_cache
from html import escape
from importlib import resources as importlib_resources
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from robin.fusionvis.analysis import analyse_segments
from robin.fusionvis.annotation import AnnotationIndex, parse_annotation
from robin.fusionvis.bam import AlignmentSegment


@lru_cache(maxsize=1)
def _annotation() -> AnnotationIndex:
    path = (
        importlib_resources.files("robin.resources")
        / "gencode.v45.basic.annotation.gff3.gz"
    )
    return parse_annotation(Path(str(path)))


def candidate_from_dataframe(
    subset: pd.DataFrame,
    gene_pair: Iterable[str],
    *,
    min_support: int = 1,
) -> dict[str, Any] | None:
    """Build a FusionVis candidate from ROBIN's selected fusion rows."""
    pair = {str(gene).strip() for gene in gene_pair}
    segments = []
    for row in subset.to_dict("records"):
        gene = str(row.get("col4") or "").strip()
        if gene not in pair:
            continue
        start = int(row.get("reference_start") or 0)
        end = int(row.get("reference_end") or start)
        query_start = int(row.get("read_start") or 0)
        query_end = int(row.get("read_end") or query_start)
        segments.append(
            AlignmentSegment(
                read_name=str(row.get("read_id") or ""),
                chrom=str(row.get("reference_id") or ""),
                start=min(start, end) + 1,
                end=max(start, end),
                strand=str(row.get("strand") or "+"),
                mapq=int(row.get("mapping_quality") or 0),
                query_start=min(query_start, query_end),
                query_end=max(query_start, query_end),
                cigar=str(row.get("cigar") or "*"),
                source=(
                    "supplementary" if bool(row.get("is_supplementary")) else "record"
                ),
                is_supplementary=bool(row.get("is_supplementary")),
                is_secondary=bool(row.get("is_secondary")),
            )
        )

    result = analyse_segments(
        segments,
        _annotation(),
        min_mapq=0,
        min_support=max(1, int(min_support)),
    )
    for candidate in result.get("candidates", []):
        names = {
            str(candidate.get("geneA", {}).get("name") or "").strip(),
            str(candidate.get("geneB", {}).get("name") or "").strip(),
        }
        if names == pair:
            return candidate
    return None


def render_candidate_svg(candidate: dict[str, Any], width: int = 1500) -> str:
    """Render a selected FusionVis candidate as a responsive inline SVG."""
    gene_a = candidate["geneA"]
    gene_b = candidate["geneB"]
    bundle = candidate.get("modelBundle") or {}
    models = list(bundle.get("models") or [])[:4]
    product_height = max(len(models), 1) * 108
    height = 405 + product_height
    left, right = 42, width - 42
    green, red = "#2f855a", "#c2413b"
    elements = [
        f'<svg viewBox="0 0 {width} {height}" role="img" '
        f'aria-label="{escape(gene_a["name"])} {escape(gene_b["name"])} fusion model" '
        'style="width:100%;height:auto;display:block;background:#fff">',
        _svg_text(
            left,
            34,
            f'{gene_a["name"]} - {gene_b["name"]}',
            24,
            weight=750,
        ),
        _svg_text(
            left,
            61,
            _candidate_summary(candidate, bundle),
            15,
            color="#52606d",
        ),
        '<line x1="0" y1="82" x2="100%" y2="82" stroke="#d8e0e5"/>',
    ]

    _draw_gene_track(elements, gene_a, left, right, 145, green)
    _draw_gene_track(elements, gene_b, left, right, 235, red)
    _draw_observed_join(elements, candidate, left, right, 145, 235)
    _draw_breakpoint_support(elements, candidate.get("breakpointSupport"), left, 292)

    elements.extend(
        [
            _svg_text(left, 346, "Predicted fusion product models", 18, weight=750),
            _svg_text(
                left,
                369,
                "Exon-scale candidate molecules; spacing is compressed and not genomic scale.",
                13,
                color="#667580",
            ),
        ]
    )
    if not models:
        elements.append(
            _svg_text(
                left,
                414,
                "No fusion product model could be resolved from these reads.",
                14,
                color="#667580",
            )
        )
    for index, model in enumerate(models):
        _draw_product_model(
            elements,
            model,
            candidate,
            left,
            right,
            420 + index * 108,
            green,
            red,
        )

    elements.append("</svg>")
    return "".join(elements)


def _candidate_summary(candidate: dict[str, Any], bundle: dict[str, Any]) -> str:
    support = int(candidate.get("support") or 0)
    model_count = len(bundle.get("models") or [])
    resolved = int(bundle.get("resolvedReads") or 0)
    breakpoint = (candidate.get("breakpointSupport") or {}).get("total") or {}
    fraction = breakpoint.get("fraction")
    support_text = ""
    if breakpoint.get("total"):
        percentage = round(float(fraction or 0) * 100)
        support_text = (
            f'; breakpoint-local support {breakpoint.get("support", 0)}/'
            f'{breakpoint["total"]} ({percentage}%)'
        )
    return (
        f"{support} split-read support; {model_count} fusion-forward model"
        f'{"s" if model_count != 1 else ""}; {resolved}/{support} reads resolved'
        f"{support_text}"
    )


def _draw_gene_track(
    output: list[str],
    gene: dict[str, Any],
    x1: float,
    x2: float,
    y: float,
    color: str,
) -> None:
    output.append(
        _svg_text(
            x1,
            y - 42,
            f'{gene["name"]}: {gene["chrom"]} {gene.get("strand", ".")} strand',
            17,
            weight=750,
        )
    )
    output.append(
        _svg_text(
            x1,
            y - 20,
            f'{gene["chrom"]}:{int(gene["start"]):,}-{int(gene["end"]):,}',
            12,
            color="#667580",
        )
    )
    output.append(
        f'<line x1="{x1}" y1="{y}" x2="{x2}" y2="{y}" '
        f'stroke="{color}" stroke-width="4"/>'
    )
    direction = 1 if gene.get("strand") != "-" else -1
    for x in range(int(x1 + 65), int(x2 - 20), 85):
        output.append(_svg_chevron(x, y, direction, color))
    for exon in gene.get("exons") or []:
        start = _scale(int(exon["start"]), gene, x1, x2)
        end = _scale(int(exon["end"]), gene, x1, x2)
        exon_x, exon_width = min(start, end), max(abs(end - start), 8)
        label = f'E{exon.get("number")}' if exon.get("number") else ""
        output.append(
            f'<rect x="{exon_x:.1f}" y="{y - 15}" width="{exon_width:.1f}" '
            f'height="30" rx="2" fill="{color}" stroke="#fff" stroke-width="2"/>'
        )
        if label and exon_width > 24:
            output.append(
                _svg_text(
                    exon_x + exon_width / 2,
                    y + 5,
                    label,
                    11,
                    color="#fff",
                    weight=700,
                    anchor="middle",
                )
            )


def _draw_observed_join(
    output: list[str],
    candidate: dict[str, Any],
    x1: float,
    x2: float,
    y1: float,
    y2: float,
) -> None:
    support = candidate.get("breakpointSupport") or {}
    points = support.get("breakpoints") or []
    if len(points) < 2:
        return
    genes = {
        candidate["geneA"]["id"]: candidate["geneA"],
        candidate["geneB"]["id"]: candidate["geneB"],
    }
    coords = []
    for point in points[:2]:
        gene = genes.get((point.get("gene") or {}).get("id"))
        if not gene:
            continue
        coords.append(
            (
                _scale(int(point.get("position") or gene["start"]), gene, x1, x2),
                y1 if gene["id"] == candidate["geneA"]["id"] else y2,
            )
        )
    if len(coords) != 2:
        return
    output.append(
        f'<line x1="{coords[0][0]:.1f}" y1="{coords[0][1] + 18}" '
        f'x2="{coords[1][0]:.1f}" y2="{coords[1][1] - 18}" '
        'stroke="#71808a" stroke-width="2" stroke-dasharray="7 6"/>'
    )
    output.append(
        _svg_text(
            (coords[0][0] + coords[1][0]) / 2,
            (coords[0][1] + coords[1][1]) / 2,
            "observed join",
            12,
            color="#52606d",
            anchor="middle",
        )
    )


def _draw_breakpoint_support(
    output: list[str], support: dict[str, Any] | None, x: float, y: float
) -> None:
    if not support or not support.get("total"):
        return
    total = support["total"]
    output.append(_svg_text(x, y, "Breakpoint-local read support", 14, weight=750))
    output.append(
        _svg_text(
            x,
            y + 22,
            _support_label(total) + " total reads spanning inferred breakpoint loci",
            12,
            color="#667580",
        )
    )
    output.append(
        _svg_text(
            x + 410,
            y + 22,
            "5'→3': " + _support_label(support.get("fiveToThree") or {}),
            12,
            color="#667580",
        )
    )
    output.append(
        _svg_text(
            x + 750,
            y + 22,
            "3'→5': " + _support_label(support.get("threeToFive") or {}),
            12,
            color="#667580",
        )
    )


def _draw_product_model(
    output: list[str],
    model: dict[str, Any],
    candidate: dict[str, Any],
    x1: float,
    x2: float,
    y: float,
    green: str,
    red: str,
) -> None:
    output.append(
        _svg_text(
            x1,
            y,
            f'{model.get("label", "Fusion model")} ({model.get("support", 0)} reads)',
            14,
            weight=750,
        )
    )
    output.append(_svg_text(x1 + 285, y + 25, "5′", 12, color="#667580"))
    output.append(_svg_text(x2, y + 25, "3′", 12, color="#667580", anchor="end"))
    parts = list(model.get("parts") or [])
    if not parts:
        return
    track_x = x1 + 285
    gap = 18
    available = x2 - track_x - gap * (len(parts) - 1)
    spans = [max(_part_span(part), 1) for part in parts]
    total_span = sum(spans)
    cursor = track_x
    for index, (part, span) in enumerate(zip(parts, spans)):
        width = max(80, available * span / total_span)
        gene_id = (part.get("gene") or {}).get("id")
        color = green if gene_id == candidate["geneA"]["id"] else red
        exons = (
            part.get("productExons")
            or part.get("retainedExons")
            or part.get("exons")
            or []
        )
        if not exons:
            output.append(
                f'<rect x="{cursor:.1f}" y="{y + 39}" width="{width:.1f}" '
                f'height="28" rx="2" fill="{color}"/>'
            )
        else:
            exon_spans = [
                max(int(exon["end"]) - int(exon["start"]) + 1, 1) for exon in exons
            ]
            exon_gap = 5
            exon_available = width - exon_gap * (len(exons) - 1)
            exon_cursor = cursor
            for exon, exon_span in zip(exons, exon_spans):
                exon_width = max(12, exon_available * exon_span / sum(exon_spans))
                output.append(
                    f'<rect x="{exon_cursor:.1f}" y="{y + 39}" '
                    f'width="{exon_width:.1f}" height="28" rx="2" '
                    f'fill="{color}" stroke="#fff" stroke-width="1.5"/>'
                )
                label = f'E{exon.get("number")}' if exon.get("number") else ""
                if label and exon_width > 24:
                    output.append(
                        _svg_text(
                            exon_cursor + exon_width / 2,
                            y + 58,
                            label,
                            10,
                            color="#fff",
                            weight=700,
                            anchor="middle",
                        )
                    )
                exon_cursor += exon_width + exon_gap
        output.append(
            _svg_text(
                cursor,
                y + 88,
                f'{(part.get("gene") or {}).get("name", "gene")} '
                f'gene-{"forward" if part.get("orientation") == "+" else "reverse"}',
                11,
                color="#52606d",
            )
        )
        cursor += width + gap


def _scale(position: int, gene: dict[str, Any], x1: float, x2: float) -> float:
    start, end = int(gene["start"]), int(gene["end"])
    if end <= start:
        return x1
    return x1 + (position - start) / (end - start) * (x2 - x1)


def _part_span(part: dict[str, Any]) -> int:
    exons = (
        part.get("productExons") or part.get("retainedExons") or part.get("exons") or []
    )
    if exons:
        return sum(max(int(exon["end"]) - int(exon["start"]) + 1, 1) for exon in exons)
    return max(int(part.get("end") or 0) - int(part.get("start") or 0), 1)


def _support_label(bucket: dict[str, Any]) -> str:
    total = int(bucket.get("total") or 0)
    support = int(bucket.get("support") or 0)
    percentage = round(float(bucket.get("fraction") or 0) * 100) if total else 0
    return f"{support}/{total} ({percentage}%)"


def _svg_chevron(x: float, y: float, direction: int, color: str) -> str:
    offset = 8 * direction
    return (
        f'<path d="M {x - offset} {y - 8} L {x} {y} '
        f'L {x - offset} {y + 8}" fill="none" stroke="{color}" '
        'stroke-width="3" stroke-linecap="round"/>'
    )


def _svg_text(
    x: float,
    y: float,
    text: str,
    size: int,
    *,
    color: str = "#18242c",
    weight: int = 500,
    anchor: str = "start",
) -> str:
    return (
        f'<text x="{x}" y="{y}" fill="{color}" font-size="{size}" '
        f'font-weight="{weight}" text-anchor="{anchor}" '
        'font-family="Inter,Arial,sans-serif">'
        f"{escape(str(text))}</text>"
    )
