"""Shared presentation rules for coverage metrics."""

from __future__ import annotations

import math
from typing import Any, Dict

# Machine-readable on/off-target measurements written to the sample-tracking TSV.
COVERAGE_READ_LENGTH_EXPORT_FIELDS: tuple[str, ...] = (
    "on_target_reads",
    "off_target_reads",
    "on_target_bases",
    "off_target_bases",
    "on_target_read_percent",
    "on_target_base_percent",
    "mean_on_target_length",
    "mean_off_target_length",
    "median_on_target_length",
    "median_off_target_length",
)


def _export_number(value: Any, *, digits: int | None = None) -> Any:
    """Format a measurement for TSV export (plain number, empty if unknown)."""
    if value is None or value == "":
        return ""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(number):
        return ""
    if digits is None:
        return int(round(number))
    return round(number, digits)


def coverage_read_length_export_fields(
    stats: Dict[str, Any] | None,
) -> Dict[str, Any]:
    """Numeric on/off-target length measurements for sample-tracking export."""
    empty = {key: "" for key in COVERAGE_READ_LENGTH_EXPORT_FIELDS}
    if not stats or not stats.get("available"):
        return empty

    on_reads = _export_number(stats.get("on_target_reads"))
    off_reads = _export_number(stats.get("off_target_reads"))
    on_bases = _export_number(stats.get("on_target_bases"))
    off_bases = _export_number(stats.get("off_target_bases"))
    on_reads_n = on_reads if isinstance(on_reads, int) else 0
    off_reads_n = off_reads if isinstance(off_reads, int) else 0
    on_bases_n = on_bases if isinstance(on_bases, int) else 0
    off_bases_n = off_bases if isinstance(off_bases, int) else 0
    total_reads = on_reads_n + off_reads_n
    total_bases = on_bases_n + off_bases_n
    return {
        "on_target_reads": on_reads,
        "off_target_reads": off_reads,
        "on_target_bases": on_bases,
        "off_target_bases": off_bases,
        "on_target_read_percent": (
            _export_number(100.0 * on_reads_n / total_reads, digits=2)
            if total_reads
            else ""
        ),
        "on_target_base_percent": (
            _export_number(100.0 * on_bases_n / total_bases, digits=2)
            if total_bases
            else ""
        ),
        "mean_on_target_length": _export_number(
            stats.get("mean_on_target_length"), digits=1
        ),
        "mean_off_target_length": _export_number(
            stats.get("mean_off_target_length"), digits=1
        ),
        "median_on_target_length": _export_number(
            stats.get("median_on_target_length"), digits=1
        ),
        "median_off_target_length": _export_number(
            stats.get("median_off_target_length"), digits=1
        ),
    }


def coverage_quality_name(target_coverage: float) -> str:
    """Return the coverage quality tier used throughout the GUI."""
    if target_coverage >= 30:
        return "Excellent"
    if target_coverage >= 20:
        return "Good"
    if target_coverage >= 10:
        return "Moderate"
    return "Insufficient"


def format_read_count(count: int | None) -> str:
    try:
        value = int(count or 0)
    except (TypeError, ValueError):
        return "0"
    return f"{value:,}"


def read_length_summary_line(stats: Dict[str, Any] | None) -> str:
    """One-line on/off-target median summary for the sample insight card."""
    from robin.analysis.bam_preprocessor import format_read_length

    if not stats or not stats.get("available"):
        return "On/off-target read length: not available"
    on_median = format_read_length(stats.get("median_on_target_length"))
    off_median = format_read_length(stats.get("median_off_target_length"))
    return f"On-target median {on_median} · Off-target median {off_median}"
