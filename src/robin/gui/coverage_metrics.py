"""Shared presentation rules for coverage metrics."""

from __future__ import annotations

from typing import Any, Dict


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
