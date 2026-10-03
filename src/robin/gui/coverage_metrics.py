"""Shared presentation rules for coverage metrics."""

from __future__ import annotations


def coverage_quality_name(target_coverage: float) -> str:
    """Return the coverage quality tier used throughout the GUI."""
    if target_coverage >= 30:
        return "Excellent"
    if target_coverage >= 20:
        return "Good"
    if target_coverage >= 10:
        return "Moderate"
    return "Insufficient"
