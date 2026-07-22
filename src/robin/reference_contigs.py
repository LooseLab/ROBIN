"""Shared rules for which reference contigs appear in coverage and CNV plots."""

from __future__ import annotations

import re
from typing import Iterable, List, Optional, Sequence, Union

import natsort

# Autosomes and sex chromosomes (chr1–22, X, Y); excludes chrM and alt/unplaced contigs.
CANONICAL_CONTIG_RE = re.compile(r"^chr(\d+|X|Y)$")
CANONICAL_WITH_M_CONTIG_RE = re.compile(r"^chr(\d+|X|Y|M)$")

REFERENCE_CONTIG_SCOPE_CANONICAL = "canonical_only"
REFERENCE_CONTIG_SCOPE_INCLUDE_M = "include_chrM"
REFERENCE_CONTIG_SCOPE_ALL = "all_reference_contigs"
REFERENCE_CONTIG_SCOPES = (
    REFERENCE_CONTIG_SCOPE_CANONICAL,
    REFERENCE_CONTIG_SCOPE_INCLUDE_M,
    REFERENCE_CONTIG_SCOPE_ALL,
)
DEFAULT_REFERENCE_CONTIG_SCOPE = REFERENCE_CONTIG_SCOPE_CANONICAL

REFERENCE_CONTIG_SCOPE_LABELS = {
    REFERENCE_CONTIG_SCOPE_CANONICAL: "Standard chromosomes only (chr1–22, X, Y)",
    REFERENCE_CONTIG_SCOPE_INCLUDE_M: "Standard chromosomes + chrM",
    REFERENCE_CONTIG_SCOPE_ALL: "All reference contigs (incl. unplaced / alt)",
}


def resolve_reference_contig_scope(scope: Optional[str]) -> str:
    """Return a valid contig scope key, falling back to the default."""
    if scope in REFERENCE_CONTIG_SCOPES:
        return scope
    return DEFAULT_REFERENCE_CONTIG_SCOPE


def is_visible_contig(name: str, scope: Optional[str] = None) -> bool:
    """Return True when a contig name should appear in plots for the given scope."""
    label = str(name).strip()
    if not label:
        return False
    resolved = resolve_reference_contig_scope(scope)
    if resolved == REFERENCE_CONTIG_SCOPE_ALL:
        return label.startswith("chr")
    if resolved == REFERENCE_CONTIG_SCOPE_INCLUDE_M:
        return bool(CANONICAL_WITH_M_CONTIG_RE.match(label))
    return bool(CANONICAL_CONTIG_RE.match(label))


def is_canonical_contig(name: str) -> bool:
    """Return True for standard autosomes/sex chromosomes (excludes chrM and alt contigs)."""
    return is_visible_contig(name, REFERENCE_CONTIG_SCOPE_CANONICAL)


def filter_contigs(
    names: Iterable[Union[str, object]],
    scope: Optional[str] = None,
) -> List[str]:
    """Natural-sort contig names that pass the scope filter."""
    resolved = resolve_reference_contig_scope(scope)
    filtered = [str(name) for name in names if is_visible_contig(str(name), resolved)]
    return natsort.natsorted(filtered, key=str)


def filter_dataframe_contigs(
    df,
    column: str,
    scope: Optional[str] = None,
):
    """Filter a DataFrame to rows whose ``column`` value is visible under ``scope``."""
    if df is None or getattr(df, "empty", True):
        return df
    mask = df[column].astype(str).map(lambda value: is_visible_contig(value, scope))
    return df.loc[mask].copy()
