"""
Shared helpers for classifying ClinVar significance from a VCF INFO field.

Both the GUI (snp_processing) and the PDF report (reporting/sections/variants)
previously had divergent implementations. This module provides one consistent
source of truth so that counts and significance flags match across surfaces.

Semantics:
- `is_pathogenic` is True for CLNSIG values that directly indicate pathogenic
  or likely pathogenic (including combined "Pathogenic/Likely_pathogenic").
- For conflicting germline classifications, CLNSIGCONF is inspected and the
  variant is called pathogenic only if pathogenic/likely-pathogenic submissions
  outnumber benign/likely-benign submissions.
- `is_oncogenic` is True for ONC values indicating oncogenic or likely
  oncogenic, with ONCCONF resolved the same way as CLNSIGCONF.
- `is_vus` is True for Uncertain_significance (VUS) on CLNSIG or ONC.
  VUS is never treated as pathogenic/oncogenic.
- `is_somatic_significant` is True for SCI Tier I (strong) or Tier II
  (potential) somatic clinical impact.
- `is_clinvar_significant` is the union of pathogenic/oncogenic/somatic
  significant tracks plus VUS (used for GUI highlighting/filtering).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Iterable, Mapping, Optional

PATHOGENIC_TERMS: tuple[str, ...] = (
    "pathogenic/likely_pathogenic",
    "pathogenic/likely pathogenic",
    "likely_pathogenic",
    "likely pathogenic",
    "pathogenic",
)

BENIGN_TERMS: tuple[str, ...] = (
    "benign/likely_benign",
    "benign/likely benign",
    "likely_benign",
    "likely benign",
    "benign",
)

CONFLICTING_GERMLINE_TERMS: tuple[str, ...] = (
    "conflicting_classifications_of_pathogenicity",
    "conflicting_interpretations_of_pathogenicity",
    "conflicting classifications of pathogenicity",
    "conflicting interpretations of pathogenicity",
)

ONCOGENIC_TERMS: tuple[str, ...] = (
    "oncogenic/likely_oncogenic",
    "oncogenic/likely oncogenic",
    "likely_oncogenic",
    "likely oncogenic",
    "oncogenic",
)

BENIGN_ONC_TERMS: tuple[str, ...] = (
    "benign/likely_benign",
    "benign/likely benign",
    "likely_benign",
    "likely benign",
    "benign",
)

CONFLICTING_ONCOGENIC_TERMS: tuple[str, ...] = (
    "conflicting_classifications_of_oncogenicity",
    "conflicting_interpretations_of_oncogenicity",
    "conflicting classifications of oncogenicity",
    "conflicting interpretations of oncogenicity",
)

VUS_TERMS: tuple[str, ...] = (
    "uncertain_significance",
    "uncertain significance",
)

SCI_SIGNIFICANT_TERMS: tuple[str, ...] = (
    "tier_i_-_strong",
    "tier_ii_-_potential",
)


@dataclass(frozen=True)
class ClinVarSignificance:
    """Structured ClinVar significance across germline, oncogenic, and somatic tracks."""

    is_clinvar_significant: bool
    is_pathogenic: bool
    is_oncogenic: bool
    is_vus: bool
    is_somatic_significant: bool
    has_conflicting_germline: bool
    has_conflicting_oncogenic: bool
    raw_clnsig: str
    raw_clnsigconf: str
    raw_onc: str
    raw_oncconf: str
    raw_sci: str

    def as_display_dict(self) -> dict:
        return {
            "is_clinvar_significant": bool(self.is_clinvar_significant),
            "is_pathogenic": bool(self.is_pathogenic),
            "is_oncogenic": bool(self.is_oncogenic),
            "is_vus": bool(self.is_vus),
            "is_somatic_significant": bool(self.is_somatic_significant),
            "has_conflicting_germline": bool(self.has_conflicting_germline),
            "has_conflicting_oncogenic": bool(self.has_conflicting_oncogenic),
            "CLNSIG": self.raw_clnsig,
            "CLNSIGCONF": self.raw_clnsigconf,
            "ONC": self.raw_onc,
            "ONCCONF": self.raw_oncconf,
            "SCI": self.raw_sci,
        }


@dataclass(frozen=True)
class PathogenicityClassification:
    """Germline-only ClinVar pathogenicity result (backward compatible)."""

    is_pathogenic: bool
    has_conflicting: bool
    raw_clnsig: str
    raw_clnsigconf: str

    def as_display_dict(self) -> dict:
        return {
            "is_pathogenic": bool(self.is_pathogenic),
            "has_conflicting_classifications": bool(self.has_conflicting),
            "CLNSIG": self.raw_clnsig,
            "CLNSIGCONF": self.raw_clnsigconf,
        }


def _extract_info_field(info_str: str, key: str) -> str:
    """Return the raw string value for `key=...` in a VCF INFO string."""
    if not info_str or "=" not in info_str:
        return ""
    target = key + "="
    for field in info_str.split(";"):
        if field.startswith(target):
            return field[len(target) :]
    return ""


def _contains_any(haystack: str, needles: Iterable[str]) -> bool:
    low = haystack.lower()
    return any(n in low for n in needles)


_CONF_COUNT_RE = re.compile(r"^(?P<label>[^()]+?)\((?P<count>-?\d+)\)\s*$")


def _score_conflicting_conf(
    conf_value: str,
    *,
    significant_terms: tuple[str, ...],
    benign_terms: tuple[str, ...],
) -> tuple[int, int]:
    """
    Parse a CLNSIGCONF/ONCCONF field such as
    "Pathogenic(2)|Likely_pathogenic(1)|Benign(3)" into (significant, benign)
    submission counts.
    """
    significant = 0
    benign = 0
    if not conf_value:
        return 0, 0

    for raw_part in conf_value.split("|"):
        part = raw_part.strip()
        if not part:
            continue
        match = _CONF_COUNT_RE.match(part)
        if match:
            label = match.group("label").strip().lower()
            try:
                count = int(match.group("count"))
            except (TypeError, ValueError):
                count = 0
        else:
            label = part.lower()
            count = 1
        if count <= 0:
            continue

        if _contains_any(label, significant_terms):
            significant += count
        elif _contains_any(label, benign_terms):
            benign += count

    return significant, benign


def _classify_track(
    raw_value: str,
    raw_conf: str,
    *,
    significant_terms: tuple[str, ...],
    benign_terms: tuple[str, ...],
    conflicting_terms: tuple[str, ...],
) -> tuple[bool, bool]:
    """Return (is_significant, has_conflicting) for a single ClinVar track."""
    if not raw_value:
        return False, False

    value_lc = raw_value.lower()
    has_conflicting = any(term in value_lc for term in conflicting_terms)

    if has_conflicting:
        if raw_conf:
            significant_count, benign_count = _score_conflicting_conf(
                raw_conf,
                significant_terms=significant_terms,
                benign_terms=benign_terms,
            )
            if significant_count > benign_count:
                return True, True
        return False, True

    if _contains_any(value_lc, significant_terms):
        return True, False

    return False, False


def _classify_somatic_sci(raw_sci: str) -> bool:
    if not raw_sci:
        return False
    return _contains_any(raw_sci.lower(), SCI_SIGNIFICANT_TERMS)


def _is_vus_value(raw_value: str) -> bool:
    """Return True when a ClinVar track value is Uncertain_significance (VUS)."""
    if not raw_value:
        return False
    return _contains_any(raw_value.lower(), VUS_TERMS)


def _info_mapping_to_string(info: Any) -> str:
    """Serialize a pysam-style INFO mapping into a VCF INFO string."""
    if not info:
        return ""
    if isinstance(info, str):
        return info

    parts: list[str] = []
    if isinstance(info, Mapping):
        items = info.items()
    else:
        return ""

    for key, value in items:
        key_str = str(key)
        if value is True:
            parts.append(key_str)
            continue
        if isinstance(value, (list, tuple)):
            value_str = "|".join(str(item) for item in value)
        else:
            value_str = str(value)
        parts.append(f"{key_str}={value_str}")
    return ";".join(parts)


def classify_clinvar_significance_from_mapping(
    info: Any,
) -> ClinVarSignificance:
    """Classify significance from a pysam-style INFO mapping."""
    return classify_clinvar_significance(_info_mapping_to_string(info))


def classify_clinvar_significance(
    info_str: Optional[str],
) -> ClinVarSignificance:
    """
    Classify a VCF record from its INFO string into a ClinVarSignificance.

    Accepts the raw `INFO` text (semicolon-separated, as emitted by VCF 4.x).
    """
    if not info_str:
        return ClinVarSignificance(
            False, False, False, False, False, False, False, "", "", "", "", ""
        )

    clnsig = _extract_info_field(info_str, "CLNSIG")
    clnsigconf = _extract_info_field(info_str, "CLNSIGCONF")
    onc = _extract_info_field(info_str, "ONC")
    oncconf = _extract_info_field(info_str, "ONCCONF")
    sci = _extract_info_field(info_str, "SCI")

    is_pathogenic, has_conflicting_germline = _classify_track(
        clnsig,
        clnsigconf,
        significant_terms=PATHOGENIC_TERMS,
        benign_terms=BENIGN_TERMS,
        conflicting_terms=CONFLICTING_GERMLINE_TERMS,
    )
    is_oncogenic, has_conflicting_oncogenic = _classify_track(
        onc,
        oncconf,
        significant_terms=ONCOGENIC_TERMS,
        benign_terms=BENIGN_ONC_TERMS,
        conflicting_terms=CONFLICTING_ONCOGENIC_TERMS,
    )
    is_somatic_significant = _classify_somatic_sci(sci)

    # VUS is highlightable but never collapses into pathogenic/oncogenic.
    is_vus = (not is_pathogenic and _is_vus_value(clnsig)) or (
        not is_oncogenic and _is_vus_value(onc)
    )

    is_clinvar_significant = (
        is_pathogenic or is_oncogenic or is_somatic_significant or is_vus
    )

    return ClinVarSignificance(
        is_clinvar_significant=is_clinvar_significant,
        is_pathogenic=is_pathogenic,
        is_oncogenic=is_oncogenic,
        is_vus=is_vus,
        is_somatic_significant=is_somatic_significant,
        has_conflicting_germline=has_conflicting_germline,
        has_conflicting_oncogenic=has_conflicting_oncogenic,
        raw_clnsig=clnsig,
        raw_clnsigconf=clnsigconf,
        raw_onc=onc,
        raw_oncconf=oncconf,
        raw_sci=sci,
    )


def classify_clinvar(info_str: Optional[str]) -> PathogenicityClassification:
    """Germline-only classification (backward compatible wrapper)."""
    result = classify_clinvar_significance(info_str)
    return PathogenicityClassification(
        is_pathogenic=result.is_pathogenic,
        has_conflicting=result.has_conflicting_germline,
        raw_clnsig=result.raw_clnsig,
        raw_clnsigconf=result.raw_clnsigconf,
    )


def is_pathogenic_from_info(info_str: Optional[str]) -> bool:
    """Convenience wrapper returning only the germline `is_pathogenic` boolean."""
    return classify_clinvar(info_str).is_pathogenic


def is_clinvar_significant_from_info(info_str: Optional[str]) -> bool:
    """Convenience wrapper returning the combined ClinVar significance boolean."""
    return classify_clinvar_significance(info_str).is_clinvar_significant


def is_clinvar_significant_from_mapping(info: Any) -> bool:
    """Convenience wrapper for pysam-style INFO mappings."""
    return classify_clinvar_significance_from_mapping(info).is_clinvar_significant
