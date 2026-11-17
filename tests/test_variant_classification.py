"""Tests for ClinVar significance classification."""

from robin.analysis.variant_classification import (
    classify_clinvar_significance,
    is_clinvar_significant_from_info,
)

H3C2_ONCOGENIC_INFO = (
    "ALLELEID=3881014;CLNHGVS=NC_000006.12:g.26031978T>A;CLNVC=single_nucleotide_variant;"
    "GENEINFO=H3C2:8358;MC=SO:0001583|missense_variant;"
    "ONC=Likely_oncogenic;ONCDN=Neoplasm;"
    "SCI=Tier_I_-_Strong;SCIDN=Diffuse_midline_glioma,_H3_K27M-mutant"
)


def test_oncogenic_and_somatic_without_clnsig_is_significant():
    result = classify_clinvar_significance(H3C2_ONCOGENIC_INFO)
    assert result.is_clinvar_significant
    assert not result.is_pathogenic
    assert result.is_oncogenic
    assert result.is_somatic_significant
    assert result.raw_onc == "Likely_oncogenic"
    assert result.raw_sci == "Tier_I_-_Strong"


def test_germline_pathogenic_is_significant():
    info = "CLNSIG=Pathogenic;CLNDN=Some_disease"
    result = classify_clinvar_significance(info)
    assert result.is_clinvar_significant
    assert result.is_pathogenic
    assert not result.is_oncogenic
    assert not result.is_somatic_significant


def test_conflicting_germline_resolved_by_clnsigconf():
    info = (
        "CLNSIG=Conflicting_classifications_of_pathogenicity;"
        "CLNSIGCONF=Pathogenic(2)|Benign(1)"
    )
    result = classify_clinvar_significance(info)
    assert result.is_clinvar_significant
    assert result.is_pathogenic
    assert result.has_conflicting_germline


def test_conflicting_germline_not_significant_when_benign_wins():
    info = (
        "CLNSIG=Conflicting_classifications_of_pathogenicity;"
        "CLNSIGCONF=Pathogenic(1)|Benign(3)"
    )
    result = classify_clinvar_significance(info)
    assert not result.is_clinvar_significant
    assert not result.is_pathogenic
    assert result.has_conflicting_germline


def test_conflicting_label_does_not_false_positive_pathogenic():
    info = "CLNSIG=Conflicting_classifications_of_pathogenicity"
    result = classify_clinvar_significance(info)
    assert not result.is_pathogenic
    assert not result.is_clinvar_significant


def test_oncogenic_only():
    info = "ONC=Oncogenic"
    assert is_clinvar_significant_from_info(info)


def test_sci_tier_ii_is_significant():
    info = "SCI=Tier_II_-_Potential"
    result = classify_clinvar_significance(info)
    assert result.is_somatic_significant
    assert result.is_clinvar_significant


def test_sci_tier_iii_not_significant():
    info = "SCI=Tier_III_-_Unknown"
    result = classify_clinvar_significance(info)
    assert not result.is_somatic_significant
    assert not result.is_clinvar_significant


def test_uncertain_onc_is_significant_vus():
    info = "ONC=Uncertain_significance"
    result = classify_clinvar_significance(info)
    assert result.is_clinvar_significant
    assert result.is_vus
    assert not result.is_pathogenic
    assert not result.is_oncogenic


def test_germline_vus_is_significant_but_not_pathogenic():
    info = "CLNSIG=Uncertain_significance;CLNDN=Some_disease"
    result = classify_clinvar_significance(info)
    assert result.is_clinvar_significant
    assert result.is_vus
    assert not result.is_pathogenic
    assert not result.is_oncogenic


def test_pathogenic_not_also_vus():
    info = "CLNSIG=Pathogenic/Uncertain_significance"
    result = classify_clinvar_significance(info)
    assert result.is_clinvar_significant
    assert result.is_pathogenic
    assert not result.is_vus


def test_empty_info_not_significant():
    assert not is_clinvar_significant_from_info("")
    assert not is_clinvar_significant_from_info(None)
