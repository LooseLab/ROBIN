import pandas as pd

from robin.coverage_metrics import calculate_coverage_metrics


def test_coverage_metrics_apply_same_contig_scope_to_all_values():
    coverage = pd.DataFrame(
        {
            "#rname": ["chr1", "chrM", "chrUn_gl0001"],
            "covbases": [800, 200, 100],
            "endpos": [1000, 100, 100],
        }
    )
    targets = pd.DataFrame(
        {
            "chrom": ["chr1", "chrM", "chrUn_gl0001"],
            "bases": [400, 100, 50],
            "length": [500, 100, 100],
        }
    )

    metrics = calculate_coverage_metrics(coverage, targets, "canonical_only")

    assert metrics["global_coverage"] == 0.8
    assert metrics["target_coverage"] == 0.8
    assert metrics["enrichment"] == 1.0
    assert list(metrics["coverage_df"]["#rname"]) == ["chr1"]
    assert list(metrics["target_df"]["chrom"]) == ["chr1"]


def test_coverage_metrics_include_chr_m_when_configured():
    coverage = pd.DataFrame(
        {"#rname": ["chr1", "chrM"], "covbases": [800, 200], "endpos": [1000, 100]}
    )
    targets = pd.DataFrame(
        {"chrom": ["chr1", "chrM"], "bases": [400, 100], "length": [500, 100]}
    )

    metrics = calculate_coverage_metrics(coverage, targets, "include_chrM")

    assert metrics["global_coverage"] == 1000 / 1100
    assert metrics["target_coverage"] == 500 / 600
