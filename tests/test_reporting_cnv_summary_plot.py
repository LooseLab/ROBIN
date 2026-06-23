from __future__ import annotations

import numpy as np
import pandas as pd

from robin.reporting.plotting import create_CNV_plot


class _FakeResult:
    def __init__(self, cnv: dict) -> None:
        self.cnv = cnv


def _sample_cnv() -> dict:
    return {
        "chr1": np.array([2.0, 2.1, 1.9, 2.0]),
        "chr2": np.array([2.0, 2.2, 2.1, 2.0]),
    }


def test_create_cnv_plot_ploidy_mode_returns_jpeg() -> None:
    result = _FakeResult(_sample_cnv())
    buf = create_CNV_plot(result, {"bin_width": 1_000_000})
    data = buf.getvalue()
    assert data[:2] == b"\xff\xd8"


def test_log2_linear_axis_limits_uses_percentile_with_minimum_span() -> None:
    from robin.reporting.plotting import _log2_linear_axis_limits

    # Tight cluster should still show at least +/- 2.
    tight = np.array([0.0, 0.05, -0.05, 0.1, -0.1])
    y_min, y_max = _log2_linear_axis_limits(tight)
    assert y_min <= -2.0
    assert y_max >= 2.0

    # Outliers can expand up to the upper cap.
    wide = np.array([0.0] * 95 + [2.5, -2.5, 3.0, -3.0, 0.1])
    y_min, y_max = _log2_linear_axis_limits(wide)
    assert y_max <= 4.0 + 0.01
    assert y_min >= -4.0 - 0.01
    assert y_min < 0 < y_max


def test_report_log2_from_ploidy_matches_gui() -> None:
    """Report CNV log2 mode uses the same transform as the GUI."""
    from robin.analysis.cnv_analysis import compute_cnv_log2_from_ploidy

    cnv = {"chr5": np.array([3.5, 3.5, 3.5])}
    log2 = compute_cnv_log2_from_ploidy(cnv, "Male")
    np.testing.assert_allclose(log2["chr5"], np.log2(3.5 / 2.0), rtol=1e-6)


def test_create_cnv_plot_log2_ratio_mode_returns_jpeg() -> None:
    result = _FakeResult(_sample_cnv())
    log2_ratios = {
        "chr1": np.array([0.1, -0.2, 0.0, 0.3]),
        "chr2": np.array([-0.1, 0.4, -0.3, 0.0]),
    }
    buf = create_CNV_plot(
        result,
        {"bin_width": 1_000_000},
        normalized_cnv=log2_ratios,
        use_normalized_difference=True,
        plot_bin_width=1_000_000,
    )
    data = buf.getvalue()
    assert data[:2] == b"\xff\xd8"


def test_downsample_cnv_for_plot_groups_values() -> None:
    from robin.analysis.cnv_analysis import downsample_cnv_for_plot

    values = np.array([1.0, 3.0, 5.0, 7.0], dtype=float)
    x_bp, out = downsample_cnv_for_plot(values, analysis_bin_width=12_000, plot_bin_width=24_000)
    assert len(out) == 2
    assert out[0] == 2.0
    assert out[1] == 6.0
    assert x_bp[0] == 12_000
    assert x_bp[1] == 36_000


def test_downsample_cnv_chromosome_track_keeps_full_x_axis() -> None:
    from robin.analysis.cnv_analysis import downsample_cnv_chromosome_track

    analysis_bw = 12_000
    n_bins = 1000
    values = np.linspace(0.0, 1.0, n_bins)
    x_mb, out, x_max_mb = downsample_cnv_chromosome_track(
        values, analysis_bw, plot_bin_width=500_000,
    )
    assert x_max_mb == n_bins * analysis_bw / 1_000_000
    assert len(out) < n_bins
    assert len(x_mb) == len(out)
    assert float(x_mb[-1]) < x_max_mb


def test_create_cnv_plot_per_chromosome_log2_mode_returns_jpeg() -> None:
    from robin.reporting.plotting import create_CNV_plot_per_chromosome

    sample = {
        "chr1": np.linspace(-0.2, 0.3, 40),
        "chr2": np.linspace(0.1, -0.1, 30),
    }
    log2 = {
        "chr1": np.linspace(-0.2, 0.3, 40),
        "chr2": np.linspace(0.1, -0.1, 30),
    }
    panel_genes = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start_pos": [0],
            "end_pos": [20_000_000],
            "gene": ["TEST1"],
        }
    )
    target_coverage = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "startpos": [1_000_000],
            "endpos": [2_000_000],
            "name": ["TEST1"],
            "length": [1_000_000],
            "coverage": [25.0],
            "bases": [25_000_000],
        }
    )

    class _FakeResult:
        cnv = sample

    plots = create_CNV_plot_per_chromosome(
        _FakeResult(),
        {"bin_width": 1_000_000},
        chromosomes=["chr1"],
        panel_genes_df=panel_genes,
        normalized_cnv=log2,
        target_coverage_df=target_coverage,
        use_log2_ratio=True,
    )
    assert len(plots) == 1
    assert plots[0][0] == "chr1"
    assert plots[0][1].getvalue()[:2] == b"\xff\xd8"
