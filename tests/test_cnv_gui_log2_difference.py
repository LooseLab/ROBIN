from __future__ import annotations

import numpy as np

from robin.analysis.cnv_analysis import (
    compute_cnv_log2_from_ploidy,
    expected_cnv_ploidy_baseline,
)
from robin.gui.components.cnv import (
    _build_cnv_track_scatter_series,
    _cnv_genome_x_extent_bp,
    _cnv_plot_bin_bp_from_ui,
    _cnv_plot_bin_key_from_ui,
    _cnv_sample_relative_stats,
    _cnv_split_points_by_zscore,
    _recompute_cnv_log2_state,
)


def test_log2_from_ploidy_gain_three_vs_two_copies() -> None:
    log2 = compute_cnv_log2_from_ploidy({"chr1": np.array([3.0, 3.0, 3.0])}, "Male")
    np.testing.assert_allclose(log2["chr1"], np.log2(1.5), rtol=1e-6)


def test_recompute_cnv_log2_state_from_ploidy_track() -> None:
    state = {
        "cnv": {"chr5": np.array([3.5, 3.5, 3.5])},
        "xy": "Male",
    }
    _recompute_cnv_log2_state(state)
    np.testing.assert_allclose(
        state["cnv_log2"]["chr5"],
        np.log2(3.5 / 2.0),
        rtol=1e-6,
    )


def test_expected_baseline_chrx_male() -> None:
    assert expected_cnv_ploidy_baseline("chrX", "Male") == 1.0
    assert expected_cnv_ploidy_baseline("chr1", "Male") == 2.0


def test_build_cnv_track_scatter_series_filters_non_finite() -> None:
    track = {"chr1": np.array([0.0, np.nan, 1.0])}
    series = _build_cnv_track_scatter_series(
        track,
        selected="chr1",
        binw_analysis=1_000_000,
        plot_bin_width=1_000_000,
        chrom_palette=["#000"],
        filter_finite=True,
    )
    assert len(series) == 1
    ys = [pt[1] for pt in series[0]["data"]]
    assert all(np.isfinite(y) for y in ys)
    assert len(ys) == 2


def test_genome_x_extent_uses_analysis_bins_not_plot_bins() -> None:
    track = {
        "chr1": np.zeros(250),
        "chr2": np.zeros(242),
    }
    assert _cnv_genome_x_extent_bp(track, 1_000_000, "All") == 492_000_000
    assert _cnv_genome_x_extent_bp(track, 1_000_000, "chr1") == 250_000_000


def test_plot_bin_ui_value_resolves_numeric_select_index() -> None:
    """NiceGUI may expose select values as numeric indices; 1 Mb is index 2."""
    assert _cnv_plot_bin_key_from_ui(2) == "1 Mb"
    assert _cnv_plot_bin_bp_from_ui(2) == 1_000_000
    assert _cnv_plot_bin_bp_from_ui("1 Mb") == 1_000_000
    assert _cnv_plot_bin_bp_from_ui(0) is None


def test_sample_relative_stats_use_autosomes_only() -> None:
    cnv = {
        "chr1": np.array([2.0, 2.0]),
        "chr2": np.array([4.0, 4.0]),
        "chrX": np.array([20.0, 20.0]),
    }
    mean_val, std_val = _cnv_sample_relative_stats(cnv, None, use_log=False)
    assert mean_val == 3.0
    assert std_val == 1.0


def test_split_points_by_zscore_uses_sample_relative_thresholds() -> None:
    pts = [[0.0, 2.0], [1.0, 3.0], [2.0, 4.0]]
    high, low, norm = _cnv_split_points_by_zscore(pts, mean_val=3.0, std_val=1.0)
    assert low == [[0.0, 2.0]]
    assert norm == [[1.0, 3.0]]
    assert high == [[2.0, 4.0]]


def test_downsampled_genome_series_do_not_overlap_between_chromosomes() -> None:
    track = {
        "chr1": np.full(250, 2.0),
        "chr2": np.full(242, 2.0),
        "chr3": np.full(198, 2.0),
    }
    series = _build_cnv_track_scatter_series(
        track,
        selected="All",
        binw_analysis=1_000_000,
        plot_bin_width=10_000_000,
        chrom_palette=["#000", "#111", "#222"],
    )
    ranges = []
    for s in series:
        xs = [pt[0] for pt in s["data"]]
        ranges.append((min(xs), max(xs)))
    for (_, end), (start, _) in zip(ranges, ranges[1:]):
        assert start >= end
