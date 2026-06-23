from __future__ import annotations

import numpy as np

from robin.analysis.cnv_analysis import (
    compute_cnv_log2_from_ploidy,
    expected_cnv_ploidy_baseline,
)
from robin.gui.components.cnv import (
    _build_cnv_track_scatter_series,
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
