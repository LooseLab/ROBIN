"""Tests for log2-based arm/whole-chromosome CNV calling."""

from __future__ import annotations

import math

import numpy as np

from robin.analysis.cnv_analysis import (
    CNV_CALLING_MIN_BIN_WIDTH,
    coarsen_cnv_track_values,
    prepare_cnv_calling_track,
    resolve_cnv_calling_bin_width,
    resolve_cnv_calling_track,
)
from robin.classification_config import (
    get_cnv_thresholds,
    is_arm_event,
    is_whole_chromosome_event,
)


def test_resolve_cnv_calling_track_matches_log2_ploidy() -> None:
    cnv = {"chr1": np.array([3.0, 3.0, 3.0])}
    track = resolve_cnv_calling_track(cnv, "Male")
    np.testing.assert_allclose(track["chr1"], np.log2(1.5), rtol=1e-6)


def test_cnv_thresholds_are_log2_call_levels() -> None:
    gain, loss = get_cnv_thresholds("chr1", "Male")
    assert gain == 0.3
    assert loss == -0.3


def test_whole_chromosome_gain_on_log2_trisomy_signal() -> None:
    gain_thr, loss_thr = get_cnv_thresholds("chr1", "Male")
    log2_trisomy = math.log2(3 / 2)
    is_event, event_type = is_whole_chromosome_event(
        log2_trisomy,
        log2_trisomy,
        0.85,
        0.0,
        0.85,
        0.0,
        gain_thr,
        loss_thr,
    )
    assert is_event is True
    assert event_type == "GAIN"


def test_whole_chromosome_gain_on_moderate_visible_shift() -> None:
    """~0.5 log2 shift on the plot should call with ±0.3 thresholds."""
    gain_thr, loss_thr = get_cnv_thresholds("chr1", "Male")
    is_event, event_type = is_whole_chromosome_event(
        0.5,
        0.5,
        1.0,
        0.0,
        1.0,
        0.0,
        gain_thr,
        loss_thr,
    )
    assert is_event is True
    assert event_type == "GAIN"


def test_whole_chromosome_loss_on_moderate_visible_shift() -> None:
    gain_thr, loss_thr = get_cnv_thresholds("chr1", "Male")
    is_event, event_type = is_whole_chromosome_event(
        -0.5,
        -0.5,
        0.0,
        1.0,
        0.0,
        1.0,
        gain_thr,
        loss_thr,
    )
    assert is_event is True
    assert event_type == "LOSS"


def test_arm_gain_requires_directional_proportion() -> None:
    gain_thr, loss_thr = get_cnv_thresholds("chr1", "Male")
    is_event, event_type = is_arm_event(
        math.log2(3 / 2),
        0.5,
        0.0,
        gain_thr,
        loss_thr,
    )
    assert is_event is True
    assert event_type == "GAIN"


def test_arm_event_not_called_below_proportion_threshold() -> None:
    gain_thr, loss_thr = get_cnv_thresholds("chr1", "Male")
    is_event, _ = is_arm_event(math.log2(3 / 2), 0.3, 0.0, gain_thr, loss_thr)
    assert is_event is False


def test_resolve_cnv_calling_bin_width_minimum_one_megabase() -> None:
    assert resolve_cnv_calling_bin_width(100_000) == CNV_CALLING_MIN_BIN_WIDTH
    assert resolve_cnv_calling_bin_width(2_000_000) == 2_000_000


def test_coarsen_cnv_track_values_averages_to_one_megabase() -> None:
    values = np.arange(20, dtype=float)
    coarsened = coarsen_cnv_track_values(values, 100_000, 1_000_000)
    assert len(coarsened) == 2
    np.testing.assert_allclose(coarsened[0], np.mean(values[:10]))
    np.testing.assert_allclose(coarsened[1], np.mean(values[10:20]))


def test_prepare_cnv_calling_track_coarsens_log2_track() -> None:
    cnv = {"chr1": np.full(20, 3.0)}
    track, calling_bw = prepare_cnv_calling_track(cnv, 100_000, "Male")
    assert calling_bw == CNV_CALLING_MIN_BIN_WIDTH
    assert len(track["chr1"]) == 2
    np.testing.assert_allclose(track["chr1"], np.log2(1.5), rtol=1e-6)


def test_whole_chromosome_event_suppresses_arm_events() -> None:
    """Arm GAIN/LOSS rows should not appear when WHOLE_CHR_* is called."""
    import pandas as pd

    from robin.analysis.cnv_classification import detect_cnv_events

    bin_width = 1_000_000
    n_bins = 100
    cytobands = pd.DataFrame(
        [
            {"chrom": "chr7", "name": "p22", "start_pos": 0, "end_pos": 40_000_000},
            {"chrom": "chr7", "name": "q31", "start_pos": 40_000_000, "end_pos": 100_000_000},
        ]
    )
    cnv_data = {"chr7": np.full(n_bins, 0.5, dtype=float)}

    events = detect_cnv_events(
        cnv_data=cnv_data,
        bin_width=bin_width,
        sex_estimate="Male",
        cytobands_df=cytobands,
    )

    assert len(events) == 1
    assert events[0].event_type == "WHOLE_CHR_GAIN"
    assert events[0].chromosome == "chr7"


def test_detect_cnv_events_for_sample_reads_on_disk_outputs(tmp_path) -> None:
    import pickle

    from robin.analysis.cnv_classification import (
        detect_cnv_events_for_sample,
        format_cnv_events_card_lines,
    )

    cnv_map = {"chr7": np.full(100, 3.0, dtype=float)}
    np.save(tmp_path / "CNV.npy", cnv_map)
    np.save(
        tmp_path / "CNV_dict.npy",
        {"bin_width": 1_000_000, "variance": 1.0},
    )
    with (tmp_path / "XYestimate.pkl").open("wb") as handle:
        pickle.dump("Female", handle)

    events = detect_cnv_events_for_sample(tmp_path)
    whole_text, arm_text = format_cnv_events_card_lines(events)

    assert len(events) == 1
    assert events[0].event_type == "WHOLE_CHR_GAIN"
    assert "chr7 GAIN" in whole_text
    assert arm_text.startswith("Arm-level: none detected")
