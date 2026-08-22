"""Tests for Tucan asset helpers and scores CSV formatting (no model download)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from robin.analysis.tucan_analysis import (
    SCORE_META_COLUMNS,
    _top_prediction,
    append_tucan_scores,
    binarize_methylation_calls,
)
from robin.reporting.sections.classification import drop_classifier_score_meta_columns
from robin.utils.tucan_manager import (
    DEFAULT_NUM_CPGS,
    DEFAULT_PROBE_MARGIN,
    get_default_num_cpgs,
    get_probe_margin,
    write_mapping_probes_bed,
)


def test_write_mapping_probes_bed(tmp_path: Path) -> None:
    tucan_probe = tmp_path / "probe.bed"
    tucan_probe.write_text(
        "chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n"
        "chr1\t10\t11\tcg0001\t0\t+\n"
        "chr2\t20\t21\tcg0002\t1\t-\n",
        encoding="utf-8",
    )
    out = tmp_path / "mapping.bed"
    write_mapping_probes_bed(tucan_probe, out)
    mapped = pd.read_csv(out, sep="\t")
    assert list(mapped.columns) == ["chr", "start", "end", "ID_REF"]
    assert str(mapped.iloc[0]["chr"]) == "1"
    assert int(mapped.iloc[0]["start"]) == 10
    assert int(mapped.iloc[0]["end"]) == 11
    assert mapped.iloc[0]["ID_REF"] == "cg0001"
    assert str(mapped.iloc[1]["chr"]) == "2"
    assert mapped.iloc[1]["ID_REF"] == "cg0002"


def test_binarize_methylation_calls() -> None:
    df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2"],
            "chromStart": [1, 2, 3],
            "chromEnd": [2, 3, 4],
            "methylation_call": [0, 3, 1],
            "probe_id": ["a", "b", "c"],
        }
    )
    out = binarize_methylation_calls(df)
    assert out["methylation_call"].tolist() == [0, 1, 1]


def test_append_tucan_scores_renames_probes(tmp_path: Path) -> None:
    scores = tmp_path / "tucan_scores.csv"
    prediction = pd.DataFrame(
        [
            {
                "ACC": 0.1,
                "SCHW": 0.9,
                "probes": 10000,
            }
        ]
    )
    append_tucan_scores(str(scores), prediction, timestamp_ms=12345.0)
    written = pd.read_csv(scores)
    assert "timestamp" in written.columns
    assert "number_probes" in written.columns
    assert "covered_cpgs" in written.columns
    assert written.iloc[0]["SCHW"] == pytest.approx(0.9)
    assert int(written.iloc[0]["number_probes"]) == 10000
    assert int(written.iloc[0]["timestamp"]) == 12345


def test_top_prediction_skips_meta() -> None:
    df = pd.DataFrame(
        [
            {
                "timestamp": 1,
                "number_probes": 50,
                "covered_cpgs": 50,
                "probes": 50,
                "A": 0.2,
                "B": 0.8,
            }
        ]
    )
    top_class, top_score, covered = _top_prediction(df)
    assert top_class == "B"
    assert top_score == pytest.approx(0.8)
    assert covered == 50
    assert "probes" in SCORE_META_COLUMNS


def test_drop_meta_includes_probes() -> None:
    df = pd.DataFrame(
        [
            {
                "timestamp": 1,
                "probes": 10,
                "number_probes": 10,
                "SCHW": 0.95,
            }
        ]
    )
    cleaned = drop_classifier_score_meta_columns(df)
    assert list(cleaned.columns) == ["SCHW"]


def test_default_num_cpgs(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("ROBIN_TUCAN_NUM_CPGS", raising=False)
    assert get_default_num_cpgs() == DEFAULT_NUM_CPGS
    monkeypatch.setenv("ROBIN_TUCAN_NUM_CPGS", "15000")
    assert get_default_num_cpgs() == 15000


def test_probe_margin_defaults_to_exact_match(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("ROBIN_TUCAN_PROBE_MARGIN", raising=False)
    assert DEFAULT_PROBE_MARGIN == 0
    assert get_probe_margin() == 0
    monkeypatch.setenv("ROBIN_TUCAN_PROBE_MARGIN", "25")
    assert get_probe_margin() == 25
    monkeypatch.setenv("ROBIN_TUCAN_PROBE_MARGIN", "-1")
    with pytest.raises(ValueError, match="must be >= 0"):
        get_probe_margin()
