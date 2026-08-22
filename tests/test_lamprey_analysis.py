"""Tests for Lamprey integration helpers (no HuggingFace download required)."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from robin.analysis.lamprey_analysis import (
    LampreyPredictor,
    append_lamprey_scores,
    bedmethyl_to_probe_calls,
    build_feature_vector,
    load_probe_names,
    load_probe_position_map,
)
from robin.classification_config import CLASSIFIER_CONFIDENCE_THRESHOLDS
from robin.utils.lamprey_manager import (
    DEFAULT_GENOME_BUILD,
    LampreyResearchAckRequired,
    get_feature_order_bed_path,
    get_probes_bed_path,
    require_research_ack,
    research_ack_given,
)


class FakeSession:
    def __init__(self, logits):
        self.logits = np.asarray(logits, dtype=np.float32)
        self.last_feeds = None

    def run(self, output_names, input_feed):
        self.last_feeds = input_feed
        return [self.logits]


def test_sturgeon_matched_confidence_tiers():
    assert (
        CLASSIFIER_CONFIDENCE_THRESHOLDS["lamprey"]
        == CLASSIFIER_CONFIDENCE_THRESHOLDS["sturgeon"]
    )


def test_research_ack_gate(monkeypatch):
    monkeypatch.delenv("ROBIN_LAMPREY_RESEARCH_ACK", raising=False)
    assert research_ack_given() is False
    with pytest.raises(LampreyResearchAckRequired):
        require_research_ack()
    monkeypatch.setenv("ROBIN_LAMPREY_RESEARCH_ACK", "1")
    assert research_ack_given() is True
    require_research_ack()


def test_hg38_only_probe_bed():
    path = get_probes_bed_path(DEFAULT_GENOME_BUILD)
    assert path.name == "probe_hg38.bed"
    with pytest.raises(ValueError, match="hg38"):
        get_probes_bed_path("hg19")


def test_probe_assets_from_installed_lamprey():
    probes = load_probe_names(get_feature_order_bed_path())
    assert len(probes) == 353232
    mapping = load_probe_position_map(get_probes_bed_path("hg38"))
    assert mapping[("chr1", 69590)] == "cg21870274"
    assert mapping[("chr1", 69591)] == "cg21870274"


def test_bedmethyl_to_probe_calls_binarization():
    position_map = {
        ("chr1", 10): "cg1",
        ("chr1", 11): "cg1",
        ("chr2", 20): "cg2",
        ("chr3", 30): "cg3",
    }
    df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2", "chr3"],
            "chromStart": [10, 11, 20, 30],
            "percent_modified": [80.0, 90.0, 10.0, 50.0],
        }
    )
    calls = bedmethyl_to_probe_calls(df, position_map)
    assert calls["cg1"] == 1
    assert calls["cg2"] == -1
    assert "cg3" not in calls  # exactly 50% dropped


def test_build_feature_vector_and_predict(tmp_path: Path):
    probe_names = ["cg1", "cg2", "cg3"]
    vector, n_used = build_feature_vector(probe_names, {"cg1": 1, "cg3": -1})
    assert n_used == 2
    assert vector.tolist() == [1.0, 0.0, -1.0]

    # Fake logits favoring class_b after temperature/softmax
    session = FakeSession([[0.1, 2.0, 0.05]])
    predictor = LampreyPredictor(
        session=session,
        probe_names=probe_names,
        class_names=["class_a", "class_b", "class_c"],
        bin_centers=np.array([0.0, 100.0]),
        temps=np.array([0.0, 0.0]),
        input_name="input",
        output_name="output",
    )
    prediction = predictor.predict_from_probe_calls({"cg1": 1, "cg3": -1})
    assert prediction.top_class == "class_b"
    assert prediction.covered_cpgs == 2
    assert abs(sum(prediction.scores.values()) - 1.0) < 1e-5

    scores_path = tmp_path / "lamprey_scores.csv"
    append_lamprey_scores(str(scores_path), prediction, timestamp_ms=1.0)
    with scores_path.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    assert row["covered_cpgs"] == "2"
    assert "class_b" in row
