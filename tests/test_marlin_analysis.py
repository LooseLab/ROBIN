"""Tests for MARLIN asset resolution and probe mapping (no TensorFlow required)."""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
import pytest

from robin.analysis.marlin_analysis import (
    append_marlin_scores,
    bedmethyl_to_probe_values,
    load_probe_position_map,
)
from robin.analysis.marlin_python.predictor import (
    LiveMARLINPredictor,
    MARLINPredictor,
    extract_probe_names_from_rdata,
    read_class_names_from_xlsx,
)
from robin.utils.marlin_manager import (
    DEFAULT_GENOME_BUILD,
    get_annotations_path,
    get_features_path,
    get_marlin_resources_dir,
    get_probes_bed_path,
    resolve_model_path,
)


class FakeModel:
    def __init__(self, result):
        self.result = result
        self.last_batch = None

    def __call__(self, batch, training=False):
        self.last_batch = batch
        return [self.result]

    def predict(self, batch, verbose=0):
        self.last_batch = batch
        return [self.result]


def _batch_as_list(batch):
    if hasattr(batch, "tolist"):
        return batch.tolist()
    return list(batch)


def test_marlin_resources_are_packaged():
    resources = get_marlin_resources_dir()
    assert resources.is_dir()
    assert get_features_path().is_file()
    assert get_annotations_path().is_file()
    assert get_probes_bed_path(DEFAULT_GENOME_BUILD).is_file()
    assert get_probes_bed_path("hg19").is_file()
    assert get_probes_bed_path("t2t").is_file()


def test_unsupported_genome_build_raises():
    with pytest.raises(ValueError, match="Unsupported MARLIN genome build"):
        get_probes_bed_path("hg37")


def test_extract_probe_names_and_classes():
    probes = extract_probe_names_from_rdata(get_features_path())
    assert len(probes) == 357340
    assert probes[0] == "cg18478105"
    assert probes[-1] == "cg12623625"

    classes = read_class_names_from_xlsx(get_annotations_path())
    assert len(classes) == 42
    assert classes[0] == "AMKL_mixed"
    assert classes[-1] == "ZNF384-r"


def test_load_probe_position_map_hg38():
    mapping = load_probe_position_map(get_probes_bed_path("hg38"))
    assert len(mapping) > 100000
    # Spot-check first probe rows from the hg38 BED.
    assert mapping[("chr1", 69590)] == "cg21870274"
    assert mapping[("chr1", 69591)] == "cg21870274"


def test_bedmethyl_to_probe_values_averages_strands():
    position_map = {
        ("chr1", 10): "cg1",
        ("chr1", 11): "cg1",
        ("chr2", 20): "cg2",
    }
    df = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr2"],
            "chromStart": [10, 11, 20],
            "percent_modified": [80.0, 60.0, 10.0],
        }
    )
    values = bedmethyl_to_probe_values(df, position_map)
    assert values["cg1"] == pytest.approx(0.7)
    assert values["cg2"] == pytest.approx(0.1)


def test_predict_from_probe_values_binarization():
    model = FakeModel([0.2, 0.8])
    predictor = MARLINPredictor(
        model=model,
        probe_names=["cg1", "cg2", "cg3", "cg4"],
        class_names=["class_a", "class_b"],
    )
    prediction = predictor.predict_from_probe_values(
        {"cg1": 0.9, "cg2": 0.1, "cg4": None, "cg_missing": 0.7}
    )
    assert _batch_as_list(model.last_batch) == [[1.0, -1.0, 0.0, 0.0]]
    assert prediction.covered_cpgs == 2
    assert prediction.top_class == "class_b"
    assert prediction.top_score == pytest.approx(0.8)


def test_live_predictor_and_scores_csv(tmp_path: Path):
    model = FakeModel([0.4, 0.6])
    predictor = MARLINPredictor(
        model=model,
        probe_names=["cg1", "cg2"],
        class_names=["low", "high"],
    )
    live = LiveMARLINPredictor(predictor)
    live.update_from_probe_values({"cg1": 0.9})
    prediction = live.predict_current()

    scores_path = tmp_path / "marlin_scores.csv"
    append_marlin_scores(str(scores_path), prediction, timestamp_ms=123456.0)
    append_marlin_scores(str(scores_path), prediction, timestamp_ms=234567.0)

    with scores_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2
    assert rows[0]["timestamp"] == "123456.0"
    assert rows[0]["covered_cpgs"] == "1"
    assert rows[0]["number_probes"] == "1"
    assert float(rows[0]["high"]) == pytest.approx(0.6)


def test_resolve_model_path_uses_explicit_file(tmp_path: Path):
    model = tmp_path / "marlin_v1.model.hdf5"
    model.write_bytes(b"fake-model")
    resolved = resolve_model_path(model_path=model, download_if_missing=False)
    assert resolved == model.resolve()


def test_resolve_model_path_missing_without_download(tmp_path: Path, monkeypatch):
    missing = tmp_path / "missing.hdf5"
    monkeypatch.delenv("ROBIN_MARLIN_MODEL_PATH", raising=False)
    with pytest.raises(FileNotFoundError):
        resolve_model_path(model_path=missing, download_if_missing=False)
