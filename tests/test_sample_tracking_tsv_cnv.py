from __future__ import annotations

import csv
from types import SimpleNamespace

import numpy as np
import pandas as pd

from robin.gui_launcher import GUILauncher


def _stub(tmp_path):
    return SimpleNamespace(
        monitored_directory=str(tmp_path),
        samples_table=None,
        _last_samples_rows=[{"sample_id": "sampleA"}],
        _mnpflex_results_dir_for_sample=lambda *_args, **_kwargs: None,
    )


def _read_row(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return next(
            csv.DictReader(
                (line for line in handle if not line.startswith("#")),
                delimiter="\t",
            )
        )


def test_tracking_tsv_derives_broad_cnv_counts_without_legacy_csv(tmp_path, monkeypatch):
    sample_dir = tmp_path / "sampleA"
    sample_dir.mkdir()
    np.save(sample_dir / "CNV.npy", {"chr1": np.zeros(4)})
    np.save(sample_dir / "CNV_dict.npy", {"bin_width": 1_000_000})
    pd.DataFrame([{"counter_mapped_bases": 100}]).to_csv(
        sample_dir / "master.csv", index=False
    )
    monkeypatch.setattr(
        "robin.analysis.cnv_classification.detect_cnv_events_for_sample",
        lambda _sample_dir: [
            SimpleNamespace(event_type="WHOLE_CHR_GAIN"),
            SimpleNamespace(event_type="LOSS"),
        ],
    )
    export_dir = tmp_path / "export"
    export_dir.mkdir()
    monkeypatch.setattr("tempfile.gettempdir", lambda: str(export_dir))

    row = _read_row(
        GUILauncher._build_sample_tracking_tsv_export(_stub(tmp_path), ["sampleA"])
    )

    assert row["analysis_cnv_broad_status"] == "analysed"
    assert row["analysis_cnv_broad_whole_chromosome_events"] == "1"
    assert row["analysis_cnv_broad_arm_events"] == "1"
    assert row["analysis_cnv_broad_broad_gain_events"] == "1"
    assert row["analysis_cnv_broad_broad_loss_events"] == "1"


def test_tracking_tsv_marks_cnv_not_analyzed_without_outputs(tmp_path, monkeypatch):
    sample_dir = tmp_path / "sampleA"
    sample_dir.mkdir()
    pd.DataFrame([{"counter_mapped_bases": 100}]).to_csv(
        sample_dir / "master.csv", index=False
    )
    export_dir = tmp_path / "export"
    export_dir.mkdir()
    monkeypatch.setattr("tempfile.gettempdir", lambda: str(export_dir))

    row = _read_row(
        GUILauncher._build_sample_tracking_tsv_export(_stub(tmp_path), ["sampleA"])
    )

    assert row["analysis_cnv_broad_status"] == "not_analyzed"
    assert row["analysis_cnv_broad_whole_chromosome_events"] == ""
