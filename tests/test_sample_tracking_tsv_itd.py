from __future__ import annotations

import csv
from types import SimpleNamespace

import pandas as pd

from robin.gui_launcher import GUILauncher


def test_sample_tracking_tsv_includes_itd_results(tmp_path, monkeypatch) -> None:
    sample_dir = tmp_path / "sampleA"
    sample_dir.mkdir()
    pd.DataFrame(
        [
            {
                "gene": "NPM1",
                "chrom": "chr5",
                "position": 170837543,
                "length": 4,
                "support": 12,
                "coverage": 80,
                "vaf": 0.15,
                "label": "ITD",
            }
        ]
    ).to_csv(sample_dir / "itd_events.csv", index=False)
    pd.DataFrame(
        [{"gene": "NPM1", "n_events": 1, "max_support": 12, "max_vaf": 0.15}]
    ).to_csv(sample_dir / "itd_summary.csv", index=False)

    stub = SimpleNamespace(
        monitored_directory=str(tmp_path),
        samples_table=None,
        _last_samples_rows=[{"sample_id": "sampleA"}],
        _mnpflex_results_dir_for_sample=lambda *_args, **_kwargs: None,
    )
    export_dir = tmp_path / "export"
    export_dir.mkdir()
    monkeypatch.setattr("tempfile.gettempdir", lambda: str(export_dir))

    out_path = GUILauncher._build_sample_tracking_tsv_export(stub, ["sampleA"])
    assert out_path
    with open(out_path, newline="", encoding="utf-8") as fh:
        rows = list(
            csv.DictReader(
                (row for row in fh if not row.startswith("#")),
                delimiter="\t",
            )
        )

    row = rows[0]
    assert row["analysis_itd_status"] == "Events detected"
    assert row["analysis_itd_event_count"] == "1"
    assert row["analysis_itd_called_genes"] == "1"
    assert row["analysis_itd_genes"] == "NPM1"
    assert row["analysis_itd_event_list"] == (
        "NPM1|chr5:170837543|4bp|support=12|vaf=0.1500"
    )
