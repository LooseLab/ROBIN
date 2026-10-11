from __future__ import annotations

import csv
from types import SimpleNamespace

import pandas as pd

from robin.gui_launcher import GUILauncher


def test_sample_tracking_tsv_exports_numeric_unit_fields(tmp_path, monkeypatch) -> None:
    sample_dir = tmp_path / "sampleA"
    sample_dir.mkdir()
    (sample_dir / "coverage_main.csv").write_text(
        "#rname,covbases,endpos\nchr1,728,100\n"
    )
    (sample_dir / "bed_coverage_main.csv").write_text(
        "chrom,bases,length\nchr1,391,100\n"
    )
    (sample_dir / "final_mgmt.csv").write_text(
        "status,pred,average\nMethylated,87.0,42.5\n"
    )
    (sample_dir / "final_mgmt.bed").write_text("chr10\t1\t2\n")

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
        row = next(
            csv.DictReader(
                (line for line in fh if not line.startswith("#")),
                delimiter="\t",
            )
        )

    assert row["analysis_coverage_global_coverage_x"] == "7.28"
    assert row["analysis_coverage_target_coverage_x"] == "3.91"
    assert row["analysis_coverage_enrichment_x"] == "0.54"
    assert row["analysis_mgmt_average_methylation_percent"] == "42.5"
    assert row["analysis_mgmt_prediction_score_percent"] == "87.0"
