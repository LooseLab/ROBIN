from pathlib import Path

from robin.gui.components.summary import _extract_mgmt_data


def test_extract_mgmt_data_counts_final_mgmt_bed(tmp_path: Path) -> None:
    (tmp_path / "final_mgmt.csv").write_text(
        "status,pred,average\nMethylated,0.87,42.5\n"
    )
    (tmp_path / "final_mgmt.bed").write_text(
        "chr10\t1\t2\nchr10\t3\t4\nchr10\t5\t6\n"
    )

    data = _extract_mgmt_data(tmp_path)

    assert data["status"] == "Methylated"
    assert data["average_methylation"] == "42.5%"
    assert data["prediction_score"] == "0.9%"
    assert data["cpg_sites"] == 3
