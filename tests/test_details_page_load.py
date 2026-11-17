"""Helpers that keep the More Details page from blocking on first paint."""

from __future__ import annotations

from pathlib import Path

from robin.gui.components.coverage import load_target_gene_table_rows
from robin.gui.components.fusion import build_fusion_pairs_table_rows


def test_fusion_pairs_builder_returns_empty_without_pickle(tmp_path: Path) -> None:
    rows = build_fusion_pairs_table_rows(tmp_path, show_fusion_target=True)
    assert rows == []


def test_load_target_gene_table_rows_vectorized(tmp_path: Path) -> None:
    csv_path = tmp_path / "target_coverage.csv"
    csv_path.write_text(
        "chrom,startpos,endpos,name,coverage\n"
        "chr1,100,200,GENEA,12.5\n"
        "chr2,1000,2500,GENEB,40\n",
        encoding="utf-8",
    )
    rows, regions = load_target_gene_table_rows(csv_path)
    assert [row["name"] for row in rows] == ["GENEA", "GENEB"]
    assert rows[0]["startpos"] == "100"
    assert rows[1]["endpos"] == "2,500"
    assert regions["GENEA"] == "chr1:1-10200"
    assert rows[0]["__row_idx"] == 0
    assert rows[1]["coverage"] == 40.0


def test_load_target_gene_table_rows_computes_coverage(tmp_path: Path) -> None:
    csv_path = tmp_path / "bed_coverage_main.csv"
    csv_path.write_text(
        "chrom,startpos,endpos,name,length,bases\n"
        "chr7,1,10,EGFR,10,25\n",
        encoding="utf-8",
    )
    rows, regions = load_target_gene_table_rows(csv_path)
    assert len(rows) == 1
    assert rows[0]["name"] == "EGFR"
    assert rows[0]["coverage"] == 2.5
    assert "EGFR" in regions
