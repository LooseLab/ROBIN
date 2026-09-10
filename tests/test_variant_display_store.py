"""Tests for Parquet-backed SNP/INDEL display paging."""

from __future__ import annotations

import json
from pathlib import Path

from robin.analysis.snp_processing import (
    INDEL_DISPLAY_JSON,
    INDEL_DISPLAY_PARQUET,
    SNP_DISPLAY_JSON,
    SNP_DISPLAY_PARQUET,
    VARIANT_SIDECAR_FORMAT,
    VariantDisplayStore,
    VariantTableFilters,
    build_snp_display_data,
    write_clair_variant_display_files,
    write_variant_display_bundle,
)
from robin.gui.theme import clamp_qtable_server_pagination


def _write_ann_vcf(path: Path) -> None:
    path.write_text(
        (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            "chr1\t1000\t.\tA\tG\t80\tPASS\t"
            "ANN=G|missense_variant|MODERATE|TP53|ENSG0001|transcript|ENST1|"
            "protein_coding|5/11|c.100A>G|p.Thr34Ala|100/1200|100/1100|34/366||;"
            "CLNSIG=Pathogenic;DP=40\tGT\t0/1\n"
            "chr1\t2000\t.\tC\tT\t12\tLowQual\t"
            "ANN=T|synonymous_variant|LOW|EGFR|ENSG0002|transcript|ENST2|"
            "protein_coding|2/8|c.20C>T|p.Ser7Ser|20/900|20/800|7/266||;"
            "DP=8\tGT\t0/1\n"
            "chr2\t500\t.\tG\tA\t55\tPASS\t"
            "ANN=A|missense_variant|MODERATE|BRAF|ENSG0003|transcript|ENST3|"
            "protein_coding|15/18|c.1799T>A|p.Val600Glu|1799/2000|1799/1900|600/766||;"
            "ONC=Oncogenic;DP=22\tGT\t0/1\n"
        ),
        encoding="utf-8",
    )


def test_write_bundle_omits_rows_from_sidecar(tmp_path: Path) -> None:
    vcf = tmp_path / "snpsift_output.vcf"
    _write_ann_vcf(vcf)
    display = build_snp_display_data(vcf)
    assert display is not None
    json_path = tmp_path / SNP_DISPLAY_JSON
    parquet_path = tmp_path / SNP_DISPLAY_PARQUET
    write_variant_display_bundle(display, json_path, parquet_path)

    sidecar = json.loads(json_path.read_text(encoding="utf-8"))
    assert sidecar["format"] == VARIANT_SIDECAR_FORMAT
    assert "rows_all" not in sidecar
    assert "rows_pathogenic" not in sidecar
    assert sidecar["summary"]["total_variants"] == 3
    assert parquet_path.is_file()
    assert sidecar["parquet"] == SNP_DISPLAY_PARQUET


def test_parquet_store_pages_and_filters(tmp_path: Path) -> None:
    vcf = tmp_path / "snpsift_output.vcf"
    _write_ann_vcf(vcf)
    write_clair_variant_display_files(tmp_path)

    store = VariantDisplayStore.open_snp(tmp_path)
    assert store is not None
    assert store.uses_parquet
    assert store.total_variants == 3

    page, total = store.page(VariantTableFilters(), offset=0, limit=2)
    assert total == 3
    assert len(page) == 2
    assert page[0]["__row_idx"] == 0

    pass_page, pass_total = store.page(VariantTableFilters(pass_only=True), offset=0, limit=10)
    assert pass_total == 2
    assert {row["FILTER"] for row in pass_page} == {"PASS"}

    sig_page, sig_total = store.page(
        VariantTableFilters(significant_only=True), offset=0, limit=10
    )
    assert sig_total >= 1
    assert all(
        str(row.get("is_clinvar_significant", "")).upper() in {"YES", "TRUE", "1", "PATHOGENIC"}
        for row in sig_page
    )

    gene_page, gene_total = store.page(
        VariantTableFilters(search_text="EGFR", search_fields=("Gene_Name", "CHROM", "POS")),
        offset=0,
        limit=10,
    )
    assert gene_total == 1
    assert "EGFR" in str(gene_page[0].get("Gene_Name", ""))

    low_qual, low_total = store.page(VariantTableFilters(min_qual=50), offset=0, limit=10)
    assert low_total == 2
    assert all(float(row["QUAL"]) >= 50 for row in low_qual)

    fetched = store.fetch_row(int(page[1]["__row_idx"]))
    assert fetched is not None
    assert fetched["POS"] == page[1]["POS"]


def test_open_falls_back_to_legacy_json_rows(tmp_path: Path) -> None:
    vcf = tmp_path / "snpsift_output.vcf"
    _write_ann_vcf(vcf)
    display = build_snp_display_data(vcf)
    assert display is not None
    json_path = tmp_path / SNP_DISPLAY_JSON
    json_path.write_text(json.dumps(display), encoding="utf-8")

    store = VariantDisplayStore.open_snp(tmp_path)
    assert store is not None
    assert not store.uses_parquet
    page, total = store.page(VariantTableFilters(pass_only=True), offset=0, limit=10)
    assert total == 2
    assert len(page) == 2


def test_indel_bundle_round_trip(tmp_path: Path) -> None:
    vcf = tmp_path / "snpsift_indel_output.vcf"
    _write_ann_vcf(vcf)
    write_clair_variant_display_files(tmp_path)
    assert (tmp_path / INDEL_DISPLAY_JSON).is_file()
    assert (tmp_path / INDEL_DISPLAY_PARQUET).is_file()
    store = VariantDisplayStore.open_indel(tmp_path, vcf_fallback=False)
    assert store is not None
    assert store.uses_parquet
    _rows, total = store.page(VariantTableFilters(), offset=0, limit=10)
    assert total == 3


def test_clamp_caps_quasar_all(tmp_path: Path) -> None:
    pag = clamp_qtable_server_pagination(
        {"page": 1, "rowsPerPage": 0, "rowsNumber": 50_000},
        rows_number=50_000,
        rows_per_page_default=100,
        rows_per_page_max=250,
    )
    assert pag["rowsPerPage"] == 250
    assert pag["page"] == 1


def test_unfiltered_page_uses_sidecar_total(tmp_path: Path) -> None:
    vcf = tmp_path / "snpsift_output.vcf"
    _write_ann_vcf(vcf)
    write_clair_variant_display_files(tmp_path)
    store = VariantDisplayStore.open_snp(tmp_path)
    assert store is not None
    filters = VariantTableFilters()
    assert filters.is_active is False
    page, total = store.page(filters, offset=0, limit=2)
    assert total == 3
    assert len(page) == 2
    empty, total_only = store.page(filters, offset=0, limit=0)
    assert empty == []
    assert total_only == 3
    assert VariantTableFilters(pass_only=True).is_active is True
