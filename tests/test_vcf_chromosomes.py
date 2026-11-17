"""Tests for VCF chromosome normalization before SnpSift."""

from pathlib import Path

from robin.analysis.utilities.vcf_chromosomes import (
    chrom_to_clinvar,
    chrom_to_ucsc,
    count_vcf_variants,
    normalize_vcf_chromosomes_for_clinvar,
    restore_vcf_chromosomes,
)


def test_chrom_to_clinvar_and_ucsc_roundtrip() -> None:
    assert chrom_to_clinvar("chr2") == "2"
    assert chrom_to_clinvar("chrX") == "X"
    assert chrom_to_clinvar("chrM") == "MT"
    assert chrom_to_clinvar("2") == "2"
    assert chrom_to_ucsc("2") == "chr2"
    assert chrom_to_ucsc("MT") == "chrM"


def test_normalize_and_restore_vcf(tmp_path: Path) -> None:
    src = tmp_path / "in.vcf"
    normalized = tmp_path / "norm.vcf"
    raw = tmp_path / "raw.vcf"
    restored = tmp_path / "out.vcf"

    src.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr2,length=242193529>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "chr2\t208248388\t.\tC\tT\t68.8\tPASS\tANN=missense",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    mapping = normalize_vcf_chromosomes_for_clinvar(src, normalized)
    assert mapping == {"2": "chr2"}
    assert count_vcf_variants(normalized) == 1
    norm_text = normalized.read_text(encoding="utf-8")
    assert "\n2\t208248388\t" in norm_text
    assert "##contig=<ID=2," in norm_text

    raw.write_text(
        norm_text.replace("PASS", "PASS;CLNSIG=Pathogenic;ONC=Oncogenic"),
        encoding="utf-8",
    )
    restore_vcf_chromosomes(raw, restored, mapping)
    out_text = restored.read_text(encoding="utf-8")
    assert "\nchr2\t208248388\t" in out_text
    assert "CLNSIG=Pathogenic" in out_text
