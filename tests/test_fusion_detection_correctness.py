from pathlib import Path
from unittest.mock import Mock, patch

import pandas as pd
import pysam

from robin.analysis.fusion_analysis import process_multiple_files
from robin.analysis.fusion_work import (
    FusionMetadata,
    GeneRegion,
    _generate_output_files,
    process_bam_single_pass,
)


def _write_split_read_bam(
    path: Path,
    include_sa_tag: bool,
    include_supplementary_record: bool = False,
) -> None:
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [
            {"SN": "chr1", "LN": 10_000},
            {"SN": "chr2", "LN": 10_000},
        ],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        primary = pysam.AlignedSegment(bam.header)
        primary.query_name = "read-1"
        primary.query_sequence = "A" * 600
        primary.flag = 0
        primary.reference_id = 0
        primary.reference_start = 100
        primary.mapping_quality = 60
        primary.cigarstring = "300M300S"
        if include_sa_tag:
            primary.set_tag("SA", "chr2,501,+,300S300M,60,0;")
        bam.write(primary)

        if include_supplementary_record or not include_sa_tag:
            supplementary = pysam.AlignedSegment(bam.header)
            supplementary.query_name = "read-1"
            supplementary.query_sequence = "A" * 600
            supplementary.flag = 2048
            supplementary.reference_id = 1
            supplementary.reference_start = 500
            supplementary.mapping_quality = 60
            supplementary.cigarstring = "300S300M"
            bam.write(supplementary)


def _scan_with_test_regions(bam_path: Path):
    target_regions = {
        "chr1": [GeneRegion(50, 450, "GENE1")],
        "chr2": [GeneRegion(450, 850, "GENE2")],
    }
    starts = {"chr1": [50], "chr2": [450]}
    with (
        patch("robin.analysis.fusion_work._ensure_gene_regions_loaded"),
        patch.dict(
            "robin.analysis.fusion_work._gene_regions_cache",
            {"test": target_regions},
            clear=True,
        ),
        patch.dict(
            "robin.analysis.fusion_work._all_gene_regions_cache",
            {"shared": target_regions},
            clear=True,
        ),
        patch.dict(
            "robin.analysis.fusion_work._gene_region_starts_cache",
            {"test": starts},
            clear=True,
        ),
        patch.dict(
            "robin.analysis.fusion_work._all_gene_region_starts_cache",
            {"shared": starts},
            clear=True,
        ),
        patch.dict(
            "robin.analysis.fusion_work._gene_region_ncls_cache",
            {"test": {}},
            clear=True,
        ),
        patch.dict(
            "robin.analysis.fusion_work._all_gene_region_ncls_cache",
            {"shared": {}},
            clear=True,
        ),
        patch.dict(
            "robin.analysis.fusion_work._combined_gene_region_ncls_cache",
            {"test": {}},
            clear=True,
        ),
    ):
        return process_bam_single_pass(
            str(bam_path),
            "test",
            supplementary_read_ids=["read-1"],
            supplementary_read_ids_complete=True,
        )


def test_sa_only_alignment_is_gene_annotated(tmp_path: Path):
    bam_path = tmp_path / "sa-only.bam"
    _write_split_read_bam(bam_path, include_sa_tag=True)

    target, genome, _ = _scan_with_test_regions(bam_path)

    assert set(target["col4"]) == {"GENE1", "GENE2"}
    assert set(genome["col4"]) == {"GENE1", "GENE2"}


def test_supplementary_record_without_sa_tag_is_gene_annotated(tmp_path: Path):
    bam_path = tmp_path / "supplementary-record.bam"
    _write_split_read_bam(bam_path, include_sa_tag=False)

    target, genome, _ = _scan_with_test_regions(bam_path)

    assert set(target["col4"]) == {"GENE1", "GENE2"}
    assert set(genome["col4"]) == {"GENE1", "GENE2"}


def test_sa_and_supplementary_record_are_one_canonical_alignment(tmp_path: Path):
    bam_path = tmp_path / "sa-and-supplementary-record.bam"
    _write_split_read_bam(
        bam_path,
        include_sa_tag=True,
        include_supplementary_record=True,
    )

    target, genome, _ = _scan_with_test_regions(bam_path)

    assert list(target["col4"]).count("GENE1") == 1
    assert list(target["col4"]).count("GENE2") == 1
    assert list(genome["col4"]).count("GENE1") == 1
    assert list(genome["col4"]).count("GENE2") == 1


def _candidate_rows(read_id: str, gene: str, chromosome: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "col1": [chromosome],
            "col2": [100],
            "col3": [500],
            "col4": [gene],
            "reference_id": [chromosome],
            "reference_start": [150],
            "reference_end": [450],
            "read_id": [read_id],
            "mapping_quality": [60],
            "strand": ["+"],
            "read_start": [0],
            "read_end": [300],
            "is_secondary": [False],
            "is_supplementary": [False],
            "mapping_span": [300],
        }
    )


def test_sample_wide_calling_combines_evidence_across_parquet_parts(tmp_path: Path):
    sample_id = "S1"
    dataset = tmp_path / sample_id / "target_candidates_dataset"
    dataset.mkdir(parents=True)

    parts = [
        _candidate_rows("read-1", "GENE1", "chr1"),
        pd.concat(
            [
                _candidate_rows("read-1", "GENE2", "chr2"),
                _candidate_rows("read-2", "GENE1", "chr1"),
                _candidate_rows("read-2", "GENE2", "chr2"),
                _candidate_rows("read-3", "GENE1", "chr1"),
                _candidate_rows("read-3", "GENE2", "chr2"),
            ],
            ignore_index=True,
        ),
    ]
    for index, part in enumerate(parts):
        part.to_parquet(dataset / f"part_{index}.parquet", index=False)

    metadata = FusionMetadata(
        sample_id=sample_id,
        file_path="accumulated",
        analysis_timestamp=0,
        target_panel="test",
    )
    with patch(
        "robin.analysis.fusion_work._generate_fusion_breakpoint_bed",
        return_value=False,
    ):
        _generate_output_files(
            sample_id,
            {},
            metadata,
            str(tmp_path),
            generate_master_bed=False,
        )

    result = pd.read_csv(tmp_path / sample_id / "fusion_candidates_master.csv")
    assert set(result["read_id"]) == {"read-1", "read-2", "read-3"}
    assert set(result["tag"]) == {"GENE1,GENE2"}


def test_multi_file_result_reports_real_counts_and_paths(tmp_path: Path):
    metadata = [{"sample_id": "S1", "has_supplementary_reads": True}]
    with (
        patch("robin.analysis.fusion_work._ensure_gene_regions_loaded"),
        patch(
            "robin.analysis.fusion_analysis.process_bam_with_staging",
            return_value=({}, False),
        ),
        patch(
            "robin.analysis.fusion_analysis.accumulate_fusion_candidates",
            return_value={
                "status": "success",
                "target_candidates": 7,
                "genome_wide_candidates": 11,
            },
        ),
    ):
        result = process_multiple_files(
            ["input.bam"],
            metadata,
            str(tmp_path),
            Mock(),
            target_panel="test",
        )

    assert result["fusion_data"]["target_candidates_count"] == 7
    assert result["fusion_data"]["genome_wide_candidates_count"] == 11
    assert result["target_fusion_path"].endswith("fusion_candidates_master.csv")
    assert result["genome_wide_fusion_path"].endswith("fusion_candidates_all.csv")
