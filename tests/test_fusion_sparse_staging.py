from pathlib import Path
from unittest.mock import patch

import pandas as pd

from robin.analysis.fusion_work import (
    FusionMetadata,
    _get_pending_count,
    accumulate_fusion_candidates,
    process_bam_with_staging,
)


def _metadata(sample_id: str) -> FusionMetadata:
    return FusionMetadata(
        sample_id=sample_id,
        file_path="input.bam",
        analysis_timestamp=0.0,
        target_panel="rCNS2",
    )


def test_sparse_staging_writes_only_non_empty_candidate_types(tmp_path: Path):
    sample_id = "S1"
    genome_candidates = pd.DataFrame(
        {
            "read_id": ["read-1", "read-1"],
            "col4": ["GENE1", "GENE2"],
            "reference_id": ["chr1", "chr2"],
            "reference_start": [100, 300],
            "reference_end": [200, 400],
        }
    )

    with patch(
        "robin.analysis.fusion_work.process_bam_single_pass",
        return_value=(None, genome_candidates, None),
    ):
        _, should_accumulate = process_bam_with_staging(
            "input.bam",
            str(tmp_path),
            {},
            _metadata(sample_id),
            "rCNS2",
            has_supplementary=True,
            work_dir=str(tmp_path),
            batch_size=1,
        )

    staging_dir = tmp_path / sample_id / "_fusion_staging"
    assert should_accumulate
    assert list(staging_dir.glob("target_*.parquet")) == []
    assert len(list(staging_dir.glob("genome_*.parquet"))) == 1
    assert list(staging_dir.glob("master_bed_*.parquet")) == []
    assert _get_pending_count(str(tmp_path), sample_id) == 1

    with patch("robin.analysis.fusion_work._generate_output_files"):
        result = accumulate_fusion_candidates(
            str(tmp_path),
            sample_id,
            "rCNS2",
            force=True,
            batch_size=1,
        )

    assert result["status"] == "success"
    assert result["files_processed"] == 1
    assert result["genome_wide_candidates"] == 2
    assert _get_pending_count(str(tmp_path), sample_id) == 0
    assert (
        len(
            list(
                (tmp_path / sample_id / "genome_wide_candidates_dataset").glob(
                    "*.parquet"
                )
            )
        )
        == 1
    )


def test_fully_empty_candidate_result_needs_no_placeholder_parquet(tmp_path: Path):
    sample_id = "S1"
    with patch(
        "robin.analysis.fusion_work.process_bam_single_pass",
        return_value=(None, None, None),
    ):
        _, should_accumulate = process_bam_with_staging(
            "input.bam",
            str(tmp_path),
            {},
            _metadata(sample_id),
            "rCNS2",
            has_supplementary=True,
            work_dir=str(tmp_path),
            batch_size=1,
        )

    staging_dir = tmp_path / sample_id / "_fusion_staging"
    assert should_accumulate
    assert list(staging_dir.glob("*.parquet")) == []
    assert _get_pending_count(str(tmp_path), sample_id) == 1

    with patch("robin.analysis.fusion_work._generate_output_files") as generate_outputs:
        result = accumulate_fusion_candidates(
            str(tmp_path),
            sample_id,
            "rCNS2",
            force=True,
            batch_size=1,
        )

    assert result["status"] == "success"
    assert result["files_processed"] == 1
    assert result["target_candidates"] == 0
    assert result["genome_wide_candidates"] == 0
    assert result["master_bed_candidates"] == 0
    assert _get_pending_count(str(tmp_path), sample_id) == 0
    generate_outputs.assert_not_called()
