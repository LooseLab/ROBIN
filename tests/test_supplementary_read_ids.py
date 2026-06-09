from pathlib import Path
from unittest.mock import Mock, patch

from robin.analysis.bam_preprocessor import (
    BamMetadata,
    _persist_supplementary_read_ids,
    calculate_bam_summary,
)
from robin.analysis.fusion_analysis import _load_supplementary_read_ids
from robin.analysis.fusion_work import FusionMetadata, process_bam_with_staging


def _metadata(read_ids):
    return BamMetadata(
        file_path="input.bam",
        file_size=1,
        creation_time=0.0,
        extracted_data={
            "sample_id": "S1",
            "supplementary_read_ids": list(read_ids),
        },
    )


def test_supplementary_ids_are_persisted_per_bam(tmp_path: Path):
    first = _metadata(["read-b", "read-a"])
    second = _metadata(["read-c"])

    _persist_supplementary_read_ids(
        first,
        str(tmp_path / "run-1" / "chunk.bam"),
        str(tmp_path),
    )
    _persist_supplementary_read_ids(
        second,
        str(tmp_path / "run-2" / "chunk.bam"),
        str(tmp_path),
    )

    first_path = Path(first.extracted_data["supplementary_read_ids_path"])
    second_path = Path(second.extracted_data["supplementary_read_ids_path"])
    assert first_path != second_path
    assert first_path.parent.name == "_supplementary_read_ids"
    assert first_path.read_text().splitlines() == ["read-a", "read-b"]
    assert second_path.read_text().splitlines() == ["read-c"]
    assert first.extracted_data["supplementary_read_ids_complete"] is True
    assert first.extracted_data["supplementary_read_ids_count"] == 2
    assert "supplementary_read_ids" not in first.extracted_data


def test_valid_per_bam_ids_remain_marked_complete(tmp_path: Path):
    metadata = _metadata(["read-a", "read-b"])
    _persist_supplementary_read_ids(
        metadata,
        str(tmp_path / "chunk.bam"),
        str(tmp_path),
    )

    loaded = _load_supplementary_read_ids(metadata.extracted_data, Mock())

    assert loaded == ["read-a", "read-b"]
    assert metadata.extracted_data["supplementary_read_ids_complete"] is True


def test_missing_or_incomplete_id_file_falls_back_to_sa_scanning(tmp_path: Path):
    missing = {
        "supplementary_read_ids_path": str(tmp_path / "missing.txt"),
        "supplementary_read_ids_count": 2,
        "supplementary_read_ids_complete": True,
    }
    assert _load_supplementary_read_ids(missing, Mock()) == []
    assert missing["supplementary_read_ids_complete"] is False

    mismatched_path = tmp_path / "mismatched.txt"
    mismatched_path.write_text("read-a\n")
    mismatched = {
        "supplementary_read_ids_path": str(mismatched_path),
        "supplementary_read_ids_count": 2,
        "supplementary_read_ids_complete": True,
    }
    assert _load_supplementary_read_ids(mismatched, Mock()) == ["read-a"]
    assert mismatched["supplementary_read_ids_complete"] is False


def test_bam_summary_preserves_supplementary_id_completeness():
    summary = calculate_bam_summary(
        {
            "has_supplementary_reads": True,
            "supplementary_read_ids": ["read-a"],
            "supplementary_read_ids_complete": True,
        }
    )

    assert summary["supplementary_read_ids_complete"] is True


def test_complete_flag_reaches_single_pass_scanner(tmp_path: Path):
    fusion_metadata = FusionMetadata(
        sample_id="S1",
        file_path="input.bam",
        analysis_timestamp=0.0,
        target_panel="rCNS2",
    )
    with patch(
        "robin.analysis.fusion_work.process_bam_single_pass",
        return_value=(None, None, None),
    ) as scan:
        process_bam_with_staging(
            "input.bam",
            str(tmp_path),
            {"supplementary_read_ids_complete": True},
            fusion_metadata,
            "rCNS2",
            has_supplementary=True,
            supplementary_read_ids=["read-a"],
            work_dir=str(tmp_path),
            batch_size=10,
        )

    assert scan.call_args.args[4] == ["read-a"]
    assert scan.call_args.args[5] is True
