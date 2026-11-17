"""Tests for incremental target.bam folding during accumulation/finalize."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pysam
import pytest

from robin.analysis.target_analysis import (
    _bam_has_any_alignment,
    finalize_accumulation_for_sample,
    fold_bam_into_target,
    sample_needs_target_bam_finalize,
)


def _write_tiny_bam(path: Path, *, start: int, name: str, chrom: str = "chr1") -> None:
    header = {
        "HD": {"VN": "1.0", "SO": "coordinate"},
        "SQ": [{"LN": 1000, "SN": chrom}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        aln = pysam.AlignedSegment()
        aln.query_name = name
        aln.query_sequence = "ACGT"
        aln.flag = 0
        aln.reference_id = 0
        aln.reference_start = start
        aln.mapping_quality = 20
        aln.cigar = ((0, 4),)
        aln.query_qualities = pysam.qualitystring_to_array("IIII")
        out.write(aln)
    pysam.index(str(path))


def _write_empty_bam(path: Path, *, chrom: str = "chr1") -> None:
    header = {
        "HD": {"VN": "1.0", "SO": "coordinate"},
        "SQ": [{"LN": 1000, "SN": chrom}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header):
        pass


def _read_names(bam_path: Path) -> list[str]:
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        return [aln.query_name for aln in bam.fetch(until_eof=True)]


def test_fold_first_batch_copies_and_indexes(tmp_path: Path) -> None:
    incoming = tmp_path / "batch_1.bam"
    target = tmp_path / "target.bam"
    _write_tiny_bam(incoming, start=10, name="read-a")

    assert fold_bam_into_target(str(target), str(incoming)) == "copied"
    assert target.is_file()
    assert Path(f"{target}.bai").is_file()
    assert _read_names(target) == ["read-a"]


def test_fold_second_batch_is_two_way_merge(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    merge_input_counts: list[int] = []
    real_merge = pysam.merge

    def _spy_merge(*args):
        args_list = list(args)
        out_idx = args_list.index("-o")
        merge_input_counts.append(len(args_list[out_idx + 2 :]))
        return real_merge(*args)

    monkeypatch.setattr("robin.analysis.target_analysis.pysam.merge", _spy_merge)

    target = tmp_path / "target.bam"
    first = tmp_path / "batch_1.bam"
    second = tmp_path / "batch_2.bam"
    _write_tiny_bam(first, start=10, name="read-a")
    _write_tiny_bam(second, start=50, name="read-b")

    assert fold_bam_into_target(str(target), str(first)) == "copied"
    assert fold_bam_into_target(str(target), str(second)) == "merged"
    assert merge_input_counts == [2]
    assert set(_read_names(target)) == {"read-a", "read-b"}
    assert Path(f"{target}.bai").is_file()
    assert not (tmp_path / ".target_fold.tmp.bam").exists()


def test_fold_skips_empty_incoming(tmp_path: Path) -> None:
    target = tmp_path / "target.bam"
    empty = tmp_path / "empty.bam"
    _write_empty_bam(empty)

    assert not _bam_has_any_alignment(str(empty))
    assert fold_bam_into_target(str(target), str(empty)) == "skipped_empty"
    assert not target.exists()


def test_sample_needs_finalize_for_batches_or_staging(tmp_path: Path) -> None:
    sample_dir = tmp_path / "S1"
    sample_dir.mkdir()
    assert sample_needs_target_bam_finalize(str(sample_dir)) is False

    batch = sample_dir / "batch_1.bam"
    batch.write_bytes(b"BAM")
    assert sample_needs_target_bam_finalize(str(sample_dir)) is True

    batch.unlink()
    staging = sample_dir / "_staging"
    staging.mkdir()
    (staging / "coverage_000001.parquet").write_bytes(b"x")
    assert sample_needs_target_bam_finalize(str(sample_dir)) is True


def test_finalize_folds_leftover_batches_sequentially(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sample_id = "S1"
    sample_dir = tmp_path / sample_id
    sample_dir.mkdir()
    first = sample_dir / "batch_1.bam"
    second = sample_dir / "batch_2.bam"
    _write_tiny_bam(first, start=10, name="read-a")
    _write_tiny_bam(second, start=50, name="read-b")

    merge_input_counts: list[int] = []
    real_merge = pysam.merge

    def _spy_merge(*args):
        args_list = list(args)
        out_idx = args_list.index("-o")
        merge_input_counts.append(len(args_list[out_idx + 2 :]))
        return real_merge(*args)

    monkeypatch.setattr("robin.analysis.target_analysis.pysam.merge", _spy_merge)
    monkeypatch.setattr(
        "robin.analysis.target_analysis.TargetAnalysis._get_pending_count",
        lambda self, sid: 0,
    )

    result = finalize_accumulation_for_sample(
        sample_id=sample_id,
        work_dir=str(tmp_path),
        target_panel="rCNS2",
    )

    assert result["final_merge"] == "success"
    assert result["batch_files_merged"] == 2
    assert merge_input_counts == [2]
    assert not first.exists()
    assert not second.exists()
    target = sample_dir / "target.bam"
    assert target.is_file()
    assert set(_read_names(target)) == {"read-a", "read-b"}


def test_finalize_already_folded_when_target_exists(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sample_id = "S1"
    sample_dir = tmp_path / sample_id
    sample_dir.mkdir()
    target = sample_dir / "target.bam"
    _write_tiny_bam(target, start=10, name="read-a")

    monkeypatch.setattr(
        "robin.analysis.target_analysis.TargetAnalysis._get_pending_count",
        lambda self, sid: 0,
    )

    result = finalize_accumulation_for_sample(
        sample_id=sample_id,
        work_dir=str(tmp_path),
        target_panel="rCNS2",
    )

    assert result["final_merge"] == "already_folded"
    assert result["batch_files_merged"] == 0
    assert _read_names(target) == ["read-a"]


def test_failed_fold_leaves_batch_for_finalize(tmp_path: Path) -> None:
    incoming = tmp_path / "batch_1.bam"
    target = tmp_path / "target.bam"
    header = {
        "HD": {"VN": "1.0", "SO": "coordinate"},
        "SQ": [{"LN": 1000, "SN": "chr1"}],
    }
    with pysam.AlignmentFile(str(incoming), "wb", header=header) as out:
        aln = pysam.AlignedSegment()
        aln.query_name = "read-a"
        aln.query_sequence = "ACGT"
        aln.flag = 0
        aln.reference_id = 0
        aln.reference_start = 10
        aln.mapping_quality = 20
        aln.cigar = ((0, 4),)
        aln.query_qualities = pysam.qualitystring_to_array("IIII")
        out.write(aln)

    with patch(
        "robin.analysis.target_analysis.pysam.index",
        side_effect=RuntimeError("index boom"),
    ):
        with pytest.raises(RuntimeError, match="index boom"):
            fold_bam_into_target(str(target), str(incoming))

    assert incoming.is_file()
