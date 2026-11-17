"""Issue #150: the 50,000-read BAM limit must count reads, not alignment records.

MinKNOW rolls BAM files over by read count. One read can carry several
alignment records (its primary plus supplementary pieces, or secondary
placements), so the limit counts only records that are neither secondary nor
supplementary: exactly one per read.

These tests build their own small BAMs, so they run without the binary
fixtures that tests/test_bam_preprocessor.py skips on.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import pysam
import pytest

import robin.analysis.bam_preprocessor as bam_preprocessor
from robin.analysis.bam_preprocessor import (
    bam_preprocessing_handler,
    extract_bam_metadata,
    process_bam_reads,
)
from robin.logging_config import JobLogger

_READ_LIMIT = 50_000

_PRIMARY = 0x0
_UNMAPPED = 0x4
_SECONDARY = 0x100
_SUPPLEMENTARY = 0x800

_HEADER = {
    "HD": {"VN": "1.6", "SO": "unsorted"},
    "SQ": [{"SN": "chr1", "LN": 10_000_000}, {"SN": "chr2", "LN": 10_000_000}],
    "RG": [
        {
            "ID": "rg1",
            "LB": "Sample150",
            # A CpG modbase model keeps the handler from raising a GUI warning.
            "DS": "runid=run150 basecall_model=model modbase_models=model_5mCG_5hmCG@v1",
        }
    ],
}

# (read name, flag, contig index, start)
Record = Tuple[str, int, int, int]

# Three reads, seven records.
_MIXED_RECORDS: List[Record] = [
    ("readA", _PRIMARY, 0, 100),
    ("readA", _SUPPLEMENTARY, 1, 5_000),
    ("readA", _SUPPLEMENTARY, 0, 9_000),
    ("readB", _PRIMARY, 0, 200),
    ("readB", _SECONDARY, 0, 800),
    ("readB", _SECONDARY | _SUPPLEMENTARY, 0, 7_000),
    ("readC", _UNMAPPED, -1, -1),
]


def _write_bam(path: Path, records: Iterable[Record]) -> str:
    seq = "ACGT" * 5
    quals = pysam.qualitystring_to_array("I" * len(seq))
    with pysam.AlignmentFile(str(path), "wb", header=_HEADER) as bam:
        for name, flag, contig, start in records:
            seg = pysam.AlignedSegment()
            seg.query_name = name
            seg.flag = flag
            seg.query_sequence = seq
            seg.query_qualities = quals
            if flag & _UNMAPPED:
                seg.reference_id = -1
                seg.reference_start = -1
            else:
                seg.reference_id = contig
                seg.reference_start = start
                seg.cigarstring = f"{len(seq)}M"
                seg.mapping_quality = 60
            bam.write(seg)
    return str(path)


def _primaries(count: int) -> List[Record]:
    return [(f"read{i}", _PRIMARY, 0, (i * 100) % 9_000_000) for i in range(count)]


def _supplementaries(count: int) -> List[Record]:
    """Supplementary pieces of the first ``count`` reads, placed on chr2."""
    return [(f"read{i}", _SUPPLEMENTARY, 1, i * 100) for i in range(count)]


def test_process_bam_reads_counts_each_read_once(tmp_path: Path) -> None:
    out = process_bam_reads(_write_bam(tmp_path / "mixed.bam", _MIXED_RECORDS))

    assert out is not None
    # readA, readB and readC each have one primary or unmapped record.
    assert out["primary_read_count"] == 3
    # Alignment-record statistics are unchanged: supplementary records still count.
    assert out["mapped_reads"] == 4
    assert out["unmapped_reads"] == 1
    assert out["supplementary_reads"] == 3


def test_primary_read_count_reaches_extracted_metadata(tmp_path: Path) -> None:
    """The handler reads extracted_data, which keeps only keys that
    calculate_bam_summary passes through."""
    meta = extract_bam_metadata(_write_bam(tmp_path / "mixed.bam", _MIXED_RECORDS))

    assert meta.extracted_data["primary_read_count"] == 3


class _Context:
    """The parts of WorkflowContext that bam_preprocessing_handler uses."""

    def __init__(self, filepath: str, work_dir: Path) -> None:
        self.filepath = filepath
        self.metadata: Dict[str, Any] = {"work_dir": str(work_dir)}
        self.results: Dict[str, Any] = {}
        self.errors: List[Dict[str, Any]] = []

    def add_metadata(self, key: str, value: Any) -> None:
        self.metadata[key] = value

    def add_result(self, job_type: str, result: Any) -> None:
        self.results[job_type] = result

    def add_error(self, job_type: str, error: str) -> None:
        self.errors.append({"job_type": job_type, "error": error})


class _Job:
    def __init__(self, filepath: str, work_dir: Path) -> None:
        self.job_id = 150
        self.job_type = "preprocessing"
        self.context = _Context(filepath, work_dir)


@pytest.fixture
def run_handler(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    # get_job_logger's first call writes robin.log into the working directory
    # and replaces the root logger's handlers; a plain JobLogger avoids both.
    monkeypatch.setattr(
        bam_preprocessor,
        "get_job_logger",
        lambda job_id, job_type, filepath=None: JobLogger(
            logging.getLogger("robin.preprocessing"), job_id, job_type, filepath
        ),
    )
    # The limit applies only when large-BAM mode is off, whatever the shell sets.
    monkeypatch.setattr(bam_preprocessor, "_PROCESS_LARGE_BAMS_INDIVIDUALLY", False)

    def _run(records: Iterable[Record]) -> _Job:
        job = _Job(_write_bam(tmp_path / "chunk.bam", records), tmp_path / "work")
        bam_preprocessing_handler(job)
        # The handler turns exceptions into context errors; surface them here.
        assert job.context.errors == []
        return job

    return _run


def test_handler_admits_bam_whose_extra_records_are_supplementary(run_handler) -> None:
    # The reporter's MinKNOW setting was 49,998 reads per file. Five chimeric
    # reads add supplementary records, making 50,003 alignment records.
    job = run_handler(_primaries(49_998) + _supplementaries(5))

    metadata = job.context.metadata
    assert metadata.get("skip_downstream") is None, (
        f"skipped: {metadata.get('skip_reason')}, {metadata.get('skip_read_count')} reads"
    )
    assert job.context.results["preprocessing"]["status"] == "success"


def test_handler_still_skips_bam_with_more_than_50000_reads(run_handler) -> None:
    job = run_handler(_primaries(_READ_LIMIT + 1))

    assert job.context.metadata["skip_downstream"] is True
    assert job.context.metadata["skip_reason"] == "too_many_reads"
    assert job.context.metadata["skip_read_count"] == _READ_LIMIT + 1
    assert job.context.results["preprocessing"]["status"] == "skipped"
