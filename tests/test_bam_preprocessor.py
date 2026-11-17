"""Tests for `robin.analysis.bam_preprocessor` using committed tiny BAM fixtures."""

from __future__ import annotations

import os
from pathlib import Path

import pysam
import pytest

from robin.analysis.bam_preprocessor import (
    BamMetadata,
    calculate_bam_summary,
    extract_bam_metadata,
    get_rg_tags_from_bam,
    process_bam_reads,
    _extract_sample_id_from_bam,
    _get_modbase_model_warning,
    _get_modbase_model_warning_level,
)

_FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "bam"
_PASS_BAM = _FIXTURE_DIR / "test_pass.bam"
_FAIL_BAM = _FIXTURE_DIR / "test_fail.bam"


pytestmark = pytest.mark.skipif(
    not _PASS_BAM.is_file() or not _FAIL_BAM.is_file(),
    reason="BAM fixtures missing: tests/fixtures/bam/test_pass.bam and test_fail.bam",
)


def _primary_mapped_unmapped(bam_path: str) -> tuple[int, int]:
    """
    Mirror `process_bam_reads`: only non-secondary alignments contribute to
    mapped/unmapped totals.
    """
    mapped = 0
    unmapped = 0
    with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary:
                continue
            if read.is_unmapped:
                unmapped += 1
            else:
                mapped += 1
    return mapped, unmapped


def test_fixture_paths_use_pass_vs_fail_substring() -> None:
    """`state` comes from `'pass' in bam_file` (full path string); fail fixture must not contain ``pass``."""
    assert "pass" in str(_PASS_BAM).lower()
    assert "pass" not in str(_FAIL_BAM).lower()


def test_get_rg_tags_from_bam_matches_header() -> None:
    with pysam.AlignmentFile(str(_PASS_BAM), "rb", check_sq=False) as bam:
        rg = get_rg_tags_from_bam(bam)
        assert rg is not None
        hdr0 = bam.header.get("RG", [])[0]
        assert rg[4] == hdr0.get("LB")
        assert rg[3] == hdr0.get("DS", "").split()[0].removeprefix("runid=")


def test_get_rg_tags_from_bam_none_when_no_rg() -> None:
    mock = type("M", (), {})()
    mock.header = {"RG": []}
    assert get_rg_tags_from_bam(mock) is None
    assert get_rg_tags_from_bam(None) is None


def test_get_rg_tags_from_bam_extracts_modbase_models_by_key() -> None:
    mock = type("M", (), {})()
    mock.header = {
        "RG": [
            {
                "ID": "run_model",
                "DS": (
                    "runid=run basecall_model=model "
                    "experiment_id=sample "
                    "modbase_models=model_5mCG_5hmCG@v1"
                ),
            }
        ]
    }

    rg = get_rg_tags_from_bam(mock)

    assert rg is not None
    assert rg[2] == "model"
    assert rg[3] == "run"
    assert rg[9] == "model_5mCG_5hmCG@v1"


@pytest.mark.parametrize(
    ("model", "warning_fragment", "level"),
    [
        ("model_5mCG_5hmCG@v1", None, None),
        (
            "dna_r10.4.1_e8.2_400bps_hac@v6.0.0_5mC_5hmC@v1",
            "incorrect and slower than expected",
            "warning",
        ),
        ("model_6mA@v1", "without 5mCG_5hmCG", "warning"),
        (
            "modbase_model_version_id",
            "does not record which modbase model was used",
            "info",
        ),
        (None, "does not report modbase_models", "warning"),
    ],
)
def test_modbase_model_warning(model, warning_fragment, level) -> None:
    warning = _get_modbase_model_warning(model)
    if warning_fragment is None:
        assert warning is None
        assert level is None
    else:
        assert warning_fragment in warning
        assert _get_modbase_model_warning_level(model) == level


def test_process_bam_reads_pass_fixture() -> None:
    out = process_bam_reads(str(_PASS_BAM))
    assert out is not None
    assert out["state"] == "pass"
    assert out["sample_id"] == "Sample_104"
    assert out["mapped_reads"] == 1
    assert out["unmapped_reads"] == 0
    assert out["pass_mapped_reads_num"] == 1
    assert out["fail_mapped_reads_num"] == 0


def test_process_bam_reads_fail_fixture() -> None:
    out = process_bam_reads(str(_FAIL_BAM))
    assert out is not None
    assert out["state"] == "fail"
    assert out["sample_id"] == "Sample_104"
    assert out["mapped_reads"] == 2
    assert out["unmapped_reads"] == 0
    assert out["pass_mapped_reads_num"] == 0
    assert out["fail_mapped_reads_num"] == 2


def test_process_bam_reads_counts_match_independent_enumerator() -> None:
    for path in (str(_PASS_BAM), str(_FAIL_BAM)):
        out = process_bam_reads(path)
        assert out is not None
        m, u = _primary_mapped_unmapped(path)
        assert out["mapped_reads"] == m
        assert out["unmapped_reads"] == u


def test_calculate_bam_summary_means_from_process_output() -> None:
    out = process_bam_reads(str(_PASS_BAM))
    assert out is not None
    summary = calculate_bam_summary(out)
    if out["mapped_reads_num"] > 0:
        assert summary["mean_mapped_length"] == pytest.approx(
            out["mapped_bases"] / out["mapped_reads_num"]
        )
    if out["unmapped_reads_num"] > 0:
        assert summary["mean_unmapped_length"] == pytest.approx(
            out["unmapped_bases"] / out["unmapped_reads_num"]
        )
    assert summary["state"] == out["state"]
    assert summary["has_mgmt_reads"] == out.get("has_mgmt_reads", False)


def test_extract_bam_metadata_pass_fixture() -> None:
    meta = extract_bam_metadata(str(_PASS_BAM))
    assert isinstance(meta, BamMetadata)
    assert meta.file_path == str(_PASS_BAM)
    assert meta.file_size == os.stat(_PASS_BAM).st_size
    assert "bam_preprocessing_complete" in meta.processing_steps
    assert meta.extracted_data["sample_id"] == "Sample_104"
    assert meta.extracted_data["state"] == "pass"
    assert meta.extracted_data.get("mapped_reads") == 1


def test_extract_bam_metadata_fail_fixture() -> None:
    meta = extract_bam_metadata(str(_FAIL_BAM))
    assert meta.extracted_data["state"] == "fail"
    assert meta.extracted_data.get("mapped_reads") == 2


def test_extract_sample_id_from_bam_fixtures() -> None:
    assert _extract_sample_id_from_bam(str(_PASS_BAM)) == "Sample_104"
    assert _extract_sample_id_from_bam(str(_FAIL_BAM)) == "Sample_104"


def test_bam_metadata_dataclass_post_init() -> None:
    m = BamMetadata(file_path="/x.bam", file_size=1, creation_time=0.0, extracted_data=None, processing_steps=None)
    assert m.extracted_data == {}
    assert m.processing_steps == []
