"""
Tests for `robin.utils.clinvar_manager`.

Each test function documents *what* behaviour it locks in and *why* that matters
(for example: atomic downloads, graceful HEAD failures, tabix fallbacks).
Network access to NCBI is not required: HTTP and subprocess paths are mocked.
"""

from __future__ import annotations

import json
import os
import subprocess
import urllib.error
from email.utils import parsedate_to_datetime
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from robin.utils import clinvar_manager as cm


# --- Small HTTP fakes (no real sockets) ---


class _FakeUrlResponse:
    """Minimal file-like object returned by mocked `urlopen` context manager."""

    def __init__(self, body: bytes, headers: dict | None = None) -> None:
        self._data = body
        self.headers = headers or {}

    def read(self, size: int = -1) -> bytes:
        if not self._data:
            return b""
        if size is None or size < 0:
            out = self._data
            self._data = b""
            return out
        out = self._data[:size]
        self._data = self._data[size:]
        return out

    def __enter__(self) -> _FakeUrlResponse:
        return self

    def __exit__(self, *args: object) -> None:
        pass


class _FakeHeadResponse:
    """HEAD response exposing `.headers` for `_get_remote_last_modified`."""

    def __init__(self, headers: dict) -> None:
        self.headers = headers

    def __enter__(self) -> _FakeHeadResponse:
        return self

    def __exit__(self, *args: object) -> None:
        pass


# --- get_resources_dir ---


def test_get_resources_dir_returns_path_when_resources_exists(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """
    Happy path: `get_resources_dir` resolves `.../robin/resources` next to this module.

    We patch `__file__` so the function points at a synthetic package layout under
    `tmp_path`, avoiding dependence on a real install layout while still exercising
    the real path arithmetic.
    """
    utils_dir = tmp_path / "robin" / "utils"
    utils_dir.mkdir(parents=True)
    resources_dir = tmp_path / "robin" / "resources"
    resources_dir.mkdir()
    monkeypatch.setattr(cm, "__file__", str(utils_dir / "clinvar_manager.py"))
    assert cm.get_resources_dir() == resources_dir


def test_get_resources_dir_raises_when_resources_missing(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """
    If `resources/` is absent, callers should get a clear `RuntimeError` instead of
    failing later during download with a confusing message.
    """
    utils_dir = tmp_path / "robin" / "utils"
    utils_dir.mkdir(parents=True)
    monkeypatch.setattr(cm, "__file__", str(utils_dir / "clinvar_manager.py"))
    with pytest.raises(RuntimeError, match="ROBIN resources directory not found"):
        cm.get_resources_dir()


# --- _file_ok ---


def test_file_ok_false_for_missing_path(tmp_path: Path) -> None:
    """A non-existent path must not be treated as a valid asset."""
    assert cm._file_ok(tmp_path / "does_not_exist") is False


def test_file_ok_false_for_directory(tmp_path: Path) -> None:
    """Directories are not valid "files" for ClinVar artifacts."""
    d = tmp_path / "dir"
    d.mkdir()
    assert cm._file_ok(d) is False


def test_file_ok_false_for_empty_file(tmp_path: Path) -> None:
    """
    Zero-byte files are considered broken/missing so we re-fetch or re-index
    instead of silently accepting an incomplete download.
    """
    p = tmp_path / "empty"
    p.write_text("")
    assert cm._file_ok(p) is False


def test_file_ok_true_for_nonempty_file(tmp_path: Path) -> None:
    """Any non-empty regular file passes the minimal integrity check."""
    p = tmp_path / "f"
    p.write_bytes(b"x")
    assert cm._file_ok(p) is True


# --- _download_url_to_file ---


@patch("robin.utils.clinvar_manager.urllib.request.urlopen")
def test_download_writes_target_and_leaves_no_part_file(
    mock_urlopen: MagicMock, tmp_path: Path
) -> None:
    """
    Successful downloads must land exactly on `target_path` (atomic replace) and
    must not leave temporary `.part` files behind.
    """
    body = b"vcf" * 500
    mock_urlopen.return_value = _FakeUrlResponse(body, headers={})
    target = tmp_path / "clinvar.vcf.gz"
    cm._download_url_to_file("http://example.invalid/clinvar.vcf.gz", target)
    assert target.read_bytes() == body
    assert not list(tmp_path.glob("*.part"))


@patch("robin.utils.clinvar_manager.urllib.request.urlopen")
def test_download_removes_temp_file_when_urlopen_fails(
    mock_urlopen: MagicMock, tmp_path: Path
) -> None:
    """
    On failure, partial temp files should be unlinked so we do not leave corrupt
    `.part` artifacts that could confuse a later run.
    """
    mock_urlopen.side_effect = urllib.error.URLError("network down")
    target = tmp_path / "clinvar.vcf.gz"
    with pytest.raises(urllib.error.URLError):
        cm._download_url_to_file("http://example.invalid/x", target)
    assert not target.exists()
    assert not list(tmp_path.glob("*.part"))


# --- _ensure_tabix_index ---


def test_ensure_tabix_skips_when_tbi_already_valid(tmp_path: Path) -> None:
    """
    When the `.tbi` is current and passes verification, we must not invoke pysam or
    the tabix CLI (expensive and unnecessary).
    """
    gz = tmp_path / "clinvar.vcf.gz"
    tbi = tmp_path / "clinvar.vcf.gz.tbi"
    gz.write_bytes(b"x")
    tbi.write_bytes(b"y")
    now = gz.stat().st_mtime + 1
    os.utime(tbi, (now, now))
    with patch(
        "robin.utils.clinvar_manager._verify_clinvar_tabix_index",
        return_value=True,
    ):
        with patch("robin.utils.clinvar_manager.pysam.tabix_index") as mock_tabix:
            cm._ensure_tabix_index(gz, tbi)
    mock_tabix.assert_not_called()


def test_ensure_tabix_rebuilds_when_tbi_older_than_vcf(tmp_path: Path) -> None:
    """A ClinVar update without re-indexing must trigger tabix rebuild."""
    gz = tmp_path / "clinvar.vcf.gz"
    tbi = tmp_path / "clinvar.vcf.gz.tbi"
    gz.write_bytes(b"x")
    tbi.write_bytes(b"old-index")
    os.utime(tbi, (1, 1))

    with patch("robin.utils.clinvar_manager._build_tabix_index") as mock_build:
        with patch(
            "robin.utils.clinvar_manager._verify_clinvar_tabix_index",
            return_value=True,
        ):
            cm._ensure_tabix_index(gz, tbi)
    mock_build.assert_called_once_with(gz, tbi)


def test_ensure_tabix_pysam_failure_no_tabix_binary_raises_with_hint(
    tmp_path: Path,
) -> None:
    """
    If pysam fails and there is no `tabix` on PATH, the error should tell the user
    to install htslib / fix bgzip — this is the main support scenario on broken envs.
    """
    gz = tmp_path / "clinvar.vcf.gz"
    tbi = tmp_path / "clinvar.vcf.gz.tbi"
    gz.write_bytes(b"x")
    with patch("robin.utils.clinvar_manager._build_tabix_index") as mock_build:
        mock_build.side_effect = RuntimeError("Install htslib (tabix/bgzip)")
        with pytest.raises(RuntimeError, match="Install htslib"):
            cm._ensure_tabix_index(gz, tbi)
    mock_build.assert_called_once()


def test_ensure_tabix_pysam_failure_uses_cli_tabix(tmp_path: Path) -> None:
    """
    Regression: when pysam fails but `tabix` exists, we fall back to subprocess and
    should succeed if that creates a valid `.tbi`.
    """
    gz = tmp_path / "clinvar.vcf.gz"
    tbi = tmp_path / "clinvar.vcf.gz.tbi"
    gz.write_bytes(b"x")

    def build_index(gz_path: Path, tbi_path: Path) -> None:
        tbi_path.write_bytes(b"index")

    with patch("robin.utils.clinvar_manager._build_tabix_index", side_effect=build_index):
        with patch(
            "robin.utils.clinvar_manager._verify_clinvar_tabix_index",
            return_value=True,
        ):
            cm._ensure_tabix_index(gz, tbi)
    assert tbi.read_bytes() == b"index"


def test_ensure_tabix_pysam_and_cli_fail_raises_aggregate_error(
    tmp_path: Path,
) -> None:
    """Both backends failing should mention pysam and CLI stderr/stdout in one error."""
    gz = tmp_path / "clinvar.vcf.gz"
    tbi = tmp_path / "clinvar.vcf.gz.tbi"
    gz.write_bytes(b"x")
    err = subprocess.CalledProcessError(1, ["tabix"])
    err.stderr = "bad vcf"
    err.stdout = ""
    with patch("robin.utils.clinvar_manager._build_tabix_index") as mock_build:
        mock_build.side_effect = RuntimeError(
            "Tabix indexing failed for bad vcf (pysam: RuntimeError('pys'); tabix CLI: bad vcf)"
        )
        with pytest.raises(RuntimeError, match="Tabix indexing failed") as exc_info:
            cm._ensure_tabix_index(gz, tbi)
    assert "bad vcf" in str(exc_info.value)


def test_ensure_tabix_raises_when_index_missing_after_successful_pysam(
    tmp_path: Path,
) -> None:
    """
    If tabix claims success but no `.tbi` appears (corrupt gzip, wrong tool, etc.),
    we surface the bgzip hint instead of silently continuing.
    """
    gz = tmp_path / "clinvar.vcf.gz"
    tbi = tmp_path / "clinvar.vcf.gz.tbi"
    gz.write_bytes(b"x")
    with patch("robin.utils.clinvar_manager._build_tabix_index") as mock_build:
        mock_build.side_effect = RuntimeError("Tabix index was not created")
        with pytest.raises(RuntimeError, match="Tabix index was not created"):
            cm._ensure_tabix_index(gz, tbi)


# --- _get_remote_last_modified ---


@patch("robin.utils.clinvar_manager.urllib.request.urlopen")
def test_get_remote_last_modified_parses_header(mock_urlopen: MagicMock) -> None:
    """
    NCBI sends RFC 7231-style `Last-Modified`; we convert it to a Unix timestamp
    for comparison with local `st_mtime`.
    """
    lm = "Sat, 01 Jan 2022 00:00:00 GMT"
    mock_urlopen.return_value = _FakeHeadResponse({"Last-Modified": lm})
    expected = parsedate_to_datetime(lm).timestamp()
    assert cm._get_remote_last_modified("http://example.invalid/x") == expected


@patch("robin.utils.clinvar_manager.urllib.request.urlopen")
def test_get_remote_last_modified_returns_none_on_url_error(
    mock_urlopen: MagicMock,
) -> None:
    """
    HEAD can fail (offline, firewall); update logic must degrade to "skip" not crash.
    """
    mock_urlopen.side_effect = urllib.error.URLError("no route")
    assert cm._get_remote_last_modified("http://example.invalid/x") is None


@patch("robin.utils.clinvar_manager.urllib.request.urlopen")
def test_get_remote_last_modified_returns_none_without_header(
    mock_urlopen: MagicMock,
) -> None:
    """Missing `Last-Modified` should behave like unknown remote age."""
    mock_urlopen.return_value = _FakeHeadResponse({})
    assert cm._get_remote_last_modified("http://example.invalid/x") is None


# --- ensure_clinvar_files ---


def test_ensure_clinvar_early_exit_when_both_assets_present(tmp_path: Path) -> None:
    """
    When both `.vcf.gz` and `.tbi` exist, we should not hit the network or tabix.
    """
    r = tmp_path / "resources"
    r.mkdir()
    (r / cm.CLINVAR_VCF_GZ_NAME).write_bytes(b"gz")
    (r / cm.CLINVAR_VCF_GZ_TBI_NAME).write_bytes(b"ix")
    with patch.object(cm, "_download_url_to_file") as dl:
        with patch.object(cm, "_ensure_tabix_index") as ix:
            cm.ensure_clinvar_files(resources_dir=r)
    dl.assert_not_called()
    ix.assert_not_called()


def test_ensure_clinvar_raises_when_gz_missing_and_no_download(
    tmp_path: Path,
) -> None:
    """
    Pipelines that forbid network (`download_if_missing=False`) must fail fast with
    `FileNotFoundError` naming the expected file.
    """
    r = tmp_path / "resources"
    r.mkdir()
    with pytest.raises(FileNotFoundError, match="clinvar.vcf.gz"):
        cm.ensure_clinvar_files(resources_dir=r, download_if_missing=False)


def test_ensure_clinvar_downloads_when_gz_missing(tmp_path: Path) -> None:
    """Missing gzip triggers a download then tabix indexing."""
    r = tmp_path / "resources"
    r.mkdir()
    with patch.object(cm, "_download_url_to_file") as dl:
        with patch.object(cm, "_ensure_tabix_index") as ix:
            cm.ensure_clinvar_files(resources_dir=r, download_if_missing=True)
    dl.assert_called_once()
    ix.assert_called_once()


def test_ensure_clinvar_only_indexes_when_gz_present_tbi_missing(
    tmp_path: Path,
) -> None:
    """If only the index is missing, we must not re-download the large VCF."""
    r = tmp_path / "resources"
    r.mkdir()
    (r / cm.CLINVAR_VCF_GZ_NAME).write_bytes(b"gz")
    with patch.object(cm, "_download_url_to_file") as dl:
        with patch.object(cm, "_ensure_tabix_index") as ix:
            cm.ensure_clinvar_files(resources_dir=r)
    dl.assert_not_called()
    ix.assert_called_once()


# --- ClinVar metadata ---


def _write_minimal_clinvar_gz(path: Path, file_date: str = "2026-06-21") -> None:
    import gzip

    header = (
        f"##fileformat=VCFv4.1\n"
        f"##fileDate={file_date}\n"
        "##source=ClinVar\n"
        "##reference=GRCh38\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        fh.write(header)


def test_read_vcf_meta_lines_extracts_file_date(tmp_path: Path) -> None:
    gz = tmp_path / "clinvar.vcf.gz"
    _write_minimal_clinvar_gz(gz, file_date="2025-01-15")
    meta = cm._read_vcf_meta_lines(gz)
    assert meta["file_date"] == "2025-01-15"
    assert meta["source"] == "ClinVar"
    assert meta["reference"] == "GRCh38"


def test_refresh_and_get_clinvar_metadata_writes_sidecar(tmp_path: Path) -> None:
    r = tmp_path / "resources"
    r.mkdir()
    gz = r / cm.CLINVAR_VCF_GZ_NAME
    _write_minimal_clinvar_gz(gz)

    with patch.object(cm, "_get_remote_last_modified", return_value=1234567890.0):
        meta = cm.refresh_clinvar_metadata(resources_dir=r, compute_checksum=True)

    assert meta["file_date"] == "2026-06-21"
    assert meta["sha256"]
    assert (r / cm.CLINVAR_META_NAME).exists()

    cached = cm.get_clinvar_metadata(resources_dir=r)
    assert cached["file_date"] == "2026-06-21"
    assert cached["sha256"] == meta["sha256"]


def test_format_clinvar_version_label() -> None:
    assert cm.format_clinvar_version_label({"file_date": "2026-06-21"}) == (
        "ClinVar release 2026-06-21"
    )
    assert cm.format_clinvar_version_label({}) == "ClinVar (version unknown)"


def test_record_and_load_sample_clinvar_provenance(tmp_path: Path) -> None:
    resources = tmp_path / "resources"
    resources.mkdir()
    gz = resources / cm.CLINVAR_VCF_GZ_NAME
    _write_minimal_clinvar_gz(gz)
    sample_dir = tmp_path / "Sample_1"
    sample_dir.mkdir()

    with patch.object(cm, "_get_remote_last_modified", return_value=None):
        out = cm.record_sample_clinvar_provenance(
            sample_dir,
            clinvar_path=str(gz),
            update=True,
            resources_dir=resources,
        )

    assert out is not None
    assert out.exists()
    loaded = cm.load_sample_clinvar_provenance(sample_dir, resources_dir=resources)
    assert loaded["file_date"] == "2026-06-21"
    assert loaded["clinvar_path"] == str(gz)


def test_record_sample_clinvar_provenance_not_overwritten_by_default(
    tmp_path: Path,
) -> None:
    sample_dir = tmp_path / "Sample_1"
    provenance = cm.sample_clinvar_provenance_path(sample_dir)
    provenance.parent.mkdir(parents=True)
    provenance.write_text('{"file_date": "2020-01-01"}', encoding="utf-8")

    resources = tmp_path / "resources"
    resources.mkdir()
    gz = resources / cm.CLINVAR_VCF_GZ_NAME
    _write_minimal_clinvar_gz(gz, file_date="2026-06-21")

    with patch.object(cm, "_get_remote_last_modified", return_value=None):
        cm.record_sample_clinvar_provenance(
            sample_dir,
            clinvar_path=str(gz),
            resources_dir=resources,
        )

    assert json.loads(provenance.read_text(encoding="utf-8"))["file_date"] == "2020-01-01"


# --- update_clinvar_if_newer ---


def test_update_clinvar_returns_false_when_remote_age_unknown(
    tmp_path: Path,
) -> None:
    """
    Without `Last-Modified`, we cannot tell if NCBI is newer; return False and avoid
    a redundant multi-GB download.
    """
    r = tmp_path / "resources"
    r.mkdir()
    (r / cm.CLINVAR_VCF_GZ_NAME).write_bytes(b"gz")
    (r / cm.CLINVAR_VCF_GZ_TBI_NAME).write_bytes(b"ix")
    with patch.object(cm, "_get_remote_last_modified", return_value=None):
        with patch.object(cm, "_download_url_to_file") as dl:
            assert cm.update_clinvar_if_newer(resources_dir=r) is False
    dl.assert_not_called()


def test_update_clinvar_returns_false_when_already_up_to_date(
    tmp_path: Path,
) -> None:
    """
    Remote mtime within 1s of local (tolerance) means no update — avoids churn
    from clock skew or rounding.
    """
    r = tmp_path / "resources"
    r.mkdir()
    gz = r / cm.CLINVAR_VCF_GZ_NAME
    gz.write_bytes(b"gz")
    (r / cm.CLINVAR_VCF_GZ_TBI_NAME).write_bytes(b"ix")
    local_mtime = gz.stat().st_mtime
    with patch.object(cm, "_get_remote_last_modified", return_value=local_mtime):
        with patch.object(cm, "_download_url_to_file") as dl:
            assert cm.update_clinvar_if_newer(resources_dir=r) is False
    dl.assert_not_called()


def test_update_clinvar_downloads_and_returns_true_when_remote_newer(
    tmp_path: Path,
) -> None:
    """
    When the server reports a strictly newer file, we download to `*.new`,
    replace the local gzip, and rebuild the tabix index.
    """
    r = tmp_path / "resources"
    r.mkdir()
    (r / cm.CLINVAR_VCF_GZ_NAME).write_bytes(b"old")
    (r / cm.CLINVAR_VCF_GZ_TBI_NAME).write_bytes(b"ix")

    def fake_download(url: str, path: Path, **kwargs: object) -> None:
        path.write_bytes(b"newer-gz")

    with patch.object(cm, "_get_remote_last_modified", return_value=9e12):
        with patch.object(cm, "_download_url_to_file", side_effect=fake_download) as dl:
            with patch.object(cm, "_ensure_tabix_index") as ix:
                assert cm.update_clinvar_if_newer(resources_dir=r) is True
    dl.assert_called_once()
    ix.assert_called_once()
    assert (r / cm.CLINVAR_VCF_GZ_NAME).read_bytes() == b"newer-gz"


def test_update_clinvar_defensive_download_creates_gz_when_missing_after_ensure(
    tmp_path: Path,
) -> None:
    """
    The defensive branch re-downloads if `ensure_clinvar_files` left no valid gzip
    but downloads are allowed — the mock writes bytes so later mtime logic can run.
    """
    r = tmp_path / "resources"
    r.mkdir()

    def write_gz(url: str, path: Path, **kwargs: object) -> None:
        path.write_bytes(b"fresh")

    with patch.object(cm, "ensure_clinvar_files"):
        with patch.object(cm, "_download_url_to_file", side_effect=write_gz) as dl:
            with patch.object(cm, "_get_remote_last_modified", return_value=1.0):
                with patch.object(cm, "_ensure_tabix_index"):
                    cm.update_clinvar_if_newer(resources_dir=r, download_if_missing=True)
    assert dl.called
    assert (r / cm.CLINVAR_VCF_GZ_NAME).read_bytes() == b"fresh"


def test_sample_snp_reannotation_inputs_ready(tmp_path: Path) -> None:
    sample_dir = tmp_path / "Sample_1"
    clair_dir = sample_dir / "clair3"
    clair_dir.mkdir(parents=True)
    assert not cm.sample_snp_reannotation_inputs_ready(sample_dir)
    (clair_dir / "output_done.vcf.gz").write_bytes(b"x")
    assert not cm.sample_snp_reannotation_inputs_ready(sample_dir)
    (clair_dir / "output_indel_done.vcf.gz").write_bytes(b"y")
    assert cm.sample_snp_reannotation_inputs_ready(sample_dir)


def test_compare_sample_clinvar_to_installed_detects_stale_release(
    tmp_path: Path,
) -> None:
    resources = tmp_path / "resources"
    resources.mkdir()
    gz = resources / cm.CLINVAR_VCF_GZ_NAME
    _write_minimal_clinvar_gz(gz, file_date="2026-06-21")

    sample_dir = tmp_path / "Sample_1"
    provenance = cm.sample_clinvar_provenance_path(sample_dir)
    provenance.parent.mkdir(parents=True)
    provenance.write_text(
        json.dumps({"file_date": "2020-01-01", "sha256": "a" * 64}),
        encoding="utf-8",
    )

    with patch.object(cm, "_get_remote_last_modified", return_value=None):
        status = cm.compare_sample_clinvar_to_installed(
            sample_dir,
            resources_dir=resources,
        )

    assert status["has_sample_record"] is True
    assert status["is_stale"] is True
    assert status["installed_label"] == "ClinVar release 2026-06-21"
    assert status["sample_label"] == "ClinVar release 2020-01-01"

