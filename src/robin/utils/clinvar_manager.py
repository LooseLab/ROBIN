from __future__ import annotations

import gzip
import hashlib
import json
import logging
import os
import shutil
import subprocess
import tempfile
import urllib.error
import urllib.request
from datetime import datetime, timezone
from email.utils import format_datetime, parsedate_to_datetime
from pathlib import Path
from typing import Any, Optional
import pysam

logger = logging.getLogger("robin.clinvar")

CLINVAR_VCF_GZ_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
)

CLINVAR_VCF_GZ_NAME = "clinvar.vcf.gz"
CLINVAR_VCF_GZ_TBI_NAME = "clinvar.vcf.gz.tbi"
CLINVAR_META_NAME = "clinvar.vcf.gz.meta.json"
SAMPLE_CLINVAR_PROVENANCE_NAME = "clinvar_provenance.json"


def get_resources_dir() -> Path:
    """
    Return the on-disk directory backing the `robin.resources` package.

    Notes:
    - This is expected to be writable in the typical ROBIN dev / workstation setup.
    - If it is not writable, the downloader will raise.
    """

    # Avoid importing `robin` at runtime.
    # `robin/__init__.py` imports analysis modules and can fail if optional
    # dependencies are not installed in minimal environments.
    this_file = Path(__file__).resolve()
    resources_dir = this_file.parent.parent / "resources"
    if not resources_dir.exists():
        raise RuntimeError(f"ROBIN resources directory not found: {resources_dir}")
    return resources_dir


def _file_ok(path: Path) -> bool:
    return path.exists() and path.is_file() and path.stat().st_size > 0


def _read_vcf_meta_lines(gz_path: Path) -> dict[str, str]:
    """Read ClinVar release fields from the VCF header (`##fileDate`, etc.)."""
    fields: dict[str, str] = {}
    try:
        with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.startswith("##"):
                    break
                if line.startswith("##fileDate="):
                    fields["file_date"] = line.split("=", 1)[1].strip()
                elif line.startswith("##source="):
                    fields["source"] = line.split("=", 1)[1].strip()
                elif line.startswith("##reference="):
                    fields["reference"] = line.split("=", 1)[1].strip()
    except (OSError, gzip.BadGzipFile) as exc:
        logger.warning("Could not read ClinVar VCF header from %s: %s", gz_path, exc)
    return fields


def _file_sha256(path: Path, *, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _meta_path(resources_dir: Path) -> Path:
    return resources_dir / CLINVAR_META_NAME


def format_clinvar_version_label(metadata: Optional[dict[str, Any]]) -> str:
    """Human-readable ClinVar release label for CLI output and reports."""
    if not metadata:
        return "ClinVar (version unknown)"
    file_date = str(metadata.get("file_date") or "").strip()
    if file_date:
        return f"ClinVar release {file_date}"
    return "ClinVar (version unknown)"


def refresh_clinvar_metadata(
    *,
    resources_dir: Optional[Path] = None,
    url: str = CLINVAR_VCF_GZ_URL,
    compute_checksum: bool = True,
) -> dict[str, Any]:
    """
    Read ClinVar release metadata from the local VCF and write a sidecar JSON file.

    The canonical release identifier is the VCF `##fileDate` header.
    """
    resources_dir = resources_dir or get_resources_dir()
    gz_path = resources_dir / CLINVAR_VCF_GZ_NAME
    if not _file_ok(gz_path):
        return {}

    header = _read_vcf_meta_lines(gz_path)
    stat = gz_path.stat()
    remote_mtime = _get_remote_last_modified(url)

    metadata: dict[str, Any] = {
        "file_date": header.get("file_date", ""),
        "source": header.get("source", ""),
        "reference": header.get("reference", ""),
        "url": url,
        "local_path": str(gz_path),
        "file_size_bytes": stat.st_size,
        "file_mtime_unix": stat.st_mtime,
        "recorded_at": datetime.now(timezone.utc).isoformat(),
    }
    if remote_mtime is not None:
        metadata["remote_last_modified_unix"] = remote_mtime
        try:
            metadata["remote_last_modified"] = format_datetime(
                datetime.fromtimestamp(remote_mtime, tz=timezone.utc)
            ).replace("+0000", "GMT")
        except (OverflowError, OSError, ValueError):
            logger.debug(
                "Could not format remote Last-Modified timestamp: %s", remote_mtime
            )
    if compute_checksum:
        metadata["sha256"] = _file_sha256(gz_path)

    _meta_path(resources_dir).write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )
    logger.info(
        "ClinVar metadata recorded: %s",
        format_clinvar_version_label(metadata),
    )
    return metadata


def get_clinvar_metadata(
    *,
    resources_dir: Optional[Path] = None,
    url: str = CLINVAR_VCF_GZ_URL,
    compute_checksum: bool = False,
) -> dict[str, Any]:
    """
    Return ClinVar metadata for the installed resource file.

    Uses the sidecar JSON when it matches the current VCF mtime; otherwise
    refreshes from the VCF header.
    """
    resources_dir = resources_dir or get_resources_dir()
    gz_path = resources_dir / CLINVAR_VCF_GZ_NAME
    meta_path = _meta_path(resources_dir)

    if not _file_ok(gz_path):
        return {}

    current_mtime = gz_path.stat().st_mtime
    if meta_path.exists():
        try:
            cached = json.loads(meta_path.read_text(encoding="utf-8"))
            if cached.get("file_mtime_unix") == current_mtime:
                return cached
        except (OSError, json.JSONDecodeError, TypeError):
            pass

    return refresh_clinvar_metadata(
        resources_dir=resources_dir,
        url=url,
        compute_checksum=compute_checksum,
    )


def sample_clinvar_provenance_path(sample_dir: Path | str) -> Path:
    """Per-sample ClinVar provenance written at SNP annotation time."""
    return Path(sample_dir) / "clair3" / SAMPLE_CLINVAR_PROVENANCE_NAME


def record_sample_clinvar_provenance(
    sample_dir: Path | str,
    *,
    clinvar_path: str,
    used_for: str = "snp_annotation",
    update: bool = False,
    resources_dir: Optional[Path] = None,
) -> Optional[Path]:
    """
    Persist which ClinVar release was used for variant annotation on a sample.

    By default this is written once and not overwritten, so later ClinVar updates
    do not erase the release that produced existing annotations.
    """
    out_path = sample_clinvar_provenance_path(sample_dir)
    if out_path.exists() and not update:
        return out_path

    metadata = get_clinvar_metadata(
        resources_dir=resources_dir,
        compute_checksum=True,
    )
    if not metadata:
        return None

    metadata = dict(metadata)
    metadata["clinvar_path"] = clinvar_path
    metadata["used_for"] = used_for
    metadata["recorded_at"] = datetime.now(timezone.utc).isoformat()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    logger.info(
        "Recorded ClinVar provenance for %s: %s",
        sample_dir,
        format_clinvar_version_label(metadata),
    )
    return out_path


def load_sample_clinvar_provenance(
    sample_dir: Path | str,
    *,
    resources_dir: Optional[Path] = None,
) -> dict[str, Any]:
    """
    Load ClinVar provenance for a sample.

    Prefers the per-sample file written during SNP annotation; falls back to the
    currently installed ClinVar resource metadata.
    """
    provenance_path = sample_clinvar_provenance_path(sample_dir)
    if provenance_path.exists():
        try:
            return json.loads(provenance_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError, TypeError):
            logger.warning("Could not read ClinVar provenance at %s", provenance_path)

    return get_clinvar_metadata(resources_dir=resources_dir, compute_checksum=False)


def sample_snp_reannotation_inputs_ready(sample_dir: Path | str) -> bool:
    """Return True when Clair3 outputs exist for annotation-only reruns."""
    clair_dir = Path(sample_dir) / "clair3"
    return (
        (clair_dir / "output_done.vcf.gz").is_file()
        and (clair_dir / "output_indel_done.vcf.gz").is_file()
    )


def compare_sample_clinvar_to_installed(
    sample_dir: Path | str,
    *,
    resources_dir: Optional[Path] = None,
) -> dict[str, Any]:
    """
    Compare the ClinVar release used to annotate a sample with the installed resource.

    Returns labels for the GUI plus ``is_stale`` when the installed ClinVar differs
    from the per-sample provenance recorded at annotation time.
    """
    installed = get_clinvar_metadata(
        resources_dir=resources_dir,
        compute_checksum=True,
    )
    provenance_path = sample_clinvar_provenance_path(sample_dir)
    has_sample_record = provenance_path.exists()
    sample_meta: dict[str, Any] = {}
    if has_sample_record:
        try:
            sample_meta = json.loads(provenance_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError, TypeError):
            logger.warning("Could not read ClinVar provenance at %s", provenance_path)
            has_sample_record = False

    installed_date = str(installed.get("file_date") or "").strip()
    sample_date = str(sample_meta.get("file_date") or "").strip()
    installed_sha = str(installed.get("sha256") or "").strip()
    sample_sha = str(sample_meta.get("sha256") or "").strip()

    is_stale = False
    if has_sample_record and installed:
        if installed_sha and sample_sha:
            is_stale = installed_sha != sample_sha
        elif installed_date and sample_date:
            is_stale = installed_date != sample_date

    return {
        "installed_metadata": installed,
        "sample_metadata": sample_meta,
        "has_sample_record": has_sample_record,
        "installed_label": format_clinvar_version_label(installed),
        "sample_label": (
            format_clinvar_version_label(sample_meta)
            if has_sample_record
            else "Annotation release not recorded"
        ),
        "is_stale": is_stale,
        "can_reannotate": sample_snp_reannotation_inputs_ready(sample_dir),
    }


def _download_url_to_file(url: str, target_path: Path, *, timeout_s: int = 600) -> None:
    """
    Download `url` to `target_path` atomically (temp file + rename).
    """

    target_path.parent.mkdir(parents=True, exist_ok=True)

    # Use a temp file in the same directory so rename is atomic.
    with tempfile.NamedTemporaryFile(
        mode="wb", suffix=".part", prefix=target_path.name + ".", dir=str(target_path.parent), delete=False
    ) as tmp:
        tmp_path = Path(tmp.name)
    try:
        msg = f"Downloading ClinVar from {url}"
        print(msg)
        logger.info("Downloading ClinVar from %s", url)
        req = urllib.request.Request(url, headers={"User-Agent": "robin-clinvar/1.0"})
        with urllib.request.urlopen(req, timeout=timeout_s) as resp:
            total: Optional[int] = None
            try:
                cl = resp.headers.get("Content-Length")
                if cl:
                    total = int(cl)
            except Exception:
                total = None

            # Import click lazily to keep module usable in non-CLI contexts.
            try:
                import click  # type: ignore
            except Exception:
                click = None  # type: ignore

            chunk_size = 1024 * 1024  # 1 MiB
            downloaded = 0

            if click is not None and total and total > 0:
                with click.progressbar(
                    length=total,
                    label="Downloading ClinVar",
                    show_eta=True,
                    show_percent=True,
                ) as bar:
                    with open(tmp_path, "wb") as f:
                        while True:
                            chunk = resp.read(chunk_size)
                            if not chunk:
                                break
                            f.write(chunk)
                            downloaded += len(chunk)
                            bar.update(len(chunk))
            else:
                # Fallback: stream and emit sparse size updates.
                last_reported_mb = -1
                with open(tmp_path, "wb") as f:
                    while True:
                        chunk = resp.read(chunk_size)
                        if not chunk:
                            break
                        f.write(chunk)
                        downloaded += len(chunk)
                        if click is not None:
                            mb = downloaded // (1024 * 1024)
                            if mb // 16 != last_reported_mb // 16:
                                last_reported_mb = mb
                                click.echo(f"Downloading ClinVar: {mb} MiB downloaded...")

        tmp_path.replace(target_path)
        print(f"ClinVar download complete: {target_path}")
        logger.info("Downloaded to %s", target_path)
    except Exception:
        try:
            if tmp_path.exists():
                tmp_path.unlink()
        except OSError:
            pass
        raise


def _verify_clinvar_tabix_index(gz_path: Path) -> bool:
    """
    Smoke-test that the ClinVar tabix index can fetch a known GRCh38 site.

    Catches stale or corrupt indices that still exist on disk (SnpSift then fails
    with ``Invalid GZIP header`` while exiting 0).
    """
    try:
        pos = 208248388
        with pysam.TabixFile(str(gz_path)) as tabix_file:
            rows = list(tabix_file.fetch("2", pos - 1, pos))
        return bool(rows)
    except Exception as exc:
        logger.debug("ClinVar tabix index verification failed for %s: %s", gz_path, exc)
        return False


def _tabix_index_needs_rebuild(gz_path: Path, tbi_path: Path) -> bool:
    """True when the index is missing or older than the bgzipped VCF."""
    if not _file_ok(tbi_path):
        return True
    if not _file_ok(gz_path):
        return False
    return tbi_path.stat().st_mtime < gz_path.stat().st_mtime


def _build_tabix_index(gz_path: Path, tbi_path: Path) -> None:
    """Create or overwrite a tabix index for a bgzipped VCF."""
    print(f"Creating tabix index for {gz_path.name}")
    logger.info("Creating tabix index for %s", gz_path)
    last_err: Optional[BaseException] = None
    try:
        pysam.tabix_index(
            str(gz_path),
            preset="vcf",
            force=True,
            keep_original=True,
        )
    except Exception as exc:
        last_err = exc
        logger.warning("pysam.tabix_index failed for %s: %s", gz_path, exc)
        tabix_bin = shutil.which("tabix")
        if tabix_bin:
            try:
                subprocess.run(
                    [tabix_bin, "-p", "vcf", str(gz_path)],
                    check=True,
                    capture_output=True,
                    text=True,
                )
            except subprocess.CalledProcessError as sub_exc:
                stderr = (sub_exc.stderr or "").strip()
                stdout = (sub_exc.stdout or "").strip()
                detail = stderr or stdout or str(sub_exc)
                raise RuntimeError(
                    f"Tabix indexing failed for {gz_path} (pysam: {exc!r}; "
                    f"tabix CLI: {detail})"
                ) from sub_exc
        else:
            raise RuntimeError(
                f"Tabix indexing failed for {gz_path}: {exc!r}. "
                "Install htslib (tabix/bgzip) or fix the VCF: it must be "
                "block-gzipped (bgzip), not plain gzip, and a valid sorted VCF."
            ) from exc

    if not _file_ok(tbi_path):
        hint = (
            f"Tabix index was not created at {tbi_path}. "
            "If the VCF was recompressed with gzip instead of bgzip, run: "
            f"bgzip -d {gz_path.name} && bgzip {gz_path.with_suffix('').name} "
            f"&& tabix -p vcf {gz_path.name}"
        )
        if last_err:
            raise RuntimeError(hint) from last_err
        raise RuntimeError(hint)


def _ensure_tabix_index(gz_path: Path, tbi_path: Path) -> None:
    """
    Ensure a tabix index exists for a bgzipped VCF.

    Rebuilds when the index is missing, older than the VCF, or fails a smoke query.
    """
    if _file_ok(tbi_path) and _file_ok(gz_path):
        if not _tabix_index_needs_rebuild(gz_path, tbi_path):
            if _verify_clinvar_tabix_index(gz_path):
                return
            logger.warning(
                "ClinVar tabix index at %s failed verification; rebuilding.",
                tbi_path,
            )
        else:
            logger.info(
                "ClinVar tabix index %s is older than %s; rebuilding.",
                tbi_path.name,
                gz_path.name,
            )
    elif _file_ok(tbi_path):
        if _verify_clinvar_tabix_index(gz_path):
            return
        logger.warning(
            "ClinVar tabix index at %s failed verification; rebuilding.",
            tbi_path,
        )

    if _file_ok(tbi_path):
        try:
            tbi_path.unlink()
        except OSError as exc:
            logger.warning("Could not remove stale ClinVar tabix index %s: %s", tbi_path, exc)

    _build_tabix_index(gz_path, tbi_path)

    if not _verify_clinvar_tabix_index(gz_path):
        raise RuntimeError(
            f"ClinVar tabix index at {tbi_path} was rebuilt but verification still failed."
        )


def ensure_clinvar_files(
    *,
    resources_dir: Optional[Path] = None,
    download_if_missing: bool = True,
    url: str = CLINVAR_VCF_GZ_URL,
) -> None:
    """
    Ensure required ClinVar files exist in `robin.resources`:
    - `clinvar.vcf.gz`
    - `clinvar.vcf.gz.tbi`

    Strategy:
    - Download `clinvar.vcf.gz` when missing.
    - Create tabix index (`.tbi`) when missing.
    """

    resources_dir = resources_dir or get_resources_dir()
    gz_path = resources_dir / CLINVAR_VCF_GZ_NAME
    tbi_path = resources_dir / CLINVAR_VCF_GZ_TBI_NAME

    # If everything is already present, do nothing.
    if _file_ok(gz_path) and _file_ok(tbi_path):
        return

    gz_ok = _file_ok(gz_path)
    tbi_ok = _file_ok(tbi_path)

    print(
        "ClinVar status check: "
        f"{gz_path.name}={'OK' if gz_ok else 'MISSING'}; "
        f"{tbi_path.name}={'OK' if tbi_ok else 'MISSING'}"
    )

    if not gz_ok:
        if not download_if_missing:
            raise FileNotFoundError(
                f"Missing ClinVar file in {resources_dir}: expected {gz_path.name}"
            )
        _download_url_to_file(url, gz_path)
    _ensure_tabix_index(gz_path, tbi_path)
    if _file_ok(gz_path):
        refresh_clinvar_metadata(
            resources_dir=resources_dir,
            url=url,
            compute_checksum=True,
        )


def _get_remote_last_modified(url: str, *, timeout_s: int = 30) -> Optional[float]:
    """
    Best-effort: fetch HTTP Last-Modified as a unix timestamp (seconds).
    Returns None when header is unavailable or request fails.
    """

    req = urllib.request.Request(url, headers={"User-Agent": "robin-clinvar/1.0"})
    req.method = "HEAD"
    try:
        with urllib.request.urlopen(req, timeout=timeout_s) as resp:
            last_modified = resp.headers.get("Last-Modified")
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
        return None

    if not last_modified:
        return None

    try:
        dt = parsedate_to_datetime(last_modified)
        return dt.timestamp()
    except Exception:
        return None


def update_clinvar_if_newer(
    *,
    resources_dir: Optional[Path] = None,
    url: str = CLINVAR_VCF_GZ_URL,
    download_if_missing: bool = True,
    timeout_s: int = 600,
) -> bool:
    """
    Update ClinVar if NCBI reports a newer `Last-Modified` than the local file.

    Returns:
        True if an update was performed, else False.
    """

    resources_dir = resources_dir or get_resources_dir()
    gz_path = resources_dir / CLINVAR_VCF_GZ_NAME
    tbi_path = resources_dir / CLINVAR_VCF_GZ_TBI_NAME

    # Ensure local ClinVar artifacts exist.
    ensure_clinvar_files(
        resources_dir=resources_dir,
        download_if_missing=download_if_missing,
        url=url,
    )

    if not _file_ok(gz_path):
        # Should not happen due to ensure_clinvar_files, but keep defensive.
        if not download_if_missing:
            return False
        _download_url_to_file(url, gz_path, timeout_s=timeout_s)

    local_mtime = gz_path.stat().st_mtime
    remote_mtime = _get_remote_last_modified(url)

    if remote_mtime is None:
        logger.info(
            "ClinVar update skipped: remote Last-Modified unavailable for %s", url
        )
        print("ClinVar update skipped: remote Last-Modified unavailable.")
        refresh_clinvar_metadata(
            resources_dir=resources_dir,
            url=url,
            compute_checksum=False,
        )
        return False

    # Add a small tolerance to avoid re-downloading due to timestamp rounding.
    if remote_mtime <= local_mtime + 1:
        logger.info("ClinVar already up to date (local=%s, remote=%s).", local_mtime, remote_mtime)
        print("ClinVar already up to date.")
        refresh_clinvar_metadata(
            resources_dir=resources_dir,
            url=url,
            compute_checksum=False,
        )
        return False

    # Download newer gz to temp then replace.
    tmp_path = Path(str(gz_path) + ".new")
    try:
        _download_url_to_file(url, tmp_path, timeout_s=timeout_s)
        tmp_path.replace(gz_path)
        _ensure_tabix_index(gz_path, tbi_path)
    finally:
        try:
            if tmp_path.exists():
                tmp_path.unlink()
        except OSError:
            pass

    logger.info("ClinVar updated successfully.")
    refresh_clinvar_metadata(
        resources_dir=resources_dir,
        url=url,
        compute_checksum=True,
    )
    return True

