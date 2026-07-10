"""Sample delete and archive helpers for admin GUI operations."""

from __future__ import annotations

import tarfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, Tuple

from robin.minknow.sample_id import (
    SAMPLE_IDENTIFIER_MANIFEST_FILENAME,
    validate_custom_sample_id,
)
from robin.utils.docker_fs import chown_tree_to_host_user, remove_tree

# Directories excluded from archives (housekeeping / staging).
HOUSEKEEPING_DIR_NAMES = frozenset({"_locks", "_fusion_staging"})


@dataclass(frozen=True)
class SampleRemovalAssessment:
    """Whether a sample may be deleted or archived."""

    removable: bool
    reason: str = ""


@dataclass(frozen=True)
class SampleLifecycleResult:
    """Outcome of delete or archive."""

    sample_id: str
    action: str
    bytes_processed: int
    archive_path: str = ""


def is_housekeeping_dir(name: str) -> bool:
    """Return True when a directory name should be omitted from archives."""
    return name in HOUSEKEEPING_DIR_NAMES or name.startswith("_")


def assess_sample_for_removal(
    *,
    origin: str,
    active_jobs: int,
    pending_jobs: int,
) -> SampleRemovalAssessment:
    """Return whether delete/archive is allowed for this sample."""
    if origin == "Live":
        return SampleRemovalAssessment(
            removable=False,
            reason="Live runs cannot be deleted or archived.",
        )
    if active_jobs > 0:
        return SampleRemovalAssessment(
            removable=False,
            reason=f"Sample has {active_jobs} active job(s).",
        )
    if pending_jobs > 0:
        return SampleRemovalAssessment(
            removable=False,
            reason=f"Sample has {pending_jobs} pending job(s).",
        )
    return SampleRemovalAssessment(removable=True)


def resolve_sample_dir(work_dir: Path, sample_id: str) -> Path:
    """Resolve and validate a sample directory under work_dir."""
    validated_id = validate_custom_sample_id(sample_id)
    base = work_dir.expanduser().resolve()
    sample_dir = (base / validated_id).resolve()
    try:
        sample_dir.relative_to(base)
    except ValueError as exc:
        raise ValueError(f"Invalid sample path for {validated_id!r}") from exc
    return sample_dir


def _path_is_within(child: Path, parent: Path) -> bool:
    try:
        child.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def validate_archive_destination(destination: Path, work_dir: Path) -> Path:
    """Validate archive destination is an existing directory outside work_dir."""
    dest = destination.expanduser().resolve()
    base = work_dir.expanduser().resolve()

    if not dest.exists():
        raise ValueError("Archive destination does not exist.")
    if not dest.is_dir():
        raise ValueError("Archive destination must be a folder.")
    if not os_access_writable(dest):
        raise ValueError("Archive destination is not writable.")

    if _path_is_within(dest, base):
        raise ValueError("Archive destination must be outside the work directory.")
    if _path_is_within(base, dest):
        raise ValueError("Archive destination cannot contain the work directory.")

    return dest


def os_access_writable(path: Path) -> bool:
    """Best-effort writable check without creating files."""
    import os

    return os.access(path, os.W_OK | os.X_OK)


def iter_archive_members(sample_dir: Path) -> Iterator[Tuple[Path, str]]:
    """Yield (absolute_path, arcname) pairs for files included in an archive."""
    if not sample_dir.is_dir():
        return

    for path in sorted(sample_dir.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(sample_dir)
        if any(is_housekeeping_dir(part) for part in rel.parts[:-1]):
            continue
        arcname = str(Path(sample_dir.name) / rel)
        yield path, arcname


def _directory_size(paths: Iterable[Path]) -> int:
    total = 0
    seen: set[Path] = set()
    for path in paths:
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        try:
            total += resolved.stat().st_size
        except OSError:
            continue
    return total


def _archive_filename(sample_id: str) -> str:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    return f"{sample_id}_{stamp}.tar.gz"


def archive_sample_data(
    work_dir: Path,
    sample_id: str,
    destination_dir: Path,
) -> SampleLifecycleResult:
    """Archive sample folder to destination as tar.gz, then remove the source folder."""
    sample_dir = resolve_sample_dir(work_dir, sample_id)
    if not sample_dir.is_dir():
        raise FileNotFoundError(f"Sample folder not found: {sample_dir}")

    if not chown_tree_to_host_user(sample_dir):
        raise PermissionError(
            f"Cannot read all files under {sample_dir} for archiving. "
            "Root-owned Docker outputs (for example Clair3) could not be "
            "reassigned to the current user. Ensure Docker is available and "
            "retry, or fix ownership manually."
        )

    members = list(iter_archive_members(sample_dir))
    if not members:
        raise ValueError(f"No archivable files found for sample {sample_id!r}.")

    dest = validate_archive_destination(destination_dir, work_dir)
    bytes_before = _directory_size(path for path, _ in members)
    archive_path = dest / _archive_filename(sample_id)

    added_arcnames: set[str] = set()
    with tarfile.open(archive_path, "w:gz") as tar:
        for path, arcname in members:
            if arcname in added_arcnames:
                continue
            tar.add(path, arcname=arcname, recursive=False)
            added_arcnames.add(arcname)

    if not archive_path.is_file() or archive_path.stat().st_size == 0:
        archive_path.unlink(missing_ok=True)
        raise OSError("Archive file was not created successfully.")

    with tarfile.open(archive_path, "r:gz") as tar:
        names = tar.getnames()
        if not names:
            archive_path.unlink(missing_ok=True)
            raise OSError("Archive verification failed: archive is empty.")
        manifest_arc = f"{sample_dir.name}/{SAMPLE_IDENTIFIER_MANIFEST_FILENAME}"
        if (sample_dir / SAMPLE_IDENTIFIER_MANIFEST_FILENAME).is_file():
            if manifest_arc not in names:
                archive_path.unlink(missing_ok=True)
                raise OSError(
                    f"Archive verification failed: missing {SAMPLE_IDENTIFIER_MANIFEST_FILENAME}."
                )

    remove_tree(sample_dir)
    return SampleLifecycleResult(
        sample_id=sample_id,
        action="archive",
        bytes_processed=bytes_before,
        archive_path=str(archive_path),
    )


def delete_sample_data(work_dir: Path, sample_id: str) -> SampleLifecycleResult:
    """Permanently delete a sample folder under work_dir."""
    sample_dir = resolve_sample_dir(work_dir, sample_id)
    if not sample_dir.is_dir():
        raise FileNotFoundError(f"Sample folder not found: {sample_dir}")

    bytes_before = 0
    for path in sample_dir.rglob("*"):
        if path.is_file():
            try:
                bytes_before += path.stat().st_size
            except OSError:
                pass

    chown_tree_to_host_user(sample_dir)
    remove_tree(sample_dir)
    return SampleLifecycleResult(
        sample_id=sample_id,
        action="delete",
        bytes_processed=bytes_before,
    )


def result_to_audit_details(result: SampleLifecycleResult) -> Dict[str, Any]:
    """Serialize a lifecycle result for audit logging."""
    details: Dict[str, Any] = {
        "action": result.action,
        "bytes_processed": result.bytes_processed,
    }
    if result.archive_path:
        details["archive_path"] = result.archive_path
    return details
