"""Push live readfish target updates when master BED files change."""

from __future__ import annotations

import logging
import os
import re
import threading
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import tomllib
import tomli_w

LOGGER = logging.getLogger(__name__)

_MASTER_BED_FILENAME_RE = re.compile(r"^master_\d{3}\.bed$")


@dataclass(frozen=True)
class ReadfishLiveSession:
    """Active readfish run registered for live target updates."""

    sample_id: str
    base_toml_path: str
    region_name: str = "robin_panel"
    last_master_bed_path: Optional[str] = None


def live_toml_path(base_toml: str | Path) -> Path:
    """Return the readfish live-update path for a base experiment TOML."""
    return Path(f"{Path(base_toml).expanduser()}{'_live'}")


def find_latest_master_bed(*, work_dir: str | Path, sample_id: str) -> Optional[Path]:
    """Return the newest ``master_NNN.bed`` under the sample bed directory."""
    bed_dir = Path(work_dir).expanduser() / sample_id / "bed_files"
    if not bed_dir.is_dir():
        return None

    candidates: list[tuple[int, Path]] = []
    for bed_file in bed_dir.glob("master_*.bed"):
        if not _MASTER_BED_FILENAME_RE.match(bed_file.name):
            continue
        if not bed_file.is_file():
            continue
        try:
            counter = int(bed_file.stem.rsplit("_", 1)[-1])
        except ValueError:
            continue
        candidates.append((counter, bed_file))

    if not candidates:
        return None
    return max(candidates, key=lambda item: item[0])[1]


def write_live_readfish_toml(
    *,
    base_toml_path: str | Path,
    targets_bed: str | Path,
    region_name: str = "robin_panel",
) -> Path:
    """Write ``{base_toml}_live`` with an updated region ``targets`` path."""
    base_path = Path(base_toml_path).expanduser()
    if not base_path.is_file():
        raise FileNotFoundError(f"Base readfish TOML not found: {base_path}")

    targets_path = Path(targets_bed).expanduser()
    if not targets_path.is_file():
        raise FileNotFoundError(f"Targets BED not found: {targets_path}")

    with base_path.open("rb") as handle:
        document = tomllib.load(handle)

    regions = document.get("regions")
    if not isinstance(regions, list) or not regions:
        raise ValueError(f"Base readfish TOML has no [[regions]] table: {base_path}")

    updated = False
    for region in regions:
        if not isinstance(region, dict):
            continue
        if region.get("name") == region_name:
            region["targets"] = str(targets_path)
            updated = True
            break

    if not updated:
        first = regions[0]
        if not isinstance(first, dict):
            raise ValueError(f"Invalid [[regions]] entry in {base_path}")
        first["targets"] = str(targets_path)

    destination = live_toml_path(base_path)
    temporary_path = Path(f"{destination}.tmp")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with temporary_path.open("wb") as handle:
        tomli_w.dump(document, handle)
    os.replace(temporary_path, destination)
    return destination


class ReadfishLiveRegistry:
    """Thread-safe registry of readfish runs awaiting live target updates."""

    _lock = threading.Lock()
    _sessions: dict[str, ReadfishLiveSession] = {}

    @classmethod
    def register(cls, session: ReadfishLiveSession) -> None:
        with cls._lock:
            cls._sessions[session.sample_id] = session
        LOGGER.info(
            "Registered readfish live updates for sample %s (base TOML: %s)",
            session.sample_id,
            session.base_toml_path,
        )

    @classmethod
    def unregister(cls, sample_id: str) -> None:
        with cls._lock:
            cls._sessions.pop(sample_id, None)

    @classmethod
    def get(cls, sample_id: str) -> Optional[ReadfishLiveSession]:
        with cls._lock:
            return cls._sessions.get(sample_id)

    @classmethod
    def notify_master_bed(
        cls,
        *,
        sample_id: str,
        master_bed_path: str | Path,
    ) -> Optional[Path]:
        """Write a live TOML when a new master BED is available for a registered sample."""
        master_path = Path(master_bed_path).expanduser()
        if not master_path.is_file():
            LOGGER.debug("Master BED path does not exist yet: %s", master_path)
            return None

        with cls._lock:
            session = cls._sessions.get(sample_id)
            if session is None:
                return None
            if session.last_master_bed_path == str(master_path):
                return None
            updated_session = ReadfishLiveSession(
                sample_id=session.sample_id,
                base_toml_path=session.base_toml_path,
                region_name=session.region_name,
                last_master_bed_path=str(master_path),
            )
            cls._sessions[sample_id] = updated_session

        try:
            live_path = write_live_readfish_toml(
                base_toml_path=updated_session.base_toml_path,
                targets_bed=master_path,
                region_name=updated_session.region_name,
            )
        except Exception:
            LOGGER.exception(
                "Failed to write readfish live TOML for sample %s", sample_id
            )
            return None

        LOGGER.info(
            "Updated readfish live targets for sample %s: %s -> %s",
            sample_id,
            master_path,
            live_path,
        )
        return live_path


def notify_readfish_live_targets(
    *,
    sample_id: str,
    master_bed_path: str | Path | None = None,
    work_dir: str | Path | None = None,
) -> Optional[Path]:
    """Notify the live updater that a master BED file is ready.

    When ``master_bed_path`` is omitted and ``work_dir`` is provided, the latest
    ``master_NNN.bed`` for the sample is resolved automatically.
    """
    if master_bed_path in (None, ""):
        if work_dir is None:
            return None
        resolved = find_latest_master_bed(work_dir=work_dir, sample_id=sample_id)
        if resolved is None:
            return None
        master_bed_path = resolved

    return ReadfishLiveRegistry.notify_master_bed(
        sample_id=sample_id,
        master_bed_path=master_bed_path,
    )
