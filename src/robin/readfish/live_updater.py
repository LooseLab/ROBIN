"""Push live readfish target updates when master BED files change."""

from __future__ import annotations

import json
import logging
import os
import re
import threading
import tomllib
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

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


SAMPLE_READFISH_TOML_NAME = "readfish.toml"


def expected_readfish_toml_name(sample_id: str) -> str:
    """Return the legacy experiment TOML filename for a sample (cwd layout)."""
    safe = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in sample_id)
    return f"readfish_{safe or 'run'}.toml"


def resolve_sample_readfish_dir(
    *,
    sample_id: str,
    work_directory: str | Path | None = None,
    output_dir: str | Path | None = None,
) -> Path:
    """Return the directory where readfish TOML/logs should live.

    Preference order:
    1. Explicit ``output_dir``
    2. ``{work_directory}/{sample_id}`` (ROBIN sample output folder)
    3. Current working directory (legacy fallback)
    """
    if output_dir is not None and str(output_dir).strip():
        return Path(output_dir).expanduser()
    if work_directory is not None and str(work_directory).strip():
        return Path(work_directory).expanduser() / sample_id
    return Path.cwd()


def resolve_readfish_toml_path(
    *,
    sample_id: str,
    work_directory: str | Path | None = None,
    output_dir: str | Path | None = None,
) -> Path:
    """Return the path for the base readfish experiment TOML."""
    directory = resolve_sample_readfish_dir(
        sample_id=sample_id,
        work_directory=work_directory,
        output_dir=output_dir,
    )
    using_sample_folder = bool(
        work_directory is not None
        and str(work_directory).strip()
        and (output_dir is None or not str(output_dir).strip())
    ) or (
        output_dir is not None
        and str(output_dir).strip()
        and directory.name == sample_id
    )
    if using_sample_folder:
        return directory / SAMPLE_READFISH_TOML_NAME
    return directory / expected_readfish_toml_name(sample_id)


def discover_readfish_toml(
    sample_id: str,
    *,
    work_dir: str | Path | None = None,
    master_bed_path: str | Path | None = None,
    region_name: str = "robin_panel",
) -> Optional[Path]:
    """Locate the base readfish TOML near the run without a live session.

    Prefers ``{work_dir}/{sample}/readfish.toml``, then legacy
    ``readfish_<sample>.toml`` under cwd / parents of the master BED.
    """
    del region_name  # reserved for future region-aware discovery
    legacy_name = expected_readfish_toml_name(sample_id)
    names = (SAMPLE_READFISH_TOML_NAME, legacy_name)
    roots: list[Path] = []

    env_dir = os.environ.get("ROBIN_READFISH_TOML_DIR")
    if env_dir:
        roots.append(Path(env_dir).expanduser())

    roots.append(Path.cwd())

    if work_dir is not None:
        work = Path(work_dir).expanduser()
        roots.extend([work / sample_id, work, work.parent])

    if master_bed_path is not None:
        master = Path(master_bed_path).expanduser()
        roots.extend(
            [
                master.parent.parent,  # sample dir
                master.parent,  # bed_files
                master.parent.parent.parent,  # work_dir
                master.parent.parent.parent.parent,
            ]
        )

    seen: set[str] = set()
    for root in roots:
        try:
            resolved_root = root.resolve()
        except OSError:
            resolved_root = root
        key = str(resolved_root)
        if key in seen:
            continue
        seen.add(key)
        for name in names:
            candidate = resolved_root / name
            if candidate.is_file():
                _announce(
                    f"Discovered base TOML for sample {sample_id!r} via fallback: "
                    f"{candidate}"
                )
                return candidate
    return None


def live_session_dir() -> Path:
    """Directory used to persist live-update sessions across processes."""
    override = os.environ.get("ROBIN_READFISH_LIVE_DIR")
    if override:
        return Path(override).expanduser()
    return Path.home() / ".robin" / "readfish_live"


def live_session_path(sample_id: str) -> Path:
    """Return the on-disk session path for ``sample_id``."""
    safe = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in sample_id)
    return live_session_dir() / f"{safe}.json"


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

    stamp_path = Path(f"{destination}.stamp")
    stamp_path.write_text(
        f"master_bed={targets_path}\n"
        f"base_toml={base_path}\n"
        f"region={region_name}\n",
        encoding="utf-8",
    )
    return destination


def _announce(message: str) -> None:
    text = f"[readfish] {message}"
    print(text, flush=True)
    LOGGER.info("%s", text)


class ReadfishLiveRegistry:
    """Thread-safe registry of readfish runs awaiting live target updates.

    Sessions are kept in memory and mirrored to disk so Ray workers / other
    processes can still write ``*_live`` when a master BED appears.
    """

    _lock = threading.Lock()
    _sessions: dict[str, ReadfishLiveSession] = {}

    @classmethod
    def register(cls, session: ReadfishLiveSession) -> None:
        with cls._lock:
            cls._sessions[session.sample_id] = session
            _write_session_file(session)
        _announce(
            f"Registered live updates for sample {session.sample_id!r} "
            f"(base TOML: {session.base_toml_path}; "
            f"session file: {live_session_path(session.sample_id)})"
        )

    @classmethod
    def unregister(cls, sample_id: str) -> None:
        with cls._lock:
            cls._sessions.pop(sample_id, None)
            path = live_session_path(sample_id)
            try:
                path.unlink(missing_ok=True)
            except OSError:
                LOGGER.debug(
                    "Could not remove live session file %s", path, exc_info=True
                )

    @classmethod
    def get(cls, sample_id: str) -> Optional[ReadfishLiveSession]:
        with cls._lock:
            session = cls._sessions.get(sample_id)
            if session is not None:
                return session
            loaded = _read_session_file(sample_id)
            if loaded is not None:
                cls._sessions[sample_id] = loaded
            return loaded

    @classmethod
    def list_sessions(cls) -> list[ReadfishLiveSession]:
        """Return in-memory sessions plus any persisted session files."""
        with cls._lock:
            by_id = dict(cls._sessions)
            session_dir = live_session_dir()
            if session_dir.is_dir():
                for path in session_dir.glob("*.json"):
                    try:
                        data = json.loads(path.read_text(encoding="utf-8"))
                        session = ReadfishLiveSession(
                            sample_id=str(data["sample_id"]),
                            base_toml_path=str(data["base_toml_path"]),
                            region_name=str(data.get("region_name") or "robin_panel"),
                            last_master_bed_path=data.get("last_master_bed_path"),
                        )
                    except (
                        OSError,
                        KeyError,
                        TypeError,
                        ValueError,
                        json.JSONDecodeError,
                    ):
                        continue
                    by_id.setdefault(session.sample_id, session)
            return sorted(by_id.values(), key=lambda item: item.sample_id)

    @classmethod
    def notify_master_bed(
        cls,
        *,
        sample_id: str,
        master_bed_path: str | Path,
        work_dir: str | Path | None = None,
        region_name: str = "robin_panel",
    ) -> Optional[Path]:
        """Write a live TOML when a new master BED is available for a sample."""
        master_path = Path(master_bed_path).expanduser()
        if not master_path.is_file():
            _announce(
                f"Live update skipped for sample {sample_id!r}: "
                f"master BED not found ({master_path})"
            )
            return None

        with cls._lock:
            session = cls._sessions.get(sample_id) or _read_session_file(sample_id)
            if session is None:
                discovered = discover_readfish_toml(
                    sample_id,
                    work_dir=work_dir,
                    master_bed_path=master_path,
                    region_name=region_name,
                )
                if discovered is None:
                    _announce(
                        f"Live update skipped for sample {sample_id!r}: "
                        "no registered readfish session and could not find "
                        f"{SAMPLE_READFISH_TOML_NAME} or "
                        f"{expected_readfish_toml_name(sample_id)} "
                        f"(looked for {live_session_path(sample_id)}; "
                        f"cwd={Path.cwd()})"
                    )
                    return None
                session = ReadfishLiveSession(
                    sample_id=sample_id,
                    base_toml_path=str(discovered),
                    region_name=region_name,
                )
                cls._sessions[sample_id] = session
                _write_session_file(session)
                _announce(
                    f"Auto-registered live session for sample {sample_id!r} "
                    f"from discovered TOML {discovered}"
                )

            if session.last_master_bed_path == str(master_path):
                _announce(
                    f"Live update skipped for sample {sample_id!r}: "
                    f"already applied {master_path}"
                )
                return None
            updated_session = ReadfishLiveSession(
                sample_id=session.sample_id,
                base_toml_path=session.base_toml_path,
                region_name=session.region_name,
                last_master_bed_path=str(master_path),
            )
            cls._sessions[sample_id] = updated_session
            _write_session_file(updated_session)

        try:
            live_path = write_live_readfish_toml(
                base_toml_path=updated_session.base_toml_path,
                targets_bed=master_path,
                region_name=updated_session.region_name,
            )
        except Exception as exc:
            LOGGER.exception(
                "Failed to write readfish live TOML for sample %s", sample_id
            )
            _announce(f"Live update FAILED for sample {sample_id!r}: {exc}")
            return None

        _announce(
            f"Wrote live TOML for sample {sample_id!r}: {live_path} "
            f"(targets={master_path})"
        )
        _announce(f"Verify stamp: {live_path}.stamp")
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

    If no live session was registered (common when master BED generation runs in
    a Ray worker), falls back to discovering ``readfish_<sample>.toml`` near the
    cwd / work_dir / master BED path.
    """
    if master_bed_path in (None, ""):
        if work_dir is None:
            return None
        resolved = find_latest_master_bed(work_dir=work_dir, sample_id=sample_id)
        if resolved is None:
            _announce(
                f"Live update skipped for sample {sample_id!r}: "
                f"no master_NNN.bed under {work_dir}/{sample_id}/bed_files"
            )
            return None
        master_bed_path = resolved

    return ReadfishLiveRegistry.notify_master_bed(
        sample_id=sample_id,
        master_bed_path=master_bed_path,
        work_dir=work_dir,
    )


def _write_session_file(session: ReadfishLiveSession) -> None:
    path = live_session_path(session.sample_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(f"{path}.tmp")
    temporary.write_text(json.dumps(asdict(session), indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def _read_session_file(sample_id: str) -> Optional[ReadfishLiveSession]:
    path = live_session_path(sample_id)
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return ReadfishLiveSession(
            sample_id=str(data["sample_id"]),
            base_toml_path=str(data["base_toml_path"]),
            region_name=str(data.get("region_name") or "robin_panel"),
            last_master_bed_path=data.get("last_master_bed_path"),
        )
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
        LOGGER.debug("Could not read live session file %s", path, exc_info=True)
        return None
