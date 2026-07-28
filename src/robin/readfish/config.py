"""readfish settings loaded from workflow TOML."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Optional

DEFAULT_DORADO_ADDRESS = "ipc:///tmp/.guppy/5555"


@dataclass(frozen=True)
class ReadfishConfig:
    """Optional overrides for readfish experiment TOML generation and launch."""

    dorado_address: str = DEFAULT_DORADO_ADDRESS
    dorado_config: Optional[str] = None
    minimap2_index: Optional[str] = None
    log_dir: Optional[str] = None
    # Default True: PromethION / ROBIN runs require mapper_settings.mappy_rs.
    prom: bool = True
    mappy_rs_threads: int = 4
    min_chunks: int = 1
    max_chunks: int = 4
    validate_on_start: bool = True
    readfish_executable: str = "readfish"
    live_updates_enabled: bool = True
    live_region_name: str = "robin_panel"
    # Wait for MinKNOW acquisition before launching readfish.
    start_wait_timeout_seconds: float = 600.0
    start_wait_poll_seconds: float = 10.0

    @classmethod
    def from_mapping(cls, data: Optional[Mapping[str, Any]]) -> Optional[ReadfishConfig]:
        """Return a config when a ``[readfish]`` table is present."""
        if not isinstance(data, Mapping):
            return None
        if not data:
            return cls()

        dorado_address = _optional_str(data.get("dorado_address")) or DEFAULT_DORADO_ADDRESS
        log_dir = _optional_str(data.get("log_dir"))
        minimap2_index = _optional_str(data.get("minimap2_index"))
        dorado_config = _optional_str(data.get("dorado_config"))
        readfish_executable = _optional_str(data.get("readfish_executable")) or "readfish"
        live_region_name = _optional_str(data.get("live_region_name")) or "robin_panel"
        mappy_rs_threads = max(4, int(data.get("mappy_rs_threads", 4)))

        return cls(
            dorado_address=dorado_address,
            dorado_config=dorado_config,
            minimap2_index=minimap2_index,
            log_dir=log_dir,
            prom=bool(data.get("prom", True)),
            mappy_rs_threads=mappy_rs_threads,
            min_chunks=int(data.get("min_chunks", 1)),
            max_chunks=int(data.get("max_chunks", 4)),
            validate_on_start=bool(data.get("validate_on_start", True)),
            readfish_executable=readfish_executable,
            live_updates_enabled=bool(data.get("live_updates_enabled", True)),
            live_region_name=live_region_name,
            start_wait_timeout_seconds=float(
                data.get("start_wait_timeout_seconds", 600.0)
            ),
            start_wait_poll_seconds=float(data.get("start_wait_poll_seconds", 10.0)),
        )

    def resolve_log_file(self, *, sample_id: str, output_dir: Path) -> Path:
        if self.log_dir:
            base = Path(self.log_dir).expanduser()
        else:
            base = output_dir
        base.mkdir(parents=True, exist_ok=True)
        safe_sample = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in sample_id)
        return base / f"readfish_{safe_sample}.log"


def resolve_minimap2_index(
    *,
    alignment_reference: str,
    explicit_index: Optional[str] = None,
) -> str:
    """Choose the mapper reference path for readfish."""
    if explicit_index:
        return str(Path(explicit_index).expanduser())
    reference = Path(alignment_reference).expanduser()
    mmi = reference.with_suffix(".mmi")
    if mmi.is_file():
        return str(mmi)
    return str(reference)


DEFAULT_DORADO_VERSION = "v5.2.0"


def dorado_config_name(
    basecall_simplex_model: str,
    *,
    prefer_fast: bool = True,
) -> str:
    """Normalize a Dorado model name for readfish ``caller_settings.dorado.config``.

    Dorado servers expect the full model id including ``@version`` and a trailing
    ``||`` (e.g. ``dna_r10.4.1_e8.2_400bps_fast@v5.2.0||``). When deriving from a
    MinKNOW HAC/SUP simplex preset, prefer the matching fast model for adaptive
    decisions unless ``prefer_fast`` is false.
    """
    text = basecall_simplex_model.strip().rstrip("|").strip()
    if not text:
        raise ValueError("Dorado config model name is empty")

    if prefer_fast:
        for tier in ("_sup", "_hac"):
            if tier in text:
                text = text.replace(tier, "_fast", 1)
                break

    if "@" not in text:
        text = f"{text}@{DEFAULT_DORADO_VERSION}"

    return f"{text}||"


def _optional_str(value: Any) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    return text or None
