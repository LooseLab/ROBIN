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
    prom: bool = False
    mappy_rs_threads: int = 4
    min_chunks: int = 1
    max_chunks: int = 4
    validate_on_start: bool = True
    readfish_executable: str = "readfish"
    live_updates_enabled: bool = True
    live_region_name: str = "robin_panel"

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

        return cls(
            dorado_address=dorado_address,
            dorado_config=dorado_config,
            minimap2_index=minimap2_index,
            log_dir=log_dir,
            prom=bool(data.get("prom", False)),
            mappy_rs_threads=int(data.get("mappy_rs_threads", 4)),
            min_chunks=int(data.get("min_chunks", 1)),
            max_chunks=int(data.get("max_chunks", 4)),
            validate_on_start=bool(data.get("validate_on_start", True)),
            readfish_executable=readfish_executable,
            live_updates_enabled=bool(data.get("live_updates_enabled", True)),
            live_region_name=live_region_name,
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


def dorado_config_name(basecall_simplex_model: str) -> str:
    """Strip the Dorado version suffix from a MinKNOW simplex model name."""
    return basecall_simplex_model.split("@", 1)[0].strip()


def _optional_str(value: Any) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    return text or None
