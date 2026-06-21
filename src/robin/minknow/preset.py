"""ROBIN-compliant MinKNOW run presets."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

ROBIN_DEFAULT_BAM_READS_PER_FILE = 50_000
ROBIN_DEFAULT_BAM_BATCH_DURATION = 0
ROBIN_DEFAULT_EXPERIMENT_DURATION_HOURS = 24.0


@dataclass(frozen=True)
class RobinRunPreset:
    """MinKNOW protocol start settings matching ROBIN sequencing requirements."""

    kit: str
    basecall_simplex_model: str
    alignment_reference: str
    modified_models: tuple[str, ...] = ()
    bed_file: Optional[str] = None
    read_until_filter: Optional[str] = None
    read_until_reference: Optional[str] = None
    read_until_bed_file: Optional[str] = None
    bam_reads_per_file: int = ROBIN_DEFAULT_BAM_READS_PER_FILE
    bam_batch_duration: int = ROBIN_DEFAULT_BAM_BATCH_DURATION
    experiment_duration_hours: float = ROBIN_DEFAULT_EXPERIMENT_DURATION_HOURS
    experiment_group: str = "ROBIN_RUN"
    position: Optional[str] = None
    product_code: Optional[str] = None
    config_name: Optional[str] = None
    mux_scan_period: float = 1.5
    enable_basecalling: bool = True
    enable_bam: bool = True

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> RobinRunPreset:
        """Load from a ``[minknow.preset]`` table or flat preset mapping."""
        nested = data.get("preset")
        if isinstance(nested, Mapping):
            source: Mapping[str, Any] = nested
        else:
            source = data

        kit = _required_str(source, "kit")
        basecall_simplex_model = _required_str(
            source,
            "basecall_simplex_model",
            aliases=("basecall_simplex", "simplex_model"),
        )
        alignment_reference = (
            _optional_str(source.get("alignment_reference"))
            or _optional_str(source.get("reference"))
            or _optional_str(source.get("alignment_ref"))
            or ""
        )

        modified_models = _optional_str_list(
            source.get("modified_models", source.get("basecall_modified_models"))
        )
        bed_file = _optional_str(source.get("bed_file"))
        read_until_filter = _optional_str(source.get("read_until_filter"))
        read_until_reference = _optional_str(
            source.get("read_until_reference", source.get("read_until_ref"))
        )
        read_until_bed_file = _optional_str(source.get("read_until_bed_file"))

        bam_reads = source.get("bam_reads_per_file", ROBIN_DEFAULT_BAM_READS_PER_FILE)
        bam_batch = source.get("bam_batch_duration", ROBIN_DEFAULT_BAM_BATCH_DURATION)
        duration = source.get(
            "experiment_duration_hours",
            source.get("experiment_duration", ROBIN_DEFAULT_EXPERIMENT_DURATION_HOURS),
        )

        return cls(
            kit=kit,
            basecall_simplex_model=basecall_simplex_model,
            alignment_reference=alignment_reference,
            modified_models=modified_models,
            bed_file=bed_file,
            read_until_filter=read_until_filter,
            read_until_reference=read_until_reference,
            read_until_bed_file=read_until_bed_file,
            bam_reads_per_file=int(bam_reads),
            bam_batch_duration=int(bam_batch),
            experiment_duration_hours=float(duration),
            experiment_group=_optional_str(source.get("experiment_group"))
            or _optional_str(source.get("protocol_group_id"))
            or "ROBIN_RUN",
            position=_optional_str(source.get("position")),
            product_code=_optional_str(source.get("product_code")),
            config_name=_optional_str(source.get("config_name")),
            mux_scan_period=float(source.get("mux_scan_period", 1.5)),
            enable_basecalling=bool(source.get("enable_basecalling", True)),
            enable_bam=bool(source.get("enable_bam", True)),
        )

    def with_overrides(self, **kwargs: Any) -> RobinRunPreset:
        allowed = {f.name for f in self.__dataclass_fields__.values()}
        filtered = {key: value for key, value in kwargs.items() if value is not None}
        unknown = set(filtered) - allowed
        if unknown:
            raise ValueError(f"Unknown preset override(s): {', '.join(sorted(unknown))}")
        return replace(self, **filtered)

    def validate(self, *, check_paths: bool = False) -> list[str]:
        """Return a list of validation errors (empty if valid).

        ``check_paths`` verifies reference/BED paths on the **local** filesystem
        (the machine running ROBIN). Preset paths normally live on the MinKNOW
        host, so path checking is opt-in.
        """
        errors: list[str] = []

        if not self.kit:
            errors.append("kit is required")
        if self.enable_basecalling and not self.basecall_simplex_model:
            errors.append("basecall_simplex_model is required when basecalling is enabled")
        if self.enable_basecalling and not self.alignment_reference:
            errors.append("alignment_reference is required when basecalling is enabled")
        if self.bed_file and not self.alignment_reference:
            errors.append("bed_file requires alignment_reference")
        if self.read_until_filter and not self.effective_read_until_reference():
            errors.append("read_until_filter requires read_until_reference or alignment_reference")
        if self.read_until_bed_file and not self.effective_read_until_reference():
            errors.append("read_until_bed_file requires read_until_reference or alignment_reference")
        if self.read_until_filter and self.read_until_filter not in {"enrich", "deplete"}:
            errors.append("read_until_filter must be 'enrich' or 'deplete'")
        if self.bam_reads_per_file <= 0:
            errors.append("bam_reads_per_file must be positive")
        if self.bam_batch_duration < 0:
            errors.append("bam_batch_duration must be >= 0")
        if self.experiment_duration_hours <= 0:
            errors.append("experiment_duration_hours must be positive")
        if not self.enable_bam:
            errors.append("ROBIN requires BAM output (enable_bam must be true)")

        if check_paths:
            errors.extend(self._path_errors())

        return errors

    def _path_errors(self) -> list[str]:
        errors: list[str] = []
        for label, path in (
            ("alignment_reference", self.alignment_reference),
            ("bed_file", self.bed_file),
            ("read_until_reference", self.effective_read_until_reference()),
            ("read_until_bed_file", self.effective_read_until_bed_file()),
        ):
            if not path:
                continue
            if not Path(path).expanduser().exists():
                errors.append(
                    f"{label} not found on this machine at {path} "
                    "(preset paths must exist on the MinKNOW host; "
                    "local check is optional)"
                )
        return errors

    def effective_read_until_reference(self) -> Optional[str]:
        return self.read_until_reference or self.alignment_reference

    def effective_read_until_bed_file(self) -> Optional[str]:
        return self.read_until_bed_file or self.bed_file

    def adaptive_sampling_enabled(self) -> bool:
        return bool(self.read_until_filter and self.effective_read_until_reference())

    def summary_lines(self) -> list[str]:
        lines = [
            f"Kit: {self.kit}",
            f"Simplex model: {self.basecall_simplex_model}",
        ]
        if self.modified_models:
            lines.append(f"Modified models: {', '.join(self.modified_models)}")
        lines.append(f"Alignment reference: {self.alignment_reference}")
        if self.bed_file:
            lines.append(f"BED file: {self.bed_file}")
        lines.append(
            f"BAM rollover: {self.bam_reads_per_file} reads, "
            f"batch duration {self.bam_batch_duration}s"
        )
        lines.append(f"Duration: {self.experiment_duration_hours} h")
        if self.adaptive_sampling_enabled():
            lines.append(
                f"Adaptive sampling: {self.read_until_filter} "
                f"({self.effective_read_until_reference()})"
            )
        return lines


def _required_str(
    data: Mapping[str, Any],
    key: str,
    *,
    aliases: Sequence[str] = (),
) -> str:
    for candidate in (key, *aliases):
        value = _optional_str(data.get(candidate))
        if value:
            return value
    raise ValueError(f"Preset missing required key: {key}")


def _optional_str(value: Any) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _optional_str_list(value: Any) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        text = value.strip()
        return (text,) if text else ()
    if isinstance(value, Sequence):
        return tuple(str(item).strip() for item in value if str(item).strip())
    return (str(value).strip(),)
