"""MNP-Flex integration configuration (API vs local Docker)."""

from __future__ import annotations

import os
import shlex
from dataclasses import dataclass
from typing import List, Literal, Optional

MNPFlexBackend = Literal["docker", "api", "disabled"]
MNPFlexDockerInput = Literal["full", "subset"]


@dataclass(frozen=True)
class MNPFlexConfig:
    backend: MNPFlexBackend
    docker_image: Optional[str]
    docker_timeout_s: int
    docker_binary: str
    docker_extra_args: List[str]
    docker_cmd_args: List[str]
    docker_input: MNPFlexDockerInput
    docker_technology: Optional[str]
    docker_impute: Optional[str]
    docker_json: Optional[str]
    docker_mode: Optional[str]
    docker_sample_type: Optional[str]
    docker_threshold: Optional[str]
    docker_columns: Optional[str]
    username: Optional[str]
    password: Optional[str]
    base_url: str
    workflow_id: int
    client_id: str
    client_secret: str
    scope: str

    def is_enabled(self) -> bool:
        return self.backend in ("docker", "api")

    def docker_cli_options(self) -> List[tuple[str, str]]:
        """Container flags after the image name (1.1+ / preview CLI)."""
        pairs = (
            ("--technology", self.docker_technology),
            ("--columns", self.docker_columns),
            ("--mode", self.docker_mode),
            ("--sample_type", self.docker_sample_type),
            ("--threshold", self.docker_threshold),
            ("--impute", self.docker_impute),
            ("--json", self.docker_json),
        )
        return [(flag, value) for flag, value in pairs if value]

    def docker_options_payload(self) -> dict:
        return {
            "technology": self.docker_technology,
            "columns": self.docker_columns,
            "mode": self.docker_mode,
            "sample_type": self.docker_sample_type,
            "threshold": self.docker_threshold,
            "impute": self.docker_impute,
            "json": self.docker_json,
        }

    def describe_backend(self) -> str:
        if self.backend == "docker":
            parts = [self.docker_image or "image not set"]
            if self.docker_mode:
                mode = self.docker_mode
                if self.docker_mode == "adaptive" and self.docker_threshold:
                    mode = f"{self.docker_mode}/{self.docker_threshold}"
                parts.append(mode)
            if self.docker_impute:
                parts.append(f"impute={self.docker_impute}")
            return f"Docker ({', '.join(parts)})"
        if self.backend == "api":
            return f"Epignostix API ({self.base_url})"
        return "disabled"

    def validation_error(self) -> Optional[str]:
        if self.backend == "disabled":
            return None
        if self.backend == "docker":
            if not (self.docker_image or "").strip():
                return (
                    "MNPFLEX_BACKEND=docker requires MNPFLEX_DOCKER_IMAGE "
                    "(full image name and tag, e.g. mnpflex-synnovis:1.0.0)."
                )
            return None
        if self.backend == "api":
            if not self.username or not self.password:
                return (
                    "MNPFLEX_BACKEND=api requires MNPFLEX_USERNAME and "
                    "MNPFLEX_PASSWORD."
                )
            return None
        return f"Unsupported MNPFLEX_BACKEND value: {self.backend!r}"


def _resolve_backend(raw: str) -> MNPFlexBackend:
    value = (raw or "").strip().lower()
    if value in ("docker", "api", "disabled"):
        return value  # type: ignore[return-value]
    if value in ("", "none", "off"):
        return "disabled"
    raise ValueError(
        f"Invalid MNPFLEX_BACKEND={raw!r}. Expected docker, api, or disabled."
    )


def _optional_env(name: str) -> Optional[str]:
    raw = os.getenv(name)
    if raw is None:
        return None
    value = str(raw).strip()
    return value or None


def _uses_extended_docker_cli(image: Optional[str]) -> bool:
    """Official 1.1+/preview images accept --technology, --impute, --mode, etc."""
    name = (image or "").strip().lower()
    if not name or name.startswith("mnpflex-synnovis"):
        return False
    return name.startswith("mnpflex:") or "preview" in name


def _default_docker_technology(image: Optional[str]) -> Optional[str]:
    """ROBIN feeds 18-column bedMethyl; official 1.1+ images need --technology nanopore."""
    if _uses_extended_docker_cli(image):
        return "nanopore"
    return None


def _env_or_default(name: str, default: Optional[str]) -> Optional[str]:
    value = _optional_env(name)
    if value is None:
        return default
    if value.lower() in {"off", "none"}:
        return None
    return value


def _infer_backend_when_unset(
    *,
    docker_image: Optional[str],
    username: Optional[str],
    password: Optional[str],
) -> MNPFlexBackend:
    """Preserve legacy behaviour: credentials alone enable the API backend."""
    if docker_image and not (username and password):
        return "docker"
    if username and password:
        return "api"
    return "disabled"


def load_mnpflex_config() -> MNPFlexConfig:
    """Load MNP-Flex settings from the server environment."""
    docker_image = (os.getenv("MNPFLEX_DOCKER_IMAGE") or "").strip() or None
    username = os.getenv("MNPFLEX_USERNAME") or os.getenv("EPIGNOSTIX_USERNAME")
    password = os.getenv("MNPFLEX_PASSWORD") or os.getenv("EPIGNOSTIX_PASSWORD")

    backend_raw = os.getenv("MNPFLEX_BACKEND")
    if backend_raw is None or not str(backend_raw).strip():
        backend = _infer_backend_when_unset(
            docker_image=docker_image,
            username=username,
            password=password,
        )
    else:
        backend = _resolve_backend(str(backend_raw))

    workflow_id_env = os.getenv("MNPFLEX_WORKFLOW_ID", "18")
    try:
        workflow_id = int(workflow_id_env)
    except ValueError:
        workflow_id = 18

    docker_input_raw = (os.getenv("MNPFLEX_DOCKER_INPUT") or "full").strip().lower()
    if docker_input_raw not in ("full", "subset"):
        docker_input_raw = "full"
    docker_input: MNPFlexDockerInput = docker_input_raw  # type: ignore[assignment]

    extra_raw = os.getenv("MNPFLEX_DOCKER_EXTRA_ARGS", "")
    docker_extra_args = shlex.split(extra_raw) if extra_raw.strip() else []
    cmd_raw = os.getenv("MNPFLEX_DOCKER_CMD_ARGS", "")
    docker_cmd_args = shlex.split(cmd_raw) if cmd_raw.strip() else []

    try:
        docker_timeout_s = int(os.getenv("MNPFLEX_DOCKER_TIMEOUT", "3600"))
    except ValueError:
        docker_timeout_s = 3600

    extended = _uses_extended_docker_cli(docker_image)
    if extended:
        technology = _env_or_default(
            "MNPFLEX_DOCKER_TECHNOLOGY",
            _default_docker_technology(docker_image),
        )
        impute = _env_or_default("MNPFLEX_DOCKER_IMPUTE", "true")
        json_report = _env_or_default("MNPFLEX_DOCKER_JSON", "true")
        mode = _env_or_default("MNPFLEX_DOCKER_MODE", "adaptive")
        sample_type = _env_or_default("MNPFLEX_DOCKER_SAMPLE_TYPE", None)
        threshold = _env_or_default(
            "MNPFLEX_DOCKER_THRESHOLD",
            "lookup" if mode == "adaptive" else None,
        )
        columns = _optional_env("MNPFLEX_DOCKER_COLUMNS")
    else:
        technology = None
        impute = None
        json_report = None
        mode = None
        sample_type = None
        threshold = None
        columns = None

    return MNPFlexConfig(
        backend=backend,
        docker_image=docker_image,
        docker_timeout_s=docker_timeout_s,
        docker_binary=(os.getenv("MNPFLEX_DOCKER_BINARY") or "docker").strip(),
        docker_extra_args=docker_extra_args,
        docker_cmd_args=docker_cmd_args,
        docker_input=docker_input,
        docker_technology=technology,
        docker_impute=impute,
        docker_json=json_report,
        docker_mode=mode,
        docker_sample_type=sample_type,
        docker_threshold=threshold,
        docker_columns=columns,
        username=username,
        password=password,
        base_url=os.getenv("MNPFLEX_BASE_URL", "https://app.epignostix.com"),
        workflow_id=workflow_id,
        client_id=os.getenv("MNPFLEX_CLIENT_ID", "ROBIN"),
        client_secret=os.getenv("MNPFLEX_CLIENT_SECRET", "SECRET"),
        scope=os.getenv("MNPFLEX_SCOPE", ""),
    )


def is_mnpflex_enabled() -> bool:
    config = load_mnpflex_config()
    return config.is_enabled() and config.validation_error() is None
