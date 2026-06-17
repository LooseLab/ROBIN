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
    docker_input: MNPFlexDockerInput
    username: Optional[str]
    password: Optional[str]
    base_url: str
    workflow_id: int
    client_id: str
    client_secret: str
    scope: str

    def is_enabled(self) -> bool:
        return self.backend in ("docker", "api")

    def describe_backend(self) -> str:
        if self.backend == "docker":
            return f"Docker ({self.docker_image or 'image not set'})"
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

    try:
        docker_timeout_s = int(os.getenv("MNPFLEX_DOCKER_TIMEOUT", "3600"))
    except ValueError:
        docker_timeout_s = 3600

    return MNPFlexConfig(
        backend=backend,
        docker_image=docker_image,
        docker_timeout_s=docker_timeout_s,
        docker_binary=(os.getenv("MNPFLEX_DOCKER_BINARY") or "docker").strip(),
        docker_extra_args=docker_extra_args,
        docker_input=docker_input,
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
