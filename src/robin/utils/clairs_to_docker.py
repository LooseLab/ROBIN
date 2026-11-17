"""Prepare a ClairS-TO image whose bundled resources are readable by any UID.

Upstream ``hkubal/clairs-to`` images ship CNA/Verdict loci files (and similar
resources) with restrictive archive permissions such as ``0640`` owned by a
non-root UID. ROBIN runs ClairS-TO as the host ``uid:gid`` so bind-mounted
outputs are not root-owned. That combination makes the loci files exist but
not be readable, and Verdict then reports them as missing.

ROBIN keeps ``--user <host_uid>:<host_gid>``. When the upstream image is not
already world-readable, we derive a local image with ``chmod -R a+rX`` on
the bundled resource trees (equivalent to a tiny ``FROM hkubal/clairs-to``
Dockerfile) and run that instead.
"""

from __future__ import annotations

import logging
from typing import Any, Mapping, Optional, Sequence

logger = logging.getLogger(__name__)

CLAIRS_TO_READABLE_IMAGE = "robin/clairs-to:readable"
SOURCE_IMAGE_LABEL = "robin.clairs_to.source_image"
SOURCE_ID_LABEL = "robin.clairs_to.source_id"
# Layout 1 committed the chmod container's ENTRYPOINT ["/bin/sh"], so
# `/opt/bin/run_clairs_to` was executed as a shell script. Layout 2 restores
# the upstream Entrypoint/Cmd on commit.
LAYOUT_LABEL = "robin.clairs_to.layout"
LAYOUT_VERSION = "2"

CLAIRS_TO_RESOURCE_PATHS: tuple[str, ...] = (
    "/opt/micromamba/envs/clairs-to/bin/clairs-to_models",
    "/opt/micromamba/envs/clairs-to/bin/clairs-to_databases",
    "/opt/micromamba/envs/clairs-to/bin/clairs-to_cna_data",
    "/opt/bin",
)

CNA_LOCI_PROBE = (
    "/opt/micromamba/envs/clairs-to/bin/clairs-to_cna_data/"
    "reference_files/loci_files/G1000_loci_hg38_chr1.txt"
)


class ClairsToImageError(RuntimeError):
    """ClairS-TO image is missing, could not be patched, or is unreadable."""


def _image_labels(image: Any) -> dict[str, str]:
    labels: Mapping[str, Any] | None = getattr(image, "labels", None)
    if labels:
        return {str(k): str(v) for k, v in labels.items() if v is not None}
    attrs = getattr(image, "attrs", None) or {}
    config = attrs.get("Config") or attrs.get("ContainerConfig") or {}
    raw = config.get("Labels") or {}
    return {str(k): str(v) for k, v in raw.items() if v is not None}


def _image_id(image: Any) -> str:
    image_id = getattr(image, "id", None) or (getattr(image, "attrs", {}) or {}).get("Id")
    if not image_id:
        raise ClairsToImageError("Could not determine Docker image ID")
    return str(image_id)


def _get_image(client: Any, name: str) -> Any:
    try:
        return client.images.get(name)
    except Exception:
        return None


def _image_config(image: Any) -> dict[str, Any]:
    attrs = getattr(image, "attrs", None) or {}
    config = attrs.get("Config") or attrs.get("ContainerConfig") or {}
    return dict(config) if isinstance(config, dict) else {}


def _commit_config_from_source(
    source_image_obj: Any,
    *,
    source_image: str,
    source_id: str,
) -> dict[str, Any]:
    """Keep upstream Entrypoint/Cmd; docker commit otherwise inherits /bin/sh."""
    config = _image_config(source_image_obj)
    labels = {
        str(k): str(v)
        for k, v in (config.get("Labels") or {}).items()
        if v is not None
    }
    labels[SOURCE_IMAGE_LABEL] = source_image
    labels[SOURCE_ID_LABEL] = source_id
    labels[LAYOUT_LABEL] = LAYOUT_VERSION
    # An omitted/null Entrypoint on commit keeps the chmod container's
    # ENTRYPOINT ["/bin/sh"]. Empty list restores "no entrypoint".
    entrypoint = config.get("Entrypoint") or []
    return {
        "Entrypoint": entrypoint,
        "Cmd": config.get("Cmd"),
        "WorkingDir": config.get("WorkingDir") or "",
        "Env": list(config.get("Env") or []),
        "User": config.get("User") or "",
        "Labels": labels,
    }


def _run_in_image(
    client: Any,
    image: str,
    script: str,
    *,
    user: Optional[str] = None,
) -> tuple[int, str]:
    """Run ``script`` with ``/bin/sh -c``, overriding any image ENTRYPOINT."""
    container: dict[str, Any] | None = None
    try:
        kwargs: dict[str, Any] = {
            "image": image,
            "command": ["-c", script],
            "entrypoint": ["/bin/sh"],
        }
        if user:
            kwargs["user"] = user
        container = client.api.create_container(**kwargs)
        container_id = container.get("Id")
        if not container_id:
            raise ClairsToImageError("Docker did not return a container ID")
        client.api.start(container=container_id)
        result = client.api.wait(container=container_id)
        status = int(result.get("StatusCode", 1))
        logs = client.api.logs(container=container_id, stdout=True, stderr=True)
        text = logs.decode("utf-8", errors="replace") if isinstance(logs, (bytes, bytearray)) else str(logs or "")
        return status, text.strip()
    except ClairsToImageError:
        raise
    except Exception as exc:
        raise ClairsToImageError(
            f"Failed to run command in ClairS-TO image {image!r}: {exc}"
        ) from exc
    finally:
        if container and container.get("Id"):
            try:
                client.api.remove_container(container=container["Id"], force=True)
            except Exception:
                logger.debug("Could not remove temporary ClairS-TO container", exc_info=True)


def _chmod_script(paths: Sequence[str] = CLAIRS_TO_RESOURCE_PATHS) -> str:
    quoted = " ".join(f"'{path}'" for path in paths)
    return (
        "set -eu; "
        f"for d in {quoted}; do "
        'if [ -e "$d" ]; then chmod -R a+rX "$d"; fi; '
        "done"
    )


def _probe_script(path: str = CNA_LOCI_PROBE) -> str:
    # Missing resources (older/layout-changed images) are not a hard failure.
    # Existing but unreadable files are: that is the Verdict "does not appear
    # to exist" failure mode.
    return (
        f'f="{path}"; '
        'if [ -e "$f" ]; then '
        'if [ -r "$f" ]; then echo READABLE; exit 0; fi; '
        'echo UNREADABLE; exit 1; '
        "fi; "
        "echo MISSING; exit 0"
    )


def probe_clairs_to_resource_access(
    client: Any,
    image: str,
    *,
    user: Optional[str],
) -> tuple[str, str]:
    """Return ``(status, logs)`` where status is readable/unreadable/missing/error."""
    try:
        code, logs = _run_in_image(client, image, _probe_script(), user=user)
    except ClairsToImageError as exc:
        return "error", str(exc)
    if code == 0:
        token = "missing" if "MISSING" in logs.upper() else "readable"
        return token, logs
    if "UNREADABLE" in logs.upper() or code == 1:
        return "unreadable", logs
    return "error", logs or f"probe exited with status {code}"


def derived_image_matches_source(
    client: Any,
    source_image: str,
    source_id: str,
    *,
    derived_image: str = CLAIRS_TO_READABLE_IMAGE,
) -> bool:
    image = _get_image(client, derived_image)
    if image is None:
        return False
    labels = _image_labels(image)
    return (
        labels.get(SOURCE_ID_LABEL) == source_id
        and labels.get(SOURCE_IMAGE_LABEL) == source_image
        and labels.get(LAYOUT_LABEL) == LAYOUT_VERSION
    )


def _commit_readable_image(
    client: Any,
    source_image: str,
    source_id: str,
    source_image_obj: Any,
    *,
    derived_image: str = CLAIRS_TO_READABLE_IMAGE,
) -> str:
    logger.info(
        "Deriving host-user-readable ClairS-TO image %s from %s "
        "(chmod -R a+rX on bundled models/databases/CNA resources)",
        derived_image,
        source_image,
    )
    container: dict[str, Any] | None = None
    try:
        container = client.api.create_container(
            image=source_image,
            command=["-c", _chmod_script()],
            entrypoint=["/bin/sh"],
            user="0:0",
        )
        container_id = container.get("Id")
        if not container_id:
            raise ClairsToImageError("Docker did not return a container ID for chmod")
        client.api.start(container=container_id)
        result = client.api.wait(container=container_id)
        status = int(result.get("StatusCode", 1))
        if status != 0:
            logs = client.api.logs(container=container_id, stdout=True, stderr=True)
            text = (
                logs.decode("utf-8", errors="replace")
                if isinstance(logs, (bytes, bytearray))
                else str(logs or "")
            )
            raise ClairsToImageError(
                f"Failed to chmod ClairS-TO resources in {source_image}: "
                f"exit {status}: {text.strip()}"
            )
        repo, tag = derived_image.split(":", 1)
        commit_conf = _commit_config_from_source(
            source_image_obj, source_image=source_image, source_id=source_id
        )
        client.api.commit(
            container=container_id,
            repository=repo,
            tag=tag,
            conf=commit_conf,
        )
        logger.info("Derived ClairS-TO image ready: %s", derived_image)
        return derived_image
    except ClairsToImageError:
        raise
    except Exception as exc:
        raise ClairsToImageError(
            f"Failed to derive readable ClairS-TO image from {source_image}: {exc}"
        ) from exc
    finally:
        if container and container.get("Id"):
            try:
                client.api.remove_container(container=container["Id"], force=True)
            except Exception:
                logger.debug("Could not remove ClairS-TO chmod container", exc_info=True)


def ensure_readable_clairs_to_image(
    client: Any,
    source_image: str,
    *,
    user: Optional[str],
    derived_image: str = CLAIRS_TO_READABLE_IMAGE,
) -> str:
    """Return an image name that is readable by ``user``.

    If ``user`` is None (root), the upstream image is used. Otherwise a local
    derived image is reused when it still matches the source image ID, or built
    by chmod'ing bundled resources and committing the result.
    """
    if not user:
        logger.info(
            "ClairS-TO will run as root; using upstream image %s without deriving",
            source_image,
        )
        return source_image

    source = _get_image(client, source_image)
    if source is None:
        raise ClairsToImageError(
            f"ClairS-TO image {source_image!r} is not available locally"
        )
    source_id = _image_id(source)

    if derived_image_matches_source(
        client, source_image, source_id, derived_image=derived_image
    ):
        logger.info(
            "Reusing derived ClairS-TO image %s for source %s (%s)",
            derived_image,
            source_image,
            source_id[:19],
        )
        return derived_image

    logger.warning(
        "Deriving a local ClairS-TO image with world-readable bundled resources "
        "so SNP calling can run as host user %s (upstream=%s).",
        user,
        source_image,
    )
    return _commit_readable_image(
        client,
        source_image,
        source_id,
        source,
        derived_image=derived_image,
    )


def assert_clairs_to_resources_readable(
    client: Any,
    image: str,
    *,
    user: Optional[str],
) -> None:
    """Raise if CNA loci (when present) are not readable by ``user``."""
    if not user:
        return
    status, logs = probe_clairs_to_resource_access(client, image, user=user)
    if status in {"readable", "missing"}:
        logger.info(
            "ClairS-TO resource preflight for %s as user %s: %s",
            image,
            user,
            status,
        )
        return
    raise ClairsToImageError(
        f"ClairS-TO image {image!r} contains bundled resources that are not "
        f"readable by runtime user {user}. Verdict's alleleCounter will then "
        f"report loci files as missing even though they are present. "
        f"ROBIN tried to derive a world-readable image ({CLAIRS_TO_READABLE_IMAGE}) "
        f"while still running as the host UID/GID. Probe: {logs or status}. "
        "Do not work around this by running ClairS-TO as root."
    )
