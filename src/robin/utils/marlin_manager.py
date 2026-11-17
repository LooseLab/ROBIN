"""Ensure MARLIN model and reference assets are available locally.

Reference files (probe list, class annotations, probe BEDs) ship with ROBIN under
``robin.resources.marlin``. The trained Keras model (~1.1 GiB) is downloaded from
Zenodo on first use.
"""

from __future__ import annotations

import logging
import os
import tempfile
import urllib.error
import urllib.request
from pathlib import Path
from typing import Optional

logger = logging.getLogger("robin.marlin")

MARLIN_MODEL_NAME = "marlin_v1.model.hdf5"
MARLIN_FEATURES_NAME = "marlin_v1.features.RData"
MARLIN_ANNOTATIONS_NAME = "marlin_v1.class_annotations.xlsx"
DEFAULT_GENOME_BUILD = "hg38"
SUPPORTED_GENOME_BUILDS = ("hg19", "hg38", "t2t")

# Public Zenodo record for the MARLIN v1 Keras model.
DEFAULT_MARLIN_MODEL_URL = (
    "https://zenodo.org/api/records/15565404/files/marlin_v1.model.hdf5/content"
)

# Env overrides for model path / download URL.
ENV_MODEL_PATH = "ROBIN_MARLIN_MODEL_PATH"
ENV_MODEL_URL = "ROBIN_MARLIN_MODEL_URL"
ENV_CACHE_DIR = "ROBIN_MARLIN_CACHE_DIR"


def get_marlin_resources_dir() -> Path:
    """Return the packaged MARLIN reference directory."""
    this_file = Path(__file__).resolve()
    resources_dir = this_file.parent.parent / "resources" / "marlin"
    if not resources_dir.is_dir():
        raise RuntimeError(f"MARLIN resources directory not found: {resources_dir}")
    return resources_dir


def get_marlin_cache_dir() -> Path:
    """Writable cache for the downloaded MARLIN model weights."""
    override = os.environ.get(ENV_CACHE_DIR)
    if override:
        path = Path(override).expanduser().resolve()
    else:
        path = Path.home() / ".cache" / "robin" / "marlin"
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_features_path() -> Path:
    path = get_marlin_resources_dir() / MARLIN_FEATURES_NAME
    if not path.is_file():
        raise FileNotFoundError(f"Missing MARLIN features file: {path}")
    return path


def get_annotations_path() -> Path:
    path = get_marlin_resources_dir() / MARLIN_ANNOTATIONS_NAME
    if not path.is_file():
        raise FileNotFoundError(f"Missing MARLIN annotations file: {path}")
    return path


def get_probes_bed_path(genome_build: str = DEFAULT_GENOME_BUILD) -> Path:
    build = (genome_build or DEFAULT_GENOME_BUILD).lower().strip()
    if build not in SUPPORTED_GENOME_BUILDS:
        raise ValueError(
            f"Unsupported MARLIN genome build {genome_build!r}; "
            f"expected one of {SUPPORTED_GENOME_BUILDS}"
        )
    path = get_marlin_resources_dir() / f"marlin_v1.probes_{build}.bed.gz"
    if not path.is_file():
        raise FileNotFoundError(f"Missing MARLIN probe BED: {path}")
    return path


def _file_ok(path: Path) -> bool:
    return path.exists() and path.is_file() and path.stat().st_size > 0


def _download_url_to_file(url: str, target_path: Path, *, timeout_s: int = 3600) -> None:
    """Download ``url`` to ``target_path`` atomically (temp file + rename)."""
    target_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(
        mode="wb",
        suffix=".part",
        prefix=target_path.name + ".",
        dir=str(target_path.parent),
        delete=False,
    ) as tmp:
        tmp_path = Path(tmp.name)

    try:
        logger.info("Downloading MARLIN model from %s", url)
        print(f"Downloading MARLIN model from {url}")
        req = urllib.request.Request(url, headers={"User-Agent": "robin-marlin/1.0"})
        with urllib.request.urlopen(req, timeout=timeout_s) as resp:
            total: Optional[int] = None
            try:
                content_length = resp.headers.get("Content-Length")
                if content_length:
                    total = int(content_length)
            except Exception:
                total = None

            try:
                import click  # type: ignore
            except Exception:
                click = None  # type: ignore

            chunk_size = 1024 * 1024
            downloaded = 0

            if click is not None and total and total > 0:
                with click.progressbar(
                    length=total,
                    label="Downloading MARLIN model",
                    show_eta=True,
                    show_percent=True,
                ) as bar:
                    with tmp_path.open("wb") as out:
                        while True:
                            chunk = resp.read(chunk_size)
                            if not chunk:
                                break
                            out.write(chunk)
                            downloaded += len(chunk)
                            bar.update(len(chunk))
            else:
                with tmp_path.open("wb") as out:
                    while True:
                        chunk = resp.read(chunk_size)
                        if not chunk:
                            break
                        out.write(chunk)
                        downloaded += len(chunk)
                        if total and downloaded % (50 * chunk_size) < chunk_size:
                            pct = 100.0 * downloaded / total
                            print(f"  MARLIN download: {pct:.1f}% ({downloaded}/{total} bytes)")

        if total is not None and downloaded != total:
            raise RuntimeError(
                f"Incomplete MARLIN download: got {downloaded} bytes, expected {total}"
            )
        if downloaded <= 0:
            raise RuntimeError("MARLIN download produced an empty file")

        tmp_path.replace(target_path)
        logger.info("MARLIN model saved to %s (%s bytes)", target_path, downloaded)
        print(f"MARLIN model saved to {target_path}")
    except Exception:
        try:
            if tmp_path.exists():
                tmp_path.unlink()
        except OSError:
            pass
        raise


def resolve_model_path(
    *,
    model_path: Optional[str | Path] = None,
    download_if_missing: bool = True,
    url: Optional[str] = None,
) -> Path:
    """
    Resolve the on-disk MARLIN Keras model path.

    Precedence:
    1. Explicit ``model_path`` argument
    2. ``ROBIN_MARLIN_MODEL_PATH`` env var
    3. Cache path under ``~/.cache/robin/marlin/`` (or ``ROBIN_MARLIN_CACHE_DIR``)

    When the file is missing and ``download_if_missing`` is True, download from
    Zenodo (or ``ROBIN_MARLIN_MODEL_URL`` / ``url``).
    """
    if model_path is not None:
        path = Path(model_path).expanduser().resolve()
    else:
        env_path = os.environ.get(ENV_MODEL_PATH)
        if env_path:
            path = Path(env_path).expanduser().resolve()
        else:
            path = get_marlin_cache_dir() / MARLIN_MODEL_NAME

    if _file_ok(path):
        return path

    if not download_if_missing:
        raise FileNotFoundError(
            f"MARLIN model not found at {path}. "
            f"Install robin[marlin], set {ENV_MODEL_PATH}, or allow auto-download."
        )

    download_url = url or os.environ.get(ENV_MODEL_URL) or DEFAULT_MARLIN_MODEL_URL
    try:
        _download_url_to_file(download_url, path)
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as exc:
        raise RuntimeError(
            f"Failed to download MARLIN model from {download_url} to {path}: {exc}"
        ) from exc

    if not _file_ok(path):
        raise RuntimeError(f"MARLIN model download completed but file is missing: {path}")
    return path


def ensure_marlin_assets(
    *,
    genome_build: str = DEFAULT_GENOME_BUILD,
    model_path: Optional[str | Path] = None,
    download_if_missing: bool = True,
) -> dict[str, Path]:
    """
    Ensure reference assets exist and the model is available.

    Returns a dict with keys: model, features, annotations, probes.
    """
    features = get_features_path()
    annotations = get_annotations_path()
    probes = get_probes_bed_path(genome_build)
    model = resolve_model_path(
        model_path=model_path, download_if_missing=download_if_missing
    )
    return {
        "model": model,
        "features": features,
        "annotations": annotations,
        "probes": probes,
    }
