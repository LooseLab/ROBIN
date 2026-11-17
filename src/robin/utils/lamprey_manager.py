"""Ensure Lamprey model assets are available (research / evaluation use only).

Lamprey itself is **not** vendored in ROBIN. Users must install the Lamprey
package separately (and accept its evaluation license), e.g.::

    pip install git+ssh://git@github.com/princessmaximacenter/lamprey.git

Probe BEDs and ``classification_system.yaml`` are loaded from the installed
``lamprey`` package. The ONNX weights (~7 GiB) are downloaded from HuggingFace
on first use into ``~/.cache/robin/lamprey/`` only after research acknowledgment.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import Optional

logger = logging.getLogger("robin.lamprey")

DEFAULT_GENOME_BUILD = "hg38"
SUPPORTED_GENOME_BUILDS = ("hg38",)  # ROBIN Lamprey integration is hg38-only
DEFAULT_HF_REPO = "tachterberg/Lamprey"
MODEL_ONNX_NAME = "model.onnx"
BIN_CENTERS_NAME = "bin_centers.npy"
TEMPS_NAME = "temps.npy"

ENV_MODEL_PATH = "ROBIN_LAMPREY_MODEL_PATH"
ENV_MODEL_REPO = "ROBIN_LAMPREY_MODEL_REPO"
ENV_CACHE_DIR = "ROBIN_LAMPREY_CACHE_DIR"
ENV_RESEARCH_ACK = "ROBIN_LAMPREY_RESEARCH_ACK"

RESEARCH_DISCLAIMER = (
    "Lamprey is licensed for internal, non-commercial research and evaluation "
    "only. It must not be used for clinical care, diagnosis, or medical "
    "decision-making. See the Lamprey LICENSE (Oncode / Cyclomics / UMCU). "
    "Set ROBIN_LAMPREY_RESEARCH_ACK=1 to acknowledge this and enable model download."
)


class LampreyNotInstalledError(RuntimeError):
    """Raised when the external Lamprey package is not importable."""


class LampreyResearchAckRequired(RuntimeError):
    """Raised when research acknowledgment is required but missing."""


def research_ack_given() -> bool:
    value = (os.environ.get(ENV_RESEARCH_ACK) or "").strip().lower()
    return value in {"1", "true", "yes", "y", "on"}


def require_research_ack(*, for_download: bool = False) -> None:
    if research_ack_given():
        return
    action = "download the Lamprey model" if for_download else "run Lamprey"
    raise LampreyResearchAckRequired(
        f"Cannot {action} without research acknowledgment. {RESEARCH_DISCLAIMER}"
    )


def get_lamprey_package_dir() -> Path:
    """Return the installed ``lamprey`` package directory."""
    try:
        import lamprey
    except ModuleNotFoundError as exc:
        raise LampreyNotInstalledError(
            "The Lamprey package is not installed. Install it separately "
            "(research/evaluation license applies), e.g.\n"
            "  pip install git+ssh://git@github.com/princessmaximacenter/lamprey.git\n"
            "or from a local clone:\n"
            "  pip install -e /path/to/Lamprey\n"
            f"{RESEARCH_DISCLAIMER}"
        ) from exc
    return Path(lamprey.__file__).resolve().parent


def get_lamprey_files_dir() -> Path:
    path = get_lamprey_package_dir() / "files"
    if not path.is_dir():
        raise FileNotFoundError(f"Lamprey package files directory missing: {path}")
    return path


def get_classification_yaml_path() -> Path:
    path = get_lamprey_files_dir() / "classification_system.yaml"
    if not path.is_file():
        raise FileNotFoundError(f"Missing Lamprey classification YAML: {path}")
    return path


def get_probes_bed_path(genome_build: str = DEFAULT_GENOME_BUILD) -> Path:
    """Return the genome-coordinate probe BED (hg38 only in ROBIN)."""
    build = (genome_build or DEFAULT_GENOME_BUILD).lower().strip()
    if build not in SUPPORTED_GENOME_BUILDS:
        raise ValueError(
            f"ROBIN Lamprey integration supports only {SUPPORTED_GENOME_BUILDS}; "
            f"got {genome_build!r}"
        )
    path = get_lamprey_files_dir() / f"probe_{build}.bed"
    if not path.is_file():
        raise FileNotFoundError(f"Missing Lamprey probe BED: {path}")
    return path


def get_feature_order_bed_path() -> Path:
    """
    Return the probe BED that defines ONNX feature-vector order.

    Upstream Lamprey ``cli.predict`` always loads ``probe_t2t.bed`` for probe
    name order regardless of the preprocessing genome. We mirror that here.
    """
    path = get_lamprey_files_dir() / "probe_t2t.bed"
    if not path.is_file():
        raise FileNotFoundError(f"Missing Lamprey feature-order BED: {path}")
    return path


def get_lamprey_cache_dir() -> Path:
    override = os.environ.get(ENV_CACHE_DIR)
    if override:
        path = Path(override).expanduser().resolve()
    else:
        path = Path.home() / ".cache" / "robin" / "lamprey"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _file_ok(path: Path) -> bool:
    return path.exists() and path.is_file() and path.stat().st_size > 0


def _model_dir_complete(model_dir: Path) -> bool:
    return (
        _file_ok(model_dir / MODEL_ONNX_NAME)
        and _file_ok(model_dir / BIN_CENTERS_NAME)
        and _file_ok(model_dir / TEMPS_NAME)
    )


def resolve_model_dir(
    *,
    model_path: Optional[str | Path] = None,
    download_if_missing: bool = True,
    hf_repo: Optional[str] = None,
) -> Path:
    """
    Resolve a local directory containing ``model.onnx`` + calibration arrays.

    Precedence:
    1. Explicit ``model_path`` (directory)
    2. ``ROBIN_LAMPREY_MODEL_PATH``
    3. Cache under ``~/.cache/robin/lamprey/<repo>``

    HuggingFace download requires ``ROBIN_LAMPREY_RESEARCH_ACK=1``.
    """
    if model_path is not None:
        path = Path(model_path).expanduser().resolve()
    else:
        env_path = os.environ.get(ENV_MODEL_PATH)
        if env_path:
            path = Path(env_path).expanduser().resolve()
        else:
            repo = hf_repo or os.environ.get(ENV_MODEL_REPO) or DEFAULT_HF_REPO
            safe_name = repo.replace("/", "__")
            path = get_lamprey_cache_dir() / safe_name

    if path.is_dir() and _model_dir_complete(path):
        return path

    if path.is_file():
        raise ValueError(
            f"ROBIN_LAMPREY_MODEL_PATH must be a directory containing "
            f"{MODEL_ONNX_NAME}, not a file: {path}"
        )

    if not download_if_missing:
        raise FileNotFoundError(
            f"Lamprey model directory incomplete or missing: {path}. "
            f"Expected {MODEL_ONNX_NAME}, {BIN_CENTERS_NAME}, {TEMPS_NAME}."
        )

    require_research_ack(for_download=True)
    repo = hf_repo or os.environ.get(ENV_MODEL_REPO) or DEFAULT_HF_REPO
    try:
        from huggingface_hub import snapshot_download
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "huggingface_hub is required to download Lamprey models. "
            "Install with: pip install 'robin[lamprey]'"
        ) from exc

    logger.info("Downloading Lamprey model from HuggingFace: %s", repo)
    print(f"Downloading Lamprey model from HuggingFace: {repo}")
    print(RESEARCH_DISCLAIMER)
    path.mkdir(parents=True, exist_ok=True)
    downloaded = snapshot_download(repo_id=repo, local_dir=str(path))
    resolved = Path(downloaded).resolve()
    if not _model_dir_complete(resolved):
        raise RuntimeError(
            f"Lamprey download from {repo} completed but required files are missing in {resolved}"
        )
    return resolved


def ensure_lamprey_assets(
    *,
    genome_build: str = DEFAULT_GENOME_BUILD,
    model_path: Optional[str | Path] = None,
    download_if_missing: bool = True,
) -> dict[str, Path]:
    """
    Ensure Lamprey package files exist and the ONNX model is available.

    Returns keys: model_dir, classification_yaml, probes_bed, feature_order_bed.
    """
    # Installing Lamprey implies the user accepted its evaluation license; still
    # require explicit ack before auto-download of weights.
    get_lamprey_package_dir()
    classification_yaml = get_classification_yaml_path()
    probes_bed = get_probes_bed_path(genome_build)
    feature_order_bed = get_feature_order_bed_path()
    model_dir = resolve_model_dir(
        model_path=model_path, download_if_missing=download_if_missing
    )
    return {
        "model_dir": model_dir,
        "classification_yaml": classification_yaml,
        "probes_bed": probes_bed,
        "feature_order_bed": feature_order_bed,
    }
