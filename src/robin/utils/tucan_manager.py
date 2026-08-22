"""Ensure Tucan model assets are available locally.

Tucan (UMCUGenetics) classifies pediatric solid tumors and lymphomas from sparse
Nanopore methylation calls. The pretrained ensemble is hosted on Hugging Face
(``MerelJongmans/model``). ROBIN downloads it on first use, packages ``model.zip``
in the format expected by ``tucan.cli.predict``, and derives a Sturgeon-style
probe BED for ``modkit_pileup_file_to_bed``.
"""

from __future__ import annotations

import logging
import os
import zipfile
from pathlib import Path
from typing import Optional

logger = logging.getLogger("robin.tucan")

DEFAULT_HF_REPO = "MerelJongmans/model"
MODEL_DIR_NAME = "model"
MODEL_ZIP_NAME = "model.zip"
PROBE_BED_NAME = "probe.bed"
MAPPING_PROBES_NAME = "probes_for_mapping.bed"
CHECKPOINTS_DIR_NAME = "checkpoints"
CLASSIFICATION_YAML_NAME = "classification_system.yaml"

ENV_MODEL_PATH = "ROBIN_TUCAN_MODEL_PATH"
ENV_MODEL_REPO = "ROBIN_TUCAN_MODEL_REPO"
ENV_CACHE_DIR = "ROBIN_TUCAN_CACHE_DIR"
ENV_NUM_CPGS = "ROBIN_TUCAN_NUM_CPGS"
ENV_NUM_SAMPLINGS = "ROBIN_TUCAN_NUM_SAMPLINGS"
ENV_PROBE_MARGIN = "ROBIN_TUCAN_PROBE_MARGIN"

DEFAULT_NUM_CPGS = 10_000
DEFAULT_NUM_SAMPLINGS = 1
# Exact chrom/position match to probe.bed (not Sturgeon's historic ±25 bp window).
DEFAULT_PROBE_MARGIN = 0


class TucanNotInstalledError(RuntimeError):
    """Raised when the external tucan package is not importable."""


def get_tucan_cache_dir() -> Path:
    """Writable cache for the downloaded Tucan model and derived assets."""
    override = os.environ.get(ENV_CACHE_DIR)
    if override:
        path = Path(override).expanduser().resolve()
    else:
        path = Path.home() / ".cache" / "robin" / "tucan"
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_default_num_cpgs() -> int:
    raw = (os.environ.get(ENV_NUM_CPGS) or "").strip()
    if not raw:
        return DEFAULT_NUM_CPGS
    try:
        value = int(raw)
    except ValueError as exc:
        raise ValueError(
            f"Invalid {ENV_NUM_CPGS}={raw!r}; expected an integer"
        ) from exc
    if value <= 0:
        raise ValueError(f"{ENV_NUM_CPGS} must be positive, got {value}")
    return value


def get_default_num_samplings() -> int:
    raw = (os.environ.get(ENV_NUM_SAMPLINGS) or "").strip()
    if not raw:
        return DEFAULT_NUM_SAMPLINGS
    try:
        value = int(raw)
    except ValueError as exc:
        raise ValueError(
            f"Invalid {ENV_NUM_SAMPLINGS}={raw!r}; expected an integer"
        ) from exc
    if value <= 0:
        raise ValueError(f"{ENV_NUM_SAMPLINGS} must be positive, got {value}")
    return value


def get_probe_margin() -> int:
    """
    Half-window (bp) when mapping modkit calls onto Tucan probes.

    Default is 0 (exact match). Set ``ROBIN_TUCAN_PROBE_MARGIN=25`` for
    Sturgeon-compatible ±25 bp aggregation.
    """
    raw = (os.environ.get(ENV_PROBE_MARGIN) or "").strip()
    if not raw:
        return DEFAULT_PROBE_MARGIN
    try:
        value = int(raw)
    except ValueError as exc:
        raise ValueError(
            f"Invalid {ENV_PROBE_MARGIN}={raw!r}; expected a non-negative integer"
        ) from exc
    if value < 0:
        raise ValueError(f"{ENV_PROBE_MARGIN} must be >= 0, got {value}")
    return value


def require_tucan_package() -> None:
    """Raise a clear error when ``tucan`` is not installed."""
    try:
        import tucan  # noqa: F401
    except ModuleNotFoundError as exc:
        raise TucanNotInstalledError(
            "The tucan package is not installed. Install the optional extra:\n"
            "  pip install 'robin[tucan]'\n"
            "or from source:\n"
            "  pip install git+https://github.com/UMCUGenetics/tucan.git"
        ) from exc


def _file_ok(path: Path) -> bool:
    return path.exists() and path.is_file() and path.stat().st_size > 0


def _dir_has_model_assets(model_dir: Path) -> bool:
    if not model_dir.is_dir():
        return False
    probe = model_dir / PROBE_BED_NAME
    yaml_path = model_dir / CLASSIFICATION_YAML_NAME
    checkpoints = model_dir / CHECKPOINTS_DIR_NAME
    if not (_file_ok(probe) and _file_ok(yaml_path) and checkpoints.is_dir()):
        return False
    # Ensemble of four checkpoints as used by tucan.cli.predict
    for i in range(4):
        if not _file_ok(checkpoints / f"checkpoint_{i}.pt"):
            return False
    return True


def _download_hf_model(model_dir: Path, *, repo_id: Optional[str] = None) -> Path:
    try:
        from huggingface_hub import snapshot_download
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "huggingface_hub is required to download the Tucan model. "
            "Install with: pip install 'robin[tucan]'"
        ) from exc

    repo = repo_id or os.environ.get(ENV_MODEL_REPO) or DEFAULT_HF_REPO
    model_dir.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Downloading Tucan model from Hugging Face repo %s", repo)
    print(f"Downloading Tucan model from Hugging Face ({repo})…")
    path = snapshot_download(
        repo_id=repo,
        repo_type="model",
        local_dir=str(model_dir),
        local_dir_use_symlinks=False,
    )
    resolved = Path(path)
    if not _dir_has_model_assets(resolved):
        raise RuntimeError(
            f"Tucan model download completed but expected assets are missing under {resolved}"
        )
    return resolved


def _create_model_zip(model_dir: Path, zip_path: Path) -> Path:
    """Create ``model.zip`` with a top-level ``model/`` folder (Tucan CLI convention)."""
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = zip_path.with_suffix(zip_path.suffix + ".partial")
    if tmp_path.exists():
        tmp_path.unlink()
    with zipfile.ZipFile(tmp_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for root, _dirs, files in os.walk(model_dir):
            for name in files:
                full = Path(root) / name
                # Archive paths must be model/<relative>
                arcname = Path(MODEL_DIR_NAME) / full.relative_to(model_dir)
                zf.write(full, arcname.as_posix())
    tmp_path.replace(zip_path)
    logger.info("Tucan model zip written to %s", zip_path)
    return zip_path


def write_mapping_probes_bed(tucan_probe_bed: Path, output_path: Path) -> Path:
    """
    Convert Tucan ``probe.bed`` (chrom/chromStart/chromEnd/name) to the
    Sturgeon-style probe table expected by ``modkit_pileup_file_to_bed``
    (chr/start/end/ID_REF).
    """
    import pandas as pd

    probes = pd.read_csv(tucan_probe_bed, sep="\t")
    rename = {}
    if "chrom" in probes.columns:
        rename["chrom"] = "chr"
    if "chromStart" in probes.columns:
        rename["chromStart"] = "start"
    if "chromEnd" in probes.columns:
        rename["chromEnd"] = "end"
    if "name" in probes.columns:
        rename["name"] = "ID_REF"
    if "probe_id" in probes.columns and "ID_REF" not in rename.values():
        rename["probe_id"] = "ID_REF"
    probes = probes.rename(columns=rename)

    required = {"chr", "start", "end", "ID_REF"}
    missing = required - set(probes.columns)
    if missing:
        raise ValueError(
            f"Tucan probe BED {tucan_probe_bed} missing columns {sorted(missing)}; "
            f"found {list(probes.columns)}"
        )

    out = probes[["chr", "start", "end", "ID_REF"]].copy()
    out["chr"] = out["chr"].astype(str).str.removeprefix("chr")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_path, sep="\t", index=False)
    return output_path


def resolve_model_zip(
    *,
    model_path: Optional[str | Path] = None,
    download_if_missing: bool = True,
) -> Path:
    """
    Resolve the on-disk Tucan ``model.zip`` path.

    Precedence:
    1. Explicit ``model_path`` argument
    2. ``ROBIN_TUCAN_MODEL_PATH`` env var
    3. Cache path under ``~/.cache/robin/tucan/model.zip``
    """
    if model_path is not None:
        path = Path(model_path).expanduser().resolve()
    else:
        env_path = os.environ.get(ENV_MODEL_PATH)
        if env_path:
            path = Path(env_path).expanduser().resolve()
        else:
            path = get_tucan_cache_dir() / MODEL_ZIP_NAME

    if _file_ok(path):
        return path

    if not download_if_missing:
        raise FileNotFoundError(
            f"Tucan model zip not found at {path}. "
            f"Install robin[tucan], set {ENV_MODEL_PATH}, or allow auto-download."
        )

    cache = get_tucan_cache_dir()
    model_dir = cache / MODEL_DIR_NAME
    if not _dir_has_model_assets(model_dir):
        _download_hf_model(model_dir)
    _create_model_zip(
        model_dir, path if path.parent == cache else cache / MODEL_ZIP_NAME
    )
    # If caller requested a custom path outside cache, copy/create there.
    final = path if path.suffix.lower() == ".zip" else cache / MODEL_ZIP_NAME
    if final != (cache / MODEL_ZIP_NAME) and _file_ok(cache / MODEL_ZIP_NAME):
        if not _file_ok(final):
            import shutil

            final.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(cache / MODEL_ZIP_NAME, final)
    if not _file_ok(final):
        raise RuntimeError(f"Tucan model zip missing after download: {final}")
    return final


def ensure_tucan_assets(
    *,
    model_path: Optional[str | Path] = None,
    download_if_missing: bool = True,
) -> dict[str, Path]:
    """
    Ensure model zip + mapping probe BED exist.

    Returns keys: ``model_zip``, ``model_dir``, ``probe_bed``, ``mapping_probes``.
    """
    cache = get_tucan_cache_dir()
    model_zip = resolve_model_zip(
        model_path=model_path, download_if_missing=download_if_missing
    )

    model_dir = cache / MODEL_DIR_NAME
    if not _dir_has_model_assets(model_dir):
        # Model zip may have been provided externally; extract into cache.
        extract_root = cache / "_extract"
        extract_root.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(model_zip, "r") as zf:
            zf.extractall(extract_root)
        candidate = extract_root / MODEL_DIR_NAME
        if not candidate.is_dir():
            # Zip may already be rooted at model contents
            candidate = extract_root
        if not _dir_has_model_assets(candidate):
            raise RuntimeError(
                f"Could not locate Tucan model assets after extracting {model_zip}"
            )
        if model_dir.exists() and model_dir.resolve() != candidate.resolve():
            import shutil

            if model_dir.is_dir():
                shutil.rmtree(model_dir)
            shutil.copytree(candidate, model_dir)
        elif not model_dir.exists():
            import shutil

            shutil.copytree(candidate, model_dir)

    probe_bed = model_dir / PROBE_BED_NAME
    mapping_probes = cache / MAPPING_PROBES_NAME
    if not _file_ok(mapping_probes):
        write_mapping_probes_bed(probe_bed, mapping_probes)

    return {
        "model_zip": model_zip,
        "model_dir": model_dir,
        "probe_bed": probe_bed,
        "mapping_probes": mapping_probes,
    }
