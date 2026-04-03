"""Copy panel BED and reference genome for Nanopore sequencing setup."""

from __future__ import annotations

import shutil
import urllib.error
import urllib.request
from pathlib import Path
from typing import Optional
from urllib.parse import urlparse

# GRCh38 primary assembly (no alt), UCSC-style contig names — common default for alignment.
DEFAULT_GRCH38_REFERENCE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/"
    "GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/"
    "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
)


def panel_bed_filename(panel: str) -> str:
    return f"{panel}_panel_name_uniq.bed"


def panel_source_filename(panel: str) -> str:
    """Unprocessed panel BED packaged as ``{panel}_panel_source.bed`` (ship or ``robin add-panel``)."""
    return f"{panel}_panel_source.bed"


def panel_source_available(panel: str) -> bool:
    """True if ``{panel}_panel_source.bed`` exists in ``robin.resources``."""
    import importlib.resources as importlib_resources

    res = importlib_resources.files("robin.resources").joinpath(panel_source_filename(panel))
    try:
        return res.is_file()
    except Exception:
        return False


def copy_panel_source_bed_if_present(panel: str, dest_dir: Path) -> Optional[Path]:
    """Copy stored original BED if present; returns destination path or None."""
    import importlib.resources as importlib_resources
    from importlib.resources import as_file

    name = panel_source_filename(panel)
    res = importlib_resources.files("robin.resources").joinpath(name)
    if not res.is_file():
        return None
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / name
    with as_file(res) as src_path:
        shutil.copy2(src_path, dest)
    return dest


def copy_panel_bed_to(panel: str, dest_dir: Path) -> Path:
    """Copy the packaged panel BED into dest_dir; returns the destination path."""
    import importlib.resources as importlib_resources
    from importlib.resources import as_file

    name = panel_bed_filename(panel)
    res = importlib_resources.files("robin.resources").joinpath(name)
    if not res.is_file():
        raise FileNotFoundError(
            f"No packaged panel BED for '{panel}' (expected '{name}' in robin.resources)."
        )
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / name
    with as_file(res) as bed_path:
        shutil.copy2(bed_path, dest)
    return dest


def reference_looks_like_url(reference: str) -> bool:
    ref = reference.strip()
    return ref.startswith("http://") or ref.startswith("https://")


def _download_url(url: str, target_path: Path, *, label: str = "Downloading reference") -> None:
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "robin-sequencing-files/1.0"},
    )
    chunk_size = 1024 * 1024
    with urllib.request.urlopen(req) as resp:
        target_path.parent.mkdir(parents=True, exist_ok=True)
        total: Optional[int] = None
        try:
            cl = resp.headers.get("Content-Length")
            if cl:
                total = int(cl)
        except Exception:
            total = None

        try:
            import click  # type: ignore
        except Exception:
            click = None  # type: ignore

        downloaded = 0
        if click is not None and total and total > 0:
            with click.progressbar(length=total, label=label, show_eta=True, show_percent=True) as bar:
                with target_path.open("wb") as f:
                    while True:
                        chunk = resp.read(chunk_size)
                        if not chunk:
                            break
                        f.write(chunk)
                        downloaded += len(chunk)
                        bar.update(len(chunk))
        else:
            with target_path.open("wb") as f:
                while True:
                    chunk = resp.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)


def describe_reference_action(reference: str, dest_dir: Path) -> str:
    """Human-readable description of what will happen for the reference genome."""
    ref = reference.strip()
    if reference_looks_like_url(ref):
        parsed = urlparse(ref)
        basename = Path(parsed.path).name or "reference.fa.gz"
        return f"Download reference genome from:\n    {ref}\n  Save as: {dest_dir / basename}"
    path = Path(ref).expanduser().resolve()
    return f"Copy reference genome from:\n    {path}\n  To: {dest_dir / path.name}"


def materialize_reference(reference: str, dest_dir: Path) -> Path:
    """Download URL or copy local file into dest_dir; returns final path."""
    ref = reference.strip()
    dest_dir.mkdir(parents=True, exist_ok=True)

    if reference_looks_like_url(ref):
        parsed = urlparse(ref)
        basename = Path(parsed.path).name or "reference.fa.gz"
        dest = dest_dir / basename
        try:
            _download_url(ref, dest)
        except urllib.error.URLError as e:
            raise RuntimeError(f"Failed to download reference: {e}") from e
        return dest

    path = Path(ref).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Reference file not found: {path}")
    dest = dest_dir / path.name
    shutil.copy2(path, dest)
    return dest
