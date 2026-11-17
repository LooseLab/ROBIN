"""Copy panel BED and reference genome for Nanopore sequencing setup."""

from __future__ import annotations

import shutil
import urllib.error
import urllib.request
from pathlib import Path
from typing import Optional
from urllib.parse import urlparse

# NCBI publishes the GRCh38 analysis FASTA as plain gzip. samtools faidx and
# other random-access readers need block gzip (bgzip) instead.
_BGZIP_NOTE = (
    "Plain gzip is re-encoded as bgzip so the FASTA can be indexed "
    "(samtools faidx / random access)."
)
_BGZIP_PROGRESS_MIN_BYTES = 8 * 1024 * 1024

# GRCh38 primary assembly (no alt), UCSC-style contig names — common default for alignment.
DEFAULT_GRCH38_REFERENCE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/"
    "GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/"
    "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
)


def panel_bed_filename(panel: str) -> str:
    return f"{panel}_panel_name_uniq.bed"


def resolve_panel_bed_path(panel: str) -> Optional[Path]:
    """Return the packaged panel BED path installed with ROBIN, if present."""
    return _resolve_packaged_panel_bed(panel, panel_bed_filename(panel))


def panel_stranded_bed_filename(panel: str) -> str:
    """Stranded panel BED for MinKNOW adaptive sampling (``{panel}_panel_source.bed``)."""
    return panel_source_filename(panel)


def resolve_panel_stranded_bed_path(
    panel: str,
    *,
    reference_path: Optional[str | Path] = None,
) -> Optional[Path]:
    """Resolve the stranded panel BED used for MinKNOW / adaptive sampling.

    Prefers ``{panel}_panel_source.bed`` next to the workflow reference FASTA
    (the layout produced by ``robin utils sequencing-files``), then the packaged
    copy shipped with ROBIN.
    """
    name = panel_stranded_bed_filename(panel)
    if reference_path not in (None, ""):
        candidate = Path(reference_path).expanduser().resolve().parent / name
        if candidate.is_file():
            return candidate
    return _resolve_packaged_panel_bed(panel, name)


def _resolve_packaged_panel_bed(panel: str, name: str) -> Optional[Path]:
    import importlib.resources as importlib_resources

    res = importlib_resources.files("robin.resources").joinpath(name)
    try:
        if not res.is_file():
            return None
    except Exception:
        return None

    try:
        from robin import resources

        bed_path = Path(resources.__file__).resolve().parent / name
        if bed_path.is_file():
            return bed_path
    except Exception:
        pass

    try:
        from importlib.resources import as_file

        with as_file(res) as bed_path:
            return Path(bed_path)
    except Exception:
        return None


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
        lines = [
            "Download reference genome from:",
            f"    {ref}",
            f"  Save as: {dest_dir / basename}",
        ]
        if basename.endswith(".gz"):
            lines.append(f"  {_BGZIP_NOTE}")
        return "\n".join(lines)
    path = Path(ref).expanduser().resolve()
    lines = [
        "Copy reference genome from:",
        f"    {path}",
        f"  To: {dest_dir / path.name}",
    ]
    if path.name.endswith(".gz"):
        lines.append(f"  {_BGZIP_NOTE}")
    return "\n".join(lines)


def _is_gzip(path: Path) -> bool:
    try:
        with path.open("rb") as handle:
            return handle.read(2) == b"\x1f\x8b"
    except OSError:
        return False


def _is_bgzf(path: Path) -> bool:
    """True when the file is block-gzip (bgzip), not plain gzip."""
    try:
        with path.open("rb") as handle:
            header = handle.read(12)
            if len(header) < 12 or header[:3] != b"\x1f\x8b\x08":
                return False
            # FEXTRA must be set; the extra field starts immediately after the header.
            if not (header[3] & 0x04):
                return False
            extra_len = int.from_bytes(header[10:12], "little")
            extra = handle.read(extra_len)
    except OSError:
        return False

    offset = 0
    while offset + 4 <= len(extra):
        subfield_id = extra[offset : offset + 2]
        subfield_len = int.from_bytes(extra[offset + 2 : offset + 4], "little")
        if subfield_id == b"BC" and subfield_len == 2:
            return True
        offset += 4 + subfield_len
    return False


class _CountingReader:
    """File wrapper that advances a click progress bar by compressed bytes read."""

    def __init__(self, raw, progress) -> None:
        self._raw = raw
        self._progress = progress
        length = getattr(progress, "length", None)
        self._remaining = length if isinstance(length, int) and length > 0 else None

    def read(self, n: int = -1) -> bytes:
        data = self._raw.read(n)
        if data and self._progress is not None and self._remaining:
            step = min(len(data), self._remaining)
            self._remaining -= step
            if step:
                self._progress.update(step)
        return data


def _recompress_gzip_as_bgzf(path: Path) -> None:
    """Replace a plain-gzip file with a bgzip file of the same decompressed bytes."""
    import gzip
    from contextlib import nullcontext

    import pysam

    tmp = path.with_name(path.name + ".bgzf.tmp")
    if tmp.exists():
        try:
            tmp.unlink()
        except OSError:
            pass

    src_size = path.stat().st_size
    try:
        click = None
        try:
            import click as click_mod

            click = click_mod
        except Exception:
            click = None

        with path.open("rb") as raw:
            show_progress = click is not None and src_size >= _BGZIP_PROGRESS_MIN_BYTES
            progress_cm = (
                click.progressbar(
                    length=src_size,
                    label="Re-encoding reference as bgzip",
                    show_eta=True,
                    show_percent=True,
                )
                if show_progress
                else nullcontext(None)
            )
            with progress_cm as progress:
                reader = _CountingReader(raw, progress)
                with gzip.GzipFile(fileobj=reader, mode="rb") as src, pysam.BGZFile(
                    str(tmp), "wb"
                ) as dst:
                    shutil.copyfileobj(src, dst, length=1024 * 1024)
        if not tmp.is_file() or tmp.stat().st_size == 0 or not _is_bgzf(tmp):
            raise RuntimeError("re-encoded file is not valid bgzip")
        tmp.replace(path)
    except Exception as exc:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass
        raise RuntimeError(
            f"Could not re-encode {path.name} as bgzip: {exc}. "
            "The original gzip file is unchanged."
        ) from exc


def ensure_reference_bgzf(path: Path) -> Path:
    """Rewrite a plain-gzip reference as bgzip in place.

    Uncompressed files and files that are already bgzip are left unchanged.
    The filename is kept (bgzip still uses a ``.gz`` suffix).
    """
    if not path.is_file() or not _is_gzip(path) or _is_bgzf(path):
        return path
    _recompress_gzip_as_bgzf(path)
    return path


def materialize_reference(reference: str, dest_dir: Path) -> Path:
    """Download URL or copy local file into dest_dir; returns final path.

    Gzipped FASTA that is not already bgzip is re-encoded as bgzip so tools
    such as ``samtools faidx`` can index it. The output name is unchanged.
    """
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
        return ensure_reference_bgzf(dest)

    path = Path(ref).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Reference file not found: {path}")
    dest = dest_dir / path.name
    shutil.copy2(path, dest)
    return ensure_reference_bgzf(dest)
