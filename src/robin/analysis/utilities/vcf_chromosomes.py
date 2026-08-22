"""
Helpers for aligning VCF contig names with ClinVar (GRCh38) before SnpSift.

ClinVar's VCF uses ``1``, ``2``, … ``22``, ``X``, ``Y``, ``MT`` while Clair3 /
snpEff output in ROBIN typically uses UCSC-style ``chr1``, ``chr2``, ….
SnpSift tabix lookups fail when names do not match, which drops ClinVar INFO
fields and (with some SnpSift builds) can discard variants entirely.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Callable, TextIO

logger = logging.getLogger(__name__)

_CONTIG_ID_RE = re.compile(r"(##contig=<ID=)([^,>]+)")


def chrom_to_clinvar(chrom: str) -> str:
    """Map a sample VCF contig name to ClinVar / NCBI GRCh38 naming."""
    name = str(chrom or "").strip()
    if not name:
        return name
    if name.lower().startswith("chr"):
        name = name[3:]
    upper = name.upper()
    if upper in {"M", "MT"}:
        return "MT"
    return name


def chrom_to_ucsc(chrom: str) -> str:
    """Map a ClinVar-style contig back to UCSC ``chr*`` naming."""
    name = str(chrom or "").strip()
    if not name:
        return name
    if name.lower().startswith("chr"):
        return name if name.startswith("chr") else f"chr{name[3:]}"
    upper = name.upper()
    if upper == "MT":
        return "chrM"
    if upper in {"X", "Y"} or name.isdigit():
        return f"chr{name}"
    return name


def _rewrite_contig_header(line: str, transform: Callable[[str], str]) -> str:
    match = _CONTIG_ID_RE.search(line)
    if not match:
        return line
    new_id = transform(match.group(2))
    return f"{match.group(1)}{new_id}{line[match.end(2):]}"


def rewrite_vcf_chromosomes(
    src: Path | str,
    dst: Path | str,
    transform: Callable[[str], str],
) -> int:
    """
    Stream ``src`` to ``dst``, rewriting contig names on data lines and ##contig headers.

    Returns the number of variant records written.
    """
    src_path = Path(src)
    dst_path = Path(dst)
    dst_path.parent.mkdir(parents=True, exist_ok=True)

    variant_count = 0
    with (
        src_path.open("r", encoding="utf-8", errors="replace") as fin,
        dst_path.open("w", encoding="utf-8") as fout,
    ):
        variant_count = _stream_rewrite_vcf(fin, fout, transform)
    return variant_count


def _stream_rewrite_vcf(
    fin: TextIO,
    fout: TextIO,
    transform: Callable[[str], str],
) -> int:
    variant_count = 0
    for line in fin:
        if line.startswith("##contig"):
            fout.write(_rewrite_contig_header(line.rstrip("\n"), transform))
            fout.write("\n")
            continue
        if line.startswith("#"):
            fout.write(line)
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            fout.write(line)
            continue

        fields[0] = transform(fields[0])
        fout.write("\t".join(fields))
        fout.write("\n")
        variant_count += 1
    return variant_count


def normalize_vcf_chromosomes_for_clinvar(
    src: Path | str,
    dst: Path | str,
) -> dict[str, str]:
    """
    Rewrite ``src`` to ``dst`` using ClinVar contig names.

    Returns ``clinvar_chrom -> original_chrom`` for each contig observed.
    """
    mapping: dict[str, str] = {}

    def _transform(chrom: str) -> str:
        clinvar = chrom_to_clinvar(chrom)
        mapping.setdefault(clinvar, chrom)
        return clinvar

    count = rewrite_vcf_chromosomes(src, dst, _transform)
    logger.info(
        "Normalized VCF chromosomes for ClinVar: %s -> %s (%d variants, %d contigs)",
        src,
        dst,
        count,
        len(mapping),
    )
    return mapping


def restore_vcf_chromosomes(
    src: Path | str,
    dst: Path | str,
    clinvar_to_original: dict[str, str],
) -> int:
    """Restore sample contig naming after SnpSift using the normalize mapping."""

    def _transform(chrom: str) -> str:
        if chrom in clinvar_to_original:
            return clinvar_to_original[chrom]
        return chrom_to_ucsc(chrom)

    count = rewrite_vcf_chromosomes(src, dst, _transform)
    logger.info(
        "Restored VCF chromosomes after ClinVar annotation: %s -> %s (%d variants)",
        src,
        dst,
        count,
    )
    return count


def count_vcf_variants(path: Path | str) -> int:
    """Count non-header records in a VCF."""
    total = 0
    with Path(path).open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line and not line.startswith("#"):
                total += 1
    return total
