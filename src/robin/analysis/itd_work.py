"""
ITD / insertion calling in curated hotspots (nasvar-style CIGAR indel counts).

Mirrors the fusion stage → accumulate → finalize lifecycle, but extracts CIGAR
insertions inside panel-overlapping hotspots (FLT3, NPM1, …) rather than SA
chimeras.
"""

from __future__ import annotations

import json
import logging
import os
import shutil
from collections import defaultdict
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import pandas as pd
import pysam

from robin.analysis.master_bed_generator import FileLock
from robin.utils.sequencing_files import resolve_panel_bed_path

logger = logging.getLogger(__name__)

DEFAULT_HOTSPOTS_NAME = "itd_hotspots.hg38.json"

ITD_COUNT_COLUMNS = [
    "gene",
    "chrom",
    "position",
    "length",
    "support",
    "coverage",
    "label",
    "bam_path",
]

ITD_EVENT_COLUMNS = [
    "gene",
    "chrom",
    "position",
    "length",
    "support",
    # VAF denominator: mean gene depth from target coverage when available,
    # otherwise local hotspot spanning depth. Do not also emit hotspot_coverage /
    # gene_coverage — those duplicated support or coverage in the UI.
    "coverage",
    "vaf",
    "label",
    "exon_number",
    "transcript_id",
    "exon_id",
]

# Columns shown in GUI / PDF (also used to strip legacy duplicate fields).
ITD_EVENT_DISPLAY_COLUMNS = list(ITD_EVENT_COLUMNS)

ITD_EVENT_COLUMN_LABELS = {
    "gene": "Gene",
    "chrom": "Chrom",
    "position": "Position",
    "length": "Length",
    "support": "Support",
    "coverage": "Coverage",
    "vaf": "VAF",
    "label": "Label",
    "exon_number": "Exon",
    "transcript_id": "Transcript",
    "exon_id": "Exon ID",
}

_LEGACY_EVENT_COVERAGE_COLUMNS = ("hotspot_coverage", "gene_coverage")


def normalize_itd_events_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Return events with a single coverage column for display / reporting.

    Older CSVs may include ``hotspot_coverage`` and ``gene_coverage`` alongside
    ``coverage``; prefer gene mean depth when present, then drop the duplicates.
    """
    if df is None or df.empty:
        return pd.DataFrame(columns=ITD_EVENT_DISPLAY_COLUMNS)

    work = df.copy()
    if "gene_coverage" in work.columns:
        if "coverage" not in work.columns:
            work["coverage"] = work["gene_coverage"]
        else:
            gene_cov = pd.to_numeric(work["gene_coverage"], errors="coerce")
            coverage = pd.to_numeric(work["coverage"], errors="coerce")
            work["coverage"] = gene_cov.fillna(coverage)

    drop_cols = [c for c in _LEGACY_EVENT_COVERAGE_COLUMNS if c in work.columns]
    if drop_cols:
        work = work.drop(columns=drop_cols)

    ordered = [c for c in ITD_EVENT_DISPLAY_COLUMNS if c in work.columns]
    extras = [c for c in work.columns if c not in ordered and c != "__row_key"]
    return work[ordered + extras]


@dataclass(frozen=True)
class ExonFeature:
    """Canonical-transcript exon used for panel scan / event annotation."""

    start: int
    end: int
    number: Optional[str] = None
    transcript_id: Optional[str] = None
    exon_id: Optional[str] = None

    def contains(self, position: int) -> bool:
        return self.start <= position <= self.end


@dataclass(frozen=True)
class ItdHotspot:
    """Curated ITD / insertion window (nasvar ItdRegion equivalent)."""

    gene: str
    chrom: str
    start: int
    end: int
    min_length: int = 3
    min_frequency: float = 0.05
    min_supporting_reads: int = 2
    label: str = "ITD"
    transcript: Optional[str] = None
    # Optional sparse scan windows (inclusive). When set, BAM extraction only
    # fetches these intervals instead of the full [start, end] span.
    scan_intervals: Optional[Tuple[Tuple[int, int], ...]] = None
    # Canonical exons overlapping this window (for event annotation).
    exons: Optional[Tuple[ExonFeature, ...]] = None

    def contains(self, position: int) -> bool:
        if self.scan_intervals:
            return any(s <= position <= e for s, e in self.scan_intervals)
        return self.start <= position <= self.end

    def iter_scan_intervals(self) -> Tuple[Tuple[int, int], ...]:
        if self.scan_intervals:
            return self.scan_intervals
        return ((self.start, self.end),)

    def exon_at(self, position: int) -> Optional[ExonFeature]:
        """Return the overlapping exon feature at ``position``, if any."""
        if not self.exons:
            return None
        hits = [exon for exon in self.exons if exon.contains(position)]
        if not hits:
            return None
        # Prefer the tightest exon if multiple overlap (rare).
        return min(hits, key=lambda exon: (exon.end - exon.start, exon.start))


def default_hotspots_path() -> Path:
    """Return the packaged hg38 ITD hotspot JSON path."""
    try:
        from robin import resources

        path = Path(resources.__file__).resolve().parent / DEFAULT_HOTSPOTS_NAME
        if path.is_file():
            return path
    except Exception:
        pass
    return Path(__file__).resolve().parents[1] / "resources" / DEFAULT_HOTSPOTS_NAME


def load_itd_hotspots(path: Optional[str | Path] = None) -> Dict[str, ItdHotspot]:
    """Load gene → hotspot mapping from JSON."""
    hotspot_path = Path(path) if path else default_hotspots_path()
    with open(hotspot_path, encoding="utf-8") as handle:
        raw = json.load(handle)
    if not isinstance(raw, dict):
        raise ValueError(f"ITD hotspot config must be a JSON object: {hotspot_path}")

    hotspots: Dict[str, ItdHotspot] = {}
    for gene, entry in raw.items():
        if not isinstance(entry, Mapping):
            raise ValueError(f"Invalid hotspot entry for {gene!r}")
        hotspots[str(gene)] = ItdHotspot(
            gene=str(gene),
            chrom=str(entry["chrom"]),
            start=int(entry["start"]),
            end=int(entry["end"]),
            min_length=int(entry.get("min_length", 3)),
            min_frequency=float(entry.get("min_frequency", 0.05)),
            min_supporting_reads=int(entry.get("min_supporting_reads", 2)),
            label=str(entry.get("label", "ITD")),
            transcript=(
                str(entry["transcript"]) if entry.get("transcript") is not None else None
            ),
        )
    return hotspots


def gene_symbol_from_panel_name(name: str) -> str:
    """Extract a gene symbol from panel BED name fields (e.g. FLT3_NM_...)."""
    token = str(name).strip()
    if not token:
        return ""
    # Prefer the first underscore-delimited token when a transcript id is appended.
    return token.split("_", 1)[0].strip()


def _mean_coverage_by_gene_from_target_df(df: pd.DataFrame) -> Dict[str, float]:
    """Mean depth per gene symbol from ``target_coverage.csv`` rows."""
    if df.empty or "name" not in df.columns:
        return {}
    work = df.copy()
    if "coverage" not in work.columns:
        if "bases" in work.columns and "length" in work.columns:
            length = work["length"].replace(0, pd.NA)
            work["coverage"] = work["bases"] / length
        else:
            return {}
    out: Dict[str, list[float]] = defaultdict(list)
    for row in work.itertuples(index=False):
        gene = gene_symbol_from_panel_name(str(row.name)).upper()
        cov = getattr(row, "coverage", None)
        if gene and cov is not None and pd.notna(cov) and float(cov) > 0:
            out[gene].append(float(cov))
    return {gene: float(sum(vals) / len(vals)) for gene, vals in out.items()}


def load_gene_target_coverage(sample_dir: str | Path) -> Dict[str, float]:
    """
    Load mean per-gene depth from the target workflow ``target_coverage.csv``.

    This is the same mean coverage shown in the target GUI (e.g. FLT3 12.63x)
    and is used as the VAF denominator proxy for ITD events. Returns an empty
    dict when target analysis has not written coverage yet.
    """
    sample_dir = Path(sample_dir)
    coverage_csv = sample_dir / "target_coverage.csv"
    if not coverage_csv.is_file():
        return {}
    try:
        return _mean_coverage_by_gene_from_target_df(pd.read_csv(coverage_csv))
    except Exception as exc:
        logger.debug("Could not read %s: %s", coverage_csv, exc)
        return {}


def panel_gene_symbols(panel: str) -> set[str]:
    """Return gene symbols present in the packaged panel BED."""
    bed_path = resolve_panel_bed_path(panel)
    if bed_path is None or not bed_path.is_file():
        raise FileNotFoundError(f"Panel BED not found for target_panel={panel!r}")

    symbols: set[str] = set()
    with open(bed_path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            symbol = gene_symbol_from_panel_name(parts[3])
            if symbol:
                symbols.add(symbol.upper())
    return symbols


# Defaults for panel-wide scan windows (stricter than curated hotspots).
PANEL_SCAN_DEFAULTS = {
    "min_length": 4,
    "min_frequency": 0.05,
    "min_supporting_reads": 3,
    "label": "insertion",
}

VALID_ITD_REGION_MODES = frozenset({"hotspots", "panel", "both"})
DEFAULT_ANNOTATION_NAME = "gencode.v45.basic.annotation.gff3.gz"
_PANEL_EXON_CACHE: Dict[str, Dict[str, Tuple[ExonFeature, ...]]] = {}


def default_annotation_path() -> Path:
    """Return the packaged GENCODE basic annotation path."""
    try:
        from robin import resources

        path = Path(resources.__file__).resolve().parent / DEFAULT_ANNOTATION_NAME
        if path.is_file():
            return path
    except Exception:
        pass
    return Path(__file__).resolve().parents[1] / "resources" / DEFAULT_ANNOTATION_NAME


def _merge_intervals(
    intervals: Sequence[Tuple[int, int]],
    *,
    pad: int = 0,
) -> Tuple[Tuple[int, int], ...]:
    """Merge inclusive intervals with optional padding."""
    if not intervals:
        return ()
    ordered = sorted(
        (max(0, int(start) - pad), int(end) + pad) for start, end in intervals
    )
    merged: List[List[int]] = [[ordered[0][0], ordered[0][1]]]
    for start, end in ordered[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return tuple((start, end) for start, end in merged)


def _parse_gff_attrs(attrs_s: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for attr in attrs_s.strip().split(";"):
        if not attr or "=" not in attr:
            continue
        key, value = attr.split("=", 1)
        out[key] = value
    return out


def load_canonical_exon_features_for_panel(
    panel: str,
    *,
    annotation_path: Optional[str | Path] = None,
    exon_pad: int = 1,
) -> Dict[str, Tuple[ExonFeature, ...]]:
    """
    Load Ensembl-canonical exon features for genes on ``panel``.

    GENCODE GFF coordinates are converted from 1-based inclusive to 0-based
    inclusive so they match pysam CIGAR insertion anchors. Features keep
    exon_number / transcript_id for event annotation. An optional ``exon_pad``
    (default ±1) softens exon boundaries after conversion. Results are cached
    per panel in-process.
    """
    cache_key = f"{panel}|{exon_pad}|{annotation_path or ''}|0based"
    cached = _PANEL_EXON_CACHE.get(cache_key)
    if cached is not None:
        return cached

    import gzip

    symbols = {s.upper() for s in panel_gene_symbols(panel)}
    ann_path = Path(annotation_path) if annotation_path else default_annotation_path()
    if not ann_path.is_file():
        logger.warning(
            "Annotation not found at %s; panel exon filter disabled", ann_path
        )
        _PANEL_EXON_CACHE[cache_key] = {}
        return {}

    open_fn = gzip.open if str(ann_path).endswith(".gz") else open
    canonical_tx: set[str] = set()
    exons_by_gene: Dict[str, List[ExonFeature]] = defaultdict(list)

    with open_fn(ann_path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            attrs = _parse_gff_attrs(fields[8])
            if feature == "transcript":
                if "Ensembl_canonical" not in attrs.get("tag", ""):
                    continue
                gene_name = attrs.get("gene_name")
                transcript_id = attrs.get("transcript_id")
                if gene_name and transcript_id and gene_name.upper() in symbols:
                    canonical_tx.add(transcript_id)
            elif feature == "exon":
                gene_name = attrs.get("gene_name")
                transcript_id = attrs.get("transcript_id")
                if (
                    transcript_id in canonical_tx
                    and gene_name
                    and gene_name.upper() in symbols
                ):
                    # GFF3 = 1-based inclusive; pysam ref positions = 0-based.
                    start = int(fields[3]) - 1
                    end = int(fields[4]) - 1
                    if exon_pad:
                        start = max(0, start - int(exon_pad))
                        end = end + int(exon_pad)
                    exons_by_gene[gene_name.upper()].append(
                        ExonFeature(
                            start=start,
                            end=end,
                            number=attrs.get("exon_number"),
                            transcript_id=transcript_id,
                            exon_id=attrs.get("exon_id"),
                        )
                    )

    result: Dict[str, Tuple[ExonFeature, ...]] = {
        gene: tuple(sorted(features, key=lambda exon: (exon.start, exon.end)))
        for gene, features in exons_by_gene.items()
    }
    _PANEL_EXON_CACHE[cache_key] = result
    logger.info(
        "Loaded canonical exons for %d/%d panel genes from %s (0-based, pad=%d)",
        len(result),
        len(symbols),
        ann_path.name,
        exon_pad,
    )
    return result


def load_canonical_exon_intervals_for_panel(
    panel: str,
    *,
    annotation_path: Optional[str | Path] = None,
    exon_pad: int = 1,
) -> Dict[str, Tuple[Tuple[int, int], ...]]:
    """Compatibility wrapper: exon features collapsed to (start, end) intervals."""
    features = load_canonical_exon_features_for_panel(
        panel, annotation_path=annotation_path, exon_pad=exon_pad
    )
    return {
        gene: _merge_intervals([(exon.start, exon.end) for exon in exons], pad=0)
        for gene, exons in features.items()
    }


def _attach_exon_features(
    hotspots: Mapping[str, ItdHotspot],
    panel: str,
    *,
    exon_pad: int = 1,
) -> Dict[str, ItdHotspot]:
    """Attach overlapping canonical exon features onto hotspots (for annotation)."""
    try:
        feature_map = load_canonical_exon_features_for_panel(panel, exon_pad=exon_pad)
    except Exception as exc:
        logger.debug("Could not attach exon features: %s", exc)
        return dict(hotspots)

    out: Dict[str, ItdHotspot] = {}
    for gene, hotspot in hotspots.items():
        features = feature_map.get(gene.upper())
        if not features:
            out[gene] = hotspot
            continue
        clipped = tuple(
            exon
            for exon in features
            if exon.end >= hotspot.start and exon.start <= hotspot.end
        )
        out[gene] = replace(hotspot, exons=clipped or None)
    return out


def hotspots_from_panel(
    panel: str,
    *,
    min_length: int = PANEL_SCAN_DEFAULTS["min_length"],
    min_frequency: float = PANEL_SCAN_DEFAULTS["min_frequency"],
    min_supporting_reads: int = PANEL_SCAN_DEFAULTS["min_supporting_reads"],
    label: str = PANEL_SCAN_DEFAULTS["label"],
    intersect_exons: bool = True,
    exon_pad: int = 1,
) -> Dict[str, ItdHotspot]:
    """
    Build one ITD scan window per gene from the panel BED.

    When a gene appears in multiple BED rows, intervals are merged to the
    bounding box (min start, max end) on the first observed contig.

    When ``intersect_exons`` is true, scan windows are restricted to
    Ensembl-canonical exons (GFF→0-based, default ``exon_pad=1``) that fall
    inside the panel target interval, reducing intron/intergenic miscalls.
    """
    bed_path = resolve_panel_bed_path(panel)
    if bed_path is None or not bed_path.is_file():
        raise FileNotFoundError(f"Panel BED not found for target_panel={panel!r}")

    display_names: Dict[str, str] = {}
    merged_intervals: Dict[str, Tuple[str, int, int]] = {}
    with open(bed_path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0].strip()
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            if end < start:
                continue
            gene = gene_symbol_from_panel_name(parts[3])
            if not gene:
                continue
            key = gene.upper()
            display_names.setdefault(key, gene)
            if key in merged_intervals:
                prev_chrom, prev_start, prev_end = merged_intervals[key]
                if prev_chrom != chrom:
                    logger.debug(
                        "Skipping %s interval on %s; already merged on %s",
                        gene,
                        chrom,
                        prev_chrom,
                    )
                    continue
                merged_intervals[key] = (
                    prev_chrom,
                    min(prev_start, start),
                    max(prev_end, end),
                )
            else:
                merged_intervals[key] = (chrom, start, end)

    exon_features: Dict[str, Tuple[ExonFeature, ...]] = {}
    if intersect_exons:
        try:
            exon_features = load_canonical_exon_features_for_panel(
                panel, exon_pad=exon_pad
            )
        except Exception as exc:
            logger.warning(
                "Failed loading exon annotation for panel %s: %s; "
                "falling back to full gene windows",
                panel,
                exc,
            )
            exon_features = {}

    hotspots: Dict[str, ItdHotspot] = {}
    exon_restricted = 0
    for key, (chrom, start, end) in merged_intervals.items():
        gene = display_names.get(key, key)
        scan_intervals: Optional[Tuple[Tuple[int, int], ...]] = None
        clipped_exons: Optional[Tuple[ExonFeature, ...]] = None
        features = exon_features.get(key)
        if features:
            clipped_list = [
                ExonFeature(
                    start=max(start, exon.start),
                    end=min(end, exon.end),
                    number=exon.number,
                    transcript_id=exon.transcript_id,
                    exon_id=exon.exon_id,
                )
                for exon in features
                if exon.end >= start and exon.start <= end
            ]
            clipped_list = [exon for exon in clipped_list if exon.end >= exon.start]
            if clipped_list:
                clipped_exons = tuple(clipped_list)
                scan_intervals = _merge_intervals(
                    [(exon.start, exon.end) for exon in clipped_exons], pad=0
                )
                exon_restricted += 1
        hotspots[gene] = ItdHotspot(
            gene=gene,
            chrom=chrom,
            start=start,
            end=end,
            min_length=int(min_length),
            min_frequency=float(min_frequency),
            min_supporting_reads=int(min_supporting_reads),
            label=str(label),
            scan_intervals=scan_intervals,
            exons=clipped_exons,
        )
    total_scan_bp = sum(
        sum(e - s + 1 for s, e in h.iter_scan_intervals()) for h in hotspots.values()
    )
    logger.info(
        "ITD panel scan windows for %s: %d genes (%.1f Mb gene span, %.1f Mb scan, "
        "%d exon-restricted)",
        panel,
        len(hotspots),
        sum(h.end - h.start + 1 for h in hotspots.values()) / 1e6,
        total_scan_bp / 1e6,
        exon_restricted,
    )
    return hotspots


def resolve_itd_hotspots(
    panel: str,
    *,
    region_mode: str = "hotspots",
    hotspots_path: Optional[str | Path] = None,
    panel_min_length: int = PANEL_SCAN_DEFAULTS["min_length"],
    panel_min_frequency: float = PANEL_SCAN_DEFAULTS["min_frequency"],
    panel_min_supporting_reads: int = PANEL_SCAN_DEFAULTS["min_supporting_reads"],
) -> Dict[str, ItdHotspot]:
    """
    Resolve active ITD scan windows for a panel.

    Modes:
      - ``hotspots``: curated JSON windows filtered to genes on the panel
      - ``panel``: every gene interval from the panel BED
      - ``both``: panel genes, with curated windows overriding matching genes
        (keeps FLT3/NPM1/… at the known hotspot rather than the full gene span)
    """
    mode = str(region_mode or "hotspots").strip().lower()
    if mode not in VALID_ITD_REGION_MODES:
        raise ValueError(
            f"Invalid itd_region_mode={region_mode!r}; "
            f"expected one of {sorted(VALID_ITD_REGION_MODES)}"
        )

    curated = filter_hotspots_by_panel(load_itd_hotspots(hotspots_path), panel)

    if mode == "hotspots":
        return _attach_exon_features(curated, panel)

    panel_windows = hotspots_from_panel(
        panel,
        min_length=panel_min_length,
        min_frequency=panel_min_frequency,
        min_supporting_reads=panel_min_supporting_reads,
    )
    if mode == "panel":
        return panel_windows

    # both: curated overrides same-named genes (case-insensitive).
    by_upper = {g.upper(): h for g, h in panel_windows.items()}
    for gene, hotspot in curated.items():
        by_upper[gene.upper()] = hotspot
    # Prefer curated gene-key casing when overridden.
    out: Dict[str, ItdHotspot] = {}
    for gene, hotspot in panel_windows.items():
        key = gene.upper()
        chosen = by_upper[key]
        out[chosen.gene] = chosen
    for gene, hotspot in curated.items():
        out[hotspot.gene] = hotspot
    logger.info(
        "ITD region mode=both: %d panel windows, %d curated overrides → %d active",
        len(panel_windows),
        len(curated),
        len(out),
    )
    # Re-attach exons so curated overrides still carry exon annotation.
    return _attach_exon_features(out, panel)


def filter_hotspots_by_panel(
    hotspots: Mapping[str, ItdHotspot],
    panel: str,
) -> Dict[str, ItdHotspot]:
    """Keep hotspots whose gene symbol appears in the active target panel."""
    symbols = panel_gene_symbols(panel)
    filtered = {
        gene: hotspot
        for gene, hotspot in hotspots.items()
        if gene.upper() in symbols
    }
    logger.info(
        "ITD hotspots after panel filter (%s): %s (from %d configured)",
        panel,
        sorted(filtered),
        len(hotspots),
    )
    return filtered


def resolve_bam_contig(bam: pysam.AlignmentFile, chrom: str) -> Optional[str]:
    """Map a hotspot chrom name onto a contig present in the BAM header."""
    references = set(bam.references)
    if chrom in references:
        return chrom
    if chrom.startswith("chr"):
        alt = chrom[3:]
        if alt in references:
            return alt
    else:
        alt = f"chr{chrom}"
        if alt in references:
            return alt
    return None


def extract_indels_in_region(
    bam: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
) -> Tuple[Dict[int, Dict[int, int]], Dict[int, int]]:
    """
    Count CIGAR insertions and spanning coverage in [start, end].

    Coverage is tracked sparsely (only at indel anchor positions) so this
    scales to full gene / panel intervals, not just small curated hotspots.

    Returns:
        indels: position → {signed_length → read count}
            (+len = insertion, -len = deletion; deletions are retained for
            completeness but not used in ITD event calling)
        coverage: position → unique spanning read count at indel anchors
    """
    if end < start:
        raise ValueError(f"Invalid region {chrom}:{start}-{end}")

    indel_reads: Dict[int, Dict[int, set[str]]] = defaultdict(lambda: defaultdict(set))
    # Clipped ref-coverage intervals for primary reads overlapping the region.
    read_intervals: List[Tuple[int, int, str]] = []

    bam_chrom = resolve_bam_contig(bam, chrom)
    if bam_chrom is None:
        logger.warning("Contig %s not found in BAM; skipping region", chrom)
        return {}, {}

    for read in bam.fetch(bam_chrom, start, end + 1):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.cigartuples is None or read.reference_start is None:
            continue

        qname = read.query_name or f"unnamed:{id(read)}"
        ref_pos = int(read.reference_start)
        for op, length in read.cigartuples:
            # pysam: M=0, I=1, D=2, N=3, S=4, H=5, P=6, =7, X=8
            if op in (0, 7, 8):  # match / sequence match / mismatch
                seg_start = max(start, ref_pos)
                seg_end = min(end, ref_pos + length - 1)
                if seg_end >= seg_start:
                    read_intervals.append((seg_start, seg_end, qname))
                ref_pos += length
            elif op == 1:  # insertion
                if start <= ref_pos <= end:
                    indel_reads[ref_pos][length].add(qname)
            elif op in (2, 3):  # deletion / skip
                seg_start = max(start, ref_pos)
                seg_end = min(end, ref_pos + length - 1)
                if seg_end >= seg_start:
                    read_intervals.append((seg_start, seg_end, qname))
                if start <= ref_pos <= end:
                    indel_reads[ref_pos][-length].add(qname)
                ref_pos += length
            elif op == 6:  # pad
                ref_pos += length
            # soft/hard clip: ignore

    indels = {
        pos: {length: len(names) for length, names in lengths.items()}
        for pos, lengths in indel_reads.items()
    }

    # Spanning depth only at indel anchors (O(indels × reads), not O(region length)).
    coverage: Dict[int, int] = {}
    for pos, lengths in indel_reads.items():
        names = set()
        for s, e, qn in read_intervals:
            if s <= pos <= e:
                names.add(qn)
        # Insertion-only anchors still count their carriers.
        for carriers in lengths.values():
            names.update(carriers)
        coverage[pos] = len(names)
    return indels, coverage


def _depth_at(
    coverage: Mapping[int, int] | Sequence[int],
    position: int,
    region_start: int,
) -> int:
    """Look up spanning depth from sparse dict or dense list coverage."""
    if isinstance(coverage, Mapping):
        return int(coverage.get(position, 0))
    cov_idx = position - region_start
    if 0 <= cov_idx < len(coverage):
        return int(coverage[cov_idx])
    return 0


def _add_spanning_reads(
    spanning: List[set[str]],
    region_start: int,
    region_end: int,
    ref_pos: int,
    length: int,
    qname: str,
) -> None:
    """Legacy helper retained for tests that build dense coverage lists."""
    seg_start = max(region_start, ref_pos)
    seg_end = min(region_end, ref_pos + length - 1)
    if seg_end < seg_start:
        return
    for pos in range(seg_start, seg_end + 1):
        spanning[pos - region_start].add(qname)


def call_events_from_counts(
    hotspot: ItdHotspot,
    indels: Mapping[int, Mapping[int, int]],
    coverage: Mapping[int, int] | Sequence[int],
    *,
    use_nasvar_frequency: bool = False,
    gene_mean_coverage: Optional[float] = None,
) -> List[Dict[str, Any]]:
    """
    Threshold, cluster nearby similar-length insertions, and compute VAF.

    Only positive lengths (insertions) are reported.

    VAF prefers mean gene depth from the target workflow
    (``gene_mean_coverage``, e.g. FLT3 12.63x) as a proxy for the number of
    reads expected to span the insertion locus; that value is stored as
    ``coverage``. Local hotspot spanning counts are only used when gene mean
    depth is unavailable. Nearby alleles (similar length / position) are
    collapsed into one event; support is the sum of member length-bin counts
    (each bin once).
    """
    region_start = hotspot.start

    # Seed candidates that pass per-bin thresholds.
    candidates: List[Dict[str, Any]] = []
    for position in sorted(indels):
        if not hotspot.contains(position):
            continue
        depth = _depth_at(coverage, position, region_start)
        for length, support in sorted(indels[position].items()):
            if length < hotspot.min_length:
                continue
            if support < hotspot.min_supporting_reads:
                continue
            candidates.append(
                {
                    "gene": hotspot.gene,
                    "chrom": hotspot.chrom,
                    "position": int(position),
                    "length": int(length),
                    "support": int(support),
                    "coverage": int(depth),
                    "label": hotspot.label,
                }
            )

    if not candidates:
        return []

    # Cluster nearby / similar-length alleles (nasvar proximity window).
    candidates.sort(key=lambda row: (-row["support"], row["position"], -row["length"]))
    used = [False] * len(candidates)
    clustered: List[Dict[str, Any]] = []

    for i, seed in enumerate(candidates):
        if used[i]:
            continue
        members = [seed]
        used[i] = True
        changed = True
        while changed:
            changed = False
            for j, other in enumerate(candidates):
                if used[j]:
                    continue
                if any(_alleles_mergeable(member, other) for member in members):
                    members.append(other)
                    used[j] = True
                    changed = True

        # Representative = highest-support member; support = sum of unique bins.
        representative = max(members, key=lambda row: (row["support"], -row["position"]))
        support = sum(member["support"] for member in members)
        # Spanning depth at the representative anchor (includes ref + alt reads).
        depth = int(representative["coverage"])
        clustered.append(
            {
                "gene": representative["gene"],
                "chrom": representative["chrom"],
                "position": representative["position"],
                "length": representative["length"],
                "support": int(support),
                "coverage": int(depth),
                "label": representative["label"],
            }
        )

    filtered: List[Dict[str, Any]] = []
    for event in clustered:
        support = float(event["support"])
        hotspot_depth = float(event["coverage"])
        # Prefer mean gene depth for VAF; fall back to local hotspot depth.
        if gene_mean_coverage is not None and gene_mean_coverage > 0:
            depth = float(gene_mean_coverage)
        else:
            depth = hotspot_depth
        event["coverage"] = round(depth, 4)
        if depth <= 0:
            continue
        if use_nasvar_frequency:
            # nasvar: support / (coverage - support); only defined when coverage > support.
            denom = depth - support
            if denom <= 0:
                continue
            vaf = support / denom
        else:
            vaf = support / depth
        # Safety clamp if support exceeds the depth proxy.
        vaf = min(vaf, 1.0)
        if hotspot.min_frequency > 0 and vaf < hotspot.min_frequency:
            continue
        event["vaf"] = round(vaf, 6)
        exon = hotspot.exon_at(int(event["position"]))
        event["exon_number"] = exon.number if exon else None
        event["transcript_id"] = (
            (exon.transcript_id if exon else None) or hotspot.transcript
        )
        event["exon_id"] = exon.exon_id if exon else None
        filtered.append(event)

    filtered.sort(key=lambda row: (-row["support"], row["position"], -row["length"]))
    return filtered


def _alleles_mergeable(a: Mapping[str, Any], b: Mapping[str, Any]) -> bool:
    """True when two insertion alleles fall in the nasvar nearby-merge window."""
    if a["position"] == b["position"] and a["length"] == b["length"]:
        return False
    length_a = float(a["length"])
    length_b = float(b["length"])
    if abs(int(a["position"]) - int(b["position"])) >= int(length_a * 1.1):
        return False
    # Lengths within ±10% of each other (symmetric check).
    return (0.9 * length_b < length_a < 1.1 * length_b) and (
        0.9 * length_a < length_b < 1.1 * length_a
    )


def process_bam_itd_pass(
    bam_path: str,
    hotspots: Mapping[str, ItdHotspot],
) -> pd.DataFrame:
    """Scan one BAM over active hotspots; return per-(gene,pos,len) count rows."""
    rows: List[Dict[str, Any]] = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for gene, hotspot in hotspots.items():
            # Merge indel/coverage across sparse scan intervals for this gene.
            merged_indels: Dict[int, Dict[int, int]] = defaultdict(dict)
            merged_coverage: Dict[int, int] = {}
            for seg_start, seg_end in hotspot.iter_scan_intervals():
                indels, coverage = extract_indels_in_region(
                    bam, hotspot.chrom, seg_start, seg_end
                )
                for position, lengths in indels.items():
                    if not hotspot.contains(position):
                        continue
                    for length, support in lengths.items():
                        if length <= 0:
                            continue
                        prev = merged_indels[position].get(length, 0)
                        merged_indels[position][length] = max(prev, int(support))
                    merged_coverage[position] = max(
                        merged_coverage.get(position, 0),
                        int(coverage.get(position, 0)),
                    )
            for position, lengths in merged_indels.items():
                depth = int(merged_coverage.get(position, 0))
                for length, support in lengths.items():
                    rows.append(
                        {
                            "gene": gene,
                            "chrom": hotspot.chrom,
                            "position": int(position),
                            "length": int(length),
                            "support": int(support),
                            "coverage": depth,
                            "label": hotspot.label,
                            "bam_path": os.path.basename(bam_path),
                        }
                    )
    if not rows:
        return pd.DataFrame(columns=ITD_COUNT_COLUMNS)
    return pd.DataFrame(rows, columns=ITD_COUNT_COLUMNS)


def _sample_dir(work_dir: str, sample_id: str) -> str:
    path = os.path.join(work_dir, sample_id)
    os.makedirs(path, exist_ok=True)
    return path


def _staging_dir(work_dir: str, sample_id: str) -> str:
    path = os.path.join(_sample_dir(work_dir, sample_id), "_itd_staging")
    os.makedirs(path, exist_ok=True)
    return path


def _dataset_dir(work_dir: str, sample_id: str) -> str:
    path = os.path.join(_sample_dir(work_dir, sample_id), "itd_candidates_dataset")
    os.makedirs(path, exist_ok=True)
    return path


def _atomic_counter_increment(work_dir: str, sample_id: str) -> int:
    sample_dir = _sample_dir(work_dir, sample_id)
    lock_dir = os.path.join(sample_dir, "_locks")
    os.makedirs(lock_dir, exist_ok=True)
    lock_file = os.path.join(lock_dir, "itd_counter.lock")
    counter_file = os.path.join(sample_dir, "itd_analysis_counter.txt")

    with FileLock(lock_file, timeout=30.0):
        counter = 0
        if os.path.exists(counter_file):
            try:
                with open(counter_file, encoding="utf-8") as handle:
                    counter = int(handle.read().strip() or "0")
            except (ValueError, OSError):
                counter = 0
        with open(counter_file, "w", encoding="utf-8") as handle:
            handle.write(str(counter + 1))
        return counter


def _increment_pending_count(work_dir: str, sample_id: str, delta: int = 1) -> int:
    sample_dir = _sample_dir(work_dir, sample_id)
    lock_dir = os.path.join(sample_dir, "_locks")
    os.makedirs(lock_dir, exist_ok=True)
    lock_file = os.path.join(lock_dir, "itd_pending.lock")
    count_file = os.path.join(sample_dir, "itd_pending_count.txt")

    with FileLock(lock_file, timeout=30.0):
        current = 0
        if os.path.exists(count_file):
            try:
                with open(count_file, encoding="utf-8") as handle:
                    current = int(handle.read().strip() or "0")
            except (ValueError, OSError):
                current = 0
        new_count = max(0, current + int(delta))
        with open(count_file, "w", encoding="utf-8") as handle:
            handle.write(str(new_count))
        return new_count


def process_bam_itd_with_staging(
    bam_path: str,
    work_dir: str,
    sample_id: str,
    hotspots: Mapping[str, ItdHotspot],
    *,
    batch_size: int = 20,
) -> Tuple[Dict[str, Any], bool]:
    """Process one BAM into ITD staging; return (stats, should_accumulate)."""
    counts = process_bam_itd_pass(bam_path, hotspots)
    counter = _atomic_counter_increment(work_dir, sample_id)
    staging_dir = _staging_dir(work_dir, sample_id)
    staging_path = os.path.join(staging_dir, f"itd_{counter:06d}.parquet")

    if not counts.empty:
        counts.to_parquet(staging_path, index=False, engine="pyarrow", compression="snappy")
        logger.info(
            "ITD staging: wrote %d indel count rows from %s",
            len(counts),
            os.path.basename(bam_path),
        )
    else:
        # Touch an empty marker so pending accounting stays aligned with fusion.
        counts.to_parquet(staging_path, index=False, engine="pyarrow", compression="snappy")
        logger.debug("ITD staging: no insertions in hotspots for %s", bam_path)

    pending = _increment_pending_count(work_dir, sample_id, delta=1)
    should_accumulate = pending >= batch_size
    return {
        "rows": int(len(counts)),
        "staging_path": staging_path,
        "pending": pending,
    }, should_accumulate


def _aggregate_count_frames(frames: Iterable[pd.DataFrame]) -> pd.DataFrame:
    """Sum support and coverage across identical (gene, chrom, pos, length) bins.

    Coverage is summed across BAM parts (pooled depth). Per-position coverage is
    later floored by insertion support when calling events so VAF stays ≤ 1.
    """
    pieces = [frame for frame in frames if frame is not None and not frame.empty]
    if not pieces:
        return pd.DataFrame(columns=ITD_COUNT_COLUMNS)

    combined = pd.concat(pieces, ignore_index=True)
    grouped = (
        combined.groupby(["gene", "chrom", "position", "length", "label"], as_index=False)
        .agg(support=("support", "sum"), coverage=("coverage", "sum"))
    )
    grouped["bam_path"] = ""
    return grouped[ITD_COUNT_COLUMNS]


def _events_from_aggregated_counts(
    counts: pd.DataFrame,
    hotspots: Mapping[str, ItdHotspot],
    gene_coverage: Optional[Mapping[str, float]] = None,
) -> pd.DataFrame:
    if counts.empty:
        return pd.DataFrame(columns=ITD_EVENT_COLUMNS)

    event_rows: List[Dict[str, Any]] = []
    for gene, group in counts.groupby("gene"):
        hotspot = hotspots.get(str(gene))
        if hotspot is None:
            continue
        indels: Dict[int, Dict[int, int]] = defaultdict(dict)
        coverage: Dict[int, int] = {}
        for row in group.itertuples(index=False):
            position = int(row.position)
            length = int(row.length)
            support = int(row.support)
            depth = int(row.coverage)
            # Same (pos, length) should already be aggregated; keep max if duplicated.
            indels[position][length] = max(indels[position].get(length, 0), support)
            # Pooled depth across length bins at this position should be shared.
            coverage[position] = max(coverage.get(position, 0), depth)
        gene_depth: Optional[float] = None
        if gene_coverage:
            gene_depth = gene_coverage.get(str(gene).upper()) or gene_coverage.get(
                str(gene)
            )
        event_rows.extend(
            call_events_from_counts(
                hotspot, indels, coverage, gene_mean_coverage=gene_depth
            )
        )

    if not event_rows:
        return pd.DataFrame(columns=ITD_EVENT_COLUMNS)
    return pd.DataFrame(event_rows, columns=ITD_EVENT_COLUMNS)


def accumulate_itd_candidates(
    work_dir: str,
    sample_id: str,
    hotspots: Mapping[str, ItdHotspot],
    *,
    force: bool = False,
    batch_size: int = 20,
) -> Dict[str, Any]:
    """Move staged ITD count parts into the append-only dataset and rebuild reports."""
    sample_dir = _sample_dir(work_dir, sample_id)
    lock_dir = os.path.join(sample_dir, "_locks")
    os.makedirs(lock_dir, exist_ok=True)
    lock_file = os.path.join(lock_dir, "itd_accumulate.lock")

    with FileLock(lock_file, timeout=60.0):
        pending_file = os.path.join(sample_dir, "itd_pending_count.txt")
        pending = 0
        if os.path.exists(pending_file):
            try:
                with open(pending_file, encoding="utf-8") as handle:
                    pending = int(handle.read().strip() or "0")
            except (ValueError, OSError):
                pending = 0

        if not force and pending < batch_size:
            return {"accumulated": False, "pending": pending, "events": 0}

        staging_dir = _staging_dir(work_dir, sample_id)
        staging_files = sorted(
            str(path)
            for path in Path(staging_dir).glob("itd_*.parquet")
            if path.is_file()
        )
        if not staging_files and not force:
            return {"accumulated": False, "pending": pending, "events": 0}

        frames = []
        for path in staging_files:
            try:
                frames.append(pd.read_parquet(path))
            except Exception as exc:
                logger.warning("Could not read ITD staging file %s: %s", path, exc)

        batch = _aggregate_count_frames(frames)
        dataset_dir = _dataset_dir(work_dir, sample_id)
        if not batch.empty:
            part_id = len(list(Path(dataset_dir).glob("part_*.parquet")))
            part_path = os.path.join(dataset_dir, f"part_{part_id:06d}.parquet")
            batch.to_parquet(part_path, index=False, engine="pyarrow", compression="snappy")

        for path in staging_files:
            try:
                os.remove(path)
            except OSError:
                pass

        with open(pending_file, "w", encoding="utf-8") as handle:
            handle.write("0")

        # Rebuild from full dataset
        dataset_parts = sorted(Path(dataset_dir).glob("part_*.parquet"))
        dataset_frames = []
        for path in dataset_parts:
            try:
                dataset_frames.append(pd.read_parquet(path))
            except Exception as exc:
                logger.warning("Could not read ITD dataset part %s: %s", path, exc)

        aggregated = _aggregate_count_frames(dataset_frames)
        gene_coverage_map = load_gene_target_coverage(sample_dir)
        events = _events_from_aggregated_counts(
            aggregated, hotspots, gene_coverage=gene_coverage_map
        )
        events_path = os.path.join(sample_dir, "itd_events.csv")
        summary_path = os.path.join(sample_dir, "itd_summary.csv")
        events.to_csv(events_path, index=False)

        summary_rows = []
        for gene, hotspot in hotspots.items():
            gene_events = events[events["gene"] == gene] if not events.empty else events
            summary_rows.append(
                {
                    "gene": gene,
                    "label": hotspot.label,
                    "chrom": hotspot.chrom,
                    "start": hotspot.start,
                    "end": hotspot.end,
                    "n_events": int(len(gene_events)),
                    "max_support": int(gene_events["support"].max()) if len(gene_events) else 0,
                    "max_vaf": float(gene_events["vaf"].max()) if len(gene_events) else 0.0,
                }
            )
        pd.DataFrame(summary_rows).to_csv(summary_path, index=False)

        # Compact hotspot config used for this sample (debug / provenance).
        # Omit exon feature lists — they can be megabytes in panel mode.
        config_path = os.path.join(sample_dir, "itd_hotspots_used.json")
        with open(config_path, "w", encoding="utf-8") as handle:
            payload = {}
            for gene, hotspot in hotspots.items():
                entry = asdict(hotspot)
                entry.pop("exons", None)
                if entry.get("scan_intervals") and len(entry["scan_intervals"]) > 20:
                    entry["scan_intervals"] = (
                        f"<{len(hotspot.scan_intervals)} intervals>"
                    )
                payload[gene] = entry
            json.dump(payload, handle, indent=2)

        logger.info(
            "ITD accumulation for %s: %d events written to %s",
            sample_id,
            len(events),
            events_path,
        )

        qc_paths: Dict[str, Any] = {}
        if not events.empty:
            try:
                bam_paths = sorted(
                    str(path)
                    for path in Path(sample_dir).glob("batch_*.bam")
                    if path.is_file()
                )
                if not bam_paths:
                    bam_paths = sorted(
                        str(path)
                        for path in Path(sample_dir).glob("*.bam")
                        if path.is_file() and not path.name.endswith(".bai")
                    )
                if bam_paths:
                    qc_paths = write_itd_read_qc(
                        sample_dir, events, bam_paths=bam_paths
                    )
            except Exception as exc:
                logger.warning("ITD read QC failed for %s: %s", sample_id, exc)

        return {
            "accumulated": True,
            "pending": 0,
            "events": int(len(events)),
            "events_path": events_path,
            "summary_path": summary_path,
            "parts_ingested": len(staging_files),
            **qc_paths,
        }


def _softclip_bases(read: pysam.AlignedSegment) -> Tuple[int, int]:
    left = right = 0
    cigars = read.cigartuples or []
    if cigars and cigars[0][0] == 4:
        left = int(cigars[0][1])
    if cigars and cigars[-1][0] == 4:
        right = int(cigars[-1][1])
    return left, right


def _insertion_base_qualities(
    read: pysam.AlignedSegment,
    anchor_pos: int,
    ins_len: int,
) -> List[int]:
    """Return base qualities for the insertion of ``ins_len`` anchored at ``anchor_pos``."""
    if read.cigartuples is None or read.reference_start is None:
        return []
    quals = read.query_qualities
    if quals is None:
        return []
    ref_pos = int(read.reference_start)
    query_pos = 0
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            ref_pos += length
            query_pos += length
        elif op == 1:
            if ref_pos == anchor_pos and length == ins_len:
                return [int(quals[query_pos + i]) for i in range(length)]
            query_pos += length
        elif op in (2, 3):
            ref_pos += length
        elif op == 4:  # soft clip
            query_pos += length
        elif op == 6:  # pad
            ref_pos += length
        # hard clip (5): neither advances
    return []


def _length_compatible(observed: int, target: int) -> bool:
    """True when insertion lengths fall in the nasvar ±10% merge window."""
    if observed == target:
        return True
    if observed <= 0 or target <= 0:
        return False
    a = float(observed)
    b = float(target)
    return (0.9 * b < a < 1.1 * b) and (0.9 * a < b < 1.1 * a)


def collect_supporting_read_qc_for_event(
    bam: pysam.AlignmentFile,
    *,
    chrom: str,
    position: int,
    length: int,
    gene: str = "",
    bam_path: str = "",
) -> List[Dict[str, Any]]:
    """Collect per-read QC metrics for insertions supporting one called event."""
    bam_chrom = resolve_bam_contig(bam, chrom)
    if bam_chrom is None:
        return []

    rows: List[Dict[str, Any]] = []
    # Tiny window around the anchor.
    fetch_start = max(0, int(position) - 5)
    fetch_end = int(position) + 5
    for read in bam.fetch(bam_chrom, fetch_start, fetch_end + 1):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.cigartuples is None or read.reference_start is None:
            continue

        ref_pos = int(read.reference_start)
        matched_len: Optional[int] = None
        for op, op_len in read.cigartuples:
            if op in (0, 7, 8):
                ref_pos += op_len
            elif op == 1:
                if ref_pos == int(position) and _length_compatible(op_len, int(length)):
                    matched_len = int(op_len)
                    break
            elif op in (2, 3, 6):
                ref_pos += op_len
        if matched_len is None:
            continue

        soft_l, soft_r = _softclip_bases(read)
        query_len = int(read.query_length or 0)
        soft_total = soft_l + soft_r
        soft_frac = (soft_total / query_len) if query_len else 0.0
        ins_quals = _insertion_base_qualities(read, int(position), matched_len)
        mean_ins_q = (
            float(sum(ins_quals) / len(ins_quals)) if ins_quals else float("nan")
        )
        min_ins_q = float(min(ins_quals)) if ins_quals else float("nan")
        nm = read.get_tag("NM") if read.has_tag("NM") else None
        aln_len = int(read.reference_length or 0)
        rows.append(
            {
                "gene": gene,
                "chrom": chrom,
                "position": int(position),
                "event_length": int(length),
                "observed_length": matched_len,
                "read_name": read.query_name,
                "bam_path": os.path.basename(bam_path) if bam_path else "",
                "mapq": int(read.mapping_quality),
                "strand": "-" if read.is_reverse else "+",
                "query_length": query_len,
                "aligned_length": aln_len,
                "softclip_left": soft_l,
                "softclip_right": soft_r,
                "softclip_frac": round(soft_frac, 4),
                "mean_ins_baseq": None if mean_ins_q != mean_ins_q else round(mean_ins_q, 2),
                "min_ins_baseq": None if min_ins_q != min_ins_q else int(min_ins_q),
                "nm": int(nm) if nm is not None else None,
            }
        )
    return rows


def collect_itd_read_qc(
    events: pd.DataFrame,
    bam_paths: Sequence[str],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build read-level and event-level QC tables for called ITD events.

    Returns:
        read_qc: one row per supporting read
        event_qc: aggregated metrics per called event
    """
    if events is None or events.empty or not bam_paths:
        empty_read = pd.DataFrame()
        empty_event = pd.DataFrame()
        return empty_read, empty_event

    read_rows: List[Dict[str, Any]] = []
    for bam_path in bam_paths:
        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                for event in events.itertuples(index=False):
                    read_rows.extend(
                        collect_supporting_read_qc_for_event(
                            bam,
                            chrom=str(event.chrom),
                            position=int(event.position),
                            length=int(event.length),
                            gene=str(getattr(event, "gene", "")),
                            bam_path=str(bam_path),
                        )
                    )
        except Exception as exc:
            logger.warning("ITD QC: could not scan %s: %s", bam_path, exc)

    read_qc = pd.DataFrame(read_rows)
    if read_qc.empty:
        return read_qc, pd.DataFrame()

    # Deduplicate same read seen in multiple overlapping BAM slices (rare).
    read_qc = read_qc.drop_duplicates(
        subset=["gene", "chrom", "position", "event_length", "read_name"],
        keep="first",
    )

    agg_rows: List[Dict[str, Any]] = []
    group_cols = ["gene", "chrom", "position", "event_length"]
    for keys, group in read_qc.groupby(group_cols, sort=False):
        gene, chrom, position, event_length = keys
        n = len(group)
        fwd = int((group["strand"] == "+").sum())
        agg_rows.append(
            {
                "gene": gene,
                "chrom": chrom,
                "position": int(position),
                "length": int(event_length),
                "n_qc_reads": n,
                "mapq_median": float(group["mapq"].median()),
                "mapq_min": int(group["mapq"].min()),
                "mapq_lt20_frac": float((group["mapq"] < 20).mean()),
                "mapq_lt40_frac": float((group["mapq"] < 40).mean()),
                "strand_forward_frac": float(fwd / n) if n else 0.0,
                "query_length_median": float(group["query_length"].median()),
                "aligned_length_median": float(group["aligned_length"].median()),
                "softclip_frac_median": float(group["softclip_frac"].median()),
                "mean_ins_baseq_median": (
                    float(group["mean_ins_baseq"].dropna().median())
                    if group["mean_ins_baseq"].notna().any()
                    else None
                ),
                "min_ins_baseq_median": (
                    float(group["min_ins_baseq"].dropna().median())
                    if group["min_ins_baseq"].notna().any()
                    else None
                ),
                "nm_median": (
                    float(group["nm"].dropna().median())
                    if group["nm"].notna().any()
                    else None
                ),
            }
        )
    event_qc = pd.DataFrame(agg_rows)
    if not event_qc.empty and not events.empty:
        event_qc = event_qc.merge(
            events[
                [
                    c
                    for c in (
                        "gene",
                        "chrom",
                        "position",
                        "length",
                        "support",
                        "vaf",
                        "exon_number",
                        "transcript_id",
                    )
                    if c in events.columns
                ]
            ],
            on=["gene", "chrom", "position", "length"],
            how="left",
        )
    return read_qc, event_qc


def write_itd_read_qc(
    sample_dir: str | Path,
    events: pd.DataFrame,
    *,
    bam_paths: Optional[Sequence[str]] = None,
) -> Dict[str, str]:
    """Write ``itd_read_qc.csv`` and ``itd_event_qc.csv`` under ``sample_dir``."""
    sample_dir = Path(sample_dir)
    if bam_paths is None:
        bam_paths = sorted(
            str(path)
            for path in sample_dir.glob("batch_*.bam")
            if path.is_file()
        )
    read_qc, event_qc = collect_itd_read_qc(events, list(bam_paths or []))
    read_path = sample_dir / "itd_read_qc.csv"
    event_path = sample_dir / "itd_event_qc.csv"
    read_qc.to_csv(read_path, index=False)
    event_qc.to_csv(event_path, index=False)
    logger.info(
        "ITD QC: %d supporting reads → %s; %d event summaries → %s",
        len(read_qc),
        read_path,
        len(event_qc),
        event_path,
    )
    return {"read_qc_path": str(read_path), "event_qc_path": str(event_path)}


def call_itds_for_bam(
    bam_path: str,
    *,
    target_panel: str,
    hotspots_path: Optional[str | Path] = None,
    region_mode: str = "hotspots",
) -> pd.DataFrame:
    """One-shot helper: resolve ITD windows and call events on a single BAM."""
    hotspots = resolve_itd_hotspots(
        target_panel, region_mode=region_mode, hotspots_path=hotspots_path
    )
    if not hotspots:
        return pd.DataFrame(columns=ITD_EVENT_COLUMNS)
    counts = process_bam_itd_pass(bam_path, hotspots)
    return _events_from_aggregated_counts(counts, hotspots)


def clear_itd_staging(work_dir: str, sample_id: str) -> None:
    """Remove staging directory contents (tests / recovery)."""
    staging = Path(_staging_dir(work_dir, sample_id))
    if staging.is_dir():
        shutil.rmtree(staging, ignore_errors=True)
