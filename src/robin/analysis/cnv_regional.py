"""Shared regional CNV analysis helpers for reports and GUI."""

from __future__ import annotations

import logging
import os
import re

import numpy as np
import pandas as pd

from robin import resources
from robin.analysis.cnv_classification import CNVEvent
from robin.utils.sequencing_files import panel_bed_filename

logger = logging.getLogger(__name__)

from robin.reference_contigs import CANONICAL_CONTIG_RE, is_canonical_contig

REPORTABLE_CHROMOSOME_RE = CANONICAL_CONTIG_RE
SIGNIFICANT_CNV_STATES = {"GAIN", "LOSS", "HIGH_GAIN", "DEEP_LOSS"}


def is_reportable_chromosome(chromosome: str) -> bool:
    """Return True for standard autosomes/sex chromosomes used in CNV reporting."""
    return is_canonical_contig(chromosome)


def load_panel_gene_bed(output_dir: str) -> tuple[str | None, pd.DataFrame]:
    """Load the target panel gene BED for the sample's analysis panel."""
    panel = None
    master_csv = os.path.join(output_dir, "master.csv")
    if os.path.exists(master_csv):
        try:
            master_df = pd.read_csv(master_csv)
            if not master_df.empty and "analysis_panel" in master_df.columns:
                panel_val = master_df.iloc[0]["analysis_panel"]
                if panel_val is not None and str(panel_val).strip().lower() not in ("", "nan"):
                    panel = str(panel_val).strip()
        except Exception as exc:
            logger.debug("Could not read analysis panel from master.csv: %s", exc)

    empty = pd.DataFrame(columns=["chrom", "start_pos", "end_pos", "gene"])
    if not panel:
        return None, empty

    bed_path = os.path.join(
        os.path.dirname(os.path.abspath(resources.__file__)),
        panel_bed_filename(panel),
    )
    if not os.path.exists(bed_path):
        logger.warning("Target panel BED not found for panel '%s' at %s", panel, bed_path)
        return panel, empty

    return panel, pd.read_csv(
        bed_path,
        sep="\t",
        header=None,
        names=["chrom", "start_pos", "end_pos", "gene"],
    )


def load_target_coverage_df(output_dir: str) -> pd.DataFrame:
    """Load per-target coverage for panel gene coverage overlays."""
    path = os.path.join(output_dir, "target_coverage.csv")
    empty = pd.DataFrame(
        columns=["chrom", "startpos", "endpos", "name", "length", "coverage", "bases"]
    )
    if not os.path.exists(path):
        return empty
    try:
        return pd.read_csv(path)
    except Exception as exc:
        logger.debug("Could not read target coverage from %s: %s", path, exc)
        return empty


def format_chromosome_cnv_status(
    chromosome: str,
    events: list[CNVEvent],
    cytoband_analysis: pd.DataFrame,
) -> str:
    """Summarize detected CNV changes for an individual chromosome plot caption."""
    chrom_events = [event for event in events if event.chromosome == chromosome]
    if chrom_events:
        parts = []
        for event in chrom_events:
            if event.event_type.startswith("WHOLE_CHR_"):
                parts.append(
                    f"Whole chromosome {event.event_type.replace('WHOLE_CHR_', '')}"
                )
            elif event.arm:
                parts.append(f"{event.arm}-arm {event.event_type}")
            else:
                parts.append(event.event_type)
        return "; ".join(parts)

    if not cytoband_analysis.empty:
        regional = [
            f"{row['cnv_state']}: {row['name']}"
            for _, row in cytoband_analysis.iterrows()
            if row["cnv_state"] in SIGNIFICANT_CNV_STATES
        ]
        if regional:
            if len(regional) > 3:
                return "; ".join(regional[:3]) + f"; +{len(regional) - 3} more"
            return "; ".join(regional)

    return "No significant CNV change"


def build_significant_regions(cytoband_analysis: pd.DataFrame) -> list[dict]:
    """Convert cytoband analysis rows into plot highlight regions."""
    regions = []
    if cytoband_analysis.empty:
        return regions

    for _, row in cytoband_analysis.iterrows():
        if row["cnv_state"] not in SIGNIFICANT_CNV_STATES:
            continue
        regions.append(
            {
                "start_pos": int(row["start_pos"]),
                "end_pos": int(row["end_pos"]),
                "type": row["cnv_state"],
                "name": row.get("name", ""),
            }
        )
    return regions


def panel_genes_in_region(
    panel_genes_df: pd.DataFrame,
    chrom: str,
    start_pos: int,
    end_pos: int,
) -> list[str]:
    """Return sorted unique panel gene labels overlapping a genomic interval."""
    if panel_genes_df.empty:
        return []

    hits = panel_genes_df[
        (panel_genes_df["chrom"] == chrom)
        & (panel_genes_df["start_pos"] <= end_pos)
        & (panel_genes_df["end_pos"] >= start_pos)
    ]
    labels: list[str] = []
    seen: set[str] = set()
    for gene_value in hits["gene"].astype(str):
        for part in str(gene_value).split(","):
            label = part.strip()
            if label and label.lower() != "nan" and label not in seen:
                seen.add(label)
                labels.append(label)
    return sorted(labels)


def format_panel_genes_for_table(panel_genes: list[str]) -> str:
    """Format panel gene labels for report tables."""
    return ", ".join(panel_genes) if panel_genes else "—"


def build_regional_cnv_events(
    cytoband_analysis: pd.DataFrame,
    panel_genes_df: pd.DataFrame,
) -> list[dict]:
    """Build focal cytoband CNV events with overlapping panel target genes."""
    regional_events: list[dict] = []
    if cytoband_analysis.empty:
        return regional_events

    for _, row in cytoband_analysis.iterrows():
        if row["cnv_state"] not in SIGNIFICANT_CNV_STATES:
            continue

        chrom = str(row["chrom"])
        start_pos = int(row["start_pos"])
        end_pos = int(row["end_pos"])
        regional_events.append(
            {
                "chrom": chrom.replace("chr", ""),
                "chromosome": chrom,
                "region": str(row.get("name", "")),
                "start_pos": start_pos,
                "end_pos": end_pos,
                "start_mb": start_pos / 1_000_000,
                "end_mb": end_pos / 1_000_000,
                "length_mb": (end_pos - start_pos) / 1_000_000,
                "mean_cnv": float(row["mean_cnv"]),
                "state": str(row["cnv_state"]),
                "panel_genes": panel_genes_in_region(
                    panel_genes_df, chrom, start_pos, end_pos,
                ),
            }
        )

    return regional_events


def format_regional_event_table_row(event: dict) -> dict:
    """Format a regional CNV event dict for NiceGUI table rows."""
    return {
        "chrom": event["chrom"],
        "region": event["region"],
        "start_mb": f"{event['start_mb']:.2f}",
        "end_mb": f"{event['end_mb']:.2f}",
        "length_mb": f"{event['length_mb']:.2f}",
        "mean_cnv": f"{event['mean_cnv']:.3f}",
        "state": event["state"],
        "panel_genes": format_panel_genes_for_table(event["panel_genes"]),
    }


def analyze_cytoband_cnv(
    cnv_data: dict,
    chromosome: str,
    cnv_dict: dict,
    cytobands_bed: pd.DataFrame,
    centromere_bed: pd.DataFrame,
    gene_bed: pd.DataFrame,
    sex_estimate: str,
) -> pd.DataFrame:
    """
    Analyze CNV values within each cytoband to detect duplications and deletions.

    Expects ``cnv_data`` on the log2(ploidy / expected) scale (0 = normal).
    Gain/loss cut-offs are the same fixed thresholds as arm / whole-chromosome
    calling (``get_cnv_thresholds``, default ±0.30). Non-finite bins (no
    coverage) are excluded from band means; bands without usable data are
    marked ``NO_DATA``.
    """
    from robin.classification_config import get_cnv_thresholds

    logger.debug(f"\n{'='*50}")
    logger.debug(f"Starting CNV analysis for {chromosome}")
    logger.debug(f"CNV data keys: {list(cnv_data.keys())}")

    if "bin_width" not in cnv_dict:
        logger.debug("No cnv_dict or bin_width available")
        return pd.DataFrame()

    logger.debug(f"Bin width: {cnv_dict['bin_width']}")

    if cnv_dict["bin_width"] > 10_000_000:
        logger.debug("Resolution insufficient for CNV calling")
        return pd.DataFrame()

    if chromosome not in cnv_data:
        return pd.DataFrame()

    bin_width = int(cnv_dict["bin_width"])
    chromosome_cytobands = cytobands_bed[cytobands_bed["chrom"] == chromosome].copy()
    logger.debug(f"Number of cytobands for {chromosome}: {len(chromosome_cytobands)}")
    if chromosome_cytobands.empty:
        return pd.DataFrame()

    chrom_arr = np.asarray(cnv_data[chromosome], dtype=float)
    if not np.any(np.isfinite(chrom_arr)):
        logger.debug(f"No finite CNV bins for {chromosome}")
        return pd.DataFrame()

    # Same absolute cut-offs as arm / whole-chromosome event detection.
    cytoband_gain_threshold, cytoband_loss_threshold = get_cnv_thresholds(
        chromosome, sex_estimate
    )
    logger.debug(
        "Thresholds - Cytoband gain: %.3f, loss: %.3f (arm-rule cut-offs)",
        cytoband_gain_threshold,
        cytoband_loss_threshold,
    )

    centromere = centromere_bed[centromere_bed["chrom"] == chromosome]
    cen_start = int(centromere["start_pos"].iloc[0]) if not centromere.empty else None
    cen_end = int(centromere["end_pos"].iloc[0]) if not centromere.empty else None

    max_expected_cytobands = len(chromosome_cytobands)
    merged_cytobands = [None] * max_expected_cytobands
    merged_idx = 0
    current_group = None

    for _, cytoband in chromosome_cytobands.iterrows():
        start_bin = int(cytoband["start_pos"] / bin_width)
        end_bin = int(cytoband["end_pos"] / bin_width)
        start_bin = max(0, start_bin)
        end_bin = min(chrom_arr.size - 1, end_bin)

        # Centromeric / satellite bands are not assessed (same as prior regional logic).
        band_overlaps_centromere = (
            cen_start is not None
            and cen_end is not None
            and int(cytoband["start_pos"]) < cen_end
            and int(cytoband["end_pos"]) > cen_start
        )

        if band_overlaps_centromere:
            mean_cnv = float("nan")
            state = "NO_DATA"
        elif start_bin < chrom_arr.size and end_bin >= start_bin:
            region_cnv = chrom_arr[start_bin : end_bin + 1]
            finite = region_cnv[np.isfinite(region_cnv)]
            if finite.size == 0:
                mean_cnv = float("nan")
                state = "NO_DATA"
            else:
                mean_cnv = float(np.mean(finite))
                if mean_cnv >= cytoband_gain_threshold:
                    state = "GAIN"
                elif mean_cnv <= cytoband_loss_threshold:
                    state = "LOSS"
                else:
                    state = "NORMAL"
        else:
            mean_cnv = float("nan")
            state = "NO_DATA"

        if current_group is None:
            current_group = {
                "chrom": cytoband["chrom"],
                "start_pos": cytoband["start_pos"],
                "end_pos": cytoband["end_pos"],
                "name": cytoband["name"],
                "mean_cnv": [mean_cnv],
                "cnv_state": state,
                "bands": [cytoband["name"]],
                "length": cytoband["end_pos"] - cytoband["start_pos"],
                "genes": [],
            }
        elif state == current_group["cnv_state"]:
            current_group["end_pos"] = cytoband["end_pos"]
            current_group["mean_cnv"].append(mean_cnv)
            current_group["bands"].append(cytoband["name"])
        else:
            if current_group["cnv_state"] in SIGNIFICANT_CNV_STATES:
                genes_in_region = gene_bed[
                    (gene_bed["chrom"] == current_group["chrom"])
                    & (gene_bed["start_pos"] <= current_group["end_pos"])
                    & (gene_bed["end_pos"] >= current_group["start_pos"])
                ]["gene"].tolist()
                current_group["genes"] = genes_in_region

            current_group["name"] = (
                f"{current_group['chrom']} {current_group['bands'][0]}-{current_group['bands'][-1]}"
            )
            finite_means = [
                v for v in current_group["mean_cnv"] if np.isfinite(v)
            ]
            current_group["mean_cnv"] = (
                float(np.mean(finite_means)) if finite_means else float("nan")
            )
            current_group["length"] = (
                current_group["end_pos"] - current_group["start_pos"]
            )

            if current_group["cnv_state"] not in ("NORMAL", "NO_DATA"):
                merged_cytobands[merged_idx] = current_group
                merged_idx += 1

            current_group = {
                "chrom": cytoband["chrom"],
                "start_pos": cytoband["start_pos"],
                "end_pos": cytoband["end_pos"],
                "name": cytoband["name"],
                "mean_cnv": [mean_cnv],
                "cnv_state": state,
                "bands": [cytoband["name"]],
                "length": cytoband["end_pos"] - cytoband["start_pos"],
                "genes": [],
            }

    if current_group is not None:
        if current_group["cnv_state"] in SIGNIFICANT_CNV_STATES:
            genes_in_region = gene_bed[
                (gene_bed["chrom"] == current_group["chrom"])
                & (gene_bed["start_pos"] <= current_group["end_pos"])
                & (gene_bed["end_pos"] >= current_group["start_pos"])
            ]["gene"].tolist()
            current_group["genes"] = genes_in_region

        current_group["name"] = (
            f"{current_group['chrom']} {current_group['bands'][0]}-{current_group['bands'][-1]}"
        )
        finite_means = [v for v in current_group["mean_cnv"] if np.isfinite(v)]
        current_group["mean_cnv"] = (
            float(np.mean(finite_means)) if finite_means else float("nan")
        )
        current_group["length"] = (
            current_group["end_pos"] - current_group["start_pos"]
        )

        if current_group["cnv_state"] not in ("NORMAL", "NO_DATA"):
            merged_cytobands[merged_idx] = current_group
            merged_idx += 1

    merged_cytobands = merged_cytobands[:merged_idx]

    merged_df = pd.DataFrame(merged_cytobands)
    if not merged_df.empty:
        merged_df = merged_df.sort_values("start_pos")

    return merged_df
