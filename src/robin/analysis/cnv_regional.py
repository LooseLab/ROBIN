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

REPORTABLE_CHROMOSOME_RE = re.compile(r"^chr(\d+|X|Y)$")
SIGNIFICANT_CNV_STATES = {"GAIN", "LOSS", "HIGH_GAIN", "DEEP_LOSS"}


def is_reportable_chromosome(chromosome: str) -> bool:
    """Return True for standard autosomes/sex chromosomes used in CNV reporting."""
    return chromosome != "chrM" and bool(REPORTABLE_CHROMOSOME_RE.match(chromosome))


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
    Uses dynamic thresholds based on data variation for more robust detection.
    """
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

    bin_width = cnv_dict["bin_width"]
    chromosome_cytobands = cytobands_bed[cytobands_bed["chrom"] == chromosome].copy()
    logger.debug(f"Number of cytobands for {chromosome}: {len(chromosome_cytobands)}")

    max_expected_cytobands = len(chromosome_cytobands)
    merged_cytobands = [None] * max_expected_cytobands
    merged_idx = 0
    whole_chr_event = False
    whole_chr_state = "NORMAL"

    if chromosome in cnv_data:
        logger.debug(
            f"\nAnalyzing chromosome {chromosome} for whole chromosome events:"
        )

        mask = np.ones(len(cnv_data[chromosome]), dtype=bool)
        centromere = centromere_bed[centromere_bed["chrom"] == chromosome]
        if not centromere.empty:
            cent_start_bin = int(centromere["start_pos"].iloc[0] / bin_width)
            cent_end_bin = int(centromere["end_pos"].iloc[0] / bin_width)
            mask[cent_start_bin:cent_end_bin] = False
            logger.debug(f"Excluded centromere region: {cent_start_bin}-{cent_end_bin}")
        chr_cnv = cnv_data[chromosome][mask]

        chr_mean = np.mean(chr_cnv)
        chr_std = np.std(chr_cnv)
        logger.debug(f"Chromosome-wide mean: {chr_mean:.3f}, std: {chr_std:.3f}")

        chromosome_means = []
        for chrom in cnv_data:
            if chrom.startswith("chr") and chrom[3:].isdigit():
                mask = np.ones(len(cnv_data[chrom]), dtype=bool)
                cent = centromere_bed[centromere_bed["chrom"] == chrom]
                if not cent.empty:
                    cent_start = int(cent["start_pos"].iloc[0] / bin_width)
                    cent_end = int(cent["end_pos"].iloc[0] / bin_width)
                    mask[cent_start:cent_end] = False
                chrom_data = cnv_data[chrom][mask]
                if len(chrom_data) > 0:
                    chromosome_means.append(np.mean(chrom_data))

        means_std = np.std(chromosome_means)
        means_mean = np.mean(chromosome_means)
        logger.debug(
            f"Mean of chromosome means: {means_mean:.3f}, std of means: {means_std:.3f}"
        )

        if chromosome.startswith("chr") and chromosome[3:].isdigit():
            gain_threshold = means_mean + (1.0 * means_std)
            loss_threshold = means_mean - (1.0 * means_std)
            cytoband_gain_threshold = chr_mean + (1.0 * chr_std)
            cytoband_loss_threshold = chr_mean - (1.0 * chr_std)
        elif chromosome == "chrX":
            gain_threshold = means_mean + (1.0 * means_std)
            loss_threshold = means_mean - (1.0 * means_std)
            cytoband_gain_threshold = chr_mean + (1.0 * chr_std)
            cytoband_loss_threshold = chr_mean - (1.0 * chr_std)
        elif chromosome == "chrY":
            if sex_estimate in ("Male", "XY"):
                gain_threshold = means_mean + (1.0 * means_std)
                loss_threshold = means_mean - (1.0 * means_std)
                cytoband_gain_threshold = chr_mean + (1.0 * chr_std)
                cytoband_loss_threshold = chr_mean - (1.0 * chr_std)
            else:
                gain_threshold = means_mean + (1.2 * means_std)
                loss_threshold = means_mean - (1.2 * means_std)
                cytoband_gain_threshold = chr_mean + (1.2 * chr_std)
                cytoband_loss_threshold = chr_mean - (1.2 * chr_std)
        else:
            gain_threshold = means_mean + (1.0 * means_std)
            loss_threshold = means_mean - (1.0 * means_std)
            cytoband_gain_threshold = chr_mean + (1.0 * chr_std)
            cytoband_loss_threshold = chr_mean - (1.0 * chr_std)

        logger.debug(
            f"Thresholds - Whole chr gain: {gain_threshold:.3f}, loss: {loss_threshold:.3f}"
        )
        logger.debug(
            f"Thresholds - Cytoband gain: {cytoband_gain_threshold:.3f}, loss: {cytoband_loss_threshold:.3f}"
        )

        bins_above_gain = np.sum(chr_cnv > gain_threshold) / len(chr_cnv)
        bins_below_loss = np.sum(chr_cnv < loss_threshold) / len(chr_cnv)

        logger.debug(
            f"Proportion of bins - Above gain: {bins_above_gain:.3f}, Below loss: {bins_below_loss:.3f}"
        )

        min_proportion = 0.7
        if bins_above_gain > min_proportion:
            whole_chr_event = True
            whole_chr_state = "GAIN"
            logger.debug(f"WHOLE CHROMOSOME EVENT DETECTED: {chromosome} GAIN")
        elif bins_below_loss > min_proportion:
            whole_chr_event = True
            whole_chr_state = "LOSS"
            logger.debug(f"WHOLE CHROMOSOME EVENT DETECTED: {chromosome} LOSS")

        if whole_chr_event:
            genes_in_chr = gene_bed[gene_bed["chrom"] == chromosome]["gene"].tolist()

            merged_cytobands[merged_idx] = {
                "chrom": chromosome,
                "start_pos": chromosome_cytobands["start_pos"].min(),
                "end_pos": chromosome_cytobands["end_pos"].max(),
                "name": f"{chromosome} WHOLE CHROMOSOME {whole_chr_state}",
                "mean_cnv": chr_mean,
                "cnv_state": whole_chr_state,
                "length": chromosome_cytobands["end_pos"].max()
                - chromosome_cytobands["start_pos"].min(),
                "genes": genes_in_chr,
            }
            merged_idx += 1

        current_group = None

        for _, cytoband in chromosome_cytobands.iterrows():
            start_bin = int(cytoband["start_pos"] / bin_width)
            end_bin = int(cytoband["end_pos"] / bin_width)

            if start_bin < len(cnv_data[chromosome]):
                region_cnv = cnv_data[chromosome][start_bin : end_bin + 1]
                mean_cnv = np.mean(region_cnv) if len(region_cnv) > 0 else 0

                if mean_cnv > cytoband_gain_threshold:
                    state = "GAIN"
                elif mean_cnv < cytoband_loss_threshold:
                    state = "LOSS"
                else:
                    state = "NORMAL"
            else:
                mean_cnv = 0
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
                current_group["mean_cnv"] = np.mean(current_group["mean_cnv"])
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
            current_group["mean_cnv"] = np.mean(current_group["mean_cnv"])
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
