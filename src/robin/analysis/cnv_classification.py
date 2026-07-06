"""
Centralized CNV classification and analysis functions.

This module provides CNV analysis functions that use the centralized
classification rules from classification_config.py.
"""

import logging
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import natsort
from typing import Dict, List, Tuple, Optional, Any

try:
    from importlib import resources as importlib_resources
except ImportError:  # pragma: no cover
    import importlib_resources  # type: ignore

from robin.classification_config import (
    get_cnv_thresholds,
    is_whole_chromosome_event,
    is_arm_event,
    is_resolution_sufficient
)

logger = logging.getLogger(__name__)


class CNVEvent:
    """Represents a CNV event with metadata."""
    
    def __init__(
        self,
        chromosome: str,
        event_type: str,
        mean_cnv: float,
        start_pos: int,
        end_pos: int,
        length: int,
        genes: List[str] = None,
        confidence: str = "Unknown",
        arm: Optional[str] = None,
        proportion_affected: float = 0.0
    ):
        self.chromosome = chromosome
        self.event_type = event_type  # 'GAIN', 'LOSS', 'WHOLE_CHR_GAIN', 'WHOLE_CHR_LOSS'
        self.mean_cnv = mean_cnv
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.length = length
        self.genes = genes or []
        self.confidence = confidence
        self.arm = arm
        self.proportion_affected = proportion_affected
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for GUI/reporting."""
        return {
            "chromosome": self.chromosome,
            "event_type": self.event_type,
            "mean_cnv": self.mean_cnv,
            "start_pos": self.start_pos,
            "end_pos": self.end_pos,
            "length": self.length,
            "genes": self.genes,
            "confidence": self.confidence,
            "arm": self.arm,
            "proportion_affected": self.proportion_affected,
            "start_pos_mb": f"{self.start_pos/1e6:.2f}",
            "end_pos_mb": f"{self.end_pos/1e6:.2f}",
            "length_mb": f"{self.length/1e6:.2f}",
            "mean_cnv_str": f"{self.mean_cnv:.3f}",
            "genes_str": ", ".join(self.genes) if self.genes else "",
        }


def _finite_arm_values(values: List[float]) -> List[float]:
    return [float(v) for v in values if np.isfinite(v)]


def _arm_bin_proportions(
    values: List[float], gain_threshold: float, loss_threshold: float
) -> Tuple[float, float]:
    if not values:
        return 0.0, 0.0
    n = len(values)
    gain_prop = sum(1 for v in values if v > gain_threshold) / n
    loss_prop = sum(1 for v in values if v < loss_threshold) / n
    return gain_prop, loss_prop


def analyze_chromosome_arms(
    cnv_data: Dict[str, np.ndarray],
    chromosome: str,
    bin_width: int,
    sex_estimate: str,
    cytobands_df: pd.DataFrame
) -> Tuple[
    Optional[float],
    Optional[float],
    float,
    float,
    float,
    float,
]:
    """
    Analyze p and q arms of a chromosome for CNV events.

    Returns:
        Tuple of (
            p_arm_mean, q_arm_mean,
            p_arm_proportion_gain, p_arm_proportion_loss,
            q_arm_proportion_gain, q_arm_proportion_loss,
        )
    """
    if chromosome not in cnv_data:
        return None, None, 0.0, 0.0, 0.0, 0.0

    gain_threshold, loss_threshold = get_cnv_thresholds(chromosome, sex_estimate)

    chr_cytobands = cytobands_df[cytobands_df["chrom"] == chromosome]
    if chr_cytobands.empty:
        return None, None, 0.0, 0.0, 0.0, 0.0

    logger.debug(f"{chromosome} cytoband names: {chr_cytobands['name'].tolist()}")

    p_arm_cytobands = chr_cytobands[
        chr_cytobands["name"].str.startswith("p", na=False)
    ]

    p_arm_mean = None
    p_arm_proportion_gain = 0.0
    p_arm_proportion_loss = 0.0

    if not p_arm_cytobands.empty:
        p_arm_values: List[float] = []
        for _, band in p_arm_cytobands.iterrows():
            start_pos_bin = max(0, int(band["start_pos"] // bin_width))
            end_pos_bin = min(len(cnv_data[chromosome]) - 1, int(band["end_pos"] // bin_width))
            if end_pos_bin >= start_pos_bin:
                region_values = cnv_data[chromosome][start_pos_bin:end_pos_bin + 1]
                p_arm_values.extend(region_values)
        p_arm_values = _finite_arm_values(p_arm_values)
        if p_arm_values:
            p_arm_mean = float(np.mean(p_arm_values))
            p_arm_proportion_gain, p_arm_proportion_loss = _arm_bin_proportions(
                p_arm_values, gain_threshold, loss_threshold
            )

    logger.debug(
        f"{chromosome} p-arm: {len(p_arm_cytobands)} bands found, mean={p_arm_mean}, "
        f"gain_prop={p_arm_proportion_gain:.3f}, loss_prop={p_arm_proportion_loss:.3f}"
    )

    q_arm_cytobands = chr_cytobands[
        chr_cytobands["name"].str.startswith("q", na=False)
    ]

    q_arm_mean = None
    q_arm_proportion_gain = 0.0
    q_arm_proportion_loss = 0.0

    if not q_arm_cytobands.empty:
        q_arm_values: List[float] = []
        for _, band in q_arm_cytobands.iterrows():
            start_pos_bin = max(0, int(band["start_pos"] // bin_width))
            end_pos_bin = min(len(cnv_data[chromosome]) - 1, int(band["end_pos"] // bin_width))
            if end_pos_bin >= start_pos_bin:
                region_values = cnv_data[chromosome][start_pos_bin:end_pos_bin + 1]
                q_arm_values.extend(region_values)
        q_arm_values = _finite_arm_values(q_arm_values)
        if q_arm_values:
            q_arm_mean = float(np.mean(q_arm_values))
            q_arm_proportion_gain, q_arm_proportion_loss = _arm_bin_proportions(
                q_arm_values, gain_threshold, loss_threshold
            )

    logger.debug(
        f"{chromosome} q-arm: {len(q_arm_cytobands)} bands found, mean={q_arm_mean}, "
        f"gain_prop={q_arm_proportion_gain:.3f}, loss_prop={q_arm_proportion_loss:.3f}"
    )

    return (
        p_arm_mean,
        q_arm_mean,
        p_arm_proportion_gain,
        p_arm_proportion_loss,
        q_arm_proportion_gain,
        q_arm_proportion_loss,
    )


def detect_cnv_events(
    cnv_data: Dict[str, np.ndarray],
    bin_width: int,
    sex_estimate: str,
    cytobands_df: pd.DataFrame,
    gene_df: Optional[pd.DataFrame] = None
) -> List[CNVEvent]:
    """
    Detect CNV events using centralized classification rules.

    Expects ``cnv_data`` on the log2(ploidy / expected copy number) scale —
    use ``prepare_cnv_calling_track()`` from ``cnv_analysis`` to build the
    track from absolute ploidy (``CNV.npy``), coarsened to at least 1 Mb bins.

    Args:
        cnv_data: Per-chromosome log2 ratio arrays
        bin_width: Bin width in base pairs
        sex_estimate: Sex estimate
        cytobands_df: Cytobands dataframe
        gene_df: Optional gene dataframe

    Returns:
        List of CNVEvent objects
    """
    events = []
    
    logger.debug(f"Detecting CNV events with sex_estimate='{sex_estimate}', bin_width={bin_width}")
    
    # Check if resolution is sufficient
    if not is_resolution_sufficient(bin_width):
        logger.warning(f"Resolution insufficient for CNV calling: bin_width={bin_width}")
        return events
    
    # Analyze each chromosome
    logger.debug(f"Available chromosomes: {list(cnv_data.keys())}")
    for chromosome in natsort.natsorted(cnv_data.keys()):
        if chromosome == "chrM" or not chromosome.startswith("chr"):
            continue
        
        # Skip Y chromosome for male samples (expected absence)
        if chromosome == "chrY" and sex_estimate.upper() in ("XY", "MALE"):
            logger.debug(f"Skipping {chromosome} for male sample")
            continue
        
        logger.debug(f"Analyzing chromosome {chromosome} for CNV events")
        
        # Get thresholds
        gain_threshold, loss_threshold = get_cnv_thresholds(chromosome, sex_estimate)
        
        (
            p_arm_mean,
            q_arm_mean,
            p_arm_proportion_gain,
            p_arm_proportion_loss,
            q_arm_proportion_gain,
            q_arm_proportion_loss,
        ) = analyze_chromosome_arms(
            cnv_data, chromosome, bin_width, sex_estimate, cytobands_df
        )
        
        # Check for whole chromosome events
        if p_arm_mean is not None and q_arm_mean is not None:
            logger.debug(
                f"{chromosome}: p_mean={p_arm_mean:.3f}, q_mean={q_arm_mean:.3f}, "
                f"p_gain={p_arm_proportion_gain:.3f}, p_loss={p_arm_proportion_loss:.3f}, "
                f"q_gain={q_arm_proportion_gain:.3f}, q_loss={q_arm_proportion_loss:.3f}"
            )
            is_whole_chr, event_type = is_whole_chromosome_event(
                p_arm_mean,
                q_arm_mean,
                p_arm_proportion_gain,
                p_arm_proportion_loss,
                q_arm_proportion_gain,
                q_arm_proportion_loss,
                gain_threshold,
                loss_threshold,
            )
            
            if is_whole_chr:
                # Create whole chromosome event
                chr_cytobands = cytobands_df[cytobands_df["chrom"] == chromosome]
                if not chr_cytobands.empty:
                    start_pos = int(chr_cytobands["start_pos"].min())
                    end_pos = int(chr_cytobands["end_pos"].max())
                    length = end_pos - start_pos
                    
                    # Get genes in chromosome
                    genes = []
                    if gene_df is not None:
                        genes = gene_df[gene_df["chrom"] == chromosome]["gene"].astype(str).tolist()
                    
                    chr_vals = cnv_data[chromosome]
                    chr_mean = float(np.nanmean(np.asarray(chr_vals, dtype=float)))
                    if event_type == "GAIN":
                        arm_prop = max(p_arm_proportion_gain, q_arm_proportion_gain)
                    else:
                        arm_prop = max(p_arm_proportion_loss, q_arm_proportion_loss)
                    
                    event = CNVEvent(
                        chromosome=chromosome,
                        event_type=f"WHOLE_CHR_{event_type}",
                        mean_cnv=chr_mean,
                        start_pos=start_pos,
                        end_pos=end_pos,
                        length=length,
                        genes=genes,
                        confidence="High" if arm_prop > 0.8 else "Medium",
                        proportion_affected=arm_prop,
                    )
                    events.append(event)
                    logger.info(f"Detected whole chromosome {event_type} for {chromosome}")

            elif chromosome != "chrY":
                # Arm-level events only when no whole-chromosome call on this chromosome
                if p_arm_mean is not None:
                    is_p_event, p_event_type = is_arm_event(
                        p_arm_mean,
                        p_arm_proportion_gain,
                        p_arm_proportion_loss,
                        gain_threshold,
                        loss_threshold,
                    )
                    if is_p_event:
                        chr_cytobands = cytobands_df[cytobands_df["chrom"] == chromosome]
                        p_bands = chr_cytobands[chr_cytobands["name"].str.startswith("p", na=False)]
                        if not p_bands.empty:
                            start_pos = int(p_bands["start_pos"].min())
                            end_pos = int(p_bands["end_pos"].max())
                            length = end_pos - start_pos
                            
                            # Get genes in p arm
                            genes = []
                            if gene_df is not None:
                                genes = gene_df[
                                    (gene_df["chrom"] == chromosome) &
                                    (gene_df["start_pos"] <= end_pos) &
                                    (gene_df["end_pos"] >= start_pos)
                                ]["gene"].astype(str).tolist()
                            
                            p_prop = (
                                p_arm_proportion_gain
                                if p_event_type == "GAIN"
                                else p_arm_proportion_loss
                            )
                            event = CNVEvent(
                                chromosome=chromosome,
                                event_type=p_event_type,
                                mean_cnv=p_arm_mean,
                                start_pos=start_pos,
                                end_pos=end_pos,
                                length=length,
                                genes=genes,
                                confidence="High" if p_prop > 0.8 else "Medium",
                                arm="p",
                                proportion_affected=p_prop,
                            )
                            events.append(event)
                            logger.info(f"Detected p-arm {p_event_type} for {chromosome}")
                
                # Check q arm
                if q_arm_mean is not None:
                    is_q_event, q_event_type = is_arm_event(
                        q_arm_mean,
                        q_arm_proportion_gain,
                        q_arm_proportion_loss,
                        gain_threshold,
                        loss_threshold,
                    )
                    if is_q_event:
                        chr_cytobands = cytobands_df[cytobands_df["chrom"] == chromosome]
                        q_bands = chr_cytobands[chr_cytobands["name"].str.startswith("q", na=False)]
                        if not q_bands.empty:
                            start_pos = int(q_bands["start_pos"].min())
                            end_pos = int(q_bands["end_pos"].max())
                            length = end_pos - start_pos
                            
                            # Get genes in q arm
                            genes = []
                            if gene_df is not None:
                                genes = gene_df[
                                    (gene_df["chrom"] == chromosome) &
                                    (gene_df["start_pos"] <= end_pos) &
                                    (gene_df["end_pos"] >= start_pos)
                                ]["gene"].astype(str).tolist()
                            
                            q_prop = (
                                q_arm_proportion_gain
                                if q_event_type == "GAIN"
                                else q_arm_proportion_loss
                            )
                            event = CNVEvent(
                                chromosome=chromosome,
                                event_type=q_event_type,
                                mean_cnv=q_arm_mean,
                                start_pos=start_pos,
                                end_pos=end_pos,
                                length=length,
                                genes=genes,
                                confidence="High" if q_prop > 0.8 else "Medium",
                                arm="q",
                                proportion_affected=q_prop,
                            )
                            events.append(event)
                            logger.info(f"Detected q-arm {q_event_type} for {chromosome}")
        else:
            # Single arm chromosome - use stricter threshold
            # Only use this logic if we truly have only one arm (like chrY in some cases)
            # For chromosomes that should have both arms, this indicates a cytoband parsing issue
            logger.warning(f"{chromosome}: Only one arm detected - this may indicate a cytoband parsing issue")
            
            whole_chr_mean = float(np.mean(cnv_data[chromosome]))
            single_arm_multiplier = 1.5  # From CNV_EVENT_RULES
            
            if abs(whole_chr_mean) > abs(gain_threshold) * single_arm_multiplier:
                chr_cytobands = cytobands_df[cytobands_df["chrom"] == chromosome]
                if not chr_cytobands.empty:
                    start_pos = int(chr_cytobands["start_pos"].min())
                    end_pos = int(chr_cytobands["end_pos"].max())
                    length = end_pos - start_pos
                    
                    # Get genes in chromosome
                    genes = []
                    if gene_df is not None:
                        genes = gene_df[gene_df["chrom"] == chromosome]["gene"].astype(str).tolist()
                    
                    event_type = "GAIN" if whole_chr_mean > gain_threshold else "LOSS"
                    event = CNVEvent(
                        chromosome=chromosome,
                        event_type=f"WHOLE_CHR_{event_type}",
                        mean_cnv=whole_chr_mean,
                        start_pos=start_pos,
                        end_pos=end_pos,
                        length=length,
                        genes=genes,
                        confidence="Medium",  # Single arm events are less certain
                        proportion_affected=1.0  # Entire chromosome
                    )
                    events.append(event)
                    logger.info(f"Detected single-arm whole chromosome {event_type} for {chromosome}")
    
    return events


def get_cnv_summary(events: List[CNVEvent]) -> Dict[str, Any]:
    """
    Generate a summary of CNV events.
    
    Args:
        events: List of CNVEvent objects
    
    Returns:
        Dictionary with summary statistics
    """
    summary = {
        "total_events": len(events),
        "whole_chromosome_events": [],
        "arm_events": [],
        "gene_containing_events": [],
        "total_genes_affected": set(),
        "high_confidence_events": 0,
        "medium_confidence_events": 0,
    }
    
    for event in events:
        if event.event_type.startswith("WHOLE_CHR_"):
            summary["whole_chromosome_events"].append(event)
        else:
            summary["arm_events"].append(event)
        
        if event.genes:
            summary["gene_containing_events"].append(event)
            summary["total_genes_affected"].update(event.genes)
        
        if event.confidence == "High":
            summary["high_confidence_events"] += 1
        elif event.confidence == "Medium":
            summary["medium_confidence_events"] += 1
    
    summary["total_genes_affected"] = len(summary["total_genes_affected"])
    
    return summary


def format_cnv_event_short_label(event: CNVEvent) -> str:
    """Compact label for a single threshold-triggered CNV event."""
    chrom = event.chromosome.replace("chr", "")
    if event.event_type.startswith("WHOLE_CHR_"):
        direction = event.event_type.replace("WHOLE_CHR_", "")
        return f"chr{chrom} {direction}"
    arm = event.arm or ""
    return f"chr{chrom}{arm} {event.event_type}"


def format_cnv_events_card_lines(events: List[CNVEvent]) -> Tuple[str, str]:
    """Whole-chromosome and arm-level summary lines for the CNV insight card."""
    whole_chr = [e for e in events if e.event_type.startswith("WHOLE_CHR_")]
    arm_level = [e for e in events if not e.event_type.startswith("WHOLE_CHR_")]

    if whole_chr:
        whole_text = "Whole chromosome: " + ", ".join(
            format_cnv_event_short_label(e) for e in whole_chr
        )
    else:
        whole_text = "Whole chromosome: none detected"

    if arm_level:
        arm_text = "Arm-level: " + ", ".join(
            format_cnv_event_short_label(e) for e in arm_level
        )
    else:
        arm_text = "Arm-level: none detected"

    return whole_text, arm_text


def format_cnv_events_section_summary(events: List[CNVEvent]) -> str:
    """One-line summary above the arm / whole-chromosome events table."""
    summary = get_cnv_summary(events)
    if summary["total_events"] <= 0:
        return "No threshold triggered CNV events detected"

    summary_text = f"Detected {summary['total_events']} CNV events: "
    parts: List[str] = []
    if summary["whole_chromosome_events"]:
        parts.append(f"{len(summary['whole_chromosome_events'])} whole chromosome")
    if summary["arm_events"]:
        parts.append(f"{len(summary['arm_events'])} arm-specific")
    if summary["gene_containing_events"]:
        parts.append(f"{summary['total_genes_affected']} genes affected")
    return summary_text + ", ".join(parts)


def _normalize_sex_label(xy_val: Any) -> str:
    try:
        s = str(xy_val).strip().upper()
        if s in ("MALE", "XY"):
            return "Male"
        if s in ("FEMALE", "XX"):
            return "Female"
    except Exception:
        pass
    return "Unknown"


def load_cytobands_df() -> pd.DataFrame:
    """Load UCSC cytoband definitions bundled with ROBIN."""
    try:
        res_path = importlib_resources.files("robin.resources") / "cytoBand.txt"
        return pd.read_csv(
            res_path,
            sep="\t",
            header=None,
            names=["chrom", "start_pos", "end_pos", "name", "stain"],
        )
    except Exception:
        return pd.DataFrame(columns=["chrom", "start_pos", "end_pos", "name", "stain"])


def load_gene_bed_for_sample(sample_dir: Path) -> pd.DataFrame:
    """Load panel gene BED for a sample (same resolution order as the CNV GUI)."""
    empty = pd.DataFrame(columns=["chrom", "start_pos", "end_pos", "gene"])
    try:
        panel = ""
        master_csv_path = sample_dir / "master.csv"
        if master_csv_path.exists():
            df = pd.read_csv(master_csv_path)
            if not df.empty and "analysis_panel" in df.columns:
                panel_val = df.iloc[0]["analysis_panel"]
                if panel_val and str(panel_val).strip():
                    panel = str(panel_val).strip()

        if not panel:
            bed_filename = "unique_genes.bed"
        elif panel == "rCNS2":
            bed_filename = "rCNS2_panel_name_uniq.bed"
        elif panel == "AML":
            bed_filename = "AML_panel_name_uniq.bed"
        else:
            bed_filename = f"{panel}_panel_name_uniq.bed"

        for name in (bed_filename, "unique_genes.bed"):
            try:
                res_path = importlib_resources.files("robin.resources") / name
                if res_path.exists():
                    return pd.read_csv(
                        res_path,
                        sep="\t",
                        header=None,
                        names=["chrom", "start_pos", "end_pos", "gene"],
                    )
            except Exception:
                continue
    except Exception:
        pass
    return empty


def detect_cnv_events_for_sample(sample_dir: Path) -> List[CNVEvent]:
    """Detect arm / whole-chromosome CNV events from on-disk sample outputs."""
    from robin.analysis.cnv_analysis import prepare_cnv_calling_track

    cnv_npy = sample_dir / "CNV.npy"
    cnv_dict_npy = sample_dir / "CNV_dict.npy"
    if not cnv_npy.exists() or not cnv_dict_npy.exists():
        return []

    try:
        cnv_map = np.load(cnv_npy, allow_pickle=True).item()
        cnv_dict = np.load(cnv_dict_npy, allow_pickle=True).item()
        if not isinstance(cnv_map, dict) or not isinstance(cnv_dict, dict):
            return []

        sex_estimate = "Unknown"
        xy_pkl = sample_dir / "XYestimate.pkl"
        if xy_pkl.exists():
            with xy_pkl.open("rb") as handle:
                loaded = pickle.load(handle)
            if loaded:
                sex_estimate = _normalize_sex_label(loaded)

        analysis_binw = int(cnv_dict.get("bin_width", 1_000_000))
        calling_cnv, calling_binw = prepare_cnv_calling_track(
            cnv_map, analysis_binw, sex_estimate
        )
        if not calling_cnv:
            return []

        return detect_cnv_events(
            cnv_data=calling_cnv,
            bin_width=int(calling_binw),
            sex_estimate=sex_estimate,
            cytobands_df=load_cytobands_df(),
            gene_df=load_gene_bed_for_sample(sample_dir),
        )
    except Exception as exc:
        logger.debug("CNV event detection failed for %s: %s", sample_dir, exc)
        return []
