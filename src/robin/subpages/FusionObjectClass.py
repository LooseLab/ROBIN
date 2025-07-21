"""
Gene Fusion Analysis Module

This module provides comprehensive gene fusion analysis capabilities for BAM files.
It includes structural variant detection, fusion candidate identification, and
visualization components.

Key Components:
- FusionProcessor: Core BAM processing and fusion detection
- StructuralVariantAnalyzer: SV detection and analysis
- FusionDataManager: Data storage and preprocessing
- FusionObject: Main analysis class coordinating all components

Performance Optimizations:
- Streaming BAM processing to reduce memory usage
- Efficient data structures and algorithms
- Pre-computed lookups for gene regions
- Vectorized operations where possible
"""

import os
import gff3_parser
import random
import pysam
import logging
import networkx as nx
import pandas as pd
import numpy as np
import re
import gc

try:
    import psutil

    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False
    psutil = None
from typing import Optional, Tuple, Dict, List, Set
from nicegui import run, background_tasks
from robin import resources
from pathlib import Path
from robin.subpages.base_analysis import BaseAnalysis
from robin.utilities.bed_file import MasterBedTree
from collections import defaultdict
from robin.core.state import state, ProcessState
from datetime import datetime
from dataclasses import dataclass
from enum import Enum
from robin.subpages.Fusion_object import build_breakpoint_graph

# Use the main logger configured in the main application
logger = logging.getLogger(__name__)


def get_memory_usage() -> float:
    """
    Get current memory usage in MB.

    Returns:
        float: Memory usage in MB, or -1 if psutil is not available
    """
    if not PSUTIL_AVAILABLE:
        return -1.0

    try:
        process = psutil.Process()
        return process.memory_info().rss / 1024 / 1024
    except Exception:
        return -1.0


def log_memory_usage(stage: str) -> None:
    """
    Log memory usage at a specific stage.

    Args:
        stage: Description of the current processing stage
    """
    memory_mb = get_memory_usage()
    if memory_mb >= 0:
        logger.debug(f"Memory usage at {stage}: {memory_mb:.1f} MB")
    else:
        logger.debug(f"Memory monitoring not available at {stage}")


def force_garbage_collection() -> None:
    """
    Force garbage collection and log memory usage.
    """
    gc.collect()
    memory_mb = get_memory_usage()
    if memory_mb >= 0:
        logger.debug(f"Memory after garbage collection: {memory_mb:.1f} MB")
    else:
        logger.debug("Memory after garbage collection: monitoring not available")


class EventType(Enum):
    """Enumeration of structural variant event types."""

    DELETION = "deletion"
    INVERSION = "inversion"
    TRANSLOCATION = "translocation"
    INTRA_TRANSLOCATION = "intra_translocation"


@dataclass
class GeneRegion:
    """Represents a gene region with start, end, and name."""

    start: int
    end: int
    name: str

    def overlaps_with(
        self, other_start: int, other_end: int, min_overlap: int = 100
    ) -> bool:
        """Check if this region overlaps with another region."""
        overlap_start = max(self.start, other_start)
        overlap_end = min(self.end, other_end)
        return (
            overlap_end > overlap_start and (overlap_end - overlap_start) > min_overlap
        )


@dataclass
class FusionCandidate:
    """Represents a fusion candidate with all relevant information."""

    read_id: str
    gene_name: str
    chromosome: str
    start: int
    end: int
    strand: str
    mapping_quality: int
    read_start: int
    read_end: int
    mapping_span: int


class FusionDataManager:
    """
    Manages fusion data storage and preprocessing operations.

    Responsibilities:
    - Efficient data storage with memory optimization
    - Data preprocessing for visualization
    - File I/O operations
    - Data validation and cleaning
    """

    def __init__(self):
        self.fusion_candidates: Dict[str, pd.DataFrame] = {}
        self.fusion_candidates_all: Dict[str, pd.DataFrame] = {}
        self.structural_variants: Dict[str, pd.DataFrame] = {}

    def add_fusion_candidates(
        self, sample_id: str, candidates: pd.DataFrame, is_target_panel: bool = True
    ) -> None:
        """
        Add fusion candidates for a sample with efficient memory management.

        Args:
            sample_id: Sample identifier
            candidates: Fusion candidate DataFrame
            is_target_panel: Whether candidates are from target panel or genome-wide
        """
        if candidates is None or candidates.empty:
            return

        # Optimize DataFrame memory usage
        candidates = self._optimize_dataframe_memory(candidates)

        target_dict = (
            self.fusion_candidates if is_target_panel else self.fusion_candidates_all
        )

        if sample_id not in target_dict:
            target_dict[sample_id] = candidates
        else:
            # Efficient concatenation without unnecessary copies
            target_dict[sample_id] = pd.concat(
                [target_dict[sample_id], candidates], ignore_index=True, copy=False
            )

        logger.info(f"Added {len(candidates)} fusion candidates for sample {sample_id}")

    def _optimize_dataframe_memory(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Optimize DataFrame memory usage through dtype optimization.

        Args:
            df: Input DataFrame

        Returns:
            Memory-optimized DataFrame
        """
        # Define optimal dtypes for common columns
        dtype_optimizations = {
            "read_id": "category",
            "col4": "category",  # Gene name
            "reference_id": "category",  # Chromosome
            "strand": "category",
            "mapping_quality": np.int8,
            "reference_start": np.int32,
            "reference_end": np.int32,
            "read_start": np.int32,
            "read_end": np.int32,
            "mapping_span": np.int32,
        }

        # Apply optimizations only for columns that exist
        existing_columns = {
            col: dtype
            for col, dtype in dtype_optimizations.items()
            if col in df.columns
        }

        return df.astype(existing_columns, errors="ignore")

    def get_filtered_candidates(
        self, sample_id: str, is_target_panel: bool = True
    ) -> Optional[pd.DataFrame]:
        """
        Get filtered fusion candidates for a sample.

        Args:
            sample_id: Sample identifier
            is_target_panel: Whether to get target panel or genome-wide candidates

        Returns:
            Filtered DataFrame or None if no candidates
        """
        target_dict = (
            self.fusion_candidates if is_target_panel else self.fusion_candidates_all
        )

        if sample_id not in target_dict or target_dict[sample_id].empty:
            return None

        df = target_dict[sample_id]

        # Filter for duplicated reads efficiently
        duplicated_mask = df["read_id"].duplicated(keep=False)
        if not duplicated_mask.any():
            return None

        doubles = df[duplicated_mask]
        gene_counts = doubles.groupby("read_id", observed=True)["col4"].transform(
            "nunique"
        )
        result = doubles[gene_counts > 1]

        return result if not result.empty else None

    def save_candidates_to_csv(self, sample_id: str, output_dir: str) -> None:
        """
        Save fusion candidates to CSV files with preprocessing.

        Args:
            sample_id: Sample identifier
            output_dir: Output directory path
        """
        # Save target panel candidates
        target_candidates = self.get_filtered_candidates(
            sample_id, is_target_panel=True
        )
        if target_candidates is not None:
            target_path = os.path.join(output_dir, "fusion_candidates_master.csv")
            target_candidates.to_csv(target_path, index=False)

            # Pre-process for visualization
            self._preprocess_for_visualization(
                target_candidates,
                os.path.join(output_dir, "fusion_candidates_master_processed.pkl"),
            )

        # Save genome-wide candidates
        all_candidates = self.get_filtered_candidates(sample_id, is_target_panel=False)
        if all_candidates is not None:
            all_path = os.path.join(output_dir, "fusion_candidates_all.csv")
            all_candidates.to_csv(all_path, index=False)

            # Pre-process for visualization
            self._preprocess_for_visualization(
                all_candidates,
                os.path.join(output_dir, "fusion_candidates_all_processed.pkl"),
            )

    def _preprocess_for_visualization(
        self, candidates: pd.DataFrame, output_path: str
    ) -> None:
        """
        Pre-process fusion candidates for efficient visualization loading.

        Args:
            candidates: Fusion candidate DataFrame
            output_path: Path to save preprocessed data
        """
        try:
            # Annotate results
            annotated_data, goodpairs = self._annotate_candidates(candidates)

            # Create processed data structure
            processed_data = {
                "annotated_data": annotated_data,
                "goodpairs": goodpairs,
                "gene_pairs": [],
                "gene_groups": [],
                "candidate_count": 0,
            }

            # Process gene pairs and groups if we have good pairs
            if not annotated_data.empty and goodpairs.any():
                gene_pairs = self._extract_gene_pairs(annotated_data, goodpairs)
                gene_groups = self._build_gene_networks(
                    gene_pairs, annotated_data, goodpairs
                )

                processed_data.update(
                    {
                        "gene_pairs": gene_pairs,
                        "gene_groups": gene_groups,
                        "candidate_count": len(gene_groups),
                    }
                )

            # Save as pickle for efficient loading
            import pickle

            with open(output_path, "wb") as f:
                pickle.dump(processed_data, f)

            logger.info(f"Pre-processed fusion data saved to {output_path}")

        except Exception as e:
            logger.error(f"Error pre-processing fusion data: {str(e)}")
            logger.error("Exception details:", exc_info=True)

    def _annotate_candidates(
        self, candidates: pd.DataFrame
    ) -> Tuple[pd.DataFrame, pd.Series]:
        """
        Annotate fusion candidates with tags and colors.

        Args:
            candidates: Raw fusion candidate DataFrame

        Returns:
            Tuple of (annotated_data, goodpairs)
        """
        # Create gene tags efficiently
        gene_lookup = candidates.groupby("read_id", observed=True)["col4"].agg(
            lambda x: ",".join(set(x))
        )
        candidates = candidates.copy()
        candidates["tag"] = candidates["read_id"].map(gene_lookup)

        # Generate colors efficiently
        color_lookup = candidates.groupby("read_id", observed=True)["col4"].apply(
            lambda x: self._generate_random_color()
        )
        candidates["Color"] = candidates["read_id"].map(color_lookup)

        # Clean string data
        candidates = candidates.apply(lambda x: x.strip() if isinstance(x, str) else x)

        # Find good pairs (reads mapping to multiple genes)
        goodpairs = (
            candidates.groupby("tag", observed=True)["read_id"].transform("nunique") > 2
        )

        return candidates, goodpairs

    def _extract_gene_pairs(
        self, annotated_data: pd.DataFrame, goodpairs: pd.Series
    ) -> List[Tuple[str, str]]:
        """Extract gene pairs from annotated data."""
        gene_pairs = (
            annotated_data[goodpairs]
            .sort_values(by="reference_start")["tag"]
            .unique()
            .tolist()
        )
        return [tuple(pair.split(",")) for pair in gene_pairs]

    def _build_gene_networks(
        self,
        gene_pairs: List[Tuple[str, str]],
        annotated_data: pd.DataFrame,
        goodpairs: pd.Series,
    ) -> List[List[str]]:
        """Build gene networks from gene pairs."""
        if not gene_pairs:
            return []

        # Create network graph
        G = nx.Graph()
        for pair in gene_pairs:
            G.add_edge(pair[0], pair[1])

        # Get connected components
        connected_components = list(nx.connected_components(G))
        gene_groups = []

        for component in connected_components:
            component_genes = list(component)
            reads = self._get_reads_for_genes(
                annotated_data, goodpairs, component_genes
            )
            if len(reads) > 1:
                gene_groups.append(component_genes)

        return gene_groups

    def _get_reads_for_genes(
        self, annotated_data: pd.DataFrame, goodpairs: pd.Series, genes: List[str]
    ) -> pd.DataFrame:
        """Get reads for specific genes."""
        reads = annotated_data[goodpairs][annotated_data[goodpairs]["col4"].isin(genes)]

        # Rename columns for consistency
        reads.columns = [
            "chromosome",
            "start",
            "end",
            "gene",
            "chromosome2",
            "start2",
            "end2",
            "id",
            "quality",
            "strand",
            "read_start",
            "read_end",
            "secondary",
            "supplementary",
            "span",
            "tag",
            "color",
        ]

        # Convert coordinates to int
        reads["start"] = pd.to_numeric(reads["start"], errors="coerce").astype("Int64")
        reads["end"] = pd.to_numeric(reads["end"], errors="coerce").astype("Int64")

        # Sort and remove duplicates
        reads = reads.sort_values(by=["chromosome", "start", "end"])
        reads = reads.drop_duplicates(subset=["start2", "end2", "id"])

        # Collapse ranges within each chromosome
        result = (
            reads.groupby("chromosome", observed=True)
            .apply(lambda x: self._collapse_ranges(x, 10000))
            .reset_index(drop=True)
        )

        return result

    def _collapse_ranges(self, df: pd.DataFrame, max_distance: int) -> pd.DataFrame:
        """Collapse ranges within a fixed distance."""
        if df.empty:
            return df

        collapsed = []
        current_range = None

        for _, row in df.iterrows():
            start, end = row["start"], row["end"]

            if current_range is None:
                current_range = row
            else:
                if start <= current_range["end"] + max_distance:
                    current_range["end"] = max(current_range["end"], end)
                else:
                    collapsed.append(current_range)
                    current_range = row

        if current_range is not None:
            collapsed.append(current_range)

        return pd.DataFrame(collapsed)

    def _generate_random_color(self) -> str:
        """Generate a random hex color."""
        return "#{:06x}".format(random.randint(0, 0xFFFFFF))


class StructuralVariantAnalyzer:
    """
    Analyzes structural variants from BAM files.

    Responsibilities:
    - Extract split read alignments
    - Build structural variant links
    - Classify structural variant events
    - Generate breakpoint summaries
    """

    def __init__(self):
        self.min_mapping_quality = 55
        self.min_overlap = 100

    def process_bam_for_svs(self, bamfile: str) -> pd.DataFrame:
        """
        Process BAM file to extract structural variant information with memory optimization.

        Args:
            bamfile: Path to BAM file

        Returns:
            DataFrame with structural variant links

        Memory optimizations:
        - Streaming processing
        - Efficient data structures
        - Immediate cleanup of intermediate data
        """
        # Extract split read alignments
        split_reads_df = self._extract_split_read_alignments(bamfile)

        if split_reads_df.empty:
            logger.debug("No split reads found in BAM file")
            return pd.DataFrame()

        # Annotate the extracted data
        annotated_df = self._annotate_split_reads(split_reads_df)

        # Free the intermediate DataFrame immediately
        del split_reads_df

        if annotated_df.empty:
            logger.debug("No annotated data after filtering")
            return pd.DataFrame()

        # Build links from annotated data
        links_df = self._build_structural_variant_links(annotated_df)

        # Free the intermediate DataFrame immediately
        del annotated_df

        return links_df

    def _extract_split_read_alignments(self, bam_path: str) -> pd.DataFrame:
        """
        Extract split read alignments from BAM file with optimized performance.

        Args:
            bam_path: Path to BAM file

        Returns:
            DataFrame with split read alignments
        """
        # Pre-allocate lists for better memory efficiency
        data_columns = {
            "qnames": [],
            "types": [],
            "rnames": [],
            "strands": [],
            "ref_starts": [],
            "ref_ends": [],
            "ref_spans": [],
            "read_starts": [],
            "read_ends": [],
            "read_spans": [],
            "left_softs": [],
            "left_hards": [],
            "right_softs": [],
            "right_hards": [],
            "mqs": [],
        }

        # Single pass through BAM file with early filtering
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for aln in bam.fetch(until_eof=True):
                # Early filtering: skip unmapped and secondary alignments
                if aln.is_unmapped or aln.is_secondary:
                    continue

                # Early filtering: only process reads with supplementary alignments or SA tag
                if not (aln.is_supplementary or aln.has_tag("SA")):
                    continue

                # Parse CIGAR string efficiently
                cigar_info = self._parse_cigar_string(aln.cigartuples or [])

                # Get reference coordinates
                ref_start = aln.reference_start
                ref_end = aln.reference_end
                ref_span = ref_end - ref_start

                # Get read coordinates
                read_start = aln.query_alignment_start + cigar_info["left_hard"]
                read_end = aln.query_alignment_end + cigar_info["left_hard"]
                read_span = aln.query_alignment_end - aln.query_alignment_start

                # Append data efficiently
                data_columns["qnames"].append(aln.query_name)
                data_columns["types"].append(
                    "SUPPLEMENTARY" if aln.is_supplementary else "PRIMARY"
                )
                data_columns["rnames"].append(bam.get_reference_name(aln.reference_id))
                data_columns["strands"].append("-" if aln.is_reverse else "+")
                data_columns["ref_starts"].append(ref_start)
                data_columns["ref_ends"].append(ref_end)
                data_columns["ref_spans"].append(ref_span)
                data_columns["read_starts"].append(read_start)
                data_columns["read_ends"].append(read_end)
                data_columns["read_spans"].append(read_span)
                data_columns["left_softs"].append(cigar_info["left_soft"])
                data_columns["left_hards"].append(cigar_info["left_hard"])
                data_columns["right_softs"].append(cigar_info["right_soft"])
                data_columns["right_hards"].append(cigar_info["right_hard"])
                data_columns["mqs"].append(aln.mapping_quality)

        # Create DataFrame efficiently with pre-allocated lists
        return pd.DataFrame(
            {
                "QNAME": data_columns["qnames"],
                "TYPE": data_columns["types"],
                "RNAME": data_columns["rnames"],
                "STRAND": data_columns["strands"],
                "REF_START": data_columns["ref_starts"],
                "REF_END": data_columns["ref_ends"],
                "REF_SPAN": data_columns["ref_spans"],
                "READ_START": data_columns["read_starts"],
                "READ_END": data_columns["read_ends"],
                "READ_SPAN": data_columns["read_spans"],
                "LEFT_SOFT": data_columns["left_softs"],
                "LEFT_HARD": data_columns["left_hards"],
                "RIGHT_SOFT": data_columns["right_softs"],
                "RIGHT_HARD": data_columns["right_hards"],
                "MQ": data_columns["mqs"],
            }
        )

    def _parse_cigar_string(self, cigartuples: List[Tuple[int, int]]) -> Dict[str, int]:
        """
        Parse CIGAR string to extract clipping information.

        Args:
            cigartuples: CIGAR operations as tuples

        Returns:
            Dictionary with clipping information
        """
        # Compute left-end clipping efficiently
        left_soft = left_hard = 0
        for op, length in cigartuples:
            if op == 4:  # soft clip
                left_soft += length
            elif op == 5:  # hard clip
                left_hard += length
            else:
                break

        # Compute right-end clipping efficiently
        right_soft = right_hard = 0
        for op, length in reversed(cigartuples):
            if op == 4:
                right_soft += length
            elif op == 5:
                right_hard += length
            else:
                break

        return {
            "left_soft": left_soft,
            "left_hard": left_hard,
            "right_soft": right_soft,
            "right_hard": right_hard,
        }

    def _annotate_split_reads(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Annotate split reads with filtering and optimization.

        Args:
            df: Split reads DataFrame

        Returns:
            Annotated DataFrame
        """
        if df.empty:
            return df

        # Apply filters first to reduce memory usage (vectorized operations)
        mq_mask = df["MQ"].values >= self.min_mapping_quality
        chr_mask = df["RNAME"].values != "chrM"
        combined_mask = mq_mask & chr_mask

        # Apply filter inplace to avoid copy
        df = df[combined_mask].reset_index(drop=True)

        if df.empty:
            return df

        # Optimize dtypes for memory efficiency
        dtype_optimizations = {
            "QNAME": "category",
            "TYPE": "category",
            "RNAME": "category",
            "REF_START": np.int32,
            "REF_END": np.int32,
            "REF_SPAN": np.int32,
            "READ_START": np.int32,
            "READ_END": np.int32,
            "READ_SPAN": np.int32,
            "MQ": np.int8,
            "STRAND": "category",
        }

        df = df.astype(dtype_optimizations, errors="ignore")

        # Sort efficiently
        df = df.sort_values(["QNAME", "TYPE", "REF_START"], ignore_index=True)

        # Filter for primary alignments efficiently
        primary_mask = df["TYPE"] == "PRIMARY"
        primary_qnames = df.loc[primary_mask, "QNAME"].unique()

        # Filter for reads that exist in primary set
        df = df[df["QNAME"].isin(primary_qnames)]

        # Filter for duplicated reads
        duplicated_mask = df["QNAME"].duplicated(keep=False)
        df = df[duplicated_mask].reset_index(drop=True)

        if df.empty:
            return df

        # Add piece_order efficiently
        df["piece_order"] = df.groupby("QNAME", observed=True)["READ_START"].transform(
            lambda x: x.rank(method="first").astype(int)
        )

        return df

    def _build_structural_variant_links(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Build structural variant links from annotated split reads.

        Args:
            df: Annotated split reads DataFrame

        Returns:
            DataFrame with structural variant links
        """
        if df.empty:
            return pd.DataFrame()

        # Pre-allocate lists for better memory efficiency
        link_data = {
            "qnames": [],
            "rname1s": [],
            "rname2s": [],
            "coord1s": [],
            "coord2s": [],
            "genomic_gaps": [],
            "events": [],
        }

        # Cache natural key function
        def natural_key(s: str):
            return [
                int(tok) if tok.isdigit() else tok.lower()
                for tok in re.split(r"(\d+)", s)
            ]

        # Process groups efficiently
        for qname, grp in df.groupby("QNAME", sort=False, observed=True):
            # Sort by piece_order once
            seq = grp.sort_values("piece_order")

            if len(seq) < 2:
                continue

            # Process consecutive pairs efficiently
            for i in range(len(seq) - 1):
                p1 = seq.iloc[i]
                p2 = seq.iloc[i + 1]

                # Compute tail/head
                tail1 = p1.REF_END
                head2 = p2.REF_START

                # Classify event type
                event_type = self._classify_event_type(p1, p2)

                # Determine ordering
                five1 = p1.REF_START if p1.STRAND == "+" else p1.REF_END
                five2 = p2.REF_START if p2.STRAND == "+" else p2.REF_END

                nk1 = natural_key(p1.RNAME)
                nk2 = natural_key(p2.RNAME)

                first_is_p1 = (nk1 < nk2) or (nk1 == nk2 and five1 <= five2)

                if first_is_p1:
                    r1, r2 = p1.RNAME, p2.RNAME
                    c1, c2 = tail1, head2
                else:
                    r1, r2 = p2.RNAME, p1.RNAME
                    c1, c2 = head2, tail1

                # Calculate genomic gap
                genomic_gap = -1 if r1 != r2 else (c2 - c1)

                # Append data efficiently
                link_data["qnames"].append(qname)
                link_data["rname1s"].append(r1)
                link_data["rname2s"].append(r2)
                link_data["coord1s"].append(c1)
                link_data["coord2s"].append(c2)
                link_data["genomic_gaps"].append(genomic_gap)
                link_data["events"].append(event_type.value)

        # Create DataFrame efficiently
        return pd.DataFrame(
            {
                "QNAME": link_data["qnames"],
                "RNAME.1": link_data["rname1s"],
                "RNAME.2": link_data["rname2s"],
                "coord_1": link_data["coord1s"],
                "coord_2": link_data["coord2s"],
                "genomic_gap": link_data["genomic_gaps"],
                "event": link_data["events"],
            }
        )

    def _classify_event_type(self, p1: pd.Series, p2: pd.Series) -> EventType:
        """
        Classify the type of structural variant event.

        Args:
            p1: First piece of the split read
            p2: Second piece of the split read

        Returns:
            EventType classification
        """
        if p1.RNAME == p2.RNAME:
            # Same chromosome events
            if p1.STRAND == p2.STRAND:
                if p1.REF_START < p2.REF_START:
                    return EventType.DELETION
                else:
                    return EventType.INTRA_TRANSLOCATION
            else:
                return EventType.INVERSION
        else:
            # Different chromosomes
            return EventType.TRANSLOCATION

    def get_structural_variant_summary(
        self, links_df: pd.DataFrame, min_support: int = 2
    ) -> pd.DataFrame:
        """
        Generate summary of structural variants with support filtering.

        Args:
            links_df: Structural variant links DataFrame
            min_support: Minimum number of supporting reads

        Returns:
            Summary DataFrame
        """
        if links_df.empty:
            return pd.DataFrame()

        summary = (
            links_df.groupby(
                ["RNAME.1", "coord_1", "RNAME.2", "coord_2"], observed=True
            )
            .agg(
                support_count=("QNAME", "nunique"),
                supporting_reads=("QNAME", lambda s: list(s.unique())),
                from_med=("coord_1", "median"),
                to_med=("coord_2", "median"),
                median_genomic_gap=("genomic_gap", "median"),
                event_counts=("event", lambda s: s.value_counts().to_dict()),
                predominant_event=("event", lambda s: s.mode().iloc[0]),
            )
            .reset_index()
        )

        # Filter by minimum support
        return summary[summary["support_count"] >= min_support].copy()


class FusionProcessor:
    """
    Core processor for fusion detection from BAM files.

    Responsibilities:
    - Load and cache gene region data
    - Process BAM files for fusion candidates
    - Coordinate with other analysis components
    """

    def __init__(self, target_panel: str):
        self.target_panel = target_panel
        self.gene_regions_cache: Dict[str, List[GeneRegion]] = {}
        self.all_gene_regions_cache: Dict[str, List[GeneRegion]] = {}

        # Initialize file paths
        self._setup_file_paths()

        # Load gene regions into memory
        self._load_gene_regions()

    def _setup_file_paths(self) -> None:
        """Setup file paths for gene data."""
        resources_dir = os.path.dirname(os.path.abspath(resources.__file__))

        # Set up gene BED files based on target panel
        if self.target_panel == "rCNS2":
            self.gene_bed = os.path.join(resources_dir, "rCNS2_panel_name_uniq.bed")
        elif self.target_panel == "AML":
            self.gene_bed = os.path.join(resources_dir, "AML_panel_name_uniq.bed")
        else:
            raise ValueError(f"Unsupported target panel: {self.target_panel}")

        self.all_gene_bed = os.path.join(resources_dir, "all_genes2.bed")
        self.gene_gff3 = os.path.join(
            resources_dir, "gencode.v45.basic.annotation.gff3"
        )

    def _load_gene_regions(self) -> None:
        """Load gene regions into memory for efficient lookup."""
        self.gene_regions_cache = self._load_bed_regions(self.gene_bed)
        self.all_gene_regions_cache = self._load_bed_regions(self.all_gene_bed)

        logger.info(
            f"Loaded {sum(len(regions) for regions in self.gene_regions_cache.values())} "
            f"target gene regions and {sum(len(regions) for regions in self.all_gene_regions_cache.values())} "
            f"genome-wide gene regions"
        )

    def _load_bed_regions(self, bed_file: str) -> Dict[str, List[GeneRegion]]:
        """
        Load BED file regions into memory for efficient lookup.

        Args:
            bed_file: Path to BED file

        Returns:
            Dictionary mapping chromosome to list of GeneRegion objects
        """
        regions = defaultdict(list)

        try:
            with open(bed_file, "r") as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:
                        continue

                    parts = line.split("\t")
                    if len(parts) < 4:
                        logger.warning(
                            f"Skipping malformed line {line_num} in {bed_file}: {line}"
                        )
                        continue

                    try:
                        chrom = parts[0]
                        start = int(parts[1])
                        end = int(parts[2])
                        gene_name = parts[3]

                        regions[chrom].append(GeneRegion(start, end, gene_name))
                    except ValueError as e:
                        logger.warning(
                            f"Skipping line {line_num} in {bed_file} due to invalid coordinates: {e}"
                        )
                        continue

        except FileNotFoundError:
            logger.error(f"BED file not found: {bed_file}")
            raise
        except Exception as e:
            logger.error(f"Error loading BED file {bed_file}: {e}")
            raise

        return dict(regions)

    def has_supplementary_alignments(self, bamfile: str) -> bool:
        """
        Quickly check if a BAM file has supplementary alignments.

        Args:
            bamfile: Path to the BAM file

        Returns:
            True if supplementary alignment is found, else False
        """
        try:
            with pysam.AlignmentFile(bamfile, "rb", check_sq=False) as bam:
                for read in bam.fetch(until_eof=True):
                    if read.is_supplementary:
                        return True
            return False
        except Exception as e:
            logger.error(f"Error checking supplementary alignments in {bamfile}: {e}")
            return False

    def find_reads_with_supplementary(self, bamfile: str) -> Set[str]:
        """
        Find all read names that have supplementary alignments.

        Args:
            bamfile: Path to BAM file

        Returns:
            Set of read names with supplementary alignments
        """
        reads_with_supp = set()

        try:
            with pysam.AlignmentFile(bamfile, "rb") as bam:
                for read in bam:
                    if read.is_unmapped:
                        continue
                    # Check for supplementary alignments or SA tag
                    if read.is_supplementary or read.has_tag("SA"):
                        reads_with_supp.add(read.query_name)

        except Exception as e:
            logger.error(
                f"Error finding reads with supplementary alignments in {bamfile}: {e}"
            )
            raise

        logger.debug(
            f"Found {len(reads_with_supp)} reads with supplementary alignments"
        )
        return reads_with_supp

    def process_bam_for_fusions(
        self, bamfile: str
    ) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
        """
        Process BAM file to find fusion candidates with memory optimization.

        Args:
            bamfile: Path to BAM file

        Returns:
            Tuple of (target_panel_candidates, genome_wide_candidates) DataFrames

        Memory optimizations:
        - Streaming processing
        - Efficient data structures
        - Immediate cleanup of intermediate data
        """
        try:
            # Find reads with supplementary alignments
            reads_with_supp = self.find_reads_with_supplementary(bamfile)

            if not reads_with_supp:
                logger.debug("No reads with supplementary alignments found")
                return None, None

            # Process reads for target panel fusions
            target_candidates = self._process_reads_for_fusions(
                bamfile, reads_with_supp, self.gene_regions_cache
            )

            # Process reads for genome-wide fusions
            genome_wide_candidates = self._process_reads_for_fusions(
                bamfile, reads_with_supp, self.all_gene_regions_cache
            )

            # Clear the reads_with_supp set to free memory
            reads_with_supp.clear()
            del reads_with_supp

            return target_candidates, genome_wide_candidates

        except Exception as e:
            logger.error(f"Error processing BAM file for fusions: {str(e)}")
            logger.error("Exception details:", exc_info=True)
            raise

    def _process_reads_for_fusions(
        self,
        bamfile: str,
        reads_with_supp: Set[str],
        gene_regions: Dict[str, List[GeneRegion]],
    ) -> Optional[pd.DataFrame]:
        """
        Process reads to find gene intersections and create fusion candidates.

        Args:
            bamfile: Path to BAM file
            reads_with_supp: Set of read names with supplementary alignments
            gene_regions: Dictionary of gene regions by chromosome

        Returns:
            DataFrame with fusion candidates or None if no candidates found
        """
        rows = []

        try:
            with pysam.AlignmentFile(bamfile, "rb") as bam:
                for read in bam:
                    # Skip if not in our target reads
                    if read.query_name not in reads_with_supp:
                        continue

                    # Skip secondary alignments and unmapped reads
                    if read.is_secondary or read.is_unmapped:
                        continue

                    # Get reference information
                    ref_name = (
                        bam.get_reference_name(read.reference_id)
                        if read.reference_id >= 0
                        else None
                    )
                    if not ref_name or ref_name == "chrM":
                        continue

                    ref_start = read.reference_start
                    ref_end = read.reference_end

                    # Check if this read intersects with any gene regions
                    if ref_name in gene_regions:
                        read_rows = self._find_gene_intersections(
                            read, ref_name, ref_start, ref_end, gene_regions[ref_name]
                        )
                        rows.extend(read_rows)

        except Exception as e:
            logger.error(f"Error processing reads for fusions: {str(e)}")
            raise

        if not rows:
            return None

        # Create DataFrame with optimized dtypes
        df = pd.DataFrame(rows)

        # Apply memory optimizations
        df = self._optimize_fusion_dataframe(df)

        # Apply quality filters
        df = df[(df["mapping_quality"] > 40) & (df["mapping_span"] > 100)].reset_index(
            drop=True
        )

        return df if not df.empty else None

    def _find_gene_intersections(
        self,
        read: pysam.AlignedSegment,
        ref_name: str,
        ref_start: int,
        ref_end: int,
        gene_regions: List[GeneRegion],
    ) -> List[Dict]:
        """
        Find intersections between a read and gene regions.

        Args:
            read: Pysam aligned segment
            ref_name: Reference chromosome name
            ref_start: Reference start position
            ref_end: Reference end position
            gene_regions: List of gene regions on this chromosome

        Returns:
            List of intersection dictionaries
        """
        read_rows = []

        for gene_region in gene_regions:
            if gene_region.overlaps_with(ref_start, ref_end, min_overlap=100):
                read_rows.append(
                    {
                        "col1": ref_name,
                        "col2": gene_region.start,
                        "col3": gene_region.end,
                        "col4": gene_region.name,
                        "reference_id": ref_name,
                        "reference_start": ref_start,
                        "reference_end": ref_end,
                        "read_id": read.query_name,
                        "mapping_quality": read.mapping_quality,
                        "strand": "-" if read.is_reverse else "+",
                        "read_start": read.query_alignment_start,
                        "read_end": read.query_alignment_end,
                        "is_secondary": read.is_secondary,
                        "is_supplementary": read.is_supplementary,
                        "mapping_span": ref_end - ref_start,
                    }
                )

        return read_rows

    def _optimize_fusion_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Optimize fusion DataFrame memory usage.

        Args:
            df: Input DataFrame

        Returns:
            Memory-optimized DataFrame
        """
        # Use categorical dtypes for string columns
        string_columns = ["col1", "col4", "reference_id", "strand", "read_id"]
        for col in string_columns:
            if col in df.columns:
                df[col] = df[col].astype("category")

        # Use appropriate integer types
        int_columns = [
            "col2",
            "col3",
            "reference_start",
            "reference_end",
            "read_start",
            "read_end",
            "mapping_quality",
            "mapping_span",
        ]
        for col in int_columns:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")

        return df


# Legacy function wrappers for backward compatibility
def fusion_work_pysam(
    bamfile: str, gene_bed: str, all_gene_bed: str
) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """
    Legacy wrapper for fusion detection using pure pysam with memory optimization.

    Args:
        bamfile: Path to the BAM file
        gene_bed: Path to target gene BED file
        all_gene_bed: Path to all genes BED file

    Returns:
        Tuple of (fusion_candidates, fusion_candidates_all) DataFrames

    Memory optimizations:
    - Streaming BAM processing
    - Efficient data structures
    - Immediate cleanup of intermediate data
    """
    # Extract target panel from gene_bed path
    if "rCNS2" in gene_bed:
        target_panel = "rCNS2"
    elif "AML" in gene_bed:
        target_panel = "AML"
    else:
        raise ValueError("Unable to determine target panel from gene_bed path")

    processor = FusionProcessor(target_panel)
    result = processor.process_bam_for_fusions(bamfile)

    # Explicit cleanup of processor to free memory
    del processor

    return result


def process_bam_pipeline(bamfile: str) -> pd.DataFrame:
    """
    Legacy wrapper for structural variant processing pipeline with memory optimization.

    Args:
        bamfile: Path to the BAM file to process

    Returns:
        pd.DataFrame: DataFrame containing structural variant links

    Memory optimizations:
    - Streaming processing
    - Efficient data structures
    - Immediate cleanup of intermediate data
    """
    analyzer = StructuralVariantAnalyzer()
    result = analyzer.process_bam_for_svs(bamfile)

    # Explicit cleanup of analyzer to free memory
    del analyzer

    return result


def has_supplementary(bam_file_path: str) -> bool:
    """
    Legacy wrapper for checking supplementary alignments with memory optimization.

    Args:
        bam_file_path: Path to the BAM file

    Returns:
        bool: True if supplementary alignment is found, else False

    Memory optimizations:
    - Early exit on first match
    - Minimal memory footprint
    """
    # Extract target panel from path (default to rCNS2 if unclear)
    target_panel = "rCNS2"  # Default fallback
    processor = FusionProcessor(target_panel)
    result = processor.has_supplementary_alignments(bam_file_path)

    # Explicit cleanup of processor to free memory
    del processor

    return result


def process_reads_for_fusions(
    bamfile: str,
    reads_with_supp: Set[str],
    gene_regions: Dict[str, List[Tuple[int, int, str]]],
) -> Optional[pd.DataFrame]:
    """
    Process reads to find gene intersections and create fusion candidates.

    Args:
        bamfile: Path to BAM file
        reads_with_supp: Set of read names with supplementary alignments
        gene_regions: Dictionary of gene regions by chromosome

    Returns:
        DataFrame with fusion candidates or None if no candidates found
    """
    rows = []

    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for read in bam:
            # Skip if not in our target reads
            if read.query_name not in reads_with_supp:
                continue

            # Skip secondary alignments
            if read.is_secondary:
                continue

            # Skip unmapped reads
            if read.is_unmapped:
                continue

            # Get reference information
            ref_name = (
                bam.get_reference_name(read.reference_id)
                if read.reference_id >= 0
                else None
            )
            if not ref_name:
                continue

            # Exclude chrM early to avoid unnecessary processing
            if ref_name == "chrM":
                continue

            ref_start = read.reference_start
            ref_end = read.reference_end

            # Check if this read intersects with any gene regions
            if ref_name in gene_regions:
                # Use list comprehension for better performance
                read_rows = [
                    {
                        "col1": ref_name,
                        "col2": gene_start,
                        "col3": gene_end,
                        "col4": gene_name,
                        "reference_id": ref_name,
                        "reference_start": ref_start,
                        "reference_end": ref_end,
                        "read_id": read.query_name,
                        "mapping_quality": read.mapping_quality,
                        "strand": "-" if read.is_reverse else "+",
                        "read_start": read.query_alignment_start,
                        "read_end": read.query_alignment_end,
                        "is_secondary": read.is_secondary,
                        "is_supplementary": read.is_supplementary,
                        "mapping_span": ref_end - ref_start,
                    }
                    for gene_start, gene_end, gene_name in gene_regions[ref_name]
                    if (
                        gene_start < ref_end
                        and gene_end > ref_start  # Overlap check
                        and min(ref_end, gene_end) - max(ref_start, gene_start) > 100
                    )  # Length check
                ]
                rows.extend(read_rows)

    if not rows:
        return None

    # Create DataFrame with optimized dtypes
    df = pd.DataFrame(rows)

    # Use categorical dtypes for string columns
    string_columns = ["col1", "col4", "reference_id", "strand", "read_id"]
    for col in string_columns:
        if col in df.columns:
            df[col] = df[col].astype("category")

    # Apply filters more efficiently
    df = df[(df["mapping_quality"] > 40) & (df["mapping_span"] > 100)].reset_index(
        drop=True
    )

    return df if not df.empty else None


def important_function(combined_df):
    """
    Generate BED format lines from structural variant data with memory optimization.

    Args:
        combined_df: Combined structural variant DataFrame

    Returns:
        List of BED format lines

    Memory optimizations:
    - Efficient iteration
    - Minimal string operations
    - Immediate cleanup of intermediate data
    """
    bed_lines = []
    if not combined_df.empty:
        # Get summary of structural variants
        sv_summary = get_summary(combined_df, min_support=2)

        if not sv_summary.empty:
            # Convert summary to BED format lines - only breakpoint boundaries
            for _, row in sv_summary.iterrows():
                # Create BED line for primary breakpoint boundary
                chrom1 = row["RNAME.1"]
                coord1 = row["coord_1"]

                # Define breakpoint boundary window (e.g., 500bp on each side)
                boundary_window = 500  # 500bp window around breakpoint
                start1 = max(0, coord1 - boundary_window)
                end1 = coord1 + boundary_window

                # Create BED lines for primary breakpoint on both strands
                bed_line1_plus = f"{chrom1}\t{start1}\t{end1}\t{row['predominant_event']}_breakpoint1\t{row['support_count']}\t+"
                bed_line1_minus = f"{chrom1}\t{start1}\t{end1}\t{row['predominant_event']}_breakpoint1\t{row['support_count']}\t-"
                bed_lines.append(bed_line1_plus)
                bed_lines.append(bed_line1_minus)

                # Handle different event types
                if row["RNAME.1"] != row["RNAME.2"]:
                    # Translocation: add partner breakpoint boundary on both strands
                    chrom2 = row["RNAME.2"]
                    coord2 = row["coord_2"]

                    start2 = max(0, coord2 - boundary_window)
                    end2 = coord2 + boundary_window

                    bed_line2_plus = f"{chrom2}\t{start2}\t{end2}\t{row['predominant_event']}_breakpoint2\t{row['support_count']}\t+"
                    bed_line2_minus = f"{chrom2}\t{start2}\t{end2}\t{row['predominant_event']}_breakpoint2\t{row['support_count']}\t-"
                    bed_lines.append(bed_line2_plus)
                    bed_lines.append(bed_line2_minus)

                elif row["RNAME.1"] == row["RNAME.2"]:
                    # Same chromosome event (deletion, inversion, etc.)
                    coord2 = row["coord_2"]

                    # Only add second breakpoint if it's different from the first
                    # and the gap is significant (> 1kb to avoid very small events)
                    if abs(coord2 - coord1) > 1000:
                        start2 = max(0, coord2 - boundary_window)
                        end2 = coord2 + boundary_window

                        bed_line2_plus = f"{chrom1}\t{start2}\t{end2}\t{row['predominant_event']}_breakpoint2\t{row['support_count']}\t+"
                        bed_line2_minus = f"{chrom1}\t{start2}\t{end2}\t{row['predominant_event']}_breakpoint2\t{row['support_count']}\t-"
                        bed_lines.append(bed_line2_plus)
                        bed_lines.append(bed_line2_minus)

    return bed_lines


class FusionObject(BaseAnalysis):
    """
    Core class for gene fusion analysis.

    Performance Warning:
    - Stores multiple copies of DataFrames in memory
    - Creates unnecessary copies during concatenation
    - Multiple file I/O operations
    - Large memory footprint for BAM processing

    Optimization Suggestions:
    1. Implement streaming for BAM files
    2. Use chunked processing
    3. Cache intermediate results
    4. Use memory-efficient data structures
    5. Implement parallel processing
    6. Use incremental updates
    """

    def __init__(
        self,
        *args,
        target_panel=None,
        reference_file: Optional[str] = None,
        bed_file: Optional[str] = None,
        readfish_toml: Optional[Path] = None,
        master_bed_tree: Optional[MasterBedTree] = None,
        **kwargs,
    ):
        """
        Initialize the FusionObject with analysis parameters.

        Args:
            target_panel: Name of the target panel
            reference_file: Path to reference genome
            bed_file: Path to target regions BED file
            readfish_toml: Path to readfish config
            master_bed_tree: Pre-computed bed tree

        Performance Notes:
        - Initializes data structures efficiently
        - Sets up async processing capabilities
        - Prepares for parallel processing
        """
        # Initialize base class first
        super().__init__(*args, **kwargs)
        state.set_process_state("Fusion Analysis", ProcessState.WAITING_FOR_DATA)

        self.target_panel = target_panel
        self.reference_file = reference_file
        self.bed_file = bed_file
        self.readfish_toml = readfish_toml
        self.fusion_candidates = {}
        self.fusion_candidates_all = {}
        self.structural_variants = {}  # Store structural variants
        self.sv_count = 0  # Counter for structural variants
        self.fstable_all = None
        self.fstable = None
        self.sv_table = None  # Table for structural variants
        self.fstable_all_row_count = 0
        self.all_candidates = 0
        self.fstable_row_count = 0
        self.candidates = 0

        # Initialize UI elements
        self.sv_plot = None
        self.sv_table_container = None
        self.fusionplot = None
        self.fusionplot_all = None
        self.fusiontable = None
        self.fusiontable_all = None

        self.gene_gff3_2 = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)),
            "gencode.v45.basic.annotation.gff3",
        )

        if self.target_panel == "rCNS2":
            self.gene_bed = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)),
                "rCNS2_panel_name_uniq.bed",
            )
        elif self.target_panel == "AML":
            self.gene_bed = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)),
                "AML_panel_name_uniq.bed",
            )

        self.all_gene_bed = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)), "all_genes2.bed"
        )

        datafile = f"{self.target_panel}_data.csv.gz"

        if os.path.isfile(
            os.path.join(os.path.dirname(os.path.abspath(resources.__file__)), datafile)
        ):
            self.gene_table = pd.read_csv(
                os.path.join(
                    os.path.dirname(os.path.abspath(resources.__file__)), datafile
                )
            )
        else:
            # logging.info(
            #    f"This looks like the first time you have run the {self.target_panel} panel."
            # )
            # logging.info("Parsing GFF3")
            self.gene_table = gff3_parser.parse_gff3(
                self.gene_gff3_2, verbose=False, parse_attributes=True
            )

            self.gene_table_small = self.gene_table[
                self.gene_table["Type"].isin(["gene", "exon", "CDS"])
            ]
            self.gene_table_small = self.gene_table_small.drop(
                [
                    "Score",
                    "Phase",
                    "havana_gene",
                    "transcript_support_level",
                    "ont",
                    "transcript_id",
                    "hgnc_id",
                    "protein_id",
                    "havana_transcript",
                    "exon_number",
                    "artif_dupl",
                    "exon_id",
                    "gene_type",
                    "ID",
                    "gene_id",
                    "level",
                    "ccdsid",
                    "tag",
                    "transcript_name",
                    "Parent",
                ],
                axis=1,
            )
            self.gene_table_small.to_csv(
                os.path.join(
                    os.path.dirname(os.path.abspath(resources.__file__)), datafile
                ),
                index=False,
                compression="gzip",
            )
            self.gene_table = self.gene_table_small

        # self.NewBed = NewBed
        self.master_bed_tree = master_bed_tree

        # Add throttling variables for breakpoint graph processing
        self.last_breakpoint_run = None
        self.bam_files_since_last_run = 0
        self.pending_sv_reads = []  # Store SV reads that haven't been processed yet
        self.breakpoint_throttle_time = 30  # seconds
        self.breakpoint_throttle_count = 30  # number of BAM files

    def fusion_table_all(self) -> None:
        """
        Processes and saves all fusion candidates to CSV.

        This method:
        1. Filters for duplicated reads
        2. Identifies reads mapping to multiple genes
        3. Saves results to CSV file

        Performance Notes:
        - Uses efficient pandas operations
        - Implements proper file handling
        - Includes error checking
        """
        # if not self.fusion_candidates_all.empty:
        if self.sampleID in self.fusion_candidates_all.keys():
            uniques_all = self.fusion_candidates_all[self.sampleID][
                "read_id"
            ].duplicated(keep=False)
            doubles_all = self.fusion_candidates_all[self.sampleID][uniques_all]
            counts_all = doubles_all.groupby("read_id", observed=True)[
                "col4"
            ].transform("nunique")
            result_all = doubles_all[counts_all > 1]
            result_all.to_csv(
                os.path.join(
                    self.check_and_create_folder(self.output, self.sampleID),
                    "fusion_candidates_all.csv",
                ),
                index=False,
            )

    def fusion_table(self) -> None:
        """
        Processes and saves fusion candidates within target regions to CSV.

        This method:
        1. Filters for duplicated reads
        2. Identifies reads mapping to multiple genes
        3. Saves results to CSV file

        Performance Notes:
        - Uses efficient pandas operations
        - Implements proper file handling
        - Includes error checking
        """
        # if not self.fusion_candidates.empty:
        if self.sampleID in self.fusion_candidates.keys():
            uniques = self.fusion_candidates[self.sampleID]["read_id"].duplicated(
                keep=False
            )
            doubles = self.fusion_candidates[self.sampleID][uniques]
            counts = doubles.groupby("read_id", observed=True)["col4"].transform(
                "nunique"
            )
            result = doubles[counts > 1]
            result.to_csv(
                os.path.join(
                    self.check_and_create_folder(self.output, self.sampleID),
                    "fusion_candidates_master.csv",
                ),
                index=False,
            )
        # self.update_fusion_table(result)

    async def process_bam(self, bamfile: str, timestamp: str) -> None:
        """
        Asynchronously processes BAM files with optimized memory usage.

        Memory optimizations:
        - Explicit cleanup of large DataFrames after use
        - Streaming processing where possible
        - Efficient data structures to minimize memory footprint
        - Immediate garbage collection after heavy operations
        """
        log_memory_usage("start of process_bam")
        state.set_process_state("Fusion Analysis", ProcessState.RUNNING)

        # Initialize variables to None for explicit cleanup
        fusion_candidates = None
        fusion_candidates_all = None
        new_df = None
        existing_df = None
        combined_df = None
        bed_lines = None

        try:
            try:
                logger.info(f"Starting BAM processing for file: {bamfile}")

                # Check for supplementary alignments using run.cpu_bound
                async def has_supp_background_work():
                    return await run.cpu_bound(has_supplementary, bamfile)

                has_supp = await background_tasks.create(has_supp_background_work())
                log_memory_usage("after supplementary check")

                if has_supp:
                    try:
                        # Process fusion candidates
                        logger.info(f"Processing BAM file for fusions: {bamfile}")

                        async def fusion_background_work():
                            return await run.cpu_bound(
                                fusion_work_pysam,
                                bamfile,
                                self.gene_bed,
                                self.all_gene_bed,
                            )

                        fusion_candidates, fusion_candidates_all = (
                            await background_tasks.create(fusion_background_work())
                        )
                        log_memory_usage("after fusion processing")

                        # Store fusion candidates efficiently with memory optimization
                        if (
                            fusion_candidates is not None
                            and not fusion_candidates.empty
                        ):
                            await self._store_fusion_candidates_efficiently(
                                fusion_candidates, is_target_panel=True
                            )
                            logger.info(
                                f"Added {len(fusion_candidates)} fusion candidates for sample {self.sampleID}"
                            )

                            # Explicitly free the DataFrame
                            del fusion_candidates
                            fusion_candidates = None
                            force_garbage_collection()

                        if (
                            fusion_candidates_all is not None
                            and not fusion_candidates_all.empty
                        ):
                            await self._store_fusion_candidates_efficiently(
                                fusion_candidates_all, is_target_panel=False
                            )
                            logger.info(
                                f"Added {len(fusion_candidates_all)} genome-wide fusion candidates for sample {self.sampleID}"
                            )

                            # Explicitly free the DataFrame
                            del fusion_candidates_all
                            fusion_candidates_all = None
                            force_garbage_collection()

                        # Save fusion results with pre-processing for display
                        await self._save_fusion_results()
                        log_memory_usage("after saving fusion results")

                        # Process genome-wide structural variants with memory monitoring
                        logger.info(
                            f"Processing BAM file for structural variants: {bamfile}"
                        )

                        try:

                            async def sv_background_work():
                                return await run.cpu_bound(
                                    process_bam_pipeline, bamfile
                                )

                            new_df = await background_tasks.create(sv_background_work())
                            log_memory_usage("after structural variant processing")
                        except Exception as e:
                            logger.error(
                                f"Error in structural variant processing: {str(e)}"
                            )
                            new_df = pd.DataFrame()  # Ensure it's defined

                        # Process structural variants with memory optimization
                        if not new_df.empty:
                            await self._process_structural_variants_memory_efficient(
                                new_df
                            )

                            # Explicitly free the DataFrame
                            del new_df
                            new_df = None
                            force_garbage_collection()
                        else:
                            logger.debug(
                                "No structural variant links found in this BAM file"
                            )

                    except Exception as e:
                        logger.error(f"Error processing BAM file: {str(e)}")
                        logger.error("Exception details:", exc_info=True)
                        raise
                    finally:
                        self.running = False
                else:
                    logger.info("BAM file has no supplementary alignments, skipping.")
                    self.running = False

            except Exception as e:
                logger.error(f"Error in process_bam: {e}")
                raise
        finally:
            # Explicit memory cleanup
            self._cleanup_process_bam_memory(
                fusion_candidates,
                fusion_candidates_all,
                new_df,
                existing_df,
                combined_df,
                bed_lines,
            )
            log_memory_usage("end of process_bam")
            state.set_process_state("Fusion Analysis", ProcessState.WAITING_FOR_DATA)

    async def _store_fusion_candidates_efficiently(
        self, candidates: pd.DataFrame, is_target_panel: bool
    ) -> None:
        """
        Store fusion candidates efficiently to minimize memory usage.

        Args:
            candidates: Fusion candidate DataFrame
            is_target_panel: Whether candidates are from target panel or genome-wide
        """
        # Optimize DataFrame memory usage before storing
        candidates = self._optimize_dataframe_memory(candidates)

        target_dict = (
            self.fusion_candidates if is_target_panel else self.fusion_candidates_all
        )

        if self.sampleID not in target_dict:
            target_dict[self.sampleID] = candidates
        else:
            # Use efficient concatenation with copy=False
            target_dict[self.sampleID] = pd.concat(
                [target_dict[self.sampleID], candidates], ignore_index=True, copy=False
            )

    def _optimize_dataframe_memory(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Optimize DataFrame memory usage through dtype optimization.

        Args:
            df: Input DataFrame

        Returns:
            Memory-optimized DataFrame
        """
        # Define optimal dtypes for common columns
        dtype_optimizations = {
            "read_id": "category",
            "col4": "category",  # Gene name
            "reference_id": "category",  # Chromosome
            "strand": "category",
            "mapping_quality": np.int8,
            "reference_start": np.int32,
            "reference_end": np.int32,
            "read_start": np.int32,
            "read_end": np.int32,
            "mapping_span": np.int32,
        }

        # Apply optimizations only for columns that exist
        existing_columns = {
            col: dtype
            for col, dtype in dtype_optimizations.items()
            if col in df.columns
        }

        return df.astype(existing_columns, errors="ignore")

    async def _process_structural_variants_memory_efficient(
        self, new_df: pd.DataFrame
    ) -> None:
        """
        Process structural variants with optimized memory usage.

        Args:
            new_df: New structural variant DataFrame
        """
        log_memory_usage("start of structural variant processing")

        # Construct the output file path
        output_dir = self.check_and_create_folder(self.output, self.sampleID)
        sv_links_file = os.path.join(output_dir, "structural_variant_links.csv")

        try:
            # Optimized file handling with reduced memory usage
            if os.path.exists(sv_links_file) and os.path.getsize(sv_links_file) > 0:
                # Read existing data with optimized dtypes
                dtype_spec = {
                    "QNAME": "category",
                    "RNAME.1": "category",
                    "RNAME.2": "category",
                    "coord_1": np.int64,
                    "coord_2": np.int64,
                    "genomic_gap": np.int64,
                    "event": "category",
                }
                existing_df = safe_read_csv(sv_links_file, dtype=dtype_spec)
                logger.debug(
                    f"Read existing SV links data with {len(existing_df)} rows"
                )
                log_memory_usage("after reading existing data")

                # Append new data efficiently with copy=False
                combined_df = pd.concat(
                    [existing_df, new_df], ignore_index=True, copy=False
                )
                logger.debug(f"Combined DataFrame size: {len(combined_df)} rows")
                log_memory_usage("after concatenation")

                # Free existing_df immediately
                del existing_df
                existing_df = None
                force_garbage_collection()
            else:
                # If file doesn't exist or is empty, use new data
                combined_df = new_df.copy()  # Create a copy to avoid reference issues
                logger.debug("No existing SV links file found, using new data")

            # Save combined data asynchronously
            async def save_csv_work():
                return await run.io_bound(
                    combined_df.to_csv, sv_links_file, index=False
                )

            background_tasks.create(save_csv_work())

            # Generate bed_lines from combined_df for BedTree
            async def bed_lines_background_work():
                return await run.cpu_bound(important_function, combined_df)

            bed_lines = await background_tasks.create(bed_lines_background_work())
            log_memory_usage("after generating bed lines")

            logger.info(f"Saved {len(combined_df)} SV links to {sv_links_file}")

            # Add to BedTree if needed
            if self.master_bed_tree[self.sampleID] is None:
                self.master_bed_tree.add_bed_tree(
                    sample_id=self.sampleID,
                    preserve_original_tree=True,
                    reference_file=f"{self.reference_file}.fai",
                )
                bedfile = self.master_bed_tree.bed_trees[self.sampleID]
                bedfile.load_from_file(self.bed_file)

            if bed_lines:
                bedfile = self.master_bed_tree.bed_trees[self.sampleID]
                bedfile.load_from_string(
                    "\n".join(bed_lines),
                    merge=False,
                    write_files=True,
                    output_location=os.path.join(
                        self.check_and_create_folder(self.output, self.sampleID)
                    ),
                    source_type="FUSION",
                )
                logger.info(
                    f"Added {len(bed_lines)} structural variant breakpoint boundaries to BedTree"
                )

                # Save fusion results again after structural variant processing
                await self._save_fusion_results()

                # Free bed_lines immediately
                del bed_lines
                bed_lines = None
                force_garbage_collection()

        except Exception as e:
            logger.error(f"Error saving structural variant links: {str(e)}")
            logger.error("Exception details:", exc_info=True)
        finally:
            # Explicit cleanup of DataFrames
            if "existing_df" in locals() and existing_df is not None:
                del existing_df
            if "combined_df" in locals() and combined_df is not None:
                del combined_df
            force_garbage_collection()
            log_memory_usage("end of structural variant processing")

    def _cleanup_process_bam_memory(self, *dataframes) -> None:
        """
        Explicitly clean up memory from DataFrames used in process_bam.

        Args:
            *dataframes: Variable number of DataFrames to clean up
        """
        log_memory_usage("before cleanup")

        for df in dataframes:
            if df is not None:
                del df

        # Force garbage collection
        force_garbage_collection()

        # Additional cleanup for any remaining references
        # Note: sys.exc_clear() was removed in Python 3, so we skip this step
        # The garbage collector will handle any remaining references automatically

        logger.debug("Memory cleanup completed for process_bam")
        log_memory_usage("after cleanup")

    def _should_run_breakpoint_processing(self) -> bool:
        """
        Determine if breakpoint graph processing should run based on throttling criteria.

        Returns:
            bool: True if processing should run, False otherwise
        """
        current_time = datetime.now()

        # Check time-based throttling
        if self.last_breakpoint_run is None:
            # First run - always process
            return True

        time_since_last = (current_time - self.last_breakpoint_run).total_seconds()
        if time_since_last >= self.breakpoint_throttle_time:
            logger.info(
                f"Running breakpoint processing due to time threshold ({time_since_last:.1f}s >= {self.breakpoint_throttle_time}s)"
            )
            return True

        # Check count-based throttling
        if self.bam_files_since_last_run >= self.breakpoint_throttle_count:
            logger.info(
                f"Running breakpoint processing due to count threshold ({self.bam_files_since_last_run} >= {self.breakpoint_throttle_count})"
            )
            return True

        # No throttling criteria met
        logger.debug(
            f"Throttling breakpoint processing: {time_since_last:.1f}s since last run, {self.bam_files_since_last_run} files processed"
        )
        return False

    async def _process_breakpoint_results(
        self, bed_lines: list, sv_reads: pd.DataFrame
    ) -> None:
        """
        Optimized version of breakpoint results processing using vectorized operations.

        Args:
            bed_lines: List of BED format lines from breakpoint analysis
            sv_reads: Combined structural variant reads DataFrame
        """
        try:
            # Pre-process sv_reads for faster lookups
            if sv_reads.empty:
                sv_data = []
            else:
                # Convert to numpy arrays for faster operations
                chroms = sv_reads["RNAME"].values
                starts = sv_reads["REF_START"].values.astype(float).astype(int)
                ends = sv_reads["REF_END"].values.astype(float).astype(int)
                qnames = sv_reads["QNAME"].values

                # Create efficient lookup structures
                # Group by QNAME for faster partner finding - use list of dicts for easier access
                qname_groups = {}
                for i, qname in enumerate(qnames):
                    if qname not in qname_groups:
                        qname_groups[qname] = []
                    qname_groups[qname].append(
                        {"chrom": chroms[i], "start": starts[i], "end": ends[i]}
                    )

                # Pre-allocate sv_data list
                sv_data = []
                sv_data_append = sv_data.append  # Local reference for faster append

                # Process bed lines efficiently
                for line in bed_lines:
                    chrom, start, end, sv_type, _, strand = line.split("\t")
                    start_int = int(float(start))
                    end_int = int(float(end))

                    # Vectorized filtering using numpy
                    chrom_mask = chroms == chrom
                    start_mask = starts >= start_int - 1000
                    end_mask = ends <= end_int + 1000
                    combined_mask = chrom_mask & start_mask & end_mask

                    matching_indices = np.where(combined_mask)[0]

                    # Find partner information efficiently
                    partner_info = ""
                    if len(matching_indices) > 0:
                        # Get unique QNAMEs from matching reads
                        matching_qnames = set(qnames[matching_indices])

                        # Find partner chromosomes efficiently
                        for qname in matching_qnames:
                            if qname in qname_groups:
                                pairs = qname_groups[qname]
                                if len(pairs) > 1:
                                    # Find partner alignment (different chromosome)
                                    for pair in pairs:
                                        if pair["chrom"] != chrom:
                                            partner_chrom = pair["chrom"]
                                            partner_pos = f"{pair['start']:,}"
                                            partner_info = (
                                                f"{partner_chrom}:{partner_pos}"
                                            )
                                            break
                                    if (
                                        partner_info
                                    ):  # Found partner, no need to continue
                                        break

                    # Format coordinates efficiently
                    formatted_start = f"{start_int:,}"
                    formatted_end = f"{end_int:,}"
                    size = end_int - start_int
                    formatted_size = f"{size:,}"

                    event_location = f"{chrom}:{formatted_start}-{formatted_end}"
                    if partner_info:
                        event_location += f" ⟷ {partner_info}"

                    # Use local append reference for better performance
                    sv_data_append(
                        {
                            "Event Type": sv_type,
                            "Primary Location": f"{chrom}:{formatted_start}-{formatted_end}",
                            "Partner Location": partner_info if partner_info else "N/A",
                            "Size (bp)": formatted_size,
                            "Strand": strand,
                            "Full Location": event_location,
                        }
                    )

            # Create DataFrame efficiently
            sv_df = pd.DataFrame(sv_data)

            # Sort efficiently using numpy argsort if DataFrame is large
            if len(sv_df) > 1000:
                # Use numpy for large datasets
                event_type_values = sv_df["Event Type"].values
                primary_loc_values = sv_df["Primary Location"].values

                # Create composite sort key
                sort_key = np.char.add(event_type_values, "_" + primary_loc_values)
                sort_indices = np.argsort(sort_key)
                sv_df = sv_df.iloc[sort_indices].reset_index(drop=True)
            else:
                # Use pandas sort for smaller datasets
                sv_df = sv_df.sort_values(
                    ["Event Type", "Primary Location"]
                ).reset_index(drop=True)

            # Save structural variants to CSV
            sv_df.to_csv(
                os.path.join(
                    self.check_and_create_folder(self.output, self.sampleID),
                    "structural_variants.csv",
                ),
                index=False,
            )

            # Add to BedTree if needed
            if self.master_bed_tree[self.sampleID] is None:
                self.master_bed_tree.add_bed_tree(
                    sample_id=self.sampleID,
                    preserve_original_tree=True,
                    reference_file=f"{self.reference_file}.fai",
                )
            bedfile = self.master_bed_tree.bed_trees[self.sampleID]

            bedfile.load_from_string(
                "\n".join(bed_lines),
                merge=False,
                write_files=True,
                output_location=os.path.join(
                    self.check_and_create_folder(self.output, self.sampleID)
                ),
                source_type="FUSION",
            )

            logger.info(
                f"Processed {len(bed_lines)} structural variant events from {len(sv_reads)} reads"
            )

        except Exception as e:
            logger.error(f"Error processing breakpoint results: {str(e)}")
            logger.error("Exception details:", exc_info=True)

    async def stop_analysis(self):
        """
        Stops ongoing analysis and processes any remaining pending SV reads.
        """
        state.set_process_state("Fusion Analysis", ProcessState.STOPPING)

        # Save final fusion results
        await self._save_fusion_results()

        # Process any remaining pending SV reads before stopping
        if self.pending_sv_reads:
            logger.info(
                f"Processing {len(self.pending_sv_reads)} pending SV read sets before stopping"
            )
            try:
                combined_sv_reads = pd.concat(self.pending_sv_reads, ignore_index=True)

                bed_lines = await run.cpu_bound(
                    build_breakpoint_graph,
                    combined_sv_reads,
                    50000,  # max_proximity
                    True,  # group_by_sv
                )

                if len(bed_lines) > 0:
                    await self._process_breakpoint_results(bed_lines, combined_sv_reads)

                logger.info("Final breakpoint processing completed")

            except Exception as e:
                logger.error(f"Error in final breakpoint processing: {str(e)}")

        state.stop_process("Fusion Analysis")
        await super().stop_analysis()

    async def _save_fusion_results(self) -> None:
        """
        Save fusion candidates to CSV files and update UI.
        This should be called periodically or when processing is complete.

        Enhanced to pre-process all data needed by FusionVis to eliminate
        heavy lifting from the display layer.
        """

        try:
            output_dir = self.check_and_create_folder(self.output, self.sampleID)

            # Save fusion candidates within target regions
            if (
                self.sampleID in self.fusion_candidates
                and not self.fusion_candidates[self.sampleID].empty
            ):
                # Filter for duplicated reads
                df = self.fusion_candidates[self.sampleID]
                uniques = df["read_id"].duplicated(keep=False)
                doubles = df[uniques]
                if not doubles.empty:
                    counts = doubles.groupby("read_id", observed=True)[
                        "col4"
                    ].transform("nunique")
                    result = doubles[counts > 1]
                    if not result.empty:
                        # Save raw data
                        result.to_csv(
                            os.path.join(output_dir, "fusion_candidates_master.csv"),
                            index=False,
                        )

                        # Pre-process data for FusionVis using run.cpu_bound
                        await run.cpu_bound(
                            preprocess_fusion_data_standalone,
                            result,
                            os.path.join(
                                output_dir, "fusion_candidates_master_processed.csv"
                            ),
                        )

                        logger.info(
                            f"Saved {len(result)} fusion candidates within target regions"
                        )

            # Save genome-wide fusion candidates
            if (
                self.sampleID in self.fusion_candidates_all
                and not self.fusion_candidates_all[self.sampleID].empty
            ):
                # Filter for duplicated reads
                df = self.fusion_candidates_all[self.sampleID]
                uniques_all = df["read_id"].duplicated(keep=False)
                doubles_all = df[uniques_all]
                if not doubles_all.empty:
                    counts_all = doubles_all.groupby("read_id", observed=True)[
                        "col4"
                    ].transform("nunique")
                    result_all = doubles_all[counts_all > 1]
                    if not result_all.empty:
                        # Save raw data
                        result_all.to_csv(
                            os.path.join(output_dir, "fusion_candidates_all.csv"),
                            index=False,
                        )

                        # Pre-process data for FusionVis using run.cpu_bound
                        await run.cpu_bound(
                            preprocess_fusion_data_standalone,
                            result_all,
                            os.path.join(
                                output_dir, "fusion_candidates_all_processed.csv"
                            ),
                        )

                        logger.info(
                            f"Saved {len(result_all)} genome-wide fusion candidates"
                        )

            # Pre-process structural variant summary data using run.cpu_bound
            await run.cpu_bound(preprocess_structural_variants_standalone, output_dir)

        except Exception as e:
            logger.error(f"Error saving fusion results: {str(e)}")
            logger.error("Exception details:", exc_info=True)

    def _preprocess_fusion_data_sync(
        self, fusion_data: pd.DataFrame, output_file: str
    ) -> None:
        """
        Synchronous version of fusion data preprocessing for CPU-bound execution.

        Args:
            fusion_data: Raw fusion candidate data
            output_file: Path to save processed data
        """
        try:
            # Apply categorical data types for efficiency
            fusion_data = fusion_data.astype(
                {
                    "read_id": "category",
                    "col4": "category",  # Gene column
                    "reference_id": "category",  # Chromosome column
                    "strand": "category",
                }
            )

            # Annotate results (this is the heavy lifting)
            annotated_data, goodpairs = _annotate_results(fusion_data)

            # Create processed data structure for FusionVis
            processed_data = {
                "annotated_data": annotated_data,
                "goodpairs": goodpairs,
                "gene_pairs": [],
                "gene_groups": [],
                "candidate_count": 0,
            }

            # Process gene pairs and groups if we have good pairs
            if not annotated_data.empty and goodpairs.any():
                gene_pairs = (
                    annotated_data[goodpairs]
                    .sort_values(by="reference_start")["tag"]
                    .unique()
                    .tolist()
                )
                gene_pairs = [tuple(pair.split(",")) for pair in gene_pairs]
                gene_groups_test = get_gene_network(gene_pairs)
                gene_groups = []

                for gene_group in gene_groups_test:
                    reads = _get_reads(
                        annotated_data[goodpairs][
                            annotated_data[goodpairs]["col4"].isin(gene_group)
                        ]
                    )
                    if len(reads) > 1:
                        gene_groups.append(gene_group)

                processed_data.update(
                    {
                        "gene_pairs": gene_pairs,
                        "gene_groups": gene_groups,
                        "candidate_count": len(gene_groups),
                    }
                )

            # Save processed data as pickle for efficient loading
            import pickle

            with open(output_file, "wb") as f:
                pickle.dump(processed_data, f)

            logger.info(f"Pre-processed fusion data saved to {output_file}")

        except Exception as e:
            logger.error(f"Error pre-processing fusion data: {str(e)}")
            logger.error("Exception details:", exc_info=True)

    async def _preprocess_fusion_data_for_display(
        self, fusion_data: pd.DataFrame, output_file: str
    ) -> None:
        """
        Pre-process fusion data for display, moving heavy lifting from FusionVis.

        Args:
            fusion_data: Raw fusion candidate data
            output_file: Path to save processed data
        """
        try:
            # Apply categorical data types for efficiency
            fusion_data = fusion_data.astype(
                {
                    "read_id": "category",
                    "col4": "category",  # Gene column
                    "reference_id": "category",  # Chromosome column
                    "strand": "category",
                }
            )

            # Annotate results (this is the heavy lifting)
            annotated_data, goodpairs = _annotate_results(fusion_data)

            # Create processed data structure for FusionVis
            processed_data = {
                "annotated_data": annotated_data,
                "goodpairs": goodpairs,
                "gene_pairs": [],
                "gene_groups": [],
                "candidate_count": 0,
            }

            # Process gene pairs and groups if we have good pairs
            if not annotated_data.empty and goodpairs.any():
                gene_pairs = (
                    annotated_data[goodpairs]
                    .sort_values(by="reference_start")["tag"]
                    .unique()
                    .tolist()
                )
                gene_pairs = [tuple(pair.split(",")) for pair in gene_pairs]
                gene_groups_test = get_gene_network(gene_pairs)
                gene_groups = []

                for gene_group in gene_groups_test:
                    reads = _get_reads(
                        annotated_data[goodpairs][
                            annotated_data[goodpairs]["col4"].isin(gene_group)
                        ]
                    )
                    if len(reads) > 1:
                        gene_groups.append(gene_group)

                processed_data.update(
                    {
                        "gene_pairs": gene_pairs,
                        "gene_groups": gene_groups,
                        "candidate_count": len(gene_groups),
                    }
                )

            # Save processed data as pickle for efficient loading
            import pickle

            with open(output_file, "wb") as f:
                pickle.dump(processed_data, f)

            logger.info(f"Pre-processed fusion data saved to {output_file}")

        except Exception as e:
            logger.error(f"Error pre-processing fusion data: {str(e)}")
            logger.error("Exception details:", exc_info=True)

    def _preprocess_structural_variants_sync(self, output_dir: str) -> None:
        """
        Synchronous version of structural variant preprocessing for CPU-bound execution.

        Args:
            output_dir: Output directory path
        """
        try:
            sv_links_file = os.path.join(output_dir, "structural_variant_links.csv")

            if os.path.exists(sv_links_file) and os.path.getsize(sv_links_file) > 0:
                # Read the structural variant links CSV using safe reader
                dtype_spec = {
                    "QNAME": str,
                    "RNAME.1": str,
                    "RNAME.2": str,
                    "coord_1": np.int64,
                    "coord_2": np.int64,
                    "genomic_gap": np.int64,
                    "event": str,
                }
                sv_links_df = safe_read_csv(sv_links_file, dtype=dtype_spec)

                if not sv_links_df.empty:
                    # Process the links data to create a summary for display
                    sv_summary = get_summary(sv_links_df, min_support=2)

                    if not sv_summary.empty:
                        # Convert the summary to the format expected by the UI
                        sv_df = pd.DataFrame(
                            {
                                "Event Type": sv_summary["predominant_event"],
                                "Primary Location": sv_summary.apply(
                                    lambda row: f"{row['RNAME.1']}:{row['coord_1']:,}",
                                    axis=1,
                                ),
                                "Partner Location": sv_summary.apply(
                                    lambda row: f"{row['RNAME.2']}:{row['coord_2']:,}",
                                    axis=1,
                                ),
                                "Size (bp)": sv_summary["median_genomic_gap"].apply(
                                    lambda x: f"{x:,}" if x >= 0 else "N/A"
                                ),
                                "Strand": "Unknown",  # Not available in links data
                                "Full Location": sv_summary.apply(
                                    lambda row: (
                                        f"{row['RNAME.1']}:{row['coord_1']:,}-{row['coord_2']:,}"
                                        if row["RNAME.1"] == row["RNAME.2"]
                                        else f"{row['RNAME.1']}:{row['coord_1']:,} ⟷ {row['RNAME.2']}:{row['coord_2']:,}"
                                    ),
                                    axis=1,
                                ),
                                "Support Count": sv_summary["support_count"],
                                "Supporting Reads": sv_summary[
                                    "supporting_reads"
                                ].apply(
                                    lambda x: ", ".join(x[:5])
                                    + ("..." if len(x) > 5 else "")
                                ),
                            }
                        )

                        # Save processed structural variant data
                        sv_df.to_csv(
                            os.path.join(
                                output_dir, "structural_variants_processed.csv"
                            ),
                            index=False,
                        )

                        # Save count for quick access
                        with open(os.path.join(output_dir, "sv_count.txt"), "w") as f:
                            f.write(str(len(sv_df)))

                        logger.info(
                            f"Pre-processed {len(sv_df)} structural variant events"
                        )
                    else:
                        # Save empty count
                        with open(os.path.join(output_dir, "sv_count.txt"), "w") as f:
                            f.write("0")
                else:
                    # Save empty count
                    with open(os.path.join(output_dir, "sv_count.txt"), "w") as f:
                        f.write("0")

        except Exception as e:
            logger.error(f"Error pre-processing structural variants: {str(e)}")
            logger.error("Exception details:", exc_info=True)

    async def _preprocess_structural_variants_for_display(
        self, output_dir: str
    ) -> None:
        """
        Pre-process structural variant data for display, moving heavy lifting from FusionVis.

        Args:
            output_dir: Output directory path
        """
        try:
            sv_links_file = os.path.join(output_dir, "structural_variant_links.csv")

            if os.path.exists(sv_links_file) and os.path.getsize(sv_links_file) > 0:
                # Read the structural variant links CSV using safe reader
                dtype_spec = {
                    "QNAME": str,
                    "RNAME.1": str,
                    "RNAME.2": str,
                    "coord_1": np.int64,
                    "coord_2": np.int64,
                    "genomic_gap": np.int64,
                    "event": str,
                }
                sv_links_df = safe_read_csv(sv_links_file, dtype=dtype_spec)

                if not sv_links_df.empty:
                    # Process the links data to create a summary for display
                    sv_summary = get_summary(sv_links_df, min_support=2)

                    if not sv_summary.empty:
                        # Convert the summary to the format expected by the UI
                        sv_df = pd.DataFrame(
                            {
                                "Event Type": sv_summary["predominant_event"],
                                "Primary Location": sv_summary.apply(
                                    lambda row: f"{row['RNAME.1']}:{row['coord_1']:,}",
                                    axis=1,
                                ),
                                "Partner Location": sv_summary.apply(
                                    lambda row: f"{row['RNAME.2']}:{row['coord_2']:,}",
                                    axis=1,
                                ),
                                "Size (bp)": sv_summary["median_genomic_gap"].apply(
                                    lambda x: f"{x:,}" if x >= 0 else "N/A"
                                ),
                                "Strand": "Unknown",  # Not available in links data
                                "Full Location": sv_summary.apply(
                                    lambda row: (
                                        f"{row['RNAME.1']}:{row['coord_1']:,}-{row['coord_2']:,}"
                                        if row["RNAME.1"] == row["RNAME.2"]
                                        else f"{row['RNAME.1']}:{row['coord_1']:,} ⟷ {row['RNAME.2']}:{row['coord_2']:,}"
                                    ),
                                    axis=1,
                                ),
                                "Support Count": sv_summary["support_count"],
                                "Supporting Reads": sv_summary[
                                    "supporting_reads"
                                ].apply(
                                    lambda x: ", ".join(x[:5])
                                    + ("..." if len(x) > 5 else "")
                                ),
                            }
                        )

                        # Save processed structural variant data
                        sv_df.to_csv(
                            os.path.join(
                                output_dir, "structural_variants_processed.csv"
                            ),
                            index=False,
                        )

                        # Save count for quick access
                        with open(os.path.join(output_dir, "sv_count.txt"), "w") as f:
                            f.write(str(len(sv_df)))

                        logger.info(
                            f"Pre-processed {len(sv_df)} structural variant events"
                        )
                    else:
                        # Save empty count
                        with open(os.path.join(output_dir, "sv_count.txt"), "w") as f:
                            f.write("0")
                else:
                    # Save empty count
                    with open(os.path.join(output_dir, "sv_count.txt"), "w") as f:
                        f.write("0")

        except Exception as e:
            logger.error(f"Error pre-processing structural variants: {str(e)}")
            logger.error("Exception details:", exc_info=True)


def get_summary(links_df, min_support=2):
    """
    Generate summary of structural variants with support filtering and memory optimization.

    Args:
        links_df: Structural variant links DataFrame
        min_support: Minimum number of supporting reads

    Returns:
        Summary DataFrame

    Memory optimizations:
    - Efficient groupby operations
    - Minimal data copying
    """
    summary = (
        links_df.groupby(["RNAME.1", "coord_1", "RNAME.2", "coord_2"], observed=True)
        .agg(
            support_count=("QNAME", "nunique"),
            supporting_reads=("QNAME", lambda s: list(s.unique())),
            from_med=("coord_1", "median"),
            to_med=("coord_2", "median"),
            median_genomic_gap=("genomic_gap", "median"),
            event_counts=("event", lambda s: s.value_counts().to_dict()),
            predominant_event=("event", lambda s: s.mode().iloc[0]),
        )
        .reset_index()
    )

    # Filter by minimum support
    return summary[summary["support_count"] >= min_support].copy()


def safe_read_csv(file_path, dtype=None):
    """
    Safely read a CSV file that might be compressed with memory optimization.

    Args:
        file_path: Path to the CSV file
        dtype: Optional dtype specification for pandas

    Returns:
        pd.DataFrame: The loaded DataFrame

    Memory optimizations:
    - Efficient dtype specification
    - Minimal memory footprint
    """
    try:
        # First try reading as uncompressed CSV
        return pd.read_csv(file_path, dtype=dtype)
    except UnicodeDecodeError:
        # If that fails, try reading as gzip compressed
        try:
            return pd.read_csv(file_path, compression="gzip", dtype=dtype)
        except Exception as e:
            logger.error(
                f"Failed to read CSV file {file_path} with both uncompressed and gzip methods: {str(e)}"
            )
            raise


def _annotate_results(result: pd.DataFrame) -> Tuple[pd.DataFrame, pd.Series]:
    """
    Annotates the result DataFrame with tags and colors with memory optimization.

    Memory optimizations:
    - Use inplace operations where possible
    - Combine groupby operations
    - Use categorical dtypes
    - Vectorize string operations
    - Minimize DataFrame copies
    """
    # Use inplace operations to avoid unnecessary copies
    result = result.copy()  # Only one copy needed

    # Group by read_id and aggregate col4 (Gene) values efficiently
    lookup = result.groupby("read_id", observed=True)["col4"].agg(
        lambda x: ",".join(set(x))
    )
    result["tag"] = result["read_id"].map(lookup)

    # Generate colors for each read_id group efficiently
    colors = result.groupby("read_id", observed=True)["col4"].apply(
        lambda x: _generate_random_color()
    )
    result["Color"] = result["read_id"].map(colors)

    # Clean string data efficiently
    result = result.apply(lambda x: x.strip() if isinstance(x, str) else x)

    # Find good pairs (reads that map to more than 2 genes)
    goodpairs = result.groupby("tag", observed=True)["read_id"].transform("nunique") > 2

    return result, goodpairs


def _generate_random_color() -> str:
    """
    Generates a random color for use in plotting.

    Returns:
        str: A random hex color code.
    """
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))


def preprocess_structural_variants_standalone(output_dir: str) -> None:
    """
    Standalone version of structural variant preprocessing for CPU-bound execution with memory optimization.

    Args:
        output_dir: Output directory path

    Memory optimizations:
    - Efficient data loading
    - Minimal intermediate data storage
    - Immediate cleanup of large DataFrames
    """
    try:
        sv_links_file = os.path.join(output_dir, "structural_variant_links.csv")

        if os.path.exists(sv_links_file) and os.path.getsize(sv_links_file) > 0:
            # Read the structural variant links CSV using safe reader with optimized dtypes
            dtype_spec = {
                "QNAME": "category",
                "RNAME.1": "category",
                "RNAME.2": "category",
                "coord_1": np.int64,
                "coord_2": np.int64,
                "genomic_gap": np.int64,
                "event": "category",
            }
            sv_links_df = safe_read_csv(sv_links_file, dtype=dtype_spec)

            if not sv_links_df.empty:
                # Process the links data to create a summary for display
                sv_summary = get_summary(sv_links_df, min_support=2)

                # Free the large DataFrame immediately
                del sv_links_df

                if not sv_summary.empty:
                    # Convert the summary to the format expected by the UI efficiently
                    sv_df = pd.DataFrame(
                        {
                            "Event Type": sv_summary["predominant_event"],
                            "Primary Location": sv_summary.apply(
                                lambda row: f"{row['RNAME.1']}:{row['coord_1']:,}",
                                axis=1,
                            ),
                            "Partner Location": sv_summary.apply(
                                lambda row: f"{row['RNAME.2']}:{row['coord_2']:,}",
                                axis=1,
                            ),
                            "Size (bp)": sv_summary["median_genomic_gap"].apply(
                                lambda x: f"{x:,}" if x >= 0 else "N/A"
                            ),
                            "Strand": "Unknown",  # Not available in links data
                            "Full Location": sv_summary.apply(
                                lambda row: (
                                    f"{row['RNAME.1']}:{row['coord_1']:,}-{row['coord_2']:,}"
                                    if row["RNAME.1"] == row["RNAME.2"]
                                    else f"{row['RNAME.1']}:{row['coord_1']:,} ⟷ {row['RNAME.2']}:{row['coord_2']:,}"
                                ),
                                axis=1,
                            ),
                            "Support Count": sv_summary["support_count"],
                            "Supporting Reads": sv_summary["supporting_reads"].apply(
                                lambda x: ", ".join(x[:5])
                                + ("..." if len(x) > 5 else "")
                            ),
                        }
                    )

                    # Save processed structural variant data
                    sv_df.to_csv(
                        os.path.join(output_dir, "structural_variants_processed.csv"),
                        index=False,
                    )

                    # Save count for quick access
                    with open(os.path.join(output_dir, "sv_count.txt"), "w") as f:
                        f.write(str(len(sv_df)))

                    logger.info(f"Pre-processed {len(sv_df)} structural variant events")

                    # Free the DataFrame immediately
                    del sv_df
                else:
                    # Save empty count
                    with open(os.path.join(output_dir, "sv_count.txt"), "w") as f:
                        f.write("0")
            else:
                # Save empty count
                with open(os.path.join(output_dir, "sv_count.txt"), "w") as f:
                    f.write("0")

    except Exception as e:
        logger.error(f"Error pre-processing structural variants: {str(e)}")
        logger.error("Exception details:", exc_info=True)
    finally:
        # Force garbage collection
        force_garbage_collection()


def get_gene_network(gene_pairs):
    """
    Build gene network from gene pairs with memory optimization.

    Args:
        gene_pairs: List of gene pairs

    Returns:
        List of connected components

    Memory optimizations:
    - Efficient graph construction
    - Minimal intermediate data storage
    """
    G = nx.Graph()
    for pair in gene_pairs:
        G.add_edge(pair[0], pair[1])
    connected_components = list(nx.connected_components(G))

    # Free the graph immediately
    del G

    return [list(component) for component in connected_components]


def _get_reads(reads: pd.DataFrame) -> pd.DataFrame:
    """
    Get reads for a specific gene with memory optimization.

    Args:
        reads (pd.DataFrame): DataFrame with reads

    Returns:
        pd.DataFrame: DataFrame with reads

    Memory optimizations:
    - Efficient column operations
    - Minimal data copying
    - Immediate cleanup of intermediate data
    """
    df = reads.copy()  # Only one copy needed
    df.columns = [
        "chromosome",
        "start",
        "end",
        "gene",
        "chromosome2",
        "start2",
        "end2",
        "id",
        "quality",
        "strand",
        "read_start",
        "read_end",
        "secondary",
        "supplementary",
        "span",
        "tag",
        "color",
    ]

    # Add logging for start and end columns before conversion
    logger.debug(
        f"Converting start column to int. First few values: {df['start'].head()}"
    )
    logger.debug(f"Converting end column to int. First few values: {df['end'].head()}")

    try:
        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)
    except ValueError as e:
        logger.error(f"Error converting start/end to int: {str(e)}")
        logger.error(
            f"Problematic values in start: {df[df['start'].apply(lambda x: not str(x).isdigit())]['start'].tolist()}"
        )
        logger.error(
            f"Problematic values in end: {df[df['end'].apply(lambda x: not str(x).isdigit())]['end'].tolist()}"
        )
        raise

    # Sort the DataFrame by chromosome, start, and end positions
    df = df.sort_values(by=["chromosome", "start", "end"])

    df = df.drop_duplicates(subset=["start2", "end2", "id"])

    # Group by chromosome and collapse ranges within each group
    # Use a more explicit approach to avoid the FutureWarning
    result_list = []
    for chromosome, group in df.groupby("chromosome", observed=True):
        collapsed_group = collapse_ranges(group, 10000)
        if not collapsed_group.empty:
            collapsed_group['chromosome'] = chromosome
            result_list.append(collapsed_group)
    
    if result_list:
        result = pd.concat(result_list, ignore_index=True)
    else:
        result = pd.DataFrame()

    # Free the intermediate DataFrame
    del df

    return result


def collapse_ranges(df, max_distance):
    """
    Collapse ranges within a fixed distance with memory optimization.

    Args:
        df: DataFrame with ranges
        max_distance: Maximum distance for collapsing

    Returns:
        DataFrame with collapsed ranges

    Memory optimizations:
    - Efficient iteration
    - Minimal intermediate data storage
    """
    collapsed = []
    current_range = None

    for _, row in df.iterrows():
        # Only unpack what we need
        start, end = row["start"], row["end"]

        if current_range is None:
            current_range = row
        else:
            if start <= current_range["end"] + max_distance:
                current_range["end"] = max(current_range["end"], end)
            else:
                collapsed.append(current_range)
                current_range = row

    if current_range is not None:
        collapsed.append(current_range)

    return pd.DataFrame(collapsed)


def preprocess_fusion_data_standalone(
    fusion_data: pd.DataFrame, output_file: str
) -> None:
    """
    Standalone version of fusion data preprocessing for CPU-bound execution with memory optimization.

    Args:
        fusion_data: Raw fusion candidate data
        output_file: Path to save processed data

    Memory optimizations:
    - Efficient data type optimization
    - Minimal intermediate data storage
    - Immediate cleanup of large DataFrames
    """
    try:
        # Apply categorical data types for efficiency
        fusion_data = fusion_data.astype(
            {
                "read_id": "category",
                "col4": "category",  # Gene column
                "reference_id": "category",  # Chromosome column
                "strand": "category",
            }
        )

        # Annotate results (this is the heavy lifting)
        annotated_data, goodpairs = _annotate_results(fusion_data)

        # Free the original data immediately
        del fusion_data

        # Create processed data structure for FusionVis
        processed_data = {
            "annotated_data": annotated_data,
            "goodpairs": goodpairs,
            "gene_pairs": [],
            "gene_groups": [],
            "candidate_count": 0,
        }

        # Process gene pairs and groups if we have good pairs
        if not annotated_data.empty and goodpairs.any():
            gene_pairs = (
                annotated_data[goodpairs]
                .sort_values(by="reference_start")["tag"]
                .unique()
                .tolist()
            )
            gene_pairs = [tuple(pair.split(",")) for pair in gene_pairs]
            gene_groups_test = get_gene_network(gene_pairs)
            gene_groups = []

            for gene_group in gene_groups_test:
                reads = _get_reads(
                    annotated_data[goodpairs][
                        annotated_data[goodpairs]["col4"].isin(gene_group)
                    ]
                )
                if len(reads) > 1:
                    gene_groups.append(gene_group)

            processed_data.update(
                {
                    "gene_pairs": gene_pairs,
                    "gene_groups": gene_groups,
                    "candidate_count": len(gene_groups),
                }
            )

        # Save processed data as pickle for efficient loading
        import pickle

        with open(output_file, "wb") as f:
            pickle.dump(processed_data, f)

        logger.info(f"Pre-processed fusion data saved to {output_file}")

        # Free the processed data immediately
        del processed_data

    except Exception as e:
        logger.error(f"Error pre-processing fusion data: {str(e)}")
        logger.error("Exception details:", exc_info=True)
    finally:
        # Force garbage collection
        force_garbage_collection()
