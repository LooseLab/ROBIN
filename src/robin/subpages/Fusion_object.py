"""
Fusion Gene Identification Module

This module provides functionality for identifying gene fusion candidates from BAM files.
The primary components include running external tools to extract and process gene fusion
candidates, visualizing the results using a GUI, and managing the process asynchronously.

Performance Considerations:
- Avoid unnecessary DataFrame copies
- Use inplace operations where possible
- Implement efficient memory management
- Use streaming for large files
- Cache intermediate results
- Use parallel processing where appropriate

Memory Optimization Notes:
1. DataFrame Operations:
   - Avoid .copy() unless necessary
   - Use inplace=True for modifications
   - Use categorical dtypes for string columns
   - Implement chunked processing for large files
   - Use memory-efficient data structures

2. BAM Processing:
   - Stream BAM files instead of loading entirely
   - Process reads in chunks
   - Use efficient data structures for alignments
   - Implement parallel processing
   - Cache intermediate results

3. File I/O:
   - Use streaming for large files
   - Implement incremental processing
   - Use compression for storage
   - Cache frequently accessed data
   - Implement lazy loading

4. Visualization:
   - Use efficient data structures
   - Implement pagination
   - Cache plot generation
   - Use lazy loading for large datasets
   - Implement virtual scrolling
"""

import os
import sys
import subprocess
import gff3_parser
import random
import pysam
import logging
import networkx as nx
import pandas as pd
import numpy as np
import re
import click
import asyncio
from typing import Optional, Tuple, Dict, List, Set
from nicegui import ui, run, app, background_tasks
from robin import theme, resources
from dna_features_viewer import GraphicFeature, GraphicRecord
from pathlib import Path
import matplotlib
from matplotlib.ticker import FuncFormatter
from robin.subpages.base_analysis import BaseAnalysis, BaseVis
from robin.utilities.decompress import decompress_gzip_file
from robin.utilities.bed_file import MasterBedTree
from collections import Counter, defaultdict
from robin.core.state import state, ProcessState
from datetime import datetime
from itertools import combinations

matplotlib.use("agg")
from matplotlib import pyplot as plt

# Use the main logger configured in the main application
logger = logging.getLogger(__name__)

# Configure matplotlib font settings
plt.rcParams["font.family"] = ["sans-serif"]
plt.rcParams["font.sans-serif"] = [
    "SF Pro Text",
    "-apple-system",
    "BlinkMacSystemFont",
    "Helvetica",
    "Arial",
    "sans-serif",
]

# Override DNA Features Viewer default font settings
# GraphicRecord.default_font_family = "sans-serif"
# GraphicFeature.default_font_family = "sans-serif"

os.environ["CI"] = "1"
STRAND = {"+": 1, "-": -1}


decompress_gzip_file(
    os.path.join(
        os.path.dirname(os.path.abspath(resources.__file__)),
        "gencode.v45.basic.annotation.gff3.gz",
    )
)


def get_aligned_read_coords(aln):
    """
    Return (read_start, read_end) on the original read,
    *unsorted*—we'll sort them downstream.
    """
    ct = aln.cigartuples or []
    # 1) leading S/H
    rs = 0
    for op, l in ct:
        if op in (4, 5):  # S or H
            rs += l
        else:
            break
    # 2) trailing S/H
    te = 0
    for op, l in reversed(ct):
        if op in (4, 5):
            te += l
        else:
            break
    re = aln.query_length - te
    return rs, re


def extract_split_read_alignments(bam_path):
    """
    Optimized version of split read alignment extraction with early filtering and memory efficiency.
    
    Performance improvements:
    - Single pass through BAM file instead of two passes
    - Early filtering to reduce memory usage
    - Pre-allocated lists for better memory efficiency
    - Optimized CIGAR parsing
    """
    # Pre-allocate lists with estimated capacity to avoid resizing
    estimated_capacity = 10000  # Conservative estimate
    qnames = []
    types = []
    rnames = []
    strands = []
    ref_starts = []
    ref_ends = []
    ref_spans = []
    read_starts = []
    read_ends = []
    read_spans = []
    left_softs = []
    left_hards = []
    right_softs = []
    right_hards = []
    mqs = []
    
    # Single pass through BAM file with early filtering
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for aln in bam.fetch(until_eof=True):
            # Early filtering: skip unmapped and secondary alignments
            if aln.is_unmapped or aln.is_secondary:
                continue
                
            # Early filtering: only process reads with supplementary alignments or SA tag
            if not (aln.is_supplementary or aln.has_tag("SA")):
                continue

            # Optimized CIGAR parsing
            ct = aln.cigartuples or []
            
            # Compute left-end clipping efficiently
            left_soft = left_hard = 0
            for op, l in ct:
                if op == 4:  # soft clip
                    left_soft += l
                elif op == 5:  # hard clip
                    left_hard += l
                else:
                    break

            # Compute right-end clipping efficiently
            right_soft = right_hard = 0
            for op, l in reversed(ct):
                if op == 4:
                    right_soft += l
                elif op == 5:
                    right_hard += l
                else:
                    break

            # Get reference coordinates
            ref_start = aln.reference_start
            ref_end = aln.reference_end
            ref_span = ref_end - ref_start

            # Get read coordinates
            read_start = aln.query_alignment_start + left_hard
            read_end = aln.query_alignment_end + left_hard
            read_span = aln.query_alignment_end - aln.query_alignment_start

            # Append data efficiently
            qnames.append(aln.query_name)
            types.append("SUPPLEMENTARY" if aln.is_supplementary else "PRIMARY")
            rnames.append(bam.get_reference_name(aln.reference_id))
            strands.append("-" if aln.is_reverse else "+")
            ref_starts.append(ref_start)
            ref_ends.append(ref_end)
            ref_spans.append(ref_span)
            read_starts.append(read_start)
            read_ends.append(read_end)
            read_spans.append(read_span)
            left_softs.append(left_soft)
            left_hards.append(left_hard)
            right_softs.append(right_soft)
            right_hards.append(right_hard)
            mqs.append(aln.mapping_quality)

    # Create DataFrame efficiently with pre-allocated lists
    return pd.DataFrame({
        "QNAME": qnames,
        "TYPE": types,
        "RNAME": rnames,
        "STRAND": strands,
        "REF_START": ref_starts,
        "REF_END": ref_ends,
        "REF_SPAN": ref_spans,
        "READ_START": read_starts,
        "READ_END": read_ends,
        "READ_SPAN": read_spans,
        "LEFT_SOFT": left_softs,
        "LEFT_HARD": left_hards,
        "RIGHT_SOFT": right_softs,
        "RIGHT_HARD": right_hards,
        "MQ": mqs,
    })


def build_links_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Optimized version of links DataFrame construction with improved performance.
    
    From a DataFrame of split‐read pieces (with columns QNAME, RNAME, STRAND,
    REF_START, REF_END, READ_START, READ_END, piece_order, etc.), build a
    links table of consecutive piece‐pairs.  Classify each link as per your
    logic (deletion / intra_translocation / other), but *sort* each pair so
    that:
      - RNAME.1 is the lexicographically ("naturally") smallest chromosome, or
        if on the same chr, the piece whose 5′ mapping is more 5′.
      - coord_1 corresponds to RNAME.1 (the "tail" for that piece),
        coord_2 to RNAME.2 (the "head" for that piece).
    genomic_gap is set to –1 if the two pieces land on different chromosomes,
    else head2–tail1.
    
    Performance improvements:
    - Pre-allocated lists for better memory efficiency
    - Optimized event classification logic
    - Reduced function calls in loops
    - Early exits for edge cases
    """
    
    # Early exit if DataFrame is empty
    if df.empty:
        return pd.DataFrame()
    
    # Pre-allocate lists with estimated capacity
    estimated_capacity = len(df) // 2  # Conservative estimate
    qnames = []
    rname1s = []
    rname2s = []
    coord1s = []
    coord2s = []
    genomic_gaps = []
    events = []
    
    # Cache natural_key function to avoid repeated compilation
    def natural_key(s: str):
        # split digits vs letters, so e.g. "chr2" < "chr10"
        return [
            int(tok) if tok.isdigit() else tok.lower() for tok in re.split(r"(\d+)", s)
        ]

    # Process groups efficiently
    for qname, grp in df.groupby("QNAME", sort=False, observed=True):
        # Sort by piece_order once
        seq = grp.sort_values("piece_order")
        
        # Early exit if only one piece
        if len(seq) < 2:
            continue
            
        # Process consecutive pairs efficiently
        for i in range(len(seq) - 1):
            p1 = seq.iloc[i]
            p2 = seq.iloc[i + 1]
            
            # compute tail/head
            tail1 = p1.REF_END
            head2 = p2.REF_START

            # Optimized event classification logic
            if p1.RNAME == p2.RNAME:
                # Same chromosome events
                same_strand = p1.STRAND == p2.STRAND
                if same_strand:
                    if p1.STRAND == "+":
                        ev = "deletion" if p1.REF_START < p2.REF_START else "intra_translocation"
                    else:  # p1.STRAND == "-"
                        ev = "deletion" if p1.REF_START < p2.REF_START else "intra_translocation"
                else:
                    ev = "inversion"
            else:
                # Different chromosomes
                ev = "translocation"

            # Optimized ordering logic
            five1 = p1.REF_START if p1.STRAND == "+" else p1.REF_END
            five2 = p2.REF_START if p2.STRAND == "+" else p2.REF_END

            # Cache natural key comparisons
            nk1 = natural_key(p1.RNAME)
            nk2 = natural_key(p2.RNAME)
            
            first_is_p1 = (nk1 < nk2) or (nk1 == nk2 and five1 <= five2)

            if first_is_p1:
                r1, r2 = p1.RNAME, p2.RNAME
                c1, c2 = tail1, head2
            else:
                r1, r2 = p2.RNAME, p1.RNAME
                c1, c2 = head2, tail1

            # genomic_gap per spec
            genomic_gap = -1 if r1 != r2 else (c2 - c1)

            # Append data efficiently
            qnames.append(qname)
            rname1s.append(r1)
            rname2s.append(r2)
            coord1s.append(c1)
            coord2s.append(c2)
            genomic_gaps.append(genomic_gap)
            events.append(ev)

    # Create DataFrame efficiently with pre-allocated lists
    return pd.DataFrame({
        "QNAME": qnames,
        "RNAME.1": rname1s,
        "RNAME.2": rname2s,
        "coord_1": coord1s,
        "coord_2": coord2s,
        "genomic_gap": genomic_gaps,
        "event": events,
    })


def annotate_df(df):
    """
    Optimized version of DataFrame annotation with improved performance.
    
    Performance improvements:
    - Early filtering to reduce memory usage
    - Vectorized operations where possible
    - Efficient dtype optimization
    - Reduced DataFrame copies
    """
    # Early exit if DataFrame is empty
    if df.empty:
        return df
    
    # Apply filters first to reduce memory usage (vectorized operations)
    mq_mask = df["MQ"].values >= 55
    chr_mask = df["RNAME"].values != "chrM"
    combined_mask = mq_mask & chr_mask
    
    # Apply filter inplace to avoid copy
    df = df[combined_mask].reset_index(drop=True)
    
    # Early exit if no data after filtering
    if df.empty:
        return df

    # Optimize dtypes using astype() with dictionary for memory efficiency
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

    # Apply dtypes efficiently
    df = df.astype(dtype_optimizations, errors="ignore")

    # Sort efficiently with optimized parameters
    df = df.sort_values(["QNAME", "TYPE", "REF_START"], ignore_index=True)

    # Filter for primary alignments efficiently
    primary_mask = df["TYPE"] == "PRIMARY"
    primary_qnames = df.loc[primary_mask, "QNAME"].unique()
    
    # Filter for reads that exist in primary set
    df = df[df["QNAME"].isin(primary_qnames)]
    
    # Filter for duplicated reads (keep=False means mark all duplicates)
    duplicated_mask = df["QNAME"].duplicated(keep=False)
    df = df[duplicated_mask].reset_index(drop=True)
    
    # Early exit if no duplicated reads
    if df.empty:
        return df
    
    # Add piece_order efficiently using groupby with observed=True
    df["piece_order"] = df.groupby("QNAME", observed=True)["READ_START"].transform(
        lambda x: x.rank(method="first").astype(int)
    )
    
    return df


def get_summary(links_df, min_support=2):
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

    # then e.g. pick deletions with enough support
    return summary[
        # (summary['predominant_event'] == 'deletion') &
        (summary["support_count"] >= min_support)
    ].copy()


def has_supplementary(bam_file_path):
    """
    Quickly checks if a BAM file has any supplementary alignments.

    Parameters:
        bam_file_path (str): Path to the BAM file.

    Returns:
        bool: True if supplementary alignment is found, else False.
    """
    import time
    start = time.time()
    with pysam.AlignmentFile(bam_file_path, "rb", check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_supplementary:
                check_time = time.time() - start
                print(f"      has_supplementary (FOUND): {check_time:.3f}s")
                return True
    check_time = time.time() - start
    print(f"      has_supplementary (NOT FOUND): {check_time:.3f}s")
    return False

def merge_overlapping_intervals(group: pd.DataFrame) -> pd.DataFrame:
    """
    Merge overlapping intervals within a single group (chrom, strand, sv_type).
    For any rows that overlap or touch, collapse them into one row with:
       min_bin = the smallest min_bin
       max_bin = the largest max_bin
    Returns a DataFrame of merged rows for this group.
    """
    # Sort by min_bin so we can process in ascending order
    group = group.sort_values(by="min_bin")

    merged = []
    current_min = None
    current_max = None

    for row in group.itertuples():
        # For clarity, row has attributes like row.min_bin, row.max_bin
        if current_min is None:
            # first interval in this group
            current_min = row.min_bin
            current_max = row.max_bin
        else:
            # check if this interval overlaps the current merge
            if row.min_bin <= current_max:
                # Overlap => expand the current_max if needed
                current_max = max(current_max, row.max_bin)
            else:
                # No overlap => finalize the previous merged interval
                merged.append(
                    {
                        "chrom": row.chrom,
                        "strand": row.strand,
                        "sv_type": row.sv_type,
                        "min_bin": current_min,
                        "max_bin": current_max,
                    }
                )
                # start a new merge interval
                current_min = row.min_bin
                current_max = row.max_bin

    # Finalize the last interval if it exists
    if current_min is not None:
        merged.append(
            {
                "chrom": group["chrom"].iloc[0],  # same for entire group
                "strand": group["strand"].iloc[0],  # same for entire group
                "sv_type": group["sv_type"].iloc[0],
                "min_bin": current_min,
                "max_bin": current_max,
            }
        )

    return pd.DataFrame(merged)


def get_gene_network(gene_pairs):
    G = nx.Graph()
    for pair in gene_pairs:
        G.add_edge(pair[0], pair[1])
    connected_components = list(nx.connected_components(G))
    return [list(component) for component in connected_components]


def collapse_overlaps(df: pd.DataFrame) -> pd.DataFrame:
    """
    For each (chrom, strand, sv_type), merge overlapping [min_bin, max_bin] intervals.
    Returns a new DataFrame with the collapsed intervals.
    """
    grouped = df.groupby(["chrom", "strand", "sv_type"], group_keys=False)
    merged_df = grouped.apply(merge_overlapping_intervals)
    return merged_df.reset_index(drop=True)


def build_breakpoint_graph(df, max_proximity=500000000, group_by_sv=False):
    """
    Ultra-optimized version with minimal memory overhead and maximum performance.

    This function identifies structural variant breakpoints by:
    1. Connecting reads that belong to the same read (QNAME)
    2. Connecting reads that are within proximity on the same chromosome/strand
    3. Finding connected components (clusters) of related breakpoints
    4. Converting clusters to BED format lines for visualization

    Args:
        df: DataFrame containing structural variant reads with columns:
            - RNAME: chromosome name
            - REF_START/REF_END: reference coordinates
            - STRAND: strand (+ or -)
            - SV_TYPE: structural variant type
            - QNAME: read name
        max_proximity: Maximum distance (bp) to consider reads as proximal
        group_by_sv: Whether to group by SV type in addition to chrom/strand

    Returns:
        List of BED format strings representing breakpoint regions
    """
    # Early exit if no data to process
    if df.empty:
        return []

    # Step 1: Data Preparation - Convert DataFrame columns to numpy arrays for performance
    # Extract chromosome names and ensure they're strings
    chroms = df["RNAME"].values.astype(str)  # Ensure string dtype
    # Extract start and end positions as 32-bit integers (memory efficient)
    starts = df["REF_START"].values.astype(np.int32)
    ends = df["REF_END"].values.astype(np.int32)
    # Calculate center positions of each alignment (midpoint)
    positions = ((starts + ends) // 2).astype(np.int32)
    # Extract strand information (+ or -)
    strands = df["STRAND"].values.astype(str)  # Ensure string dtype
    # Extract structural variant types (DELETION, INSERTION, etc.)
    sv_types = df["SV_TYPE"].values.astype(str)  # Ensure string dtype
    # Extract read names (QNAME)
    qnames = df["QNAME"].values.astype(str)  # Ensure string dtype

    # Step 2: Edge Set Initialization
    # Pre-allocate edges set with conservative size estimate to avoid resizing
    edges = set()  # Will store pairs of indices representing connected reads

    # Step 3: Same Read Connections - Connect all alignments from the same read
    # Find unique read names and their inverse mapping for efficient grouping
    unique_qnames, qname_inverse = np.unique(qnames, return_inverse=True)
    # Iterate through each unique read name
    for qname_idx in range(len(unique_qnames)):
        # Find all row indices that belong to this read
        read_indices = np.where(qname_inverse == qname_idx)[0]
        # If read has multiple alignments, connect them all pairwise
        if len(read_indices) > 1:
            # Use itertools.combinations for faster pairwise generation
            # This creates all possible pairs of indices for this read
            for i, j in combinations(read_indices, 2):
                edges.add((i, j))  # Add edge between these two alignments

    # Step 4: Proximity Connections - Connect reads that are close to each other
    # Create grouping keys based on chromosome, strand, and optionally SV type
    if group_by_sv:
        # Group by chromosome, strand, AND structural variant type
        group_keys = [
            f"{chrom}_{strand}_{sv_type}"
            for chrom, strand, sv_type in zip(chroms, strands, sv_types)
        ]
    else:
        # Group only by chromosome and strand
        group_keys = [f"{chrom}_{strand}" for chrom, strand in zip(chroms, strands)]

    # Create temporary DataFrame for efficient grouping operations
    df_temp = pd.DataFrame(
        {
            "group_key": group_keys,  # Grouping key for each row
            "position": positions,  # Center position of each alignment
            "index": np.arange(len(df)),  # Original row index
        }
    )

    # Process each group (same chromosome, strand, SV type) separately
    for group_key, group in df_temp.groupby("group_key"):
        # Skip groups with only one alignment (no proximity connections possible)
        if len(group) < 2:
            continue

        # Sort by position for efficient proximity search (O(n log n) vs O(n²))
        group_sorted = group.sort_values("position")
        group_indices = group_sorted["index"].values  # Original row indices
        group_positions = group_sorted["position"].values  # Sorted positions

        # For each position, find all subsequent positions within max_proximity
        for i in range(len(group_positions)):
            # Vectorized proximity search: find all positions within max_proximity
            # This is much faster than nested loops
            proximity_mask = (
                group_positions[i + 1 :] - group_positions[i]
            ) <= max_proximity
            if np.any(proximity_mask):
                # Add edges to all positions within proximity
                # Calculate how many positions are within proximity
                num_proximal = np.sum(proximity_mask)
                # Add edges to all proximal positions
                for j in range(i + 1, i + 1 + num_proximal):
                    edges.add((group_indices[i], group_indices[j]))

    # Step 5: Connected Components Detection using Union-Find Algorithm
    # Define Union-Find data structure for efficient connected component detection
    class UnionFind:
        def __init__(self, n):
            # Initialize parent array: each element is its own parent initially
            self.parent = list(range(n))
            # Initialize rank array for union-by-rank optimization
            self.rank = [0] * n

        def find(self, x):
            # Path compression: make all nodes on path point directly to root
            if self.parent[x] != x:
                self.parent[x] = self.find(self.parent[x])
            return self.parent[x]

        def union(self, x, y):
            # Union by rank: attach smaller tree to root of larger tree
            px, py = self.find(x), self.find(y)
            if px == py:  # Already in same component
                return
            if self.rank[px] < self.rank[py]:
                self.parent[px] = py
            elif self.rank[px] > self.rank[py]:
                self.parent[py] = px
            else:
                self.parent[py] = px
                self.rank[px] += 1  # Increase rank when trees have same height

    # Build connected components using Union-Find
    uf = UnionFind(len(df))  # Initialize Union-Find with number of rows
    # Union all connected pairs
    for edge in edges:
        uf.union(edge[0], edge[1])

    # Step 6: Group nodes by their connected component
    components = {}  # Dictionary: root -> list of indices in this component
    for i in range(len(df)):
        root = uf.find(i)  # Find the root of this node
        if root not in components:
            components[root] = []
        components[root].append(i)  # Add this node to its component

    # Step 7: Process components and generate BED format lines
    bed_lines = []  # Will store BED format strings

    # Process each connected component
    for component_indices in components.values():
        # Filter components: must have 3-100 alignments (avoid noise and huge clusters)
        if len(component_indices) < 3 or len(component_indices) > 100:
            continue

        # Extract data for this component efficiently using numpy indexing
        comp_chroms = chroms[component_indices]  # Chromosomes in this component
        comp_positions = positions[component_indices]  # Positions in this component
        comp_strands = strands[component_indices]  # Strands in this component
        comp_sv_types = sv_types[component_indices]  # SV types in this component
        comp_qnames = qnames[component_indices]  # Read names in this component

        # Check component quality criteria
        unique_chroms = set(comp_chroms)  # Number of different chromosomes
        unique_reads = set(comp_qnames)  # Number of different reads

        # Quality filter: must have >2 unique reads and <4 chromosomes
        # This ensures we have multiple supporting reads but not too many chromosomes
        if len(unique_reads) > 2 and len(unique_chroms) < 4:
            # Calculate BED coordinates using vectorized operations
            bin_size = 1000  # 1kb bins for coordinate calculation
            bins = comp_positions // bin_size  # Convert positions to bin numbers

            # Create strand-specific masks for different bin calculations
            strand_mask = comp_strands == "+"  # Boolean mask for positive strand

            # Calculate min/max bins differently for + and - strands
            # Positive strand: extend 1 bin upstream, 5 bins downstream
            # Negative strand: extend 5 bins upstream, 1 bin downstream
            min_bins = np.where(strand_mask, bins - 1, bins - 5) * bin_size
            max_bins = np.where(strand_mask, bins + 5, bins + 1) * bin_size

            # Create BED lines directly without intermediate DataFrame
            for i in range(len(component_indices)):
                chrom = comp_chroms[i]  # Chromosome name
                strand = comp_strands[i]  # Strand (+ or -)
                # Use SV type if grouping by SV, otherwise "UNKNOWN"
                sv_type = comp_sv_types[i] if group_by_sv else "UNKNOWN"
                min_bin = min_bins[i]  # Start coordinate
                max_bin = max_bins[i]  # End coordinate

                # Format as BED line: chrom start end name score strand
                bed_line = f"{chrom}\t{min_bin}\t{max_bin}\t{sv_type}\t.\t{strand}"
                bed_lines.append(bed_line)

    # Step 8: Post-processing
    # Remove duplicates efficiently while preserving order
    bed_lines = list(dict.fromkeys(bed_lines))  # Preserves order, removes duplicates

    # Limit output size for performance (prevent memory issues with huge datasets)
    if len(bed_lines) > 1000:
        bed_lines = bed_lines[:1000]
        logger.warning("Limited bed lines output to 1000 entries for performance")

    return bed_lines


def _get_reads(reads: pd.DataFrame) -> pd.DataFrame:
    """
    Get reads for a specific gene.

    Args:
        reads (pd.DataFrame): DataFrame with reads

    Returns:
        pd.DataFrame: DataFrame with reads
    """
    df = reads
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
    result = (
        df.groupby("chromosome")
        .apply(lambda x: collapse_ranges(x, 10000))
        .reset_index(drop=True)
    )
    return result


def _annotate_results(result: pd.DataFrame) -> Tuple[pd.DataFrame, pd.Series]:
    """
    Annotates the result DataFrame with tags and colors.

    Performance Warning:
    - Creates unnecessary copy of entire DataFrame (result_copy = result.copy())
    - Consider using inplace operations or views where possible
    - Multiple groupby operations could be combined
    - String operations could be vectorized

    Optimization Suggestions:
    1. Use inplace operations where possible
    2. Combine groupby operations
    3. Use categorical dtypes
    4. Vectorize string operations
    5. Consider using numpy for numerical operations
    """
    # TODO: Optimize by removing unnecessary copy and using inplace operations
    result_copy = result.copy()  # Unnecessary copy of entire DataFrame
    # Group by read_id and aggregate col4 (Gene) values
    lookup = result_copy.groupby("read_id", observed=True)["col4"].agg(
        lambda x: ",".join(set(x))
    )
    tags = result_copy["read_id"].map(lookup.get)
    result_copy.loc[:, "tag"] = tags
    result = result_copy
    # Generate colors for each read_id group
    colors = result.groupby("read_id", observed=True).apply(
        lambda x: _generate_random_color()
    )
    result = result.map(lambda x: x.strip() if isinstance(x, str) else x)
    result["Color"] = result["read_id"].map(colors.get)
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


def safe_read_csv(file_path, dtype=None):
    """
    Safely read a CSV file that might be compressed.
    
    Args:
        file_path: Path to the CSV file
        dtype: Optional dtype specification for pandas
        
    Returns:
        pd.DataFrame: The loaded DataFrame
    """
    try:
        # First try reading as uncompressed CSV
        return pd.read_csv(file_path, dtype=dtype)
    except UnicodeDecodeError:
        # If that fails, try reading as gzip compressed
        try:
            return pd.read_csv(file_path, compression='gzip', dtype=dtype)
        except Exception as e:
            logger.error(f"Failed to read CSV file {file_path} with both uncompressed and gzip methods: {str(e)}")
            raise


def preprocess_fusion_data_standalone(fusion_data: pd.DataFrame, output_file: str) -> None:
    """
    Standalone version of fusion data preprocessing for CPU-bound execution.
    
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
            "candidate_count": 0
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
            
            processed_data.update({
                "gene_pairs": gene_pairs,
                "gene_groups": gene_groups,
                "candidate_count": len(gene_groups)
            })
        
        # Save processed data as pickle for efficient loading
        import pickle
        with open(output_file, 'wb') as f:
            pickle.dump(processed_data, f)
            
        logger.info(f"Pre-processed fusion data saved to {output_file}")
        
    except Exception as e:
        logger.error(f"Error pre-processing fusion data: {str(e)}")
        logger.error("Exception details:", exc_info=True)


def preprocess_structural_variants_standalone(output_dir: str) -> None:
    """
    Standalone version of structural variant preprocessing for CPU-bound execution.
    
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
                        os.path.join(output_dir, "structural_variants_processed.csv"),
                        index=False
                    )
                    
                    # Save count for quick access
                    with open(os.path.join(output_dir, "sv_count.txt"), 'w') as f:
                        f.write(str(len(sv_df)))
                    
                    logger.info(
                        f"Pre-processed {len(sv_df)} structural variant events"
                    )
                else:
                    # Save empty count
                    with open(os.path.join(output_dir, "sv_count.txt"), 'w') as f:
                        f.write("0")
            else:
                # Save empty count
                with open(os.path.join(output_dir, "sv_count.txt"), 'w') as f:
                    f.write("0")
                    
    except Exception as e:
        logger.error(f"Error pre-processing structural variants: {str(e)}")
        logger.error("Exception details:", exc_info=True)


def process_bam_pipeline(bamfile):
    """
    Complete pipeline for processing BAM files to extract structural variant links.
    
    This function combines all steps in a single CPU-bound task to reduce async overhead:
    1. Extract split read alignments from BAM file
    2. Annotate the extracted data
    3. Build links from annotated data
    
    Args:
        bamfile: Path to the BAM file to process
        
    Returns:
        pd.DataFrame: DataFrame containing structural variant links
    """
    import time
    start_total = time.time()
    print(f"    process_bam_pipeline START: {bamfile}")
    
    # Step 1: Extract split read alignments from BAM file
    start_extract = time.time()
    split_reads_df = extract_split_read_alignments(bamfile)
    extract_time = time.time() - start_extract
    print(f"      Extract split reads: {extract_time:.3f}s")
    
    # Early exit if no split reads found
    if split_reads_df.empty:
        logger.debug("No split reads found in BAM file")
        total_time = time.time() - start_total
        print(f"    process_bam_pipeline END (no split reads): {total_time:.3f}s")
        return pd.DataFrame()
    
    # Step 2: Annotate the extracted data
    start_annotate = time.time()
    annotated_df = annotate_df(split_reads_df)
    annotate_time = time.time() - start_annotate
    print(f"      Annotate data: {annotate_time:.3f}s")
    
    # Early exit if no annotated data
    if annotated_df.empty:
        logger.debug("No annotated data after filtering")
        total_time = time.time() - start_total
        print(f"    process_bam_pipeline END (no annotated data): {total_time:.3f}s")
        return pd.DataFrame()
    
    # Step 3: Build links from annotated data
    start_build_links = time.time()
    new_df = build_links_df(annotated_df)
    build_links_time = time.time() - start_build_links
    print(f"      Build links: {build_links_time:.3f}s")
    
    total_time = time.time() - start_total
    print(f"    process_bam_pipeline END: {total_time:.3f}s")
    return new_df


def extract_bam_info(bam_file):
    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Initialize lists to store information
    read_ids = []
    reference_ids = []
    reference_starts = []
    reference_ends = []
    read_starts = []
    read_ends = []
    strands = []
    mapping_qualities = []
    is_secondary = []
    is_supplementary = []

    # Iterate over each read in the BAM file
    for read in bam:
        # Append the read information to the lists
        read_ids.append(read.query_name)
        reference_ids.append(bam.get_reference_name(read.reference_id))
        reference_starts.append(read.reference_start)
        reference_ends.append(read.reference_end)
        read_starts.append(read.query_alignment_start)
        read_ends.append(read.query_alignment_end)
        strands.append("-" if read.is_reverse else "+")
        mapping_qualities.append(read.mapping_quality)
        is_secondary.append(read.is_secondary)
        is_supplementary.append(read.is_supplementary)

    # Close the BAM file
    bam.close()

    # Create a DataFrame
    df = pd.DataFrame(
        {
            "read_id": read_ids,
            "reference_id": reference_ids,
            "reference_start": reference_starts,
            "reference_end": reference_ends,
            "read_start": read_starts,
            "read_end": read_ends,
            "strand": strands,
            "mapping_quality": mapping_qualities,
            "is_secondary": is_secondary,
            "is_supplementary": is_supplementary,
        }
    )

    return df


# Function to collapse ranges within a fixed distance
def collapse_ranges(df, max_distance):
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


def run_command(command: str) -> None:
    """
    Run a shell command and handle exceptions.

    Args:
        command (str): The shell command to run.
    """
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        # logging.error(f"Command failed: {command}")
        raise e


def fusion_work_pysam(
    bamfile: str,
    gene_bed: str,
    all_gene_bed: str,
) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """
    Simplified fusion candidate detection using pure pysam instead of external tools.

    Args:
        bamfile: Path to the BAM file
        gene_bed: Path to target gene BED file
        all_gene_bed: Path to all genes BED file

    Returns:
        Tuple of (fusion_candidates, fusion_candidates_all) DataFrames
    """
    import time
    start_total = time.time()
    print(f"    fusion_work_pysam START: {bamfile}")
    
    fusion_candidates: Optional[pd.DataFrame] = None
    fusion_candidates_all: Optional[pd.DataFrame] = None

    try:
        # Load BED files into memory for efficient lookup
        start_bed_load = time.time()
        gene_regions = load_bed_regions(gene_bed)
        all_gene_regions = load_bed_regions(all_gene_bed)
        bed_load_time = time.time() - start_bed_load
        print(f"      Load BED regions: {bed_load_time:.3f}s")

        # Find reads with supplementary alignments
        start_find_supp = time.time()
        reads_with_supp = find_reads_with_supplementary(bamfile)
        find_supp_time = time.time() - start_find_supp
        print(f"      Find supplementary reads: {find_supp_time:.3f}s")

        if not reads_with_supp:
            logger.debug("No reads with supplementary alignments found")
            total_time = time.time() - start_total
            print(f"    fusion_work_pysam END (no supp): {total_time:.3f}s")
            return None, None

        # Process reads and find gene intersections
        start_process_target = time.time()
        fusion_candidates = process_reads_for_fusions(
            bamfile, reads_with_supp, gene_regions
        )
        process_target_time = time.time() - start_process_target
        print(f"      Process target fusions: {process_target_time:.3f}s")
        
        start_process_all = time.time()
        fusion_candidates_all = process_reads_for_fusions(
            bamfile, reads_with_supp, all_gene_regions
        )
        process_all_time = time.time() - start_process_all
        print(f"      Process all fusions: {process_all_time:.3f}s")

    except Exception as e:
        logger.error(f"Error in fusion_work_pysam: {str(e)}")
        logger.error("Exception details:", exc_info=True)
        raise

    total_time = time.time() - start_total
    print(f"    fusion_work_pysam END: {total_time:.3f}s")
    return fusion_candidates, fusion_candidates_all


def load_bed_regions(bed_file: str) -> Dict[str, List[Tuple[int, int, str]]]:
    """
    Load BED file regions into memory for efficient lookup.

    Args:
        bed_file: Path to BED file

    Returns:
        Dictionary mapping chromosome to list of (start, end, gene_name) tuples
    """
    regions = defaultdict(list)

    with open(bed_file, "r") as f:
        for line in f:
            if line.strip():
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    gene_name = parts[3]
                    regions[chrom].append((start, end, gene_name))

    return dict(regions)


def find_reads_with_supplementary(bamfile: str) -> Set[str]:
    """
    Find all read names that have supplementary alignments.

    Args:
        bamfile: Path to BAM file

    Returns:
        Set of read names with supplementary alignments
    """
    reads_with_supp = set()

    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                continue
            # Check for supplementary alignments or SA tag
            if read.is_supplementary or read.has_tag("SA"):
                reads_with_supp.add(read.query_name)

    logger.debug(f"Found {len(reads_with_supp)} reads with supplementary alignments")
    return reads_with_supp


def process_reads_for_fusions(
    bamfile: str,
    reads_with_supp: Set[str],
    gene_regions: Dict[str, List[Tuple[int, int, str]]],
) -> Optional[pd.DataFrame]:
    """
    Process reads to find gene intersections and create fusion candidates.

    Args:
        bamfile: Path to BAM file
        reads_with_supp: Set of read names with supplementary alignmentsxx
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
    string_columns = ["BAM_FILE", "QNAME", "TYPE", "RNAME", "STRAND"]
    for col in string_columns:
        if col in df.columns:
            df[col] = df[col].astype("category")

    # Apply filters more efficiently
    df = df[(df["mapping_quality"] > 40) & (df["mapping_span"] > 100)].reset_index(
        drop=True
    )

    return df


class FusionVis(BaseVis):
    """
    FusionVis handles the visualization and UI components for gene fusion analysis.

    Performance Warning:
    - Stores multiple copies of DataFrames in memory
    - Creates unnecessary copies during table updates
    - Multiple file reads for gene annotations

    Optimization Suggestions:
    1. Implement lazy loading for large datasets
    2. Use memory-efficient data structures
    3. Cache gene annotations
    4. Implement pagination for tables
    5. Use virtual scrolling for large datasets
    6. Stream large files instead of loading entirely
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
        Initialize the FusionVis class with necessary parameters and UI components.

        Args:
            target_panel: Name of the target panel (e.g., 'rCNS2', 'AML')
            reference_file: Path to reference genome file
            bed_file: Path to BED file with target regions
            readfish_toml: Path to readfish configuration file
            master_bed_tree: Pre-computed bed tree for efficient region queries

        Performance Notes:
        - Initializes UI components lazily to reduce startup time
        - Caches gene annotations in memory
        - Uses efficient data structures for gene lookups
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

    async def setup_ui(self) -> None:
        """
        Sets up the user interface for the Fusion Panel and Structural Variant Analysis.
        """
        #app.config.request_timeout = 10  # Increase timeout to 10 seconds

        if self.summary:
            with self.summary:
                with ui.card().classes("w-full p-4 mb-4"):
                    with ui.row().classes("w-full items-center justify-between"):
                        # Left side - Fusion Status
                        with ui.column().classes("gap-2"):
                            ui.label("Gene Fusion Analysis").classes(
                                "text-lg font-medium"
                            )
                            with ui.row().classes("items-center gap-2"):
                                ui.label(f"Panel: {self.target_panel}").classes(
                                    "text-gray-600 font-medium"
                                )
                                ui.label("0").bind_text_from(
                                    self,
                                    "candidates",
                                    backward=lambda n: f"{n} between target fusions",
                                ).classes("px-2 py-1 rounded bg-blue-100 text-blue-600")

                        # Right side - Additional metrics
                        with ui.column().classes("gap-2 text-right"):
                            ui.label("Analysis Details").classes("font-medium")
                            ui.label("0").bind_text_from(
                                self,
                                "all_candidates",
                                backward=lambda n: f"{n} genome wide fusions",
                            ).classes("text-gray-600")

                    # Bottom row - Information
                    with ui.row().classes(
                        "w-full mt-4 text-sm text-gray-500 justify-center"
                    ):
                        ui.label(
                            "Fusion candidates identified from reads with supplementary alignments"
                        )

        with ui.card().style("width: 100%"):
            ui.label("Gene Fusion and Structural Variant Analysis").classes(
                "text-sky-600 dark:text-white"
            ).style("font-size: 150%; font-weight: 300").tailwind(
                "drop-shadow", "font-bold"
            )
            ui.label(
                "This panel identifies gene fusion candidates and structural variants from the input bam files. "
                "The panel is split into three tabs: gene fusions within the target panel, genome wide fusions, and structural variants. "
                "Events are identified on a streaming basis derived from reads with supplementary alignments. "
                "The plots are indicative of the presence of events and should be interpreted with care."
            ).style("font-size: 125%; font-weight: 300")

            with ui.tabs().classes("w-full") as tabs:
                one = ui.tab("Within Target Fusions").style(
                    "font-size: 125%; font-weight: 300"
                )
                with one:
                    self.badge_one = (
                        ui.badge("0", color="red")
                        .bind_text_from(self, "candidates", backward=lambda n: f"{n}")
                        .props("floating rounded outline")
                    )
                two = ui.tab("Genome Wide Fusions").style(
                    "font-size: 125%; font-weight: 300"
                )
                with two:
                    self.badge_two = (
                        ui.badge("0", color="red")
                        .bind_text_from(
                            self, "all_candidates", backward=lambda n: f"{n}"
                        )
                        .props("floating rounded outline")
                    )
                three = ui.tab("Structural Variants").style(
                    "font-size: 125%; font-weight: 300"
                )
                with three:
                    self.badge_three = (
                        ui.badge("0", color="red")
                        .bind_text_from(self, "sv_count", backward=lambda n: f"{n}")
                        .props("floating rounded outline")
                    )

            with ui.tab_panels(tabs, value=one).classes("w-full"):
                with ui.tab_panel(one):
                    with ui.card().style("width: 100%"):
                        ui.label("Fusion Candidates (within targets)").style(
                            "font-size: 125%; font-weight: 300"
                        ).tailwind("drop-shadow", "font-bold")
                        ui.separator()
                        with ui.expansion(
                            "Methods",
                            caption="A description of the methods used to identify gene fusions within the target panel",
                        ).classes("w-full"):
                            ui.restructured_text(
                                """
                                Gene fusion analysis within the target panel is performed through a multi-step process:

                                1. **Read Extraction**
                                Reads with supplementary alignments are identified and extracted from the BAM file.

                                2. **Target Gene Mapping**
                                Extracted reads are intersected with the target panel gene coordinates using bedtools.

                                3. **Fusion Candidate Identification**
                                Reads mapping to multiple genes within the panel are identified as fusion candidates.

                                The analysis requires:

                                * Minimum mapping quality > 40
                                * Minimum alignment length > 100bp
                                * Multiple high-quality alignments to different genes

                                For each fusion candidate, the following is provided:

                                * Gene pairs involved
                                * Precise genomic coordinates
                                * Read alignment details
                                * Mapping quality scores

                                Results are visualized in plots showing read alignments to each gene region, with individual reads color-coded for tracking.
                            """
                            ).style("font-size: 100%; font-weight: 300")
                        self.fusionplot = ui.row()
                        with self.fusionplot.classes("w-full"):
                            with ui.column().classes("gap-2"):
                                ui.label("Awaiting Fusion Plot").classes(
                                    "text-lg font-medium"
                                )
                                ui.label(
                                    "Plot will be displayed here when fusion data is available"
                                ).classes("text-gray-600")
                        self.fusiontable = ui.row().classes("w-full")
                        with self.fusiontable:
                            with ui.column().classes("gap-2"):
                                ui.label("Awaiting Fusion Data").classes(
                                    "text-lg font-medium"
                                )
                                ui.label(
                                    "Table will be displayed here when fusion data is available"
                                ).classes("text-gray-600")

                with ui.tab_panel(two):
                    with ui.card().style("width: 100%"):
                        ui.label("Fusion Candidates (genome wide)").style(
                            "font-size: 125%; font-weight: 300"
                        ).tailwind("drop-shadow", "font-bold")
                        ui.separator()
                        with ui.expansion(
                            "Methods",
                            caption="A description of the methods used to identify genome-wide fusions",
                        ).classes("w-full"):
                            ui.restructured_text(
                                """
                                Genome-wide fusion analysis extends the search beyond the target panel through the following steps:

                                1. **Initial Screening**
                                Reads with supplementary alignments are identified from the target panel regions.

                                2. **Extended Mapping**
                                Selected reads are then mapped against all annotated genes in the genome.

                                3. **Partner Gene Identification**
                                Additional gene partners are identified from genome-wide alignments.

                                The analysis maintains strict criteria:

                                * Minimum mapping quality > 40
                                * Minimum alignment length > 100bp
                                * At least one gene from the target panel
                                * Multiple high-quality alignments

                                Each fusion event includes:

                                * Target panel gene
                                * Partner gene(s) from anywhere in the genome
                                * Precise genomic coordinates
                                * Read alignment details
                                * Mapping quality scores

                                Results are displayed with gene-specific read alignment plots and comprehensive mapping details.
                            """
                            ).style("font-size: 100%; font-weight: 300")
                        self.fusionplot_all = ui.row()
                        with self.fusionplot_all.classes("w-full"):
                            with ui.column().classes("gap-2"):
                                ui.label("Awaiting Genome-wide Fusion Plot").classes(
                                    "text-lg font-medium"
                                )
                                ui.label(
                                    "Plot will be displayed here when genome-wide fusion data is available"
                                ).classes("text-gray-600")
                        self.fusiontable_all = ui.row().classes("w-full")
                        with self.fusiontable_all:
                            with ui.column().classes("gap-2"):
                                ui.label("Awaiting Genome-wide Fusion Data").classes(
                                    "text-lg font-weight: 300"
                                )
                                ui.label(
                                    "Table will be displayed here when genome-wide fusion data is available"
                                ).classes("text-gray-600")

                with ui.tab_panel(three):
                    with ui.card().style("width: 100%"):
                        ui.label("Structural Variants").style(
                            "font-size: 125%; font-weight: 300"
                        ).tailwind("drop-shadow", "font-bold")
                        ui.separator()
                        with ui.expansion(
                            "Methods",
                            caption="A description of the methods used to identify structural variants",
                        ).classes("w-full"):
                            ui.restructured_text(
                                """
                                The analysis identifies structural variants through a multi-step process:

                                1. **Breakpoint Detection**
                                Reads with supplementary alignments are extracted and analyzed to identify potential breakpoints.

                                2. **Event Clustering**
                                Breakpoints are clustered within 50kb to identify event boundaries.

                                3. **Partner Detection**
                                Read pairs mapping to different chromosomes are used to identify translocation partners.

                                The analysis identifies various structural variants including:

                                * Deletions
                                * Insertions  
                                * Inversions
                                * Translocations

                                For each event, the following information is provided:

                                * Precise genomic coordinates
                                * Event type and size
                                * For translocations: both primary and partner chromosomal locations

                                The events are summarized in plots below and detailed in the table.
                            """
                            ).style("font-size: 100%; font-weight: 300")

                        # self.sv_plot = ui.row().classes("w-full")
                        self.sv_table_container = ui.row().classes("w-full")
        #await ui.context.client.connected()
        if self.browse:
            self.show_previous_data()
        else:
            # Timer now calls lightweight show_previous_data method
            ui.timer(30, lambda: self.show_previous_data())

    def update_fusion_table_all(self, result_all: pd.DataFrame) -> None:
        """
        Updates the UI table with all fusion candidates (genome-wide).

        Args:
            result_all: DataFrame containing all fusion candidates

        Performance Optimizations:
        - Uses categorical data types for string columns
        - Pre-sorts data for efficient updates
        - Implements efficient filtering
        - Uses pagination to handle large datasets

        Potential Improvements:
        - Implement virtual scrolling for large tables
        - Add data compression for large datasets
        - Cache frequently accessed data
        - Use parallel processing for data transformations
        """
        # Add debug logging to see what columns we actually have
        logger.debug(
            f"DataFrame columns in update_fusion_table_all: {result_all.columns.tolist()}"
        )
        logger.debug(f"DataFrame head in update_fusion_table_all:\n{result_all.head()}")

        # Pre-sort and use categorical data types
        result_all = result_all.astype(
            {
                "read_id": "category",
                "col4": "category",  # This is the Gene column
                "reference_id": "category",  # This is the chrom column
                "strand": "category",
            }
        )

        if result_all.shape[0] > self.fstable_all_row_count:
            self.fstable_all_row_count = result_all.shape[0]
            if not self.fstable_all:
                self.fusiontable_all.clear()
                with self.fusiontable_all:
                    self.fstable_all = (
                        ui.table.from_pandas(
                            result_all.sort_values(by="reference_start").rename(
                                columns={
                                    "col1": "chromBED",
                                    "col2": "BS",
                                    "col3": "BE",
                                    "col4": "Gene",
                                    "reference_id": "chrom",
                                    "reference_start": "mS",
                                    "reference_end": "mE",
                                    "read_id": "readID",
                                    "mapping_quality": "mapQ",
                                    "strand": "strand",
                                    "read_start": "Read Map Start",
                                    "read_end": "Read Map End",
                                    "is_secondary": "Secondary",
                                    "is_supplementary": "Supplementary",
                                    "mapping_span": "mapping span",
                                }
                            ),
                            pagination=25,
                        )
                        .props("dense")
                        .classes("w-full")
                        .style("height: 900px")
                        .style("font-size: 100%; font-weight: 300")
                    )
                    for col in self.fstable_all.columns:
                        col["sortable"] = True

                    with self.fstable_all.add_slot("top-right"):
                        with ui.input(placeholder="Search").props(
                            "type=search"
                        ).bind_value(self.fstable_all, "filter").add_slot("append"):
                            ui.icon("search")

                self.fusionplot_all.clear()
            else:
                self.fstable_all.update_rows(
                    result_all.sort_values(by="reference_start")
                    .rename(
                        columns={
                            "col1": "chromBED",
                            "col2": "BS",
                            "col3": "BE",
                            "col4": "Gene",
                            "reference_id": "chrom",
                            "reference_start": "mS",
                            "reference_end": "mE",
                            "read_id": "readID",
                            "mapping_quality": "mapQ",
                            "strand": "strand",
                            "read_start": "Read Map Start",
                            "read_end": "Read Map End",
                            "is_secondary": "Secondary",
                            "is_supplementary": "Supplementary",
                            "mapping_span": "mapping span",
                        }
                    )
                    .to_dict("records")
                )

                self.fstable_all.update()
                self.fusionplot_all.clear()

            result_all, goodpairs = _annotate_results(result_all)
            self.all_candidates = 0

            if not result_all.empty:
                with self.fusionplot_all.classes("w-full"):
                    gene_pairs = (
                        result_all[goodpairs]
                        .sort_values(by="reference_start")["tag"]
                        .unique()
                        .tolist()
                    )
                    gene_pairs = [pair.split(", ") for pair in gene_pairs]
                    gene_groups_test = get_gene_network(gene_pairs)
                    gene_groups = []
                    for gene_group in gene_groups_test:
                        reads = _get_reads(
                            result_all[goodpairs][
                                result_all[goodpairs]["col4"].isin(gene_group)
                            ]
                        )
                        if len(reads) > 1:
                            gene_groups.append(gene_group)
                    self.all_candidates = len(gene_groups)
                    with ui.row().classes("w-full"):
                        ui.select(
                            options=gene_groups,
                            with_input=True,
                            on_change=lambda e: show_gene_pair(
                                e.value, result_all, goodpairs
                            ),
                        ).classes("w-40")
                    with ui.row().classes("w-full"):
                        self.all_card = ui.card()
                        with self.all_card:
                            ui.label("Select gene pair to see results.").tailwind(
                                "drop-shadow", "font-bold"
                            )

                    def show_gene_pair(gene_group: str, result_all, goodpairs) -> None:
                        ui.notify(gene_group)
                        self.all_card.clear()
                        with self.all_card:
                            with ui.row():
                                ui.label(f"{gene_group}").tailwind(
                                    "drop-shadow", "font-bold"
                                )
                            with ui.row():
                                reads = result_all[goodpairs][
                                    result_all[goodpairs]["col4"].isin(gene_group)
                                ]
                                self.create_fusion_plot(gene_group, reads)

    def update_fusion_table(self, result: pd.DataFrame) -> None:
        """
        Updates the UI table with fusion candidates within target regions.

        Args:
            result: DataFrame containing fusion candidates

        Performance Optimizations:
        - Uses categorical data types
        - Implements efficient indexing
        - Pre-sorts data
        - Uses sets for unique operations

        Potential Improvements:
        - Implement incremental updates
        - Add data caching
        - Use parallel processing for data transformations
        - Implement virtual scrolling
        """
        # Add debug logging to see what columns we actually have
        logger.debug(
            f"DataFrame columns in update_fusion_table: {result.columns.tolist()}"
        )
        logger.debug(f"DataFrame head in update_fusion_table:\n{result.head()}")

        # Pre-sort and use categorical data types
        result = result.astype(
            {
                "read_id": "category",
                "col4": "category",  # This is the Gene column
                "reference_id": "category",  # This is the chrom column
                "strand": "category",
            }
        )

        # Annotate results and get goodpairs before using them
        result, goodpairs = _annotate_results(result)

        # Use more efficient operations
        if not hasattr(self, "_sorted_result"):
            self._sorted_result = result.sort_values(by="reference_start")

        # Use sets for unique operations - convert lists to tuples to make them hashable
        gene_pairs = set(
            tuple(pair.split(",")) for pair in result[goodpairs]["tag"].unique()
        )

        # Implement better indexing
        if not hasattr(self, "_gene_index"):
            self._gene_index = result.set_index(
                "col4"
            )  # Using col4 which is the Gene column

        if result.shape[0] > self.fstable_row_count:
            self.fstable_row_count = result.shape[0]
            if not self.fstable:
                self.fusiontable.clear()
                with self.fusiontable:
                    self.fstable = (
                        ui.table.from_pandas(
                            result.sort_values(by="reference_start").rename(
                                columns={
                                    "col1": "chromBED",
                                    "col2": "BS",
                                    "col3": "BE",
                                    "col4": "Gene",
                                    "reference_id": "chrom",
                                    "reference_start": "mS",
                                    "reference_end": "mE",
                                    "read_id": "readID",
                                    "mapping_quality": "mapQ",
                                    "strand": "strand",
                                    "read_start": "Read Map Start",
                                    "read_end": "Read Map End",
                                    "is_secondary": "Secondary",
                                    "is_supplementary": "Supplementary",
                                    "mapping_span": "mapping span",
                                }
                            ),
                            pagination=25,
                        )
                        .props("dense")
                        .classes("w-full")
                        .style("height: 900px")
                        .style("font-size: 100%; font-weight: 300")
                    )
                    for col in self.fstable.columns:
                        col["sortable"] = True

                    with self.fstable.add_slot("top-right"):
                        with ui.input(placeholder="Search").props(
                            "type=search"
                        ).bind_value(self.fstable, "filter").add_slot("append"):
                            ui.icon("search")

                self.fusionplot.clear()
            else:
                self.fstable.update_rows(
                    result.sort_values(by="reference_start")
                    .rename(
                        columns={
                            "col1": "chromBED",
                            "col2": "BS",
                            "col3": "BE",
                            "col4": "Gene",
                            "reference_id": "chrom",
                            "reference_start": "mS",
                            "reference_end": "mE",
                            "read_id": "readID",
                            "mapping_quality": "mapQ",
                            "strand": "strand",
                            "read_start": "Read Map Start",
                            "read_end": "Read Map End",
                            "is_secondary": "Secondary",
                            "is_supplementary": "Supplementary",
                            "mapping_span": "mapping span",
                        }
                    )
                    .to_dict("records")
                )

                self.fusionplot.clear()

            self.candidates = 0

            if not (result.empty):
                with self.fusionplot.classes("w-full"):
                    gene_pairs = (
                        result[goodpairs]
                        .sort_values(by="reference_start")["tag"]
                        .unique()
                        .tolist()
                    )
                    # Convert to tuples for network analysis
                    gene_pairs = [tuple(pair.split(",")) for pair in gene_pairs]
                    gene_groups_test = get_gene_network(gene_pairs)
                    gene_groups = []
                    for gene_group in gene_groups_test:
                        reads = _get_reads(
                            result[goodpairs][
                                result[goodpairs]["col4"].isin(gene_group)
                            ]
                        )
                        if len(reads) > 1:
                            gene_groups.append(gene_group)
                    self.candidates = len(gene_groups)
                    with ui.row().classes("w-full"):
                        ui.select(
                            options=gene_groups,
                            with_input=True,
                            on_change=lambda e: show_gene_pair(
                                e.value, result, goodpairs
                            ),
                        ).classes("w-40")
                    with ui.row().classes("w-full"):
                        self.card = ui.card()
                        with self.card:
                            ui.label("Select gene pair to see results.").tailwind(
                                "drop-shadow", "font-bold"
                            )

                    def show_gene_pair(gene_group: str, result, goodpairs) -> None:
                        ui.notify(gene_group)
                        self.card.clear()
                        with self.card:
                            with ui.row():
                                ui.label(f"{gene_group}").tailwind(
                                    "drop-shadow", "font-bold"
                                )
                            with ui.row():
                                reads = result[goodpairs][
                                    result[goodpairs]["col4"].isin(gene_group)
                                ]
                                self.create_fusion_plot(gene_group, reads)

    def create_fusion_plot(self, title: str, reads: pd.DataFrame) -> None:
        """
        Creates an interactive plot showing read alignments for fusion candidates.

        Args:
            title: Plot title
            reads: DataFrame containing read alignments

        Performance Considerations:
        - Plot generation can be CPU-intensive
        - Large datasets may cause memory issues
        - Complex visualizations may impact UI responsiveness

        Potential Improvements:
        - Implement plot caching
        - Use parallel processing for plot generation
        - Add progressive loading for large datasets
        - Implement plot simplification for large datasets
        """
        with ui.card().classes("w-full no-shadow border-[2px]"):
            result = _get_reads(reads)
            df = reads

            # Function to rank overlapping ranges
            def rank_overlaps(df, start_col, end_col):
                # Sort by start and end columns
                df = df.sort_values(by=["gene", start_col, end_col]).reset_index(
                    drop=True
                )
                ranks = []
                current_rank = 0
                current_end = -1

                for _, row in df.iterrows():
                    if row[start_col] > current_end:
                        current_rank = 0
                    ranks.append(current_rank)
                    current_end = max(current_end, row[end_col])
                    current_rank += 1

                return ranks

            # Add logging for start2 and end2 columns before conversion
            logger.debug(
                f"Converting start2 column to int. First few values: {df['start2'].head()}"
            )
            logger.debug(
                f"Converting end2 column to int. First few values: {df['end2'].head()}"
            )

            try:
                df["start2"] = df["start2"].astype(int)
                df["end2"] = df["end2"].astype(int)
            except ValueError as e:
                logger.error(f"Error converting start2/end2 to int: {str(e)}")
                logger.error(
                    f"Problematic values in start2: {df[df['start2'].apply(lambda x: not str(x).isdigit())]['start2'].tolist()}"
                )
                logger.error(
                    f"Problematic values in end2: {df[df['end2'].apply(lambda x: not str(x).isdigit())]['end2'].tolist()}"
                )
                raise

            df = df.sort_values(by=["gene", "start2", "end2"])

            df["rank"] = rank_overlaps(df, "start2", "end2")

            lines = df.sort_values(by=["id", "read_start"]).reset_index(drop=True)

            # Function to assign occurrence number
            def assign_occurrence(group):
                group["Occurrence"] = range(1, len(group) + 1)
                return group

            # Apply the function to each group
            lines = (
                lines.groupby("id", observed=True)
                .apply(assign_occurrence)
                .reset_index(drop=True)
            )

            # Function to find join coordinates
            def find_joins(group):
                group = group.reset_index(drop=True)
                for i in range(len(group)):
                    if i < len(group) - 1:
                        if group.loc[i, "Occurrence"] == 1:
                            group.loc[i, "Join_Gene"] = group.loc[i + 1, "gene"]
                            group.loc[i, "Join_Start"] = group.loc[i, "start2"]
                            group.loc[i, "Join_Chromosome"] = group.loc[
                                i + 1, "chromosome2"
                            ]
                            group.loc[i, "spanB"] = group.loc[i + 1, "span"]
                            group.loc[i, "start3"] = group.loc[i + 1, "start2"]
                            group.loc[i, "end3"] = group.loc[i + 1, "end2"]
                            group.loc[i, "rankB"] = group.loc[i + 1, "rank"]
                            if group.loc[i + 1, "strand"] == "+":
                                group.loc[i, "Join_End"] = group.loc[i + 1, "end2"]
                            else:
                                group.loc[i, "Join_End"] = group.loc[i + 1, "start2"]
                return group

            # Initialize columns for the coordinates where the read joins the next read
            lines["Join_Start"] = None
            lines["Join_End"] = None
            lines["Join_Chromosome"] = None
            lines["Join_Gene"] = None
            lines["spanB"] = None
            lines["start3"] = None
            lines["end3"] = None
            lines["rankB"] = None

            # Apply the function to each group
            lines = (
                lines.groupby("id", observed=True)
                .apply(find_joins)
                .reset_index(drop=True)
            )
            # Remove rows containing NA values
            lines = lines.dropna()

            if len(result) > 1:
                gene_table = self.gene_table
                with ui.pyplot(figsize=(19, 5)).classes("w-full"):
                    # Create figure with tight layout
                    plt.rcParams["figure.constrained_layout.use"] = True
                    plt.rcParams["figure.constrained_layout.h_pad"] = 0.05
                    plt.rcParams["figure.constrained_layout.w_pad"] = 0.05

                    num_plots = 2 * len(result)
                    num_cols = len(result)
                    num_rows = (num_plots + num_cols - 1) // num_cols

                    for i, ax in enumerate(range(num_plots), start=1):
                        ax = plt.subplot(num_rows, num_cols, i)
                        # Remove extra padding around subplot
                        ax.set_position(
                            [
                                ax.get_position().x0,
                                ax.get_position().y0,
                                ax.get_position().width * 1.1,
                                ax.get_position().height * 1.1,
                            ]
                        )

                        row, col = divmod(i - 1, num_cols)
                        data = result.iloc[col]

                        chrom = data["chromosome"]
                        start = data["start"]
                        end = data["end"]

                        def human_readable_format(x, pos):
                            return f"{x / 1e6:.2f}"

                        if row == 1:
                            features = []
                            for index, row in gene_table[
                                gene_table["Seqid"].eq(chrom)
                                & gene_table["Start"].le(end)
                                & gene_table["End"].ge(start)
                            ].iterrows():
                                if row["Type"] == "gene":
                                    features.append(
                                        GraphicFeature(
                                            start=int(row["Start"]),
                                            end=int(row["End"]),
                                            strand=STRAND[row["Strand"]],
                                            thickness=8,
                                            color="#ffd700",
                                            label=row["gene_name"],
                                            fontdict={
                                                "color": "black",
                                                "fontsize": 8,
                                            },
                                        )
                                    )

                            for index, row in (
                                gene_table[
                                    gene_table["gene_name"].eq(data["gene"])
                                    & gene_table["Source"].eq("HAVANA")
                                    & gene_table["Type"].eq("exon")
                                ]
                                .groupby(["Seqid", "Start", "End", "Type", "Strand"])
                                .count()
                                .reset_index()
                                .iterrows()
                            ):
                                features.append(
                                    GraphicFeature(
                                        start=int(row["Start"]),
                                        end=int(row["End"]),
                                        strand=STRAND[row["Strand"]],
                                        thickness=4,
                                        color="#C0C0C0",
                                    )
                                )

                            record = GraphicRecord(
                                sequence_length=end - start,
                                first_index=start,
                                features=features,
                            )
                            ax = plt.gca()
                            record.plot(
                                ax=ax,
                                with_ruler=False,
                                draw_line=True,
                                strand_in_label_threshold=4,
                            )
                            # Adjust subplot spacing
                            ax.margins(x=0.02, y=0.1)

                        else:
                            features2 = []
                            for index, row in (
                                df[df["chromosome"].eq(chrom)]
                                .sort_values(by="id")
                                .iterrows()
                            ):
                                features2.append(
                                    GraphicFeature(
                                        start=int(row["start2"]),
                                        end=int(row["end2"]),
                                        strand=STRAND[row["strand"]],
                                        color=row["color"],
                                    )
                                )

                            record2 = GraphicRecord(
                                sequence_length=end - start,
                                first_index=start,
                                features=features2,
                            )
                            ax = plt.gca()
                            record2.plot(ax=ax)
                            ax.xaxis.set_major_formatter(
                                FuncFormatter(human_readable_format)
                            )
                            ax.tick_params(axis="x", labelsize=8)
                            ax.set_xlabel(
                                f'Position (Mb) - {chrom} - {data["gene"]}', fontsize=10
                            )
                            ax.set_title(f'{data["gene"]}', pad=2)
                            # Adjust subplot spacing
                            ax.margins(x=0.02, y=0.1)

                    # Adjust overall figure layout
                    plt.subplots_adjust(hspace=0.4, wspace=0.3)
                    gene_counter = Counter()

                    for index, row in lines.iterrows():
                        gene_counter[row["gene"]] += 1
                        gene_counter[row["Join_Gene"]] += 1

    def show_previous_data(self) -> None:
        """
        Loads and displays previously analyzed data.
        
        This method now only loads pre-processed data from FusionObject,
        eliminating all heavy computational work from the display layer.
        """
        if not self.browse:
            for item in app.storage.general[self.mainuuid]:
                if item == "sample_ids":
                    for sample in app.storage.general[self.mainuuid][item]:
                        self.sampleID = sample
            output = self.output
        if self.browse:
            output = self.check_and_create_folder(self.output, self.sampleID)

        # Load pre-processed fusion data (lightweight operation)
        self._load_preprocessed_fusion_data(output)
        
        # Load pre-processed structural variant data (lightweight operation)
        self._load_preprocessed_structural_variants(output)

    def _load_preprocessed_fusion_data(self, output: str) -> None:
        """Load pre-processed fusion data for display."""
        try:
            # Load target panel fusion data
            master_processed_file = os.path.join(output, "fusion_candidates_master_processed.csv")
            if os.path.exists(master_processed_file):
                import pickle
                with open(master_processed_file, 'rb') as f:
                    processed_data = pickle.load(f)
                
                # Update UI with pre-processed data
                self._update_fusion_ui_from_processed_data(processed_data, is_target_panel=True)
            
            # Load genome-wide fusion data
            all_processed_file = os.path.join(output, "fusion_candidates_all_processed.csv")
            if os.path.exists(all_processed_file):
                import pickle
                with open(all_processed_file, 'rb') as f:
                    processed_data = pickle.load(f)
                
                # Update UI with pre-processed data
                self._update_fusion_ui_from_processed_data(processed_data, is_target_panel=False)
                
        except Exception as e:
            logger.error(f"Error loading pre-processed fusion data: {str(e)}")

    def _load_preprocessed_structural_variants(self, output: str) -> None:
        """Load pre-processed structural variant data for display."""
        try:
            # Load structural variant count
            sv_count_file = os.path.join(output, "sv_count.txt")
            if os.path.exists(sv_count_file):
                with open(sv_count_file, 'r') as f:
                    self.sv_count = int(f.read().strip())
            
            # Load processed structural variant data
            sv_processed_file = os.path.join(output, "structural_variants_processed.csv")
            if os.path.exists(sv_processed_file):
                sv_df = pd.read_csv(sv_processed_file)
                
                # Update UI elements if they exist
                if hasattr(self, "sv_plot") and self.sv_plot is not None:
                    self.sv_plot.clear()
                    with self.sv_plot:
                        self.create_sv_plot(sv_df, self.sampleID)
                if (
                    hasattr(self, "sv_table_container")
                    and self.sv_table_container is not None
                ):
                    self.update_sv_table(sv_df)
                    
                logger.info(f"Loaded {len(sv_df)} pre-processed structural variant events")
            else:
                self.sv_count = 0
                
        except Exception as e:
            logger.error(f"Error loading pre-processed structural variants: {str(e)}")
            self.sv_count = 0

    def _update_fusion_ui_from_processed_data(self, processed_data: dict, is_target_panel: bool) -> None:
        """Update UI with pre-processed fusion data."""
        try:
            annotated_data = processed_data["annotated_data"]
            goodpairs = processed_data["goodpairs"]
            gene_groups = processed_data["gene_groups"]
            candidate_count = processed_data["candidate_count"]
            
            if is_target_panel:
                # Update target panel fusion table and UI
                self._update_fusion_table_from_processed(annotated_data, goodpairs, gene_groups, candidate_count)
            else:
                # Update genome-wide fusion table and UI
                self._update_fusion_table_all_from_processed(annotated_data, goodpairs, gene_groups, candidate_count)
                
        except Exception as e:
            logger.error(f"Error updating fusion UI from processed data: {str(e)}")

    def _update_fusion_table_from_processed(self, annotated_data: pd.DataFrame, goodpairs: pd.Series, gene_groups: list, candidate_count: int) -> None:
        """Update target panel fusion table with pre-processed data."""
        if annotated_data.shape[0] > self.fstable_row_count:
            self.fstable_row_count = annotated_data.shape[0]
            if not self.fstable:
                self.fusiontable.clear()
                with self.fusiontable:
                    self.fstable = (
                        ui.table.from_pandas(
                            annotated_data.sort_values(by="reference_start").rename(
                                columns={
                                    "col1": "chromBED",
                                    "col2": "BS",
                                    "col3": "BE",
                                    "col4": "Gene",
                                    "reference_id": "chrom",
                                    "reference_start": "mS",
                                    "reference_end": "mE",
                                    "read_id": "readID",
                                    "mapping_quality": "mapQ",
                                    "strand": "strand",
                                    "read_start": "Read Map Start",
                                    "read_end": "Read Map End",
                                    "is_secondary": "Secondary",
                                    "is_supplementary": "Supplementary",
                                    "mapping_span": "mapping span",
                                }
                            ),
                            pagination=25,
                        )
                        .props("dense")
                        .classes("w-full")
                        .style("height: 900px")
                        .style("font-size: 100%; font-weight: 300")
                    )
                    for col in self.fstable.columns:
                        col["sortable"] = True

                    with self.fstable.add_slot("top-right"):
                        with ui.input(placeholder="Search").props(
                            "type=search"
                        ).bind_value(self.fstable, "filter").add_slot("append"):
                            ui.icon("search")

                self.fusionplot.clear()
            else:
                self.fstable.update_rows(
                    annotated_data.sort_values(by="reference_start")
                    .rename(
                        columns={
                            "col1": "chromBED",
                            "col2": "BS",
                            "col3": "BE",
                            "col4": "Gene",
                            "reference_id": "chrom",
                            "reference_start": "mS",
                            "reference_end": "mE",
                            "read_id": "readID",
                            "mapping_quality": "mapQ",
                            "strand": "strand",
                            "read_start": "Read Map Start",
                            "read_end": "Read Map End",
                            "is_secondary": "Secondary",
                            "is_supplementary": "Supplementary",
                            "mapping_span": "mapping span",
                        }
                    )
                    .to_dict("records")
                )

                self.fusionplot.clear()

            self.candidates = candidate_count

            if not annotated_data.empty and candidate_count > 0:
                with self.fusionplot.classes("w-full"):
                    with ui.row().classes("w-full"):
                        ui.select(
                            options=gene_groups,
                            with_input=True,
                            on_change=lambda e: self._show_gene_pair_from_processed(
                                e.value, annotated_data, goodpairs, is_target_panel=True
                            ),
                        ).classes("w-40")
                    with ui.row().classes("w-full"):
                        self.card = ui.card()
                        with self.card:
                            ui.label("Select gene pair to see results.").tailwind(
                                "drop-shadow", "font-bold"
                            )

    def _update_fusion_table_all_from_processed(self, annotated_data: pd.DataFrame, goodpairs: pd.Series, gene_groups: list, candidate_count: int) -> None:
        """Update genome-wide fusion table with pre-processed data."""
        if annotated_data.shape[0] > self.fstable_all_row_count:
            self.fstable_all_row_count = annotated_data.shape[0]
            if not self.fstable_all:
                self.fusiontable_all.clear()
                with self.fusiontable_all:
                    self.fstable_all = (
                        ui.table.from_pandas(
                            annotated_data.sort_values(by="reference_start").rename(
                                columns={
                                    "col1": "chromBED",
                                    "col2": "BS",
                                    "col3": "BE",
                                    "col4": "Gene",
                                    "reference_id": "chrom",
                                    "reference_start": "mS",
                                    "reference_end": "mE",
                                    "read_id": "readID",
                                    "mapping_quality": "mapQ",
                                    "strand": "strand",
                                    "read_start": "Read Map Start",
                                    "read_end": "Read Map End",
                                    "is_secondary": "Secondary",
                                    "is_supplementary": "Supplementary",
                                    "mapping_span": "mapping span",
                                }
                            ),
                            pagination=25,
                        )
                        .props("dense")
                        .classes("w-full")
                        .style("height: 900px")
                        .style("font-size: 100%; font-weight: 300")
                    )
                    for col in self.fstable_all.columns:
                        col["sortable"] = True

                    with self.fstable_all.add_slot("top-right"):
                        with ui.input(placeholder="Search").props(
                            "type=search"
                        ).bind_value(self.fstable_all, "filter").add_slot("append"):
                            ui.icon("search")

                self.fusionplot_all.clear()
            else:
                self.fstable_all.update_rows(
                    annotated_data.sort_values(by="reference_start")
                    .rename(
                        columns={
                            "col1": "chromBED",
                            "col2": "BS",
                            "col3": "BE",
                            "col4": "Gene",
                            "reference_id": "chrom",
                            "reference_start": "mS",
                            "reference_end": "mE",
                            "read_id": "readID",
                            "mapping_quality": "mapQ",
                            "strand": "strand",
                            "read_start": "Read Map Start",
                            "read_end": "Read Map End",
                            "is_secondary": "Secondary",
                            "is_supplementary": "Supplementary",
                            "mapping_span": "mapping span",
                        }
                    )
                    .to_dict("records")
                )

                self.fusionplot_all.clear()

            self.all_candidates = candidate_count

            if not annotated_data.empty and candidate_count > 0:
                with self.fusionplot_all.classes("w-full"):
                    with ui.row().classes("w-full"):
                        ui.select(
                            options=gene_groups,
                            with_input=True,
                            on_change=lambda e: self._show_gene_pair_from_processed(
                                e.value, annotated_data, goodpairs, is_target_panel=False
                            ),
                        ).classes("w-40")
                    with ui.row().classes("w-full"):
                        self.all_card = ui.card()
                        with self.all_card:
                            ui.label("Select gene pair to see results.").tailwind(
                                "drop-shadow", "font-bold"
                            )

    def _show_gene_pair_from_processed(self, gene_group: str, annotated_data: pd.DataFrame, goodpairs: pd.Series, is_target_panel: bool) -> None:
        """Show gene pair plot from pre-processed data."""
        ui.notify(gene_group)
        
        if is_target_panel:
            self.card.clear()
            with self.card:
                with ui.row():
                    ui.label(f"{gene_group}").tailwind("drop-shadow", "font-bold")
                with ui.row():
                    reads = annotated_data[goodpairs][
                        annotated_data[goodpairs]["col4"].isin(gene_group)
                    ]
                    self.create_fusion_plot(gene_group, reads)
        else:
            self.all_card.clear()
            with self.all_card:
                with ui.row():
                    ui.label(f"{gene_group}").tailwind("drop-shadow", "font-bold")
                with ui.row():
                    reads = annotated_data[goodpairs][
                        annotated_data[goodpairs]["col4"].isin(gene_group)
                    ]
                    self.create_fusion_plot(gene_group, reads)

    def create_sv_plot(self, sv_df: pd.DataFrame, sample_id: str) -> None:
        """
        Creates a plot showing structural variant events.

        Args:
            sv_df: DataFrame containing structural variant data
            sample_id: Sample identifier
        """
        if sv_df.empty:
            with ui.column().classes("gap-2"):
                ui.label("No Structural Variants Found").classes("text-lg font-medium")
                ui.label(
                    "No structural variant events with sufficient support were identified"
                ).classes("text-gray-600")
            return

        with ui.pyplot(figsize=(12, 8)).classes("w-full"):
            # Create figure with tight layout
            plt.rcParams["figure.constrained_layout.use"] = True

            # Create subplots
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

            # Plot 1: Event type distribution
            event_counts = sv_df["Event Type"].value_counts()
            ax1.bar(event_counts.index, event_counts.values, color="skyblue", alpha=0.7)
            ax1.set_title(
                f"Structural Variant Events by Type - {sample_id}",
                fontsize=14,
                fontweight="bold",
            )
            ax1.set_ylabel("Number of Events", fontsize=12)
            ax1.tick_params(axis="x", rotation=45)

            # Add value labels on bars
            for i, v in enumerate(event_counts.values):
                ax1.text(
                    i, v + 0.1, str(v), ha="center", va="bottom", fontweight="bold"
                )

            # Plot 2: Support count distribution
            support_counts = sv_df["Support Count"].value_counts().sort_index()
            ax2.bar(
                support_counts.index.astype(str),
                support_counts.values,
                color="lightcoral",
                alpha=0.7,
            )
            ax2.set_title("Events by Support Count", fontsize=14, fontweight="bold")
            ax2.set_xlabel("Support Count", fontsize=12)
            ax2.set_ylabel("Number of Events", fontsize=12)

            # Add value labels on bars
            for i, v in enumerate(support_counts.values):
                ax2.text(
                    i, v + 0.1, str(v), ha="center", va="bottom", fontweight="bold"
                )

            # plt.tight_layout()

    def update_sv_table(self, sv_df: pd.DataFrame) -> None:
        """
        Updates the structural variant table in the UI.

        Args:
            sv_df: DataFrame containing structural variant data
        """
        if sv_df.empty:
            with ui.column().classes("gap-2"):
                ui.label("No Structural Variants Found").classes("text-lg font-medium")
                ui.label(
                    "No structural variant events with sufficient support were identified"
                ).classes("text-gray-600")
            return

        # Clear existing table
        self.sv_table_container.clear()

        with self.sv_table_container:
            # Create table with structural variant data
            self.sv_table = (
                ui.table.from_pandas(
                    sv_df,
                    pagination=25,
                )
                .props("dense")
                .classes("w-full")
                .style("height: 600px")
                .style("font-size: 100%; font-weight: 300")
            )

            # Make columns sortable
            for col in self.sv_table.columns:
                col["sortable"] = True

            # Add search functionality
            with self.sv_table.add_slot("top-right"):
                with ui.input(placeholder="Search structural variants").props(
                    "type=search"
                ).bind_value(self.sv_table, "filter").add_slot("append"):
                    ui.icon("search")


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
            counts_all = doubles_all.groupby("read_id")["col4"].transform("nunique")
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
            counts = doubles.groupby("read_id")["col4"].transform("nunique")
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
        Asynchronously processes BAM files with throttled breakpoint graph processing.
        """
        import time
        start_total = time.time()
        print(f"\n=== FUSION PROCESS_BAM START: {bamfile} ===")
        
        state.set_process_state("Fusion Analysis", ProcessState.RUNNING)
        try:
            try:
                logger.info(f"Starting BAM processing for file: {bamfile}")

                # Check for supplementary alignments
                start_supp_check = time.time()
                has_supp = has_supplementary(bamfile)
                supp_check_time = time.time() - start_supp_check
                print(f"  Supplementary check: {supp_check_time:.3f}s")

                if has_supp:
                    try:
                        # Process fusion candidates
                        start_fusion = time.time()
                        logger.info(f"Processing BAM file for fusions: {bamfile}")

                        # Execute fusion processing in background task
                        try:
                            async def fusion_background_work():
                                loop = asyncio.get_event_loop()
                                return await loop.run_in_executor(
                                    None, 
                                    fusion_work_pysam,
                                    bamfile,
                                    self.gene_bed,
                                    self.all_gene_bed,
                                )
                            
                            fusion_candidates, fusion_candidates_all = await background_tasks.create(
                                fusion_background_work()
                            )
                            fusion_time = time.time() - start_fusion
                            print(f"  Fusion processing: {fusion_time:.3f}s")
                        except Exception as e:
                            logger.error(f"Error in fusion processing: {str(e)}")
                            fusion_candidates, fusion_candidates_all = None, None
                            fusion_time = time.time() - start_fusion
                            print(f"  Fusion processing (ERROR): {fusion_time:.3f}s")
                        
                        # Store fusion candidates in the class dictionaries
                        start_store_fusion = time.time()
                        if (
                            fusion_candidates is not None
                            and not fusion_candidates.empty
                        ):
                            if self.sampleID not in self.fusion_candidates:
                                self.fusion_candidates[self.sampleID] = (
                                    fusion_candidates
                                )
                            else:
                                # Append to existing data
                                self.fusion_candidates[self.sampleID] = pd.concat(
                                    [
                                        self.fusion_candidates[self.sampleID],
                                        fusion_candidates,
                                    ],
                                    ignore_index=True,
                                )
                                
                            logger.info(
                                f"Added {len(fusion_candidates)} fusion candidates for sample {self.sampleID}"
                            )
                        
                        if (
                            fusion_candidates_all is not None
                            and not fusion_candidates_all.empty
                        ):
                            if self.sampleID not in self.fusion_candidates_all:
                                self.fusion_candidates_all[self.sampleID] = (
                                    fusion_candidates_all
                                )
                            else:
                                # Append to existing data
                                self.fusion_candidates_all[self.sampleID] = pd.concat(
                                    [
                                        self.fusion_candidates_all[self.sampleID],
                                        fusion_candidates_all,
                                    ],
                                    ignore_index=True,
                                )
                            logger.info(
                                f"Added {len(fusion_candidates_all)} genome-wide fusion candidates for sample {self.sampleID}"
                            )
                        store_fusion_time = time.time() - start_store_fusion
                        print(f"  Store fusion candidates: {store_fusion_time:.3f}s")

                        # Save fusion results with pre-processing for display
                        start_save_fusion = time.time()
                        await self._save_fusion_results()
                        save_fusion_time = time.time() - start_save_fusion
                        print(f"  Save fusion results: {save_fusion_time:.3f}s")
                        
                        # Process genome-wide structural variants with performance monitoring
                        start_sv = time.time()
                        logger.info(
                            f"Processing BAM file for structural variants: {bamfile}"
                        )

                        # Execute structural variant processing in background task
                        try:
                            async def sv_background_work():
                                loop = asyncio.get_event_loop()
                                return await loop.run_in_executor(None, process_bam_pipeline, bamfile)
                            
                            new_df = await background_tasks.create(sv_background_work())
                            sv_time = time.time() - start_sv
                            print(f"  Structural variant processing: {sv_time:.3f}s")
                            logger.info(f"Structural variant processing completed in {sv_time:.2f} seconds")
                        except Exception as e:
                            logger.error(f"Error in structural variant processing: {str(e)}")
                            new_df = pd.DataFrame()
                            sv_time = time.time() - start_sv
                            print(f"  Structural variant processing (ERROR): {sv_time:.3f}s")
                        
                        
                        # Save the new_df to a file, appending to existing data with optimized I/O
                        if not new_df.empty:
                            start_sv_save = time.time()
                            # Construct the output file path
                            output_dir = self.check_and_create_folder(
                                self.output, self.sampleID
                            )
                            sv_links_file = os.path.join(
                                output_dir, "structural_variant_links.csv"
                            )

                            try:
                                # Optimized file handling with reduced memory usage
                                if (
                                    os.path.exists(sv_links_file)
                                    and os.path.getsize(sv_links_file) > 0
                                ):
                                    # Read existing data with optimized dtypes using safe reader
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

                                    # Append new data efficiently
                                    combined_df = pd.concat(
                                        [existing_df, new_df], ignore_index=True, copy=False
                                    )
                                    logger.debug(
                                        f"Combined DataFrame size: {len(combined_df)} rows"
                                    )
                                else:
                                    # If file doesn't exist or is empty, use new data
                                    combined_df = new_df
                                    logger.debug(
                                        "No existing SV links file found, using new data"
                                    )

                                # Save combined data without compression for compatibility
                                combined_df.to_csv(sv_links_file, index=False)
                                sv_save_time = time.time() - start_sv_save
                                print(f"    Save SV links to CSV: {sv_save_time:.3f}s")

                                # Generate bed_lines from combined_df for BedTree - focus on breakpoint boundaries only
                                start_bed_lines = time.time()
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
                                            boundary_window = (
                                                500  # 500bp window around breakpoint
                                            )
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

                                                start2 = max(
                                                    0, coord2 - boundary_window
                                                )
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
                                                    start2 = max(
                                                        0, coord2 - boundary_window
                                                    )
                                                    end2 = coord2 + boundary_window

                                                    bed_line2_plus = f"{chrom1}\t{start2}\t{end2}\t{row['predominant_event']}_breakpoint2\t{row['support_count']}\t+"
                                                    bed_line2_minus = f"{chrom1}\t{start2}\t{end2}\t{row['predominant_event']}_breakpoint2\t{row['support_count']}\t-"
                                                    bed_lines.append(bed_line2_plus)
                                                    bed_lines.append(bed_line2_minus)

                                bed_lines_time = time.time() - start_bed_lines
                                print(f"    Generate BED lines: {bed_lines_time:.3f}s")
                                logger.info(
                                    f"Saved {len(combined_df)} SV links to {sv_links_file}"
                                )
                                logger.info(
                                    f"Generated {len(bed_lines)} breakpoint boundary BED lines for BedTree"
                                )

                                # Add to BedTree if needed
                                start_bedtree = time.time()
                                if (
                                    bed_lines
                                    and self.master_bed_tree[self.sampleID] is None
                                ):
                                    self.master_bed_tree.add_bed_tree(
                                        sample_id=self.sampleID,
                                        preserve_original_tree=True,
                                        reference_file=f"{self.reference_file}.fai",
                                    )
                                    bedfile = self.master_bed_tree.bed_trees[
                                        self.sampleID
                                    ]
                                    bedfile.load_from_file(self.bed_file)

                                if bed_lines:
                                    bedfile = self.master_bed_tree.bed_trees[
                                        self.sampleID
                                    ]
                                    bedfile.load_from_string(
                                        "\n".join(bed_lines),
                                        merge=False,
                                        write_files=True,
                                        output_location=os.path.join(
                                            self.check_and_create_folder(
                                                self.output, self.sampleID
                                            )
                                        ),
                                        source_type="FUSION",
                                    )
                                    bedtree_time = time.time() - start_bedtree
                                    print(f"    BedTree operations: {bedtree_time:.3f}s")
                                    logger.info(
                                        f"Added {len(bed_lines)} structural variant breakpoint boundaries to BedTree"
                                    )

                                    # Save fusion results again after structural variant processing
                                    start_save_final = time.time()
                                    await self._save_fusion_results()
                                    save_final_time = time.time() - start_save_final
                                    print(f"    Final save fusion results: {save_final_time:.3f}s")

                            except Exception as e:
                                logger.error(
                                    f"Error saving structural variant links: {str(e)}"
                                )
                                logger.error("Exception details:", exc_info=True)
                                # Continue processing even if save fails
                        else:
                            print(f"    No SV data to save")
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
                    print(f"  No supplementary alignments found, skipping")
                    logger.info("BAM file has no supplementary alignments, skipping.")
                    self.running = False

            except Exception as e:
                logger.error(f"Error in process_bam: {e}")
                raise
        finally:
            total_time = time.time() - start_total
            print(f"=== FUSION PROCESS_BAM END: {total_time:.3f}s total ===\n")
            state.set_process_state("Fusion Analysis", ProcessState.WAITING_FOR_DATA)

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

                async def final_bed_background_work(combined_sv_reads):
                    loop = asyncio.get_event_loop()
                    return await loop.run_in_executor(
                        None,
                        build_breakpoint_graph,
                        combined_sv_reads,
                        50000,  # max_proximity
                        True,   # group_by_sv
                    )

                bed_lines = await background_tasks.create(
                    final_bed_background_work(combined_sv_reads)
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
        import time
        start_total = time.time()
        print(f"      _save_fusion_results START")
        
        try:
            output_dir = self.check_and_create_folder(self.output, self.sampleID)
            
            # Save fusion candidates within target regions
            start_target = time.time()
            if (
                self.sampleID in self.fusion_candidates
                and not self.fusion_candidates[self.sampleID].empty
            ):
                # Filter for duplicated reads
                df = self.fusion_candidates[self.sampleID]
                uniques = df["read_id"].duplicated(keep=False)
                doubles = df[uniques]
                if not doubles.empty:
                    counts = doubles.groupby("read_id")["col4"].transform("nunique")
                    result = doubles[counts > 1]
                    if not result.empty:
                        # Save raw data
                        result.to_csv(
                            os.path.join(output_dir, "fusion_candidates_master.csv"),
                            index=False,
                        )
                        
                        # Pre-process data for FusionVis in background task
                        async def preprocess_master_background_work():
                            loop = asyncio.get_event_loop()
                            return await loop.run_in_executor(
                                None,
                                preprocess_fusion_data_standalone,
                                result, 
                                os.path.join(output_dir, "fusion_candidates_master_processed.csv")
                            )
                        
                        background_tasks.create(preprocess_master_background_work())
                        
                        logger.info(
                            f"Saved {len(result)} fusion candidates within target regions"
                        )
            target_time = time.time() - start_target
            print(f"        Target fusion processing: {target_time:.3f}s")

            # Save genome-wide fusion candidates
            start_all = time.time()
            if (
                self.sampleID in self.fusion_candidates_all
                and not self.fusion_candidates_all[self.sampleID].empty
            ):
                # Filter for duplicated reads
                df = self.fusion_candidates_all[self.sampleID]
                uniques_all = df["read_id"].duplicated(keep=False)
                doubles_all = df[uniques_all]
                if not doubles_all.empty:
                    counts_all = doubles_all.groupby("read_id")["col4"].transform(
                        "nunique"
                    )
                    result_all = doubles_all[counts_all > 1]
                    if not result_all.empty:
                        # Save raw data
                        result_all.to_csv(
                            os.path.join(output_dir, "fusion_candidates_all.csv"),
                            index=False,
                        )
                        
                        # Pre-process data for FusionVis in background task
                        async def preprocess_all_background_work():
                            loop = asyncio.get_event_loop()
                            return await loop.run_in_executor(
                                None,
                                preprocess_fusion_data_standalone,
                                result_all, 
                                os.path.join(output_dir, "fusion_candidates_all_processed.csv")
                            )
                        
                        background_tasks.create(preprocess_all_background_work())
                        
                        logger.info(
                            f"Saved {len(result_all)} genome-wide fusion candidates"
                        )
            all_time = time.time() - start_all
            print(f"        All fusion processing: {all_time:.3f}s")

            # Pre-process structural variant summary data in background task
            start_sv_preprocess = time.time()
            async def preprocess_sv_background_work():
                loop = asyncio.get_event_loop()
                return await loop.run_in_executor(None, preprocess_structural_variants_standalone, output_dir)
            
            background_tasks.create(preprocess_sv_background_work())
            sv_preprocess_time = time.time() - start_sv_preprocess
            print(f"        SV preprocess background task: {sv_preprocess_time:.3f}s")

        except Exception as e:
            logger.error(f"Error saving fusion results: {str(e)}")
            logger.error("Exception details:", exc_info=True)
        finally:
            total_time = time.time() - start_total
            print(f"      _save_fusion_results END: {total_time:.3f}s")

    def _preprocess_fusion_data_sync(
        self, 
        fusion_data: pd.DataFrame, 
        output_file: str
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
                "candidate_count": 0
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
                
                processed_data.update({
                    "gene_pairs": gene_pairs,
                    "gene_groups": gene_groups,
                    "candidate_count": len(gene_groups)
                })
            
            # Save processed data as pickle for efficient loading
            import pickle
            with open(output_file, 'wb') as f:
                pickle.dump(processed_data, f)
                
            logger.info(f"Pre-processed fusion data saved to {output_file}")
            
        except Exception as e:
            logger.error(f"Error pre-processing fusion data: {str(e)}")
            logger.error("Exception details:", exc_info=True)

    async def _preprocess_fusion_data_for_display(
        self, 
        fusion_data: pd.DataFrame, 
        output_file: str
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
                "candidate_count": 0
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
                
                processed_data.update({
                    "gene_pairs": gene_pairs,
                    "gene_groups": gene_groups,
                    "candidate_count": len(gene_groups)
                })
            
            # Save processed data as pickle for efficient loading
            import pickle
            with open(output_file, 'wb') as f:
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
                            os.path.join(output_dir, "structural_variants_processed.csv"),
                            index=False
                        )
                        
                        # Save count for quick access
                        with open(os.path.join(output_dir, "sv_count.txt"), 'w') as f:
                            f.write(str(len(sv_df)))
                        
                        logger.info(
                            f"Pre-processed {len(sv_df)} structural variant events"
                        )
                    else:
                        # Save empty count
                        with open(os.path.join(output_dir, "sv_count.txt"), 'w') as f:
                            f.write("0")
                else:
                    # Save empty count
                    with open(os.path.join(output_dir, "sv_count.txt"), 'w') as f:
                        f.write("0")
                        
        except Exception as e:
            logger.error(f"Error pre-processing structural variants: {str(e)}")
            logger.error("Exception details:", exc_info=True)

    async def _preprocess_structural_variants_for_display(self, output_dir: str) -> None:
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
                            os.path.join(output_dir, "structural_variants_processed.csv"),
                            index=False
                        )
                        
                        # Save count for quick access
                        with open(os.path.join(output_dir, "sv_count.txt"), 'w') as f:
                            f.write(str(len(sv_df)))
                        
                        logger.info(
                            f"Pre-processed {len(sv_df)} structural variant events"
                        )
                    else:
                        # Save empty count
                        with open(os.path.join(output_dir, "sv_count.txt"), 'w') as f:
                            f.write("0")
                else:
                    # Save empty count
                    with open(os.path.join(output_dir, "sv_count.txt"), 'w') as f:
                        f.write("0")
                        
        except Exception as e:
            logger.error(f"Error pre-processing structural variants: {str(e)}")
            logger.error("Exception details:", exc_info=True)


def test_me(
    port: int,
    threads: int,
    watchfolder: str,
    output: str,
    reload: bool = False,
    browse: bool = False,
) -> None:
    """
    Sets up and runs the gene fusion identification application.

    Args:
        port (int): Port number for the server.
        threads (int): Number of threads to use for processing.
        watchfolder (str): Path to the folder to watch for new BAM files.
        output (str): Path to the output directory.
        reload (bool): Flag to reload the application on changes.
        browse (bool): Flag to enable browsing historic data.
    """
    my_connection = None
    with theme.frame("Fusion Gene Identification.", my_connection):
        TestObject = FusionObject(threads, output, progress=True)
    if not browse:
        path = watchfolder
        searchdirectory = os.fsencode(path)
        for root, d_names, f_names in os.walk(searchdirectory):
            directory = os.fsdecode(root)
            for f in f_names:
                filename = os.fsdecode(f)
                if filename.endswith(".bam"):
                    TestObject.add_bam(os.path.join(directory, filename))
    else:
        TestObject.progress_trackers.visible = False
        TestObject.show_previous_data(output)
    ui.run(port=port, reload=reload)


@click.command()
@click.option(
    "--port",
    default=12345,
    help="Port for GUI",
)
@click.option("--threads", default=4, help="Number of threads available.")
@click.argument(
    "watchfolder",
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, resolve_path=True, path_type=Path
    ),
    required=False,
)
@click.argument(
    "output",
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, resolve_path=True, path_type=Path
    ),
    required=False,
)
@click.option(
    "--browse",
    is_flag=True,
    show_default=True,
    default=False,
    help="Browse Historic Data.",
)
def main(port, threads, watchfolder, output, browse) -> None:
    """
    CLI entry point for running the gene fusion identification app.

    Args:
        port (int): The port to serve the app on.
        threads (int): Number of threads available for processing.
        watchfolder (str): Directory to watch for new BAM files.
        output (str): Directory to save output files.
        browse (bool): Enable browsing historic data.
    """
    if browse:
        click.echo("Browse mode is enabled. Only the output folder is required.")
        test_me(
            port=port,
            reload=False,
            threads=threads,
            watchfolder=None,
            output=watchfolder,
            browse=browse,
        )
    else:
        click.echo(f"Watchfolder: {watchfolder}, Output: {output}")
        if watchfolder is None or output is None:
            click.echo("Watchfolder and output are required when --browse is not set.")
            sys.exit(1)
        test_me(
            port=port,
            reload=False,
            threads=threads,
            watchfolder=watchfolder,
            output=output,
            browse=browse,
        )


if __name__ in {"__main__", "__mp_main__"}:
    print("GUI launched by auto-reload function.")
    main()
