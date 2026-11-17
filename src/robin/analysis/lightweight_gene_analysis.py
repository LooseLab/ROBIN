q#!/usr/bin/env python3
"""
Lightweight Gene Analysis Module for robin

This module provides efficient, targeted analysis of genes of interest by:
1. Intersecting target BAM files with known pathogenic SNPs from ClinVar
2. Using samtools pileup for fast, targeted coverage analysis
3. Focusing on specific genes within the target panel
4. Providing lightweight reporting without heavy variant calling

Features:
- Fast targeted analysis using samtools pileup
- ClinVar pathogenic SNP intersection
- Gene-specific reporting
- Lightweight output formats
- Integration with robin's workflow system

Classes
-------
LightweightGeneAnalysis
    Main analysis class for lightweight gene analysis.

GeneReport
    Container for gene-specific analysis results.

Dependencies
-----------
- pysam: BAM file processing
- pandas: Data manipulation and analysis
- numpy: Numerical computations
- logging: Logging for debugging and monitoring
- typing: Type hints
- tempfile: Temporary file creation
- pathlib: File system paths
- os: Operating system interface
- subprocess: External command execution (samtools)

Example Usage
-----------
.. code-block:: python

    from robin.analysis.lightweight_gene_analysis import LightweightGeneAnalysis

    # Initialize analysis
    analysis = LightweightGeneAnalysis(
        work_dir="output/",
        clinvar_path="resources/clinvar.vcf.gz"
    )

    # Analyze specific genes
    results = analysis.analyze_genes("sample.bam", ["BRCA1", "TP53"])

Notes
-----
This module is designed to be much faster than full variant calling pipelines
while still providing valuable insights into genes of interest.
"""

import os
import logging
from typing import Dict, Any, Optional, List, Tuple
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pandas as pd
import pysam
from robin.logging_config import get_job_logger
from robin.analysis.variant_classification import is_clinvar_significant_from_mapping


@dataclass
class GeneReport:
    """Container for gene-specific analysis results."""

    gene_name: str
    gene_id: str
    chromosome: str
    start_position: int
    end_position: int
    total_coverage: int
    mean_coverage: float
    coverage_std: float
    pathogenic_snps_found: int
    pathogenic_snps_details: List[Dict[str, Any]]
    coverage_profile: List[Tuple[int, int]]  # (position, coverage)
    analysis_quality: str  # "high", "medium", "low"
    pileup_data: List[Dict[str, Any]] = None  # Added for backward compatibility


class LightweightGeneAnalysis:
    """
    Lightweight gene analysis using targeted samtools pileup and ClinVar intersection.

    This class provides efficient analysis of genes of interest by:
    1. Loading and indexing ClinVar pathogenic variants
    2. Performing targeted BAM pileup at gene regions
    3. Intersecting coverage with known pathogenic SNPs
    4. Generating lightweight reports
    """

    def __init__(
        self, work_dir: str, clinvar_path: str, reference: Optional[str] = None
    ):
        """
        Initialize the lightweight gene analysis.

        Args:
            work_dir: Working directory for output files
            clinvar_path: Path to ClinVar VCF file
            reference: Optional reference genome path
        """
        self.work_dir = Path(work_dir)
        self.clinvar_path = Path(clinvar_path)
        self.reference = reference
        self.logger = logging.getLogger("robin.lightweight_gene")

        # Ensure working directory exists
        self.work_dir.mkdir(parents=True, exist_ok=True)

        # Initialize ClinVar data
        self.clinvar_variants = {}
        self.gene_variants = {}

        # Load ClinVar data
        self._load_clinvar_data()

    def _load_clinvar_data(self) -> None:
        """Load and index ClinVar pathogenic variants."""
        print(f"Loading ClinVar pathogenic variants from: {self.clinvar_path}")
        self.logger.info("Loading ClinVar pathogenic variants...")

        if not self.clinvar_path.exists():
            print(f"ClinVar file not found: {self.clinvar_path}")
            self.logger.warning(f"ClinVar file not found: {self.clinvar_path}")
            return

        try:
            print("Parsing ClinVar VCF file...")
            # Handle compressed VCF files properly
            if str(self.clinvar_path).endswith(".gz"):
                import gzip

                with gzip.open(str(self.clinvar_path), "rt") as vcf_file:
                    self._parse_clinvar_vcf(vcf_file)
            else:
                with open(str(self.clinvar_path), "r") as vcf_file:
                    self._parse_clinvar_vcf(vcf_file)

            print(
                f" Successfully loaded {len(self.clinvar_variants)} pathogenic variants from ClinVar"
            )
            print(f"Variants found for {len(self.gene_variants)} genes")
            self.logger.info(
                f"Loaded {len(self.clinvar_variants)} pathogenic variants from ClinVar"
            )
            self.logger.info(f"Variants found for {len(self.gene_variants)} genes")

        except Exception as e:
            print(f"Error loading ClinVar data: {e}")
            self.logger.error(f"Error loading ClinVar data: {e}")
            import traceback

            self.logger.error(f"Traceback: {traceback.format_exc()}")

    def _parse_clinvar_vcf(self, vcf_file) -> None:
        """Parse ClinVar VCF file content."""
        pathogenic_count = 0
        total_lines = 0

        print("   Parsing VCF file line by line...")

        for line_num, line in enumerate(vcf_file, 1):
            line = line.strip()

            # Skip header lines
            if line.startswith("#"):
                continue

            total_lines += 1

            # Show progress every 10000 lines
            if total_lines % 10000 == 0:
                print(
                    f"    Processed {total_lines:,} lines, found {pathogenic_count} pathogenic variants..."
                )

            try:
                # Parse VCF line manually
                parts = line.split("\t")
                if (
                    len(parts) < 8
                ):  # Need at least CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO
                    continue

                chrom = parts[0]
                pos = int(parts[1])
                variant_id = parts[2]
                ref = parts[3]
                alt = parts[4]

                # Parse INFO field
                info_parts = parts[7].split(";")
                info_dict = {}
                for info_part in info_parts:
                    if "=" in info_part:
                        key, value = info_part.split("=", 1)
                        info_dict[key] = value
                    else:
                        info_dict[info_part] = True

                # Create a mock record object for compatibility
                class MockRecord:
                    def __init__(self, chrom, pos, variant_id, ref, alt, info_dict):
                        self.chrom = chrom
                        self.pos = pos
                        self.id = variant_id
                        self.ref = ref
                        self.alts = [alt] if alt != "." else []
                        self.info = info_dict

                record = MockRecord(chrom, pos, variant_id, ref, alt, info_dict)

                # Check if pathogenic
                if self._is_pathogenic(record):
                    pathogenic_count += 1
                    variant_info = self._extract_variant_info(record)

                    # Store by position for fast lookup
                    key = (record.chrom, record.pos)
                    self.clinvar_variants[key] = variant_info

                    # Store by gene for gene-specific analysis
                    if (
                        "gene_info" in variant_info
                        and variant_info["gene_info"] != "unknown"
                    ):
                        gene_name = variant_info["gene_info"].split(":")[0]
                        if gene_name not in self.gene_variants:
                            self.gene_variants[gene_name] = []
                        self.gene_variants[gene_name].append(variant_info)

            except Exception as e:
                if line_num < 100:  # Only log first 100 errors to avoid spam
                    print(f"   Warning parsing line {line_num}: {e}")
                    self.logger.warning(f"Error parsing line {line_num}: {e}")
                continue

        print(
            f"    VCF parsing completed: {total_lines:,} total lines, {pathogenic_count} pathogenic variants found"
        )

    def _is_pathogenic(self, record) -> bool:
        """Check if a ClinVar variant is clinically significant."""
        try:
            return is_clinvar_significant_from_mapping(record.info)
        except Exception:
            return False

    def _extract_variant_info(self, record) -> Dict[str, Any]:
        """Extract relevant information from a ClinVar variant record."""
        info = {
            "chromosome": record.chrom,
            "position": record.pos,
            "reference": record.ref,
            "alternate": record.alts[0] if record.alts else None,
            "variant_id": record.id,
            "clinical_significance": record.info.get("CLNSIG", "unknown"),
            "disease_name": record.info.get("CLNDN", "unknown"),
            "gene_info": record.info.get("GENEINFO", "unknown"),
            "molecular_consequence": record.info.get("MC", "unknown"),
            "review_status": record.info.get("CLNREVSTAT", "unknown"),
        }

        # Add oncogenicity if available
        if "ONC" in record.info:
            info["oncogenicity"] = record.info["ONC"]
        if "SCI" in record.info:
            info["somatic_impact"] = record.info["SCI"]

        return info

    def analyze_bed_regions(
        self, bam_path: str, bed_path: str
    ) -> Dict[str, GeneReport]:
        """
        Analyze BED regions by finding intersections with ClinVar pathogenic variants.

        Args:
            bam_path: Path to target BAM file
            bed_path: Path to BED file defining regions of interest

        Returns:
            Dictionary mapping region identifiers to GeneReport objects
        """
        print("Starting lightweight analysis of BED regions...")
        self.logger.info("Starting lightweight analysis of BED regions")

        if not os.path.exists(bam_path):
            raise FileNotFoundError(f"BAM file not found: {bam_path}")

        if not os.path.exists(bed_path):
            raise FileNotFoundError(f"BED file not found: {bed_path}")

        print(f"BAM file: {bam_path}")
        print(f"BED file: {bed_path}")

        # Load BED regions
        print("Loading BED regions...")
        bed_regions = self._load_bed_regions(bed_path)
        print(f"Loaded {len(bed_regions)} BED regions")
        self.logger.info(f"Loaded {len(bed_regions)} BED regions")

        # Find ClinVar variants that intersect with BED regions
        print("Finding ClinVar variants intersecting with BED regions...")
        intersecting_variants = self._find_bed_clinvar_intersections(bed_regions)
        print(
            f" Found {len(intersecting_variants)} ClinVar variants intersecting with BED regions"
        )
        self.logger.info(
            f"Found {len(intersecting_variants)} ClinVar variants intersecting with BED regions"
        )

        # Group variants by region for analysis
        print("Grouping variants by region...")
        region_variants = self._group_variants_by_region(
            intersecting_variants, bed_regions
        )

        results = {}
        total_regions = len([r for r in region_variants.values() if r])

        print(f"Starting analysis of {total_regions} regions with variants...")

        for i, (region_id, variants) in enumerate(region_variants.items(), 1):
            if not variants:
                continue

            print(
                f" [{i}/{total_regions}] Analyzing region '{region_id}' with {len(variants)} variants..."
            )
            self.logger.info(
                f"Analyzing region {region_id} with {len(variants)} variants"
            )

            try:
                # Analyze region coverage and variants
                region_report = self._analyze_single_region(
                    bam_path, region_id, variants
                )
                results[region_id] = region_report
                print(f"Region '{region_id}' completed successfully")

            except Exception as e:
                print(f"Error analyzing region '{region_id}': {e}")
                self.logger.error(f"Error analyzing region {region_id}: {e}")
                continue

        print(f"Analysis completed! Successfully analyzed {len(results)} regions")
        self.logger.info(f"Completed analysis of {len(results)} regions")
        return results

    def analyze_genes(
        self, bam_path: str, gene_names: List[str]
    ) -> Dict[str, GeneReport]:
        """
        Analyze specific genes in a target BAM file.

        Args:
            bam_path: Path to target BAM file
            gene_names: List of gene names to analyze

        Returns:
            Dictionary mapping gene names to GeneReport objects
        """
        self.logger.info(f"Starting lightweight analysis of {len(gene_names)} genes")

        if not os.path.exists(bam_path):
            raise FileNotFoundError(f"BAM file not found: {bam_path}")

        results = {}

        for gene_name in gene_names:
            self.logger.info(f"Analyzing gene: {gene_name}")

            try:
                # Get gene variants from ClinVar
                gene_variants = self.gene_variants.get(gene_name, [])

                if not gene_variants:
                    self.logger.warning(
                        f"No pathogenic variants found for gene: {gene_name}"
                    )
                    continue

                # Analyze gene coverage and variants
                gene_report = self._analyze_single_gene(
                    bam_path, gene_name, gene_variants
                )
                results[gene_name] = gene_report

            except Exception as e:
                self.logger.error(f"Error analyzing gene {gene_name}: {e}")
                continue

        self.logger.info(f"Completed analysis of {len(results)} genes")
        return results

    def _analyze_single_gene(
        self, bam_path: str, gene_name: str, gene_variants: List[Dict]
    ) -> GeneReport:
        """Analyze a single gene for coverage and pathogenic variants."""

        # Get gene coordinates from variants
        gene_coords = self._get_gene_coordinates(gene_variants)

        if not gene_coords:
            raise ValueError(f"Could not determine coordinates for gene: {gene_name}")

        # Perform targeted pileup using pysam
        pileup_data = self._perform_targeted_pileup(bam_path, gene_coords)

        # Extract coverage data from pileup for backward compatibility
        coverage_data = [
            (entry["position"], entry["coverage"]) for entry in pileup_data
        ]

        # Intersect with pathogenic variants
        intersecting_variants = self._find_intersecting_variants(
            gene_coords, gene_variants, pileup_data
        )

        # Calculate coverage statistics
        coverage_stats = self._calculate_coverage_statistics(pileup_data)

        # Determine analysis quality
        analysis_quality = self._determine_analysis_quality(pileup_data, gene_coords)

        # Extract gene ID safely
        gene_id = "unknown"
        if gene_variants and "gene_info" in gene_variants[0]:
            gene_info = gene_variants[0]["gene_info"]
            if gene_info != "unknown" and ":" in gene_info:
                gene_id = gene_info.split(":")[1]

        # Create gene report
        gene_report = GeneReport(
            gene_name=gene_name,
            gene_id=gene_id,
            chromosome=gene_coords["chromosome"],
            start_position=gene_coords["start"],
            end_position=gene_coords["end"],
            total_coverage=coverage_stats["total_coverage"],
            mean_coverage=coverage_stats["mean_coverage"],
            coverage_std=coverage_stats["coverage_std"],
            pathogenic_snps_found=len(intersecting_variants),
            pathogenic_snps_details=intersecting_variants,
            coverage_profile=coverage_data,
            pileup_data=pileup_data,
            analysis_quality=analysis_quality,
        )

        return gene_report

    def _load_bed_regions(self, bed_path: str) -> List[Dict[str, Any]]:
        """Load BED regions from file."""
        regions = []
        try:
            print(f"   Reading BED file: {bed_path}")
            with open(bed_path, "r") as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue

                    parts = line.split("\t")
                    if len(parts) >= 3:
                        region = {
                            "chromosome": parts[0],
                            "start": int(parts[1]),
                            "end": int(parts[2]),
                            "name": (
                                parts[3] if len(parts) > 3 else f"region_{line_num}"
                            ),
                            "line_num": line_num,
                        }
                        regions.append(region)

                        # Show first few regions for debugging
                        if line_num <= 5:
                            print(
                                f"   Region {line_num}: {region['chromosome']}:{region['start']}-{region['end']} '{region['name']}'"
                            )
                        elif line_num == 6:
                            print(f"   ... and {len(regions) - 5} more regions")

            print(f"   Successfully loaded {len(regions)} BED regions")
            self.logger.info(f"Loaded {len(regions)} BED regions")
            return regions

        except Exception as e:
            print(f"   Error loading BED file: {e}")
            self.logger.error(f"Error loading BED file: {e}")
            return []

    def _find_bed_clinvar_intersections(
        self, bed_regions: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """Find ClinVar variants that intersect with BED regions."""
        intersecting_variants = []
        total_variants = len(self.clinvar_variants)

        print(
            f"    Checking {total_variants} ClinVar variants against {len(bed_regions)} BED regions..."
        )

        # Show sample of BED regions for debugging
        print("   Sample BED regions:")
        for i, region in enumerate(bed_regions[:3]):
            print(
                f"      Region {i+1}: {region['chromosome']}:{region['start']:,}-{region['end']:,} '{region['name']}'"
            )
        if len(bed_regions) > 3:
            print(f"      ... and {len(bed_regions) - 3} more regions")

        # Show sample of ClinVar variants for debugging
        print("   🧬 Sample ClinVar variants:")
        variant_samples = list(self.clinvar_variants.items())[:3]
        for i, (variant_key, variant_info) in enumerate(variant_samples):
            chrom, pos = variant_key
            print(
                f"      Variant {i+1}: {chrom}:{pos:,} ({variant_info.get('gene_info', 'unknown')})"
            )
        if total_variants > 3:
            print(f"      ... and {total_variants - 3} more variants")

        for i, (variant_key, variant_info) in enumerate(
            self.clinvar_variants.items(), 1
        ):
            chrom, pos = variant_key

            # Show progress every 1000 variants
            if i % 1000 == 0:
                print(f"   Processed {i}/{total_variants} ClinVar variants...")

            # Check if this variant intersects with any BED region
            for region in bed_regions:
                # Handle coordinate system differences (chr1 vs 1)
                region_chrom = region["chromosome"]
                variant_chrom = chrom

                # Try to match chromosomes with and without 'chr' prefix
                chrom_match = (
                    region_chrom == variant_chrom  # Exact match
                    or (
                        region_chrom.startswith("chr")
                        and region_chrom[3:] == variant_chrom
                    )  # chr1 vs 1
                    or (
                        variant_chrom.startswith("chr")
                        and variant_chrom[3:] == region_chrom
                    )  # 1 vs chr1
                    or (
                        not region_chrom.startswith("chr")
                        and not variant_chrom.startswith("chr")
                        and region_chrom == variant_chrom
                    )  # Both without chr
                )

                if chrom_match and pos >= region["start"] and pos <= region["end"]:

                    # Create a copy of variant info with region information
                    variant_copy = variant_info.copy()
                    variant_copy["intersecting_region"] = region
                    variant_copy["bed_region_id"] = region["name"]
                    intersecting_variants.append(variant_copy)

                    # Show first few intersections for debugging
                    if len(intersecting_variants) <= 5:
                        print(
                            f"   Found intersection: {chrom}:{pos:,} in region '{region['name']}' ({variant_info.get('gene_info', 'unknown')})"
                        )
                        print(f"      Matched: {variant_chrom} ↔ {region_chrom}")
                    elif len(intersecting_variants) == 6:
                        print(
                            f"   ... and {len(intersecting_variants) - 5} more intersections"
                        )

                    break  # Found intersection, move to next variant

        if not intersecting_variants:
            print("   No intersections found! This might indicate:")
            print("      • Coordinate system mismatch (e.g., 'chr1' vs '1')")
            print("      • Different genome builds (e.g., hg38 vs hg19)")
            print("      • BED regions outside ClinVar variant ranges")
            print("      • ClinVar file doesn't contain variants for your regions")

            # Show coordinate system analysis
            bed_chroms = set(region["chromosome"] for region in bed_regions)
            clinvar_chroms = set(chrom for chrom, _ in self.clinvar_variants.keys())

            print("   Coordinate system analysis:")
            print(f"      BED chromosomes: {sorted(list(bed_chroms))[:10]}...")
            print(f"      ClinVar chromosomes: {sorted(list(clinvar_chroms))[:10]}...")

            # Check for potential matches
            potential_matches = []
            for bed_chrom in list(bed_chroms)[:5]:  # Check first 5 BED chromosomes
                for clinvar_chrom in list(clinvar_chroms)[
                    :5
                ]:  # Check first 5 ClinVar chromosomes
                    if (
                        bed_chrom.startswith("chr") and bed_chrom[3:] == clinvar_chrom
                    ) or (
                        clinvar_chrom.startswith("chr")
                        and clinvar_chrom[3:] == bed_chrom
                    ):
                        potential_matches.append((bed_chrom, clinvar_chrom))

            if potential_matches:
                print(f"      Potential coordinate matches: {potential_matches}")
            else:
                print("      No obvious coordinate system matches found")

        print(
            f"    Found {len(intersecting_variants)} ClinVar variants intersecting with BED regions"
        )
        self.logger.info(
            f"Found {len(intersecting_variants)} ClinVar variants intersecting with BED regions"
        )
        return intersecting_variants

    def _group_variants_by_region(
        self,
        intersecting_variants: List[Dict[str, Any]],
        bed_regions: List[Dict[str, Any]],
    ) -> Dict[str, List[Dict[str, Any]]]:
        """Group intersecting variants by BED region."""
        region_variants = {}

        # Initialize all regions
        for region in bed_regions:
            region_variants[region["name"]] = []

        # Group variants by region
        for variant in intersecting_variants:
            region_id = variant["bed_region_id"]
            if region_id in region_variants:
                region_variants[region_id].append(variant)

        return region_variants

    def _analyze_single_region(
        self, bam_path: str, region_id: str, variants: List[Dict[str, Any]]
    ) -> GeneReport:
        """Analyze a single BED region for coverage and pathogenic variants."""

        # Get region coordinates from the first variant's intersecting region
        if not variants:
            raise ValueError(f"No variants found for region: {region_id}")

        region_info = variants[0]["intersecting_region"]
        region_coords = {
            "chromosome": region_info["chromosome"],
            "start": region_info["start"],
            "end": region_info["end"],
        }

        # Perform targeted pileup using pysam

        pileup_data = self._perform_targeted_pileup(bam_path, region_coords)

        # Extract coverage data from pileup for backward compatibility
        coverage_data = [
            (entry["position"], entry["coverage"]) for entry in pileup_data
        ]

        # Intersect with pathogenic variants
        intersecting_variants = self._find_intersecting_variants(
            region_coords, variants, pileup_data
        )

        # Calculate coverage statistics
        coverage_stats = self._calculate_coverage_statistics(pileup_data)

        # Determine analysis quality
        analysis_quality = self._determine_analysis_quality(pileup_data, region_coords)

        # Create gene report (reusing the same structure for regions)
        region_report = GeneReport(
            gene_name=region_id,
            gene_id=region_id,  # Use region name as ID
            chromosome=region_coords["chromosome"],
            start_position=region_coords["start"],
            end_position=region_coords["end"],
            total_coverage=coverage_stats["total_coverage"],
            mean_coverage=coverage_stats["mean_coverage"],
            coverage_std=coverage_stats["coverage_std"],
            pathogenic_snps_found=len(intersecting_variants),
            pathogenic_snps_details=intersecting_variants,
            coverage_profile=coverage_data,
            pileup_data=pileup_data,
            analysis_quality=analysis_quality,
        )

        return region_report

    def _get_gene_coordinates(
        self, gene_variants: List[Dict]
    ) -> Optional[Dict[str, Any]]:
        """Extract gene coordinates from variant information."""
        if not gene_variants:
            return None

        positions = [v["position"] for v in gene_variants]
        chromosomes = set(v["chromosome"] for v in gene_variants)

        return {
            "chromosome": list(chromosomes)[0],
            "start": min(positions),
            "end": max(positions),
        }

    def _perform_targeted_pileup(
        self, bam_path: str, gene_coords: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        """Perform targeted pileup using pysam for a specific gene region."""

        print(f"         Opening BAM file: {bam_path}")

        # Open BAM file
        try:
            bam_file = pysam.AlignmentFile(bam_path, "rb")
        except Exception as e:
            print(f"         Error opening BAM file: {e}")
            self.logger.error(f"Error opening BAM file {bam_path}: {e}")
            return []

        # Check if BAM file has an index
        try:
            bam_file.check_index()
            print("         BAM file has valid index")
        except Exception as e:
            print(f"         BAM file missing or invalid index: {e}")
            print(f"         Try running: samtools index {bam_path}")
            self.logger.error(f"BAM file missing index: {e}")
            bam_file.close()
            return []

        pileup_data = []

        try:
            # Get the chromosome name as it appears in the BAM file
            # Handle different chromosome naming conventions (chr1 vs 1)
            chrom_name = gene_coords["chromosome"]
            if not chrom_name.startswith("chr") and any(
                ref.startswith("chr") for ref in bam_file.references
            ):
                chrom_name = f"chr{chrom_name}"
            elif chrom_name.startswith("chr") and not any(
                ref.startswith("chr") for ref in bam_file.references
            ):
                chrom_name = chrom_name[3:]

            print(f"         🧬 Using chromosome name: {chrom_name}")

            # Check if chromosome exists in BAM
            if chrom_name not in bam_file.references:
                print(f"         Chromosome {chrom_name} not found in BAM file")
                self.logger.warning(
                    f"Chromosome {chrom_name} not found in BAM file. Available: {list(bam_file.references)[:5]}..."
                )
                return []

            # Get pileup for the specific region
            start_pos = gene_coords["start"]
            end_pos = gene_coords["end"]

            print(
                f"         Analyzing region: {chrom_name}:{start_pos:,}-{end_pos:,}"
            )
            self.logger.debug(
                f"Analyzing pileup for {chrom_name}:{start_pos}-{end_pos}"
            )

            # Use pysam's pileup method for the region
            for pileup_column in bam_file.pileup(
                chrom_name, start_pos, end_pos, truncate=True, stepper="samtools"
            ):

                position = pileup_column.pos + 1  # pysam is 0-based, convert to 1-based

                # Skip positions outside our target region
                if position < start_pos or position > end_pos:
                    continue

                # Count bases at this position
                base_counts = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0, "indel": 0}
                quality_scores = []
                read_names = []

                for read in pileup_column.pileups:
                    if read.is_del or read.is_refskip:
                        base_counts["indel"] += 1
                        continue

                    if read.query_position is not None:
                        # Get the base and quality
                        base = read.alignment.query_sequence[read.query_position]
                        quality = (
                            read.alignment.query_qualities[read.query_position]
                            if read.alignment.query_qualities
                            else 0
                        )

                        # Count the base
                        if base.upper() in base_counts:
                            base_counts[base.upper()] += 1
                        else:
                            base_counts["N"] += 1

                        quality_scores.append(quality)
                        read_names.append(read.alignment.query_name)

                # Calculate coverage (total reads at this position)
                total_coverage = sum(base_counts.values()) - base_counts["indel"]

                # Determine reference base (we'll need to get this from the BAM header or assume)
                # For now, we'll use the most common base as reference
                ref_base = max((k, v) for k, v in base_counts.items() if k != "indel")[
                    0
                ]

                # Calculate variant allele frequency if there are variants
                variant_alleles = [
                    (base, count)
                    for base, count in base_counts.items()
                    if base != ref_base and base != "indel" and count > 0
                ]

                # Create pileup entry
                pileup_entry = {
                    "position": position,
                    "reference": ref_base,
                    "coverage": total_coverage,
                    "base_counts": base_counts,
                    "variant_alleles": variant_alleles,
                    "mean_quality": np.mean(quality_scores) if quality_scores else 0,
                    "read_names": read_names[
                        :10
                    ],  # Limit to first 10 reads to avoid memory issues
                    "has_variants": len(variant_alleles) > 0,
                }

                pileup_data.append(pileup_entry)

                # Log progress for large regions
                if len(pileup_data) % 1000 == 0:
                    print(f"         Processed {len(pileup_data)} positions...")
                    self.logger.debug(f"Processed {len(pileup_data)} positions...")

            print(
                f"          Pileup analysis completed: {len(pileup_data)} positions analyzed"
            )
            self.logger.info(
                f"Generated pileup data for {len(pileup_data)} positions in region {chrom_name}:{start_pos}-{end_pos}"
            )

        except Exception as e:
            print(f"         Error during pileup analysis: {e}")
            self.logger.error(f"Error during pileup analysis: {e}")
            import traceback

            self.logger.error(f"Traceback: {traceback.format_exc()}")

        finally:
            bam_file.close()

        return pileup_data

    def _find_intersecting_variants(
        self,
        gene_coords: Dict[str, Any],
        gene_variants: List[Dict],
        pileup_data: List[Dict[str, Any]],
    ) -> List[Dict[str, Any]]:
        """Find pathogenic variants that intersect with coverage data and analyze evidence."""

        # Create position lookup for fast access
        pileup_by_position = {entry["position"]: entry for entry in pileup_data}

        intersecting_variants = []

        for variant in gene_variants:
            if variant["position"] in pileup_by_position:
                # Get the pileup data for this position
                pileup_entry = pileup_by_position[variant["position"]]

                # Create enhanced variant copy with evidence analysis
                variant_copy = variant.copy()
                variant_copy["coverage_at_variant"] = pileup_entry["coverage"]

                # Analyze if the BAM data supports this ClinVar variant
                variant_copy["evidence_analysis"] = self._analyze_variant_evidence(
                    variant, pileup_entry
                )

                intersecting_variants.append(variant_copy)

        return intersecting_variants

    def _analyze_variant_evidence(
        self, clinvar_variant: Dict[str, Any], pileup_entry: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Analyze if the BAM data supports the ClinVar variant."""

        # Get the expected reference and alternate from ClinVar
        expected_ref = clinvar_variant.get("reference", "N")
        expected_alt = clinvar_variant.get("alternate", "N")

        # Get the actual bases found in the BAM
        base_counts = pileup_entry["base_counts"]
        total_coverage = pileup_entry["coverage"]

        # Determine if this is an indel
        is_indel = len(expected_ref) != len(expected_alt)
        variant_type = "indel" if is_indel else "snv"

        # For indels, we need to check if the reference sequence matches what we expect
        if is_indel:
            
            # NOTE: Indel analysis is complex and requires:
            # 1. Reference genome sequence validation
            # 2. Multi-position pileup analysis
            # 3. CIGAR string parsing for true insertion/deletion evidence
            # 4. Context-aware sequence alignment

            if variant_type == "indel":
                if len(expected_ref) > len(expected_alt):
                    # Deletion (e.g., TG→T)
                    # Check if we see the expected reference sequence at this position
                    # For TG→T, we need to see T at position and G at position+1
                    ref_support = 0
                    alt_support = 0

                    # This is a complex case - we need to check multiple positions
                    # For now, mark as requiring manual review
                    evidence_level = "indel_requires_review"
                    vaf = 0.0

                else:
                    # Insertion (e.g., T→TG)
                    # For insertions, we need to:
                    # 1. Verify the reference base is correct
                    # 2. Check if we actually see the insertion
                    # 3. Validate the insertion context

                    # Get the reference base from the pileup
                    ref_base = pileup_entry.get("reference", "N")
                    expected_ref_base = expected_ref[0] if expected_ref else "N"

                    # Check if reference base matches what we expect
                    if ref_base != expected_ref_base:
                        # Reference mismatch - this might not be the right position
                        evidence_level = "reference_mismatch"
                        ref_support = 0
                        alt_support = 0
                        vaf = 0.0
                    else:
                        
                        # For insertions, we need to look at the actual sequence context
                        # This requires more sophisticated pileup analysis
                        # For now, mark as requiring manual review
                        evidence_level = "insertion_requires_review"
                        ref_support = base_counts.get(expected_ref_base.upper(), 0)
                        alt_support = (
                            0  # Can't determine without proper insertion analysis
                        )
                        vaf = 0.0
            else:
                # Fallback for unknown indel types
                ref_support = 0
                alt_support = 0
                evidence_level = "indel_requires_review"
                vaf = 0.0
        else:
            # SNV - use the original logic
            ref_support = base_counts.get(expected_ref.upper(), 0)
            alt_support = base_counts.get(expected_alt.upper(), 0)
            vaf = alt_support / total_coverage if total_coverage > 0 else 0

            # Determine evidence level for SNVs
            if alt_support == 0:
                evidence_level = "no_support"
            elif alt_support < 3:
                evidence_level = "low_support"
            elif vaf < 0.1:
                evidence_level = "low_frequency"
            elif vaf < 0.3:
                evidence_level = "medium_frequency"
            else:
                evidence_level = "high_frequency"

        # Check if this looks like a heterozygous or homozygous variant
        zygosity = "unknown"
        if total_coverage > 0 and not is_indel:
            if vaf > 0.8:
                zygosity = "homozygous"
            elif vaf > 0.2:
                zygosity = "heterozygous"
            else:
                zygosity = "low_frequency"
        elif is_indel:
            zygosity = "indel_requires_review"

        return {
            "reference_support": ref_support,
            "alternate_support": alt_support,
            "variant_allele_frequency": vaf,
            "evidence_level": evidence_level,
            "zygosity": zygosity,
            "total_coverage": total_coverage,
            "base_counts": base_counts,
            "mean_quality": pileup_entry.get("mean_quality", 0),
            "has_variants": pileup_entry.get("has_variants", False),
            "variant_type": variant_type,
            "is_indel": is_indel,
        }

    def _calculate_coverage_statistics(
        self, pileup_data: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """Calculate coverage statistics from pileup data."""
        if not pileup_data:
            return {"total_coverage": 0, "mean_coverage": 0.0, "coverage_std": 0.0}

        coverages = [entry["coverage"] for entry in pileup_data]

        return {
            "total_coverage": sum(coverages),
            "mean_coverage": np.mean(coverages),
            "coverage_std": np.std(coverages),
        }

    def _determine_analysis_quality(
        self, pileup_data: List[Dict[str, Any]], gene_coords: Dict[str, Any]
    ) -> str:
        """Determine the quality of the analysis based on coverage and region size."""
        if not pileup_data:
            return "low"

        region_size = gene_coords["end"] - gene_coords["start"]
        mean_coverage = np.mean([entry["coverage"] for entry in pileup_data])

        # Quality thresholds
        if mean_coverage >= 30 and region_size >= 1000:
            return "high"
        elif mean_coverage >= 10 and region_size >= 500:
            return "medium"
        else:
            return "low"

    def generate_report(
        self, gene_reports: Dict[str, GeneReport], output_path: Optional[str] = None
    ) -> str:
        """
        Generate a comprehensive report from gene analysis results.

        Args:
            gene_reports: Dictionary of gene analysis results
            output_path: Optional output path for the report

        Returns:
            Path to the generated report
        """
        if output_path is None:
            output_path = self.work_dir / "lightweight_gene_analysis_report.html"

        print(f"         Generating HTML report for {len(gene_reports)} regions...")

        # Convert to DataFrame for easy manipulation
        report_data = []
        for gene_name, report in gene_reports.items():
            report_data.append(
                {
                    "Gene": gene_name,
                    "Gene_ID": report.gene_id,
                    "Chromosome": report.chromosome,
                    "Start": report.start_position,
                    "End": report.end_position,
                    "Total_Coverage": report.total_coverage,
                    "Mean_Coverage": f"{report.mean_coverage:.1f}",
                    "Coverage_STD": f"{report.coverage_std:.1f}",
                    "Pathogenic_SNPs": report.pathogenic_snps_found,
                    "Analysis_Quality": report.analysis_quality,
                }
            )

        df = pd.DataFrame(report_data)

        # Generate HTML report
        html_content = self._generate_html_report(df, gene_reports)

        with open(output_path, "w") as f:
            f.write(html_content)

        print(f"         HTML report generated: {output_path}")
        self.logger.info(f"Generated report: {output_path}")
        return str(output_path)

    def _generate_html_report(
        self, df: pd.DataFrame, gene_reports: Dict[str, GeneReport]
    ) -> str:
        """Generate an HTML report from the analysis results."""

        html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Lightweight Gene Analysis Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
                .summary {{ margin: 20px 0; }}
                .gene-detail {{ margin: 20px 0; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }}
                .pathogenic-present {{ background-color: #ffe6e6; border-left: 4px solid #ff4444; }}
                .pathogenic-absent {{ background-color: #e6ffe6; border-left: 4px solid #44ff44; }}
                .low-coverage {{ background-color: #fff2e6; border-left: 4px solid #ffaa44; }}
                .coverage-high {{ color: #006600; font-weight: bold; }}
                .coverage-medium {{ color: #cc6600; font-weight: bold; }}
                .coverage-low {{ color: #cc0000; font-weight: bold; }}
                .variant-detail {{ margin: 10px 0; padding: 10px; background-color: #f9f9f9; border-radius: 3px; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>🧬 Lightweight Gene Analysis Report</h1>
                <p>Generated by robin Lightweight Gene Analysis Module</p>
                <p>Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            </div>
            
            <div class="summary">
                <h2> Analysis Summary</h2>
                <p><strong>Total Genes Analyzed:</strong> {len(gene_reports)}</p>
                <p><strong>Genes with Pathogenic SNPs:</strong> {sum(1 for r in gene_reports.values() if r.pathogenic_snps_found > 0)}</p>
                <p><strong>Total Pathogenic SNPs Found:</strong> {sum(r.pathogenic_snps_found for r in gene_reports.values())}</p>
                <p><strong>Genes with Low Coverage (&lt;10x):</strong> {sum(1 for r in gene_reports.values() if r.mean_coverage < 10)}</p>
            </div>
            
            <h2> Detailed Gene Analysis</h2>
        """

        # Group variants by gene for better reporting
        gene_variants = {}
        for gene_name, report in gene_reports.items():
            if report.pathogenic_snps_details:
                for variant in report.pathogenic_snps_details:
                    gene_info = variant.get("gene_info", "unknown")
                    if gene_info != "unknown" and ":" in gene_info:
                        gene_name_actual = gene_info.split(":")[0]
                        if gene_name_actual not in gene_variants:
                            gene_variants[gene_name_actual] = []
                        gene_variants[gene_name_actual].append(
                            {
                                "variant": variant,
                                "coverage": variant.get("coverage_at_variant", 0),
                                "region": gene_name,
                                "report": report,
                            }
                        )
                    else:
                        # Fallback to region name if no gene info
                        if gene_name not in gene_variants:
                            gene_variants[gene_name] = []
                        gene_variants[gene_name].append(
                            {
                                "variant": variant,
                                "coverage": variant.get("coverage_at_variant", 0),
                                "region": gene_name,
                                "report": report,
                            }
                        )

        # Add detailed gene information
        for gene_name, variants in gene_variants.items():
            # Calculate gene-level statistics
            total_variants = len(variants)
            high_coverage_variants = sum(1 for v in variants if v["coverage"] >= 10)
            low_coverage_variants = total_variants - high_coverage_variants

            # Get coverage statistics from the first report (they should be similar for the same gene)
            first_report = variants[0]["report"]

            # Determine coverage class and status
            if first_report.mean_coverage < 10:
                coverage_class = "coverage-low"
                coverage_status = " LOW COVERAGE"
                detail_class = "low-coverage"
            elif first_report.mean_coverage < 30:
                coverage_class = "coverage-medium"
                coverage_status = " MEDIUM COVERAGE"
                detail_class = "gene-detail"
            else:
                coverage_class = "coverage-high"
                coverage_status = " GOOD COVERAGE"
                detail_class = "gene-detail"

            # Determine pathogenic SNP status
            if high_coverage_variants > 0:
                pathogenic_status = f"🚨 PATHOGENIC SNPS PRESENT ({high_coverage_variants} with good coverage)"
                detail_class += " pathogenic-present"
            else:
                pathogenic_status = (
                    " NO PATHOGENIC SNPS DETECTED (or all below coverage threshold)"
                )
                detail_class += " pathogenic-absent"

            html += f"""
            <div class="{detail_class}">
                <h3>🧬 {gene_name}</h3>
                <p><strong>Status:</strong> {pathogenic_status}</p>
                <p><strong>Coverage:</strong> <span class="{coverage_class}">{coverage_status} - {first_report.mean_coverage:.1f}x mean</span></p>
                <p><strong>Total Pathogenic Variants:</strong> {total_variants} (including {low_coverage_variants} below 10x coverage)</p>
                <p><strong>Regions Analyzed:</strong> {', '.join(set(v['region'] for v in variants))}</p>
                <p><strong>Coverage Statistics:</strong></p>
                <ul>
                    <li>Total Coverage: {first_report.total_coverage:,} bases</li>
                    <li>Mean Coverage: {first_report.mean_coverage:.1f}x</li>
                    <li>Coverage Standard Deviation: {first_report.coverage_std:.1f}x</li>
                </ul>
            """

            if variants:
                html += f"<p><strong>🚨 Pathogenic Variants Found ({total_variants}):</strong></p>"

                # Group variants by coverage status
                high_coverage_group = [v for v in variants if v["coverage"] >= 10]
                low_coverage_group = [v for v in variants if v["coverage"] < 10]

                if high_coverage_group:
                    html += f"<h4> High Coverage Variants (≥10x) - {len(high_coverage_group)} variants</h4>"
                    for variant_data in high_coverage_group:
                        variant = variant_data["variant"]
                        coverage = variant_data["coverage"]
                        html += f"""
                        <div class="variant-detail">
                            <strong>Position {variant['position']:,}:</strong> {variant['reference']}→{variant['alternate']}
                            <br>Clinical Significance: {variant['clinical_significance']}
                            <br>Disease: {variant['disease_name']}
                            <br>Coverage at variant: <span class="coverage-high">{coverage}x</span>
                        """

                        # Add evidence analysis if available
                        if "evidence_analysis" in variant:
                            evidence = variant["evidence_analysis"]
                            html += f"""
                            <br><strong>Evidence Analysis:</strong>
                            <br>• Reference support: {evidence['reference_support']} reads
                            <br>• Alternate support: {evidence['alternate_support']} reads
                            <br>• Variant allele frequency: {evidence['variant_allele_frequency']:.1%}
                            <br>• Evidence level: {evidence['evidence_level'].replace('_', ' ').title()}
                            <br>• Zygosity: {evidence['zygosity'].title()}
                            <br>• Mean quality: {evidence['mean_quality']:.1f}
                            """

                        html += "</div>"

                if low_coverage_group:
                    html += f"<h4> Low Coverage Variants (<10x) - {len(low_coverage_group)} variants</h4>"
                    html += "<p><em>These variants are excluded from detailed analysis due to insufficient coverage:</em></p>"
                    for variant_data in low_coverage_group:
                        variant = variant_data["variant"]
                        coverage = variant_data["coverage"]
                        html += f"""
                        <div class="variant-detail">
                            <strong>Position {variant['position']:,}:</strong> {variant['reference']}→{variant['alternate']}
                            <br>Clinical Significance: {variant['clinical_significance']}
                            <br>Disease: {variant['disease_name']}
                            <br>Coverage at variant: <span class="coverage-low">{coverage}x (LOW COVERAGE - EXCLUDED)</span>
                        </div>
                        """
            else:
                html += "<p><strong> No pathogenic variants detected in this gene</strong></p>"

            html += "</div>"

        html += """
        </body>
        </html>
        """

        return html

    def export_json(
        self, gene_reports: Dict[str, GeneReport], output_path: Optional[str] = None
    ) -> str:
        """
        Export gene analysis results to JSON format for programmatic analysis.

        Args:
            gene_reports: Dictionary of gene analysis results
            output_path: Optional output path for the JSON file

        Returns:
            Path to the generated JSON file
        """
        if output_path is None:
            output_path = self.work_dir / "lightweight_gene_analysis_results.json"

        print(f"         Generating JSON report for {len(gene_reports)} regions...")

        # Group variants by gene for better JSON structure
        gene_variants = {}
        for gene_name, report in gene_reports.items():
            if report.pathogenic_snps_details:
                for variant in report.pathogenic_snps_details:
                    gene_info = variant.get("gene_info", "unknown")
                    if gene_info != "unknown" and ":" in gene_info:
                        gene_name_actual = gene_info.split(":")[0]
                        if gene_name_actual not in gene_variants:
                            gene_variants[gene_name_actual] = {
                                "gene_name": gene_name_actual,
                                "regions_analyzed": set(),
                                "coverage_statistics": {
                                    "total_coverage": report.total_coverage,
                                    "mean_coverage": report.mean_coverage,
                                    "coverage_std": report.coverage_std,
                                },
                                "variants": [],
                                "summary": {
                                    "total_variants": 0,
                                    "high_coverage_variants": 0,
                                    "low_coverage_variants": 0,
                                    "pathogenic_status": "unknown",
                                },
                            }
                        gene_variants[gene_name_actual]["regions_analyzed"].add(
                            gene_name
                        )

                        # Add variant with coverage information
                        coverage = variant.get("coverage_at_variant", 0)
                        variant_data = {
                            "position": variant["position"],
                            "chromosome": variant["chromosome"],
                            "reference": variant["reference"],
                            "alternate": variant["alternate"],
                            "clinical_significance": variant["clinical_significance"],
                            "disease_name": variant["disease_name"],
                            "coverage_at_variant": coverage,
                            "coverage_status": "high" if coverage >= 10 else "low",
                            "variant_id": variant.get("variant_id", "unknown"),
                            "molecular_consequence": variant.get(
                                "molecular_consequence", "unknown"
                            ),
                            "review_status": variant.get("review_status", "unknown"),
                        }

                        # Add evidence analysis if available
                        if "evidence_analysis" in variant:
                            evidence = variant["evidence_analysis"]
                            variant_data["evidence_analysis"] = {
                                "reference_support": evidence["reference_support"],
                                "alternate_support": evidence["alternate_support"],
                                "variant_allele_frequency": evidence[
                                    "variant_allele_frequency"
                                ],
                                "evidence_level": evidence["evidence_level"],
                                "zygosity": evidence["zygosity"],
                                "mean_quality": evidence["mean_quality"],
                            }

                        gene_variants[gene_name_actual]["variants"].append(variant_data)
                    else:
                        # Fallback to region name if no gene info
                        if gene_name not in gene_variants:
                            gene_variants[gene_name] = {
                                "gene_name": gene_name,
                                "regions_analyzed": set(),
                                "coverage_statistics": {
                                    "total_coverage": report.total_coverage,
                                    "mean_coverage": report.mean_coverage,
                                    "coverage_std": report.coverage_std,
                                },
                                "variants": [],
                                "summary": {
                                    "total_variants": 0,
                                    "high_coverage_variants": 0,
                                    "low_coverage_variants": 0,
                                    "pathogenic_status": "unknown",
                                },
                            }
                        gene_variants[gene_name]["regions_analyzed"].add(gene_name)

                        # Add variant data
                        coverage = variant.get("coverage_at_variant", 0)
                        variant_data = {
                            "position": variant["position"],
                            "chromosome": variant["chromosome"],
                            "reference": variant["reference"],
                            "alternate": variant["alternate"],
                            "clinical_significance": variant["clinical_significance"],
                            "disease_name": variant["disease_name"],
                            "coverage_at_variant": coverage,
                            "coverage_status": "high" if coverage >= 10 else "low",
                            "variant_id": variant.get("variant_id", "unknown"),
                            "molecular_consequence": variant.get(
                                "molecular_consequence", "unknown"
                            ),
                            "review_status": variant.get("review_status", "unknown"),
                        }

                        if "evidence_analysis" in variant:
                            evidence = variant["evidence_analysis"]
                            variant_data["evidence_analysis"] = {
                                "reference_support": evidence["reference_support"],
                                "alternate_support": evidence["alternate_support"],
                                "variant_allele_frequency": evidence[
                                    "variant_allele_frequency"
                                ],
                                "evidence_level": evidence["evidence_level"],
                                "zygosity": evidence["zygosity"],
                                "mean_quality": evidence["mean_quality"],
                            }

                        gene_variants[gene_name]["variants"].append(variant_data)

        # Calculate summary statistics and convert sets to lists
        for gene_data in gene_variants.values():
            gene_data["regions_analyzed"] = list(gene_data["regions_analyzed"])
            gene_data["summary"]["total_variants"] = len(gene_data["variants"])
            gene_data["summary"]["high_coverage_variants"] = sum(
                1 for v in gene_data["variants"] if v["coverage_status"] == "high"
            )
            gene_data["summary"]["low_coverage_variants"] = sum(
                1 for v in gene_data["variants"] if v["coverage_status"] == "low"
            )

            if gene_data["summary"]["high_coverage_variants"] > 0:
                gene_data["summary"]["pathogenic_status"] = "present_with_good_coverage"
            elif gene_data["summary"]["total_variants"] > 0:
                gene_data["summary"]["pathogenic_status"] = "present_but_low_coverage"
            else:
                gene_data["summary"]["pathogenic_status"] = "none_detected"

        # Create the final JSON structure
        json_data = {
            "metadata": {
                "analysis_type": "lightweight_gene_analysis",
                "generated_at": pd.Timestamp.now().isoformat(),
                "total_genes_analyzed": len(gene_variants),
                "total_variants_found": sum(
                    g["summary"]["total_variants"] for g in gene_variants.values()
                ),
                "genes_with_pathogenic_variants": sum(
                    1
                    for g in gene_variants.values()
                    if g["summary"]["total_variants"] > 0
                ),
                "genes_with_good_coverage_variants": sum(
                    1
                    for g in gene_variants.values()
                    if g["summary"]["high_coverage_variants"] > 0
                ),
            },
            "genes": gene_variants,
        }

        # Write JSON file with pretty formatting
        import json

        with open(output_path, "w") as f:
            json.dump(json_data, f, indent=2, default=str)

        print(f"         JSON report generated: {output_path}")
        self.logger.info(f"Exported JSON report: {output_path}")
        return str(output_path)


def lightweight_gene_analysis_handler(job, work_dir: Optional[str] = None) -> None:
    """
    Workflow handler for lightweight gene analysis.

    This handler performs fast, targeted analysis of BED regions
    by intersecting target BAM files with known pathogenic SNPs from ClinVar.

    Expected metadata:
    - bed_path: Path to BED file defining regions of interest (optional, will look for targets_exceeding_threshold.bed)
    - clinvar_path: Path to ClinVar VCF file (optional, will use default)
    - reference: Path to reference genome (optional)
    - sample_id: Sample identifier (optional, will try to extract from job context)
    """
    print(" Starting Lightweight Gene Analysis Handler...")
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)

    try:
        # Get job parameters
        bed_path = job.context.metadata.get("bed_path")
        clinvar_path = job.context.metadata.get("clinvar_path")
        reference = job.context.metadata.get("reference")

        print(" Job parameters:")
        print(f"   - BED path: {bed_path}")
        print(f"   - ClinVar path: {clinvar_path}")
        print(f"   - Reference: {reference}")

        # Determine sample directory and sample ID
        sample_dir = work_dir or job.context.metadata.get("work_dir")
        if not sample_dir:
            sample_dir = os.path.dirname(job.context.filepath)

        # Try to get sample ID from job context or metadata
        sample_id = None
        if hasattr(job.context, "get_sample_id"):
            try:
                sample_id = job.context.get_sample_id()
                print(f" Sample ID from context: {sample_id}")
            except Exception as e:
                print(f"  Could not get sample ID from context: {e}")

        # If no sample ID from context, try to extract from metadata or filepath
        if not sample_id:
            if "sample_id" in job.context.metadata:
                sample_id = job.context.metadata["sample_id"]
                print(f" Sample ID from metadata: {sample_id}")
            else:
                # Try to extract from filepath
                try:
                    filepath = job.context.filepath
                    if filepath:
                        # Extract the last part of the path as sample ID
                        sample_id = os.path.basename(filepath)
                        print(f" Sample ID extracted from filepath: {sample_id}")
                except Exception as e:
                    print(f"  Could not extract sample ID from filepath: {e}")

        print(f" Sample directory: {sample_dir}")
        print(f" Sample ID: {sample_id}")

        # Find the best available BAM file for analysis
        target_bam = None

        # ALWAYS prioritize igv_ready.bam from IGV folders first (these are properly indexed)
        print(" Searching for igv_ready.bam in IGV folders (priority 1)...")

        # Priority 1: Try igv_ready.bam from IGV folder in sample subdirectory (if sample_id is known)
        if sample_id:
            igv_bam = os.path.join(sample_dir, sample_id, "igv", "igv_ready.bam")
            if os.path.exists(igv_bam):
                target_bam = igv_bam
                print(f" Found igv_ready.bam in sample IGV folder: {target_bam}")
                print(
                    " This file should be properly indexed and ready for pileup analysis"
                )
            else:
                print(f" igv_ready.bam not found at: {igv_bam}")

        # Priority 2: Try igv_ready.bam from IGV folder in main directory
        if not target_bam:
            igv_bam = os.path.join(sample_dir, "igv", "igv_ready.bam")
            if os.path.exists(igv_bam):
                target_bam = igv_bam
                print(f" Found igv_ready.bam in main IGV folder: {target_bam}")
                print(
                    " This file should be properly indexed and ready for pileup analysis"
                )
            else:
                print(f" igv_ready.bam not found at: {igv_bam}")

        # Priority 3: Search ALL subdirectories for igv_ready.bam (IGV files are priority)
        if not target_bam:
            print(
                " No IGV igv_ready.bam found, searching all subdirectories for igv_ready.bam..."
            )
            for root, dirs, files in os.walk(sample_dir):
                if "igv_ready.bam" in files:
                    target_bam = os.path.join(root, "igv_ready.bam")
                    print(f" Found igv_ready.bam in subdirectory: {target_bam}")
                    print(
                        " This file should be properly indexed and ready for pileup analysis"
                    )
                    break

        # NO FALLBACK to target.bam - only use properly indexed IGV files
        if not target_bam:
            print(" No igv_ready.bam found in any IGV folders or subdirectories")
            print(" Searched locations:")
            if sample_id:
                print(
                    f"   - {os.path.join(sample_dir, sample_id, 'igv', 'igv_ready.bam')}"
                )
            print(f"   - {os.path.join(sample_dir, 'igv', 'igv_ready.bam')}")
            print("   - All subdirectories recursively")
            raise FileNotFoundError(
                f"No igv_ready.bam found in {sample_dir} or subdirectories. This file is required for analysis."
            )

        if not target_bam:
            raise FileNotFoundError(
                f"No suitable BAM file found in {sample_dir} or subdirectories"
            )

        # Check if BAM file has an index
        bam_index = target_bam + ".bai"
        if not os.path.exists(bam_index):
            print(f"  BAM file {target_bam} has no index (.bai file)")
            print(" Attempting to create index using samtools...")

            try:
                import subprocess

                result = subprocess.run(
                    ["samtools", "index", target_bam],
                    capture_output=True,
                    text=True,
                    timeout=300,
                )
                if result.returncode == 0:
                    print(" Successfully created BAM index")
                else:
                    print(f" Failed to create BAM index: {result.stderr}")
                    print("  Analysis may fail due to missing index")
            except Exception as e:
                print(f" Error creating BAM index: {e}")
                print("  Analysis may fail due to missing index")

        # Find BED file if not specified
        if not bed_path:
            # Try to find BED file in sample subdirectory first (if sample_id is known)
            if sample_id:
                bed_path = os.path.join(
                    sample_dir, sample_id, "targets_exceeding_threshold.bed"
                )
                if os.path.exists(bed_path):
                    print(f" Found BED file in sample subdirectory: {bed_path}")
                else:
                    print(f" BED file not found at: {bed_path}")
                    bed_path = None
            else:
                bed_path = os.path.join(sample_dir, "targets_exceeding_threshold.bed")

        if not bed_path or not os.path.exists(bed_path):
            # Try to find targets_exceeding_threshold.bed in subdirectories
            print(
                " targets_exceeding_threshold.bed not found, searching subdirectories..."
            )
            for root, dirs, files in os.walk(sample_dir):
                if "targets_exceeding_threshold.bed" in files:
                    bed_path = os.path.join(root, "targets_exceeding_threshold.bed")
                    print(f" Found BED file in subdirectory: {bed_path}")
                    break
            else:
                raise FileNotFoundError(
                    f"targets_exceeding_threshold.bed not found in {sample_dir} or any subdirectories"
                )
        else:
            print(f" Found BED file: {bed_path}")

        # Use default ClinVar path if not specified
        if not clinvar_path:
            clinvar_path = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                "..",
                "resources",
                "clinvar.vcf.gz",
            )

        print(f" ClinVar path: {clinvar_path}")

        logger.info("Starting lightweight gene analysis for BED regions")
        logger.info(f"Target BAM: {target_bam}")
        logger.info(f"BED file: {bed_path}")
        logger.info(f"ClinVar path: {clinvar_path}")

        # Initialize analysis
        print(" Initializing LightweightGeneAnalysis...")
        analysis = LightweightGeneAnalysis(
            work_dir=sample_dir, clinvar_path=clinvar_path, reference=reference
        )

        # Perform analysis using BED regions
        print(" Starting BED region analysis...")

        results = analysis.analyze_bed_regions(target_bam, bed_path)

        if not results:
            print(
                "  No regions were successfully analyzed - this may indicate no ClinVar variants intersect with your BED regions"
            )
            logger.warning(
                "No regions were successfully analyzed - this may indicate no ClinVar variants intersect with your BED regions"
            )
            # Create empty results structure
            results = {}

        # Generate JSON report
        print(" Generating JSON report...")

        # Determine the correct output directory for the sample
        if sample_id:
            output_dir = os.path.join(sample_dir, sample_id)
        else:
            output_dir = sample_dir

        # Ensure the output directory exists
        os.makedirs(output_dir, exist_ok=True)

        json_report = analysis.export_json(
            results, os.path.join(output_dir, "lightweight_gene_analysis_results.json")
        )

        print(" JSON report generated:")
        print(f"   - JSON: {json_report}")

        # Record results
        try:
            job.context.add_result(
                "lightweight_gene_analysis",
                {
                    "region_reports": len(results),
                    "json_report": json_report,
                    "regions_analyzed": list(results.keys()),
                    "total_pathogenic_snps": sum(
                        r.pathogenic_snps_found for r in results.values()
                    ),
                    "analysis_timestamp": pd.Timestamp.now().isoformat(),
                },
            )
        except Exception:
            pass

        print(" Lightweight gene analysis completed successfully!")
        print(
            f" Analyzed {len(results)} regions, found {sum(r.pathogenic_snps_found for r in results.values())} pathogenic SNPs"
        )
        logger.info("Lightweight gene analysis completed successfully")
        logger.info(
            f"Analyzed {len(results)} regions, found {sum(r.pathogenic_snps_found for r in results.values())} pathogenic SNPs"
        )

    except Exception as e:
        print(f" Lightweight gene analysis handler failed: {e}")
        job.context.add_error("lightweight_gene_analysis", str(e))
        logger.error(f"Lightweight gene analysis handler failed: {e}")
        raise
