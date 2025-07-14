"""
CNV Analysis Section for ROBIN Reports.

This module handles the Copy Number Variation (CNV) analysis section of the report.
"""

import os
import re
import pickle
import logging
import numpy as np
import pandas as pd
import natsort
from reportlab.lib.units import inch
from reportlab.platypus import PageBreak, Paragraph, Image, Spacer, Table, TableStyle
from reportlab.lib.styles import ParagraphStyle
from ..sections.base import ReportSection
from ..plotting import create_CNV_plot, create_CNV_plot_per_chromosome
from robin.subpages.CNVObjectClass import (
    Result,
    CNVAnalysis,
    CNV_Difference,
    moving_average,
    # iterate_bam_bin,
)
from robin import resources

logger = logging.getLogger(__name__)


class CNVSection(ReportSection):
    """Section containing the CNV analysis."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.current_row = []
        self.plots_per_row = 2

    def add_content(self):
        """Add the CNV analysis content to the report."""
        logger.debug("Starting CNV section processing")

        # Load CNV data and XYestimate
        XYestimate = "Unknown"  # Default value
        cnv_file = os.path.join(self.report.output, "CNV.npy")

        # Check for required files
        if not os.path.exists(cnv_file):
            logger.error("No CNV.npy file found in output directory")
            return

        # Load CNV data
        logger.debug("Loading CNV data from %s", cnv_file)
        CNVresult = np.load(cnv_file, allow_pickle="TRUE").item()
        CNVresult = Result(CNVresult)
        logger.debug("CNV data loaded with keys: %s", list(CNVresult.cnv.keys())[:5])

        cnv_dict = np.load(
            os.path.join(self.report.output, "CNV_dict.npy"), allow_pickle=True
        ).item()
        logger.debug("CNV dict loaded with keys: %s", list(cnv_dict.keys()))

        # Store cnv_dict in report for use by other methods
        self.report.cnv_dict = cnv_dict

        # Load XY estimate if available
        if os.path.exists(os.path.join(self.report.output, "XYestimate.pkl")):
            with open(os.path.join(self.report.output, "XYestimate.pkl"), "rb") as file:
                XYestimate = pickle.load(file)
                logger.debug("Loaded XY estimate: %s", XYestimate)

        # Add CNV section header
        logger.debug("Adding CNV section header")

        # Start detailed analysis section
        self.elements.append(PageBreak())
        self.elements.append(
            Paragraph(
                "Copy Number Variation Detailed Analysis",
                self.styles.styles["Heading2"],
            )
        )

        try:
            # Initialize CNVAnalysis object with the same settings as UI
            cnv_analyzer = CNVAnalysis(target_panel="rCNS2")
            cnv_analyzer.XYestimate = XYestimate

            # Load required resource files
            gene_bed_file = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)), "unique_genes.bed"
            )
            cytoband_file = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)), "cytoBand.txt"
            )
            logger.debug(
                "Resource files: gene_bed=%s, cytoband=%s", gene_bed_file, cytoband_file
            )

            # Load gene and cytoband data
            gene_bed = None
            cytobands_bed = None
            if os.path.exists(gene_bed_file):
                gene_bed = pd.read_csv(
                    gene_bed_file,
                    sep="\t",
                    names=["chrom", "start_pos", "end_pos", "gene"],
                )
                logger.debug("Loaded gene bed file with shape: %s", gene_bed.shape)
            if os.path.exists(cytoband_file):
                cytobands_bed = pd.read_csv(
                    cytoband_file,
                    sep="\t",
                    names=["chrom", "start_pos", "end_pos", "name", "stain"],
                )
                logger.debug("Loaded cytoband file with shape: %s", cytobands_bed.shape)

            # Set up CNVAnalysis object with loaded data
            cnv_analyzer.gene_bed = gene_bed
            cnv_analyzer.cytobands_bed = cytobands_bed
            cnv_analyzer.cnv_dict = cnv_dict

            # Get reference CNV data with matching bin width
            logger.debug(
                "Getting reference CNV data with bin width %s", cnv_dict["bin_width"]
            )

            r2_cnv = Result(
                np.load(
                    os.path.join(self.report.output, "CNV2.npy"), allow_pickle="TRUE"
                ).item()
            ).cnv

            # Initialize CNV_Difference object for normalized values
            result3 = CNV_Difference()

            # Calculate normalized CNV values
            logger.debug("Calculating normalized CNV values")
            for key in CNVresult.cnv.keys():
                if key != "chrM" and re.match(r"^chr(\d+|X|Y)$", key):
                    if key in r2_cnv:
                        moving_avg_data1 = moving_average(CNVresult.cnv[key])
                        moving_avg_data2 = moving_average(r2_cnv[key])
                        # Pad arrays if needed
                        if len(moving_avg_data1) != len(moving_avg_data2):
                            max_len = max(len(moving_avg_data1), len(moving_avg_data2))
                            if len(moving_avg_data1) < max_len:
                                moving_avg_data1 = np.pad(
                                    moving_avg_data1,
                                    (0, max_len - len(moving_avg_data1)),
                                )
                            if len(moving_avg_data2) < max_len:
                                moving_avg_data2 = np.pad(
                                    moving_avg_data2,
                                    (0, max_len - len(moving_avg_data2)),
                                )
                        # Calculate difference
                        result3.cnv[key] = moving_avg_data1 - moving_avg_data2

            # Set the result3 in the analyzer
            cnv_analyzer.result3 = result3

            # Calculate chromosome statistics using CNVAnalysis logic
            chromosome_stats = cnv_analyzer.calculate_chromosome_stats(
                CNVresult, r2_cnv
            )
            cnv_analyzer.chromosome_stats = chromosome_stats

            # Add gain/loss thresholds to chromosome stats
            for chrom, stats in chromosome_stats.items():
                if chrom != "global":
                    # Set default thresholds to 0.4/-0.4 for autosomes
                    stats["gain_threshold"] = 0.4  # Changed from 0.25
                    stats["loss_threshold"] = -0.4  # Changed from -0.25
                    if chrom == "chrX":
                        if XYestimate == "XY":  # Male
                            stats["gain_threshold"] = 0.3
                            stats["loss_threshold"] = -0.3
                        else:  # Female
                            stats["gain_threshold"] = 0.75
                            stats["loss_threshold"] = -0.75
                    elif chrom == "chrY":
                        if XYestimate == "XY":  # Male
                            stats["gain_threshold"] = 0.5
                            stats["loss_threshold"] = -0.5
                        else:  # Female
                            stats["gain_threshold"] = -0.2
                            stats["loss_threshold"] = -1.0

            # Add Summary Card
            logger.debug("Adding CNV summary card")
            # Create summary card table data
            summary_data = []

            # Add genetic sex row (simplified)
            summary_data.append(
                [
                    Paragraph("Genetic Sex:", self.styles.styles["Normal"]),
                    Paragraph(XYestimate, self.styles.styles["Normal"]),
                ]
            )

            # Add analysis metrics
            summary_data.append(
                [
                    Paragraph("Bin Width:", self.styles.styles["Normal"]),
                    Paragraph(
                        f"{cnv_dict['bin_width']:,}", self.styles.styles["Normal"]
                    ),
                ]
            )
            summary_data.append(
                [
                    Paragraph("Variance:", self.styles.styles["Normal"]),
                    Paragraph(
                        f"{cnv_dict.get('variance', 0):.3f}",
                        self.styles.styles["Normal"],
                    ),
                ]
            )

            # Calculate gene counts
            total_gained_genes = set()
            total_lost_genes = set()
            for chrom in natsort.natsorted(result3.cnv.keys()):
                if chrom != "chrM" and re.match(r"^chr(\d+|X|Y)$", chrom):
                    analysis = cnv_analyzer.analyze_cytoband_cnv(result3.cnv, chrom)
                    if not analysis.empty:
                        # Get genes in gained regions (including HIGH_GAIN)
                        gained = analysis[
                            analysis["cnv_state"].isin(["GAIN", "HIGH_GAIN"])
                        ]
                        for _, row in gained.iterrows():
                            if row["genes"]:
                                total_gained_genes.update(row["genes"])

                        # Get genes in lost regions (including DEEP_LOSS)
                        lost = analysis[
                            analysis["cnv_state"].isin(["LOSS", "DEEP_LOSS"])
                        ]
                        for _, row in lost.iterrows():
                            if row["genes"]:
                                total_lost_genes.update(row["genes"])

            # Add gene counts to summary
            summary_data.append(
                [
                    Paragraph("Genes in Gained Regions:", self.styles.styles["Normal"]),
                    Paragraph(
                        str(len(total_gained_genes)), self.styles.styles["Normal"]
                    ),
                ]
            )
            summary_data.append(
                [
                    Paragraph("Genes in Lost Regions:", self.styles.styles["Normal"]),
                    Paragraph(str(len(total_lost_genes)), self.styles.styles["Normal"]),
                ]
            )

            # Create summary table with styling
            if summary_data:
                formatted_summary_data = []
                for row in summary_data:
                    formatted_row = [
                        row[0].text if hasattr(row[0], "text") else str(row[0]),
                        row[1].text if hasattr(row[1], "text") else str(row[1]),
                    ]
                    formatted_summary_data.append(formatted_row)

                summary_table = self.create_table(
                    formatted_summary_data, auto_col_width=True
                )
                # Add specific styling while preserving modern table style
                summary_table.setStyle(
                    TableStyle(
                        [
                            *self.MODERN_TABLE_STYLE._cmds,
                            (
                                "ALIGN",
                                (1, 0),
                                (1, -1),
                                "RIGHT",
                            ),  # Right-align the count column
                            (
                                "FONTNAME",
                                (0, 0),
                                (-1, -1),
                                "Helvetica-Bold",
                            ),  # Bold font for all cells
                            (
                                "FONTSIZE",
                                (0, 0),
                                (-1, -1),
                                10,
                            ),  # Consistent font size
                            (
                                "LEADING",
                                (0, 0),
                                (-1, -1),
                                12,
                            ),  # Line spacing
                            (
                                "TOPPADDING",
                                (0, 0),
                                (-1, -1),
                                6,
                            ),  # Top padding
                            (
                                "BOTTOMPADDING",
                                (0, 0),
                                (-1, -1),
                                6,
                            ),  # Bottom padding
                        ]
                    )
                )
                self.elements.append(summary_table)

            # First check for arm events and whole chromosome events
            summary_whole_chr_events = []
            summary_arm_events = []

            # Analyze each chromosome
            for chrom in natsort.natsorted(CNVresult.cnv.keys()):
                if chrom != "chrM" and re.match(r"^chr(\d+|X|Y)$", chrom):
                    logger.info(f"Analyzing chromosome {chrom} for events")
                    cytoband_analysis = cnv_analyzer.analyze_cytoband_cnv(
                        result3.cnv, chrom
                    )
                    if not cytoband_analysis.empty:
                        # Skip Y chromosome for male samples
                        if chrom == "chrY" and XYestimate == "XY":
                            continue

                        gain_threshold = chromosome_stats[chrom]["gain_threshold"]
                        loss_threshold = chromosome_stats[chrom]["loss_threshold"]
                        logger.info(
                            f"{chrom} thresholds: gain={gain_threshold:.2f}, loss={loss_threshold:.2f}"
                        )

                        # Log the cytoband names for debugging
                        logger.info(
                            f"Available cytoband names: {cytoband_analysis['name'].tolist()}"
                        )

                        # Check both arms
                        p_arm_mean = None
                        q_arm_mean = None
                        p_arm_proportion = 0
                        q_arm_proportion = 0

                        # Analyze p arm
                        p_arm_cytobands = cytoband_analysis[
                            cytoband_analysis["name"].str.contains(
                                f"{chrom} p", regex=False
                            )
                        ]
                        if not p_arm_cytobands.empty:
                            p_arm_values = p_arm_cytobands["mean_cnv"].values
                            p_arm_mean = np.mean(p_arm_values)
                            logger.info(f"{chrom} p-arm mean: {p_arm_mean:.2f}")
                            # Calculate proportion affected regardless of mean
                            p_arm_proportion_gain = np.sum(
                                p_arm_values > gain_threshold
                            ) / len(p_arm_values)
                            p_arm_proportion_loss = np.sum(
                                p_arm_values < loss_threshold
                            ) / len(p_arm_values)
                            p_arm_proportion = max(
                                p_arm_proportion_gain, p_arm_proportion_loss
                            )
                            logger.info(
                                f"{chrom} p-arm proportion affected: {p_arm_proportion:.2f}"
                            )

                        # Analyze q arm
                        q_arm_cytobands = cytoband_analysis[
                            cytoband_analysis["name"].str.contains(
                                f"{chrom} q", regex=False
                            )
                        ]
                        if not q_arm_cytobands.empty:
                            q_arm_values = q_arm_cytobands["mean_cnv"].values
                            q_arm_mean = np.mean(q_arm_values)
                            logger.info(f"{chrom} q-arm mean: {q_arm_mean:.2f}")
                            # Calculate proportion affected regardless of mean
                            q_arm_proportion_gain = np.sum(
                                q_arm_values > gain_threshold
                            ) / len(q_arm_values)
                            q_arm_proportion_loss = np.sum(
                                q_arm_values < loss_threshold
                            ) / len(q_arm_values)
                            q_arm_proportion = max(
                                q_arm_proportion_gain, q_arm_proportion_loss
                            )
                            logger.info(
                                f"{chrom} q-arm proportion affected: {q_arm_proportion:.2f}"
                            )

                        # Determine if this is a whole chromosome event or arm-specific events
                        whole_chr_mean = chromosome_stats[chrom]["mean"]

                        # Check if both arms show similar changes (whole chromosome event)
                        if p_arm_mean is not None and q_arm_mean is not None:
                            # For gains, check if both arms show consistent gains
                            both_arms_gained = (
                                # Both arms must show gain above threshold
                                p_arm_mean > gain_threshold
                                and q_arm_mean > gain_threshold
                                and
                                # At least one arm must have high proportion affected
                                (p_arm_proportion > 0.7 or q_arm_proportion > 0.7)
                                and
                                # Other arm must show some effect
                                (p_arm_proportion > 0.4 and q_arm_proportion > 0.4)
                            )

                            # For losses, check if both arms show consistent losses
                            both_arms_lost = (
                                # Both arms must show loss below threshold
                                p_arm_mean < loss_threshold
                                and q_arm_mean < loss_threshold
                                and
                                # At least one arm must have high proportion affected
                                (p_arm_proportion > 0.7 or q_arm_proportion > 0.7)
                                and
                                # Other arm must show some effect
                                (p_arm_proportion > 0.4 and q_arm_proportion > 0.4)
                            )

                            # Check for whole chromosome event first
                            if both_arms_gained or both_arms_lost:
                                # Report as whole chromosome event
                                if whole_chr_mean > gain_threshold:
                                    summary_whole_chr_events.append(
                                        f"Chromosome {chrom[3:]}: GAIN (mean={whole_chr_mean:.2f})"
                                    )
                                    logger.info(
                                        f"Added whole chromosome gain for {chrom}"
                                    )
                                elif whole_chr_mean < loss_threshold:
                                    summary_whole_chr_events.append(
                                        f"Chromosome {chrom[3:]}: LOSS (mean={whole_chr_mean:.2f})"
                                    )
                                    logger.info(
                                        f"Added whole chromosome loss for {chrom}"
                                    )
                            else:
                                # If not a whole chromosome event, check for individual arm events
                                # Skip Y chromosome for arm events
                                if chrom != "chrY":
                                    if (
                                        p_arm_mean is not None
                                        and p_arm_proportion > 0.7
                                    ):
                                        if p_arm_mean > gain_threshold:
                                            summary_arm_events.append(
                                                f"Chromosome {chrom[3:]} p-arm: GAIN (mean={p_arm_mean:.2f}, {p_arm_proportion:.0%} of arm)"
                                            )
                                            logger.info(
                                                f"Added arm event: {chrom} p-arm GAIN"
                                            )
                                        elif p_arm_mean < loss_threshold:
                                            summary_arm_events.append(
                                                f"Chromosome {chrom[3:]} p-arm: LOSS (mean={p_arm_mean:.2f}, {p_arm_proportion:.0%} of arm)"
                                            )
                                            logger.info(
                                                f"Added arm event: {chrom} p-arm LOSS"
                                            )

                                    if (
                                        q_arm_mean is not None
                                        and q_arm_proportion > 0.7
                                    ):
                                        if q_arm_mean > gain_threshold:
                                            summary_arm_events.append(
                                                f"Chromosome {chrom[3:]} q-arm: GAIN (mean={q_arm_mean:.2f}, {q_arm_proportion:.0%} of arm)"
                                            )
                                            logger.info(
                                                f"Added arm event: {chrom} q-arm GAIN"
                                            )
                                        elif q_arm_mean < loss_threshold:
                                            summary_arm_events.append(
                                                f"Chromosome {chrom[3:]} q-arm: LOSS (mean={q_arm_mean:.2f}, {q_arm_proportion:.0%} of arm)"
                                            )
                                            logger.info(
                                                f"Added arm event: {chrom} q-arm LOSS"
                                            )
                        else:
                            # Only one arm has data, check for whole chromosome event
                            # For chromosomes with only one arm, use a stricter threshold
                            if abs(whole_chr_mean) > abs(gain_threshold) * 1.5:
                                if whole_chr_mean > gain_threshold:
                                    summary_whole_chr_events.append(
                                        f"Chromosome {chrom[3:]}: GAIN (mean={whole_chr_mean:.2f})"
                                    )
                                    logger.info(
                                        f"Added whole chromosome gain for {chrom}"
                                    )
                                elif whole_chr_mean < loss_threshold:
                                    summary_whole_chr_events.append(
                                        f"Chromosome {chrom[3:]}: LOSS (mean={whole_chr_mean:.2f})"
                                    )
                                    logger.info(
                                        f"Added whole chromosome loss for {chrom}"
                                    )

            # Log the final counts
            logger.info(f"Found {len(summary_arm_events)} arm events")
            logger.info(
                f"Found {len(summary_whole_chr_events)} whole chromosome events"
            )

            self.summary_elements.append(
                Paragraph(
                    "Copy Number Variation Summary",
                    ParagraphStyle(
                        "SummaryHeader",
                        parent=self.styles.styles["Heading3"],
                        fontSize=12,
                        fontName="Helvetica-Bold",
                        textColor=self.styles.COLORS["primary"],
                        spaceAfter=12,
                    ),
                )
            )

            # Add whole chromosome events to summary
            if summary_whole_chr_events:
                self.summary_elements.append(
                    Paragraph(
                        "Whole Chromosome Events:<br/> "
                        + " <br/> ".join(summary_whole_chr_events),
                        ParagraphStyle(
                            "SummaryText",
                            parent=self.styles.styles["Normal"],
                            fontSize=10,
                            fontName="Helvetica",
                            textColor=self.styles.COLORS["text"],
                            leading=14,
                            spaceAfter=12,
                        ),
                    )
                )

            # Add arm events to summary
            if summary_arm_events:
                logger.debug(f"Found {len(summary_arm_events)} arm events to report")
                self.summary_elements.append(
                    Paragraph(
                        "Chromosome Arm Events:<br/> "
                        + " <br/> ".join(summary_arm_events),
                        ParagraphStyle(
                            "SummaryText",
                            parent=self.styles.styles["Normal"],
                            fontSize=10,
                            fontName="Helvetica",
                            textColor=self.styles.COLORS["text"],
                            leading=14,
                            spaceAfter=12,
                        ),
                    )
                )

            # Generate genome-wide CNV plot
            logger.debug("Generating genome-wide CNV plot")
            img_buf = create_CNV_plot(CNVresult, cnv_dict)
            width, height = inch * 7.5, inch * 2  # A4 width minus margins
            self.summary_elements.append(Image(img_buf, width=width, height=height))
            self.summary_elements.append(
                Paragraph(
                    "Copy number variation across chromosomes",
                    ParagraphStyle(
                        "PlotCaption",
                        parent=self.styles.styles["Caption"],
                        fontSize=9,
                        fontName="Helvetica",
                        textColor=self.styles.COLORS["text"],
                        alignment=1,  # Center alignment
                        spaceBefore=6,
                        spaceAfter=12,
                    ),
                )
            )

            # Create summary of whole chromosome events and gene-containing events
            logger.debug("Creating CNV summary")

            # Add whole chromosome and arm events summary
            whole_chr_events = []
            arm_events = []
            gene_containing_events = []

            for chrom in natsort.natsorted(CNVresult.cnv.keys()):
                if chrom != "chrM" and re.match(r"^chr(\d+|X|Y)$", chrom):
                    cytoband_analysis = cnv_analyzer.analyze_cytoband_cnv(
                        result3.cnv, chrom
                    )
                    if not cytoband_analysis.empty:
                        # Find whole chromosome events
                        whole_chr = cytoband_analysis[
                            cytoband_analysis["name"].str.contains(
                                "WHOLE CHROMOSOME", na=False
                            )
                        ]
                        if not whole_chr.empty:
                            for _, row in whole_chr.iterrows():
                                whole_chr_events.append(
                                    [
                                        row["chrom"].replace("chr", ""),
                                        row["cnv_state"],
                                        f"{row['mean_cnv']:.3f}",
                                    ]
                                )

                        # Find arm-specific events, but only if there isn't a whole chromosome event
                        if whole_chr.empty:
                            for arm in ["p", "q"]:
                                arm_cytobands = cytoband_analysis[
                                    cytoband_analysis["name"].str.startswith(arm)
                                ]
                                if not arm_cytobands.empty:
                                    arm_mean = arm_cytobands["mean_cnv"].mean()
                                    arm_state = (
                                        "GAIN"
                                        if arm_mean > 0.5
                                        else "LOSS" if arm_mean < -0.5 else "NORMAL"
                                    )
                                    if arm_state != "NORMAL":
                                        # Calculate proportion of arm showing the event
                                        arm_bins_above = np.sum(
                                            arm_cytobands["mean_cnv"] > 0.5
                                        ) / len(arm_cytobands)
                                        arm_bins_below = np.sum(
                                            arm_cytobands["mean_cnv"] < -0.5
                                        ) / len(arm_cytobands)
                                        arm_proportion = max(
                                            arm_bins_above, arm_bins_below
                                        )

                                        if (
                                            arm_proportion > 0.7
                                        ):  # Require at least 70% of arm to show the event
                                            arm_events.append(
                                                [
                                                    row["chrom"].replace("chr", ""),
                                                    f"{arm}-arm",
                                                    arm_state,
                                                    f"{arm_mean:.3f}",
                                                    f"{arm_proportion:.1%}",  # Add proportion of arm affected
                                                ]
                                            )

                        # Find events with genes
                        events_with_genes = cytoband_analysis[
                            (cytoband_analysis["genes"].apply(len) > 0)
                            & (
                                ~cytoband_analysis["name"].str.contains(
                                    "WHOLE CHROMOSOME", na=False
                                )
                            )
                        ]
                        for _, row in events_with_genes.iterrows():
                            gene_containing_events.append(
                                [
                                    row["chrom"].replace("chr", ""),
                                    row["name"].replace(f"{row['chrom']} ", ""),
                                    row["cnv_state"],
                                    f"{row['mean_cnv']:.3f}",
                                    ", ".join(row["genes"]),
                                ]
                            )

            # Add whole chromosome events summary if any exist
            if whole_chr_events:
                self.elements.append(
                    Paragraph("Whole Chromosome Events", self.styles.styles["Heading4"])
                )
                whole_chr_data = [["Chr", "State", "Mean CNV"]]
                whole_chr_data.extend(whole_chr_events)
                whole_chr_table = self.create_table(
                    whole_chr_data,
                    repeat_rows=1,
                    auto_col_width=False,
                    col_widths=[inch * x for x in [0.4, 0.8, 0.8]],
                )
                whole_chr_table.setStyle(
                    TableStyle(
                        [
                            *self.MODERN_TABLE_STYLE._cmds,
                            ("ALIGN", (2, 1), (2, -1), "RIGHT"),  # Right-align mean CNV
                            ("ALIGN", (1, 1), (1, -1), "CENTER"),  # Center-align state
                        ]
                    )
                )
                self.elements.append(whole_chr_table)

            # Add chromosome arm events summary if any exist
            if arm_events:
                self.elements.append(
                    Paragraph("Chromosome Arm Events", self.styles.styles["Heading4"])
                )
                arm_data = [["Chr", "Arm", "State", "Mean CNV", "Proportion Affected"]]
                arm_data.extend(arm_events)
                arm_table = self.create_table(
                    arm_data,
                    repeat_rows=1,
                    auto_col_width=False,
                    col_widths=[inch * x for x in [0.4, 0.4, 0.8, 0.8, 1.0]],
                )
                arm_table.setStyle(
                    TableStyle(
                        [
                            *self.MODERN_TABLE_STYLE._cmds,
                            ("ALIGN", (3, 1), (3, -1), "RIGHT"),  # Right-align mean CNV
                            ("ALIGN", (2, 1), (2, -1), "CENTER"),  # Center-align state
                            (
                                "ALIGN",
                                (4, 1),
                                (4, -1),
                                "RIGHT",
                            ),  # Right-align proportion
                        ]
                    )
                )
                self.elements.append(arm_table)

            # Add gene-containing events if any exist
            if gene_containing_events:
                self.elements.append(
                    Paragraph(
                        "CNV Events Containing Genes", self.styles.styles["Heading4"]
                    )
                )
                gene_events_data = [["Chr", "Region", "State", "Mean CNV", "Genes"]]
                gene_events_data.extend(gene_containing_events)
                gene_events_table = self.create_table(
                    gene_events_data,
                    repeat_rows=1,
                    auto_col_width=False,
                    col_widths=[inch * x for x in [0.4, 1.0, 0.6, 0.6, 4.0]],
                )
                gene_events_table.setStyle(
                    TableStyle(
                        [
                            *self.MODERN_TABLE_STYLE._cmds,
                            ("ALIGN", (3, 1), (3, -1), "RIGHT"),  # Right-align mean CNV
                            ("ALIGN", (2, 1), (2, -1), "CENTER"),  # Center-align state
                        ]
                    )
                )
                self.elements.append(gene_events_table)

            # Add note about detailed view
            self.elements.append(
                Paragraph(
                    "Note: Full CNV details are available in the detailed view.",
                    ParagraphStyle(
                        "Note",
                        parent=self.styles.styles["Normal"],
                        fontSize=8,
                        textColor=self.styles.COLORS["text"],
                        spaceBefore=6,
                        spaceAfter=6,
                        italics=True,
                    ),
                )
            )

            # Initialize variables for plot layout
            self.current_row = []

            try:
                # Generate all chromosome plots at once
                logger.debug("Generating individual chromosome plots")
                chromosome_plots = create_CNV_plot_per_chromosome(CNVresult, cnv_dict)

                # Filter plots to only show chromosomes with significant changes
                significant_chromosomes = set()
                for chrom in natsort.natsorted(CNVresult.cnv.keys()):
                    if chrom != "chrM" and re.match(r"^chr(\d+|X|Y)$", chrom):
                        cytoband_analysis = cnv_analyzer.analyze_cytoband_cnv(
                            result3.cnv, chrom
                        )
                        if not cytoband_analysis.empty and any(
                            row["cnv_state"]
                            in ["GAIN", "LOSS", "HIGH_GAIN", "DEEP_LOSS"]
                            for _, row in cytoband_analysis.iterrows()
                        ):
                            significant_chromosomes.add(chrom)

                filtered_plots = [
                    (chrom, plot_buf)
                    for chrom, plot_buf in chromosome_plots
                    if chrom in significant_chromosomes
                ]

                for chrom, img_buf in filtered_plots:
                    # Add plot and its caption
                    plot_elements = [
                        Image(img_buf, width=inch * 3.5, height=inch * 1.5),
                        Paragraph(
                            f"Chromosome {chrom.replace('chr', '')}",
                            self.styles.styles["Caption"],
                        ),
                    ]
                    self.current_row.append(
                        Table(
                            [[plot_elements[0]], [plot_elements[1]]],
                            style=[("ALIGN", (0, 0), (-1, -1), "CENTER")],
                        )
                    )

                    # When row is full or it's the last plot, add the row to elements
                    if (
                        len(self.current_row) == self.plots_per_row
                        or (chrom, img_buf) == filtered_plots[-1]
                    ):
                        # If it's the last row and not full, add empty space
                        while len(self.current_row) < self.plots_per_row:
                            self.current_row.append(Spacer(inch * 3.5, inch * 1.8))

                        # Create row table and add to elements
                        plot_row = Table(
                            [self.current_row],
                            colWidths=[inch * 3.5] * self.plots_per_row,
                        )
                        plot_row.setStyle(
                            TableStyle(
                                [
                                    ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                                    ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                                ]
                            )
                        )
                        self.elements.append(plot_row)
                        self.current_row = []

                # Add detailed CNV table
                self.elements.append(Spacer(1, 12))
                self.elements.append(
                    Paragraph("Detailed CNV Events", self.styles.styles["Heading3"])
                )

                # Create detailed table
                all_cnv_events = []
                for chrom in natsort.natsorted(CNVresult.cnv.keys()):
                    if chrom != "chrM" and re.match(r"^chr(\d+|X|Y)$", chrom):
                        cytoband_analysis = cnv_analyzer.analyze_cytoband_cnv(
                            result3.cnv, chrom
                        )
                        if not cytoband_analysis.empty:
                            for _, row in cytoband_analysis.iterrows():
                                if row["cnv_state"] in [
                                    "GAIN",
                                    "LOSS",
                                    "HIGH_GAIN",
                                    "DEEP_LOSS",
                                ]:
                                    all_cnv_events.append(
                                        [
                                            row["chrom"].replace("chr", ""),
                                            row["name"].replace(f"{row['chrom']} ", ""),
                                            f"{row['start_pos']/1e6:.2f}",
                                            f"{row['end_pos']/1e6:.2f}",
                                            f"{row['length']/1e6:.2f}",
                                            f"{float(row['mean_cnv']):.3f}",
                                            row["cnv_state"],
                                            (
                                                ", ".join(row["genes"])
                                                if isinstance(row["genes"], list)
                                                else ""
                                            ),
                                        ]
                                    )

                if all_cnv_events:
                    # Convert all data to Paragraphs with proper styling
                    formatted_events = []
                    for event in all_cnv_events:
                        formatted_events.append(
                            [
                                Paragraph(
                                    event[0], self.styles.styles["Normal"]
                                ),  # Chr
                                Paragraph(
                                    event[1], self.styles.styles["Normal"]
                                ),  # Region
                                Paragraph(
                                    event[2], self.styles.styles["Normal"]
                                ),  # Start
                                Paragraph(
                                    event[3], self.styles.styles["Normal"]
                                ),  # End
                                Paragraph(
                                    event[4], self.styles.styles["Normal"]
                                ),  # Length
                                Paragraph(
                                    event[5], self.styles.styles["Normal"]
                                ),  # Mean CNV
                                Paragraph(
                                    event[6], self.styles.styles["Normal"]
                                ),  # State
                                Paragraph(
                                    event[7],
                                    ParagraphStyle(
                                        "GeneList",
                                        parent=self.styles.styles["Normal"],
                                        leading=10,  # Adjust line spacing
                                        spaceBefore=1,
                                        spaceAfter=1,
                                        wordWrap="LTR",  # Left to right word wrap
                                    ),
                                ),  # Genes
                            ]
                        )

                    # Format detailed CNV table data
                    detailed_data = [
                        [
                            "Chr",
                            "Region",
                            "Start (Mb)",
                            "End (Mb)",
                            "Length (Mb)",
                            "Mean CNV",
                            "State",
                            "Genes",
                        ]
                    ]

                    for row in formatted_events:
                        detailed_data.append(
                            [
                                row[0].text if hasattr(row[0], "text") else str(row[0]),
                                row[1].text if hasattr(row[1], "text") else str(row[1]),
                                row[2].text if hasattr(row[2], "text") else str(row[2]),
                                row[3].text if hasattr(row[3], "text") else str(row[3]),
                                row[4].text if hasattr(row[4], "text") else str(row[4]),
                                row[5].text if hasattr(row[5], "text") else str(row[5]),
                                row[6].text if hasattr(row[6], "text") else str(row[6]),
                                row[7].text if hasattr(row[7], "text") else str(row[7]),
                            ]
                        )

                    # Create detailed table
                    detailed_table = self.create_table(
                        detailed_data,
                        repeat_rows=1,
                        auto_col_width=False,
                        col_widths=[
                            inch * x for x in [0.4, 1.0, 0.6, 0.6, 0.6, 0.6, 0.6, 3.0]
                        ],
                    )

                    # Add specific styling while preserving modern table style
                    detailed_table.setStyle(
                        TableStyle(
                            [
                                *self.MODERN_TABLE_STYLE._cmds,
                                # Right-align numeric columns (Start, End, Length, Mean CNV)
                                ("ALIGN", (2, 1), (5, -1), "RIGHT"),
                                # Center-align the state column
                                ("ALIGN", (6, 1), (6, -1), "CENTER"),
                                # Left-align the remaining columns (Chr, Region, Genes)
                                ("ALIGN", (0, 1), (1, -1), "LEFT"),
                                ("ALIGN", (7, 1), (7, -1), "LEFT"),
                            ]
                        )
                    )

                    self.elements.append(detailed_table)

            except Exception as e:
                logger.error(
                    "Error processing detailed CNV analysis: %s", str(e), exc_info=True
                )
                self.elements.append(
                    Paragraph(
                        "Error processing detailed CNV analysis data",
                        self.styles.styles["Normal"],
                    )
                )

        except Exception as e:
            logger.error("Error processing CNV section: %s", str(e), exc_info=True)
            self.elements.append(
                Paragraph(
                    "Error processing CNV analysis data",
                    self.styles.styles["Normal"],
                )
            )
