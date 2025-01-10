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
from reportlab.lib.colors import HexColor, white
from reportlab.lib.styles import ParagraphStyle
from ..sections.base import ReportSection
from ..plotting import create_CNV_plot, create_CNV_plot_per_chromosome
from robin.subpages.CNV_object import (
    Result,
    CNVAnalysis,
    CNV_Difference,
    moving_average,
    iterate_bam_bin,
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
        self.elements.append(
            Paragraph("Copy Number Variation", self.styles.styles["Heading2"])
        )
        self.elements.append(Spacer(1, 12))

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

            # Load reference CNV data
            ref_file = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)),
                "HG01280_control_new.pkl",
            )
            logger.debug("Loading reference CNV data from %s", ref_file)
            with open(ref_file, "rb") as f:
                ref_cnv_dict = pickle.load(f)
            logger.debug("Reference CNV data loaded")

            # Get reference CNV data with matching bin width
            logger.debug(
                "Getting reference CNV data with bin width %s", cnv_dict["bin_width"]
            )
            r2_cnv, _, _, _ = iterate_bam_bin(
                None,
                1,
                60,
                ref_cnv_dict,
                int(logging.ERROR),
                bin_width=cnv_dict["bin_width"],
            )

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
                    stats["gain_threshold"] = 0.5  # Default gain threshold
                    stats["loss_threshold"] = -0.5  # Default loss threshold
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
                    Paragraph(
                        str(len(total_lost_genes)), self.styles.styles["Normal"]
                    ),
                ]
            )

            # Create summary table with styling
            summary_table = Table(summary_data, colWidths=[2*inch, 1.5*inch])
            summary_table.setStyle(TableStyle([
                # Header styling
                ('BACKGROUND', (0, 0), (-1, -1), HexColor("#F5F6FA")),
                ('TEXTCOLOR', (0, 0), (-1, -1), HexColor("#2C3E50")),
                ('FONTNAME', (0, 0), (0, -1), 'Helvetica'),
                ('FONTSIZE', (0, 0), (-1, -1), 8),
                ('TOPPADDING', (0, 0), (-1, -1), 6),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
                ('LEFTPADDING', (0, 0), (-1, -1), 3),
                ('RIGHTPADDING', (0, 0), (-1, -1), 3),
                # Grid styling
                ('GRID', (0, 0), (-1, -1), 0.5, HexColor("#E2E8F0")),
                ('LINEBELOW', (0, 0), (-1, -1), 0.5, HexColor("#CBD5E1")),
                # Alignment
                ('ALIGN', (0, 0), (0, -1), 'LEFT'),
                ('ALIGN', (1, 0), (1, -1), 'RIGHT'),
                ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
                # Alternating row colors
                ('ROWBACKGROUNDS', (0, 0), (-1, -1), [HexColor("#FFFFFF"), HexColor("#F8FAFC")])
            ]))

            # Add summary to summary section
            self.summary_elements.append(
                Paragraph("Copy Number Variation", self.styles.styles["Heading3"])
            )
            self.summary_elements.append(summary_table)
            self.summary_elements.append(Spacer(1, 12))

            # Add whole chromosome events summary if any exist
            whole_chr_events = []
            for chrom, stats in chromosome_stats.items():
                if chrom != "global":
                    mean_cnv = stats["mean"]
                    if mean_cnv > stats["gain_threshold"]:
                        whole_chr_events.append(f"Chromosome {chrom[3:]}: GAIN (mean={mean_cnv:.2f})")
                    elif mean_cnv < stats["loss_threshold"]:
                        whole_chr_events.append(f"Chromosome {chrom[3:]}: LOSS (mean={mean_cnv:.2f})")
            
            if whole_chr_events:
                self.summary_elements.append(
                    Paragraph(
                        "Whole Chromosome Events: " + " | ".join(whole_chr_events),
                        self.styles.styles["Normal"]
                    )
                )

            # Start detailed analysis section
            self.elements.append(PageBreak())
            self.elements.append(
                Paragraph("Copy Number Variation Analysis", self.styles.styles["Heading2"])
            )
            self.elements.append(Spacer(1, 12))

            # Add detailed content below

            # Generate genome-wide CNV plot
            logger.debug("Generating genome-wide CNV plot")
            img_buf = create_CNV_plot(CNVresult, cnv_dict)
            width, height = inch * 7.5, inch * 2.5  # A4 width minus margins
            self.elements.append(Image(img_buf, width=width, height=height))
            self.elements.append(
                Paragraph(
                    "Copy number variation across chromosomes",
                    self.styles.styles["Caption"],
                )
            )

            # Create summary of whole chromosome events and gene-containing events
            logger.debug("Creating CNV summary")

            # Define header style for tables
            header_style = ParagraphStyle(
                "TableHeader",
                parent=self.styles.styles["Normal"],
                fontSize=8,
                fontName="Helvetica-Bold",
                textColor=self.styles.COLORS["primary"],
                alignment=1,  # Center alignment
                spaceAfter=6,
                spaceBefore=6,
            )

            whole_chr_events = []
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
                                    Paragraph(
                                        row["chrom"].replace("chr", ""),
                                        self.styles.styles["Normal"],
                                    ),
                                    Paragraph(
                                        row["name"].replace(f"{row['chrom']} ", ""),
                                        self.styles.styles["Normal"],
                                    ),
                                    Paragraph(
                                        row["cnv_state"], self.styles.styles["Normal"]
                                    ),
                                    Paragraph(
                                        f"{row['mean_cnv']:.3f}",
                                        self.styles.styles["Normal"],
                                    ),
                                    Paragraph(
                                        ", ".join(row["genes"]),
                                        ParagraphStyle(
                                            "GeneList",
                                            parent=self.styles.styles["Normal"],
                                            leading=10,  # Adjust line spacing
                                            spaceBefore=1,
                                            spaceAfter=1,
                                        ),
                                    ),
                                ]
                            )

            # Add whole chromosome events summary if any exist
            if whole_chr_events:
                self.elements.append(Spacer(1, 12))
                self.elements.append(
                    Paragraph("Whole Chromosome Events", self.styles.styles["Heading4"])
                )

                # Enhance whole chromosome events with genes
                enhanced_whole_chr_events = []
                for event in whole_chr_events:
                    chrom, state, mean_cnv = event
                    # Get all genes for this chromosome from gene_bed
                    chr_genes = []
                    if gene_bed is not None:
                        chr_genes = gene_bed[gene_bed["chrom"] == f"chr{chrom}"][
                            "gene"
                        ].tolist()
                    enhanced_whole_chr_events.append(
                        [
                            Paragraph(chrom, self.styles.styles["Normal"]),
                            Paragraph(state, self.styles.styles["Normal"]),
                            Paragraph(mean_cnv, self.styles.styles["Normal"]),
                            Paragraph(
                                ", ".join(chr_genes) if chr_genes else "No genes found",
                                ParagraphStyle(
                                    "GeneList",
                                    parent=self.styles.styles["Normal"],
                                    leading=10,  # Adjust line spacing
                                    spaceBefore=1,
                                    spaceAfter=1,
                                ),
                            ),
                        ]
                    )

                whole_chr_table = Table(
                    [
                        [
                            Paragraph("Chr", header_style),
                            Paragraph("State", header_style),
                            Paragraph("Mean CNV", header_style),
                            Paragraph("Affected Genes", header_style),
                        ]
                    ]
                    + enhanced_whole_chr_events,
                    colWidths=[inch * 0.4, inch * 0.8, inch * 0.8, inch * 4.5],
                    rowHeights=None,
                )  # Allow rows to expand based on content

                whole_chr_table.setStyle(
                    TableStyle(
                        [
                            (
                                "BACKGROUND",
                                (0, 0),
                                (-1, 0),
                                self.styles.COLORS["background"],
                            ),
                            (
                                "TEXTCOLOR",
                                (0, 0),
                                (-1, 0),
                                self.styles.COLORS["primary"],
                            ),
                            ("ALIGN", (0, 0), (-1, -1), "LEFT"),
                            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                            ("FONTSIZE", (0, 0), (-1, -1), 8),
                            ("BOTTOMPADDING", (0, 0), (-1, 0), 6),
                            ("TOPPADDING", (0, 0), (-1, -1), 3),
                            ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
                            ("LEFTPADDING", (0, 0), (-1, -1), 3),
                            ("RIGHTPADDING", (0, 0), (-1, -1), 3),
                            (
                                "GRID",
                                (0, 0),
                                (-1, -1),
                                0.5,
                                self.styles.COLORS["border"],
                            ),
                            ("VALIGN", (0, 0), (-1, -1), "TOP"),
                            (
                                "ALIGN",
                                (2, 1),
                                (2, -1),
                                "RIGHT",
                            ),  # Right align mean CNV values
                        ]
                    )
                )
                self.elements.append(whole_chr_table)

            # Add gene-containing events if any exist
            if gene_containing_events:
                self.elements.append(Spacer(1, 12))
                self.elements.append(
                    Paragraph(
                        "CNV Events Containing Genes", self.styles.styles["Heading4"]
                    )
                )
                gene_events_table = Table(
                    [
                        [
                            Paragraph("Chr", header_style),
                            Paragraph("Region", header_style),
                            Paragraph("State", header_style),
                            Paragraph("Mean CNV", header_style),
                            Paragraph("Genes", header_style),
                        ]
                    ]
                    + gene_containing_events,
                    colWidths=[
                        inch * 0.4,  # Chr
                        inch * 1.0,  # Region
                        inch * 0.6,  # State
                        inch * 0.6,  # Mean CNV
                        inch * 4.0,  # Genes
                    ],
                    rowHeights=None,
                )  # Allow rows to expand based on content

                gene_events_table.setStyle(
                    TableStyle(
                        [
                            (
                                "BACKGROUND",
                                (0, 0),
                                (-1, 0),
                                self.styles.COLORS["background"],
                            ),
                            (
                                "TEXTCOLOR",
                                (0, 0),
                                (-1, 0),
                                self.styles.COLORS["primary"],
                            ),
                            ("ALIGN", (0, 0), (-1, -1), "LEFT"),
                            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                            ("FONTSIZE", (0, 0), (-1, -1), 8),
                            ("BOTTOMPADDING", (0, 0), (-1, 0), 6),
                            ("TOPPADDING", (0, 0), (-1, -1), 3),
                            ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
                            ("LEFTPADDING", (0, 0), (-1, -1), 3),
                            ("RIGHTPADDING", (0, 0), (-1, -1), 3),
                            (
                                "GRID",
                                (0, 0),
                                (-1, -1),
                                0.5,
                                self.styles.COLORS["border"],
                            ),
                            ("VALIGN", (0, 0), (-1, -1), "TOP"),
                            # Align numeric columns to the right
                            ("ALIGN", (3, 1), (3, -1), "RIGHT"),
                            # Center the state column
                            ("ALIGN", (2, 1), (2, -1), "CENTER"),
                        ]
                    )
                )
                self.elements.append(gene_events_table)

            # Add note about detailed view
            self.elements.append(Spacer(1, 12))
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

            # Add Detailed CNV Analysis Section
            self.elements.append(PageBreak())
            self.elements.append(
                Paragraph("Detailed CNV Analysis", self.styles.styles["Heading2"])
            )
            self.elements.append(Spacer(1, 12))

            # Add Analysis Summary table here
            self.elements.append(
                Paragraph("Analysis Parameters", self.styles.styles["Heading3"])
            )
            self.elements.append(Spacer(1, 6))
            self.elements.append(summary_table)
            self.elements.append(Spacer(1, 12))

            # Initialize variables for plot layout
            self.current_row = []
            plots_per_row = 2

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

                    detailed_table = Table(
                        [
                            [
                                Paragraph("Chr", header_style),
                                Paragraph("Region", header_style),
                                Paragraph("Start (Mb)", header_style),
                                Paragraph("End (Mb)", header_style),
                                Paragraph("Length (Mb)", header_style),
                                Paragraph("Mean CNV", header_style),
                                Paragraph("State", header_style),
                                Paragraph("Genes", header_style),
                            ]
                        ]
                        + formatted_events,
                        colWidths=[
                            inch * 0.4,  # Chr
                            inch * 1.0,  # Region
                            inch * 0.6,  # Start
                            inch * 0.6,  # End
                            inch * 0.6,  # Length
                            inch * 0.6,  # Mean CNV
                            inch * 0.6,  # State
                            inch * 3.0,  # Genes
                        ],
                        rowHeights=None,  # Allow rows to expand based on content
                        repeatRows=1,  # Repeat header row on new pages
                    )

                    detailed_table.setStyle(
                        TableStyle(
                            [
                                # Header styling
                                (
                                    "BACKGROUND",
                                    (0, 0),
                                    (-1, 0),
                                    self.styles.COLORS["background"],
                                ),
                                (
                                    "TEXTCOLOR",
                                    (0, 0),
                                    (-1, 0),
                                    self.styles.COLORS["primary"],
                                ),
                                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                                ("FONTSIZE", (0, 0), (-1, -1), 8),
                                # Cell padding
                                ("TOPPADDING", (0, 0), (-1, -1), 3),
                                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
                                ("LEFTPADDING", (0, 0), (-1, -1), 3),
                                ("RIGHTPADDING", (0, 0), (-1, -1), 3),
                                # Grid and alignment
                                (
                                    "GRID",
                                    (0, 0),
                                    (-1, -1),
                                    0.5,
                                    self.styles.COLORS["border"],
                                ),
                                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                                (
                                    "ALIGN",
                                    (0, 0),
                                    (-1, -1),
                                    "LEFT",
                                ),  # Default left alignment
                                # Specific column alignments
                                (
                                    "ALIGN",
                                    (2, 1),
                                    (5, -1),
                                    "RIGHT",
                                ),  # Numeric columns right-aligned
                                (
                                    "ALIGN",
                                    (6, 1),
                                    (6, -1),
                                    "CENTER",
                                ),  # State column centered
                                # Alternating row colors for better readability
                                (
                                    "ROWBACKGROUNDS",
                                    (0, 0),
                                    (-1, -1),
                                    [HexColor("#FFFFFF"), HexColor("#F8FAFC")],
                                ),
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
