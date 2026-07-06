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

# from robin.subpages.CNVObjectClass import (
#    CNVAnalysis
# )

from robin.analysis.cnv_analysis import (
    Result,
    moving_average,
    CNV_Difference,
    compute_cnv_log2_from_ploidy,
    prepare_cnv_calling_track,
    resolve_cnv_calling_bin_width,
)
from robin.analysis.cnv_classification import detect_cnv_events, get_cnv_summary, CNVEvent
from robin.analysis.cnv_regional import (
    SIGNIFICANT_CNV_STATES,
    analyze_cytoband_cnv,
    build_regional_cnv_events,
    build_significant_regions,
    format_chromosome_cnv_status,
    format_panel_genes_for_table,
    is_reportable_chromosome,
    load_panel_gene_bed,
    load_target_coverage_df,
    panel_genes_in_region,
)
from robin.classification_config import get_cnv_thresholds, is_resolution_sufficient

from robin import resources

logger = logging.getLogger(__name__)


def calculate_chromosome_stats(result, ref_result, XYestimate):
    """Calculate chromosome-wide statistics and baselines.

    Args:
        result: CNV result object for sample
        ref_result: CNV result object for reference

    Returns:
        Dictionary of chromosome statistics including means, baselines, and thresholds
    """
    stats = {}
    autosome_means = []

    # Calculate normalized values and stats for each chromosome
    for chrom in result.cnv.keys():
        if chrom != "chrM" and chrom in ref_result:
            # Calculate normalized CNV values
            sample_avg = moving_average(result.cnv[chrom])
            ref_avg = moving_average(ref_result[chrom])

            # Pad arrays if needed
            max_len = max(len(sample_avg), len(ref_avg))
            if len(sample_avg) < max_len:
                sample_avg = np.pad(sample_avg, (0, max_len - len(sample_avg)))
            if len(ref_avg) < max_len:
                ref_avg = np.pad(ref_avg, (0, max_len - len(ref_avg)))

            # Calculate normalized CNV
            normalized_cnv = sample_avg - ref_avg

            # Calculate basic statistics
            chr_mean = np.mean(normalized_cnv)
            chr_std = np.std(normalized_cnv)

            # Store autosome means for global statistics
            if chrom.startswith("chr") and chrom[3:].isdigit():
                autosome_means.append(chr_mean)

            # Set baseline and thresholds based on chromosome and sex
            if chrom == "chrX":
                if XYestimate == "XX":  # Female
                    baseline = 1.0  # Expected +1 relative to male control
                else:  # Male
                    baseline = 0.0  # Expected same as male control
            elif chrom == "chrY":
                if XYestimate == "XY":  # Male
                    baseline = 0.0
                else:  # Female
                    baseline = -1.0  # Expected absence
            else:  # Autosomes
                baseline = 0.0

            stats[chrom] = {
                "mean": chr_mean,
                "std": chr_std,
                "baseline": baseline,
                "normalized_cnv": normalized_cnv,
            }

    # Calculate global autosome statistics
    global_mean = np.mean(autosome_means)
    global_std = np.std(autosome_means)

    # Store global stats
    stats["global"] = {"mean": global_mean, "std": global_std}

    # self.chromosome_stats = stats
    return stats


class CNVSection(ReportSection):
    """Section containing the CNV analysis."""

    FULL_PLOT_WIDTH = inch * 7.5
    FULL_PLOT_HEIGHT = inch * 2.5

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
            # cnv_analyzer = CNVAnalysis(target_panel="rCNS2")
            # cnv_analyzer.XYestimate = XYestimate

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
            # cnv_analyzer.gene_bed = gene_bed
            # cnv_analyzer.cytobands_bed = cytobands_bed
            # cnv_analyzer.cnv_dict = cnv_dict

            # Get reference CNV data with matching bin width
            logger.debug(
                "Getting reference CNV data with bin width %s", cnv_dict["bin_width"]
            )

            r2_cnv = Result(
                np.load(
                    os.path.join(self.report.output, "CNV2.npy"), allow_pickle="TRUE"
                ).item()
            ).cnv

            use_normalized_summary = getattr(
                self.report, "cnv_summary_normalized", False
            )
            log2_cnv = (
                compute_cnv_log2_from_ploidy(CNVresult.cnv, XYestimate)
                if use_normalized_summary
                else None
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
            # cnv_analyzer.result3 = result3

            # Calculate chromosome statistics using CNVAnalysis logic
            chromosome_stats = calculate_chromosome_stats(CNVresult, r2_cnv, XYestimate)
            # cnv_analyzer.chromosome_stats = chromosome_stats

            # Add gain/loss thresholds to chromosome stats using centralized rules
            for chrom, stats in chromosome_stats.items():
                if chrom != "global":
                    gain_threshold, loss_threshold = get_cnv_thresholds(chrom, XYestimate)
                    stats["gain_threshold"] = gain_threshold
                    stats["loss_threshold"] = loss_threshold

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
                        f"{cnv_dict.get('variance', 0):.2f}",
                        self.styles.styles["Normal"],
                    ),
                ]
            )
            centromeres_file = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)),
                "cenSatRegions.bed",
            )

            centromere_bed = pd.read_csv(
                centromeres_file,
                usecols=[0, 1, 2, 3],
                names=["chrom", "start_pos", "end_pos", "name"],
                header=None,
                sep=r"\s+",
            )

            panel_name, panel_genes_df = load_panel_gene_bed(self.report.output)
            target_coverage_df = load_target_coverage_df(self.report.output)
            reportable_chromosomes = [
                chrom
                for chrom in natsort.natsorted(result3.cnv.keys())
                if is_reportable_chromosome(chrom)
            ]
            cytoband_analysis_by_chrom: dict[str, pd.DataFrame] = {}
            regional_cnv_events: list[dict] = []
            for chrom in reportable_chromosomes:
                cytoband_analysis = analyze_cytoband_cnv(
                    result3.cnv,
                    chrom,
                    cnv_dict,
                    cytobands_bed,
                    centromere_bed,
                    gene_bed if gene_bed is not None else pd.DataFrame(
                        columns=["chrom", "start_pos", "end_pos", "gene"]
                    ),
                    XYestimate,
                )
                cytoband_analysis_by_chrom[chrom] = cytoband_analysis
                regional_cnv_events.extend(
                    build_regional_cnv_events(cytoband_analysis, panel_genes_df)
                )

            # Calculate gene counts
            total_gained_genes = set()
            total_lost_genes = set()
            for chrom in natsort.natsorted(result3.cnv.keys()):
                if chrom != "chrM" and re.match(r"^chr(\d+|X|Y)$", chrom):
                    analysis = cytoband_analysis_by_chrom.get(chrom)
                    if analysis is None:
                        analysis = analyze_cytoband_cnv(
                            result3.cnv,
                            chrom,
                            cnv_dict,
                            cytobands_bed,
                            centromere_bed,
                            gene_bed if gene_bed is not None else pd.DataFrame(
                                columns=["chrom", "start_pos", "end_pos", "gene"]
                            ),
                            XYestimate,
                        )
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
                    formatted_summary_data,
                    auto_col_width=True,
                    compact=True,
                    font_size=9,
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
                                9,
                            ),  # Unified font size
                            (
                                "TOPPADDING",
                                (0, 0),
                                (-1, -1),
                                4,
                            ),
                            (
                                "BOTTOMPADDING",
                                (0, 0),
                                (-1, -1),
                                4,
                            ),
                        ]
                    )
                )
                self.elements.append(summary_table)

            # Detect CNV events using centralized classification rules
            logger.info("Detecting CNV events using centralized rules")
            events = []
            
            # Check if resolution is sufficient
            analysis_binw = int(cnv_dict.get("bin_width", 1000000))
            calling_binw = resolve_cnv_calling_bin_width(analysis_binw)
            if not is_resolution_sufficient(calling_binw):
                logger.warning("Resolution insufficient for CNV calling")
                summary_whole_chr_events = []
                summary_arm_events = []
            else:
                # Arm/whole-chromosome events: log2(ploidy / expected), ≥1 Mb bins
                analysis_binw = int(cnv_dict.get("bin_width", 1000000))
                calling_cnv, calling_binw = prepare_cnv_calling_track(
                    CNVresult.cnv, analysis_binw, XYestimate
                )
                events = detect_cnv_events(
                    cnv_data=calling_cnv,
                    bin_width=calling_binw,
                    sex_estimate=XYestimate,
                    cytobands_df=cytobands_bed,
                    gene_df=gene_bed
                )
                
                # Convert events to summary format
                summary_whole_chr_events = []
                summary_arm_events = []
                
                for event in events:
                    if event.event_type.startswith("WHOLE_CHR_"):
                        event_type = event.event_type.replace("WHOLE_CHR_", "")
                        summary_whole_chr_events.append(
                            f"Chromosome {event.chromosome[3:]}: {event_type} (mean={event.mean_cnv:.2f})"
                        )
                        logger.info(f"Detected whole chromosome {event_type} for {event.chromosome}")
                    else:
                        arm_label = f"{event.arm}-arm" if event.arm else "arm"
                        summary_arm_events.append(
                            f"Chromosome {event.chromosome[3:]} {arm_label}: {event.event_type} (mean={event.mean_cnv:.2f}, {event.proportion_affected:.0%} of arm)"
                        )
                        logger.info(f"Detected arm event: {event.chromosome} {arm_label} {event.event_type}")
                
                # Log the final counts
                logger.info(f"Found {len(summary_arm_events)} arm events")
                logger.info(f"Found {len(summary_whole_chr_events)} whole chromosome events")

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
                        "Chromosome Arm Events (requires visual inspection):<br/> "
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
            from robin.gui.plotting_preferences import cnv_report_plot_caption

            significant_regions: dict[str, list[dict]] = {}
            for chrom in reportable_chromosomes:
                cytoband_analysis = cytoband_analysis_by_chrom.get(chrom)
                if cytoband_analysis is not None and not cytoband_analysis.empty:
                    region_list = build_significant_regions(cytoband_analysis)
                    if region_list:
                        significant_regions[chrom] = region_list

            use_normalized_summary = getattr(
                self.report, "cnv_summary_normalized", False
            )
            img_buf = create_CNV_plot(
                CNVresult,
                cnv_dict,
                normalized_cnv=log2_cnv,
                use_normalized_difference=use_normalized_summary,
                sex_estimate=str(XYestimate),
                panel_genes_df=panel_genes_df,
                target_coverage_df=target_coverage_df,
                significant_regions=significant_regions,
            )
            summary_caption_scale = (
                "normalized_difference" if use_normalized_summary else "ploidy"
            )
            width, height = inch * 7.5, inch * 2  # A4 width minus margins
            self.summary_elements.append(Image(img_buf, width=width, height=height))
            self.summary_elements.append(
                Paragraph(
                    cnv_report_plot_caption(summary_caption_scale),
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

            # Create summary of CNV events using centralized detection
            logger.debug("Creating CNV summary using centralized events")

            # Use the events detected above
            whole_chr_events = []
            arm_events = []

            for event in events:
                if event.event_type.startswith("WHOLE_CHR_"):
                    event_type = event.event_type.replace("WHOLE_CHR_", "")
                    whole_chr_events.append([
                        event.chromosome.replace("chr", ""),
                        event_type,
                        f"{event.mean_cnv:.2f}",
                    ])
                else:
                    arm_label = f"{event.arm}-arm" if event.arm else "arm"
                    arm_events.append([
                        event.chromosome.replace("chr", ""),
                        arm_label,
                        event.event_type,
                        f"{event.mean_cnv:.2f}",
                        f"{event.proportion_affected:.1%}",
                    ])

            # Add whole chromosome events summary if any exist
            if whole_chr_events:
                self.elements.append(
                    Paragraph(
                        "Whole Chromosome Events",
                        ParagraphStyle(
                            "CNVTableTitle",
                            parent=self.styles.styles["Normal"],
                            fontSize=9,
                            fontName="Helvetica-Bold",
                            spaceAfter=4,
                        ),
                    )
                )
                self.elements.append(Spacer(1, 2))
                whole_chr_data = [["Chr", "State", "Mean CNV"]]
                whole_chr_data.extend(whole_chr_events)
                whole_chr_table = self.create_table(
                    whole_chr_data,
                    repeat_rows=1,
                    auto_col_width=False,
                    col_widths=[inch * x for x in [0.4, 0.8, 0.8]],
                    compact=True,
                    font_size=9,
                )
                whole_chr_table.setStyle(
                    TableStyle(
                        [
                            *self.MODERN_TABLE_STYLE._cmds,
                            ("FONTSIZE", (0, 0), (-1, -1), 9),
                            ("ALIGN", (2, 1), (2, -1), "RIGHT"),
                            ("ALIGN", (1, 1), (1, -1), "CENTER"),
                            ("TOPPADDING", (0, 0), (-1, -1), 4),
                            ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
                        ]
                    )
                )
                self.elements.append(whole_chr_table)
                self.elements.append(Spacer(1, 4))

            # Build arm and regional event tables. Always stack vertically (never nested
            # side-by-side) so each table can split across pages. Nested tables
            # cannot split and cause LayoutError when content exceeds frame height.
            arm_col_widths = [inch * x for x in [0.35, 0.55, 0.5, 0.5, 0.8]]
            regional_col_widths = [
                inch * x for x in [0.35, 1.1, 0.55, 0.55, 0.55, 0.55, 0.55, 1.35]
            ]

            arm_header = None
            arm_table = None
            regional_header = None
            regional_table = None

            if regional_cnv_events:
                regional_data = [[
                    "Chr",
                    "Region",
                    "Start (Mb)",
                    "End (Mb)",
                    "Length (Mb)",
                    "Mean CNV",
                    "State",
                    "Panel genes",
                ]]
                for event in regional_cnv_events:
                    panel_gene_text = format_panel_genes_for_table(event["panel_genes"])
                    regional_data.append([
                        event["chrom"],
                        event["region"],
                        f"{event['start_mb']:.2f}",
                        f"{event['end_mb']:.2f}",
                        f"{event['length_mb']:.2f}",
                        f"{event['mean_cnv']:.2f}",
                        event["state"],
                        panel_gene_text,
                    ])
                regional_table = self.create_table(
                    regional_data,
                    repeat_rows=1,
                    auto_col_width=False,
                    col_widths=regional_col_widths,
                    compact=True,
                    font_size=9,
                )
                regional_table.setStyle(
                    TableStyle(
                        [
                            *self.MODERN_TABLE_STYLE._cmds,
                            ("FONTSIZE", (0, 0), (-1, -1), 9),
                            ("ALIGN", (2, 1), (5, -1), "RIGHT"),
                            ("ALIGN", (6, 1), (6, -1), "CENTER"),
                            ("TOPPADDING", (0, 0), (-1, -1), 4),
                            ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
                        ]
                    )
                )
                regional_title = "Regional CNV Events"
                if panel_name:
                    regional_title += f" ({panel_name} panel genes)"
                regional_header = Paragraph(
                    regional_title,
                    ParagraphStyle(
                        "CNVTableTitle",
                        parent=self.styles.styles["Normal"],
                        fontSize=9,
                        fontName="Helvetica-Bold",
                        spaceAfter=4,
                    ),
                )

            if arm_events:
                arm_data = [["Chr", "Arm", "State", "Mean CNV", "Proportion Affected"]]
                arm_data.extend(arm_events)
                arm_table = self.create_table(
                    arm_data,
                    repeat_rows=1,
                    auto_col_width=False,
                    col_widths=arm_col_widths,
                    compact=True,
                    font_size=9,
                )
                arm_table.setStyle(
                    TableStyle(
                        [
                            *self.MODERN_TABLE_STYLE._cmds,
                            ("FONTSIZE", (0, 0), (-1, -1), 9),
                            ("ALIGN", (3, 1), (3, -1), "RIGHT"),
                            ("ALIGN", (2, 1), (2, -1), "CENTER"),
                            ("ALIGN", (4, 1), (4, -1), "RIGHT"),
                            ("TOPPADDING", (0, 0), (-1, -1), 4),
                            ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
                        ]
                    )
                )
                arm_header = Paragraph(
                    "Arm Events (visual inspection)",
                    ParagraphStyle(
                        "CNVTableTitle",
                        parent=self.styles.styles["Normal"],
                        fontSize=9,
                        fontName="Helvetica-Bold",
                        spaceAfter=4,
                    ),
                )

            # Stack vertically as separate flowables so each table can split across pages
            if regional_header is not None:
                self.elements.append(regional_header)
                self.elements.append(regional_table)
            if arm_header is not None:
                if regional_header is not None:
                    self.elements.append(Spacer(1, 8))
                self.elements.append(arm_header)
                self.elements.append(arm_table)

            # Add note about detailed view
            self.elements.append(
                Paragraph(
                    (
                        "Note: Regional events are derived from merged cytoband analysis. "
                        "Panel genes are targets from the sample analysis panel overlapping "
                        "each called region. Arm-level events require visual inspection."
                    ),
                    ParagraphStyle(
                        "Note",
                        parent=self.styles.styles["Normal"],
                        fontSize=9,
                        textColor=self.styles.COLORS["text"],
                        spaceBefore=3,
                        spaceAfter=3,
                        italics=True,
                    ),
                )
            )

            # Add individual chromosome plots at full page width
            try:
                chromosome_status = {}

                for chrom in reportable_chromosomes:
                    cytoband_analysis = cytoband_analysis_by_chrom[chrom]
                    chromosome_status[chrom] = format_chromosome_cnv_status(
                        chrom,
                        events,
                        cytoband_analysis,
                    )

                if panel_name and not panel_genes_df.empty:
                    panel_plot_blurb = (
                        "Scatter points show bin-level log2(ploidy / expected copy number); "
                        "the dark trace is a rolling median. Panel target lollipops (right axis) show "
                        "per-target sequencing coverage from 0 to the chromosome maximum."
                        if use_normalized_summary
                        else (
                            "Scatter points show bin-level copy number; the dark "
                            "trace is a rolling median. Panel target lollipops (right axis) show "
                            "per-target sequencing coverage from 0 to the chromosome maximum."
                        )
                    )
                    self.elements.append(
                        Paragraph(
                            (
                                f"Individual chromosome plots include lollipop markers for "
                                f"genes in the <b>{panel_name}</b> target panel "
                                f"({len(panel_genes_df)} genes). "
                                f"{panel_plot_blurb} "
                                "Panel targets are lollipops; "
                                "genes are labelled when >3 SD from the chromosome mean."
                            ),
                            ParagraphStyle(
                                "PanelGeneLegend",
                                parent=self.styles.styles["Normal"],
                                fontSize=9,
                                textColor=self.styles.COLORS["text"],
                                spaceBefore=3,
                                spaceAfter=6,
                            ),
                        )
                    )

                # Generate all chromosome plots at once
                logger.debug(
                    "Generating individual chromosome plots for %d chromosomes",
                    len(reportable_chromosomes),
                )
                chromosome_plots = create_CNV_plot_per_chromosome(
                    CNVresult,
                    cnv_dict,
                    significant_regions=significant_regions,
                    chromosomes=reportable_chromosomes,
                    panel_genes_df=panel_genes_df,
                    chromosome_status=chromosome_status,
                    normalized_cnv=log2_cnv,
                    target_coverage_df=target_coverage_df,
                    use_log2_ratio=use_normalized_summary,
                )
                plot_lookup = dict(chromosome_plots)
                plotted_chromosomes = [
                    chrom for chrom in reportable_chromosomes if chrom in plot_lookup
                ]

                for chrom in plotted_chromosomes:
                    img_buf = plot_lookup[chrom]
                    self.elements.append(
                        Image(
                            img_buf,
                            width=self.FULL_PLOT_WIDTH,
                            height=self.FULL_PLOT_HEIGHT,
                        )
                    )
                    self.elements.append(Spacer(1, 10))

                # Add detailed CNV table
                self.elements.append(Spacer(1, 6))
                self.elements.append(
                    Paragraph("Detailed CNV Events", self.styles.styles["Heading3"])
                )

                # Create detailed table from regional and arm/whole-chromosome events
                all_cnv_events = []
                for event in regional_cnv_events:
                    all_cnv_events.append([
                        event["chrom"],
                        event["region"],
                        f"{event['start_mb']:.2f}",
                        f"{event['end_mb']:.2f}",
                        f"{event['length_mb']:.2f}",
                        f"{event['mean_cnv']:.2f}",
                        event["state"],
                        format_panel_genes_for_table(event["panel_genes"]),
                    ])
                for event in events:
                    region_name = (
                        f"{event.chromosome} {event.arm}-arm"
                        if event.arm
                        else f"{event.chromosome} whole chromosome"
                    )
                    panel_genes = panel_genes_in_region(
                        panel_genes_df,
                        event.chromosome,
                        event.start_pos,
                        event.end_pos,
                    )
                    all_cnv_events.append([
                        event.chromosome.replace("chr", ""),
                        region_name.replace(f"{event.chromosome} ", ""),
                        f"{event.start_pos/1e6:.2f}",
                        f"{event.end_pos/1e6:.2f}",
                        f"{event.length/1e6:.2f}",
                        f"{event.mean_cnv:.2f}",
                        event.event_type.replace("WHOLE_CHR_", ""),
                        format_panel_genes_for_table(panel_genes),
                    ])

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
                                ),  # Panel genes
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
                            "Panel genes",
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

                    # Create detailed table (compact)
                    detailed_table = self.create_table(
                        detailed_data,
                        repeat_rows=1,
                        auto_col_width=False,
                        col_widths=[
                            inch * x for x in [0.4, 1.0, 0.6, 0.6, 0.6, 0.6, 0.6, 3.0]
                        ],
                        compact=True,
                        font_size=9,
                    )

                    # Add specific styling while preserving modern table style
                    detailed_table.setStyle(
                        TableStyle(
                            [
                                *self.MODERN_TABLE_STYLE._cmds,
                                ("ALIGN", (2, 1), (5, -1), "RIGHT"),
                                ("ALIGN", (6, 1), (6, -1), "CENTER"),
                                ("ALIGN", (0, 1), (1, -1), "LEFT"),
                                ("ALIGN", (7, 1), (7, -1), "LEFT"),
                                ("FONTSIZE", (0, 0), (-1, -1), 9),
                                ("TOPPADDING", (0, 0), (-1, -1), 4),
                                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
                            ]
                        )
                    )

                    self.elements.append(detailed_table)

                    # Build export DataFrames for CSVs
                    try:
                        # Whole chromosome events
                        if summary_whole_chr_events or whole_chr_events:
                            # whole_chr_events is already a list of lists
                            df_whole = pd.DataFrame(
                                whole_chr_events, columns=["Chr", "State", "Mean CNV"]
                            )
                            if not df_whole.empty:
                                df_whole["Mean CNV"] = pd.to_numeric(
                                    df_whole["Mean CNV"], errors="coerce"
                                )
                            self.export_frames["cnv_whole_chromosome_events"] = df_whole

                        # Chromosome arm events
                        if arm_events:
                            df_arm = pd.DataFrame(
                                arm_events,
                                columns=[
                                    "Chr",
                                    "Arm",
                                    "State",
                                    "Mean CNV",
                                    "Proportion Affected",
                                ],
                            )
                            if not df_arm.empty:
                                df_arm["Mean CNV"] = pd.to_numeric(
                                    df_arm["Mean CNV"], errors="coerce"
                                )
                                # Convert percentage strings like "70%" to 0.7
                                df_arm["Proportion Affected"] = (
                                    df_arm["Proportion Affected"]
                                    .astype(str)
                                    .str.rstrip("%")
                                )
                                df_arm["Proportion Affected"] = (
                                    pd.to_numeric(
                                        df_arm["Proportion Affected"], errors="coerce"
                                    )
                                    / 100.0
                                )
                            self.export_frames["cnv_arm_events"] = df_arm

                        # Regional CNV events
                        if regional_cnv_events:
                            df_regional = pd.DataFrame(
                                [
                                    {
                                        "Chr": event["chrom"],
                                        "Region": event["region"],
                                        "Start (Mb)": event["start_mb"],
                                        "End (Mb)": event["end_mb"],
                                        "Length (Mb)": event["length_mb"],
                                        "Mean CNV": event["mean_cnv"],
                                        "State": event["state"],
                                        "Panel genes": format_panel_genes_for_table(
                                            event["panel_genes"]
                                        ),
                                    }
                                    for event in regional_cnv_events
                                ]
                            )
                            self.export_frames["cnv_regional_events"] = df_regional

                        # Detailed CNV events
                        if all_cnv_events:
                            df_detail = pd.DataFrame(
                                all_cnv_events,
                                columns=[
                                    "Chr",
                                    "Region",
                                    "Start (Mb)",
                                    "End (Mb)",
                                    "Length (Mb)",
                                    "Mean CNV",
                                    "State",
                                    "Panel genes",
                                ],
                            )
                            if not df_detail.empty:
                                for col in [
                                    "Start (Mb)",
                                    "End (Mb)",
                                    "Length (Mb)",
                                    "Mean CNV",
                                ]:
                                    df_detail[col] = pd.to_numeric(
                                        df_detail[col], errors="coerce"
                                    )
                            self.export_frames["cnv_detailed_events"] = df_detail
                    except Exception as ex:
                        logger.error(
                            "Error building CNV export DataFrames: %s",
                            str(ex),
                            exc_info=True,
                        )

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
