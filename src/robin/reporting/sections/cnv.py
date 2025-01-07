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
from robin.subpages.CNV_object import Result, CNVAnalysis, CNV_Difference, moving_average, iterate_bam_bin
from robin import resources

logger = logging.getLogger(__name__)

class CNVSection(ReportSection):
    """Section containing the CNV analysis."""

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
                os.path.dirname(os.path.abspath(resources.__file__)), 
                "unique_genes.bed"
            )
            cytoband_file = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)), 
                "cytoBand.txt"
            )
            logger.debug("Resource files: gene_bed=%s, cytoband=%s", gene_bed_file, cytoband_file)
            
            # Load gene and cytoband data
            gene_bed = None
            cytobands_bed = None
            if os.path.exists(gene_bed_file):
                gene_bed = pd.read_csv(
                    gene_bed_file,
                    sep='\t',
                    names=['chrom', 'start_pos', 'end_pos', 'gene']
                )
                logger.debug("Loaded gene bed file with shape: %s", gene_bed.shape)
            if os.path.exists(cytoband_file):
                cytobands_bed = pd.read_csv(
                    cytoband_file,
                    sep='\t',
                    names=['chrom', 'start_pos', 'end_pos', 'name', 'stain']
                )
                logger.debug("Loaded cytoband file with shape: %s", cytobands_bed.shape)

            # Set up CNVAnalysis object with loaded data
            cnv_analyzer.gene_bed = gene_bed
            cnv_analyzer.cytobands_bed = cytobands_bed
            cnv_analyzer.cnv_dict = cnv_dict

            # Load reference CNV data
            ref_file = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)),
                "HG01280_control_new.pkl"
            )
            logger.debug("Loading reference CNV data from %s", ref_file)
            with open(ref_file, "rb") as f:
                ref_cnv_dict = pickle.load(f)
            logger.debug("Reference CNV data loaded")

            # Get reference CNV data with matching bin width
            logger.debug("Getting reference CNV data with bin width %s", cnv_dict["bin_width"])
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
                                moving_avg_data1 = np.pad(moving_avg_data1, (0, max_len - len(moving_avg_data1)))
                            if len(moving_avg_data2) < max_len:
                                moving_avg_data2 = np.pad(moving_avg_data2, (0, max_len - len(moving_avg_data2)))
                        # Calculate difference
                        result3.cnv[key] = moving_avg_data1 - moving_avg_data2

            # Set the result3 in the analyzer
            cnv_analyzer.result3 = result3

            # Calculate chromosome statistics using CNVAnalysis logic
            chromosome_stats = cnv_analyzer.calculate_chromosome_stats(CNVresult, r2_cnv)
            cnv_analyzer.chromosome_stats = chromosome_stats
            
            # Add gain/loss thresholds to chromosome stats
            for chrom, stats in chromosome_stats.items():
                if chrom != 'global':
                    stats['gain_threshold'] = 0.5  # Default gain threshold
                    stats['loss_threshold'] = -0.5  # Default loss threshold
                    if chrom == 'chrX':
                        if XYestimate == 'XY':  # Male
                            stats['gain_threshold'] = 0.3
                            stats['loss_threshold'] = -0.3
                        else:  # Female
                            stats['gain_threshold'] = 0.75
                            stats['loss_threshold'] = -0.75
                    elif chrom == 'chrY':
                        if XYestimate == 'XY':  # Male
                            stats['gain_threshold'] = 0.5
                            stats['loss_threshold'] = -0.5
                        else:  # Female
                            stats['gain_threshold'] = -0.2
                            stats['loss_threshold'] = -1.0

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

            # Add genetic sex information
            if XYestimate != "Unknown":
                self.elements.append(
                    Paragraph(
                        f"Estimated Genetic Sex: {XYestimate}",
                        ParagraphStyle(
                            'GeneticSex',
                            parent=self.styles.styles['Normal'],
                            fontSize=8,
                            textColor=self.styles.COLORS["text"],
                            spaceBefore=2,
                            spaceAfter=6
                        )
                    )
                )

            # Create detailed CNV analysis table
            logger.debug("Creating detailed CNV analysis table")
            all_cytoband_analysis = []
            for chrom in natsort.natsorted(CNVresult.cnv.keys()):
                if chrom != "chrM" and re.match(r"^chr(\d+|X|Y)$", chrom):
                    # Get CNV analysis for this chromosome using the CNVAnalysis class method
                    cytoband_analysis = cnv_analyzer.analyze_cytoband_cnv(result3.cnv, chrom)
                    if not cytoband_analysis.empty:
                        all_cytoband_analysis.append(cytoband_analysis)

            if all_cytoband_analysis:
                combined_analysis = pd.concat(all_cytoband_analysis)
                # Create table data with headers
                table_data = [["Chr", "Cytoband", "Start (Mb)", "End (Mb)", "Length (Mb)", "Mean CNV", "State", "Genes"]]
                
                # Add rows
                for _, row in combined_analysis.iterrows():
                    if row['cnv_state'] in ['GAIN', 'LOSS']:  # Only include significant events
                        genes_str = ', '.join(row['genes']) if isinstance(row['genes'], list) else ''
                        table_data.append([
                            row['chrom'].replace('chr', ''),  # Remove 'chr' prefix to save space
                            row['name'],
                            f"{row['start_pos']/1e6:.2f}",  # Convert to Mb
                            f"{row['end_pos']/1e6:.2f}",    # Convert to Mb
                            f"{row['length']/1e6:.2f}",     # Convert to Mb
                            f"{row['mean_cnv']:.3f}",
                            row['cnv_state'],
                            genes_str
                        ])

                # Create table style
                style = TableStyle([
                    ('BACKGROUND', (0, 0), (-1, 0), self.styles.COLORS["background"]),
                    ('TEXTCOLOR', (0, 0), (-1, 0), self.styles.COLORS["primary"]),
                    ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                    ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                    ('FONTSIZE', (0, 0), (-1, 0), 7),  # Smaller header font
                    ('BOTTOMPADDING', (0, 0), (-1, 0), 6),  # Reduced padding
                    ('TOPPADDING', (0, 0), (-1, 0), 6),    # Reduced padding
                    ('BACKGROUND', (0, 1), (-1, -1), white),
                    ('TEXTCOLOR', (0, 1), (-1, -1), self.styles.COLORS["text"]),
                    ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
                    ('FONTSIZE', (0, 1), (-1, -1), 7),  # Smaller body font
                    ('GRID', (0, 0), (-1, -1), 0.5, self.styles.COLORS["border"]),  # Thinner grid lines
                    # Align numeric columns to the right
                    ('ALIGN', (2, 1), (5, -1), 'RIGHT'),
                    # Center the CNV state column
                    ('ALIGN', (6, 1), (6, -1), 'CENTER'),
                    # Enable text wrapping for all cells
                    ('LEFTPADDING', (0, 0), (-1, -1), 4),   # Reduced padding
                    ('RIGHTPADDING', (0, 0), (-1, -1), 4),  # Reduced padding
                    ('TOPPADDING', (0, 1), (-1, -1), 4),    # Reduced padding
                    ('BOTTOMPADDING', (0, 1), (-1, -1), 4), # Reduced padding
                ])

                # Add color coding for gains and losses
                for i in range(1, len(table_data)):
                    if table_data[i][6] == 'GAIN':
                        style.add('TEXTCOLOR', (6, i), (6, i), self.styles.COLORS["success"])  # Green for gains
                    elif table_data[i][6] == 'LOSS':
                        style.add('TEXTCOLOR', (6, i), (6, i), self.styles.COLORS["error"])  # Red for losses

                # Create and add the table with adjusted column widths
                table = Table(
                    table_data,
                    colWidths=[
                        inch * 0.4,  # Chr (smaller)
                        inch * 0.7,  # Cytoband
                        inch * 0.7,  # Start
                        inch * 0.7,  # End
                        inch * 0.7,  # Length
                        inch * 0.7,  # Mean CNV
                        inch * 0.6,  # State
                        inch * 2.5   # Genes (wider to accommodate wrapped text)
                    ],
                    repeatRows=1  # Repeat header row on new pages
                )
                table.setStyle(style)
                
                self.elements.append(Spacer(1, 12))
                self.elements.append(
                    Paragraph("CNV Analysis Results", self.styles.styles["Heading3"])
                )
                self.elements.append(Spacer(1, 6))
                self.elements.append(table)

            # Generate individual chromosome plots
            logger.debug("Generating individual chromosome plots")
            
            # Generate chromosome events from the CNV analysis
            chromosome_events = []
            for chrom in natsort.natsorted(result3.cnv.keys()):
                if chrom != "chrM" and re.match(r"^chr(\d+|X|Y)$", chrom):
                    if chrom in chromosome_stats:
                        stats = chromosome_stats[chrom]
                        # Check for gains
                        if np.any(result3.cnv[chrom] > stats['gain_threshold']):
                            chromosome_events.append({
                                'chromosome': chrom,
                                'type': 'GAIN'
                            })
                        # Check for losses
                        if np.any(result3.cnv[chrom] < stats['loss_threshold']):
                            chromosome_events.append({
                                'chromosome': chrom,
                                'type': 'LOSS'
                            })
            
            # Generate plots for chromosomes in a grid layout
            plots_data = []  # Store plot data for grid arrangement
            for chrom in natsort.natsorted(CNVresult.cnv.keys()):
                if chrom != "chrM" and re.match(r"^chr(\d+|X|Y)$", chrom):
                    if chrom in chromosome_stats:
                        # Create a Result object with just this chromosome's data
                        single_chrom_result = Result({chrom: CNVresult.cnv[chrom]})
                        
                        # Get significant regions for this chromosome
                        significant_regions = None
                        if chromosome_events:
                            chrom_events = [event for event in chromosome_events if event['chromosome'] == chrom]
                            if chrom_events:
                                significant_regions = {
                                    chrom: [
                                        {
                                            'start_pos': 0,  # We'll use the whole chromosome length
                                            'end_pos': len(CNVresult.cnv[chrom]) * cnv_dict["bin_width"],
                                            'type': event['type']
                                        }
                                        for event in chrom_events
                                    ]
                                }
                        
                        # Create chromosome plot
                        plot_results = create_CNV_plot_per_chromosome(single_chrom_result, cnv_dict, significant_regions)
                        # Store plot data for grid arrangement
                        for plot_chrom, plot_buffer in plot_results:
                            if plot_chrom == chrom:
                                plots_data.append((chrom, plot_buffer))
                                break

            # Arrange plots in a grid (2 columns)
            if plots_data:
                PLOTS_PER_ROW = 2
                plot_width = inch * 3.6  # Adjust plot size to fit 2 per row
                plot_height = inch * 2.0
                
                # Calculate number of rows needed
                num_plots = len(plots_data)
                num_rows = (num_plots + PLOTS_PER_ROW - 1) // PLOTS_PER_ROW
                
                # Create plots grid
                for row in range(num_rows):
                    # Create a list to hold this row's plots
                    row_plots = []
                    row_captions = []
                    
                    # Get plots for this row
                    start_idx = row * PLOTS_PER_ROW
                    end_idx = min(start_idx + PLOTS_PER_ROW, num_plots)
                    
                    # Add plots to this row
                    for i in range(start_idx, end_idx):
                        chrom, plot_buffer = plots_data[i]
                        row_plots.append(Image(plot_buffer, width=plot_width, height=plot_height))
                        row_captions.append(
                            Paragraph(
                                f"CNV in {chrom}",
                                ParagraphStyle(
                                    'PlotCaption',
                                    parent=self.styles.styles['Caption'],
                                    fontSize=7,  # Slightly larger caption font
                                    leading=9,
                                    spaceBefore=2,
                                    spaceAfter=6
                                )
                            )
                        )
                    
                    # Add spacer before each row except the first
                    if row > 0:
                        self.elements.append(Spacer(1, 12))
                    
                    # Create a table for this row of plots with their captions
                    # Calculate cell widths to distribute space evenly
                    available_width = inch * 7.5  # A4 width minus margins
                    cell_width = available_width / PLOTS_PER_ROW
                    
                    # Create table data with plots and captions
                    table_data = [row_plots, row_captions]
                    
                    # Create table style
                    style = TableStyle([
                        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
                        ('LEFTPADDING', (0, 0), (-1, -1), 4),  # Slightly more padding
                        ('RIGHTPADDING', (0, 0), (-1, -1), 4),
                        ('TOPPADDING', (0, 0), (-1, -1), 2),
                        ('BOTTOMPADDING', (0, 0), (-1, -1), 2),
                    ])
                    
                    # Create and add the table
                    plot_table = Table(table_data, colWidths=[cell_width] * len(row_plots))
                    plot_table.setStyle(style)
                    self.elements.append(plot_table)

        except Exception as e:
            logger.error("Error processing CNV section: %s", str(e), exc_info=True)
            self.elements.append(
                Paragraph(
                    "Error processing CNV analysis data",
                    self.styles.styles["Normal"],
                )
            ) 