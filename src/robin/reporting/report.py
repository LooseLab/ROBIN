"""
report.py

This module contains the main function for creating the PDF report.
"""

import os
import io
import pickle
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
from PIL import Image as PILImage
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import inch
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Image,
    Table,
    TableStyle,
    PageBreak,
    CondPageBreak,
)
from reportlab.lib import colors

from .fonts import register_fonts
from .header_footer import header_footer_canvas_factory
from .plotting import (
    target_distribution_plot,
    create_CNV_plot,
    create_CNV_plot_per_chromosome,
    classification_plot,
    coverage_plot,
)
from .utils import convert_to_space_separated_string, split_text, get_target_outliers

import matplotlib
from collections import Counter

matplotlib.use("agg")
from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.font_manager import FontProperties


# Import the Result class
from robin import fonts
from robin.subpages.CNV_object import Result
from robin.subpages.Fusion_object import (
    _annotate_results,
    get_gene_network,
    _get_reads,
    STRAND,
)
from robin import resources
from dna_features_viewer import GraphicFeature, GraphicRecord
import natsort

import logging
from datetime import datetime
import textwrap

# Use the main logger configured in the main application
logger = logging.getLogger(__name__)


# Use DejaVu Sans as the default font
font_properties = FontProperties(
    fname=os.path.join(
        os.path.dirname(os.path.abspath(fonts.__file__)),
        "fira-sans-v16-latin-regular.ttf",
    )
)  # Update the path to the font file

matplotlib.rcParams["font.family"] = font_properties.get_name()


def format_number(n):
    """Format number with commas for readability"""
    try:
        return "{:,}".format(int(float(n)))
    except (ValueError, TypeError):
        return str(n)

def format_timestamp(timestamp_str):
    """Convert timestamp string to readable format"""
    # Remove brackets and timezone if present
    clean_ts = str(timestamp_str).strip("[]'\"")
    try:
        # Parse the timestamp - adjust format string as needed
        dt = datetime.strptime(clean_ts.split('+')[0], '%Y-%m-%dT%H:%M:%S.%f')
        return dt.strftime('%B %d, %Y %H:%M:%S')
    except Exception as e:
        logger.debug(f"Error formatting timestamp {timestamp_str}: {e}")
        return str(timestamp_str)


def create_auto_adjusting_table(data, style, max_width=None):
    """
    Creates a table that automatically adjusts its column widths based on content and page width.
    
    Args:
        data (List[List]): Table data as a list of rows
        style (TableStyle): Style to apply to the table
        max_width (float, optional): Maximum width for the table. Defaults to A4 width - 1 inch margins.
    
    Returns:
        Table: A reportlab Table object with optimized column widths
    """
    if not data or not data[0]:  # Check if data is empty or first row is empty
        return None
        
    # Get page width if max_width not specified
    if max_width is None:
        page_width, _ = A4
        max_width = page_width - inch  # 0.5 inch margin on each side
    
    # Calculate minimum column widths based on content
    num_cols = len(data[0])
    min_widths = [0] * num_cols
    
    for row in data:
        for i, cell in enumerate(row):
            if cell is not None:  # Handle None values
                cell_str = str(cell)
                # Estimate width based on character count (approximate)
                min_widths[i] = max(min_widths[i], len(cell_str) * 7)  # 7 points per character
    
    # Ensure minimum width of 30 points per column
    min_widths = [max(30, w) for w in min_widths]
    
    # Calculate total width
    total_width = sum(min_widths)
    
    if total_width > max_width:
        # Scale columns proportionally
        scale = max_width / total_width
        col_widths = [w * scale for w in min_widths]
    else:
        col_widths = min_widths
    
    # Create table with calculated widths
    table = Table(data, colWidths=col_widths)
    table.setStyle(style)
    
    return table

# Modern table style definition
MODERN_TABLE_STYLE = TableStyle([
    # Header styling
    ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#F5F6FA')),
    ('TEXTCOLOR', (0, 0), (-1, 0), colors.HexColor('#2C3E50')),
    ('FONTNAME', (0, 0), (-1, 0), 'FiraSans-Bold'),
    ('FONTSIZE', (0, 0), (-1, 0), 8),
    ('TOPPADDING', (0, 0), (-1, 0), 6),
    ('BOTTOMPADDING', (0, 0), (-1, 0), 6),
    
    # Body styling
    ('BACKGROUND', (0, 1), (-1, -1), colors.white),
    ('TEXTCOLOR', (0, 1), (-1, -1), colors.HexColor('#2C3E50')),
    ('FONTNAME', (0, 1), (-1, -1), 'FiraSans'),
    ('FONTSIZE', (0, 1), (-1, -1), 7),
    ('TOPPADDING', (0, 1), (-1, -1), 4),
    ('BOTTOMPADDING', (0, 1), (-1, -1), 4),
    
    # Grid styling
    ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#E2E8F0')),
    ('LINEBELOW', (0, 0), (-1, 0), 1, colors.HexColor('#CBD5E1')),
    
    # Alignment
    ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
    ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
    
    # Alternating row colors
    ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F8FAFC')]),
])

def create_pdf(filename, output):
    """
    Creates a PDF report.

    Args:
        filename (str): The filename for the PDF report.
        output (str): The directory where the output files are located.

    Returns:
        str: The filename of the created PDF report.
    """
    if filename.startswith("None"):
        final_folder = os.path.basename(os.path.normpath(output))
        filename = filename.replace("None", final_folder, 1)
        sample_id = final_folder
    else:
        sample_id = os.path.basename(os.path.normpath(output))
    logger.info(f"Creating PDF {filename} in {output} for sample {sample_id}")

    fonts_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fonts")
    # images_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "images")
    register_fonts(fonts_dir)

    styles = getSampleStyleSheet()
    
    # Update all text styles for consistency with Apple HIG
    styles["Heading1"].fontSize = 14  # Apple's recommended large title size
    styles["Heading1"].spaceAfter = 5
    styles["Heading1"].spaceBefore = 7
    styles["Heading1"].fontName = "FiraSans-Bold"
    styles["Heading1"].textColor = colors.HexColor('#000000')  # Apple prefers true black for important text

    styles["Heading2"].fontSize = 12  # Apple's recommended title size
    styles["Heading2"].spaceAfter = 4
    styles["Heading2"].spaceBefore = 6
    styles["Heading2"].fontName = "FiraSans-Bold"
    styles["Heading2"].textColor = colors.HexColor('#000000')

    styles["BodyText"].fontSize = 10  # Apple's recommended body text size
    styles["BodyText"].fontName = "FiraSans"
    styles["BodyText"].textColor = colors.HexColor('#000000')

    # Add styles for summary results
    styles.add(
        ParagraphStyle(
            name='SummaryResult',
            parent=styles['BodyText'],
            fontSize=13,
            leading=16,
            fontName="FiraSans",
            textColor=colors.HexColor('#000000'),
            spaceAfter=8,
            bulletIndent=12,
            leftIndent=24
        )
    )

    # Add style for important metrics
    styles.add(
        ParagraphStyle(
            name='Metric',
            parent=styles['BodyText'],
            fontSize=15,
            leading=18,
            fontName="FiraSans-Bold",
            textColor=colors.HexColor('#000000'),
            spaceAfter=4
        )
    )

    # Update custom styles
    styles.add(
        ParagraphStyle(
            name='Caption',
            parent=styles['Normal'],
            fontSize=8,
            leading=10,
            fontName="FiraSans",
            textColor=colors.HexColor('#666666'),
            alignment=1,  # Center alignment
            spaceAfter=6
        )
    )

    styles.add(
        ParagraphStyle(
            name='Bold',
            parent=styles['Normal'],
            fontSize=9,
            leading=11,
            fontName="FiraSans-Bold",
            textColor=colors.HexColor('#2C3E50')
        )
    )

    styles.add(
        ParagraphStyle(
            name='Smaller',
            parent=styles['Normal'],
            fontSize=8,
            leading=10,
            fontName="FiraSans",
            textColor=colors.HexColor('#2C3E50')
        )
    )

    # Add Underline style
    styles.add(
        ParagraphStyle(
            name='Underline',
            parent=styles['Heading2'],  # Based on Heading2 style
            fontSize=12,
            fontName="FiraSans-Bold",
            textColor=colors.HexColor('#2C3E50'),
            spaceAfter=4,
            spaceBefore=8
        )
    )

    # Define a consistent modern table style
    MODERN_TABLE_STYLE = TableStyle([
        # Header styling
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#F5F6FA')),  # Light blue-grey header
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.HexColor('#2C3E50')),   # Dark blue-grey text
        ('FONTNAME', (0, 0), (-1, 0), 'FiraSans-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 8),
        ('TOPPADDING', (0, 0), (-1, 0), 6),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 6),
        
        # Body styling
        ('BACKGROUND', (0, 1), (-1, -1), colors.white),
        ('TEXTCOLOR', (0, 1), (-1, -1), colors.HexColor('#2C3E50')),
        ('FONTNAME', (0, 1), (-1, -1), 'FiraSans'),
        ('FONTSIZE', (0, 1), (-1, -1), 7),
        ('TOPPADDING', (0, 1), (-1, -1), 4),
        ('BOTTOMPADDING', (0, 1), (-1, -1), 4),
        
        # Grid styling
        ('GRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#E2E8F0')),  # Light grey grid
        ('LINEBELOW', (0, 0), (-1, 0), 1, colors.HexColor('#CBD5E1')), # Slightly darker line below header
        
        # Alignment
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        
        # Alternating row colors
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#F8FAFC')]),
    ])

    # Update document margins for more compact but safe layout
    doc = SimpleDocTemplate(
        filename,
        pagesize=A4,
        rightMargin=0.5*inch,
        leftMargin=0.5*inch,
        topMargin=1.0*inch,      # Increased to 1 inch
        bottomMargin=0.5*inch,   # Reduced to 0.5 inch
    )

    # Update image handling function for more compact but safe layout
    def add_figure(elements, img, caption=None, width_scale=0.8):
        """Helper function for consistent image formatting"""
        elements.append(Spacer(1, 16))  # Increased spacing
        width, height = A4
        img_width = width * width_scale  # More compact images
        elements.append(Image(img, width=img_width, height=img_width/1.6))
        if caption:
            elements.append(Spacer(1, 4))
            elements.append(Paragraph(caption, styles["Caption"]))
        elements.append(Spacer(1, 12))  # Added space after figure
        elements.append(CondPageBreak(inch * 1))  # Add conditional page break if not enough space

    elements_summary = []
    elements = []

    masterdf = (
        pd.read_csv(os.path.join(output, "master.csv"), index_col=0, header=None)
        if os.path.exists(os.path.join(output, "master.csv"))
        else None
    )
    
    print(masterdf)
    try:
        centreID = masterdf.loc['centreID'][1]
    except KeyError:
        centreID = None  # or some default value

    elements_summary.append(Paragraph("Classification Summary", styles["Heading1"]))
    elements_summary.append(
        Paragraph(f"Sample {sample_id} has the following classifications:", styles["BodyText"])
    )
    elements_summary.append(Spacer(1, 12))

    threshold = 0.05

    try:
        # Add classification plots and details
        for name, df_name in [
            ("Sturgeon", "sturgeon_scores.csv"),
            ("NanoDX", "nanoDX_scores.csv"),
            ("PanNanoDX", "pannanodx_scores.csv"),
            ("Forest", "random_forest_scores.csv"),
        ]:
            # List files in the directory and convert them to lowercase
            files_in_directory = [f.lower() for f in os.listdir(output)]

            # Check if the lowercase version of df_name exists in the directory
            if df_name.lower() in files_in_directory:
                def find_case_insensitive_file(target_name, search_path):
                    target_name_lower = target_name.lower()
                    for path in Path(search_path).rglob('*'):
                        if path.is_file() and path.name.lower() == target_name_lower:
                            return str(path)
                    return None
                file_path = find_case_insensitive_file(df_name, output) or os.path.join(output, df_name)

                # Read the CSV file
                df_store = pd.read_csv(file_path)

                df_store2 = df_store.drop(columns=["timestamp"])

                if "number_probes" in df_store2.columns:
                    lastrow = df_store2.iloc[-1].drop("number_probes")
                else:
                    lastrow = df_store2.iloc[-1]

                lastrow_plot = lastrow.sort_values(ascending=False).head(10)
                lastrow_plot_top = lastrow.sort_values(ascending=False).head(1)
                # print (f"{name} classification: {lastrow_plot_top.index[0]} - {lastrow_plot_top.values[0]:.2f}")
                elements_summary.append(
                    Paragraph(
                        f"{name} classification: <b>{lastrow_plot_top.index[0]} - {lastrow_plot_top.values[0]:.2f}</b>",
                        styles["BodyText"],  # Changed from "Bold" to "BodyText"
                    )
                )

                img_buf = classification_plot(df_store, name, 0.05)
                img_pil = PILImage.open(img_buf)
                width_img, height_img = img_pil.size
                width, height = A4
                height = (width * 0.95) / width_img * height_img
                img = Image(img_buf, width=width * 0.95, height=height, kind="proportional")
                elements.append(img)

                elements.append(Spacer(1, 6))
                df = lastrow_plot.reset_index()
                df.columns = ["Classification", "Score"]
                df["Score"] = df["Score"].apply(lambda x: round(x, 5))
                df_transposed = df.set_index("Classification").T.reset_index()
                df_transposed.columns = [split_text(col) for col in df_transposed.columns]
                df_transposed.columns.name = None
                data = [df_transposed.columns.to_list()] + df_transposed.values.tolist()
                table = create_auto_adjusting_table(data, MODERN_TABLE_STYLE)
                elements.append(table)
                elements.append(Spacer(1, 12))
                elements.append(PageBreak())
            else:
                elements.append(
                    Paragraph(f"No {name} Classification Available", styles["BodyText"])
                )
    except Exception as e:
        logger.error(f"Error processing classification plots: {e}")
        raise

    try:
        # Add CNV plots
        if os.path.exists(os.path.join(output, "CNV.npy")):
            CNVresult = np.load(os.path.join(output, "CNV.npy"), allow_pickle="TRUE").item()
            CNVresult = Result(CNVresult)
            cnv_dict = np.load(
                os.path.join(output, "CNV_dict.npy"), allow_pickle=True
            ).item()
            file = open(os.path.join(output, "XYestimate.pkl"), "rb")
            XYestimate = pickle.load(file)
            elements_summary.append(Paragraph("Estimated Genetic Sex", styles["Heading3"]))
            elements_summary.append(
                Paragraph(
                    f"{XYestimate}",
                    styles["BodyText"],
                )
            )

            cnv_summary = create_CNV_plot(CNVresult, cnv_dict)
            img = Image(cnv_summary, width=6 * inch, height=1.5 * inch)
            elements.append(img)
            elements.append(Spacer(1, 6))

            # Create CNV per chromosome plots in a 2-column grid
            cnv_plots = create_CNV_plot_per_chromosome(CNVresult, cnv_dict)
            
            # Calculate dimensions for 2-column layout
            width, height = A4
            col_width = (width * 0.95) / 2  # 95% of page width split into 2 columns
            plot_height = col_width / 4  # Maintain aspect ratio
            
            # Process plots two at a time
            for i in range(0, len(cnv_plots), 2):
                # Create a list to hold the current row's plots
                row_plots = []
                
                # Add plots for this row (either 1 or 2 plots)
                for j in range(2):
                    if i + j < len(cnv_plots):
                        contig, img_buf = cnv_plots[i + j]
                        row_plots.append((contig, Image(img_buf, width=col_width, height=plot_height)))
                
                # Create a table for this row of plots
                plot_data = [[plot[1] for plot in row_plots]]
                if len(plot_data[0]) < 2:  # If odd number of plots, add empty cell
                    plot_data[0].append('')
                
                table = Table(plot_data, colWidths=[col_width] * 2)
                elements.append(table)
                elements.append(Spacer(1, 6))

            if XYestimate != "Unknown":
                elements.append(
                    Paragraph(f"Estimated Genetic Sex: {XYestimate}", styles["Smaller"])
                )
            elements.append(
                Paragraph(f"Current Bin Width: {cnv_dict['bin_width']}", styles["Smaller"])
            )
            elements.append(
                Paragraph(
                    f"Current Variance: {round(cnv_dict['variance'], 3)}", styles["Smaller"]
                )
            )
            elements.append(Spacer(1, 12))
            elements.append(PageBreak())
    except Exception as e:
        logger.error(f"Error processing CNV plots: {e}")
        raise

    try:
        # Add fusion gene plots and summary
        fusion_file = os.path.join(output, "fusion_candidates_master.csv")
        fusion_file_all = os.path.join(output, "fusion_candidates_all.csv")
        
        if os.path.exists(fusion_file) or os.path.exists(fusion_file_all):
            logger.info("Processing fusion gene data")
            
            # Load gene annotation data
            datafile = "rCNS2_data.csv.gz"
            gene_table_path = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)),
                datafile
            )
            
            if os.path.exists(gene_table_path):
                try:
                    gene_table = pd.read_csv(gene_table_path)
                    logger.info("Loaded gene annotation file")
                except Exception as e:
                    logger.error(f"Could not load gene annotation file: {e}")
                    gene_table = None
            else:
                logger.error("Could not find gene annotation file")
                gene_table = None

            try:
                significant_fusions = []
                significant_fusions_all = []
                total_supporting_reads = 0
                total_supporting_reads_all = 0

                # Process targeted fusions
                if os.path.exists(fusion_file):
                    fusion_candidates = pd.read_csv(
                        fusion_file,
                        dtype=str,
                        header=None,
                        skiprows=1,
                        on_bad_lines='warn'
                    )
                    logger.info(f"Loaded targeted fusion candidates with shape: {fusion_candidates.shape}")
                    
                    result, goodpairs = _annotate_results(fusion_candidates)
                    logger.info(f"Processed targeted fusion results. Good pairs: {goodpairs.sum()}")
                    
                    if not result.empty:
                        gene_pairs = result[goodpairs].sort_values(by=7)["tag"].unique().tolist()
                        stripped_list = [item.replace(" ", "") for item in gene_pairs]
                        gene_pairs = [pair.split(",") for pair in stripped_list]
                        gene_groups = get_gene_network(gene_pairs)
                        
                        for gene_group in gene_groups:
                            reads = result[goodpairs][result[goodpairs][3].isin(gene_group)]
                            supporting_reads = count_supporting_reads(reads)
                            if supporting_reads >= 3:
                                significant_fusions.append((gene_group, supporting_reads))
                                total_supporting_reads += supporting_reads

                # Process genome-wide fusions
                if os.path.exists(fusion_file_all):
                    fusion_candidates_all = pd.read_csv(
                        fusion_file_all,
                        dtype=str,
                        header=None,
                        skiprows=1,
                        on_bad_lines='warn'
                    )
                    logger.info(f"Loaded genome-wide fusion candidates with shape: {fusion_candidates_all.shape}")
                    
                    result_all, goodpairs_all = _annotate_results(fusion_candidates_all)
                    logger.info(f"Processed genome-wide fusion results. Good pairs: {goodpairs_all.sum()}")
                    
                    if not result_all.empty:
                        gene_pairs_all = result_all[goodpairs_all].sort_values(by=7)["tag"].unique().tolist()
                        stripped_list_all = [item.replace(" ", "") for item in gene_pairs_all]
                        gene_pairs_all = [pair.split(",") for pair in stripped_list_all]
                        gene_groups_all = get_gene_network(gene_pairs_all)
                        
                        for gene_group in gene_groups_all:
                            reads = result_all[goodpairs_all][result_all[goodpairs_all][3].isin(gene_group)]
                            supporting_reads = count_supporting_reads(reads)
                            if supporting_reads >= 3:
                                significant_fusions_all.append((gene_group, supporting_reads))
                                total_supporting_reads_all += supporting_reads

                # Add fusion summaries to report
                #elements_summary.append(Paragraph("Fusion Summary", styles["Heading2"]))
                
                # Targeted fusions summary
                if significant_fusions:
                    elements_summary.append(Paragraph("Targeted Gene Fusions:", styles["Heading3"]))
                    elements_summary.append(
                        Paragraph(
                            f"Total Significant Fusion Events (at least 3 reads): {len(significant_fusions)}<br/>"
                            f"Total Supporting Reads: {total_supporting_reads}",
                            styles["BodyText"],
                        )
                    )
                    fusion_list = []
                    significant_fusions.sort(key=lambda x: x[1], reverse=True)
                    for gene_group, supporting_reads in significant_fusions:
                        fusion_list.append(
                            f"• {' - '.join(gene_group)} ({supporting_reads} supporting reads)"
                        )
                    elements_summary.append(
                        Paragraph(
                            "<br/>".join(fusion_list),
                            styles["BodyText"],
                        )
                    )

                # Genome-wide fusions summary
                if significant_fusions_all:
                    elements_summary.append(Paragraph("Genome-wide Gene Fusions:", styles["Heading3"]))
                    elements_summary.append(
                        Paragraph(
                            f"Total Significant Fusion Events (at least 3 reads): {len(significant_fusions_all)}<br/>"
                            f"Total Supporting Reads: {total_supporting_reads_all}",
                            styles["BodyText"],
                        )
                    )
                    fusion_list = []
                    significant_fusions_all.sort(key=lambda x: x[1], reverse=True)
                    for gene_group, supporting_reads in significant_fusions_all:
                        fusion_list.append(
                            f"• {' - '.join(gene_group)} ({supporting_reads} supporting reads)"
                        )
                    #elements_summary.append(
                    #    Paragraph(
                    #        "<br/>".join(fusion_list),
                    #   )
                    #        styles["BodyText"],
                    #)

                if not significant_fusions and not significant_fusions_all:
                    elements_summary.append(
                        Paragraph(
                            "No significant fusion events detected (minimum 3 supporting reads required)",
                            styles["BodyText"],
                        )
                    )

                elements_summary.append(Spacer(1, 12))

                # Add fusion plots
                if gene_table is not None:
                    # Targeted fusion plots
                    if significant_fusions:
                        elements.append(Paragraph("Targeted Gene Fusion Plots", styles["Heading2"]))
                        for gene_group, supporting_reads in significant_fusions:
                            # Skip if too many genes involved
                            if len(gene_group) > 5:
                                elements.append(
                                    Paragraph(
                                        f"Gene Fusion: {' - '.join(gene_group)} ({supporting_reads} supporting reads) - "
                                        "Plot skipped due to complexity",
                                        styles["BodyText"]
                                    )
                                )
                                continue
                                
                            reads = result[goodpairs][result[goodpairs][3].isin(gene_group)]
                            fig = create_fusion_plot(reads, gene_table)
                            img_buf = io.BytesIO()
                            fig.savefig(img_buf, format='png', dpi=150, bbox_inches=None)  # Changed to 150 DPI
                            plt.close(fig)
                            img_buf.seek(0)
                            
                            elements.append(
                                Paragraph(
                                    f"Gene Fusion: {' - '.join(gene_group)} ({supporting_reads} supporting reads)", 
                                    styles["BodyText"]
                                )
                            )
                            elements.append(
                                Image(img_buf, width=5*inch, height=2.5*inch)  # Reduced from 7x3.5 to 5x2.5
                            )
                            elements.append(Spacer(1, 12))

                    # Genome-wide fusion summary table
                    if significant_fusions_all:
                        elements.append(Paragraph("Genome-wide Gene Fusions", styles["Heading2"]))
                        
                        # Add summary statistics
                        total_pairs = len(significant_fusions_all)
                        total_genes = len(set([gene for gene_group, _ in significant_fusions_all for gene in gene_group]))
                        
                        summary_text = (
                            f"Summary of genome-wide fusion analysis:\n"
                            f"• Total fusion pairs detected: {total_pairs}\n"
                            f"• Total unique genes involved: {total_genes}\n"
                            f"• Minimum supporting reads threshold: 3"
                        )
                        
                        elements.append(Paragraph(summary_text, styles["BodyText"]))
                        elements.append(Spacer(1, 12))
                        
                        # Add note about complexity
                        elements.append(
                            Paragraph(
                                "Note: Due to the complexity of genome-wide fusion events, "
                                "please refer to the interactive viewer for detailed visualization "
                                "and analysis of specific fusion pairs.",
                                styles["Italic"]
                            )
                        )
                        
                        elements.append(Spacer(1, 12))

            except pd.errors.EmptyDataError:
                logger.warning("Fusion candidates file is empty")
            except Exception as e:
                logger.error(f"Error processing fusion candidates: {e}")
                raise
    except Exception as e:
        logger.error(f"Error processing fusion plots: {e}")
        raise

    try:
        # Add coverage plots
        if os.path.exists(os.path.join(output, "coverage_main.csv")):
            cov_df_main = pd.read_csv(os.path.join(output, "coverage_main.csv"))
            bedcov_df_main = pd.read_csv(os.path.join(output, "bed_coverage_main.csv"))
            target_coverage_df = pd.read_csv(os.path.join(output, "target_coverage.csv"))
            elements_summary.append(Paragraph("Coverage Summary", styles["Heading2"]))
            elements_summary.append(
                Paragraph(
                    f"Coverage Depths - Global Estimated Coverage: {(cov_df_main['covbases'].sum() / cov_df_main['endpos'].sum()):.2f}x Targets Estimated Coverage: {(bedcov_df_main['bases'].sum() / bedcov_df_main['length'].sum()):.2f}x",
                    styles["BodyText"],
                )
            )

            if bedcov_df_main["bases"].sum() / bedcov_df_main["length"].sum() < 10:
                elements_summary.append(
                    Paragraph(
                        "Target Coverage is below the recommended 10x threshold",
                        styles["BodyText"],
                    )
                )

            outliers = get_target_outliers(target_coverage_df)
            outliers["coverage"] = outliers["coverage"].apply(lambda x: round(x, 1))
            outliers = outliers.sort_values(by="coverage", ascending=False)
            threshold = bedcov_df_main["bases"].sum() / bedcov_df_main["length"].sum()
            outliers_above_threshold = outliers[outliers["coverage"] > threshold].copy()
            outliers_below_threshold = outliers[outliers["coverage"] <= threshold].copy()
            outliers_above_threshold["name_with_coverage"] = outliers_above_threshold.apply(
                lambda row: f"{row['name']} ({row['coverage']})", axis=1
            )
            outliers_below_threshold["name_with_coverage"] = outliers_below_threshold.apply(
                lambda row: f"{row['name']} ({row['coverage']})", axis=1
            )

            gene_names = " - ".join(outliers_above_threshold["name_with_coverage"])
            elements_summary.append(Spacer(1, 6))
            elements_summary.append(
                Paragraph(
                    f"Outlier genes by coverage (high): {gene_names}", styles["Smaller"]
                )
            )
            elements_summary.append(Spacer(1, 6))
            gene_names = " - ".join(outliers_below_threshold["name_with_coverage"])
            elements_summary.append(Spacer(1, 6))
            elements_summary.append(
                Paragraph(
                    f"Outlier genes by coverage (low): {gene_names}", styles["Smaller"]
                )
            )

            # Add page break before Target Coverage section
            elements.append(PageBreak())
            elements.append(Paragraph("Target Coverage", styles["Underline"]))
            img_buf = coverage_plot(cov_df_main)
            width, height = A4
            img = Image(img_buf, width=width * 0.95, height=width, kind="proportional")
            elements.append(img)
            elements.append(Spacer(1, 12))
            elements.append(
                Paragraph(
                    "Coverage over individual targets on each chromosome. Outliers are annotated by gene name.",
                    styles["Smaller"],
                )
            )
            img_buf = target_distribution_plot(target_coverage_df)
            width, height = A4
            img = Image(img_buf, width=width * 0.9, height=width, kind="proportional")
            elements.append(img)
            elements.append(Spacer(1, 12))
            elements.append(
                Paragraph(
                    f"The following table identifies potential outliers differing significantly from the mean coverage of {(bedcov_df_main['bases'].sum() / bedcov_df_main['length'].sum()):.2f}x",
                    styles["Smaller"],
                )
            )

            outliers = get_target_outliers(target_coverage_df)
            outliers["coverage"] = outliers["coverage"].apply(lambda x: round(x, 1))
            outliers = outliers.sort_values(by="coverage", ascending=False)
            data = [outliers.columns.to_list()] + outliers.values.tolist()
            table = create_auto_adjusting_table(data, MODERN_TABLE_STYLE)
            elements.append(table)
            elements.append(Spacer(1, 12))
        else:
            elements.append(Paragraph("No Coverage Data Available", styles["BodyText"]))
    except Exception as e:
        logger.error(f"Error processing coverage plots: {e}")
        raise

    try:
        # Add run data summary with more compact spacing
        if masterdf is not None and isinstance(masterdf, pd.DataFrame):
            elements_summary.append(Paragraph("Run Data Summary", styles["Heading1"]))
            
            masterdf_dict = eval(masterdf[masterdf.index == "samples"][1]["samples"])[sample_id]
            
            # Sample Information - combine sections with less spacing
            elements_summary.append(Paragraph("Sample Information", styles["Heading2"]))
            elements_summary.append(
                Paragraph(
                    f"Sample ID: {sample_id} • "
                    f"Run Start: {format_timestamp(masterdf_dict['run_time'])} • "
                    f"Target Panel: {' '.join(masterdf.loc[(masterdf.index == 'target_panel')][1].values)}",
                    styles["BodyText"],
                )
            )
            elements_summary.append(Spacer(1, 6))
            
            # Device Details
            elements_summary.append(Paragraph("Device Details", styles["Heading2"]))
            elements_summary.append(
                Paragraph(
                    f"Sequencing Device: {convert_to_space_separated_string(masterdf_dict['devices'])} • "
                    f"Flowcell ID: {convert_to_space_separated_string(masterdf_dict['flowcell_ids'])} • "
                    f"Basecalling Model: {convert_to_space_separated_string(masterdf_dict['basecall_models'])}",
                    styles["BodyText"],
                )
            )
            elements_summary.append(Spacer(1, 6))
            
            # File Locations
            elements_summary.append(Paragraph("File Locations", styles["Heading2"]))
            elements_summary.append(
                Paragraph(
                    f"Run: {' '.join(masterdf.loc[(masterdf.index == 'watchfolder')][1].values)}<br/>"
                    f"Out: {' '.join(masterdf.loc[(masterdf.index == 'output')][1].values)}<br/>"
                    f"Ref: {' '.join(masterdf.loc[(masterdf.index == 'reference')][1].values)}",
                    styles["BodyText"],
                )
            )
            elements_summary.append(Spacer(1, 6))
            
            # Sequencing Statistics
            try:
                file_counters = eval(masterdf[masterdf.index == "samples"][1]["samples"])[sample_id]["file_counters"]
                
                elements_summary.append(Paragraph("Sequencing Statistics", styles["Heading2"]))
                elements_summary.append(
                    Paragraph(
                        f"BAM Files: {format_number(file_counters.get('bam_passed', 0))} passed, "
                        f"{format_number(file_counters.get('bam_failed', 0))} failed<br/>"
                        f"Mapped Reads: {format_number(file_counters.get('mapped_count', 0))} total "
                        f"({format_number(file_counters.get('pass_mapped_count', 0))} passed, "
                        f"{format_number(file_counters.get('fail_mapped_count', 0))} failed)<br/>"
                        f"Unmapped Reads: {format_number(file_counters.get('unmapped_count', 0))} total "
                        f"({format_number(file_counters.get('pass_unmapped_count', 0))} passed, "
                        f"{format_number(file_counters.get('fail_unmapped_count', 0))} failed)<br/>"
                        f"Total Bases: {format_number(file_counters.get('bases_count', 0))} "
                        f"({format_number(file_counters.get('pass_bases_count', 0))} passed, "
                        f"{format_number(file_counters.get('fail_bases_count', 0))} failed)",
                        styles["BodyText"],
                    )
                )
                
            except Exception as e:
                logger.info(f"Error parsing file counters: {e}")

            elements_summary.append(Spacer(1, 12))  # Final spacing before next section

    except Exception as e:
        logger.error(f"Error processing run data summary: {e}")
        raise

    try:
        # Add MGMT Promoter Methylation
        last_seen = 0
        for file in natsort.natsorted(os.listdir(output)):
            if file.endswith("_mgmt.csv"):
                count = int(file.split("_")[0])
                if count > last_seen:
                    results = pd.read_csv(os.path.join(output, file))
                    plot_out = os.path.join(output, file.replace(".csv", ".png"))
                    last_seen = count

        if last_seen > 0:
            elements.append(PageBreak())  # Add page break here
            elements.append(Paragraph("MGMT Promoter Methylation", styles["Underline"]))
            image = Image(plot_out, 6 * inch, 4 * inch)
            elements.append(image)
            data_list = [results.columns.values.tolist()]
            rounded_values = [[round_floats(val) for val in row] for row in results.values.tolist()]
            data_list.extend(rounded_values)
            table = create_auto_adjusting_table(data_list, MODERN_TABLE_STYLE)
            elements.append(table)

    except Exception as e:
        logger.error(f"Error processing MGMT Promoter Methylation: {e}")
        raise

    try:
        final_elements = elements_summary + elements
        doc.multiBuild(
            final_elements,
            canvasmaker=header_footer_canvas_factory(sample_id, centreID,styles, fonts_dir),
        )
        logger.info(f"PDF created: {filename}")
    except Exception as e:
        logger.error(f"Error finalizing PDF: {e}")
        raise

    return filename


def create_fusion_plot(reads: pd.DataFrame, gene_table: pd.DataFrame) -> plt.Figure:
    """Creates a fusion plot matching the interactive version from Fusion_object.py"""
    
    # Process reads to get result format
    result = _get_reads(reads)
    
    # Calculate reasonable width based on data
    max_span = 0
    for _, data in result.iterrows():
        span = data["end"] - data["start"]
        max_span = max(max_span, span)
    
    # Limit the figure size based on data span
    # Use log scale to prevent extremely wide plots
    width = min(6, max(4, np.log10(max_span/1000)))  # Scale based on span in kb
    height = width / 2  # Maintain aspect ratio
    
    # Create figure with moderate size and DPI
    plt.figure(figsize=(width, height), dpi=150)  # Changed to 150 DPI
    
    num_plots = 2 * len(result)
    num_cols = len(result)
    num_rows = (num_plots + num_cols - 1) // num_cols
    
    # Use a font that supports arrows
    plt.rcParams['font.family'] = ['DejaVu Sans', 'Arial', 'sans-serif']
    RIGHT_ARROW = '→'
    LEFT_ARROW = '←'
    
    # Create subplots with generous spacing
    plt.subplots_adjust(
        left=0.2,     # More left margin
        right=0.8,    # Less right margin
        bottom=0.25,  # More bottom margin
        top=0.75,     # Less top margin
        wspace=0.6,   # More space between plots
        hspace=0.5    # More vertical space
    )
    
    # Create subplots for each gene
    for i, ax in enumerate(range(num_plots), start=1):
        plt.subplot(num_rows, num_cols, i)
        row, col = divmod(i - 1, num_cols)
        data = result.iloc[col]
        
        chrom = data["chromosome"]
        start = data["start"]
        end = data["end"]
        
        if row == 1:  # Bottom row - gene structure
            features = []
            # Add gene body
            for _, gene_row in gene_table[
                gene_table["Seqid"].eq(chrom) &
                gene_table["Start"].le(end) &
                gene_table["End"].ge(start)
            ].iterrows():
                if gene_row["Type"] == "gene":
                    features.append(
                        GraphicFeature(
                            start=int(gene_row["Start"]),
                            end=int(gene_row["End"]),
                            strand=STRAND[gene_row["Strand"]],
                            thickness=8,
                            color="#ffd700",
                            label=gene_row["gene_name"],
                            fontdict={'family': 'DejaVu Sans', 'fontsize': 8}
                        )
                    )
            
            # Add exons
            for _, exon_row in gene_table[
                gene_table["gene_name"].eq(data["gene"]) &
                gene_table["Source"].eq("HAVANA") &
                gene_table["Type"].eq("exon")
            ].groupby(["Seqid", "Start", "End", "Type", "Strand"]).count().reset_index().iterrows():
                features.append(
                    GraphicFeature(
                        start=int(exon_row["Start"]),
                        end=int(exon_row["End"]),
                        strand=STRAND[exon_row["Strand"]],
                        thickness=4,
                        color="#C0C0C0",
                        arrow_style="simple"
                    )
                )
            
            record = GraphicRecord(
                sequence_length=end - start,
                first_index=start,
                features=features
            )
            ax = plt.gca()
            record.plot(ax=ax, with_ruler=False, draw_line=True, strand_in_label_threshold=4)
            
            # Adjust axis labels
            ax.set_xlabel("Position (Mb)", fontsize=8)
            ax.tick_params(axis='both', which='major', labelsize=8)
            
        else:  # Top row - read alignments
            features = []
            df = reads[reads["chromosome2"].eq(chrom)].sort_values(by="id")
            
            for _, row in df.iterrows():
                strand_marker = RIGHT_ARROW if row["strand"] == "+" else LEFT_ARROW
                features.append(
                    GraphicFeature(
                        start=int(row["start2"]),
                        end=int(row["end2"]),
                        strand=STRAND[row["strand"]],
                        color=row["color"] if "color" in row else "#ffd700",
                        label=strand_marker,
                        fontdict={'family': 'DejaVu Sans', 'fontsize': 8}
                    )
                )
            
            record = GraphicRecord(
                sequence_length=end - start,
                first_index=start,
                features=features
            )
            ax = plt.gca()
            record.plot(ax=ax, with_ruler=False, draw_line=True)
            
            # Adjust axis labels
            ax.set_xlabel("Position (Mb)", fontsize=8)
            ax.tick_params(axis='both', which='major', labelsize=8)
    
    # Don't use tight_layout - use the manual adjustments above instead
    # plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])
    
    return plt.gcf()


def count_supporting_reads(reads_df: pd.DataFrame) -> int:
    """
    Count unique supporting reads for a fusion.
    
    Args:
        reads_df (pd.DataFrame): DataFrame containing fusion reads
        
    Returns:
        int: Number of unique supporting reads
    """
    # Column 7 contains the read IDs
    return reads_df[7].nunique()

def round_floats(val):
    """Round float values to 3 decimal places, leave other types unchanged."""
    try:
        if isinstance(val, float):
            return round(val, 3)
        return val
    except:
        return val


import click

@click.command()
@click.argument('filename', type=str)
@click.argument('output', type=str)
@click.option('--debug', is_flag=True, help='Enable debug logging')
def main(filename: str, output: str, debug: bool):
    """
    Create a PDF report from ROBIN analysis results.

    Args:
        FILENAME: The filename for the PDF report
        OUTPUT: The directory containing the analysis output files
    """
    if debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    
    #try:
    pdf_file = create_pdf(filename, output)
    click.echo(f"Successfully created PDF report: {pdf_file}")
    #except Exception as e:
    #    click.echo(f"Error creating PDF report: {str(e)}", err=True)
    #    raise click.Abort()

if __name__ == '__main__':
    main()

