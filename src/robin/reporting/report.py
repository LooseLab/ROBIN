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
from reportlab.platypus import Frame, PageTemplate, BaseDocTemplate

from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Image,
    Table,
    TableStyle,
    PageBreak,
    CondPageBreak,
    HRFlowable,
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
        dt = datetime.strptime(clean_ts.split("+")[0], "%Y-%m-%dT%H:%M:%S.%f")
        return dt.strftime("%B %d, %Y %H:%M:%S")
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
                min_widths[i] = max(
                    min_widths[i], len(cell_str) * 7
                )  # 7 points per character

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


MODERN_TABLE_STYLE = TableStyle(
    [
        # Header styling
        (
            "BACKGROUND",
            (0, 0),
            (-1, 0),
            colors.HexColor("#f8f9fa"),
        ),  # Lighter background
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor("#1a237e")),
        ("FONTNAME", (0, 0), (-1, 0), "FiraSans"),
        ("FONTSIZE", (0, 0), (-1, 0), 9),
        ("TOPPADDING", (0, 0), (-1, 0), 8),
        ("BOTTOMPADDING", (0, 0), (-1, 0), 8),
        (
            "LINEBELOW",
            (0, 0),
            (-1, 0),
            1,
            colors.HexColor("#e9ecef"),
        ),  # Subtle header border
        # Body styling
        ("BACKGROUND", (0, 1), (-1, -1), colors.white),
        (
            "TEXTCOLOR",
            (0, 1),
            (-1, -1),
            colors.HexColor("#212529"),
        ),  # Darker for better contrast
        ("FONTNAME", (0, 1), (-1, -1), "FiraSans"),
        ("FONTSIZE", (0, 1), (-1, -1), 8),
        ("TOPPADDING", (0, 1), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 1), (-1, -1), 6),
        ("LEFTPADDING", (0, 0), (-1, -1), 12),  # Consistent padding
        ("RIGHTPADDING", (0, 0), (-1, -1), 12),
        # Grid styling - more subtle
        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#dee2e6")),
        # Alignment
        ("ALIGN", (0, 0), (-1, -1), "LEFT"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        # Alternating row colors - more subtle
        (
            "ROWBACKGROUNDS",
            (0, 1),
            (-1, -1),
            [colors.white, colors.HexColor("#f8f9fa")],
        ),
    ]
)


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

    # Color palette for consistent branding
    COLORS = {
        "primary": colors.HexColor("#1a237e"),
        "secondary": colors.HexColor("#283593"),
        "text": colors.HexColor("#37474f"),
        "success": colors.HexColor("#2e7d32"),
        "warning": colors.HexColor("#f57c00"),
        "error": colors.HexColor("#c62828"),
        "muted": colors.HexColor("#546e7a"),
        "border": colors.HexColor("#e9ecef"),
        "background": colors.HexColor("#f8f9fa"),
    }

    # Update all styles to use regular font weight
    for style_name in ["Title", "Heading1", "Heading2", "Normal"]:
        styles[style_name].fontName = "FiraSans"
        if style_name == "Title":
            styles[style_name].fontSize = 18
            styles[style_name].spaceAfter = 14
            styles[style_name].spaceBefore = 7
            styles[style_name].textColor = COLORS["primary"]
            styles[style_name].alignment = 1
        elif style_name == "Heading1":
            styles[style_name].fontSize = 14
            styles[style_name].spaceAfter = 7
            styles[style_name].spaceBefore = 4
            styles[style_name].textColor = COLORS["primary"]
        elif style_name == "Heading2":
            styles[style_name].fontSize = 12
            styles[style_name].spaceAfter = 6
            styles[style_name].spaceBefore = 4
            styles[style_name].textColor = COLORS["secondary"]
        else:  # Normal
            styles[style_name].fontSize = 9
            styles[style_name].leading = 14
            styles[style_name].spaceBefore = 2
            styles[style_name].spaceAfter = 2
            styles[style_name].textColor = COLORS["text"]

    # Add additional styles
    additional_styles = {
        "Smaller": {"fontSize": 7, "leading": 11, "textColor": COLORS["text"]},
        "Bold": {  # Keep name for compatibility
            "fontSize": 10,
            "leading": 12,
            "textColor": COLORS["primary"],
        },
        "Underline": {"fontSize": 10, "leading": 12, "textColor": COLORS["secondary"]},
        "SummaryCard": {
            "fontSize": 11,
            "leading": 14,
            "textColor": COLORS["primary"],
            "backColor": COLORS["background"],
            "borderColor": COLORS["border"],
            "borderWidth": 1,
            "borderPadding": 8,
            "spaceBefore": 8,
            "spaceAfter": 8,
            "bulletIndent": 0,
            "leftIndent": 8,
            "rightIndent": 8,
        },
        "Metric": {
            "fontSize": 12,
            "leading": 16,
            "textColor": COLORS["secondary"],
            "alignment": 1,
            "spaceBefore": 4,
            "spaceAfter": 4,
        },
        "Caption": {
            "fontSize": 9,
            "leading": 11,
            "textColor": COLORS["muted"],
            "alignment": 1,
            "spaceBefore": 4,
            "spaceAfter": 12,
        },
    }

    for style_name, style_props in additional_styles.items():
        style_props["parent"] = styles["Normal"]
        style_props["fontName"] = "FiraSans"
        styles.add(ParagraphStyle(name=style_name, **style_props))

    # Update table style to use regular font
    MODERN_TABLE_STYLE = TableStyle(
        [
            ("BACKGROUND", (0, 0), (-1, 0), COLORS["background"]),
            ("TEXTCOLOR", (0, 0), (-1, 0), COLORS["primary"]),
            ("FONTNAME", (0, 0), (-1, -1), "FiraSans"),
            ("FONTSIZE", (0, 0), (-1, 0), 10),
            ("TOPPADDING", (0, 0), (-1, 0), 12),
            ("BOTTOMPADDING", (0, 0), (-1, 0), 12),
            ("LINEBELOW", (0, 0), (-1, 0), 1, COLORS["border"]),
            ("BACKGROUND", (0, 1), (-1, -1), colors.white),
            ("TEXTCOLOR", (0, 1), (-1, -1), COLORS["text"]),
            ("FONTSIZE", (0, 1), (-1, -1), 9),
            ("TOPPADDING", (0, 1), (-1, -1), 8),
            ("BOTTOMPADDING", (0, 1), (-1, -1), 8),
            ("LEFTPADDING", (0, 0), (-1, -1), 16),
            ("RIGHTPADDING", (0, 0), (-1, -1), 16),
            ("GRID", (0, 0), (-1, -1), 0.5, COLORS["border"]),
            ("ALIGN", (0, 0), (-1, -1), "LEFT"),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, COLORS["background"]]),
        ]
    )

    # Helper function for section dividers
    def add_section_divider(elements):
        elements.append(
            HRFlowable(
                width="100%",
                thickness=1,
                color=COLORS["border"],
                spaceBefore=12,
                spaceAfter=12,
            )
        )

    # Helper function for summary cards - remove bold from title
    def add_summary_card(elements, title, content):
        # elements.append(Paragraph(title, styles["SummaryCard"]))  # Removed <b> tags
        elements.append(Paragraph(content, styles["Normal"]))
        elements.append(Spacer(1, 8))

    # Update document margins for better spacing
    doc = SimpleDocTemplate(
        filename,
        pagesize=A4,
        rightMargin=0.75 * inch,  # Increased margins for better readability
        leftMargin=0.75 * inch,
        topMargin=1.35 * inch,
        bottomMargin=0.75 * inch,
    )

    # Initialize lists to store document elements
    elements_summary = []
    elements = []

    # Load masterdf and get centreID first
    masterdf = (
        pd.read_csv(os.path.join(output, "master.csv"), index_col=0, header=None)
        if os.path.exists(os.path.join(output, "master.csv"))
        else None
    )

    try:
        centreID = masterdf.loc["centreID"][1] if masterdf is not None else None
    except KeyError:
        centreID = None

    def add_section_header(elements, title, level=1):
        """Helper function for consistent section headers"""
        elements.append(Spacer(1, 6))  # Reduced from 10
        elements.append(Paragraph(title, styles[f"Heading{level}"]))
        elements.append(Spacer(1, 4))  # Reduced from 6

    def add_subsection(elements, content, style="BodyText"):
        """Helper function for consistent subsection formatting"""
        elements.append(Spacer(1, 3))  # Reduced from 4
        elements.append(Paragraph(content, styles[style]))
        elements.append(Spacer(1, 3))  # Reduced from 4

    def add_figure(
        elements, img, caption=None, width_scale=0.95
    ):  # Increased width scale
        """Helper function for consistent image formatting"""
        elements.append(Spacer(1, 8))  # Reduced from 12
        width, height = A4
        img_width = width * width_scale
        elements.append(Image(img, width=img_width, height=img_width / 1.6))
        if caption:
            elements.append(Spacer(1, 2))  # Reduced from 4
            elements.append(Paragraph(caption, styles["Caption"]))
        elements.append(Spacer(1, 8))  # Reduced from 12

    # Initialize storage for classification data
    classification_data = []
    current_row = []

    # Load CNV data and XYestimate before using them
    XYestimate = "Unknown"  # Default value
    if os.path.exists(os.path.join(output, "CNV.npy")):
        CNVresult = np.load(os.path.join(output, "CNV.npy"), allow_pickle="TRUE").item()
        CNVresult = Result(CNVresult)
        cnv_dict = np.load(
            os.path.join(output, "CNV_dict.npy"), allow_pickle=True
        ).item()
        if os.path.exists(os.path.join(output, "XYestimate.pkl")):
            with open(os.path.join(output, "XYestimate.pkl"), "rb") as file:
                XYestimate = pickle.load(file)

    # Start with the report content
    # elements_summary.append(Paragraph("Analysis Summary", styles["Title"]))
    # add_section_divider(elements_summary)

    # Classification section
    elements_summary.append(Paragraph("Classification Results", styles["Heading1"]))
    # add_summary_card(
    #    elements_summary,
    #    "Sample Information",
    #    f"Sample ID: {sample_id}"
    # )

    # Add the function definition before the classification processing
    def find_case_insensitive_file(target_name, search_path):
        """Find a file regardless of case sensitivity."""
        target_name_lower = target_name.lower()
        for path in Path(search_path).rglob("*"):
            if path.is_file() and path.name.lower() == target_name_lower:
                return str(path)
        return None

    # Process classification data with new styling
    for name, df_name in [
        ("Sturgeon", "sturgeon_scores.csv"),
        ("NanoDX", "nanodx_scores.csv"),
        ("PanNanoDX", "pannanodx_scores.csv"),
        ("Forest", "random_forest_scores.csv"),
    ]:
        if df_name.lower() in [f.lower() for f in os.listdir(output)]:
            file_path = find_case_insensitive_file(df_name, output) or os.path.join(
                output, df_name
            )
            df_store = pd.read_csv(file_path)
            df_store2 = df_store.drop(columns=["timestamp"])
            lastrow = (
                df_store2.iloc[-1].drop("number_probes")
                if "number_probes" in df_store2.columns
                else df_store2.iloc[-1]
            )
            lastrow_plot_top = lastrow.sort_values(ascending=False).head(1)

            raw_confidence = float(lastrow_plot_top.values[0])
            confidence_value = (
                raw_confidence / 100.0 if name == "Forest" else raw_confidence
            )

            confidence_text = (
                "High confidence"
                if confidence_value >= 0.75
                else (
                    "Medium confidence" if confidence_value >= 0.5 else "Low confidence"
                )
            )
            confidence_color = (
                COLORS["success"]
                if confidence_value >= 0.75
                else COLORS["warning"] if confidence_value >= 0.5 else COLORS["error"]
            )

            # Create classification card
            card_content = (
                f"<b>{name} Classification</b><br/>"
                f"Result: {lastrow_plot_top.index[0]}<br/>"
                f'<font color="{confidence_color.hexval()}">{confidence_value:.1%} - {confidence_text}</font>'
            )
            if "number_probes" in df_store.columns:
                card_content += (
                    f'<br/>Features found: {int(df_store.iloc[-1]["number_probes"])}'
                )

            add_summary_card(elements_summary, name, card_content)

            # Add classification plot with improved styling
            img_buf = classification_plot(df_store, name, 0.05)
            img_pil = PILImage.open(img_buf)
            width_img, height_img = img_pil.size
            width, height = A4
            scaled_width = width * 0.85  # Slightly narrower for better presentation
            scaled_height = (scaled_width / width_img) * height_img

            elements.append(
                Paragraph(f"{name} Classification Timeline", styles["Heading2"])
            )
            elements.append(Image(img_buf, width=scaled_width, height=scaled_height))
            elements.append(
                Paragraph(
                    f"Classification confidence over time for {name}",
                    ParagraphStyle(
                        "Caption",
                        parent=styles["Normal"],
                        fontSize=9,
                        leading=11,
                        textColor=COLORS["muted"],
                        alignment=1,
                        spaceBefore=4,
                        spaceAfter=12,
                    ),
                )
            )
            add_section_divider(elements)

    # Genetic Sex Summary
    if XYestimate != "Unknown":
        elements_summary.append(Paragraph("Genetic Sex Analysis", styles["Heading2"]))

        sex_color = (
            COLORS["primary"]
            if XYestimate == "XX"
            else COLORS["secondary"] if XYestimate == "XY" else COLORS["muted"]
        )

        sex_icon = "♀" if XYestimate == "XX" else "♂" if XYestimate == "XY" else "?"

        add_summary_card(
            elements_summary,
            "Genetic Sex Estimate",
            f'<font size="16" color="{sex_color.hexval()}">{sex_icon}</font> '
            f'<font color="{sex_color.hexval()}">Estimated: {XYestimate}</font>',
        )

    # Add CNV summary if available
    if os.path.exists(os.path.join(output, "cnv.png")):
        elements_summary.append(Paragraph("Copy Number Variation", styles["Heading2"]))
        img = Image(os.path.join(output, "cnv.png"), width=width * 0.85)
        elements_summary.append(img)
        elements_summary.append(
            Paragraph(
                "Copy number variation across chromosomes",
                ParagraphStyle(
                    "Caption",
                    parent=styles["Normal"],
                    fontSize=9,
                    leading=11,
                    textColor=COLORS["muted"],
                    alignment=1,
                    spaceBefore=4,
                    spaceAfter=12,
                ),
            )
        )

    try:
        # Add fusion gene plots and summary
        fusion_file = os.path.join(output, "fusion_candidates_master.csv")
        fusion_file_all = os.path.join(output, "fusion_candidates_all.csv")

        if os.path.exists(fusion_file) or os.path.exists(fusion_file_all):
            logger.info("Processing fusion gene data")

            # Load gene annotation data
            datafile = "rCNS2_data.csv.gz"
            gene_table_path = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)), datafile
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
                        on_bad_lines="warn",
                    )
                    logger.info(
                        f"Loaded targeted fusion candidates with shape: {fusion_candidates.shape}"
                    )

                    result, goodpairs = _annotate_results(fusion_candidates)
                    logger.info(
                        f"Processed targeted fusion results. Good pairs: {goodpairs.sum()}"
                    )

                    if not result.empty:
                        gene_pairs = (
                            result[goodpairs].sort_values(by=7)["tag"].unique().tolist()
                        )
                        stripped_list = [item.replace(" ", "") for item in gene_pairs]
                        gene_pairs = [pair.split(",") for pair in stripped_list]
                        gene_groups = get_gene_network(gene_pairs)

                        for gene_group in gene_groups:
                            reads = result[goodpairs][
                                result[goodpairs][3].isin(gene_group)
                            ]
                            supporting_reads = count_supporting_reads(reads)
                            if supporting_reads >= 3:
                                significant_fusions.append(
                                    (gene_group, supporting_reads)
                                )
                                total_supporting_reads += supporting_reads

                # Process genome-wide fusions
                if os.path.exists(fusion_file_all):
                    fusion_candidates_all = pd.read_csv(
                        fusion_file_all,
                        dtype=str,
                        header=None,
                        skiprows=1,
                        on_bad_lines="warn",
                    )
                    logger.info(
                        f"Loaded genome-wide fusion candidates with shape: {fusion_candidates_all.shape}"
                    )

                    result_all, goodpairs_all = _annotate_results(fusion_candidates_all)
                    logger.info(
                        f"Processed genome-wide fusion results. Good pairs: {goodpairs_all.sum()}"
                    )

                    if not result_all.empty:
                        gene_pairs_all = (
                            result_all[goodpairs_all]
                            .sort_values(by=7)["tag"]
                            .unique()
                            .tolist()
                        )
                        stripped_list_all = [
                            item.replace(" ", "") for item in gene_pairs_all
                        ]
                        gene_pairs_all = [pair.split(",") for pair in stripped_list_all]
                        gene_groups_all = get_gene_network(gene_pairs_all)

                        for gene_group in gene_groups_all:
                            reads = result_all[goodpairs_all][
                                result_all[goodpairs_all][3].isin(gene_group)
                            ]
                            supporting_reads = count_supporting_reads(reads)
                            if supporting_reads >= 3:
                                significant_fusions_all.append(
                                    (gene_group, supporting_reads)
                                )
                                total_supporting_reads_all += supporting_reads

                # Add fusion summaries to report
                # elements_summary.append(Paragraph("Fusion Summary", styles["Heading2"]))

                # Targeted fusions summary
                if significant_fusions:
                    elements_summary.append(
                        Paragraph("Targeted Gene Fusions:", styles["Heading3"])
                    )
                    elements_summary.append(
                        Paragraph(
                            f"Total Significant Fusion Events (at least 3 reads): {len(significant_fusions)}<br/>"
                            f"Total Supporting Reads: {total_supporting_reads}",
                            styles["Normal"],
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
                            styles["Normal"],
                        )
                    )

                # Genome-wide fusions summary
                if significant_fusions_all:
                    elements_summary.append(
                        Paragraph("Genome-wide Gene Fusions:", styles["Heading3"])
                    )
                    elements_summary.append(
                        Paragraph(
                            f"Total Significant Fusion Events (at least 3 reads): {len(significant_fusions_all)}<br/>"
                            f"Total Supporting Reads: {total_supporting_reads_all}",
                            styles["Normal"],
                        )
                    )
                    fusion_list = []
                    significant_fusions_all.sort(key=lambda x: x[1], reverse=True)
                    for gene_group, supporting_reads in significant_fusions_all:
                        fusion_list.append(
                            f"• {' - '.join(gene_group)} ({supporting_reads} supporting reads)"
                        )
                    # elements_summary.append(
                    #    Paragraph(
                    #        "<br/>".join(fusion_list),
                    #   )
                    #        styles["Normal"],
                    # )

                if not significant_fusions and not significant_fusions_all:
                    elements_summary.append(
                        Paragraph(
                            "No significant fusion events detected (minimum 3 supporting reads required)",
                            styles["Normal"],
                        )
                    )

                # Add fusion plots
                if gene_table is not None:
                    # Targeted fusion plots
                    if significant_fusions:
                        elements.append(
                            Paragraph("Targeted Gene Fusion Plots", styles["Heading2"])
                        )
                        for gene_group, supporting_reads in significant_fusions:
                            # Skip if too many genes involved
                            if len(gene_group) > 5:
                                elements.append(
                                    Paragraph(
                                        f"Gene Fusion: {' - '.join(gene_group)} ({supporting_reads} supporting reads) - "
                                        "Plot skipped due to complexity",
                                        styles["Normal"],
                                    )
                                )
                                continue

                            reads = result[goodpairs][
                                result[goodpairs][3].isin(gene_group)
                            ]
                            fig = create_fusion_plot(reads, gene_table)
                            img_buf = io.BytesIO()
                            fig.savefig(
                                img_buf, format="png", dpi=150, bbox_inches=None
                            )  # Changed to 150 DPI
                            plt.close(fig)
                            img_buf.seek(0)

                            elements.append(
                                Paragraph(
                                    f"Gene Fusion: {' - '.join(gene_group)} ({supporting_reads} supporting reads)",
                                    styles["Normal"],
                                )
                            )
                            elements.append(
                                Image(
                                    img_buf, width=5 * inch, height=2.5 * inch
                                )  # Reduced from 7x3.5 to 5x2.5
                            )
                            elements.append(Spacer(1, 12))

                    # Genome-wide fusion summary table
                    if significant_fusions_all:
                        elements.append(
                            Paragraph("Genome-wide Gene Fusions", styles["Heading2"])
                        )

                        # Add summary statistics
                        total_pairs = len(significant_fusions_all)
                        total_genes = len(
                            set(
                                [
                                    gene
                                    for gene_group, _ in significant_fusions_all
                                    for gene in gene_group
                                ]
                            )
                        )

                        summary_text = (
                            f"Summary of genome-wide fusion analysis:\n"
                            f"• Total fusion pairs detected: {total_pairs}\n"
                            f"• Total unique genes involved: {total_genes}\n"
                            f"• Minimum supporting reads threshold: 3"
                        )

                        elements.append(Paragraph(summary_text, styles["Normal"]))
                        elements.append(Spacer(1, 12))

                        # Add note about complexity
                        elements.append(
                            Paragraph(
                                "Note: Due to the complexity of genome-wide fusion events, "
                                "please refer to the interactive viewer for detailed visualization "
                                "and analysis of specific fusion pairs.",
                                styles["Italic"],
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
            target_coverage_df = pd.read_csv(
                os.path.join(output, "target_coverage.csv")
            )
            elements_summary.append(Paragraph("Coverage Summary", styles["Heading2"]))
            elements_summary.append(Spacer(1, 16))
            elements_summary.append(
                Paragraph(
                    f"Coverage Depths - Global Estimated Coverage: {(cov_df_main['covbases'].sum() / cov_df_main['endpos'].sum()):.2f}x Targets Estimated Coverage: {(bedcov_df_main['bases'].sum() / bedcov_df_main['length'].sum()):.2f}x",
                    styles["Normal"],
                )
            )

            if bedcov_df_main["bases"].sum() / bedcov_df_main["length"].sum() < 10:
                elements_summary.append(
                    Paragraph(
                        "Target Coverage is below the recommended 10x threshold",
                        styles["Normal"],
                    )
                )

            outliers = get_target_outliers(target_coverage_df)
            outliers["coverage"] = outliers["coverage"].apply(lambda x: round(x, 1))
            outliers = outliers.sort_values(by="coverage", ascending=False)
            threshold = bedcov_df_main["bases"].sum() / bedcov_df_main["length"].sum()
            outliers_above_threshold = outliers[outliers["coverage"] > threshold].copy()
            outliers_below_threshold = outliers[
                outliers["coverage"] <= threshold
            ].copy()
            outliers_above_threshold["name_with_coverage"] = (
                outliers_above_threshold.apply(
                    lambda row: f"{row['name']} ({row['coverage']})", axis=1
                )
            )
            outliers_below_threshold["name_with_coverage"] = (
                outliers_below_threshold.apply(
                    lambda row: f"{row['name']} ({row['coverage']})", axis=1
                )
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
            elements.append(Paragraph("No Coverage Data Available", styles["Normal"]))

        # Add MGMT results
        last_seen = 0
        mgmt_results = None
        for file in natsort.natsorted(os.listdir(output)):
            if file.endswith("_mgmt.csv"):
                count = int(file.split("_")[0])
                if count > last_seen:
                    mgmt_results = pd.read_csv(os.path.join(output, file))
                    plot_out = os.path.join(output, file.replace(".csv", ".png"))
                    last_seen = count

        if last_seen > 0 and mgmt_results is not None:
            # Add summary to elements_summary
            elements_summary.append(
                Paragraph("MGMT Promoter Methylation", styles["Heading3"])
            )

            # Extract key metrics for summary
            try:
                methylation_status = (
                    mgmt_results["status"].iloc[0]
                    if "status" in mgmt_results.columns
                    else "Unknown"
                )
                methylation_average = (
                    mgmt_results["average"].iloc[0]
                    if "average" in mgmt_results.columns
                    else None
                )

                summary_text = f"Status: {methylation_status}"
                if methylation_average is not None:
                    summary_text += (
                        f" (Average methylation: {round(methylation_average, 3)}%)"
                    )

                elements_summary.append(
                    Paragraph(
                        summary_text,
                        styles["Normal"],
                    )
                )

                # Add detailed MGMT section to main report
                elements.append(PageBreak())
                elements.append(
                    Paragraph("MGMT Promoter Methylation Analysis", styles["Heading2"])
                )

                # Add the plot if it exists
                if os.path.exists(plot_out):
                    img = Image(plot_out, width=6 * inch, height=4 * inch)
                    elements.append(img)
                    elements.append(Spacer(1, 12))

                # Add detailed results table
                data = [["Metric", "Value"]]
                data.append(["Status", methylation_status])
                data.append(
                    ["Average Methylation (%)", f"{round(methylation_average, 3)}"]
                )
                if "pred" in mgmt_results.columns:
                    data.append(
                        [
                            "Prediction Score",
                            f"{round(mgmt_results['pred'].iloc[0], 3)}",
                        ]
                    )

                table = create_auto_adjusting_table(data, MODERN_TABLE_STYLE)
                elements.append(table)
                elements.append(Spacer(1, 12))

            except Exception as e:
                logger.error(f"Error processing MGMT results: {e}")

        # Create a separate list for the run data summary sections that will go at the end
        end_of_report_elements = []

        # Add run data summary with more compact spacing
        if masterdf is not None and isinstance(masterdf, pd.DataFrame):
            end_of_report_elements.append(PageBreak())  # Ensure it starts on a new page
            end_of_report_elements.append(
                Paragraph("Run Data Summary", styles["Heading2"])
            )

            masterdf_dict = eval(masterdf[masterdf.index == "samples"][1]["samples"])[
                sample_id
            ]

            # Sample Information - combine sections with less spacing
            end_of_report_elements.append(
                Paragraph("Sample Information", styles["Heading3"])
            )
            end_of_report_elements.append(
                Paragraph(
                    f"Sample ID: {sample_id} • "
                    f"Run Start: {format_timestamp(masterdf_dict['run_time'])} • "
                    f"Target Panel: {' '.join(masterdf.loc[(masterdf.index == 'target_panel')][1].values)}",
                    styles["Smaller"],
                )
            )
            end_of_report_elements.append(Spacer(1, 6))

            # Device Details
            end_of_report_elements.append(
                Paragraph("Device Details", styles["Heading3"])
            )
            end_of_report_elements.append(
                Paragraph(
                    f"Sequencing Device: {convert_to_space_separated_string(masterdf_dict['devices'])} • "
                    f"Flowcell ID: {convert_to_space_separated_string(masterdf_dict['flowcell_ids'])} • "
                    f"Basecalling Model: {convert_to_space_separated_string(masterdf_dict['basecall_models'])}",
                    styles["Smaller"],
                )
            )
            end_of_report_elements.append(Spacer(1, 6))

            # File Locations
            end_of_report_elements.append(
                Paragraph("File Locations", styles["Heading3"])
            )
            end_of_report_elements.append(
                Paragraph(
                    f"Run: {' '.join(masterdf.loc[(masterdf.index == 'watchfolder')][1].values)}<br/>"
                    f"Out: {' '.join(masterdf.loc[(masterdf.index == 'output')][1].values)}<br/>"
                    f"Ref: {' '.join(masterdf.loc[(masterdf.index == 'reference')][1].values)}",
                    styles["Smaller"],
                )
            )
            end_of_report_elements.append(Spacer(1, 6))

            # Sequencing Statistics
            try:
                file_counters = eval(
                    masterdf[masterdf.index == "samples"][1]["samples"]
                )[sample_id]["file_counters"]

                end_of_report_elements.append(
                    Paragraph("Sequencing Statistics", styles["Heading2"])
                )
                end_of_report_elements.append(
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
                        styles["Smaller"],
                    )
                )

            except Exception as e:
                logger.info(f"Error parsing file counters: {e}")

    except Exception as e:
        logger.error(f"Error processing run data summary: {e}")
        raise

    try:
        # Combine all elements with improved spacing
        final_elements = elements_summary + elements + end_of_report_elements

        # Build the PDF with the updated header/footer
        doc.multiBuild(
            final_elements,
            canvasmaker=header_footer_canvas_factory(
                sample_id, centreID, styles, fonts_dir
            ),
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
    width = min(6, max(4, np.log10(max_span / 1000)))  # Scale based on span in kb
    height = width / 2  # Maintain aspect ratio

    # Create figure with moderate size and DPI
    plt.figure(figsize=(width, height), dpi=150)  # Changed to 150 DPI

    num_plots = 2 * len(result)
    num_cols = len(result)
    num_rows = (num_plots + num_cols - 1) // num_cols

    # Use a font that supports arrows
    plt.rcParams["font.family"] = ["DejaVu Sans", "Arial", "sans-serif"]
    RIGHT_ARROW = "→"
    LEFT_ARROW = "←"

    # Create subplots with generous spacing
    plt.subplots_adjust(
        left=0.2,  # More left margin
        right=0.8,  # Less right margin
        bottom=0.25,  # More bottom margin
        top=0.75,  # Less top margin
        wspace=0.6,  # More space between plots
        hspace=0.5,  # More vertical space
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
                gene_table["Seqid"].eq(chrom)
                & gene_table["Start"].le(end)
                & gene_table["End"].ge(start)
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
                            fontdict={"family": "DejaVu Sans", "fontsize": 8},
                        )
                    )

            # Add exons
            for _, exon_row in (
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
                        start=int(exon_row["Start"]),
                        end=int(exon_row["End"]),
                        strand=STRAND[exon_row["Strand"]],
                        thickness=4,
                        color="#C0C0C0",
                        arrow_style="simple",
                    )
                )

            record = GraphicRecord(
                sequence_length=end - start, first_index=start, features=features
            )
            ax = plt.gca()
            record.plot(
                ax=ax, with_ruler=False, draw_line=True, strand_in_label_threshold=4
            )

            # Adjust axis labels
            ax.set_xlabel("Position (Mb)", fontsize=8)
            ax.tick_params(axis="both", which="major", labelsize=8)

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
                        fontdict={"family": "DejaVu Sans", "fontsize": 8},
                    )
                )

            record = GraphicRecord(
                sequence_length=end - start, first_index=start, features=features
            )
            ax = plt.gca()
            record.plot(ax=ax, with_ruler=False, draw_line=True)

            # Adjust axis labels
            ax.set_xlabel("Position (Mb)", fontsize=8)
            ax.tick_params(axis="both", which="major", labelsize=8)

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
@click.argument("filename", type=str)
@click.argument("output", type=str)
@click.option("--debug", is_flag=True, help="Enable debug logging")
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

    # try:
    pdf_file = create_pdf(filename, output)
    click.echo(f"Successfully created PDF report: {pdf_file}")
    # except Exception as e:
    #    click.echo(f"Error creating PDF report: {str(e)}", err=True)
    #    raise click.Abort()


if __name__ == "__main__":
    main()
