"""
report3.py

This module contains the main class and functions for creating the PDF report,
combining the styling of report2.py with the content of report.py and the
header/footer from header_footer.py.
"""

from reportlab.pdfgen import canvas
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    PageBreak,
    Image,
    Spacer,
    Table,
    TableStyle,
    Frame,
    PageTemplate,
    BaseDocTemplate,
    HRFlowable,
)
from reportlab.lib.enums import TA_LEFT, TA_RIGHT, TA_CENTER, TA_JUSTIFY
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.pagesizes import A4, inch
from reportlab.graphics.shapes import Line, LineShape, Drawing
from reportlab.lib.colors import Color, HexColor, white
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from PIL import Image as PILImage
import os
import io
import re
import pickle
import logging
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path
import matplotlib

matplotlib.use("agg")
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import natsort

from robin import images, resources, fonts
from robin.__about__ import __version__
from robin.subpages.CNV_object import Result
from robin.subpages.Fusion_object import (
    _annotate_results,
    get_gene_network,
    _get_reads,
    STRAND,
)
from .plotting import (
    target_distribution_plot,
    create_CNV_plot,
    create_CNV_plot_per_chromosome,
    classification_plot,
    coverage_plot,
)
from .utils import convert_to_space_separated_string, split_text, get_target_outliers
from dna_features_viewer import GraphicFeature, GraphicRecord

VERSION = __version__

# Use the main logger configured in the main application
logger = logging.getLogger(__name__)


class RobinReport:
    """Main class for generating ROBIN PDF reports with consistent styling."""

    def __init__(self, filename, output):
        """Initialize the report generator with filename and output directory."""
        self.filename = filename
        self.output = output
        self.sample_id = os.path.basename(os.path.normpath(output))
        if filename.startswith("None"):
            final_folder = os.path.basename(os.path.normpath(output))
            self.filename = filename.replace("None", final_folder, 1)
            self.sample_id = final_folder

        self.fonts_dir = os.path.join(os.path.dirname(os.path.abspath(fonts.__file__)))
        self.register_fonts()

        # Initialize styles
        self.styles = getSampleStyleSheet()
        self.setup_styles()

        # Initialize document elements
        self.elements_summary = []
        self.elements = []
        self.end_of_report_elements = []

        # Load master data
        self.masterdf = self.load_master_data()
        self.centreID = self.get_centre_id()

        # Create document
        self.doc = self.create_document()

    def register_fonts(self):
        """Register custom fonts for the report."""
        font_files = {
            "FiraSans": "fira-sans-v16-latin-regular.ttf",
            "FiraSans-Medium": "fira-sans-v16-latin-500.ttf",
            "FiraSans-Bold": "fira-sans-v16-latin-700.ttf",
            "FiraMono": "fira-mono-v14-latin-regular.ttf",
        }

        # Default to Helvetica if Fira Sans is not available
        self.use_default_font = True

        for font_name, font_file in font_files.items():
            font_path = os.path.join(self.fonts_dir, font_file)
            if os.path.exists(font_path):
                try:
                    pdfmetrics.registerFont(TTFont(font_name, font_path))
                    if font_name == "FiraSans":
                        self.use_default_font = False
                except Exception as e:
                    logger.warning(f"Could not register font {font_name}: {e}")
            else:
                logger.warning(f"Font file not found: {font_path}")

        # If FiraSans couldn't be registered, use Helvetica
        if self.use_default_font:
            logger.info("Using Helvetica as fallback font")

    def setup_styles(self):
        """Set up document styles with modern design."""
        # Get a fresh stylesheet and override any italic styles with regular font
        self.styles = getSampleStyleSheet()
        for style_name in self.styles.byName:
            if "Italic" in style_name:
                self.styles[style_name].fontName = (
                    "FiraSans" if not self.use_default_font else "Helvetica"
                )

        # Color palette for consistent branding
        self.COLORS = {
            "primary": HexColor("#4F9153"),  # ROBIN theme green
            "secondary": HexColor("#367AB3"),  # Blue accent
            "text": HexColor("#37474f"),
            "success": HexColor("#2e7d32"),
            "warning": HexColor("#f57c00"),
            "error": HexColor("#c62828"),
            "muted": HexColor("#546e7a"),
            "border": HexColor("#e9ecef"),
            "background": HexColor("#f8f9fa"),
        }

        # Base fonts to use
        if self.use_default_font:
            base_font = "Helvetica"
            bold_font = "Helvetica-Bold"
            medium_font = "Helvetica-Bold"  # Helvetica doesn't have medium weight
            mono_font = "Courier"
        else:
            base_font = "FiraSans"
            bold_font = "FiraSans-Bold"
            medium_font = "FiraSans-Medium"
            mono_font = "FiraMono"

        # Update base styles
        for style_name in ["Title", "Heading1", "Heading2", "Normal"]:
            if style_name == "Title":
                self.styles[style_name].fontName = bold_font
                self.styles[style_name].fontSize = 18
                self.styles[style_name].spaceAfter = 14
                self.styles[style_name].spaceBefore = 7
                self.styles[style_name].textColor = self.COLORS["primary"]
                self.styles[style_name].alignment = 1
            elif style_name == "Heading1":
                self.styles[style_name].fontName = bold_font
                self.styles[style_name].fontSize = 14
                self.styles[style_name].spaceAfter = 7
                self.styles[style_name].spaceBefore = 4
                self.styles[style_name].textColor = self.COLORS["primary"]
            elif style_name == "Heading2":
                self.styles[style_name].fontName = medium_font
                self.styles[style_name].fontSize = 12
                self.styles[style_name].spaceAfter = 6
                self.styles[style_name].spaceBefore = 4
                self.styles[style_name].textColor = self.COLORS["secondary"]
            else:  # Normal
                self.styles[style_name].fontName = base_font
                self.styles[style_name].fontSize = 9
                self.styles[style_name].leading = 14
                self.styles[style_name].spaceBefore = 2
                self.styles[style_name].spaceAfter = 2
                self.styles[style_name].textColor = self.COLORS["text"]
                # Add support for HTML-like tags
                self.styles[style_name].allowWidows = 0
                self.styles[style_name].allowOrphans = 0
                self.styles[style_name].wordWrap = 'CJK'

        # Add custom styles
        custom_styles = {
            "Smaller": {
                "fontSize": 7,
                "leading": 11,
                "textColor": self.COLORS["text"],
                "fontName": base_font,
            },
            "Bold": {
                "fontSize": 10,
                "leading": 12,
                "textColor": self.COLORS["primary"],
                "fontName": bold_font,
            },
            "Emphasis": {
                "fontSize": 10,
                "leading": 12,
                "textColor": self.COLORS["secondary"],
                "fontName": medium_font,
            },
            "SummaryCard": {
                "fontSize": 11,
                "leading": 14,
                "textColor": self.COLORS["primary"],
                "backColor": self.COLORS["background"],
                "borderColor": self.COLORS["border"],
                "borderWidth": 1,
                "borderPadding": 8,
                "spaceBefore": 8,
                "spaceAfter": 8,
                "bulletIndent": 0,
                "leftIndent": 8,
                "rightIndent": 8,
                "fontName": medium_font,
                "allowWidows": 0,
                "allowOrphans": 0,
                "wordWrap": "CJK",
                "alignment": 0,  # Left alignment
            },
            "Metric": {
                "fontSize": 12,
                "leading": 16,
                "textColor": self.COLORS["secondary"],
                "alignment": 1,
                "spaceBefore": 4,
                "spaceAfter": 4,
                "fontName": bold_font,
            },
            "Caption": {
                "fontSize": 9,
                "leading": 11,
                "textColor": self.COLORS["muted"],
                "alignment": 1,
                "spaceBefore": 4,
                "spaceAfter": 12,
                "fontName": base_font,
            },
            "MonoText": {
                "fontSize": 9,
                "leading": 12,
                "fontName": mono_font,
                "textColor": self.COLORS["text"],
                "backColor": self.COLORS["background"],
            },
        }

        for style_name, style_props in custom_styles.items():
            style_props["parent"] = self.styles["Normal"]
            # Remove existing style if it exists
            if style_name in self.styles:
                del self.styles[style_name]
            self.styles.add(ParagraphStyle(name=style_name, **style_props))

        # Table style
        self.MODERN_TABLE_STYLE = TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), self.COLORS["background"]),
                ("TEXTCOLOR", (0, 0), (-1, 0), self.COLORS["primary"]),
                ("FONTNAME", (0, 0), (-1, 0), bold_font),  # Header uses bold
                ("FONTNAME", (0, 1), (-1, -1), base_font),  # Body uses regular
                ("FONTSIZE", (0, 0), (-1, 0), 10),
                ("TOPPADDING", (0, 0), (-1, 0), 12),
                ("BOTTOMPADDING", (0, 0), (-1, 0), 12),
                ("LINEBELOW", (0, 0), (-1, 0), 1, self.COLORS["border"]),
                ("BACKGROUND", (0, 1), (-1, -1), white),
                ("TEXTCOLOR", (0, 1), (-1, -1), self.COLORS["text"]),
                ("FONTSIZE", (0, 1), (-1, -1), 9),
                ("TOPPADDING", (0, 1), (-1, -1), 8),
                ("BOTTOMPADDING", (0, 1), (-1, -1), 8),
                ("LEFTPADDING", (0, 0), (-1, -1), 16),
                ("RIGHTPADDING", (0, 0), (-1, -1), 16),
                ("GRID", (0, 0), (-1, -1), 0.5, self.COLORS["border"]),
                ("ALIGN", (0, 0), (-1, -1), "LEFT"),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                (
                    "ROWBACKGROUNDS",
                    (0, 1),
                    (-1, -1),
                    [white, self.COLORS["background"]],
                ),
            ]
        )

    def create_document(self):
        """Create the PDF document with proper margins and header/footer."""
        return SimpleDocTemplate(
            self.filename,
            pagesize=A4,
            rightMargin=0.75 * inch,
            leftMargin=0.75 * inch,
            topMargin=1.35 * inch,
            bottomMargin=0.75 * inch,
        )

    def load_master_data(self):
        """Load the master data file if it exists."""
        master_path = os.path.join(self.output, "master.csv")
        if os.path.exists(master_path):
            return pd.read_csv(master_path, index_col=0, header=None)
        return None

    def get_centre_id(self):
        """Get the centre ID from master data."""
        try:
            return (
                self.masterdf.loc["centreID"][1] if self.masterdf is not None else None
            )
        except KeyError:
            return None

    def add_section_divider(self):
        """Add a section divider to the document."""
        self.elements.append(
            HRFlowable(
                width="100%",
                thickness=1,
                color=self.COLORS["border"],
                spaceBefore=12,
                spaceAfter=12,
            )
        )

    def add_section_header(self, title, level=1):
        """Add a section header with consistent styling."""
        self.elements.append(Spacer(1, 6))
        self.elements.append(Paragraph(title, self.styles[f"Heading{level}"]))
        self.elements.append(Spacer(1, 4))

    def add_summary_card(self, content):
        """Add a summary card with consistent styling."""
        self.elements_summary.append(Paragraph(content, self.styles["SummaryCard"]))
        self.elements_summary.append(Spacer(1, 8))

    def add_figure(self, img, caption=None, width_scale=0.95):
        """Add a figure with optional caption."""
        self.elements.append(Spacer(1, 8))
        width, height = A4
        img_width = width * width_scale
        self.elements.append(Image(img, width=img_width, height=img_width / 1.6))
        if caption:
            self.elements.append(Spacer(1, 2))
            self.elements.append(Paragraph(caption, self.styles["Caption"]))
        self.elements.append(Spacer(1, 8))

    def create_auto_adjusting_table(self, data, style, max_width=None):
        """Create a table with automatically adjusted column widths."""
        if not data or not data[0]:
            return None

        if max_width is None:
            page_width, _ = A4
            max_width = page_width - inch

        num_cols = len(data[0])
        min_widths = [0] * num_cols

        for row in data:
            for i, cell in enumerate(row):
                if cell is not None:
                    cell_str = str(cell)
                    min_widths[i] = max(min_widths[i], len(cell_str) * 7)

        min_widths = [max(30, w) for w in min_widths]
        total_width = sum(min_widths)

        if total_width > max_width:
            scale = max_width / total_width
            col_widths = [w * scale for w in min_widths]
        else:
            col_widths = min_widths

        table = Table(data, colWidths=col_widths)
        table.setStyle(style)
        return table

    def generate_report(self):
        """Generate the complete PDF report."""
        try:
            logger.info("Starting report generation")
            
            # Add classification results
            logger.info("Adding classification results")
            self.add_classification_results()

            # Add CNV analysis
            logger.info("Starting CNV analysis section")
            self.add_cnv_analysis()

            # Add fusion analysis
            logger.info("Starting fusion analysis section")
            self.add_fusion_analysis()

            # Add coverage analysis
            logger.info("Starting coverage analysis section")
            self.add_coverage_analysis()

            # Add MGMT analysis
            logger.info("Starting MGMT analysis section")
            self.add_mgmt_analysis()

            # Add run data summary
            logger.info("Starting run data summary section")
            self.add_run_data_summary()

            # Add disclaimer
            logger.info("Adding disclaimer")
            self.add_disclaimer()

            # Combine all elements
            logger.info("Combining elements for final PDF")
            final_elements = (
                self.elements_summary + self.elements + self.end_of_report_elements
            )
            logger.info(f"Total elements: {len(final_elements)}")
            logger.info(f"Summary elements: {len(self.elements_summary)}")
            logger.info(f"Main elements: {len(self.elements)}")
            logger.info(f"End elements: {len(self.end_of_report_elements)}")

            # Build the PDF
            from .header_footer import header_footer_canvas_factory

            self.doc.multiBuild(
                final_elements,
                canvasmaker=header_footer_canvas_factory(
                    self.sample_id, self.centreID, self.styles, self.fonts_dir
                ),
            )

            logger.info(f"PDF created: {self.filename}")
            return self.filename

        except Exception as e:
            logger.error(f"Error generating report: {e}")
            raise

    def add_classification_results(self):
        """Add classification results section to the report."""
        self.elements_summary.append(
            Paragraph("Classification Results", self.styles["Heading1"])
        )

        # Collect all classification results first
        classification_data = []
        for name, df_name in [
            ("Sturgeon", "sturgeon_scores.csv"),
            ("NanoDX", "nanodx_scores.csv"),
            ("PanNanoDX", "pannanodx_scores.csv"),
            ("Forest", "random_forest_scores.csv"),
        ]:
            if df_name.lower() in [f.lower() for f in os.listdir(self.output)]:
                file_path = self.find_case_insensitive_file(df_name, self.output)
                if file_path:
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
                        else "Medium confidence" if confidence_value >= 0.5 else "Low confidence"
                    )
                    confidence_color = (
                        self.COLORS["success"]
                        if confidence_value >= 0.75
                        else (
                            self.COLORS["warning"]
                            if confidence_value >= 0.5
                            else self.COLORS["error"]
                        )
                    )

                    card_content = (
                        f"<b>{name} Classification</b><br/>"
                        f"Result: {lastrow_plot_top.index[0]}<br/>"
                        f'<font color="{confidence_color.hexval()}">{confidence_value:.1%} - {confidence_text}</font>'
                    )

                    if "number_probes" in df_store.columns:
                        card_content += (
                            f'<br/>Features found: {int(df_store.iloc[-1]["number_probes"])}'
                        )
                    
                    classification_data.append(Paragraph(card_content, self.styles["Normal"]))

        # Arrange results in two columns if we have multiple results
        if len(classification_data) > 1:
            # Calculate rows needed (round up division)
            num_rows = (len(classification_data) + 1) // 2
            
            # Create table data with empty cells
            table_data = []
            for i in range(num_rows):
                row = []
                # First column
                if i * 2 < len(classification_data):
                    row.append(classification_data[i * 2])
                else:
                    row.append('')
                # Second column
                if i * 2 + 1 < len(classification_data):
                    row.append(classification_data[i * 2 + 1])
                else:
                    row.append('')
                table_data.append(row)

            # Create and style the table
            table_style = TableStyle([
                ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                ('LEFTPADDING', (0, 0), (-1, -1), 20),
                ('RIGHTPADDING', (0, 0), (-1, -1), 20),
                ('TOPPADDING', (0, 0), (-1, -1), 10),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 10),
            ])

            table = Table(table_data, colWidths=['50%', '50%'])
            table.setStyle(table_style)
            self.elements_summary.append(table)
        else:
            # If only one result, add it directly
            for content in classification_data:
                self.elements_summary.append(content)
                self.elements_summary.append(Spacer(1, 8))

    def process_classification_data(self, name, file_path):
        """Process classification data for a specific model."""
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
            else "Medium confidence" if confidence_value >= 0.5 else "Low confidence"
        )
        confidence_color = (
            self.COLORS["success"]
            if confidence_value >= 0.75
            else (
                self.COLORS["warning"]
                if confidence_value >= 0.5
                else self.COLORS["error"]
            )
        )

        card_content = (
            f"<b>{name} Classification</b><br/>"
            f"Result: {lastrow_plot_top.index[0]}<br/>"
            f'<font color="{confidence_color.hexval()}">{confidence_value:.1%} - {confidence_text}</font>'
        )

        if "number_probes" in df_store.columns:
            card_content += (
                f'<br/>Features found: {int(df_store.iloc[-1]["number_probes"])}'
            )

        self.add_summary_card(card_content)

        # Add classification plot
        img_buf = classification_plot(df_store, name, 0.05)
        self.elements.append(
            Paragraph(f"{name} Classification Timeline", self.styles["Heading2"])
        )
        self.add_figure(img_buf, f"Classification confidence over time for {name}")
        self.add_section_divider()

    @staticmethod
    def find_case_insensitive_file(target_name, search_path=None):
        """Find a file regardless of case sensitivity."""
        if search_path is None:
            search_path = os.getcwd()
        target_name_lower = target_name.lower()
        for path in Path(search_path).rglob("*"):
            if path.is_file() and path.name.lower() == target_name_lower:
                return str(path)
        return None

    def add_cnv_analysis(self):
        """Add CNV analysis section to the report."""
        logger.info("Starting CNV analysis section")
        
        # Helper function for plot dimensions
        def calculate_plot_dimensions(width, is_genome_wide=False):
            """Calculate plot dimensions based on page width"""
            if is_genome_wide:
                img_width = width * 0.85
                img_height = img_width / 4  # Fixed aspect ratio for genome-wide plot
            else:
                # Reduce individual plot width to 40% of page width to allow for spacing
                img_width = width * 0.40  # Changed from (width * 0.90) / 2
                img_height = img_width * 0.5  # Maintain aspect ratio for chromosome plots
            return img_width, img_height

        # Load CNV data and XYestimate
        XYestimate = "Unknown"  # Default value
        cnv_file = os.path.join(self.output, "CNV.npy")
        
        # Check for required files
        if not os.path.exists(cnv_file):
            logger.info("No CNV.npy file found in output directory")
            return
            
        # Load CNV data
        logger.info("Found CNV.npy file")
        CNVresult = np.load(cnv_file, allow_pickle="TRUE").item()
        CNVresult = Result(CNVresult)
        logger.info(f"Loaded CNVresult with attributes: {dir(CNVresult)}")
        
        cnv_dict = np.load(
            os.path.join(self.output, "CNV_dict.npy"), allow_pickle=True
        ).item()
        logger.info(f"Loaded CNV dictionary with keys: {cnv_dict.keys()}")
        
        # Load XY estimate if available
        if os.path.exists(os.path.join(self.output, "XYestimate.pkl")):
            with open(os.path.join(self.output, "XYestimate.pkl"), "rb") as file:
                XYestimate = pickle.load(file)
        logger.info(f"Loaded XY estimate: {XYestimate}")

        # Add genetic sex summary
        if XYestimate != "Unknown":
            logger.info("Adding genetic sex analysis")
            sex_color = (
                self.COLORS["primary"]
                if XYestimate == "XX"
                else (
                    self.COLORS["secondary"]
                    if XYestimate == "XY"
                    else self.COLORS["muted"]
                )
            )

            card_content = (
                f'<font color="{sex_color.hexval()}">Estimated Genetic Sex: {XYestimate}</font>'
            )
            self.add_summary_card(card_content)

        # Add CNV section
        logger.info("Adding CNV section")
        
        # Load required resource files
        gene_bed_file = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)), 
            "unique_genes.bed"
        )
        cytoband_file = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)), 
            "cytoBand.txt"
        )
        
        # Load gene and cytoband data
        gene_bed = None
        cytobands_bed = None
        if os.path.exists(gene_bed_file):
            logger.info(f"Loading gene bed data from {gene_bed_file}")
            gene_bed = pd.read_csv(
                gene_bed_file,
                sep='\t',
                names=['chrom', 'start_pos', 'end_pos', 'gene']
            )
        if os.path.exists(cytoband_file):
            logger.info(f"Loading cytoband data from {cytoband_file}")
            cytobands_bed = pd.read_csv(
                cytoband_file,
                sep='\t',
                names=['chrom', 'start_pos', 'end_pos', 'name', 'stain']
            )
            
        # Create CNVAnalysis object with minimal required attributes
        from robin.subpages.CNV_object import CNVAnalysis
        cnv_analyzer = CNVAnalysis(target_panel="rCNS2")  # Default panel
        # Skip BedTree initialization by setting the attributes directly
        cnv_analyzer.gene_bed = gene_bed
        cnv_analyzer.cytobands_bed = cytobands_bed
        cnv_analyzer.cnv_dict = cnv_dict
        cnv_analyzer.XYestimate = XYestimate
        # Skip BedTree initialization
        cnv_analyzer.NewBed = None
        
        # Load reference CNV data
        with open(
            os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)),
                "HG01280_control_new.pkl",
            ),
            "rb",
        ) as f:
            ref_cnv_dict = pickle.load(f)
        
        # Calculate CNV differences
        from robin.subpages.CNV_object import CNV_Difference, moving_average, iterate_bam_bin
        result3 = CNV_Difference()
        
        # Get reference CNV data with matching bin width
        r2_cnv, _, _, _ = iterate_bam_bin(
            None,
            1,
            60,
            ref_cnv_dict,
            int(logging.ERROR),
            bin_width=cnv_dict["bin_width"],
        )
        
        # Calculate normalized CNV values and chromosome-specific thresholds
        chromosome_stats = {}  # Store stats for each chromosome
        result3 = CNV_Difference()
        chromosome_means = {}  # Store mean CNV value for each chromosome
        
        # First pass: calculate normalized values and chromosome means
        for key in CNVresult.cnv.keys():
            if key != "chrM" and re.match(r"^chr(\d+|X|Y)$", key):
                if key in r2_cnv:
                    moving_avg_data1 = moving_average(CNVresult.cnv[key])
                    moving_avg_data2 = moving_average(r2_cnv[key])
                    # Pad arrays to same length if needed
                    if len(moving_avg_data1) != len(moving_avg_data2):
                        max_len = max(len(moving_avg_data1), len(moving_avg_data2))
                        if len(moving_avg_data1) < max_len:
                            moving_avg_data1 = np.pad(moving_avg_data1, (0, max_len - len(moving_avg_data1)))
                        if len(moving_avg_data2) < max_len:
                            moving_avg_data2 = np.pad(moving_avg_data2, (0, max_len - len(moving_avg_data2)))
                    # Calculate difference
                    normalized_cnv = moving_avg_data1 - moving_avg_data2
                    result3.cnv[key] = normalized_cnv
                    # Calculate mean for this chromosome
                    chromosome_means[key] = np.mean(normalized_cnv)
        
        # Calculate genome-wide baseline from autosome means only
        autosome_means = [v for k, v in chromosome_means.items() if k not in ["chrX", "chrY", "chrM"]]
        genome_baseline = np.mean(autosome_means)
        genome_std = np.std(autosome_means)
        
        # Second pass: analyze each chromosome
        chromosome_wide_events = []  # Store chromosome-wide CNV events
        for key in CNVresult.cnv.keys():
            if key != "chrM" and re.match(r"^chr(\d+|X|Y)$", key):
                if key in result3.cnv:
                    normalized_cnv = result3.cnv[key]
                    chr_mean = chromosome_means[key]
                    chr_std = np.std(normalized_cnv)
                    z_score = 0  # Default value
                    
                    # Different analysis based on chromosome type
                    if key == "chrX":
                        if XYestimate == "XY":  # Male
                            # For males, X chromosome should be at half dosage
                            if chr_mean > 0.1:  # Gain from male baseline
                                chromosome_wide_events.append([
                                    key,
                                    "Whole chromosome (male X)",
                                    chr_mean,
                                    'GAIN',
                                    len(normalized_cnv) * cnv_dict["bin_width"] / 1_000_000
                                ])
                            elif chr_mean < -0.3:  # Loss from male baseline
                                chromosome_wide_events.append([
                                    key,
                                    "Whole chromosome (male X)",
                                    chr_mean,
                                    'LOSS',
                                    len(normalized_cnv) * cnv_dict["bin_width"] / 1_000_000
                                ])
                        else:  # Female or Unknown
                            z_score = (chr_mean - genome_baseline) / genome_std if genome_std > 0 else 0
                            if z_score > 1.5:
                                chromosome_wide_events.append([
                                    key,
                                    "Whole chromosome (female X)",
                                    chr_mean,
                                    'GAIN',
                                    len(normalized_cnv) * cnv_dict["bin_width"] / 1_000_000
                                ])
                            elif z_score < -1.5:
                                chromosome_wide_events.append([
                                    key,
                                    "Whole chromosome (female X)",
                                    chr_mean,
                                    'LOSS',
                                    len(normalized_cnv) * cnv_dict["bin_width"] / 1_000_000
                                ])
                    elif key == "chrY":
                        if XYestimate == "XY":  # Male
                            if chr_mean > 0.5:  # Significant gain from male baseline
                                chromosome_wide_events.append([
                                    key,
                                    "Whole chromosome (male Y)",
                                    chr_mean,
                                    'GAIN',
                                    len(normalized_cnv) * cnv_dict["bin_width"] / 1_000_000
                                ])
                            elif chr_mean < -0.5:  # Significant loss from male baseline
                                chromosome_wide_events.append([
                                    key,
                                    "Whole chromosome (male Y)",
                                    chr_mean,
                                    'LOSS',
                                    len(normalized_cnv) * cnv_dict["bin_width"] / 1_000_000
                                ])
                        else:  # Female or Unknown
                            if chr_mean > -0.2:  # Presence of Y material in female
                                chromosome_wide_events.append([
                                    key,
                                    "Y chromosome material (female)",
                                    chr_mean,
                                    'GAIN',
                                    len(normalized_cnv) * cnv_dict["bin_width"] / 1_000_000
                                ])
                    else:  # Autosomes
                        z_score = (chr_mean - genome_baseline) / genome_std if genome_std > 0 else 0
                        if z_score > 1.5:
                            chromosome_wide_events.append([
                                key,
                                "Whole chromosome",
                                chr_mean,
                                'GAIN',
                                len(normalized_cnv) * cnv_dict["bin_width"] / 1_000_000
                            ])
                        elif z_score < -1.5:
                            chromosome_wide_events.append([
                                key,
                                "Whole chromosome",
                                chr_mean,
                                'LOSS',
                                len(normalized_cnv) * cnv_dict["bin_width"] / 1_000_000
                            ])
                    
                    # Store chromosome-specific stats
                    chromosome_stats[key] = {
                        'mean': chr_mean,
                        'std': chr_std,
                        'z_score': z_score,
                        'gain_threshold': chr_mean + (1.0 * chr_std),
                        'loss_threshold': chr_mean - (1.0 * chr_std)
                    }
        
        # Update CNVAnalysis object
        cnv_analyzer.result3 = result3
        cnv_analyzer.chromosome_stats = chromosome_stats
        
        # Add CNV summary header and threshold information
        self.elements_summary.append(
            Paragraph("Copy Number Variation", self.styles["Heading2"])
        )

        # Generate genome-wide CNV plot for summary
        img_buf = create_CNV_plot(CNVresult, cnv_dict)
        width, height = A4
        img_width, img_height = calculate_plot_dimensions(width, is_genome_wide=True)
        img = Image(img_buf, width=img_width, height=img_height)
        self.elements_summary.append(img)
        self.elements_summary.append(
            Paragraph(
                "Copy number variation across chromosomes",
                self.styles["Caption"],
            )
        )
        
        # Create threshold summary table
        threshold_data = [['Chromosome', 'Mean CNV', 'Z-score', 'Gain Threshold', 'Loss Threshold']]
        for chrom in natsort.natsorted(chromosome_stats.keys()):
            stats = chromosome_stats[chrom]
            threshold_data.append([
                chrom,
                f"{stats['mean']:.3f}",
                f"{stats['z_score']:.3f}",
                f"{stats['gain_threshold']:.3f}",
                f"{stats['loss_threshold']:.3f}"
            ])
        
        # Add genome baseline to the table
        threshold_data.append([
            'Genome',
            f"{genome_baseline:.3f}",
            'baseline',
            f"{genome_baseline + genome_std:.3f}",
            f"{genome_baseline - genome_std:.3f}"
        ])
        
        # Create and style the threshold table
        threshold_table = Table(threshold_data, colWidths=[60, 70, 70, 70, 70])
        threshold_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), self.COLORS["background"]),
            ('TEXTCOLOR', (0, 0), (-1, 0), self.COLORS["primary"]),
            ('FONTNAME', (0, 0), (-1, 0), 'FiraSans-Bold' if not self.use_default_font else 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 8),
            ('FONTSIZE', (0, 1), (-1, -1), 7),
            ('TOPPADDING', (0, 0), (-1, -1), 2),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 2),
            ('LEFTPADDING', (0, 0), (-1, -1), 4),
            ('RIGHTPADDING', (0, 0), (-1, -1), 4),
            ('ALIGN', (1, 0), (-1, -1), 'RIGHT'),
            ('GRID', (0, 0), (-1, -1), 0.5, self.COLORS["border"]),
            ('BACKGROUND', (-1, -1), (-1, -1), self.COLORS["background"]),  # Highlight global row
        ]))
        
        self.elements_summary.append(Paragraph("CNV Analysis Thresholds:", self.styles["Normal"]))
        self.elements_summary.append(Spacer(1, 4))
        self.elements_summary.append(threshold_table)
        self.elements_summary.append(Spacer(1, 8))
        
        # Combine chromosome-wide events with cytoband-level events
        cnv_data = []
        
        # Add chromosome-wide events first
        for chrom, region, mean_cnv, cnv_type, length_mb in chromosome_wide_events:
            cnv_data.append([
                chrom,
                region,
                f"{length_mb:.1f}Mb",
                mean_cnv,
                cnv_type
            ])
        
        # Define list of chromosomes to analyze
        chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        
        # Then add within-chromosome events
        for chrom in chromosomes:
            if chrom in CNVresult.cnv:
                logger.info(f"Processing chromosome {chrom}")
                try:
                    # Skip within-chromosome analysis if the whole chromosome is altered
                    if not any(event[0] == chrom for event in chromosome_wide_events):
                        # Get chromosome-specific thresholds
                        stats = chromosome_stats[chrom]
                        gain_threshold = stats['gain_threshold']
                        loss_threshold = stats['loss_threshold']
                        
                        # Analyze CNV data with chromosome-specific thresholds
                        cytoband_analysis = cnv_analyzer.analyze_cytoband_cnv(result3.cnv, chrom)
                        if not cytoband_analysis.empty:
                            for _, row in cytoband_analysis.iterrows():
                                if row['mean_cnv'] > gain_threshold:
                                    length_mb = (row['end_pos'] - row['start_pos']) / 1_000_000
                                    cnv_data.append([
                                        chrom,
                                        row['name'],
                                        f"{length_mb:.1f}Mb",
                                        row['mean_cnv'],
                                        'GAIN',
                                    ])
                                elif row['mean_cnv'] < loss_threshold:
                                    length_mb = (row['end_pos'] - row['start_pos']) / 1_000_000
                                    cnv_data.append([
                                        chrom,
                                        row['name'],
                                        f"{length_mb:.1f}Mb",
                                        row['mean_cnv'],
                                        'LOSS',
                                    ])
                except Exception as e:
                    logger.error(f"Error processing {chrom}: {str(e)}")
        
        if cnv_data:
            logger.info(f"Found {len(cnv_data)} significant changes")
            # Sort data by chromosome (natural sort) and mean CNV value
            cnv_data.sort(key=lambda x: (natsort.natsort_key(x[0]), -abs(float(x[3]))))
            
            # Create table header
            table_data = [['Chr', 'Region', 'Length', 'Mean CNV', 'Type']]
            
            # Add data rows with formatting
            for row in cnv_data:
                chrom, region, length_mb, mean_cnv, cnv_type = row
                # Remove chromosome from region if it starts with it to avoid repetition
                if region.startswith(chrom):
                    region = region.replace(f"{chrom} ", "")
                
                # Add arrow icons for gains and losses
                type_with_arrow = '↑ GAIN' if cnv_type == 'GAIN' else '↓ LOSS'
                
                table_data.append([
                    chrom,
                    region,
                    length_mb,
                    f"{float(mean_cnv):.2f}",
                    type_with_arrow
                ])
            
            # Create table style with colors for gains and losses
            cnv_table_style = TableStyle([
                # Header styling
                ('BACKGROUND', (0, 0), (-1, 0), self.COLORS["background"]),
                ('TEXTCOLOR', (0, 0), (-1, 0), self.COLORS["primary"]),
                ('FONTNAME', (0, 0), (-1, 0), 'FiraSans-Bold' if not self.use_default_font else 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, 0), 8),
                ('TOPPADDING', (0, 0), (-1, 0), 4),
                ('BOTTOMPADDING', (0, 0), (-1, 0), 4),
                # Body styling
                ('FONTNAME', (0, 1), (-1, -1), 'FiraSans' if not self.use_default_font else 'Helvetica'),
                ('FONTSIZE', (0, 1), (-1, -1), 7),
                ('TOPPADDING', (0, 1), (-1, -1), 2),
                ('BOTTOMPADDING', (0, 1), (-1, -1), 2),
                ('LEFTPADDING', (0, 0), (-1, -1), 4),
                ('RIGHTPADDING', (0, 0), (-1, -1), 4),
                ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                ('ALIGN', (2, 1), (3, -1), 'RIGHT'),
                ('GRID', (0, 0), (-1, -1), 0.5, self.COLORS["border"]),
                # Add background shading for type column
                ('BACKGROUND', (4, 1), (4, -1), self.COLORS["background"]),
            ])
            
            # Add conditional formatting for GAIN/LOSS
            for i in range(1, len(table_data)):
                if '↑' in table_data[i][4]:  # GAIN
                    cnv_table_style.add('TEXTCOLOR', (4, i), (4, i), self.COLORS["success"])
                    cnv_table_style.add('BACKGROUND', (4, i), (4, i), HexColor('#e8f5e9'))  # Light green
                    cnv_table_style.add('TEXTCOLOR', (3, i), (3, i), self.COLORS["success"])
                else:  # LOSS
                    cnv_table_style.add('TEXTCOLOR', (4, i), (4, i), self.COLORS["error"])
                    cnv_table_style.add('BACKGROUND', (4, i), (4, i), HexColor('#ffebee'))  # Light red
                    cnv_table_style.add('TEXTCOLOR', (3, i), (3, i), self.COLORS["error"])
            
            # Create and add the table with adjusted column widths
            table = Table(table_data, colWidths=[45, 160, 50, 50, 45])  # Adjusted widths for new column
            table.setStyle(cnv_table_style)
            
            # Add table to summary section with less spacing
            self.elements_summary.append(
                Paragraph("CNV Summary", self.styles["Heading3"])
            )
            self.elements_summary.append(Spacer(1, 4))  # Reduced from 12
            self.elements_summary.append(table)
            self.elements_summary.append(Spacer(1, 8))  # Reduced from 12
        else:
            logger.info("No significant changes found")
            summary_color = self.COLORS["muted"]
            card_content = (
                f'<font color="{summary_color.hexval()}">'
                "<b>CNV Summary</b><br/>"
                "No significant CNV changes detected across any chromosome"
                "</font>"
            )
            self.add_summary_card(card_content)
        
        # Add per-chromosome plots in main section
        self.elements.append(PageBreak())
        self.elements.append(Paragraph("Detailed CNV Analysis", self.styles["Heading2"]))
        
        try:
            logger.info("Generating CNV plots")
            
            # Generate per-chromosome plots
            chromosome_plots = create_CNV_plot_per_chromosome(CNVresult, cnv_dict)
            
            # Get list of chromosomes with significant changes
            significant_chromosomes = set(row[0] for row in cnv_data)  # Get unique chromosomes from cnv_data
            
            # Filter chromosome plots to only show those with changes
            filtered_plots = [(chrom, plot_buf) for chrom, plot_buf in chromosome_plots if chrom in significant_chromosomes]
            
            if filtered_plots:
                # Calculate dimensions for the grid
                width, height = A4
                num_plots = len(filtered_plots)
                num_cols = 2  # Two plots per row
                num_rows = (num_plots + num_cols - 1) // num_cols  # Ceiling division
                
                # Calculate individual plot dimensions
                plot_width, plot_height = calculate_plot_dimensions(width, is_genome_wide=False)
                
                # Create table data for the grid
                table_data = []
                row_data = []
                
                for idx, (chrom, plot_buf) in enumerate(filtered_plots):
                    # Get chromosome-specific stats
                    stats = chromosome_stats[chrom]
                    mean_cnv = stats['mean']
                    std_cnv = stats['std']
                    
                    # Create a container for the plot and its caption
                    plot_elements = []
                    
                    # Add the stats as a caption first
                    plot_elements.append(
                        Paragraph(
                            f"{chrom} (Mean: {mean_cnv:.2f}, SD: {std_cnv:.2f})",
                            ParagraphStyle(
                                'PlotCaption',
                                parent=self.styles['Caption'],
                                fontSize=8,
                                leading=10,
                                alignment=1,
                                spaceBefore=0,
                                spaceAfter=5
                            )
                        )
                    )
                    
                    # Add the plot image
                    plot_elements.append(Image(plot_buf, width=plot_width, height=plot_height))
                    
                    # Add to row
                    row_data.append(plot_elements)
                    
                    # Start new row after every 2 plots
                    if len(row_data) == num_cols or idx == len(filtered_plots) - 1:
                        # Pad the last row with empty cells if needed
                        while len(row_data) < num_cols:
                            row_data.append('')
                        table_data.append(row_data)
                        row_data = []
                
                # Create the table with the plots
                plot_table = Table(table_data)
                
                # Style the table with more spacing
                plot_table.setStyle(TableStyle([
                    ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                    ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
                    ('LEFTPADDING', (0, 0), (-1, -1), 40),  # Increased horizontal padding
                    ('RIGHTPADDING', (0, 0), (-1, -1), 40), # Increased horizontal padding
                    ('TOPPADDING', (0, 0), (-1, -1), 5),    # Keep reduced vertical padding
                    ('BOTTOMPADDING', (0, 0), (-1, -1), 5), # Keep reduced vertical padding
                    ('COLWIDTH', (0, 0), (-1, -1), plot_width + 80),  # Increased column width padding
                ]))
                
                # Add less space before the table
                self.elements.append(Spacer(1, 10))  # Reduced from 20 to 10
                
                # Add the table to the document
                self.elements.append(plot_table)
                
                # Add less space before the caption
                self.elements.append(Spacer(1, 5))  # Reduced from 10 to 5
                
                self.elements.append(
                    Paragraph(
                        "Detailed copy number variation analysis for chromosomes with significant changes",
                        self.styles["Caption"],
                    )
                )
            else:
                self.elements.append(
                    Paragraph(
                        "No chromosomes showed significant copy number changes",
                        self.styles["Normal"],
                    )
                )
        except Exception as e:
            logger.error(f"Error generating CNV plots: {str(e)}", exc_info=True)

    def add_fusion_analysis(self):
        """Add fusion analysis section to the report."""
        fusion_file = os.path.join(self.output, "fusion_candidates_master.csv")
        fusion_file_all = os.path.join(self.output, "fusion_candidates_all.csv")

        if not (os.path.exists(fusion_file) or os.path.exists(fusion_file_all)):
            return

        logger.info("Processing fusion gene data")

        # Load gene annotation data
        datafile = "rCNS2_data.csv.gz"
        gene_table_path = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)), datafile
        )

        gene_table = None
        if os.path.exists(gene_table_path):
            try:
                gene_table = pd.read_csv(gene_table_path)
                logger.info("Loaded gene annotation file")
            except Exception as e:
                logger.error(f"Could not load gene annotation file: {e}")

        try:
            significant_fusions = []
            significant_fusions_all = []
            total_supporting_reads = 0
            total_supporting_reads_all = 0

            # Process targeted fusions
            if os.path.exists(fusion_file):
                fusion_candidates = pd.read_csv(
                    fusion_file, dtype=str, header=None, skiprows=1, on_bad_lines="warn"
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
                        reads = result[goodpairs][result[goodpairs][3].isin(gene_group)]
                        supporting_reads = self.count_supporting_reads(reads)
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
                        supporting_reads = self.count_supporting_reads(reads)
                        if supporting_reads >= 3:
                            significant_fusions_all.append(
                                (gene_group, supporting_reads)
                            )
                            total_supporting_reads_all += supporting_reads

            # Add fusion summaries to report
            if significant_fusions:
                self.elements_summary.append(
                    Paragraph("Targeted Gene Fusions:", self.styles["Heading3"])
                )
                self.elements_summary.append(
                    Paragraph(
                        f"Total Significant Fusion Events (at least 3 reads): {len(significant_fusions)}<br/>"
                        f"Total Supporting Reads: {total_supporting_reads}",
                        self.styles["Normal"],
                    )
                )
                fusion_list = []
                significant_fusions.sort(key=lambda x: x[1], reverse=True)
                for gene_group, supporting_reads in significant_fusions:
                    fusion_list.append(
                        f"• {' - '.join(gene_group)} ({supporting_reads} supporting reads)"
                    )
                self.elements_summary.append(
                    Paragraph(
                        "<br/>".join(fusion_list),
                        self.styles["Normal"],
                    )
                )

            # Add fusion plots if gene table is available
            if gene_table is not None and significant_fusions:
                self.elements.append(
                    Paragraph("Targeted Gene Fusion Plots", self.styles["Heading2"])
                )
                for gene_group, supporting_reads in significant_fusions:
                    if len(gene_group) > 5:
                        self.elements.append(
                            Paragraph(
                                f"Gene Fusion: {' - '.join(gene_group)} ({supporting_reads} supporting reads) - "
                                "Plot skipped due to complexity",
                                self.styles["Normal"],
                            )
                        )
                        continue

                    reads = result[goodpairs][result[goodpairs][3].isin(gene_group)]
                    fig = self.create_fusion_plot(reads, gene_table)
                    img_buf = io.BytesIO()
                    fig.savefig(img_buf, format="png", dpi=150, bbox_inches=None)
                    plt.close(fig)
                    img_buf.seek(0)

                    self.elements.append(
                        Paragraph(
                            f"Gene Fusion: {' - '.join(gene_group)} ({supporting_reads} supporting reads)",
                            self.styles["Normal"],
                        )
                    )
                    self.elements.append(
                        Image(img_buf, width=5 * inch, height=2.5 * inch)
                    )
                    self.elements.append(Spacer(1, 12))

        except pd.errors.EmptyDataError:
            logger.warning("Fusion candidates file is empty")
        except Exception as e:
            logger.error(f"Error processing fusion candidates: {e}")
            raise

    def add_coverage_analysis(self):
        """Add coverage analysis section to the report."""
        if not os.path.exists(os.path.join(self.output, "coverage_main.csv")):
            return

        cov_df_main = pd.read_csv(os.path.join(self.output, "coverage_main.csv"))
        bedcov_df_main = pd.read_csv(os.path.join(self.output, "bed_coverage_main.csv"))
        target_coverage_df = pd.read_csv(
            os.path.join(self.output, "target_coverage.csv")
        )

        self.elements_summary.append(
            Paragraph("Coverage Summary", self.styles["Heading2"])
        )
        self.elements_summary.append(Spacer(1, 16))

        global_coverage = cov_df_main["covbases"].sum() / cov_df_main["endpos"].sum()
        target_coverage = bedcov_df_main["bases"].sum() / bedcov_df_main["length"].sum()

        self.elements_summary.append(
            Paragraph(
                f"Coverage Depths - Global Estimated Coverage: {global_coverage:.2f}x "
                f"Targets Estimated Coverage: {target_coverage:.2f}x",
                self.styles["Normal"],
            )
        )

        if target_coverage < 10:
            self.elements_summary.append(
                Paragraph(
                    "Target Coverage is below the recommended 10x threshold",
                    self.styles["Normal"],
                )
            )

        # Process outliers
        outliers = get_target_outliers(target_coverage_df)
        outliers["coverage"] = outliers["coverage"].apply(lambda x: round(x, 1))
        outliers = outliers.sort_values(by="coverage", ascending=False)
        threshold = target_coverage

        outliers_above_threshold = outliers[outliers["coverage"] > threshold].copy()
        outliers_below_threshold = outliers[outliers["coverage"] <= threshold].copy()

        outliers_above_threshold["name_with_coverage"] = outliers_above_threshold.apply(
            lambda row: f"{row['name']} ({row['coverage']})", axis=1
        )
        outliers_below_threshold["name_with_coverage"] = outliers_below_threshold.apply(
            lambda row: f"{row['name']} ({row['coverage']})", axis=1
        )

        gene_names = " - ".join(outliers_above_threshold["name_with_coverage"])
        self.elements_summary.append(Spacer(1, 6))
        self.elements_summary.append(
            Paragraph(
                f"Outlier genes by coverage (high): {gene_names}",
                self.styles["Smaller"],
            )
        )

        gene_names = " - ".join(outliers_below_threshold["name_with_coverage"])
        self.elements_summary.append(Spacer(1, 6))
        self.elements_summary.append(
            Paragraph(
                f"Outlier genes by coverage (low): {gene_names}", self.styles["Smaller"]
            )
        )

        # Add coverage plots
        self.elements.append(PageBreak())
        self.elements.append(Paragraph("Target Coverage", self.styles["Emphasis"]))

        img_buf = coverage_plot(cov_df_main)
        width, height = A4
        img = Image(img_buf, width=width * 0.95, height=width, kind="proportional")
        self.elements.append(img)
        self.elements.append(Spacer(1, 12))
        self.elements.append(
            Paragraph(
                "Coverage over individual targets on each chromosome. Outliers are annotated by gene name.",
                self.styles["Smaller"],
            )
        )

        img_buf = target_distribution_plot(target_coverage_df)
        img = Image(img_buf, width=width * 0.9, height=width, kind="proportional")
        self.elements.append(img)
        self.elements.append(Spacer(1, 12))
        self.elements.append(
            Paragraph(
                f"The following table identifies potential outliers differing significantly from the mean coverage of {target_coverage:.2f}x",
                self.styles["Smaller"],
            )
        )

        # Add outliers table
        outliers = get_target_outliers(target_coverage_df)
        outliers["coverage"] = outliers["coverage"].apply(lambda x: round(x, 1))
        outliers = outliers.sort_values(by="coverage", ascending=False)
        data = [outliers.columns.to_list()] + outliers.values.tolist()
        table = self.create_auto_adjusting_table(data, self.MODERN_TABLE_STYLE)
        self.elements.append(table)
        self.elements.append(Spacer(1, 12))

    def add_mgmt_analysis(self):
        """Add MGMT analysis section to the report."""
        # Find the latest MGMT results
        last_seen = 0
        mgmt_results = None
        for file in natsort.natsorted(os.listdir(self.output)):
            if file.endswith("_mgmt.csv"):
                count = int(file.split("_")[0])
                if count > last_seen:
                    mgmt_results = pd.read_csv(os.path.join(self.output, file))
                    plot_out = os.path.join(self.output, file.replace(".csv", ".png"))
                    last_seen = count

        if last_seen == 0 or mgmt_results is None:
            return

        try:
            # Add summary to elements_summary
            self.elements_summary.append(
                Paragraph("MGMT Promoter Methylation", self.styles["Heading3"])
            )

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

            self.elements_summary.append(Paragraph(summary_text, self.styles["Normal"]))

            # Add detailed MGMT section
            self.elements.append(PageBreak())
            self.elements.append(
                Paragraph("MGMT Promoter Methylation Analysis", self.styles["Heading2"])
            )

            if os.path.exists(plot_out):
                img = Image(plot_out, width=6 * inch, height=4 * inch)
                self.elements.append(img)
                self.elements.append(Spacer(1, 12))

            # Add detailed results table
            data = [["Metric", "Value"]]
            data.append(["Status", methylation_status])
            data.append(["Average Methylation (%)", f"{round(methylation_average, 3)}"])
            if "pred" in mgmt_results.columns:
                data.append(
                    ["Prediction Score", f"{round(mgmt_results['pred'].iloc[0], 3)}"]
                )

            table = self.create_auto_adjusting_table(data, self.MODERN_TABLE_STYLE)
            self.elements.append(table)
            self.elements.append(Spacer(1, 12))

        except Exception as e:
            logger.error(f"Error processing MGMT results: {e}")

    def add_run_data_summary(self):
        """Add run data summary section to the report."""
        if not (self.masterdf is not None and isinstance(self.masterdf, pd.DataFrame)):
            return

        self.end_of_report_elements.append(PageBreak())
        self.end_of_report_elements.append(
            Paragraph("Run Data Summary", self.styles["Heading2"])
        )

        masterdf_dict = eval(
            self.masterdf[self.masterdf.index == "samples"][1]["samples"]
        )[self.sample_id]

        # Sample Information
        self.end_of_report_elements.append(
            Paragraph("Sample Information", self.styles["Heading3"])
        )
        self.end_of_report_elements.append(
            Paragraph(
                f"Sample ID: {self.sample_id} • "
                f"Run Start: {self.format_timestamp(masterdf_dict['run_time'])} • "
                f"Target Panel: {' '.join(self.masterdf.loc[(self.masterdf.index == 'target_panel')][1].values)}",
                self.styles["Smaller"],
            )
        )
        self.end_of_report_elements.append(Spacer(1, 6))

        # Device Details
        self.end_of_report_elements.append(
            Paragraph("Device Details", self.styles["Heading3"])
        )
        self.end_of_report_elements.append(
            Paragraph(
                f"Sequencing Device: {convert_to_space_separated_string(masterdf_dict['devices'])} • "
                f"Flowcell ID: {convert_to_space_separated_string(masterdf_dict['flowcell_ids'])} • "
                f"Basecalling Model: {convert_to_space_separated_string(masterdf_dict['basecall_models'])}",
                self.styles["Smaller"],
            )
        )
        self.end_of_report_elements.append(Spacer(1, 6))

        # File Locations
        self.end_of_report_elements.append(
            Paragraph("File Locations", self.styles["Heading3"])
        )
        self.end_of_report_elements.append(
            Paragraph(
                f"Run: {' '.join(self.masterdf.loc[(self.masterdf.index == 'watchfolder')][1].values)}<br/>"
                f"Out: {' '.join(self.masterdf.loc[(self.masterdf.index == 'output')][1].values)}<br/>"
                f"Ref: {' '.join(self.masterdf.loc[(self.masterdf.index == 'reference')][1].values)}",
                self.styles["Smaller"],
            )
        )
        self.end_of_report_elements.append(Spacer(1, 6))

        # Sequencing Statistics
        try:
            file_counters = eval(
                self.masterdf[self.masterdf.index == "samples"][1]["samples"]
            )[self.sample_id]["file_counters"]

            self.end_of_report_elements.append(
                Paragraph("Sequencing Statistics", self.styles["Heading2"])
            )
            self.end_of_report_elements.append(
                Paragraph(
                    f"BAM Files: {self.format_number(file_counters.get('bam_passed', 0))} passed, "
                    f"{self.format_number(file_counters.get('bam_failed', 0))} failed<br/>"
                    f"Mapped Reads: {self.format_number(file_counters.get('mapped_count', 0))} total "
                    f"({self.format_number(file_counters.get('pass_mapped_count', 0))} passed, "
                    f"{self.format_number(file_counters.get('fail_mapped_count', 0))} failed)<br/>"
                    f"Unmapped Reads: {self.format_number(file_counters.get('unmapped_count', 0))} total "
                    f"({self.format_number(file_counters.get('pass_unmapped_count', 0))} passed, "
                    f"{self.format_number(file_counters.get('fail_unmapped_count', 0))} failed)<br/>"
                    f"Total Bases: {self.format_number(file_counters.get('bases_count', 0))} "
                    f"({self.format_number(file_counters.get('pass_bases_count', 0))} passed, "
                    f"{self.format_number(file_counters.get('fail_bases_count', 0))} failed)",
                    self.styles["Smaller"],
                )
            )

        except Exception as e:
            logger.info(f"Error parsing file counters: {e}")

    def add_disclaimer(self):
        """Add the research use disclaimer to the report."""
        self.end_of_report_elements.append(PageBreak())
        self.end_of_report_elements.append(
            Paragraph("Disclaimer", self.styles["Heading2"])
        )

        disclaimer_text = (
            "This report and the data contained within it are intended for research use only and "
            "should not be used for direct diagnostic purposes. The methylation-based classification "
            "and other analyses provided here may be considered by neuropathologists as supplementary "
            "information in the context of comprehensive diagnostic assessment, which should include "
            "clinical history, radiological findings, and complete histopathological and molecular evaluation. "
            "The final interpretation and diagnosis should always be made by qualified healthcare professionals "
            "based on all available information."
        )

        self.end_of_report_elements.append(
            Paragraph(disclaimer_text, self.styles["Normal"])
        )
        self.end_of_report_elements.append(Spacer(1, 12))

    @staticmethod
    def create_fusion_plot(reads: pd.DataFrame, gene_table: pd.DataFrame) -> plt.Figure:
        """Creates a fusion plot matching the interactive version."""
        # Process reads to get result format
        result = _get_reads(reads)

        # Calculate reasonable width based on data
        max_span = 0
        for _, data in result.iterrows():
            span = data["end"] - data["start"]
            max_span = max(max_span, span)

        # Limit the figure size based on data span
        width = min(6, max(4, np.log10(max_span / 1000)))  # Scale based on span in kb
        height = width / 2  # Maintain aspect ratio

        plt.figure(figsize=(width, height), dpi=150)

        num_plots = 2 * len(result)
        num_cols = len(result)
        num_rows = (num_plots + num_cols - 1) // num_cols

        # Use a font that supports arrows
        plt.rcParams["font.family"] = ["DejaVu Sans", "Arial", "sans-serif"]
        RIGHT_ARROW = "→"
        LEFT_ARROW = "←"

        plt.subplots_adjust(
            left=0.2, right=0.8, bottom=0.25, top=0.75, wspace=0.6, hspace=0.5
        )

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

                ax.set_xlabel("Position (Mb)", fontsize=8)
                ax.tick_params(axis="both", which="major", labelsize=8)

        return plt.gcf()

    @staticmethod
    def count_supporting_reads(reads_df: pd.DataFrame) -> int:
        """Count unique supporting reads for a fusion."""
        return reads_df[7].nunique()

    @staticmethod
    def format_number(n):
        """Format number with commas for readability."""
        try:
            return "{:,}".format(int(float(n)))
        except (ValueError, TypeError):
            return str(n)

    @staticmethod
    def format_timestamp(timestamp_str):
        """Convert timestamp string to readable format."""
        clean_ts = str(timestamp_str).strip("[]'\"")
        try:
            dt = datetime.strptime(clean_ts.split("+")[0], "%Y-%m-%dT%H:%M:%S.%f")
            return dt.strftime("%B %d, %Y %H:%M:%S")
        except Exception as e:
            logger.debug(f"Error formatting timestamp {timestamp_str}: {e}")
            return str(timestamp_str)



def create_pdf(filename, output):
    """Create a PDF report from ROBIN analysis results."""
    report = RobinReport(filename, output)
    return report.generate_report()















def main():
    """Command-line interface for creating PDF reports."""
    import click

    @click.command()
    @click.argument("filename", type=str)
    @click.argument("output", type=str)
    @click.option("--debug", is_flag=True, help="Enable debug logging")
    def cli(filename, output, debug):
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

        try:
            pdf_file = create_pdf(filename, output)
            click.echo(f"Successfully created PDF report: {pdf_file}")
        except Exception as e:
            click.echo(f"Error creating PDF report: {str(e)}", err=True)
            raise click.Abort()

    cli()


if __name__ == "__main__":
    main()
