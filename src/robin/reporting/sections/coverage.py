"""
coverage.py

This module contains the coverage analysis section of the report.
"""

from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.platypus import Paragraph, Spacer, Table, TableStyle, PageBreak, Image
from reportlab.lib.styles import ParagraphStyle
from .base import ReportSection
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import io
import natsort


class CoverageSection(ReportSection):
    """Section containing the coverage analysis results."""

    def __init__(self, report):
        """Initialize the coverage section."""
        super().__init__(report)
        self._initialize_data()

    def _get_color_hex(self, color_name):
        """Convert a color from self.styles.COLORS to hex string."""
        color = self.styles.COLORS[color_name]
        if hasattr(color, 'hexval'):
            return color.hexval
        # Fallback colors if the style colors are not available
        fallback_colors = {
            "primary": "#2C3E50",
            "secondary": "#E74C3C",
            "error": "#C0392B",
            "muted": "#95A5A6"
        }
        return fallback_colors.get(color_name, "#000000")

    def _initialize_data(self):
        """Initialize coverage data from CSV files."""
        output_dir = self.report.output

        # Initialize coverage data
        self.coverage_data = {}
        self.chromosome_data = []
        self.target_data = []
        self.distribution_data = {}

        # Read coverage data if available
        if os.path.exists(os.path.join(output_dir, "coverage_main.csv")):
            self.cov_df_main = pd.read_csv(os.path.join(output_dir, "coverage_main.csv"))
            self.bedcov_df_main = pd.read_csv(os.path.join(output_dir, "bed_coverage_main.csv"))
            self.target_coverage_df = pd.read_csv(os.path.join(output_dir, "target_coverage.csv"))

            # Calculate global and target coverage
            global_coverage = self.cov_df_main['covbases'].sum() / self.cov_df_main['endpos'].sum()
            target_coverage = self.bedcov_df_main['bases'].sum() / self.bedcov_df_main['length'].sum()
            
            # Set coverage summary data
            self.coverage_data = {
                'global_coverage': global_coverage,
                'target_coverage': target_coverage
            }

            # Process chromosome-level data
            for _, row in self.cov_df_main.iterrows():
                if str(row['#rname']).startswith('chr'):  # Only process chromosome data
                    self.chromosome_data.append({
                        'name': row['#rname'],
                        'mean_coverage': row['meandepth'],
                        'covered_bases': row['covbases'],
                        'total_bases': row['endpos']
                    })

            # Process target data
            for _, row in self.target_coverage_df.iterrows():
                self.target_data.append({
                    'name': row['name'],
                    'chromosome': row['chrom'],
                    'start': row['startpos'],
                    'end': row['endpos'],
                    'coverage': row['coverage']
                })

            # Calculate distribution statistics
            coverages = self.target_coverage_df['coverage'].values
            self.distribution_data = {
                'median': np.median(coverages),
                'mean': np.mean(coverages),
                'min': np.min(coverages),
                'max': np.max(coverages),
                'above_30x': (coverages >= 30).sum() / len(coverages) * 100,
                'above_20x': (coverages >= 20).sum() / len(coverages) * 100,
                'above_10x': (coverages >= 10).sum() / len(coverages) * 100
            }

    def _create_chromosome_coverage_plot(self):
        """Create a plot showing coverage by chromosome."""
        # Filter and sort chromosomes
        pattern = r"^chr([0-9]+|X|Y)$"
        temp_covdf = self.cov_df_main[self.cov_df_main['#rname'].str.match(pattern)]
        sorteddf = temp_covdf.sort_values(
            by="#rname",
            key=lambda x: np.argsort(natsort.index_natsorted(temp_covdf["#rname"])),
        )
        sorteddf = sorteddf[sorteddf["#rname"] != "chrM"]

        # Create the plot
        plt.figure(figsize=(12, 6))
        plt.bar(sorteddf['#rname'], sorteddf['meandepth'], color="#2C3E50")
        plt.xlabel('Chromosome')
        plt.ylabel('Mean Coverage Depth')
        plt.title('Per Chromosome Coverage')
        plt.xticks(rotation=45)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()

        # Save plot to bytes buffer
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=300)
        plt.close()
        buf.seek(0)
        return buf

    def _create_target_coverage_plot(self):
        """Create a plot showing coverage distribution across targets."""
        # Calculate target statistics
        self.bedcov_df_main["length"] = self.bedcov_df_main["endpos"] - self.bedcov_df_main["startpos"] + 1
        grouped = (
            self.bedcov_df_main.groupby("chrom")
            .agg({"bases": "sum", "length": "sum"})
            .reset_index()
        )
        groupeddf = grouped.sort_values(
            by="chrom",
            key=lambda x: np.argsort(natsort.index_natsorted(grouped["chrom"])),
        )
        groupeddf = groupeddf[groupeddf["chrom"] != "chrM"]
        groupeddf["meandepth"] = groupeddf["bases"] / groupeddf["length"]

        # Filter and sort chromosomes for off-target data
        pattern = r"^chr([0-9]+|X|Y)$"
        temp_covdf = self.cov_df_main[self.cov_df_main['#rname'].str.match(pattern)]
        sorteddf = temp_covdf.sort_values(
            by="#rname",
            key=lambda x: np.argsort(natsort.index_natsorted(temp_covdf["#rname"])),
        )
        sorteddf = sorteddf[sorteddf["#rname"] != "chrM"]

        # Create the plot
        plt.figure(figsize=(12, 6))
        plt.scatter(sorteddf['#rname'], sorteddf['meandepth'], 
                   label='Off Target', color="#E74C3C", s=100)
        plt.scatter(groupeddf['chrom'], groupeddf['meandepth'], 
                   label='On Target', color="#2C3E50", s=100)
        plt.xlabel('Chromosome')
        plt.ylabel('Coverage Depth')
        plt.title('Target vs Off-Target Coverage by Chromosome')
        plt.xticks(rotation=45)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()
        plt.tight_layout()

        # Save plot to bytes buffer
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=300)
        plt.close()
        buf.seek(0)
        return buf

    def _create_target_boxplot(self):
        """Create a boxplot showing coverage distribution for targets."""
        # Prepare data
        self.target_coverage_df['coverage'] = self.target_coverage_df['coverage'].round(2)
        
        # Create figure
        plt.figure(figsize=(12, 6))
        
        # Create boxplot
        bp = plt.boxplot([self.target_coverage_df[self.target_coverage_df['chrom'] == chrom]['coverage'] 
                         for chrom in sorted(self.target_coverage_df['chrom'].unique(), key=natsort.natsort_keygen())],
                        patch_artist=True)
        
        # Style the boxplot
        plt.setp(bp['boxes'], facecolor="#2C3E50", alpha=0.6)
        plt.setp(bp['medians'], color="#E74C3C")
        plt.setp(bp['fliers'], marker='o', markerfacecolor="#C0392B")
        
        # Customize plot
        plt.xlabel('Chromosome')
        plt.ylabel('Coverage Depth')
        plt.title('Coverage Distribution by Chromosome')
        plt.xticks(range(1, len(self.target_coverage_df['chrom'].unique()) + 1),
                  sorted(self.target_coverage_df['chrom'].unique(), key=natsort.natsort_keygen()),
                  rotation=45)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.tight_layout()

        # Save plot to bytes buffer
        buf = io.BytesIO()
        plt.savefig(buf, format='png', dpi=300)
        plt.close()
        buf.seek(0)
        return buf

    def add_content(self):
        """Add the coverage analysis content to the report."""
        # Add Summary Section
        self.elements.append(
            Paragraph("Coverage Analysis", self.styles.styles["Heading2"])
        )
        self.elements.append(Spacer(1, 12))
        self.add_coverage_summary()
        
        # Add page break before detailed section
        self.elements.append(PageBreak())
        
        # Add Detailed Analysis Section
        self.elements.append(
            Paragraph("Detailed Coverage Analysis", self.styles.styles["Heading2"])
        )
        self.add_detailed_coverage()

        # Add summary to summary section
        self.summary_elements.append(
            Paragraph("Coverage Analysis", self.styles.styles["Heading3"])
        )
        self.summary_elements.append(Spacer(1, 6))
        
        if hasattr(self, 'coverage_data') and self.coverage_data:
            summary_text = []
            global_coverage = self.coverage_data.get('global_coverage', 0)
            target_coverage = self.coverage_data.get('target_coverage', 0)
            enrichment = target_coverage / global_coverage if global_coverage > 0 else 0
            
            summary_text.append(f"Target Coverage: {target_coverage:.2f}x")
            summary_text.append(f"Global Coverage: {global_coverage:.2f}x")
            summary_text.append(f"Enrichment: {enrichment:.2f}x")
            
            if hasattr(self, 'distribution_data') and self.distribution_data:
                summary_text.append(
                    f"Coverage ≥30x: {self.distribution_data.get('above_30x', 0):.1f}%"
                )
            
            self.summary_elements.append(
                Paragraph(
                    " | ".join(summary_text),
                    self.styles.styles["Normal"]
                )
            )
        else:
            self.summary_elements.append(
                Paragraph(
                    "No coverage data available",
                    self.styles.styles["Normal"]
                )
            )

    def add_coverage_summary(self):
        """Add the coverage summary section."""
        # Create styles for the summary text
        summary_style = ParagraphStyle(
            'SummaryStyle',
            parent=self.styles.styles['Normal'],
            spaceBefore=10,
            spaceAfter=10,
            leading=16
        )
        
        # Add overview text
        overview = Paragraph(
            "This section provides an analysis of sequencing coverage across the genome and targeted regions.",
            summary_style
        )
        self.elements.append(overview)
        self.elements.append(Spacer(1, 0.1 * inch))

        # Create summary table with key metrics
        if hasattr(self, 'coverage_data') and self.coverage_data is not None:
            global_coverage = self.coverage_data.get('global_coverage', 0)
            target_coverage = self.coverage_data.get('target_coverage', 0)
            enrichment = target_coverage / global_coverage if global_coverage > 0 else 0
            
            # Determine coverage quality level
            quality_level = self._get_coverage_quality(target_coverage)
            
            data = [
                ['Metric', 'Value', 'Status'],
                ['Global Coverage', f"{global_coverage:.2f}x", ''],
                ['Target Coverage', f"{target_coverage:.2f}x", quality_level],
                ['Enrichment Factor', f"{enrichment:.2f}x", '']
            ]
        else:
            data = [
                ['Metric', 'Value', 'Status'],
                ['Global Coverage', 'N/A', ''],
                ['Target Coverage', 'N/A', ''],
                ['Enrichment Factor', 'N/A', '']
            ]

        table_style = TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), self.styles.COLORS["primary"]),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 12),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
            ('BACKGROUND', (0, 1), (-1, -1), colors.white),
            ('TEXTCOLOR', (0, 1), (-1, -1), self.styles.COLORS["text"]),
            ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 1), (-1, -1), 10),
            ('GRID', (0, 0), (-1, -1), 0.5, self.styles.COLORS["border"]),
            ('ROWHEIGHTS', (0, 0), (-1, -1), 20),
        ])

        table = Table(data, colWidths=[2*inch, 2*inch, 2*inch])
        table.setStyle(table_style)
        self.elements.append(table)
        self.elements.append(Spacer(1, 0.2 * inch))

        # Add quality thresholds legend
        legend_style = ParagraphStyle(
            'LegendStyle',
            parent=self.styles.styles['Normal'],
            fontSize=8,
            textColor=self.styles.COLORS["muted"]
        )
        
        legend_text = """
        Coverage Quality Thresholds:
        • Excellent: ≥30x coverage
        • Good: ≥20x coverage
        • Moderate: ≥10x coverage
        • Insufficient: <10x coverage
        """
        legend = Paragraph(legend_text, legend_style)
        self.elements.append(legend)

    def add_detailed_coverage(self):
        """Add the detailed coverage analysis section."""
        if not hasattr(self, 'cov_df_main'):
            return

        # Add chromosome coverage plot
        title_style = ParagraphStyle(
            'TitleStyle',
            parent=self.styles.styles['Normal'],
            fontSize=12,
            spaceBefore=10,
            spaceAfter=10
        )
        
        # Add chromosome coverage plot
        self.elements.append(Paragraph("Chromosome Coverage Analysis", title_style))
        self.elements.append(Spacer(1, 0.1 * inch))
        chrom_plot = Image(self._create_chromosome_coverage_plot())
        chrom_plot.drawHeight = 4 * inch
        chrom_plot.drawWidth = 7 * inch
        self.elements.append(chrom_plot)
        self.elements.append(Spacer(1, 0.2 * inch))
        
        # Add target vs off-target plot
        self.elements.append(Paragraph("Target vs Off-Target Coverage", title_style))
        self.elements.append(Spacer(1, 0.1 * inch))
        target_plot = Image(self._create_target_coverage_plot())
        target_plot.drawHeight = 4 * inch
        target_plot.drawWidth = 7 * inch
        self.elements.append(target_plot)
        self.elements.append(Spacer(1, 0.2 * inch))
        
        # Add coverage distribution boxplot
        self.elements.append(Paragraph("Coverage Distribution by Chromosome", title_style))
        self.elements.append(Spacer(1, 0.1 * inch))
        box_plot = Image(self._create_target_boxplot())
        box_plot.drawHeight = 4 * inch
        box_plot.drawWidth = 7 * inch
        self.elements.append(box_plot)
        self.elements.append(Spacer(1, 0.2 * inch))
        
        # Add coverage distribution statistics
        self._add_coverage_distribution()

    def _add_coverage_distribution(self):
        """Add coverage distribution analysis."""
        if not hasattr(self, 'distribution_data') or self.distribution_data is None:
            return

        title_style = ParagraphStyle(
            'TitleStyle',
            parent=self.styles.styles['Normal'],
            fontSize=12,
            spaceBefore=10,
            spaceAfter=10
        )
        
        title = Paragraph("Coverage Distribution Analysis", title_style)
        self.elements.append(title)
        
        # Add coverage distribution statistics
        stats_style = ParagraphStyle(
            'StatsStyle',
            parent=self.styles.styles['Normal'],
            spaceBefore=6,
            spaceAfter=6
        )
        
        stats_text = f"""
        • Median Coverage: {self.distribution_data.get('median', 'N/A'):.2f}x
        • Mean Coverage: {self.distribution_data.get('mean', 'N/A'):.2f}x
        • Coverage Range: {self.distribution_data.get('min', 'N/A'):.2f}x - {self.distribution_data.get('max', 'N/A'):.2f}x
        • Regions ≥30x: {self.distribution_data.get('above_30x', 'N/A'):.1f}%
        • Regions ≥20x: {self.distribution_data.get('above_20x', 'N/A'):.1f}%
        • Regions ≥10x: {self.distribution_data.get('above_10x', 'N/A'):.1f}%
        """
        stats = Paragraph(stats_text, stats_style)
        self.elements.append(stats)

    def _get_coverage_quality(self, coverage):
        """Determine coverage quality level based on depth."""
        if coverage >= 30:
            return "Excellent"
        elif coverage >= 20:
            return "Good"
        elif coverage >= 10:
            return "Moderate"
        else:
            return "Insufficient"
