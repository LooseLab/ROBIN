"""
coverage.py

This module contains the coverage analysis section of the report.
"""

from reportlab.lib.units import inch
from reportlab.platypus import Paragraph, Spacer, Table, TableStyle, PageBreak, Image
from reportlab.lib.styles import ParagraphStyle
from .base import ReportSection
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Agg")
import io
import natsort
import logging

logger = logging.getLogger(__name__)

from robin.reference_contigs import is_visible_contig


class CoverageSection(ReportSection):
    """Section containing the coverage analysis results."""

    def __init__(self, report):
        """Initialize the coverage section."""
        super().__init__(report)
        self._initialize_data()

    def _get_color_hex(self, color_name):
        """Convert a color from self.styles.COLORS to hex string."""
        color = self.styles.COLORS[color_name]
        if hasattr(color, "hexval"):
            return color.hexval
        # Fallback colors if the style colors are not available
        fallback_colors = {
            "primary": "#2C3E50",
            "secondary": "#E74C3C",
            "error": "#C0392B",
            "muted": "#95A5A6",
        }
        return fallback_colors.get(color_name, "#000000")

    def _initialize_data(self):
        """Initialize coverage data from CSV files."""
        output_dir = self.report.output
        scope = getattr(self.report, "reference_contig_scope", None)
        # Initialize coverage data
        self.coverage_data = {}
        self.chromosome_data = []
        self.target_data = []
        self.distribution_data = {}

        # Read coverage data if available
        if os.path.exists(os.path.join(output_dir, "coverage_main.csv")):
            self.cov_df_main = pd.read_csv(
                os.path.join(output_dir, "coverage_main.csv")
            )
            self.bedcov_df_main = pd.read_csv(
                os.path.join(output_dir, "bed_coverage_main.csv")
            )
            self.target_coverage_df = pd.read_csv(
                os.path.join(output_dir, "target_coverage.csv")
            )

            # Apply the configured reference contig scope consistently across all
            # coverage plots and summary statistics.
            if "#rname" in self.cov_df_main.columns:
                self.cov_df_main = self.cov_df_main[
                    self.cov_df_main["#rname"]
                    .astype(str)
                    .map(lambda v: is_visible_contig(v, scope))
                ].copy()
            if "chrom" in self.bedcov_df_main.columns:
                self.bedcov_df_main = self.bedcov_df_main[
                    self.bedcov_df_main["chrom"]
                    .astype(str)
                    .map(lambda v: is_visible_contig(v, scope))
                ].copy()
            if "chrom" in self.target_coverage_df.columns:
                self.target_coverage_df = self.target_coverage_df[
                    self.target_coverage_df["chrom"]
                    .astype(str)
                    .map(lambda v: is_visible_contig(v, scope))
                ].copy()

            # Calculate global and target coverage
            global_coverage = (
                self.cov_df_main["covbases"].sum() / self.cov_df_main["endpos"].sum()
            )
            target_coverage = (
                self.bedcov_df_main["bases"].sum() / self.bedcov_df_main["length"].sum()
            )

            # Set coverage summary data
            self.coverage_data = {
                "global_coverage": global_coverage,
                "target_coverage": target_coverage,
            }

            # Process chromosome-level data
            for _, row in self.cov_df_main.iterrows():
                if is_visible_contig(str(row["#rname"]), scope):
                    self.chromosome_data.append(
                        {
                            "name": row["#rname"],
                            "mean_coverage": row["meandepth"],
                            "covered_bases": row["covbases"],
                            "total_bases": row["endpos"],
                        }
                    )

            # Process target data
            for _, row in self.target_coverage_df.iterrows():
                self.target_data.append(
                    {
                        "name": row["name"],
                        "chromosome": row["chrom"],
                        "start": row["startpos"],
                        "end": row["endpos"],
                        "coverage": row["coverage"],
                    }
                )

            # Calculate distribution statistics
            coverages = self.target_coverage_df["coverage"].values
            self.distribution_data = {
                "median": np.median(coverages),
                "mean": np.mean(coverages),
                "min": np.min(coverages),
                "max": np.max(coverages),
                "above_30x": (coverages >= 30).sum() / len(coverages) * 100,
                "above_20x": (coverages >= 20).sum() / len(coverages) * 100,
                "above_10x": (coverages >= 10).sum() / len(coverages) * 100,
            }

    def _create_chromosome_coverage_plot(self):
        """Create a bar plot showing coverage by chromosome."""
        try:
            # Check if data is available
            if not self.chromosome_data:
                plt.figure(figsize=(8, 4))
                plt.text(0.5, 0.5, "No chromosome coverage data available", 
                        ha='center', va='center', transform=plt.gca().transAxes,
                        fontsize=14, color='gray')
                plt.title("Per Chromosome Coverage")
                plt.axis('off')
            else:
                plt.figure(figsize=(8, 4))  # Reduced from default size
                chromosomes = [d["name"] for d in self.chromosome_data]
                coverages = [d["mean_coverage"] for d in self.chromosome_data]

                # Convert reportlab color to matplotlib color (hex string)
                plot_color = (
                    "#" + self.styles.COLORS["primary"].hexval()[2:]
                )  # Convert 0x... to #...
                plt.bar(chromosomes, coverages, color=plot_color)
                plt.title("Per Chromosome Coverage")
                plt.xlabel("Chromosome")
                plt.ylabel("Mean Coverage Depth")
                plt.grid(True, alpha=0.3)
                plt.xticks(rotation=45)
                plt.tight_layout()

            # Save plot to bytes buffer
            buf = io.BytesIO()
            plt.savefig(buf, format="png", dpi=300, bbox_inches='tight')
            plt.close()
            buf.seek(0)
            
            # Verify the buffer has data
            if buf.getvalue():
                return buf
            else:
                return self._create_empty_plot_buffer()
                
        except Exception as e:
            logger.warning(f"Error creating chromosome coverage plot: {e}")
            plt.close('all')  # Close any open figures
            return self._create_empty_plot_buffer()

    def _create_target_coverage_plot(self):
        """Create a plot showing coverage distribution across targets."""
        try:
            # Check if data is available
            if self.bedcov_df_main.empty or self.cov_df_main.empty:
                # Create empty plot with message
                plt.figure(figsize=(8, 4))
                plt.text(0.5, 0.5, "No coverage data available", 
                        ha='center', va='center', transform=plt.gca().transAxes,
                        fontsize=14, color='gray')
                plt.title("Target vs Off-Target Coverage by Chromosome")
                plt.axis('off')
            else:
                # Calculate target statistics
                self.bedcov_df_main["length"] = (
                    self.bedcov_df_main["endpos"] - self.bedcov_df_main["startpos"] + 1
                )
                grouped = (
                    self.bedcov_df_main.groupby("chrom")
                    .agg({"bases": "sum", "length": "sum"})
                    .reset_index()
                )
                groupeddf = grouped.sort_values(
                    by="chrom",
                    key=lambda x: np.argsort(natsort.index_natsorted(grouped["chrom"])),
                )
                groupeddf["meandepth"] = groupeddf["bases"] / groupeddf["length"]

                # Filter and sort chromosomes for off-target data
                sorteddf = self.cov_df_main.sort_values(
                    by="#rname",
                    key=lambda x: np.argsort(natsort.index_natsorted(x)),
                )

                # Create the plot
                plt.figure(figsize=(8, 4))  # Reduced size
                if not sorteddf.empty:
                    plt.scatter(
                        sorteddf["#rname"],
                        sorteddf["meandepth"],
                        label="Off Target",
                        color="#E74C3C",
                        s=50,
                    )  # Reduced marker size
                if not groupeddf.empty:
                    plt.scatter(
                        groupeddf["chrom"],
                        groupeddf["meandepth"],
                        label="On Target",
                        color="#2C3E50",
                        s=50,
                    )  # Reduced marker size
                plt.xlabel("Chromosome")
                plt.ylabel("Coverage Depth")
                plt.title("Target vs Off-Target Coverage by Chromosome")
                plt.xticks(rotation=45)
                plt.grid(True, alpha=0.3)
                plt.legend()
                plt.tight_layout()

            # Save plot to bytes buffer
            buf = io.BytesIO()
            plt.savefig(buf, format="png", dpi=300, bbox_inches='tight')
            plt.close()
            buf.seek(0)
            
            # Verify the buffer has data
            if buf.getvalue():
                return buf
            else:
                # Return a minimal valid PNG if buffer is empty
                return self._create_empty_plot_buffer()
                
        except Exception as e:
            logger.warning(f"Error creating target coverage plot: {e}")
            plt.close('all')  # Close any open figures
            return self._create_empty_plot_buffer()

    def _create_empty_plot_buffer(self):
        """Create a minimal valid PNG buffer for empty plots."""
        try:
            plt.figure(figsize=(8, 4))
            plt.text(0.5, 0.5, "No data available", 
                    ha='center', va='center', transform=plt.gca().transAxes,
                    fontsize=14, color='gray')
            plt.title("Coverage Plot")
            plt.axis('off')
            
            buf = io.BytesIO()
            plt.savefig(buf, format="png", dpi=300, bbox_inches='tight')
            plt.close()
            buf.seek(0)
            return buf
        except Exception:
            # If even this fails, return a minimal PNG
            import base64
            # Minimal 1x1 transparent PNG
            png_data = base64.b64decode(
                'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg=='
            )
            buf = io.BytesIO(png_data)
            buf.seek(0)
            return buf

    def _get_target_coverage_outliers(self):
        """Identify coverage outliers (gains/losses) per chromosome for table display."""
        if not hasattr(self, "target_coverage_df") or self.target_coverage_df.empty:
            return []
        df = self.target_coverage_df.copy()
        df["coverage"] = df["coverage"].round(2)
        sorted_chroms = sorted(
            df["chrom"].unique(), key=natsort.natsort_keygen()
        )
        rows = []
        for chrom in sorted_chroms:
            chrom_data = df[df["chrom"] == chrom]
            if len(chrom_data) < 4:
                continue
            Q1 = chrom_data["coverage"].quantile(0.25)
            Q3 = chrom_data["coverage"].quantile(0.75)
            IQR = Q3 - Q1
            lower_bound = Q1 - 1.5 * IQR
            upper_bound = Q3 + 1.5 * IQR
            outliers = chrom_data[
                (chrom_data["coverage"] < lower_bound)
                | (chrom_data["coverage"] > upper_bound)
            ]
            for _, out in outliers.iterrows():
                cov = out["coverage"]
                event_type = "Gain" if cov > upper_bound else "Loss"
                rows.append({
                    "chrom": chrom,
                    "gene": out["name"],
                    "coverage": cov,
                    "type": event_type,
                })
        return rows

    def _create_target_boxplot(self):
        """Create a boxplot showing coverage distribution for targets."""
        try:
            # Check if data is available
            if self.target_coverage_df.empty:
                plt.figure(figsize=(8, 4))
                plt.text(0.5, 0.5, "No target coverage data available", 
                        ha='center', va='center', transform=plt.gca().transAxes,
                        fontsize=14, color='gray')
                plt.title("Target Coverage Distribution")
                plt.axis('off')
            else:
                # Prepare data
                self.target_coverage_df["coverage"] = self.target_coverage_df["coverage"].round(
                    2
                )

                # Create figure and boxplot
                plt.figure(figsize=(8, 4))  # Reduced size

                # Get sorted unique chromosomes
                sorted_chroms = sorted(
                    self.target_coverage_df["chrom"].unique(), key=natsort.natsort_keygen()
                )

                # Create boxplot
                bp = plt.boxplot(
                    [
                        self.target_coverage_df[self.target_coverage_df["chrom"] == chrom][
                            "coverage"
                        ]
                        for chrom in sorted_chroms
                    ],
                    patch_artist=True,
                )

                # Style the boxplot
                plt.setp(bp["boxes"], facecolor="#2C3E50", alpha=0.6)
                plt.setp(bp["medians"], color="#E74C3C")
                plt.setp(bp["fliers"], marker="o", markerfacecolor="#C0392B")

                # Outlier labels removed - shown in table below for readability

                # Customize plot
                plt.xlabel("Chromosome")
                plt.ylabel("Coverage Depth")
                plt.title("Coverage Distribution by Chromosome")
                plt.xticks(range(1, len(sorted_chroms) + 1), sorted_chroms, rotation=45)
                plt.grid(True, alpha=0.3)

                # Adjust layout to prevent label cutoff
                plt.tight_layout()

            # Save plot to bytes buffer
            buf = io.BytesIO()
            plt.savefig(buf, format="png", dpi=300, bbox_inches="tight")
            plt.close()
            buf.seek(0)
            
            # Verify the buffer has data
            if buf.getvalue():
                return buf
            else:
                return self._create_empty_plot_buffer()
                
        except Exception as e:
            logger.warning(f"Error creating target boxplot: {e}")
            plt.close('all')  # Close any open figures
            return self._create_empty_plot_buffer()

    def add_content(self):
        """Add the coverage analysis content to the report."""
        # Add page break before detailed section
        self.elements.append(PageBreak())
        # Add Summary Section
        self.elements.append(
            Paragraph("Coverage Analysis", self.styles.styles["Heading2"])
        )
        self.elements.append(Spacer(1, 6))
        self.add_coverage_summary()

        self.add_detailed_coverage()

        # Add summary to summary section
        self.summary_elements.append(
            Paragraph("Coverage Analysis", self.styles.styles["Heading3"])
        )
        self.summary_elements.append(Spacer(1, 6))

        if hasattr(self, "coverage_data") and self.coverage_data:
            summary_text = []
            global_coverage = self.coverage_data.get("global_coverage", 0)
            target_coverage = self.coverage_data.get("target_coverage", 0)
            enrichment = target_coverage / global_coverage if global_coverage > 0 else 0

            summary_text.append(f"Target Coverage: {target_coverage:.2f}x")
            summary_text.append(f"Global Coverage: {global_coverage:.2f}x")
            summary_text.append(f"Enrichment: {enrichment:.2f}x")

            if hasattr(self, "distribution_data") and self.distribution_data:
                summary_text.append(
                    f"Coverage >=30x: {self.distribution_data.get('above_30x', 0):.1f}%"
                )

            self.summary_elements.append(
                Paragraph(" <br/> ".join(summary_text), self.styles.styles["Normal"])
            )
            # Export coverage summary and tables
            try:
                import pandas as pd

                self.export_frames["coverage_summary"] = pd.DataFrame(
                    [
                        {
                            "GlobalCoverageX": float(global_coverage),
                            "TargetCoverageX": float(target_coverage),
                            "Enrichment": float(enrichment),
                        }
                    ]
                )
                if hasattr(self, "cov_df_main"):
                    self.export_frames["coverage_chromosome"] = self.cov_df_main.copy()
                if hasattr(self, "bedcov_df_main"):
                    self.export_frames["coverage_targets_bed"] = (
                        self.bedcov_df_main.copy()
                    )
                if hasattr(self, "target_coverage_df"):
                    self.export_frames["coverage_targets"] = (
                        self.target_coverage_df.copy()
                    )
                if hasattr(self, "distribution_data") and self.distribution_data:
                    self.export_frames["coverage_distribution"] = pd.DataFrame(
                        [self.distribution_data]
                    )
            except Exception:
                pass
        else:
            self.summary_elements.append(
                Paragraph("No coverage data available", self.styles.styles["Normal"])
            )

    def add_coverage_summary(self):
        """Add the coverage summary section."""
        # Create styles for the summary text
        summary_style = ParagraphStyle(
            "SummaryStyle",
            parent=self.styles.styles["Normal"],
            fontSize=8,
            leading=10,
            spaceBefore=8,
            spaceAfter=8,
        )

        # Add overview text
        overview = Paragraph(
            "This section provides an analysis of sequencing coverage across the genome and targeted regions.",
            summary_style,
        )
        self.elements.append(overview)
        self.elements.append(Spacer(1, 4))

        # Create summary table with key metrics
        if hasattr(self, "coverage_data") and self.coverage_data is not None:
            global_coverage = self.coverage_data.get("global_coverage", 0)
            target_coverage = self.coverage_data.get("target_coverage", 0)
            enrichment = target_coverage / global_coverage if global_coverage > 0 else 0
            quality_level = self._get_coverage_quality(target_coverage)
            data = [
                ["Metric", "Value", "Status"],
                ["Global Coverage", f"{global_coverage:.2f}x", ""],
                ["Target Coverage", f"{target_coverage:.2f}x", quality_level],
                ["Enrichment Factor", f"{enrichment:.2f}x", ""],
            ]
        else:
            data = [
                ["Metric", "Value", "Status"],
                ["Global Coverage", "N/A", ""],
                ["Target Coverage", "N/A", ""],
                ["Enrichment Factor", "N/A", ""],
            ]

        # Create and style the table
        table = Table(data, colWidths=[2 * inch, 2 * inch, 2 * inch])
        table.setStyle(
            TableStyle(
                [
                    *self.MODERN_TABLE_STYLE._cmds,
                    ("ALIGN", (0, 0), (0, -1), "LEFT"),  # Metric column left-aligned
                    ("ALIGN", (1, 0), (1, -1), "RIGHT"),  # Value column right-aligned
                    ("ALIGN", (2, 0), (2, -1), "CENTER"),  # Status column centered
                ]
            )
        )
        self.elements.append(table)
        self.elements.append(Spacer(1, 6))

        # Add quality thresholds legend
        legend_style = ParagraphStyle(
            "LegendStyle",
            parent=self.styles.styles["Normal"],
            fontSize=8,
            textColor=self.styles.COLORS["muted"],
        )

        legend_text = """
        Coverage Quality Thresholds:
        • Excellent: >30x coverage
        • Good: >20x coverage
        • Moderate: >10x coverage
        • Insufficient: <10x coverage
        """
        legend = Paragraph(legend_text, legend_style)
        self.elements.append(legend)

    def add_detailed_coverage(self):
        """Add the detailed coverage analysis section."""
        if not hasattr(self, "cov_df_main"):
            return
        try:
            # Create plots and convert to Images with specified dimensions
            chrom_buf = self._create_chromosome_coverage_plot()
            target_buf = self._create_target_coverage_plot()
            box_buf = self._create_target_boxplot()
            
            # Validate buffers before creating Image objects
            if chrom_buf and chrom_buf.getvalue():
                chrom_plot = Image(chrom_buf, width=6 * inch, height=3 * inch)
                self.elements.append(chrom_plot)
            else:
                logger.warning("Skipping chromosome coverage plot - invalid buffer")
            
            if target_buf and target_buf.getvalue():
                target_plot = Image(target_buf, width=6 * inch, height=3 * inch)
                self.elements.append(target_plot)
            else:
                logger.warning("Skipping target coverage plot - invalid buffer")
            
            if box_buf and box_buf.getvalue():
                box_plot = Image(box_buf, width=6 * inch, height=3 * inch)
                self.elements.append(box_plot)
                # Add outliers table (replaces unreadable plot labels)
                outliers = self._get_target_coverage_outliers()
                if outliers:
                    try:
                        self.export_frames["coverage_outliers"] = pd.DataFrame(outliers)
                    except Exception:
                        pass
                    self.elements.append(Spacer(1, 4))
                    self.elements.append(
                        Paragraph(
                            "Coverage outliers (potential gains/losses)",
                            ParagraphStyle(
                                "CoverageOutliersTitle",
                                parent=self.styles.styles["Normal"],
                                fontSize=9,
                                fontName="Helvetica-Bold",
                                spaceAfter=4,
                            ),
                        )
                    )
                    table_data = [["Chromosome", "Gene", "Coverage", "Type"]]
                    for o in outliers:
                        table_data.append([
                            o["chrom"],
                            o["gene"],
                            f"{o['coverage']:.1f}x",
                            o["type"],
                        ])
                    outliers_table = self.create_table(
                        table_data,
                        repeat_rows=1,
                        auto_col_width=True,
                        compact=True,
                        font_size=9,
                    )
                    outliers_table.setStyle(
                        TableStyle([
                            *self.MODERN_TABLE_STYLE._cmds,
                            ("ALIGN", (2, 1), (2, -1), "RIGHT"),
                            ("ALIGN", (3, 1), (3, -1), "CENTER"),
                            ("TOPPADDING", (0, 0), (-1, -1), 4),
                            ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
                        ])
                    )
                    self.elements.append(outliers_table)
                    self.elements.append(Spacer(1, 6))
            else:
                logger.warning("Skipping box plot - invalid buffer")

            # Add coverage distribution statistics
            self._add_coverage_distribution()
        except Exception as e:
            logger.error(f"Error processing detailed coverage section: {str(e)}")
            # Continue execution - don't break the entire report

    def _add_coverage_distribution(self):
        """Add coverage distribution analysis."""
        if not hasattr(self, "distribution_data") or self.distribution_data is None:
            return

        self.elements.append(
            Paragraph("Coverage Distribution Analysis", self.styles.styles["Heading3"])
        )
        self.elements.append(Spacer(1, 4))

        # Create table data for statistics (with safe formatting for missing values)
        def _fmt(val, fmt_str: str, default: str = "N/A"):
            if val is None or (isinstance(val, str) and val == "N/A"):
                return default
            try:
                return fmt_str.format(float(val))
            except (TypeError, ValueError):
                return default

        median = self.distribution_data.get("median")
        mean = self.distribution_data.get("mean")
        min_cov = self.distribution_data.get("min")
        max_cov = self.distribution_data.get("max")
        above_30 = self.distribution_data.get("above_30x")
        above_20 = self.distribution_data.get("above_20x")
        above_10 = self.distribution_data.get("above_10x")

        stats_data = [
            ["Median Coverage", _fmt(median, "{:.2f}x")],
            ["Mean Coverage", _fmt(mean, "{:.2f}x")],
            [
                "Coverage Range",
                f"{_fmt(min_cov, '{:.2f}x')} - {_fmt(max_cov, '{:.2f}x')}",
            ],
            ["Regions >=30x", _fmt(above_30, "{:.1f}%")],
            ["Regions >=20x", _fmt(above_20, "{:.1f}%")],
            ["Regions >=10x", _fmt(above_10, "{:.1f}%")],
        ]

        # Create and style the table
        stats_table = Table(stats_data, colWidths=[2 * inch, 1.5 * inch])
        stats_table.setStyle(
            TableStyle(
                [
                    *self.MODERN_TABLE_STYLE._cmds,
                    ("ALIGN", (0, 0), (0, -1), "LEFT"),  # Metric names left-aligned
                    ("ALIGN", (1, 0), (1, -1), "RIGHT"),  # Values right-aligned
                ]
            )
        )

        self.elements.append(stats_table)
        self.elements.append(Spacer(1, 6))

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
