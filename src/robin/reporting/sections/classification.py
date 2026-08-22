"""
Classification Section for ROBIN Reports.

This module handles the methylation-based classification section of the report,
including results from Sturgeon, Random Forest, NanoDX, and PannanoDX classifiers.
"""

import logging
import os

import matplotlib
import pandas as pd
from reportlab.lib.colors import HexColor
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import Image, PageBreak, Paragraph, Spacer, Table

from ..sections.base import ReportSection

matplotlib.use("Agg")  # ensure non-interactive backend before importing pyplot
import io

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns

from robin.analysis.mnpflex_hierarchy import mnpflex_hierarchy_has_content
from robin.classification_config import get_confidence_status
from robin.reporting.mnpflex_hierarchy import (
    append_mnpflex_classifier_prediction,
    append_mnpflex_top_path_paragraph,
)

logger = logging.getLogger(__name__)

# Score CSV columns that are run statistics / metadata, not predicted classes.
# Kept as a set of lowercase names so case variants are also excluded.
CLASSIFIER_SCORE_META_COLUMNS = frozenset(
    {
        "timestamp",
        "number_probes",
        "covered_cpgs",
        "temperature",
        "diagnostic",
        "probes",
    }
)


def drop_classifier_score_meta_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Remove non-class columns from a classifier scores DataFrame."""
    if df is None or getattr(df, "empty", True):
        return df
    drop_cols = [
        col
        for col in df.columns
        if str(col).strip().lower() in CLASSIFIER_SCORE_META_COLUMNS
    ]
    if not drop_cols:
        return df
    return df.drop(columns=drop_cols)


class ClassificationSection(ReportSection):
    """Section containing the methylation classification results."""

    def _create_empty_plot_buffer(self):
        """Create a minimal valid PNG buffer for empty plots."""
        try:
            plt.figure(figsize=(8, 5), facecolor="white")
            plt.text(
                0.5,
                0.5,
                "No data available",
                ha="center",
                va="center",
                transform=plt.gca().transAxes,
                fontsize=14,
                color="gray",
            )
            plt.title("Classification Plot")
            plt.axis("off")

            buf = io.BytesIO()
            plt.savefig(
                buf, format="png", dpi=300, bbox_inches="tight", facecolor="white"
            )
            plt.close()
            buf.seek(0)

            # Validate buffer contains data
            if buf.getvalue():
                return buf
            else:
                raise ValueError("Empty buffer")
        except Exception:
            # If even this fails, return a minimal PNG
            import base64

            # Minimal 1x1 transparent PNG
            png_data = base64.b64decode(
                "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg=="
            )
            buf = io.BytesIO(png_data)
            buf.seek(0)
            return buf

    def _create_time_plot(self, df, classifier_name, figsize=(4.5, 2.75)):
        """Create a time series plot for classifier predictions.

        Args:
            df: DataFrame containing classification data
            classifier_name: Name of the classifier
            figsize: Figure size (width, height) in inches for compact layout

        Returns:
            BytesIO object containing the plot image
        """
        try:
            # Drop non-classification columns (covered_cpgs, diagnostic, etc.)
            df = drop_classifier_score_meta_columns(df)

            # Check if dataframe is empty or has no data
            if df.empty or len(df) == 0:
                logger.warning(
                    f"Empty dataframe for {classifier_name} time series plot"
                )
                return self._create_empty_plot_buffer()

            # Convert confidence values to percentages if not already
            if (
                classifier_name != "Random Forest"
            ):  # Random Forest is already in percentages
                df = df * 100

            # Get top 3 classifications (those that exceed 5% at any point) to reduce clutter
            threshold = 5
            above_threshold = df.columns[df.max() > threshold]
            if len(above_threshold) <= 3:
                top_classes = above_threshold
            else:
                final_vals = df.loc[df.index[-1], above_threshold].sort_values(
                    ascending=False
                )
                top_classes = final_vals.head(3).index

            # Check if we have any classes to plot
            if len(top_classes) == 0:
                logger.warning(
                    f"No classes above threshold for {classifier_name} time series plot"
                )
                return self._create_empty_plot_buffer()

            # Get current highest confidence for subtitle
            last_row = df.iloc[-1]
            top_prediction = last_row.sort_values(ascending=False).head(1)
            predicted_class = top_prediction.index[0]
            confidence_value = float(top_prediction.values[0])

            # Convert index to datetime
            df.index = pd.to_datetime(df.index, unit="ms")

            # Reshape to long format for seaborn
            plot_df = df[top_classes].copy().reset_index()
            # Index column may be named "timestamp" or "index" depending on source
            time_col = plot_df.columns[0]
            plot_df = plot_df.rename(columns={time_col: "Time"})
            plot_df = plot_df.melt(
                id_vars=["Time"],
                value_vars=top_classes,
                var_name="Class",
                value_name="Confidence",
            )
            # Truncate long labels for legend readability, ensuring no collisions
            # (e.g. "Embryonal - MB G3" and "Embryonal - MB G3G4 - G3" both become
            # "Embryonal - MB G3..." if truncated naively, merging distinct classes)
            unique_classes = plot_df["Class"].unique()
            truncated_map = {}
            for c in unique_classes:
                base = c if len(c) <= 20 else c[:17] + "..."
                # Ensure uniqueness: if collision, extend until distinct
                out = base
                while out in truncated_map.values() and out != c:
                    if len(c) <= len(out) + 3:
                        out = c  # Use full name to avoid collision
                        break
                    out = c[: min(len(out) + 4, len(c))] + (
                        "..." if len(c) > 20 else ""
                    )
                truncated_map[c] = out
            plot_df["Class"] = plot_df["Class"].map(truncated_map)

            # Palette must match actual number of hue levels (avoids seaborn warning)
            n_hue = plot_df["Class"].nunique()
            palette = sns.color_palette("husl", n_colors=n_hue)

            # Create plot with seaborn (use context to avoid affecting other report plots)
            with sns.axes_style("whitegrid"):
                fig, ax = plt.subplots(figsize=figsize, facecolor="white")
                ax.set_facecolor("white")

                sns.lineplot(
                    data=plot_df,
                    x="Time",
                    y="Confidence",
                    hue="Class",
                    palette=palette,
                    linewidth=2,
                    markers=True,
                    markersize=6,
                    markevery=max(1, len(df) // 15),
                    ax=ax,
                )

                # Add title and subtitle (compact for 2x2 grid)
                fig.suptitle(
                    f"{classifier_name}", y=0.98, fontsize=10, fontweight="bold"
                )
                ax.set_title(
                    f"{predicted_class} ({confidence_value:.1f}%)",
                    pad=4,
                    fontsize=8,
                )

                # Format x-axis to show HH:MM
                ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
                ax.set_xlabel("Time", fontsize=8)
                ax.set_ylabel("Confidence (%)", fontsize=8)
                ax.yaxis.set_major_formatter(
                    mticker.PercentFormatter(xmax=100, decimals=0)
                )

                # Set axis ranges
                ax.set_ylim(0, 100)

                # Legend - seaborn places it automatically; adjust for compact layout
                ax.legend(
                    bbox_to_anchor=(0, 1.08, 1, 0),
                    loc="lower left",
                    mode="expand",
                    ncol=2,
                    fontsize=8,
                    frameon=True,
                    framealpha=0.9,
                )

                # Refined spines (seaborn whitegrid already has grid)
                sns.despine(ax=ax, top=True, right=True)

                # Tick parameters
                ax.tick_params(axis="both", labelsize=7)
                plt.xticks(rotation=0)

                plt.tight_layout()

            # Save plot to bytes buffer with high DPI for crisp rendering
            buf = io.BytesIO()
            fig.savefig(
                buf, format="png", dpi=300, bbox_inches="tight", facecolor="white"
            )
            plt.close(fig)
            buf.seek(0)

            # Validate buffer contains valid image data
            if not buf.getvalue():
                logger.warning(f"Empty buffer for {classifier_name} time series plot")
                return self._create_empty_plot_buffer()

            # Try to validate it's a valid PNG by checking magic bytes
            buf_data = buf.getvalue()
            if len(buf_data) < 8 or buf_data[:8] != b"\x89PNG\r\n\x1a\n":
                logger.warning(
                    f"Invalid PNG data for {classifier_name} time series plot"
                )
                return self._create_empty_plot_buffer()

            # Reset buffer position after validation
            buf.seek(0)
            return buf

        except Exception as e:
            logger.error(f"Error creating {classifier_name} time series plot: {str(e)}")
            plt.close("all")  # Close any open figures
            return self._create_empty_plot_buffer()

    def add_content(self):
        """Add the classification content to the report."""
        # Add MNP-Flex hierarchical summary to summary section (if MNP-Flex enabled)
        try:
            from robin.reporting.sections.mnpflex import MNPFlexSection

            mnpflex = MNPFlexSection(self.report)
            results_dir = mnpflex._find_results_dir()
            if results_dir:
                summary = mnpflex._load_bundle_summary(results_dir)
            else:
                summary = None

            if summary:
                qc = summary.get("qc", {}) or {}
                mgmt = summary.get("mgmt", {}) or {}
                classifier_summary = summary.get("classifier_summary", {}) or {}
                hierarchy = classifier_summary.get("summary_hierarchical", []) or []
                has_hierarchical_summary = mnpflex_hierarchy_has_content(hierarchy)

                # Build legend text (QC, MGMT) for table caption
                legend_text = (
                    f"QC status: {qc.get('status', 'Unknown')} | "
                    f"Average coverage: {mnpflex._format_coverage_value(qc.get('avg_coverage'))} | "
                    f"Missing sites: {qc.get('missing_site_count', 'N/A')} | "
                    f"MGMT status: {mgmt.get('status', 'Unknown')} | "
                    f"MGMT average: {mnpflex._format_mgmt_value(mgmt.get('average'))}"
                )

                if not has_hierarchical_summary:
                    # No confident classification - show message + legend in same place
                    self.summary_elements.append(
                        Paragraph(
                            "MNP-Flex Hierarchical Summary",
                            self.styles.styles["Heading3"],
                        )
                    )
                    self.summary_elements.append(Spacer(1, 6))
                    self.summary_elements.append(
                        Paragraph(
                            "MNP-Flex has not returned a confident classification.",
                            self.styles.styles["Warning"],
                        )
                    )
                    self.summary_elements.append(Spacer(1, 6))
                    self.summary_elements.append(
                        Paragraph(
                            legend_text,
                            ParagraphStyle(
                                "TableLegend",
                                parent=self.styles.styles["Normal"],
                                fontSize=8,
                                leading=10,
                                textColor=HexColor("#4B5563"),
                            ),
                        )
                    )
                    self.summary_elements.append(Spacer(1, 8))
                else:
                    # Has hierarchical summary - nested tree with legend
                    self.summary_elements.append(
                        Paragraph(
                            "MNP-Flex Hierarchical Summary",
                            self.styles.styles["Heading3"],
                        )
                    )
                    self.summary_elements.append(Spacer(1, 6))
                    append_mnpflex_classifier_prediction(
                        self.summary_elements,
                        hierarchy,
                        styles=self.styles,
                        page_width=self.report.doc.width,
                        include_heading=False,
                    )
                    append_mnpflex_top_path_paragraph(
                        self.summary_elements,
                        hierarchy,
                        styles=self.styles,
                        format_score=mnpflex._format_score,
                    )
                    # Table legend
                    self.summary_elements.append(
                        Paragraph(
                            legend_text,
                            ParagraphStyle(
                                "TableLegend",
                                parent=self.styles.styles["Normal"],
                                fontSize=8,
                                leading=10,
                                textColor=HexColor("#4B5563"),
                            ),
                        )
                    )
                    self.summary_elements.append(Spacer(1, 8))
        except Exception as e:
            logger.error(
                f"Error adding MNP-Flex hierarchy to classification summary: {e}"
            )

        # Methylation Classification header - divides MNP-Flex from other classifiers
        self.summary_elements.append(
            Paragraph("Methylation Classification", self.styles.styles["Heading3"])
        )
        self.summary_elements.append(Spacer(1, 6))

        # Add section header for detailed view
        self.elements.append(
            Paragraph("Methylation Classification", self.styles.styles["Heading2"])
        )
        self.elements.append(Spacer(1, 4))

        # Dictionary of classifiers and their corresponding files
        classifiers = {
            "Sturgeon": "sturgeon_scores.csv",
            "NanoDX": "nanodx_scores.csv",
            "PanNanoDX": "pannanodx_scores.csv",
            "Random Forest": "random_forest_scores.csv",
            "MARLIN": "marlin_scores.csv",
            "Lamprey (research)": "lamprey_scores.csv",
            "Tucan": "tucan_scores.csv",
        }

        # Add summary table of all classifications
        summary_data = [["Classifier", "Predicted Class", "Confidence", "Status"]]
        plot_buffers = []  # Collect (name, buf) for 2x2 grid

        # Process each classifier
        for name, filename in classifiers.items():
            try:
                file_path = None
                for f in os.listdir(self.report.output):
                    if f.lower() == filename.lower():
                        file_path = os.path.join(self.report.output, f)
                        break

                if file_path and os.path.exists(file_path):
                    df = pd.read_csv(file_path)

                    # Create time series plot (compact for 2x2 grid)
                    try:
                        df_plot = df.copy()
                        df_plot.set_index("timestamp", inplace=True)
                        plot_buf = self._create_time_plot(df_plot, name)
                        plot_buffers.append((name, plot_buf))
                    except Exception as e:
                        logger.error(
                            f"Error creating {name} time series plot: {str(e)}"
                        )

                    # Continue with existing classification processing
                    df = drop_classifier_score_meta_columns(df)

                    # Get the last row and find top prediction
                    last_row = df.iloc[-1]
                    top_prediction = last_row.sort_values(ascending=False).head(1)
                    predicted_class = top_prediction.index[0]
                    raw_confidence = float(top_prediction.values[0])

                    # Normalize confidence value
                    confidence_value = (
                        raw_confidence / 100.0
                        if name == "Random Forest"
                        else raw_confidence
                    )

                    # Determine confidence level based on classifier using centralized config
                    classifier_key = {
                        "Sturgeon": "sturgeon",
                        "NanoDX": "nanodx",
                        "PanNanoDX": "pannanodx",
                        "Random Forest": "random_forest",
                        "MARLIN": "marlin",
                        "Lamprey (research)": "lamprey",
                        "Tucan": "tucan",
                    }.get(name, name.lower().replace(" ", "_"))
                    confidence_status, status_color = get_confidence_status(
                        classifier_key, confidence_value
                    )

                    # Add to summary table with HTML-like color formatting
                    summary_data.append(
                        [
                            name,
                            predicted_class,
                            f"{confidence_value:.1%}",
                            f'<font color="{status_color}">{confidence_status}</font>',
                        ]
                    )

            except Exception as e:
                logger.error(f"Error processing {name} classification: {str(e)}")
                continue

        # Add classification time series plots in 2x2 grid (4 per page)
        if plot_buffers:
            plot_width = 3.25 * inch
            plot_height = 2.25 * inch
            plots_per_row = 2
            # Pad to complete last row
            while len(plot_buffers) % plots_per_row != 0:
                plot_buffers.append((None, None))
            for row_start in range(0, len(plot_buffers), plots_per_row):
                row_cells = []
                for name, buf in plot_buffers[row_start : row_start + plots_per_row]:
                    if buf is not None:
                        buf.seek(0)
                        row_cells.append(
                            Image(buf, width=plot_width, height=plot_height)
                        )
                    else:
                        row_cells.append(Spacer(plot_width, plot_height))
                plot_table = Table(
                    [row_cells],
                    colWidths=[plot_width] * plots_per_row,
                )
                plot_table.setStyle(
                    [
                        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                    ]
                )
                self.elements.append(plot_table)
                self.elements.append(Spacer(1, 4))

        # Add explanation text using centralized thresholds
        from robin.classification_config import CLASSIFIER_CONFIDENCE_THRESHOLDS

        explanation_lines = [
            "Note: Classification confidence levels are defined as follows:"
        ]
        for classifier, thresholds in CLASSIFIER_CONFIDENCE_THRESHOLDS.items():
            explanation_lines.append(
                f"- {classifier.title()}: High (>={thresholds['high']:.0f}%), "
                f"Medium (>={thresholds['medium']:.0f}%), Low (<{thresholds['medium']:.0f}%)"
            )
        explanation_lines.append(
            "Multiple classifiers may provide different results based on their training data and methodology."
        )

        Explanation_text = Paragraph(
            "\n".join(explanation_lines),
            ParagraphStyle(
                "Explanation",
                parent=self.styles.styles["Normal"],
                fontSize=8,
                leading=10,
                textColor=HexColor("#4B5563"),
            ),
        )

        if len(summary_data) > 1:  # If we have any results
            # Create table with automatic column width calculation
            summary_table = self.create_table(summary_data, repeat_rows=1)
            self.summary_elements.append(summary_table)
            self.elements.append(Spacer(1, 6))

        # Export summary as DataFrame
        try:
            # Strip HTML tags from Status for CSV, keep plain text
            def _strip_html(s):
                try:
                    return (
                        str(s)
                        .replace('<font color="#059669">', "")
                        .replace('<font color="#D97706">', "")
                        .replace('<font color="#DC2626">', "")
                        .replace("</font>", "")
                    )
                except Exception:
                    return s

            rows = []
            for r in summary_data[1:]:
                rows.append(
                    {
                        "Classifier": r[0],
                        "PredictedClass": r[1],
                        "ConfidencePercent": r[2],
                        "Status": _strip_html(r[3]),
                    }
                )
            from pandas import DataFrame as _DF

            self.export_frames["classification_summary"] = _DF(rows)
        except Exception:
            pass

        # Add detailed results for each classifier
        self.elements.append(PageBreak())
        self.elements.append(
            Paragraph("Detailed Classification Results", self.styles.styles["Heading2"])
        )
        self.elements.append(Spacer(1, 4))

        # Collect classifier blocks for 2-column layout
        block_width = self.report.doc.width / 2
        table_col_widths = [block_width * 0.65, block_width * 0.35]  # Class, Score
        classifier_blocks = []

        for name, filename in classifiers.items():
            try:
                file_path = None
                for f in os.listdir(self.report.output):
                    if f.lower() == filename.lower():
                        file_path = os.path.join(self.report.output, f)
                        break

                if file_path and os.path.exists(file_path):
                    df = pd.read_csv(file_path)
                    df = drop_classifier_score_meta_columns(df)

                    last_row = df.iloc[-1]
                    top_predictions = last_row.sort_values(ascending=False).head(10)

                    # Build detailed table for top 10 predictions
                    detailed_data = [["Predicted Class", "Confidence Score"]]
                    for class_name, score in top_predictions.items():
                        confidence = score / 100.0 if name == "Random Forest" else score
                        detailed_data.append([class_name, f"{confidence:.1%}"])

                    detailed_table = self.create_table(
                        detailed_data,
                        col_widths=table_col_widths,
                        repeat_rows=1,
                        auto_col_width=False,
                        compact=True,
                    )
                    detailed_table.setStyle(
                        [
                            ("TOPPADDING", (0, 0), (-1, -1), 2),
                            ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
                            ("LEFTPADDING", (0, 0), (-1, -1), 4),
                            ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                        ]
                    )

                    # Block: header + table
                    header = Paragraph(
                        f"{name} Classification", self.styles.styles["Heading3"]
                    )
                    block = Table(
                        [[header], [detailed_table]],
                        colWidths=[block_width],
                    )
                    block.setStyle(
                        [
                            ("TOPPADDING", (0, 0), (-1, -1), 0),
                            ("BOTTOMPADDING", (0, 0), (-1, 0), 2),
                            ("BOTTOMPADDING", (0, 1), (-1, -1), 0),
                        ]
                    )
                    classifier_blocks.append((name, block, top_predictions))

                    # Export top predictions
                    try:
                        rows = []
                        for class_name, score in top_predictions.items():
                            confidence = (
                                score / 100.0 if name == "Random Forest" else score
                            )
                            rows.append(
                                {
                                    "Classifier": name,
                                    "PredictedClass": class_name,
                                    "ConfidencePercent": f"{confidence:.1%}",
                                }
                            )
                        key = f"classification_{name.lower().replace(' ', '_')}_top10"
                        from pandas import DataFrame as _DF

                        self.export_frames[key] = _DF(rows)
                    except Exception:
                        pass

            except Exception as e:
                logger.error(
                    f"Error processing detailed {name} classification: {str(e)}"
                )
                continue

        # Arrange blocks in 2-column grid
        if classifier_blocks:
            plots_per_row = 2
            for row_start in range(0, len(classifier_blocks), plots_per_row):
                row_cells = []
                for i in range(plots_per_row):
                    idx = row_start + i
                    if idx < len(classifier_blocks):
                        _, block, _ = classifier_blocks[idx]
                        row_cells.append(block)
                    else:
                        row_cells.append(Spacer(block_width, 1))

                grid_table = Table(
                    [row_cells],
                    colWidths=[block_width] * plots_per_row,
                )
                grid_table.setStyle(
                    [
                        ("VALIGN", (0, 0), (-1, -1), "TOP"),
                        ("ALIGN", (0, 0), (-1, -1), "LEFT"),
                        ("LEFTPADDING", (1, 0), (1, -1), 8),  # Gap between columns
                    ]
                )
                self.elements.append(grid_table)
                self.elements.append(Spacer(1, 4))

        # Add explanation text to detailed section
        if len(summary_data) > 1:
            self.elements.append(Explanation_text)
        else:
            self.elements.append(
                Paragraph(
                    "No classification results available.", self.styles.styles["Normal"]
                )
            )

        # Add explanation text to summary
        self.summary_elements.append(Explanation_text)
