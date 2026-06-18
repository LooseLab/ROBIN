"""
plotting.py

This module contains functions for creating plots used in the PDF report.
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import io
import matplotlib.font_manager as fm
import os
from robin.gui import fonts

import natsort

from matplotlib import gridspec

import logging
from typing import List, Optional, Tuple, Dict, Any

logger = logging.getLogger(__name__)

# Define consistent color scheme and style
MODERN_COLORS = {
    "primary": "#2C3E50",  # Dark blue-grey (matching report text)
    "secondary": "#E2E8F0",  # Light grey (matching table grid)
    "background": "#F8FAFC",  # Light background (matching table alternate rows)
    "accent": "#3498DB",  # Blue accent
    "grid": "#E2E8F0",  # Grid color
}

CNV_COLORS = {
    "points": "#4A6FA5",
    "trend": "#2C5282",
    "reference": "#94A3B8",
    "gain_fill": "#C6E6C3",
    "loss_fill": "#F5C6C6",
    "gain_edge": "#2E7D32",
    "loss_edge": "#C62828",
}


def set_modern_style():
    """Set consistent modern style for all plots"""
    # Register FiraSans font
    font_path = os.path.join(
        os.path.dirname(os.path.abspath(fonts.__file__)),
        "fira-sans-v16-latin-regular.ttf",
    )

    # Add font to matplotlib's font manager
    fm.fontManager.addfont(font_path)
    prop = fm.FontProperties(fname=font_path)
    plt.rcParams["font.family"] = prop.get_name()

    plt.style.use("seaborn-v0_8-whitegrid")
    sns.set_theme(style="whitegrid")

    plt.rcParams.update(
        {
            # Figure settings
            "figure.facecolor": "white",
            "figure.dpi": 300,
            # Axes settings
            "axes.facecolor": "white",
            "axes.edgecolor": MODERN_COLORS["primary"],
            "axes.labelcolor": MODERN_COLORS["primary"],
            "axes.titlecolor": MODERN_COLORS["primary"],
            "axes.grid": True,
            "axes.labelsize": 10,
            "axes.titlesize": 12,
            # Grid settings
            "grid.color": MODERN_COLORS["grid"],
            "grid.linestyle": "--",
            "grid.linewidth": 0.5,
            "grid.alpha": 0.5,
            # Tick settings
            "xtick.color": MODERN_COLORS["primary"],
            "ytick.color": MODERN_COLORS["primary"],
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            # Legend settings
            "legend.frameon": True,
            "legend.facecolor": "white",
            "legend.edgecolor": MODERN_COLORS["grid"],
            "legend.fontsize": 8,
            # Line settings
            "lines.linewidth": 1.5,
            "lines.markersize": 6,
        }
    )


def _apply_cnv_chromosome_axes(
    ax,
    x_max_mb: float,
    y_max: float,
    *,
    xlabel: str,
    ylabel: str,
) -> None:
    """Pin axes at the origin with no leading x padding."""
    ax.set_xlim(0, x_max_mb)
    ax.set_ylim(0, y_max)
    ax.margins(x=0, y=0)
    ax.autoscale(enable=False)
    ax.set_xlabel(xlabel, fontsize=9, color=MODERN_COLORS["primary"])
    ax.set_ylabel(ylabel, fontsize=9, color=MODERN_COLORS["primary"])
    ax.spines["left"].set_position(("data", 0))
    ax.spines["bottom"].set_position(("data", 0))
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    ax.grid(True, axis="y", color=MODERN_COLORS["grid"], linestyle="--", linewidth=0.5, alpha=0.7)
    ax.grid(False, axis="x")
    ax.tick_params(colors=MODERN_COLORS["primary"], labelsize=8)


def _apply_cnv_axes_style(ax, *, xlabel: str, ylabel: str, title: Optional[str] = None) -> None:
    """Apply consistent seaborn-inspired styling to a CNV axes."""
    ax.set_facecolor("white")
    ax.set_xlabel(xlabel, fontsize=9, color=MODERN_COLORS["primary"])
    ax.set_ylabel(ylabel, fontsize=9, color=MODERN_COLORS["primary"])
    if title:
        ax.set_title(title, fontsize=10, color=MODERN_COLORS["primary"], pad=8)
    sns.despine(ax=ax, top=True, right=True)
    ax.grid(True, axis="y", color=MODERN_COLORS["grid"], linestyle="--", linewidth=0.5, alpha=0.7)
    ax.grid(False, axis="x")
    ax.tick_params(colors=MODERN_COLORS["primary"], labelsize=8)


def _chromosome_cnv_dataframe(positions_mb, values) -> pd.DataFrame:
    """Build a tidy dataframe for seaborn CNV plotting."""
    return pd.DataFrame(
        {
            "position_mb": positions_mb,
            "ploidy": pd.Series(values, dtype=float),
        }
    )


def _plot_cnv_track(
    ax,
    cnv_df: pd.DataFrame,
    *,
    point_size: float = 4,
    point_alpha: float = 0.3,
    show_trend: bool = False,
) -> None:
    """Render CNV bin values on an axis."""
    if show_trend and "rolling_ploidy" in cnv_df.columns:
        sns.lineplot(
            data=cnv_df,
            x="position_mb",
            y="rolling_ploidy",
            ax=ax,
            color=CNV_COLORS["trend"],
            linewidth=1.2,
            alpha=0.45,
            errorbar=None,
            zorder=2,
        )

    ax.scatter(
        cnv_df["position_mb"],
        cnv_df["ploidy"],
        s=point_size,
        c=CNV_COLORS["points"],
        alpha=point_alpha,
        linewidths=0,
        edgecolors="none",
        rasterized=True,
        zorder=3,
    )


def _add_panel_gene_table(
    ax_table,
    panel_points: List[Dict[str, Any]],
    off_scale_mode: bool,
    y_max: float,
    mean_cnv: float,
) -> None:
    """Render panel target names in a readable table below the CNV track."""
    ax_table.axis("off")
    if not panel_points:
        return

    ncol = 4
    off_scale_threshold = y_max * 0.97
    sorted_points = sorted(panel_points, key=lambda item: item["mid_mb"])
    rows: list[list[str]] = []
    highlight_cells: set[tuple[int, int]] = set()

    for row_idx in range(int(np.ceil(len(sorted_points) / ncol))):
        row: list[str] = []
        for col_idx in range(ncol):
            point_idx = row_idx * ncol + col_idx
            if point_idx >= len(sorted_points):
                row.append("")
                continue
            point = sorted_points[point_idx]
            is_off_scale = off_scale_mode and point["cnv_val"] > off_scale_threshold
            is_amplified = point["cnv_val"] > max(mean_cnv + 1.0, mean_cnv * 1.35)
            ploidy_note = f" ({point['cnv_val']:.1f})"
            if is_off_scale:
                ploidy_note += " ↑"
            row.append(f"{point['label']} @ {point['mid_mb']:.0f} Mb{ploidy_note}")
            if is_off_scale or is_amplified:
                highlight_cells.add((row_idx, col_idx))
        rows.append(row)

    table = ax_table.table(
        cellText=rows,
        loc="upper center",
        cellLoc="left",
        colWidths=[0.25] * ncol,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(6.5)
    table.scale(1.0, 1.25)
    for (row_idx, col_idx), cell in table.get_celld().items():
        cell.set_edgecolor(MODERN_COLORS["grid"])
        cell.set_linewidth(0.4)
        if (row_idx, col_idx) in highlight_cells:
            cell.set_text_props(color=AMPLIFIED_COLOR, weight="bold")
        else:
            cell.set_text_props(color=MODERN_COLORS["primary"])


def _add_cnv_reference_lines(ax, mean_cnv: float, std_cnv: float, y_min: float, y_max: float) -> None:
    """Add reference ploidy guides to a CNV track."""
    ax.axhline(y=mean_cnv, color=CNV_COLORS["reference"], linestyle="--", linewidth=1.0, alpha=0.8, zorder=1)
    for offset in (std_cnv, 2 * std_cnv):
        if y_min < mean_cnv + offset <= y_max:
            ax.axhline(
                y=mean_cnv + offset,
                color=CNV_COLORS["reference"],
                linestyle=":",
                linewidth=0.8,
                alpha=0.45,
                zorder=1,
            )
        if y_min <= mean_cnv - offset < y_max:
            ax.axhline(
                y=mean_cnv - offset,
                color=CNV_COLORS["reference"],
                linestyle=":",
                linewidth=0.8,
                alpha=0.45,
                zorder=1,
            )


def target_distribution_plot(df):
    """
    Creates a target distribution plot.

    Args:
        df (pd.DataFrame): DataFrame containing the target distribution data.

    Returns:
        io.BytesIO: Buffer containing the plot image.
    """
    set_modern_style()

    df["chrom"] = pd.Categorical(
        df["chrom"], categories=natsort.natsorted(df["chrom"].unique()), ordered=True
    )
    df = df.sort_values("chrom")

    # Generate the plot
    plt.figure(figsize=(16, 8))
    boxplot = sns.boxplot(
        x="chrom",
        y="coverage",
        data=df,
        color=MODERN_COLORS["accent"],
        flierprops={"marker": "o", "markerfacecolor": MODERN_COLORS["primary"]},
    )

    plt.title(
        "Distribution of Target Coverage on Each Chromosome",
        fontsize=12,
        color=MODERN_COLORS["primary"],
        pad=20,
    )
    plt.xlabel("Chromosome", color=MODERN_COLORS["primary"])
    plt.ylabel("Coverage", color=MODERN_COLORS["primary"])
    plt.xticks(rotation=45)

    # Identify and annotate outliers
    def annotate_outliers(df, boxplot):
        # Calculate quartiles and IQR
        for chrom in df["chrom"].unique():
            chrom_data = df[df["chrom"] == chrom]
            Q1 = chrom_data["coverage"].quantile(0.25)
            Q3 = chrom_data["coverage"].quantile(0.75)
            IQR = Q3 - Q1
            lower_bound = Q1 - 1.5 * IQR
            upper_bound = Q3 + 1.5 * IQR

            # Find outliers
            outliers = chrom_data[
                (chrom_data["coverage"] < lower_bound)
                | (chrom_data["coverage"] > upper_bound)
            ]

            for idx in outliers.index:
                outlier = outliers.loc[idx]
                boxplot.annotate(
                    outlier["name"],
                    xy=(df.loc[idx, "chrom"], df.loc[idx, "coverage"]),
                    xytext=(10, 10),  # Offset the text more from the point
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=12,
                    color="red",
                )  # Increase font size

    annotate_outliers(df, boxplot)

    plt.tight_layout()
    # Save the plot as a JPG file
    buf = io.BytesIO()
    plt.savefig(buf, format="jpg", dpi=300)
    buf.seek(0)
    return buf


def _create_empty_cnv_buffer():
    """Create a minimal valid JPEG buffer for empty CNV plots."""
    try:
        plt.figure(figsize=(16, 4))
        plt.text(0.5, 0.5, "No CNV data available", 
                ha='center', va='center', transform=plt.gca().transAxes,
                fontsize=14, color='gray')
        plt.title("Copy Number Changes")
        plt.axis('off')
        
        buf = io.BytesIO()
        fig = plt.gcf()
        plt.savefig(buf, format="jpg", dpi=300, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        
        # Validate buffer contains data
        if buf.getvalue():
            return buf
        else:
            raise ValueError("Empty buffer")
    except Exception:
        # If even this fails, create a minimal valid JPEG programmatically
        try:
            # Create a minimal 1x1 white JPEG
            from PIL import Image as PILImage
            img = PILImage.new('RGB', (1, 1), color='white')
            buf = io.BytesIO()
            img.save(buf, format='JPEG')
            buf.seek(0)
            return buf
        except Exception:
            # Last resort: return a minimal valid JPEG binary directly
            # This is a valid 1x1 white JPEG
            jpeg_bytes = bytes([
                0xFF, 0xD8, 0xFF, 0xE0, 0x00, 0x10, 0x4A, 0x46, 0x49, 0x46, 0x00, 0x01,
                0x01, 0x01, 0x00, 0x48, 0x00, 0x48, 0x00, 0x00, 0xFF, 0xDB, 0x00, 0x43,
                0x00, 0x08, 0x06, 0x06, 0x07, 0x06, 0x05, 0x08, 0x07, 0x07, 0x07, 0x09,
                0x09, 0x08, 0x0A, 0x0C, 0x14, 0x0D, 0x0C, 0x0B, 0x0B, 0x0C, 0x19, 0x12,
                0x13, 0x0F, 0x14, 0x1D, 0x1A, 0x1F, 0x1E, 0x1D, 0x1A, 0x1C, 0x1C, 0x20,
                0x24, 0x2E, 0x27, 0x20, 0x22, 0x2C, 0x23, 0x1C, 0x1C, 0x28, 0x37, 0x29,
                0x2C, 0x30, 0x31, 0x34, 0x34, 0x34, 0x1F, 0x27, 0x39, 0x3D, 0x38, 0x32,
                0x3C, 0x2E, 0x33, 0x34, 0x32, 0xFF, 0xC0, 0x00, 0x0B, 0x08, 0x00, 0x01,
                0x00, 0x01, 0x01, 0x01, 0x11, 0x00, 0xFF, 0xC4, 0x00, 0x14, 0x00, 0x01,
                0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                0x00, 0x00, 0x00, 0x08, 0xFF, 0xC4, 0x00, 0x14, 0x10, 0x01, 0x00, 0x00,
                0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
                0x00, 0x00, 0xFF, 0xDA, 0x00, 0x08, 0x01, 0x01, 0x00, 0x00, 0x3F, 0x00,
                0xD2, 0xCF, 0x20, 0xFF, 0xD9
            ])
            buf = io.BytesIO(jpeg_bytes)
            buf.seek(0)
            return buf


def create_CNV_plot(result, cnv_dict):
    """
    Creates a CNV plot.

    Args:
        result (Result): CNV result object.
        cnv_dict (dict): Dictionary containing CNV data.

    Returns:
        io.BytesIO: Buffer containing the plot image.
    """
    try:
        set_modern_style()

        # Check if result has CNV data
        if not hasattr(result, 'cnv') or not result.cnv:
            logger.warning("No CNV data available for plotting")
            return _create_empty_cnv_buffer()

        # Prepare data for plotting
        plot_rows = []
        offset = 0
        contig_centers = {}
        contig_boundaries = []
        reportable = ["chr" + str(i) for i in range(0, 23)] + ["chrX", "chrY"]
        ordered_contigs = [
            contig for contig in natsort.natsorted(result.cnv.keys()) if contig in reportable
        ]

        for contig in ordered_contigs:
            values = result.cnv[contig]
            start_offset = offset
            for i, value in enumerate(values):
                plot_rows.append(
                    {
                        "contig": contig,
                        "bin_index": i + offset,
                        "position_bp": (i + offset) * cnv_dict["bin_width"],
                        "ploidy": float(value),
                    }
                )
            end_offset = offset + len(values) - 1
            contig_centers[contig] = ((start_offset + end_offset) / 2) * cnv_dict["bin_width"]
            offset += len(values)
            contig_boundaries.append(offset * cnv_dict["bin_width"])

        if not plot_rows:
            logger.warning("No plot data available for CNV plot")
            return _create_empty_cnv_buffer()

        df = pd.DataFrame(plot_rows)
        mean_value = float(df["ploidy"].mean())
        std_value = float(df["ploidy"].std())
        y_min = max(0.0, float(df["ploidy"].min()) - 0.25)
        y_max = max(mean_value + (4 * std_value), mean_value * 1.35, 2.5)

        palette = sns.color_palette("muted", n_colors=max(len(ordered_contigs), 3))
        contig_palette = dict(zip(ordered_contigs, palette))

        width = 16
        fig, ax = plt.subplots(figsize=(width, width / 4))
        sns.scatterplot(
            data=df,
            x="position_bp",
            y="ploidy",
            hue="contig",
            palette=contig_palette,
            ax=ax,
            legend=False,
            s=4,
            alpha=0.28,
            linewidth=0,
            edgecolor=None,
            rasterized=True,
        )

        for boundary in contig_boundaries[:-1]:
            ax.axvline(
                boundary,
                color=MODERN_COLORS["grid"],
                linewidth=0.6,
                linestyle="-",
                alpha=0.8,
                zorder=0,
            )

        _add_cnv_reference_lines(ax, mean_value, std_value, y_min, y_max)

        label_y = y_min + (y_max - y_min) * 0.03
        for contig, center_bp in contig_centers.items():
            ax.text(
                center_bp,
                label_y,
                contig.replace("chr", ""),
                fontsize=7,
                ha="center",
                va="bottom",
                rotation=0,
                color=MODERN_COLORS["primary"],
            )

        ax.set_ylim(y_min, y_max)
        _apply_cnv_axes_style(
            ax,
            xlabel="Genomic position (bp)",
            ylabel="Estimated ploidy",
            title="Copy number variation across chromosomes",
        )

        buf = io.BytesIO()
        fig.savefig(buf, format="jpg", dpi=300, bbox_inches="tight", pad_inches=0.08)
        plt.close(fig)
        buf.seek(0)

        if not buf.getvalue():
            logger.warning("Empty buffer for CNV plot")
            return _create_empty_cnv_buffer()

        buf_data = buf.getvalue()
        if len(buf_data) < 2 or buf_data[:2] != b'\xff\xd8':
            logger.warning("Invalid JPEG data for CNV plot")
            return _create_empty_cnv_buffer()

        buf.seek(0)
        return buf
    except Exception as e:
        logger.error(f"Error creating CNV plot: {str(e)}")
        plt.close('all')
        return _create_empty_cnv_buffer()


REPORTABLE_CHROMOSOMES = ["chr" + str(i) for i in range(0, 23)] + ["chrX", "chrY"]
CNV_CHROMOSOME_FIG_WIDTH = 7.5
PANEL_Y_EXPANSION_FACTOR = 6.0
AMPLIFIED_COLOR = "#C62828"
PANEL_GENE_COLOR = "#6A1B9A"


def _panel_target_label(gene_row) -> str:
    """Return the display label for a target-panel region."""
    for key in ("gene", "name", "target"):
        if key in gene_row.index and pd.notna(gene_row[key]):
            raw = str(gene_row[key]).strip()
            if raw and raw.lower() != "nan":
                return raw.split(",")[0].strip()
    return ""


def _collect_panel_gene_points(
    panel_genes_df: Optional[pd.DataFrame],
    contig: str,
    values,
    bin_width: int,
) -> List[Dict[str, Any]]:
    """Collect panel target positions and peak CNV across each target region."""
    if panel_genes_df is None or panel_genes_df.empty:
        return []

    values_array = np.array(values)
    genes = panel_genes_df[panel_genes_df["chrom"] == contig].sort_values("start_pos")
    points: List[Dict[str, Any]] = []

    for _, gene_row in genes.iterrows():
        label_text = _panel_target_label(gene_row)
        if not label_text:
            continue

        start_bp = float(gene_row["start_pos"])
        end_bp = float(gene_row["end_pos"])
        mid_bp = (start_bp + end_bp) / 2.0
        mid_mb = mid_bp / 1_000_000
        start_bin = max(0, int(start_bp // bin_width))
        end_bin = min(len(values_array) - 1, int(end_bp // bin_width))
        if end_bin < start_bin:
            cnv_val = 0.0
        else:
            region_vals = values_array[start_bin : end_bin + 1]
            cnv_val = float(np.max(region_vals)) if len(region_vals) else 0.0

        points.append(
            {
                "label": label_text,
                "mid_mb": mid_mb,
                "cnv_val": cnv_val,
            }
        )

    return points


def _compute_cnv_y_limits(
    values_array: np.ndarray,
    panel_points: List[Dict[str, Any]],
) -> tuple[float, float, bool]:
    """
    Choose y-axis limits that preserve bulk CNV detail while surfacing
  amplified panel targets when feasible.
    """
    mean_cnv = float(np.mean(values_array))
    std_cnv = float(np.std(values_array))
    baseline_max = max(mean_cnv + (2 * std_cnv), mean_cnv * 1.4, 2.5)

    if not panel_points:
        return 0.0, baseline_max, False

    panel_max = max(point["cnv_val"] for point in panel_points)
    if panel_max <= baseline_max * 1.05:
        return 0.0, baseline_max, False

    if panel_max <= baseline_max * PANEL_Y_EXPANSION_FACTOR:
        return 0.0, panel_max * 1.08, False

    return 0.0, baseline_max, True


def _add_panel_gene_lollipops(
    ax,
    panel_points: List[Dict[str, Any]],
    y_min: float,
    y_max: float,
    off_scale_mode: bool,
    mean_cnv: float,
) -> tuple[int, int]:
    """Draw panel target lollipop markers on the CNV track (names go in the table)."""
    if not panel_points:
        return 0, 0

    off_scale_count = 0
    off_scale_threshold = y_max * 0.97

    for point in panel_points:
        mid_mb = point["mid_mb"]
        cnv_val = point["cnv_val"]
        is_off_scale = off_scale_mode and cnv_val > off_scale_threshold
        is_amplified = cnv_val > max(mean_cnv + 1.0, mean_cnv * 1.35)

        if is_off_scale:
            stem_top = y_max * 0.96
            marker_color = AMPLIFIED_COLOR
            off_scale_count += 1
        else:
            stem_top = min(cnv_val, y_max * 0.96)
            marker_color = AMPLIFIED_COLOR if is_amplified else PANEL_GENE_COLOR

        ax.plot(
            [mid_mb, mid_mb],
            [y_min, stem_top],
            color=marker_color,
            linewidth=0.9 if is_off_scale else 0.7,
            alpha=0.75,
            zorder=4,
            linestyle="--" if is_off_scale else "-",
        )

        if is_off_scale:
            ax.scatter(
                [mid_mb],
                [y_max * 0.99],
                s=50,
                marker="^",
                color=AMPLIFIED_COLOR,
                zorder=7,
                edgecolors="none",
                clip_on=False,
            )
        else:
            ax.scatter(
                [mid_mb],
                [cnv_val],
                s=14 if is_amplified else 10,
                color=marker_color,
                zorder=5,
                edgecolors="none",
                alpha=0.95,
            )

    return len(panel_points), off_scale_count


def create_CNV_plot_per_chromosome(
    result,
    cnv_dict,
    significant_regions=None,
    chromosomes: Optional[List[str]] = None,
    panel_genes_df: Optional[pd.DataFrame] = None,
):
    """Creates CNV plots per chromosome.

    Args:
        result (Result): CNV result object.
        cnv_dict (dict): Dictionary containing CNV data.
        significant_regions (dict): Dictionary mapping chromosomes to lists of significant regions.
            Each region should be a dict with keys: 'start_pos', 'end_pos', 'type' ('GAIN' or 'LOSS')
        chromosomes (list, optional): Ordered chromosome names to plot. Defaults to reportable
            chromosomes present in the result.
        panel_genes_df (pd.DataFrame, optional): Target panel genes with columns
            chrom, start_pos, end_pos, gene.

    Returns:
        List[Tuple[str, io.BytesIO]]: List of tuples containing chromosome names and plot buffers.
    """
    plots = []
    try:
        set_modern_style()

        # Check if result has CNV data
        if not hasattr(result, 'cnv') or not result.cnv:
            logger.warning("No CNV data available for per-chromosome plotting")
            return plots

        if chromosomes is None:
            chromosomes = [
                contig
                for contig in natsort.natsorted(result.cnv.keys())
                if contig in REPORTABLE_CHROMOSOMES
            ]

        for contig in chromosomes:
            if contig not in result.cnv:
                continue
            values = result.cnv[contig]
            if contig not in REPORTABLE_CHROMOSOMES:
                continue

            # Calculate mean and standard deviation for this chromosome
            values_array = np.array(values)
            mean_cnv = np.mean(values_array)
            std_cnv = np.std(values_array)

            # Calculate positions in megabases
            positions = (
                np.arange(len(values)) * cnv_dict["bin_width"] / 1_000_000
            )  # Convert to Mb
            x_max_mb = float(positions[-1]) if len(positions) else 1.0

            panel_points = _collect_panel_gene_points(
                panel_genes_df,
                contig,
                values,
                cnv_dict["bin_width"],
            )
            y_min, y_max, off_scale_mode = _compute_cnv_y_limits(values_array, panel_points)
            y_min = 0.0
            gene_count_on_chrom = len(panel_points)

            plot_height = CNV_CHROMOSOME_FIG_WIDTH / 3.6
            table_height = 0.0
            if gene_count_on_chrom:
                table_rows = int(np.ceil(gene_count_on_chrom / 4))
                table_height = 0.35 + 0.16 * table_rows

            fig = plt.figure(figsize=(CNV_CHROMOSOME_FIG_WIDTH, plot_height + table_height))
            if gene_count_on_chrom:
                gs = fig.add_gridspec(
                    2,
                    1,
                    height_ratios=[plot_height, table_height],
                    hspace=0.22,
                )
                ax = fig.add_subplot(gs[0])
                ax_table = fig.add_subplot(gs[1])
            else:
                ax = fig.add_subplot(111)

            cnv_df = _chromosome_cnv_dataframe(positions, values)

            # Highlight significant cytoband regions beneath the CNV track
            if significant_regions and contig in significant_regions:
                for region in significant_regions[contig]:
                    start_mb = region["start_pos"] / 1_000_000
                    end_mb = region["end_pos"] / 1_000_000
                    is_gain = region["type"] in ("GAIN", "HIGH_GAIN")
                    fill_color = CNV_COLORS["gain_fill"] if is_gain else CNV_COLORS["loss_fill"]
                    edge_color = CNV_COLORS["gain_edge"] if is_gain else CNV_COLORS["loss_edge"]

                    ax.axvspan(start_mb, end_mb, color=fill_color, alpha=0.55, zorder=0)
                    mid_point = (start_mb + end_mb) / 2
                    y_pos = y_max - (0.08 * (y_max - y_min))
                    ax.text(
                        mid_point,
                        y_pos,
                        region.get("name", ""),
                        fontsize=7,
                        fontweight="bold",
                        ha="center",
                        va="center",
                        color=edge_color,
                        bbox=dict(
                            boxstyle="round,pad=0.25",
                            fc="white",
                            ec=edge_color,
                            alpha=0.92,
                            linewidth=0.9,
                        ),
                        zorder=4,
                    )

            _add_cnv_reference_lines(ax, float(mean_cnv), float(std_cnv), y_min, y_max)
            _plot_cnv_track(ax, cnv_df, point_size=4, point_alpha=0.28, show_trend=False)

            gene_count, off_scale_count = _add_panel_gene_lollipops(
                ax,
                panel_points,
                y_min,
                y_max,
                off_scale_mode,
                float(mean_cnv),
            )

            _apply_cnv_chromosome_axes(
                ax,
                x_max_mb,
                y_max,
                xlabel="Position (Mb)",
                ylabel="Estimated ploidy",
            )

            if gene_count_on_chrom:
                _add_panel_gene_table(
                    ax_table,
                    panel_points,
                    off_scale_mode,
                    y_max,
                    float(mean_cnv),
                )

            if off_scale_count:
                off_scale_labels = [
                    f"{point['label']} ({point['cnv_val']:.1f})"
                    for point in panel_points
                    if point["cnv_val"] > y_max * 0.97
                ]
                ax.text(
                    0.01,
                    0.98,
                    "Amplified (off-scale): " + ", ".join(off_scale_labels),
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=6.5,
                    color=AMPLIFIED_COLOR,
                    fontweight="bold",
                    bbox=dict(
                        boxstyle="round,pad=0.35",
                        fc="#FFEBEE",
                        ec=AMPLIFIED_COLOR,
                        linewidth=1.0,
                        alpha=0.95,
                    ),
                    zorder=9,
                )

            try:
                buf = io.BytesIO()
                fig.savefig(buf, format="jpg", dpi=300, bbox_inches="tight", pad_inches=0.06)
                plt.close(fig)
                buf.seek(0)

                # Validate buffer contains data
                if not buf.getvalue():
                    logger.warning(f"Empty buffer for chromosome {contig} CNV plot")
                    continue

                # Validate it's a valid JPEG by checking magic bytes
                buf_data = buf.getvalue()
                if len(buf_data) < 2 or buf_data[:2] != b'\xff\xd8':
                    logger.warning(f"Invalid JPEG data for chromosome {contig} CNV plot")
                    continue

                buf.seek(0)
                plots.append((contig, buf))
            except Exception as e:
                logger.error(f"Error creating CNV plot for chromosome {contig}: {str(e)}")
                plt.close('all')
                continue

    except Exception as e:
        logger.error(f"Error in create_CNV_plot_per_chromosome: {str(e)}")
        plt.close('all')

    return plots


def classification_plot(df, title, threshold):
    """
    Creates a classification plot.

    Args:
        df (pd.DataFrame): DataFrame containing the classification data.
        title (str): Title of the plot.
        threshold (float): Threshold value for filtering classifications.

    Returns:
        io.BytesIO: Buffer containing the plot image.
    """
    set_modern_style()

    df["timestamp"] = pd.to_datetime(df["timestamp"], unit="ms", utc=True)

    # Reshape the data to long format
    df_melted = df.melt(id_vars=["timestamp"], var_name="Condition", value_name="Value")
    df_melted = df_melted[df_melted["Condition"].ne("number_probes")]

    # Filter conditions that cross the threshold
    top_conditions = df_melted.groupby("Condition")["Value"].max().nlargest(10).index
    df_filtered = df_melted[df_melted["Condition"].isin(top_conditions)]

    conditions_above_threshold = df_filtered[df_filtered["Value"] > threshold][
        "Condition"
    ].unique()
    df_filtered = df_filtered[df_filtered["Condition"].isin(conditions_above_threshold)]

    # Create figure with adjusted size and margins
    fig = plt.figure(figsize=(10, 6))

    # Only create legend if we have data to plot
    if not df_filtered.empty:
        sns.lineplot(
            data=df_filtered, x="timestamp", y="Value", hue="Condition", palette="Set2"
        )

        # Move the legend below the plot with adjusted position
        plt.legend(
            title="Condition", bbox_to_anchor=(0.5, -0.3), loc="upper center", ncol=3
        )
    else:
        # If no data, create an empty plot
        plt.plot([])
        plt.text(
            0.5,
            0.5,
            "No classification data above threshold",
            horizontalalignment="center",
            verticalalignment="center",
        )

    plt.title(f"{title} Classifications over Time")
    plt.xlabel("Timestamp")
    plt.ylabel("Value")
    plt.xticks(rotation=45)

    # Format the x-axis with custom date format
    ax = plt.gca()
    ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter("%Y-%m-%d %H:%M"))

    # Adjust layout with explicit margins
    plt.subplots_adjust(bottom=0.25, left=0.1, right=0.9, top=0.9)

    # Save the plot as a JPG file with reduced DPI
    buf = io.BytesIO()
    plt.savefig(buf, format="jpg", dpi=300, bbox_inches="tight")
    plt.close(fig)  # Close the figure to free memory
    buf.seek(0)
    return buf


def coverage_plot(df):
    """
    Creates a coverage plot.

    Args:
        df (pd.DataFrame): DataFrame containing coverage data.

    Returns:
        io.BytesIO: Buffer containing the plot image.
    """
    set_modern_style()

    # df = df[df["#rname"] != "chrM"].copy()
    df = df[
        df["#rname"].isin(["chr" + str(i) for i in range(0, 23)] + ["chrX", "chrY"])
    ].copy()

    # Sort chromosomes naturally
    df["#rname"] = pd.Categorical(
        df["#rname"], categories=natsort.natsorted(df["#rname"].unique()), ordered=True
    )
    df = df.sort_values("#rname")

    # Create figure with adjusted size
    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1], hspace=0.4)

    # Font size settings
    title_fontsize = 8
    label_fontsize = 8
    tick_fontsize = 6

    # Plot number of reads per chromosome
    ax0 = plt.subplot(gs[0])
    sns.barplot(
        x="#rname", y="numreads", data=df, ax=ax0, color=MODERN_COLORS["accent"]
    )
    ax0.set_title("Number of Reads per Chromosome", fontsize=title_fontsize)
    ax0.set_xlabel("", fontsize=label_fontsize)
    ax0.set_ylabel("Number of Reads", fontsize=label_fontsize)
    ax0.tick_params(axis="x", rotation=90, labelsize=tick_fontsize)
    ax0.tick_params(axis="y", labelsize=tick_fontsize)

    # Plot number of bases per chromosome
    ax1 = plt.subplot(gs[1])
    sns.barplot(
        x="#rname", y="covbases", data=df, ax=ax1, color=MODERN_COLORS["accent"]
    )
    ax1.set_title("Number of Bases per Chromosome", fontsize=title_fontsize)
    ax1.set_xlabel("", fontsize=label_fontsize)
    ax1.set_ylabel("Number of Bases", fontsize=label_fontsize)
    ax1.tick_params(axis="x", rotation=90, labelsize=tick_fontsize)
    ax1.tick_params(axis="y", labelsize=tick_fontsize)

    # Plot mean depth per chromosome
    ax2 = plt.subplot(gs[2])
    sns.barplot(
        x="#rname", y="meandepth", data=df, ax=ax2, color=MODERN_COLORS["accent"]
    )
    ax2.set_title("Mean Depth per Chromosome", fontsize=title_fontsize)
    ax2.set_xlabel("Chromosome", fontsize=label_fontsize)
    ax2.set_ylabel("Mean Depth", fontsize=label_fontsize)
    ax2.tick_params(axis="x", rotation=90, labelsize=tick_fontsize)
    ax2.tick_params(axis="y", labelsize=tick_fontsize)

    # Adjust layout with explicit margins
    plt.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95, hspace=0.5)

    # Save the plot as a JPG file with reduced DPI
    buf = io.BytesIO()
    plt.savefig(buf, format="jpg", dpi=300, bbox_inches="tight")
    plt.close(fig)  # Close the figure to free memory
    buf.seek(0)
    return buf


def plot_classification_timeline(
    df: pd.DataFrame, classification_level: str = "class", title: Optional[str] = None
) -> Tuple[plt.Figure, plt.Axes]:
    """Plot classification changes over time.

    Args:
        df: DataFrame containing classification data
        classification_level: Level of classification to plot ('superfamily', 'family', 'class', 'subclass')
        title: Optional title for the plot

    Returns:
        Tuple of (figure, axes) objects
    """
    if df.empty:
        logger.warning("No data available for plotting")
        return None, None

    # Create figure and axes
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot classification scores
    ax.plot(
        df["timestamp"],
        df[f"{classification_level}_score"],
        marker="o",
        linestyle="-",
        label="Score",
    )

    # Add labels for each point
    for x, y, label in zip(
        df["timestamp"],
        df[f"{classification_level}_score"],
        df[f"{classification_level}_label"],
    ):
        ax.annotate(label, (x, y), xytext=(5, 5), textcoords="offset points")

    # Customize plot
    ax.set_xlabel("Time")
    ax.set_ylabel("Classification Score")
    ax.set_title(title or f"{classification_level.title()} Classification Over Time")
    ax.grid(True, alpha=0.3)

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45)

    # Adjust layout
    plt.tight_layout()

    return fig, ax


def plot_mgmt_timeline(
    df: pd.DataFrame, title: Optional[str] = None
) -> Tuple[plt.Figure, plt.Axes]:
    """Plot MGMT methylation changes over time.

    Args:
        df: DataFrame containing MGMT data
        title: Optional title for the plot

    Returns:
        Tuple of (figure, axes) objects
    """
    if df.empty:
        logger.warning("No data available for plotting")
        return None, None

    # Create figure and axes
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot methylation percentage
    ax.plot(
        df["timestamp"],
        df["mgmt_methylation"],
        marker="o",
        linestyle="-",
        label="Methylation %",
    )

    # Add status labels
    for x, y, status in zip(df["timestamp"], df["mgmt_methylation"], df["mgmt_status"]):
        ax.annotate(status, (x, y), xytext=(5, 5), textcoords="offset points")

    # Customize plot
    ax.set_xlabel("Time")
    ax.set_ylabel("Methylation Percentage")
    ax.set_title(title or "MGMT Methylation Over Time")
    ax.grid(True, alpha=0.3)

    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45)

    # Adjust layout
    plt.tight_layout()

    return fig, ax


def save_plot(fig: plt.Figure, output_path: str, dpi: int = 300):
    """Save a plot to a file.

    Args:
        fig: Figure object to save
        output_path: Path to save the plot
        dpi: DPI for the output image
    """
    try:
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        # Save plot
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)

    except Exception as e:
        logger.error(f"Error saving plot to {output_path}: {str(e)}")


def plot_to_bytes(fig: plt.Figure, format: str = "png", dpi: int = 300) -> bytes:
    """Convert a plot to bytes.

    Args:
        fig: Figure object to convert
        format: Output format ('png', 'jpg', etc.)
        dpi: DPI for the output image

    Returns:
        Bytes containing the plot image
    """
    try:
        buf = io.BytesIO()
        fig.savefig(buf, format=format, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        return buf.getvalue()

    except Exception as e:
        logger.error(f"Error converting plot to bytes: {str(e)}")
        return None
