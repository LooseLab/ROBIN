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
from matplotlib.transforms import blended_transform_factory

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
    "points": "#5B7C9D",
    "median": "#1E3A5F",
    "diploid": "#3D5A80",
    "reference": "#C5CDD6",
    "gain_fill": "#D4EDDA",
    "loss_fill": "#F8D7DA",
    "gain_edge": "#3D7A4A",
    "loss_edge": "#A94442",
    "gene": "#5C6BC0",
}

CNV_FONT = {
    "title": 10,
    "subtitle": 7.5,
    "axis": 9,
    "tick": 7,
    "annotation": 6.5,
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
    show_xlabel: bool = True,
) -> None:
    """Pin axes at the origin with no leading x padding."""
    ax.set_xlim(0, x_max_mb)
    ax.set_ylim(0, y_max)
    ax.margins(x=0, y=0)
    ax.autoscale(enable=False)
    if show_xlabel:
        ax.set_xlabel(xlabel, fontsize=CNV_FONT["axis"], color=MODERN_COLORS["primary"])
    else:
        ax.set_xlabel("")
    ax.set_ylabel(ylabel, fontsize=CNV_FONT["axis"], color=MODERN_COLORS["primary"])
    ax.spines["left"].set_position(("data", 0))
    ax.spines["bottom"].set_position(("data", 0))
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    ax.grid(True, axis="y", color=MODERN_COLORS["grid"], linestyle="--", linewidth=0.4, alpha=0.55)
    ax.grid(False, axis="x")
    ax.tick_params(colors=MODERN_COLORS["primary"], labelsize=CNV_FONT["tick"])


def _chromosome_display_name(contig: str) -> str:
    return contig.replace("chr", "")


def _plot_cnv_track(ax, cnv_df: pd.DataFrame, x_max_mb: float) -> None:
    """Render CNV bin values as scatter with a rolling median trace."""
    ax.scatter(
        cnv_df["position_mb"],
        cnv_df["ploidy"],
        s=4,
        c=CNV_COLORS["points"],
        alpha=0.35,
        linewidths=0,
        edgecolors="none",
        rasterized=True,
        zorder=2,
    )
    window = max(5, len(cnv_df) // 35)
    if len(cnv_df) >= window:
        median = cnv_df["ploidy"].rolling(window, center=True, min_periods=1).median()
        ax.plot(
            cnv_df["position_mb"],
            median,
            color=CNV_COLORS["median"],
            linewidth=1.4,
            alpha=0.95,
            solid_capstyle="round",
            zorder=4,
        )


def _is_gene_cnv_outlier(cnv_val: float, mean_cnv: float, std_cnv: float) -> bool:
    """Return True when a panel target copy number deviates >2 SD from the chromosome mean."""
    if std_cnv < 1e-6:
        return abs(cnv_val - mean_cnv) > 0.5
    return abs(cnv_val - mean_cnv) > 2 * std_cnv


SIGNIFICANT_CNV_REGION_TYPES = {"GAIN", "LOSS", "HIGH_GAIN", "DEEP_LOSS"}


def _gene_overlaps_cnv_region(mid_bp: float, regions: List[Dict[str, Any]]) -> bool:
    """Return True when a gene midpoint falls inside a called gain/loss region."""
    for region in regions:
        if region.get("type") not in SIGNIFICANT_CNV_REGION_TYPES:
            continue
        start_bp = float(region["start_pos"])
        end_bp = float(region["end_pos"])
        if start_bp <= mid_bp <= end_bp:
            return True
    return False


def _should_label_panel_gene(
    point: Dict[str, Any],
    mean_cnv: float,
    std_cnv: float,
    regions: List[Dict[str, Any]],
) -> bool:
    """Label outliers and panel genes inside called CNV regions."""
    if _is_gene_cnv_outlier(point["cnv_val"], mean_cnv, std_cnv):
        return True
    return _gene_overlaps_cnv_region(point["mid_mb"] * 1_000_000, regions)


def _is_highlighted_panel_gene(
    point: Dict[str, Any],
    mean_cnv: float,
    std_cnv: float,
    regions: List[Dict[str, Any]],
    off_scale_mode: bool,
    y_max: float,
) -> tuple[bool, bool]:
    """Return (highlight, off_scale) for a panel target."""
    cnv_val = point["cnv_val"]
    off_scale = off_scale_mode and cnv_val > y_max * 0.97
    highlight = off_scale or _should_label_panel_gene(point, mean_cnv, std_cnv, regions)
    return highlight, off_scale


def _add_cnv_regions_on_plot(ax, regions: List[Dict[str, Any]], y_max: float) -> None:
    """Shade called regions on the main track with edge guides and top brackets."""
    for region in regions:
        start_mb = region["start_pos"] / 1_000_000
        end_mb = region["end_pos"] / 1_000_000
        is_gain = region["type"] in ("GAIN", "HIGH_GAIN")
        fill_color = CNV_COLORS["gain_fill"] if is_gain else CNV_COLORS["loss_fill"]
        edge_color = CNV_COLORS["gain_edge"] if is_gain else CNV_COLORS["loss_edge"]

        ax.axvspan(start_mb, end_mb, color=fill_color, alpha=0.55, zorder=0, linewidth=0)
        ax.axvline(start_mb, color=edge_color, linestyle="--", linewidth=0.9, alpha=0.75, zorder=1)
        ax.axvline(end_mb, color=edge_color, linestyle="--", linewidth=0.9, alpha=0.75, zorder=1)
        _draw_region_bracket(
            ax, start_mb, end_mb, y_max * 0.992, edge_color, height_frac=y_max * 0.028,
        )


LOLLIPOP_LABEL_TIERS = (-0.040, -0.062, -0.084)


def _label_spacing_mb(label: str, x_max_mb: float, rotation: int) -> float:
    """Estimate horizontal label footprint in megabases."""
    char_factor = 0.0035 if rotation == 0 else 0.0024
    return max(x_max_mb * 0.011, len(label) * x_max_mb * char_factor, 1.0)


def _layout_gene_lollipop_labels(
    panel_points: List[Dict[str, Any]],
    x_max_mb: float,
    *,
    rotation: int,
) -> List[Dict[str, Any]]:
    """Assign below-axis label tiers to reduce overlap."""
    occupied: dict[int, list[tuple[float, str]]] = {
        tier: [] for tier in range(len(LOLLIPOP_LABEL_TIERS))
    }
    layouts: List[Dict[str, Any]] = []

    for point in sorted(panel_points, key=lambda item: item["mid_mb"]):
        label = point["label"]
        label_y = LOLLIPOP_LABEL_TIERS[0]
        tier = 0
        spacing = _label_spacing_mb(label, x_max_mb, rotation)
        placed = False

        for tier_idx, candidate_y in enumerate(LOLLIPOP_LABEL_TIERS):
            if all(
                abs(point["mid_mb"] - other_mb) >= max(
                    spacing,
                    _label_spacing_mb(other_label, x_max_mb, rotation),
                )
                for other_mb, other_label in occupied[tier_idx]
            ):
                tier = tier_idx
                label_y = candidate_y
                occupied[tier_idx].append((point["mid_mb"], label))
                placed = True
                break
        if not placed:
            occupied[0].append((point["mid_mb"], label))

        layouts.append({**point, "tier": tier, "label_y": label_y})

    return layouts


def _add_panel_gene_lollipops(
    ax,
    panel_points: List[Dict[str, Any]],
    x_max_mb: float,
    y_max: float,
    off_scale_mode: bool,
    mean_cnv: float,
    std_cnv: float,
    regions: List[Dict[str, Any]],
) -> None:
    """Draw panel lollipops; label outliers and genes in called CNV regions."""
    if not panel_points:
        return

    label_points = [
        point for point in panel_points
        if _should_label_panel_gene(point, mean_cnv, std_cnv, regions)
    ]
    rotation = 0 if len(label_points) <= 8 else 45
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    label_layouts = {
        (layout["label"], layout["mid_mb"]): layout
        for layout in _layout_gene_lollipop_labels(label_points, x_max_mb, rotation=rotation)
    }

    for point in panel_points:
        mid_mb = point["mid_mb"]
        cnv_val = point["cnv_val"]
        highlight, off_scale = _is_highlighted_panel_gene(
            point, mean_cnv, std_cnv, regions, off_scale_mode, y_max,
        )
        color = AMPLIFIED_COLOR if highlight else CNV_COLORS["gene"]
        y_top = min(cnv_val, y_max * 0.98)
        stem_alpha = 0.9 if highlight else 0.45
        stem_width = 1.0 if highlight else 0.6

        ax.plot(
            [mid_mb, mid_mb],
            [0.0, y_top],
            color=color,
            linewidth=stem_width,
            alpha=stem_alpha,
            zorder=5,
            solid_capstyle="round",
        )

        if off_scale:
            ax.scatter(
                [mid_mb],
                [y_max * 0.99],
                s=34,
                marker="^",
                color=AMPLIFIED_COLOR,
                zorder=7,
                edgecolors="none",
                clip_on=False,
            )
        else:
            ax.scatter(
                [mid_mb],
                [y_top],
                s=24 if highlight else 12,
                color=color,
                zorder=6,
                edgecolors="white",
                linewidths=0.35,
                alpha=0.95 if highlight else 0.7,
            )

        layout = label_layouts.get((point["label"], mid_mb))
        if layout is None:
            continue

        label = point["label"]
        if off_scale:
            label = f"{label} ({cnv_val:.1f})"
        ax.text(
            mid_mb,
            layout["label_y"],
            label,
            transform=trans,
            fontsize=CNV_FONT["annotation"],
            color=color,
            ha="center",
            va="top",
            rotation=rotation,
            rotation_mode="anchor",
            zorder=8,
            clip_on=False,
        )


def _add_cnv_ploidy_reference_lines(ax, y_max: float, x_max_mb: float) -> None:
    """Add dashed reference lines at integer ploidy values."""
    for ploidy in (1, 2, 3):
        if ploidy > y_max:
            continue
        is_diploid = ploidy == 2
        ax.axhline(
            ploidy,
            color=CNV_COLORS["diploid"] if is_diploid else CNV_COLORS["reference"],
            linestyle="--",
            linewidth=1.3 if is_diploid else 0.7,
            alpha=0.95 if is_diploid else 0.55,
            zorder=1,
        )
        ax.text(
            x_max_mb * 0.008,
            ploidy,
            str(ploidy),
            fontsize=CNV_FONT["annotation"],
            color=CNV_COLORS["diploid"] if is_diploid else CNV_COLORS["reference"],
            va="center",
            ha="left",
            zorder=1,
        )


def _draw_region_bracket(
    ax,
    start_mb: float,
    end_mb: float,
    y_level: float,
    color: str,
    height_frac: float = 0.04,
) -> None:
    """Draw a horizontal genomic bracket."""
    y_top = y_level
    y_drop = y_level - height_frac
    ax.plot(
        [start_mb, start_mb, end_mb, end_mb],
        [y_drop, y_top, y_top, y_drop],
        color=color,
        linewidth=1.1,
        solid_capstyle="butt",
        zorder=5,
        clip_on=False,
    )


def _apply_cnv_axes_style(ax, *, xlabel: str, ylabel: str, title: Optional[str] = None) -> None:
    """Apply consistent seaborn-inspired styling to a CNV axes."""
    ax.set_facecolor("white")
    ax.set_xlabel(xlabel, fontsize=CNV_FONT["axis"], color=MODERN_COLORS["primary"])
    ax.set_ylabel(ylabel, fontsize=CNV_FONT["axis"], color=MODERN_COLORS["primary"])
    if title:
        ax.set_title(title, fontsize=CNV_FONT["title"], color=MODERN_COLORS["primary"], pad=8)
    sns.despine(ax=ax, top=True, right=True)
    ax.grid(True, axis="y", color=MODERN_COLORS["grid"], linestyle="--", linewidth=0.4, alpha=0.55)
    ax.grid(False, axis="x")
    ax.tick_params(colors=MODERN_COLORS["primary"], labelsize=CNV_FONT["tick"])


def _chromosome_cnv_dataframe(positions_mb, values) -> pd.DataFrame:
    return pd.DataFrame(
        {"position_mb": positions_mb, "ploidy": pd.Series(values, dtype=float)}
    )


def _add_cnv_reference_lines(ax, mean_cnv: float, std_cnv: float, y_min: float, y_max: float) -> None:
    """Genome-wide reference guides (mean plus optional spread)."""
    ax.axhline(y=mean_cnv, color=CNV_COLORS["reference"], linestyle="--", linewidth=0.9, alpha=0.7, zorder=1)
    for offset in (std_cnv, 2 * std_cnv):
        if y_min < mean_cnv + offset <= y_max:
            ax.axhline(y=mean_cnv + offset, color=CNV_COLORS["reference"], linestyle=":", linewidth=0.6, alpha=0.4, zorder=1)
        if y_min <= mean_cnv - offset < y_max:
            ax.axhline(y=mean_cnv - offset, color=CNV_COLORS["reference"], linestyle=":", linewidth=0.6, alpha=0.4, zorder=1)


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

    merged: Dict[str, Dict[str, Any]] = {}
    for point in points:
        existing = merged.get(point["label"])
        if existing is None or point["cnv_val"] > existing["cnv_val"]:
            merged[point["label"]] = point

    return list(merged.values())


def _compute_cnv_y_limits(
    values_array: np.ndarray,
    panel_points: List[Dict[str, Any]],
) -> tuple[float, float, bool, float, float]:
    """
    Choose y-axis limits that preserve bulk CNV detail while surfacing
    amplified panel targets when feasible.
    """
    mean_cnv = float(np.mean(values_array))
    std_cnv = float(np.std(values_array))
    baseline_max = max(mean_cnv + (2 * std_cnv), mean_cnv * 1.4, 2.5)

    if not panel_points:
        return 0.0, baseline_max, False, mean_cnv, std_cnv

    panel_max = max(point["cnv_val"] for point in panel_points)
    if panel_max <= baseline_max * 1.05:
        return 0.0, baseline_max, False, mean_cnv, std_cnv

    if panel_max <= baseline_max * PANEL_Y_EXPANSION_FACTOR:
        return 0.0, panel_max * 1.08, False, mean_cnv, std_cnv

    return 0.0, baseline_max, True, mean_cnv, std_cnv


def create_CNV_plot_per_chromosome(
    result,
    cnv_dict,
    significant_regions=None,
    chromosomes: Optional[List[str]] = None,
    panel_genes_df: Optional[pd.DataFrame] = None,
    chromosome_status: Optional[Dict[str, str]] = None,
):
    """Creates CNV plots per chromosome.

    Args:
        result (Result): CNV result object.
        cnv_dict (dict): Dictionary containing CNV data.
        significant_regions (dict): Dictionary mapping chromosomes to lists of significant regions.
        chromosomes (list, optional): Ordered chromosome names to plot.
        panel_genes_df (pd.DataFrame, optional): Target panel genes.
        chromosome_status (dict, optional): Status text keyed by chromosome.

    Returns:
        List[Tuple[str, io.BytesIO]]: List of tuples containing chromosome names and plot buffers.
    """
    plots = []
    chromosome_status = chromosome_status or {}
    try:
        set_modern_style()

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

            values_array = np.array(values)
            mean_cnv = float(np.mean(values_array))
            std_cnv = float(np.std(values_array))
            positions = np.arange(len(values)) * cnv_dict["bin_width"] / 1_000_000
            x_max_mb = float(positions[-1]) if len(positions) else 1.0

            panel_points = _collect_panel_gene_points(
                panel_genes_df, contig, values, cnv_dict["bin_width"],
            )
            y_min, y_max, off_scale_mode, mean_cnv, std_cnv = _compute_cnv_y_limits(
                values_array, panel_points,
            )
            y_min = 0.0
            regions = (significant_regions or {}).get(contig, [])
            chrom_name = _chromosome_display_name(contig)
            status_text = chromosome_status.get(contig, "No significant CNV change")

            label_count = sum(
                1 for point in panel_points
                if _should_label_panel_gene(point, mean_cnv, std_cnv, regions)
            )
            total_height = CNV_CHROMOSOME_FIG_WIDTH * 0.26
            fig = plt.figure(figsize=(CNV_CHROMOSOME_FIG_WIDTH, total_height))
            ax = fig.add_subplot(111)

            bottom_margin = 0.10 if label_count == 0 else min(0.13, 0.095 + label_count * 0.004)
            fig.subplots_adjust(top=0.86, bottom=bottom_margin)
            fig.suptitle(
                f"Chromosome {chrom_name} copy-number profile",
                fontsize=CNV_FONT["title"],
                fontweight="bold",
                color=MODERN_COLORS["primary"],
                y=0.98,
            )
            if status_text and status_text != "No significant CNV change":
                fig.text(
                    0.5,
                    0.91,
                    status_text,
                    ha="center",
                    va="top",
                    fontsize=CNV_FONT["subtitle"],
                    color=MODERN_COLORS["primary"],
                )

            cnv_df = _chromosome_cnv_dataframe(positions, values)
            if regions:
                _add_cnv_regions_on_plot(ax, regions, y_max)
            _add_cnv_ploidy_reference_lines(ax, y_max, x_max_mb)
            _plot_cnv_track(ax, cnv_df, x_max_mb)
            _apply_cnv_chromosome_axes(
                ax,
                x_max_mb,
                y_max,
                xlabel=f"Position on chromosome {chrom_name} (Mb)",
                ylabel="Estimated copy number / ploidy",
                show_xlabel=True,
            )
            _add_panel_gene_lollipops(
                ax, panel_points, x_max_mb, y_max, off_scale_mode, mean_cnv, std_cnv, regions,
            )

            off_scale_outliers = [
                point for point in panel_points
                if _should_label_panel_gene(point, mean_cnv, std_cnv, regions)
                and off_scale_mode
                and point["cnv_val"] > y_max * 0.97
            ]
            if off_scale_outliers:
                amp_text = ", ".join(
                    f"{point['label']} ({point['cnv_val']:.1f})" for point in off_scale_outliers
                )
                ax.text(
                    0.01,
                    0.99,
                    f"Off-scale: {amp_text}",
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=CNV_FONT["annotation"],
                    color=AMPLIFIED_COLOR,
                    bbox=dict(
                        boxstyle="round,pad=0.25",
                        facecolor="white",
                        edgecolor=MODERN_COLORS["grid"],
                        alpha=0.92,
                    ),
                    zorder=9,
                )

            try:
                buf = io.BytesIO()
                fig.savefig(buf, format="jpg", dpi=300, bbox_inches="tight", pad_inches=0.08)
                plt.close(fig)
                buf.seek(0)

                if not buf.getvalue():
                    logger.warning(f"Empty buffer for chromosome {contig} CNV plot")
                    continue

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
