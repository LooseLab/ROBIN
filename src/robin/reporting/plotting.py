"""
plotting.py

This module contains functions for creating plots used in the PDF report.
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import io
import math
import matplotlib.font_manager as fm
import os
from robin import fonts

import natsort

from matplotlib import gridspec
from robin.subpages.CNV_object import filter_and_find_max

# Define consistent color scheme and style
MODERN_COLORS = {
    "primary": "#2C3E50",  # Dark blue-grey (matching report text)
    "secondary": "#E2E8F0",  # Light grey (matching table grid)
    "background": "#F8FAFC",  # Light background (matching table alternate rows)
    "accent": "#3498DB",  # Blue accent
    "grid": "#E2E8F0",  # Grid color
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


def create_CNV_plot(result, cnv_dict):
    """
    Creates a CNV plot.

    Args:
        result (Result): CNV result object.
        cnv_dict (dict): Dictionary containing CNV data.

    Returns:
        io.BytesIO: Buffer containing the plot image.
    """
    set_modern_style()

    # Prepare data for plotting
    plot_data = []
    offset = 0
    contig_centers = {}

    for contig, values in result.cnv.items():
        if contig in ["chr" + str(i) for i in range(0, 23)] + ["chrX", "chrY"]:
            start_offset = offset
            for i, value in enumerate(values):
                plot_data.append((contig, i + offset, value))
            end_offset = offset + len(values) - 1
            contig_centers[contig] = (
                start_offset + end_offset
            ) / 2  # Calculate central position
            offset += len(values)  # Increase the offset for the next contig

    # Convert to DataFrame
    df = pd.DataFrame(plot_data, columns=["Contig", "Position", "Value"])
    df["Position_Corrected"] = df["Position"] * cnv_dict["bin_width"]

    # Calculate the mean and standard deviation of the 'Value' column
    mean_value = df["Value"].mean()
    std_value = df["Value"].std()

    # Calculate the threshold
    threshold = mean_value + 5 * std_value

    # Create scatter plot
    width = 16
    plt.figure(figsize=(width, width / 4))
    sns.scatterplot(
        data=df,
        x="Position_Corrected",
        y="Value",
        hue="Contig",
        palette="Set2",  # Modern color palette
        legend=False,
        s=2,
        alpha=0.6,
    )

    min_y = df["Value"].min()  # Minimum y-value
    label_y_position = 4.5  # Position labels below the minimum y-value

    for contig, center in contig_centers.items():
        plt.text(
            center * cnv_dict["bin_width"],
            label_y_position,
            contig,
            fontsize=12,
            ha="center",
            va="top",
            rotation=45,
        )

    plt.ylim(min_y, threshold)  # Set the y-axis limits
    plt.title("Copy Number Changes")
    plt.xlabel("Position")
    plt.ylabel("Estimated ploidy")
    buf = io.BytesIO()
    plt.savefig(buf, format="jpg", dpi=300, bbox_inches="tight")
    buf.seek(0)
    return buf


def create_CNV_plot_per_chromosome(result, cnv_dict, significant_regions=None):
    """Creates CNV plots per chromosome.

    Args:
        result (Result): CNV result object.
        cnv_dict (dict): Dictionary containing CNV data.
        significant_regions (dict): Dictionary mapping chromosomes to lists of significant regions.
            Each region should be a dict with keys: 'start_pos', 'end_pos', 'type' ('GAIN' or 'LOSS')

    Returns:
        List[Tuple[str, io.BytesIO]]: List of tuples containing chromosome names and plot buffers.
    """
    set_modern_style()

    plots = []
    for contig, values in result.cnv.items():
        if contig in ["chr" + str(i) for i in range(0, 23)] + ["chrX", "chrY"]:
            # Calculate mean and standard deviation for this chromosome
            values_array = np.array(values)
            mean_cnv = np.mean(values_array)
            std_cnv = np.std(values_array)
            
            # Set y-axis limits based on mean and standard deviation
            y_min = 0
            y_max = mean_cnv + (2 * std_cnv)
            
            # Calculate positions in megabases
            positions = np.arange(len(values)) * cnv_dict["bin_width"] / 1_000_000  # Convert to Mb
            
            plt.figure(figsize=(4, 2))
            
            # If we have significant regions for this chromosome, highlight them
            if significant_regions and contig in significant_regions:
                for region in significant_regions[contig]:
                    start_mb = region['start_pos'] / 1_000_000  # Convert to Mb
                    end_mb = region['end_pos'] / 1_000_000  # Convert to Mb
                    
                    # Choose color based on type
                    color = '#e8f5e9' if region['type'] == 'GAIN' else '#ffebee'  # Light green for gains, light red for losses
                    alpha = 0.3
                    
                    # Add shaded region
                    plt.axvspan(start_mb, end_mb, color=color, alpha=alpha, zorder=1)
                    
                    # Add cytoband label with background
                    mid_point = (start_mb + end_mb) / 2
                    y_pos = y_max - (0.1 * (y_max - y_min))  # Position near the top
                    
                    # Create background box for text
                    bbox_props = dict(
                        boxstyle="round,pad=0.3",
                        fc='white',
                        ec=color.replace('e8', 'a5').replace('ff', 'ef'),  # Darker version of highlight color
                        alpha=0.8
                    )
                    
                    # Add text with background
                    plt.text(
                        mid_point, 
                        y_pos,
                        region.get('name', ''),  # Add cytoband name if available
                        fontsize=8,
                        fontweight='bold',
                        ha='center',
                        va='center',
                        bbox=bbox_props,
                        zorder=4
                    )

            # Plot CNV data points on top of highlighted regions
            sns.scatterplot(
                x=positions,
                y=values,
                s=2,
                color=MODERN_COLORS["accent"],
                alpha=0.6,
                zorder=2
            )

            plt.xlabel("Position (Mb)")
            plt.ylabel("Estimated Ploidy")
            plt.ylim(y_min, y_max)

            # Add horizontal lines for mean and standard deviations
            plt.axhline(y=mean_cnv, color='gray', linestyle='--', alpha=0.5, zorder=3)
            plt.axhline(y=mean_cnv + std_cnv, color='gray', linestyle=':', alpha=0.3, zorder=3)
            plt.axhline(y=mean_cnv + (2 * std_cnv), color='gray', linestyle=':', alpha=0.3, zorder=3)

            buf = io.BytesIO()
            plt.savefig(buf, format="jpg", dpi=300, bbox_inches="tight")
            buf.seek(0)
            plots.append((contig, buf))
            plt.close()

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
