"""
report.py

This module contains the main function for creating the PDF report.
"""

import os
import io
import pickle
import numpy as np
import pandas as pd
from PIL import Image as PILImage
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import inch
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer,
    Image,
    Table,
    TableStyle,
    PageBreak,
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
    for style_name in styles.byName:
        styles[style_name].fontName = "FiraSans"

    smaller_style = ParagraphStyle(name="Smaller", parent=styles["Normal"], fontSize=8)
    bold_style = ParagraphStyle(
        name="Bold", parent=styles["Normal"], fontName="FiraSans-Bold"
    )
    underline_style = ParagraphStyle(
        name="Underline",
        parent=styles["Heading1"],
        fontName="FiraSans-Bold",
        underline=True,
    )

    styles.add(smaller_style)
    styles.add(bold_style)
    styles.add(underline_style)

    doc = SimpleDocTemplate(filename, pagesize=A4)
    elements_summary = []
    elements = []

    masterdf = (
        pd.read_csv(os.path.join(output, "master.csv"), index_col=0, header=None)
        if os.path.exists(os.path.join(output, "master.csv"))
        else None
    )

    elements_summary.append(Paragraph("Classification Summary", styles["Heading1"]))
    elements_summary.append(Spacer(1, 12))
    elements_summary.append(
        Paragraph("This sample has the following classifications:", styles["BodyText"])
    )
    elements_summary.append(Spacer(1, 12))

    threshold = 0.05


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
            file_path = os.path.join(output, df_name)

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
                    f"{name} classification: {lastrow_plot_top.index[0]} - {lastrow_plot_top.values[0]:.2f}",
                    styles["Bold"],
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
            table = Table(data)
            style = TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, -1), colors.white),
                    ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),
                    ("TEXTCOLOR", (0, 1), (-1, -1), colors.black),
                    ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                    ("FONTNAME", (0, 0), (-1, 0), "FiraSans-Bold"),
                    ("FONTNAME", (0, 1), (-1, -1), "FiraSans"),
                    ("FONTSIZE", (0, 0), (-1, 0), 6),
                    ("FONTSIZE", (0, 1), (-1, -1), 5),
                    ("BOTTOMPADDING", (0, 0), (-1, 0), 1),
                    ("BACKGROUND", (0, 1), (-1, -1), colors.white),
                    ("GRID", (0, 0), (-1, -1), 1, colors.black),
                ]
            )
            table.setStyle(style)
            elements.append(table)
            elements.append(Spacer(1, 12))
            elements.append(PageBreak())
        else:
            elements.append(
                Paragraph(f"No {name} Classification Available", styles["BodyText"])
            )

    # Add CNV plots
    if os.path.exists(os.path.join(output, "CNV.npy")):
        CNVresult = np.load(os.path.join(output, "CNV.npy"), allow_pickle="TRUE").item()
        CNVresult = Result(CNVresult)
        cnv_dict = np.load(
            os.path.join(output, "CNV_dict.npy"), allow_pickle=True
        ).item()
        file = open(os.path.join(output, "XYestimate.pkl"), "rb")
        XYestimate = pickle.load(file)
        elements_summary.append(
            Paragraph(
                f"Estimated sex chromosome composition: {XYestimate}",
                styles["BodyText"],
            )
        )

        cnv_summary = create_CNV_plot(CNVresult, cnv_dict)
        img = Image(cnv_summary, width=6 * inch, height=1.5 * inch)
        elements.append(img)
        elements.append(Spacer(1, 6))

        cnv_plots = create_CNV_plot_per_chromosome(CNVresult, cnv_dict)
        for contig, img_buf in cnv_plots:
            img = Image(img_buf, width=6 * inch, height=1.5 * inch)
            elements.append(img)
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

    # Add coverage plots
    if os.path.exists(os.path.join(output, "coverage_main.csv")):
        cov_df_main = pd.read_csv(os.path.join(output, "coverage_main.csv"))
        bedcov_df_main = pd.read_csv(os.path.join(output, "bed_coverage_main.csv"))
        target_coverage_df = pd.read_csv(os.path.join(output, "target_coverage.csv"))
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

        elements.append(Paragraph("Target Coverage", styles["Underline"]))
        elements.append(
            Paragraph("This plot was generated by ROBIN.", styles["Smaller"])
        )
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
        table = Table(data)
        style = TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), colors.white),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),
                ("TEXTCOLOR", (0, 1), (-1, -1), colors.black),
                ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                ("FONTNAME", (0, 0), (-1, 0), "FiraSans-Bold"),
                ("FONTNAME", (0, 1), (-1, -1), "FiraSans"),
                ("FONTSIZE", (0, 0), (-1, 0), 6),
                ("FONTSIZE", (0, 1), (-1, -1), 5),
                ("BOTTOMPADDING", (0, 0), (-1, 0), 1),
                ("BACKGROUND", (0, 1), (-1, -1), colors.white),
                ("GRID", (0, 0), (-1, -1), 1, colors.black),
            ]
        )
        table.setStyle(style)
        elements.append(table)
        elements.append(Spacer(1, 12))
    else:
        elements.append(Paragraph("No Coverage Data Available", styles["BodyText"]))

    elements.append(PageBreak())

    fusion_results = []
    if os.path.exists(os.path.join(output, "fusion_candidates_master.csv")):
        try:
            fusion_candidates = pd.read_csv(
                os.path.join(output, "fusion_candidates_master.csv"),
                dtype=str,
                # names=["index", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, "diff"],
                # names=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, "diff"],
                header=None,
                skiprows=1,
            )
            fusion_results.append(("High", fusion_candidates))
        except pd.errors.EmptyDataError:
            pass

        # result, goodpairs = _annotate_results(fusion_candidates)
        try:
            fusion_candidates_all = pd.read_csv(
                os.path.join(output, "fusion_candidates_all.csv"),
                dtype=str,
                # names=["index", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, "diff"],
                # names=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, "diff"],
                header=None,
                skiprows=1,
            )
            fusion_results.append(("Low", fusion_candidates_all))
        except pd.errors.EmptyDataError:
            pass
        # result_all, goodpairs_all = _annotate_results(fusion_candidates_all)
        if len(fusion_results) > 0:
            for confidence, df in fusion_results:
                result, goodpairs = _annotate_results(df)
                if not result.empty:
                    logger.info(result)
                    elements.append(
                        Paragraph(
                            f"Fusion Candidates - {confidence} Confidence",
                            styles["Underline"],
                        )
                    )
                    elements.append(
                        Paragraph(
                            "This plot was generated by ROBIN.", styles["Smaller"]
                        )
                    )
                    gene_pairs = (
                        result[goodpairs].sort_values(by=7)["tag"].unique().tolist()
                    )
                    stripped_list = [item.replace(" ", "") for item in gene_pairs]

                    gene_pairs = [pair.split(",") for pair in stripped_list]

                    gene_groups_test = get_gene_network(gene_pairs)

                    gene_groups = []
                    for gene_group in gene_groups_test:
                        if (
                            len(
                                _get_reads(
                                    result[goodpairs][
                                        result[goodpairs][3].isin(gene_group)
                                    ]
                                )
                            )
                            > 1
                        ):
                            gene_groups.append(gene_group)

                    logger.info(f"{confidence}, {len(gene_groups)}")

                    for gene_group in gene_groups:
                        elements.append(Paragraph(f"{gene_group}", styles["Smaller"]))
                        reads = result[goodpairs][result[goodpairs][3].isin(gene_group)]

                        result_reads = _get_reads(reads)
                        df = reads

                        def rank_overlaps(df, start_col, end_col):
                            # Sort by start and end columns
                            df = df.sort_values(
                                by=["gene", start_col, end_col]
                            ).reset_index(drop=True)
                            ranks = []
                            current_rank = 0
                            current_end = -1

                            for _, row in df.iterrows():
                                if row[start_col] > current_end:
                                    current_rank = 0
                                ranks.append(current_rank)
                                current_end = max(current_end, row[end_col])
                                current_rank += 1

                            return ranks

                        df["start2"] = df["start2"].astype(int)
                        df["end2"] = df["end2"].astype(int)
                        df = df.sort_values(by=["gene", "start2", "end2"])

                        df["rank"] = rank_overlaps(df, "start2", "end2")

                        lines = df.sort_values(by=["id", "read_start"]).reset_index(
                            drop=True
                        )

                        # Function to assign occurrence number
                        def assign_occurrence(group):
                            group["Occurrence"] = range(1, len(group) + 1)
                            return group

                        # Apply the function to each group
                        lines = (
                            lines.groupby("id")
                            .apply(assign_occurrence)
                            .reset_index(drop=True)
                        )

                        # Function to find join coordinates
                        def find_joins(group):
                            group = group.reset_index(drop=True)
                            for i in range(len(group)):
                                if i < len(group) - 1:
                                    if (
                                        group.loc[i, "Occurrence"] == 1
                                    ):  # The first read in the sequence of read mappings
                                        group.loc[i, "Join_Gene"] = group.loc[
                                            i + 1, "gene"
                                        ]
                                        group.loc[i, "Join_Start"] = group.loc[
                                            i, "start2"
                                        ]  # Reference end of current read
                                        group.loc[i, "Join_Chromosome"] = group.loc[
                                            i + 1, "chromosome2"
                                        ]
                                        group.loc[i, "spanB"] = group.loc[i + 1, "span"]
                                        group.loc[i, "start3"] = group.loc[
                                            i + 1, "start2"
                                        ]
                                        group.loc[i, "end3"] = group.loc[i + 1, "end2"]
                                        group.loc[i, "rankB"] = group.loc[i + 1, "rank"]
                                        if group.loc[i + 1, "strand"] == "+":

                                            group.loc[i, "Join_End"] = group.loc[
                                                i + 1, "end2"
                                            ]  # Reference start of next read
                                        else:

                                            group.loc[i, "Join_End"] = group.loc[
                                                i + 1, "start2"
                                            ]  # Reference end of next read

                            return group

                        # Initialize columns for the coordinates where the read joins the next read
                        lines["Join_Start"] = None
                        lines["Join_End"] = None
                        lines["Join_Chromosome"] = None
                        lines["Join_Gene"] = None
                        lines["spanB"] = None
                        lines["start3"] = None
                        lines["end3"] = None
                        lines["rankB"] = None

                        # Apply the function to each group
                        lines = (
                            lines.groupby("id").apply(find_joins).reset_index(drop=True)
                        )
                        # Remove rows containing NA values
                        lines = lines.dropna()

                        datafile = "rCNS2_data.csv.gz"

                        if os.path.isfile(
                            os.path.join(
                                os.path.dirname(os.path.abspath(resources.__file__)),
                                datafile,
                            )
                        ):
                            gene_table = pd.read_csv(
                                os.path.join(
                                    os.path.dirname(
                                        os.path.abspath(resources.__file__)
                                    ),
                                    datafile,
                                )
                            )

                        width, height = A4

                        axdict = {}
                        if len(result_reads) > 1:
                            plt.figure(figsize=(12, 5))

                            num_plots = 2 * len(result_reads)
                            num_cols = len(
                                result_reads
                            )  # Number of columns in the subplot grid
                            num_rows = (
                                num_plots + num_cols - 1
                            ) // num_cols  # Calculate the number of rows needed

                            for i, ax in enumerate(range(num_plots), start=1):
                                plt.subplot(num_rows, num_cols, i)
                                row, col = divmod(i - 1, num_cols)
                                data = result_reads.iloc[col]

                                chrom = data["chromosome"]
                                start = data["start"]
                                end = data["end"]

                                def human_readable_format(x, pos):
                                    return f"{x / 1e6:.1f}"  # Mb

                                if row == 1:
                                    features = []
                                    for index, row in gene_table[
                                        gene_table["Seqid"].eq(chrom)
                                        & gene_table["Start"].le(end)
                                        & gene_table["End"].ge(start)
                                    ].iterrows():
                                        if row["Type"] == "gene":
                                            features.append(
                                                GraphicFeature(
                                                    start=int(row["Start"]),
                                                    end=int(row["End"]),
                                                    strand=STRAND[row["Strand"]],
                                                    thickness=8,
                                                    color="#ffd700",
                                                    label=row["gene_name"],
                                                    fontdict={
                                                        "family": "serif",
                                                        "color": "black",
                                                        "fontsize": 8,
                                                    },
                                                    # fontdict={"color": "black", "fontsize": 8},
                                                )
                                            )

                                    for index, row in (
                                        gene_table[
                                            gene_table["gene_name"].eq(data["gene"])
                                            & gene_table["Source"].eq("HAVANA")
                                            & gene_table["Type"].eq("exon")
                                        ]
                                        .groupby(
                                            ["Seqid", "Start", "End", "Type", "Strand"]
                                        )
                                        .count()
                                        .reset_index()
                                        .iterrows()
                                    ):
                                        features.append(
                                            GraphicFeature(
                                                start=int(row["Start"]),
                                                end=int(row["End"]),
                                                strand=STRAND[row["Strand"]],
                                                thickness=4,
                                                color="#C0C0C0",
                                            )
                                        )

                                    record = GraphicRecord(
                                        sequence_length=end - start,
                                        first_index=start,
                                        features=features,
                                    )
                                    ax = plt.gca()

                                    record.plot(
                                        ax=ax,
                                        with_ruler=False,
                                        draw_line=True,
                                        strand_in_label_threshold=4,
                                    )

                                else:
                                    features2 = []
                                    for index, row in (
                                        df[df["chromosome"].eq(chrom)]
                                        .sort_values(by="id")
                                        .iterrows()
                                    ):
                                        features2.append(
                                            GraphicFeature(
                                                start=int(row["start2"]),
                                                end=int(row["end2"]),
                                                strand=STRAND[row["strand"]],
                                                color=row["color"],
                                            )
                                        )

                                    record2 = GraphicRecord(
                                        sequence_length=end - start,
                                        first_index=start,
                                        features=features2,
                                    )
                                    ax = plt.gca()
                                    record2.plot(ax=ax)
                                    ax.xaxis.set_major_formatter(
                                        FuncFormatter(human_readable_format)
                                    )
                                    ax.tick_params(
                                        axis="x", labelsize=8
                                    )  # Set font size for x-axis labels
                                    ax.set_xlabel(
                                        f'Position (Mb) - {chrom} - {data["gene"]}',
                                        fontsize=10,
                                    )
                                    ax.set_title(
                                        f'{data["gene"]}', fontsize=10
                                    )  # Reduced font size for titles
                                    ax.set_xticklabels(
                                        ax.get_xticklabels(), rotation=90
                                    )

                                    ax.grid(which="both", axis="x", linestyle="")

                                    axdict[data["gene"]] = ax

                            # Adjust the space between subplots
                            plt.subplots_adjust(hspace=0.4)

                            gene_counter = Counter()

                            for index, row in lines.iterrows():
                                # Coordinates in data space of each subplot
                                # xyB = (
                                #    row["Join_Start"],
                                #    row["rank"],
                                # )  # Point in subplot 1
                                # xyA = (row["Join_End"], row["rankB"])

                                gene_counter[row["gene"]] += 1
                                gene_counter[row["Join_Gene"]] += 1
                                """
                                if row["Join_Gene"] in axdict.keys():
                                    # Create an arc connection
                                    con = ConnectionPatch(
                                        xyA=xyA,
                                        coordsA=axdict[row["Join_Gene"]].transData,
                                        xyB=xyB,
                                        coordsB=axdict[row["gene"]].transData,
                                        axesB=row["gene"],
                                        axesA=row["Join_Gene"],
                                        connectionstyle="arc3,rad=0.05",
                                        arrowstyle="->",
                                        linestyle="--",
                                        color=row["color"],
                                    )
                                """

                            buf2 = io.BytesIO()
                            plt.savefig(
                                buf2, format="jpg", dpi=300, bbox_inches="tight"
                            )
                            plt.close()
                            buf2.seek(0)

                            img = Image(
                                buf2,
                                width=width * 0.95,
                                height=width,
                                kind="proportional",
                            )
                            elements.append(img)
                            elements.append(Spacer(1, 12))

                    elements_summary.append(
                        Paragraph(
                            f"Fusion Candidates - using panel rCNS2 - {len(gene_groups)} {confidence} confidence fusions observed.",
                            styles["BodyText"],
                        )
                    )

                elements.append(PageBreak())

    else:
        elements_summary.append(
            Paragraph("No Fusion Data Available", styles["BodyText"])
        )

    # Add run data summary
    if masterdf is not None and isinstance(masterdf, pd.DataFrame):
        elements_summary.append(Paragraph("Run Data Summary", styles["Heading2"]))
        # start_time = "Placeholder"
        masterdf_dict = eval(masterdf[masterdf.index == "samples"][1]["samples"])[
            sample_id
        ]
        elements_summary.append(
            Paragraph(
                f"Sample ID: {sample_id}<br/>"
                f"Run Start: {convert_to_space_separated_string(masterdf_dict['run_time'])}<br/>"
                f"Run Folder: {' '.join(masterdf.loc[(masterdf.index == 'watchfolder')][1].values)}<br/>"
                f"Output Folder: {' '.join(masterdf.loc[(masterdf.index == 'output')][1].values)}<br/>"
                f"Target Panel: {' '.join(masterdf.loc[(masterdf.index == 'target_panel')][1].values)}<br/>"
                f"Reference: {' '.join(masterdf.loc[(masterdf.index == 'reference')][1].values)}<br/>"
                f"Sequencing Device: {convert_to_space_separated_string(masterdf_dict['devices'])}<br/>"
                f"Flowcell ID: {convert_to_space_separated_string(masterdf_dict['flowcell_ids'])}<br/>"
                f"Basecalling Model: {convert_to_space_separated_string(masterdf_dict['basecall_models'])}<br/>",
                styles["Smaller"],
            )
        )
        try:
            formatted_lines = "".join(
                (
                    f"{k}: {v}<br/>"
                    for k, v in eval(
                        masterdf[masterdf.index == "samples"][1]["samples"]
                    )[sample_id]["file_counters"].items()
                )
            )
            elements_summary.append(Paragraph(f"{formatted_lines}", styles["Smaller"]))
        except Exception as e:
            logger.info(f"Error parsing file counters: {e}")

    elements.append(Spacer(1, 12))
    elements.append(PageBreak())

    last_seen = 0
    for file in natsort.natsorted(os.listdir(output)):
        if file.endswith("_mgmt.csv"):
            count = int(file.split("_")[0])
            if count > last_seen:
                results = pd.read_csv(os.path.join(output, file))
                plot_out = os.path.join(output, file.replace(".csv", ".png"))
                last_seen = count

    if last_seen > 0:
        elements.append(Paragraph("MGMT Promoter Methylation", styles["Underline"]))
        image = Image(plot_out, 6 * inch, 4 * inch)
        elements.append(image)
        data_list = [results.columns.values.tolist()] + results.values.tolist()
        table = Table(data_list)
        style = TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), colors.white),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),
                ("TEXTCOLOR", (0, 1), (-1, -1), colors.black),
                ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                ("FONTNAME", (0, 0), (-1, 0), "FiraSans-Bold"),
                ("FONTNAME", (0, 1), (-1, -1), "FiraSans"),
                ("FONTSIZE", (0, 0), (-1, 0), 6),
                ("FONTSIZE", (0, 1), (-1, -1), 5),
                ("BOTTOMPADDING", (0, 0), (-1, 0), 1),
                ("BACKGROUND", (0, 1), (-1, -1), colors.white),
                ("GRID", (0, 0), (-1, -1), 1, colors.black),
            ]
        )
        table.setStyle(style)
        elements.append(table)

    final_elements = elements_summary + elements
    doc.multiBuild(
        final_elements,
        canvasmaker=header_footer_canvas_factory(sample_id, styles, fonts_dir),
    )
    logger.info(f"PDF created: {filename}")
    return filename
