from reportlab.lib.pagesizes import A4
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas
from reportlab.platypus import (
    SimpleDocTemplate,
    Table,
    TableStyle,
    Paragraph,
    Image,
    Spacer,
    PageBreak,
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase import pdfmetrics
import matplotlib

matplotlib.use("agg")
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties

import natsort
import seaborn as sns
import pandas as pd
import numpy as np
import io
import os
from datetime import datetime
import random
import pickle
from PIL import Image as PILImage
from typing import Tuple

from robin.subpages.CNV_object import Result

from robin import fonts
from robin import images
from robin import resources
from dna_features_viewer import GraphicFeature, GraphicRecord

STRAND = {"+": 1, "-": -1}

pdfmetrics.registerFont(
    TTFont(
        "FiraSans",
        os.path.join(
            os.path.dirname(os.path.abspath(fonts.__file__)),
            "fira-sans-v16-latin-regular.ttf",
        ),
    )
)
pdfmetrics.registerFont(
    TTFont(
        "FiraSans-Bold",
        os.path.join(
            os.path.dirname(os.path.abspath(fonts.__file__)),
            "fira-sans-v16-latin-700.ttf",
        ),
    )
)  # Assuming this is the path for the bold version

# Update styles to use the custom font
styles = getSampleStyleSheet()
for style_name in styles.byName:
    styles[style_name].fontName = "FiraSans"

# Define a smaller style for the header date
smaller_style = ParagraphStyle(name="Smaller", parent=styles["Normal"], fontSize=8)

# Define a bold style for the first header line
bold_style = ParagraphStyle(
    name="Bold", parent=styles["Normal"], fontName="FiraSans-Bold"
)

# Define an underlined style for section headings
underline_style = ParagraphStyle(
    name="Underline", parent=styles["Heading1"], underline=True
)


def convert_to_space_separated_string(array):
    try:
        import ast

        # Convert array to list and extract the string
        string_repr = array.tolist()[0]

        # Evaluate the string to convert it to an actual list
        list_repr = ast.literal_eval(string_repr)

        # Join the elements of the list into a space-separated string
        return " ".join(list_repr)
    except:
        return array


class HeaderFooterCanvas(canvas.Canvas):

    def __init__(self, sample_id, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.pages = []
        self.sample_id = sample_id  # Store the sample_id

    def showPage(self):
        self.pages.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        page_count = len(self.pages)
        for page in self.pages:
            self.__dict__.update(page)
            self.draw_canvas(page_count)
            canvas.Canvas.showPage(self)
        canvas.Canvas.save(self)

    def draw_canvas(self, page_count):
        width, height = A4

        # Add first line of the header in bold
        header1 = Paragraph("R.O.B.I.N Reports... RESEARCH USE ONLY", bold_style)
        w, h = header1.wrap(width - 2 * inch, inch)
        header1.drawOn(self, inch, height - h - inch + 36)

        # Add second line of the header
        header2 = Paragraph(f"Sample ID: {self.sample_id}", styles["Normal"])
        w, h = header2.wrap(width - 2 * inch, inch)
        header2.drawOn(self, inch, height - h - inch + 24)

        # Add third line of the header with date and time
        report_date = datetime.now().strftime(
            "R.O.B.I.N. Report Generated: %Y-%m-%d %H:%M:%S"
        )
        header3 = Paragraph(report_date, smaller_style)
        w, h = header3.wrap(width - 2 * inch, inch)
        header3.drawOn(self, inch, height - h - inch + 12)

        # Add logo to the top right corner of the header
        logo_path = os.path.join(
            os.path.dirname(os.path.abspath(images.__file__)), "ROBIN_logo_small.png"
        )
        # logo_path = "src/robin/images/Robin_logo_small.png"  # Replace with the path to your logo
        max_logo_size = 50  # Maximum width and height in pixels
        self.drawImage(
            logo_path,
            width - max_logo_size - inch,
            height - max_logo_size - inch + 36,
            width=max_logo_size,
            height=max_logo_size,
            preserveAspectRatio=True,
            mask="auto",
        )

        # Add footer
        page = f"SampleID: {self.sample_id} - Page {self._pageNumber} of {page_count}"
        x = 190
        self.saveState()
        self.setStrokeColorRGB(0, 0, 0)
        self.setLineWidth(0.5)
        self.line(66, 78, A4[0] - 66, 78)
        self.setFont("FiraSans", 7)
        self.drawString(A4[0] - x, 65, page)
        self.restoreState()


def header_footer_canvas_factory(sample_id):
    def create_canvas(*args, **kwargs):
        return HeaderFooterCanvas(sample_id, *args, **kwargs)

    return create_canvas


def target_distribution_plot(df):
    df["chrom"] = pd.Categorical(
        df["chrom"], categories=natsort.natsorted(df["chrom"].unique()), ordered=True
    )
    df = df.sort_values("chrom")

    # Generate the plot
    plt.figure(figsize=(16, 8))
    boxplot = sns.boxplot(x="chrom", y="coverage", data=df)
    plt.title("Distribution of Target Coverage on Each Chromosome")
    plt.xlabel("Chromosome")
    plt.ylabel("Coverage")
    plt.xticks(rotation=90)

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


def fusion_plot_test(title, reads):
    features = [
        GraphicFeature(
            start=0, end=20, strand=+1, color="#ffd700", label="Small feature"
        ),
        GraphicFeature(
            start=20,
            end=500,
            strand=+1,
            color="#ffcccc",
            label="Gene 1 with a very long name",
        ),
        GraphicFeature(start=400, end=700, strand=-1, color="#cffccc", label="Gene 2"),
        GraphicFeature(start=600, end=900, strand=+1, color="#ccccff", label="Gene 3"),
    ]
    record = GraphicRecord(sequence_length=1000, features=features)
    ax, _ = record.plot(figure_width=5)
    ax.figure.savefig("graphic_record_defined_by_hand.png")


def fusion_plot(title, reads):
    # ToDo: Remove hardcoding to this gene file
    datafile = "rCNS2_data.csv.gz"

    if os.path.isfile(
        os.path.join(os.path.dirname(os.path.abspath(resources.__file__)), datafile)
    ):
        gene_table = pd.read_csv(
            os.path.join(os.path.dirname(os.path.abspath(resources.__file__)), datafile)
        )

    plt.figure(figsize=(16, 2))
    ax1 = plt.gca()
    features = []
    first_index = 0
    sequence_length = 100
    x_label = ""
    for index, row in gene_table[gene_table["gene_name"].eq(title.strip())].iterrows():
        if row["Type"] == "gene":
            x_label = row["Seqid"]
            features.append(
                GraphicFeature(
                    start=int(row["Start"]),
                    end=int(row["End"]),
                    strand=STRAND[row["Strand"]],
                    thickness=4,
                    color="#ffd700",
                )
            )
            first_index = int(row["Start"]) - 1000
            sequence_length = int(row["End"]) - int(row["Start"]) + 2000
        if row["Type"] == "CDS" and row["transcript_type"] == "protein_coding":
            features.append(
                GraphicFeature(
                    start=int(row["Start"]),
                    end=int(row["End"]),
                    strand=STRAND[row["Strand"]],
                    color="#ffcccc",
                )
            )
    record = GraphicRecord(
        sequence_length=sequence_length, first_index=first_index, features=features
    )
    ax1.set_title(f"{title} - {x_label}")
    ax, _ = record.plot(ax=ax1, figure_width=5)
    buf = io.BytesIO()
    ax.figure.savefig(buf, format="jpg", dpi=300, bbox_inches="tight")
    buf.seek(0)
    plt.close()
    plt.figure(figsize=(16, 2))
    ax = plt.gca()
    features = []
    for index, row in reads.sort_values(by=7).iterrows():
        features.append(
            GraphicFeature(
                start=int(row[5]),
                end=int(row[6]),
                strand=STRAND[row[9]],
                color=row["Color"],
            )
        )
    record = GraphicRecord(
        sequence_length=sequence_length, first_index=first_index, features=features
    )
    ax, _ = record.plot(
        ax=ax, figure_width=5, with_ruler=False, draw_line=True
    )  # , figure_width=5)
    buf2 = io.BytesIO()
    ax.figure.savefig(buf2, format="jpg", dpi=300, bbox_inches="tight")
    plt.close()
    buf2.seek(0)
    return buf, buf2


def get_target_outliers(df):
    df["chrom"] = pd.Categorical(
        df["chrom"], categories=natsort.natsorted(df["chrom"].unique()), ordered=True
    )
    df = df.sort_values("chrom")
    Q1 = df["coverage"].quantile(0.25)
    Q3 = df["coverage"].quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1 * IQR
    upper_bound = Q3 + 1 * IQR
    outliers = df[(df["coverage"] < lower_bound) | (df["coverage"] > upper_bound)]
    return outliers


def create_CNV_plot_per_chromosome(result, cnv_dict):
    plots = []
    for contig, values in result.cnv.items():
        if contig != "chrM":
            plt.figure(figsize=(10, 2))
            sns.scatterplot(x=range(len(values)), y=values, s=2)
            plt.title(f"Copy Number Variation - {contig}")
            plt.xlabel("Position")
            plt.ylabel("Estimated Ploidy")
            plt.ylim(0, None)  # Adjust as needed

            buf = io.BytesIO()
            plt.savefig(buf, format="jpg", dpi=300, bbox_inches="tight")
            buf.seek(0)
            plots.append((contig, buf))
            plt.close()

    return plots


def coverage_plot(df):
    df = df[df["#rname"] != "chrM"].copy()

    # Sort chromosomes naturally
    df["#rname"] = pd.Categorical(
        df["#rname"], categories=natsort.natsorted(df["#rname"].unique()), ordered=True
    )
    df = df.sort_values("#rname")

    # Create subplots
    fig = plt.figure(figsize=(10, 6))
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1])

    # Font size settings
    title_fontsize = 8
    label_fontsize = 8
    tick_fontsize = 6

    # Plot number of reads per chromosome
    ax0 = plt.subplot(gs[0])
    sns.barplot(x="#rname", y="numreads", data=df, ax=ax0)
    ax0.set_title("Number of Reads per Chromosome", fontsize=title_fontsize)
    ax0.set_xlabel("", fontsize=label_fontsize)
    ax0.set_ylabel("Number of Reads", fontsize=label_fontsize)
    ax0.tick_params(axis="x", rotation=90, labelsize=tick_fontsize)
    ax0.tick_params(axis="y", labelsize=tick_fontsize)

    # Plot number of bases per chromosome
    ax1 = plt.subplot(gs[1])
    sns.barplot(x="#rname", y="covbases", data=df, ax=ax1)
    ax1.set_title("Number of Bases per Chromosome", fontsize=title_fontsize)
    ax1.set_xlabel("", fontsize=label_fontsize)
    ax1.set_ylabel("Number of Bases", fontsize=label_fontsize)
    ax1.tick_params(axis="x", rotation=90, labelsize=tick_fontsize)
    ax1.tick_params(axis="y", labelsize=tick_fontsize)

    # Plot mean depth per chromosome
    ax2 = plt.subplot(gs[2])
    sns.barplot(x="#rname", y="meandepth", data=df, ax=ax2)
    ax2.set_title("Mean Depth per Chromosome", fontsize=title_fontsize)
    ax2.set_xlabel("Chromosome", fontsize=label_fontsize)
    ax2.set_ylabel("Mean Depth", fontsize=label_fontsize)
    ax2.tick_params(axis="x", rotation=90, labelsize=tick_fontsize)
    ax2.tick_params(axis="y", labelsize=tick_fontsize)

    plt.tight_layout()
    # Save the plot as a JPG file
    buf = io.BytesIO()
    plt.savefig(buf, format="jpg", dpi=300, bbox_inches="tight")
    buf.seek(0)
    return buf


def _generate_random_color() -> str:
    """
    Generates a random color for use in plotting.

    Returns:
        str: A random hex color code.
    """
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))


def _annotate_results(result: pd.DataFrame) -> Tuple[pd.DataFrame, pd.Series]:
    """
    Annotates the result DataFrame with tags and colors.

    Args:
        result (pd.DataFrame): DataFrame with fusion candidates.

    Returns:
        Tuple[pd.DataFrame, pd.Series]: Annotated DataFrame and a boolean Series indicating good pairs.
    """
    result_copy = result.copy()
    lookup = result_copy.groupby(7)[3].agg(lambda x: ",".join(set(x)))
    tags = result_copy[7].map(lookup.get)
    result_copy.loc[:, "tag"] = tags
    result = result_copy
    colors = result.groupby(7).apply(
        lambda x: _generate_random_color(), include_groups=False
    )
    result["Color"] = result[7].map(colors.get)
    goodpairs = result.groupby("tag")[7].transform("nunique") > 1
    return result, goodpairs


def classification_plot(df, title, threshold):
    df["timestamp"] = pd.to_datetime(df["timestamp"], unit="ms", utc=True)

    # Reshape the data to long format
    df_melted = df.melt(id_vars=["timestamp"], var_name="Condition", value_name="Value")
    df_melted = df_melted[df_melted["Condition"].ne("number_probes")]

    # Filter conditions that cross the threshold of 0.05
    threshold = threshold

    top_conditions = df_melted.groupby("Condition")["Value"].max().nlargest(10).index
    df_filtered = df_melted[df_melted["Condition"].isin(top_conditions)]

    conditions_above_threshold = df_filtered[df_filtered["Value"] > threshold][
        "Condition"
    ].unique()
    df_filtered = df_filtered[df_filtered["Condition"].isin(conditions_above_threshold)]

    sns.set_theme(style="whitegrid", palette="colorblind")
    # Create the timecourse plot
    plt.figure(figsize=(12, 9))
    sns.lineplot(data=df_filtered, x="timestamp", y="Value", hue="Condition")
    plt.title(f"{title} Classifications over Time")
    plt.xlabel("Timestamp")
    plt.ylabel("Value")
    plt.xticks(rotation=45)

    # Format the x-axis with custom date format
    ax = plt.gca()
    ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter("%Y-%m-%d %H:%M"))

    # Move the legend below the plot
    plt.legend(
        title="Condition", bbox_to_anchor=(0.5, -0.8), loc="upper center", ncol=3
    )

    plt.tight_layout()

    # Save the plot as a JPG file
    buf = io.BytesIO()
    plt.savefig(buf, format="jpg", dpi=300, bbox_inches="tight")
    buf.seek(0)
    return buf


def create_CNV_plot(result, cnv_dict):
    # Prepare data for plotting
    plot_data = []
    offset = 0
    contig_centers = {}

    for contig, values in result.cnv.items():
        if contig != "chrM":
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
    sns.set_theme(style="whitegrid", palette="colorblind")
    plt.figure(figsize=(width, width / 4))
    sns.scatterplot(
        data=df,
        x="Position_Corrected",
        y="Value",
        hue="Contig",
        palette="colorblind",
        legend=False,
        s=1,
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


def split_text(text):
    """Helper function to split text on '-' characters."""
    return text.replace("-", "-\n")


# Create a sample DataFrame
def create_dataframe():
    data = {"A": np.random.rand(5), "B": np.random.rand(5), "C": np.random.rand(5)}
    df = pd.DataFrame(data)
    return df


# Convert DataFrame to ReportLab Table
def dataframe_to_table(df):
    data = [df.columns.values.tolist()] + df.values.tolist()
    table = Table(data)
    return table


# Create PDF
def create_pdf(filename, output):
    if filename.startswith("None"):
        final_folder = os.path.basename(os.path.normpath(output))
        filename = filename.replace("None", final_folder, 1)
        sample_id = final_folder
    else:
        sample_id = os.path.basename(os.path.normpath(output))
    print(f"Creating PDF {filename} in {output} for sample {sample_id}")

    pdfmetrics.registerFont(
        TTFont(
            "FiraSans",
            os.path.join(
                os.path.dirname(os.path.abspath(fonts.__file__)),
                "fira-sans-v16-latin-regular.ttf",
            ),
        )
    )
    pdfmetrics.registerFont(
        TTFont(
            "FiraSans-Bold",
            os.path.join(
                os.path.dirname(os.path.abspath(fonts.__file__)),
                "fira-sans-v16-latin-700.ttf",
            ),
        )
    )

    # Update styles to use the custom font
    styles = getSampleStyleSheet()
    for style_name in styles.byName:
        styles[style_name].fontName = "FiraSans"

    # Define a smaller style for the header date
    smaller_style = ParagraphStyle(name="Smaller", parent=styles["Normal"], fontSize=8)

    # Define a bold style for the first header line
    bold_style = ParagraphStyle(
        name="Bold", parent=styles["Normal"], fontName="FiraSans-Bold"
    )

    # Define an underlined style for section headings
    underline_style = ParagraphStyle(
        name="Underline", parent=styles["Heading1"], underline=True
    )

    doc = SimpleDocTemplate(filename, pagesize=A4)
    elements = []

    if os.path.exists(os.path.join(output, "master.csv")):
        masterdf = pd.read_csv(
            os.path.join(output, "master.csv"), index_col=0, header=None
        )
    else:
        masterdf = None

    # if masterdf is not None and isinstance(masterdf, pd.DataFrame):
    #    sample_id = convert_to_space_separated_string(
    #        masterdf.loc[(masterdf.index == "sample_ids")][1].values
    #    )

    # else:
    #    sample_id = None

    elements.append(Paragraph("Classification Summary", styles["Heading1"]))

    elements.append(Spacer(1, 12))

    elements.append(
        Paragraph("This sample has the following classifications:", styles["BodyText"])
    )
    elements.append(Spacer(1, 12))
    threshold = 0.05
    if os.path.exists(os.path.join(output, "sturgeon_scores.csv")):
        sturgeon_df_store = pd.read_csv(
            os.path.join(os.path.join(output, "sturgeon_scores.csv")),
            # index_col=0,
        )
        sturgeon_df_store2 = sturgeon_df_store.drop(columns=["timestamp"])
        columns_greater_than_threshold = (sturgeon_df_store2 > threshold).any()
        columns_not_greater_than_threshold = ~columns_greater_than_threshold
        Sturgeonresult = sturgeon_df_store2.columns[
            columns_not_greater_than_threshold
        ].tolist()
        Sturgeonlastrow = sturgeon_df_store2.iloc[-1].drop("number_probes")
        Sturgeonlastrow_plot = Sturgeonlastrow.sort_values(ascending=False).head(10)
        Sturgeonlastrow_plot_top = Sturgeonlastrow.sort_values(ascending=False).head(1)
        elements.append(
            Paragraph(
                f"Sturgeon classification: {Sturgeonlastrow_plot_top.index[0]} - {Sturgeonlastrow_plot_top.values[0]:.2f}",
                bold_style,
            )
        )

    else:
        elements.append(
            Paragraph("No Sturgeon Classification Available", styles["BodyText"])
        )

    if os.path.exists(os.path.join(output, "nanoDX_scores.csv")):
        nanodx_df_store = pd.read_csv(
            os.path.join(os.path.join(output, "nanoDX_scores.csv")),
            # index_col=0,
        )

        nanodx_df_store2 = nanodx_df_store.drop(columns=["timestamp"])
        columns_greater_than_threshold = (nanodx_df_store2 > threshold).any()
        columns_not_greater_than_threshold = ~columns_greater_than_threshold
        result = nanodx_df_store2.columns[columns_not_greater_than_threshold].tolist()
        NanoDXlastrow = nanodx_df_store2.iloc[-1].drop("number_probes")
        NanoDXlastrow_plot = NanoDXlastrow.sort_values(ascending=False).head(10)
        NanoDXlastrow_plot_top = NanoDXlastrow.sort_values(ascending=False).head(1)
        elements.append(
            Paragraph(
                f"NanoDX classification: {NanoDXlastrow_plot_top.index[0]} - {NanoDXlastrow_plot_top.values[0]:.2f}",
                bold_style,
            )
        )
    else:
        elements.append(
            Paragraph("No NanoDX Classification Available", styles["BodyText"])
        )

    if os.path.exists(os.path.join(output, "random_forest_scores.csv")):
        rcns2_df_store = pd.read_csv(
            os.path.join(os.path.join(output, "random_forest_scores.csv")),
            # index_col=0,
        )
        rcns2_df_store2 = rcns2_df_store.drop(columns=["timestamp"])
        columns_greater_than_threshold = (rcns2_df_store2 > threshold).any()
        columns_not_greater_than_threshold = ~columns_greater_than_threshold
        Forestresult = rcns2_df_store2.columns[
            columns_not_greater_than_threshold
        ].tolist()
        Forestlastrow = rcns2_df_store2.iloc[-1]  # .drop("number_probes")
        Forestlastrow_plot = Forestlastrow.sort_values(ascending=False).head(10)
        Forestlastrow_plot_top = Forestlastrow.sort_values(ascending=False).head(1)
        elements.append(
            Paragraph(
                f"Forest classification: {Forestlastrow_plot_top.index[0]} - {Forestlastrow_plot_top.values[0]:.2f}",
                bold_style,
            )
        )
    else:
        elements.append(
            Paragraph("No Forest Classification Available", styles["BodyText"])
        )

    if os.path.exists(os.path.join(output, "CNV.npy")):
        CNVresult = np.load(os.path.join(output, "CNV.npy"), allow_pickle="TRUE").item()
        CNVresult = Result(CNVresult)
        cnv_dict = np.load(
            os.path.join(output, "CNV_dict.npy"), allow_pickle=True
        ).item()
        file = open(os.path.join(output, "XYestimate.pkl"), "rb")
        XYestimate = pickle.load(file)
        elements.append(
            Paragraph(
                f"Estimated sex chromosome composition: {XYestimate}",
                styles["BodyText"],
            )
        )

    if os.path.exists(os.path.join(output, "coverage_main.csv")):
        cov_df_main = pd.read_csv(os.path.join(output, "coverage_main.csv"))
        bedcov_df_main = pd.read_csv(os.path.join(output, "bed_coverage_main.csv"))
        target_coverage_df = pd.read_csv(os.path.join(output, "target_coverage.csv"))
        elements.append(
            Paragraph(
                f"Coverage Depths - Global Estimated Coverage: {(cov_df_main['covbases'].sum() / cov_df_main['endpos'].sum()):.2f}x Targets Estimated Coverage: {(bedcov_df_main['bases'].sum() / bedcov_df_main['length'].sum()):.2f}x",
                styles["BodyText"],
            )
        )
        if bedcov_df_main["bases"].sum() / bedcov_df_main["length"].sum() < 10:
            elements.append(
                Paragraph(
                    "Target Coverage is below the recommended 10x threshold",
                    styles["BodyText"],
                )
            )
        # Get the outliers
        outliers = get_target_outliers(target_coverage_df)

        # Round the coverage values to 5 decimal places
        outliers["coverage"] = outliers["coverage"].apply(lambda x: round(x, 1))

        # Sort the dataframe by coverage in descending order
        outliers = outliers.sort_values(by="coverage", ascending=False)

        # Define the threshold value
        threshold = bedcov_df_main["bases"].sum() / bedcov_df_main["length"].sum()

        # Split the dataframe based on the threshold value
        outliers_above_threshold = outliers[outliers["coverage"] > threshold].copy()
        outliers_below_threshold = outliers[outliers["coverage"] <= threshold].copy()

        # Create the desired string representation using f-strings for both groups
        outliers_above_threshold["name_with_coverage"] = outliers_above_threshold.apply(
            lambda row: f"{row['name']} ({row['coverage']})", axis=1
        )
        outliers_below_threshold["name_with_coverage"] = outliers_below_threshold.apply(
            lambda row: f"{row['name']} ({row['coverage']})", axis=1
        )

        gene_names = " - ".join(outliers_above_threshold["name_with_coverage"])
        elements.append(Spacer(1, 6))
        elements.append(
            Paragraph(f"Outlier genes by coverage (high): {gene_names}", smaller_style)
        )
        elements.append(Spacer(1, 6))
        gene_names = " - ".join(outliers_below_threshold["name_with_coverage"])
        elements.append(Spacer(1, 6))
        elements.append(
            Paragraph(f"Outlier genes by coverage (low): {gene_names}", smaller_style)
        )

    else:
        elements.append(Paragraph("No Coverage Data Available", styles["BodyText"]))
    # Add summary section

    if os.path.exists(os.path.join(output, "fusion_candidates_master.csv")):
        fusion_candidates = pd.read_csv(
            os.path.join(output, "fusion_candidates_master.csv"),
            dtype=str,
            # names=["index", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, "diff"],
            names=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, "diff"],
            header=None,
            skiprows=1,
        )

        result, goodpairs = _annotate_results(fusion_candidates)

        fusion_candidates_all = pd.read_csv(
            os.path.join(output, "fusion_candidates_all.csv"),
            dtype=str,
            # names=["index", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, "diff"],
            names=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, "diff"],
            header=None,
            skiprows=1,
        )

        result_all, goodpairs_all = _annotate_results(fusion_candidates_all)

        elements.append(
            Paragraph(
                f"Fusion Candidates - using panel rCNS2 <br/> {result[goodpairs].sort_values(by=7)['tag'].nunique()} high confidence fusions observed.<br/>{result_all[goodpairs_all].sort_values(by=7)['tag'].nunique()} low confidence fusions observed.",
                styles["BodyText"],
            )
        )
    else:
        elements.append(Paragraph("No Fusion Data Available", styles["BodyText"]))

    elements.append(Spacer(1, 12))
    if masterdf is not None and isinstance(masterdf, pd.DataFrame):
        elements.append(Paragraph("Run Data Summary", styles["Heading2"]))
        start_time = "Placeholder"
        masterdf_dict = eval(masterdf[masterdf.index == "samples"][1]["samples"])[
            sample_id
        ]
        elements.append(
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
                # f"Files Seen: {masterdf.loc[(masterdf.index == 'file_counters')][1].values}<br/>",
                smaller_style,
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
            elements.append(
                Paragraph(
                    f"{formatted_lines}",
                    smaller_style,
                )
            )
        except Exception as e:
            print(f"Error parsing file counters: {e}")

    elements.append(Spacer(1, 12))

    elements.append(PageBreak())

    elements.append(Paragraph("Methylation Classifications", styles["Heading2"]))

    if os.path.exists(os.path.join(output, "sturgeon_scores.csv")):
        elements.append(
            Paragraph(
                f"Sturgeon classification: {Sturgeonlastrow_plot_top.index[0]} - {Sturgeonlastrow_plot_top.values[0]:.2f}",
                styles["Heading3"],
            )
        )
        elements.append(
            Paragraph("This plot was generated by sturgeon.", smaller_style)
        )
        img_buf = classification_plot(sturgeon_df_store, "Sturgeon", 0.05)
        # Read the image from the buffer to get its dimensions
        img_pil = PILImage.open(img_buf)
        width_img, height_img = img_pil.size

        width, height = A4

        height = (width * 0.95) / width_img * height_img

        img = Image(img_buf, width=width * 0.95, height=height, kind="proportional")
        elements.append(img)

        elements.append(Spacer(1, 6))  # Adjust spacer to minimal height

        df = Sturgeonlastrow_plot.reset_index()
        df.columns = ["Classification", "Score"]
        df["Score"] = df["Score"].apply(lambda x: round(x, 5))
        # Transpose the DataFrame
        df_transposed = df.set_index("Classification").T.reset_index()

        # Split headers on '-' characters
        df_transposed.columns = [split_text(col) for col in df_transposed.columns]

        df_transposed.columns.name = None  # Remove index name
        data = [df_transposed.columns.to_list()] + df_transposed.values.tolist()
        table = Table(data)

        # Add style to the table
        style = TableStyle(
            [
                (
                    "BACKGROUND",
                    (0, 0),
                    (-1, -1),
                    colors.white,
                ),  # Set background to white
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),  # Header text color
                ("TEXTCOLOR", (0, 1), (-1, -1), colors.black),  # Body text color
                ("ALIGN", (0, 0), (-1, -1), "CENTER"),  # Center align all cells
                ("FONTNAME", (0, 0), (-1, 0), "FiraSans-Bold"),  # Header font
                ("FONTNAME", (0, 1), (-1, -1), "FiraSans"),  # Body font
                ("FONTSIZE", (0, 0), (-1, 0), 6),  # Header font size reduced
                ("FONTSIZE", (0, 1), (-1, -1), 5),  # Body font size reduced
                ("BOTTOMPADDING", (0, 0), (-1, 0), 1),  # Header padding
                ("BACKGROUND", (0, 1), (-1, -1), colors.white),  # Body background white
                ("GRID", (0, 0), (-1, -1), 1, colors.black),  # Grid lines
            ]
        )
        table.setStyle(style)

        elements.append(table)
        # elements.append(Spacer(1, 12))
        # elements.append(PageBreak())

    if os.path.exists(os.path.join(output, "nanoDX_scores.csv")):
        elements.append(
            Paragraph(
                f"NanoDX (crossNN) classification: {NanoDXlastrow_plot_top.index[0]} - {NanoDXlastrow_plot_top.values[0]:.2f}",
                styles["Heading3"],
            )
        )
        elements.append(Paragraph("This plot was generated by NanoDX.", smaller_style))
        img_buf = classification_plot(nanodx_df_store, "NanoDX", 0.05)
        # Read the image from the buffer to get its dimensions
        img_pil = PILImage.open(img_buf)
        width_img, height_img = img_pil.size

        width, height = A4

        height = (width * 0.95) / width_img * height_img

        img = Image(img_buf, width=width * 0.95, height=height, kind="proportional")
        elements.append(img)
        df = NanoDXlastrow_plot.reset_index()
        df.columns = ["Classification", "Score"]
        df["Score"] = df["Score"].apply(lambda x: round(x, 5))

        # Transpose the DataFrame
        df_transposed = df.set_index("Classification").T.reset_index()
        df_transposed.columns.name = None  # Remove index name
        data = [df_transposed.columns.to_list()] + df_transposed.values.tolist()
        table = Table(data)

        # Add style to the table
        style = TableStyle(
            [
                (
                    "BACKGROUND",
                    (0, 0),
                    (-1, -1),
                    colors.white,
                ),  # Set background to white
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),  # Header text color
                ("TEXTCOLOR", (0, 1), (-1, -1), colors.black),  # Body text color
                ("ALIGN", (0, 0), (-1, -1), "CENTER"),  # Center align all cells
                ("FONTNAME", (0, 0), (-1, 0), "FiraSans-Bold"),  # Header font
                ("FONTNAME", (0, 1), (-1, -1), "FiraSans"),  # Body font
                ("FONTSIZE", (0, 0), (-1, 0), 6),  # Header font size reduced
                ("FONTSIZE", (0, 1), (-1, -1), 5),  # Body font size reduced
                ("BOTTOMPADDING", (0, 0), (-1, 0), 1),  # Header padding
                ("BACKGROUND", (0, 1), (-1, -1), colors.white),  # Body background white
                ("GRID", (0, 0), (-1, -1), 1, colors.black),  # Grid lines
            ]
        )

        table.setStyle(style)

        elements.append(table)
        # elements.append(Spacer(1, 12))
        # elements.append(PageBreak())

    if os.path.exists(os.path.join(output, "random_forest_scores.csv")):
        elements.append(
            Paragraph(
                f"Random Forest classification: {Forestlastrow_plot_top.index[0]} - {Forestlastrow_plot_top.values[0]:.2f}",
                styles["Heading3"],
            )
        )
        elements.append(
            Paragraph("This plot was generated by Random Forest.", smaller_style)
        )
        img_buf = classification_plot(rcns2_df_store, "Forest", 0.05)
        # Read the image from the buffer to get its dimensions
        img_pil = PILImage.open(img_buf)
        width_img, height_img = img_pil.size

        width, height = A4

        height = (width * 0.95) / width_img * height_img

        img = Image(img_buf, width=width * 0.95, height=height, kind="proportional")
        elements.append(img)
        df = Forestlastrow_plot.reset_index()
        df.columns = ["Classification", "Score"]
        df["Score"] = df["Score"].apply(lambda x: round(x, 5))
        # Transpose the DataFrame
        df_transposed = df.set_index("Classification").T.reset_index()
        df_transposed.columns.name = None  # Remove index name
        data = [df_transposed.columns.to_list()] + df_transposed.values.tolist()
        table = Table(data)

        # Add style to the table
        style = TableStyle(
            [
                (
                    "BACKGROUND",
                    (0, 0),
                    (-1, -1),
                    colors.white,
                ),  # Set background to white
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),  # Header text color
                ("TEXTCOLOR", (0, 1), (-1, -1), colors.black),  # Body text color
                ("ALIGN", (0, 0), (-1, -1), "CENTER"),  # Center align all cells
                ("FONTNAME", (0, 0), (-1, 0), "FiraSans-Bold"),  # Header font
                ("FONTNAME", (0, 1), (-1, -1), "FiraSans"),  # Body font
                ("FONTSIZE", (0, 0), (-1, 0), 6),  # Header font size reduced
                ("FONTSIZE", (0, 1), (-1, -1), 5),  # Body font size reduced
                ("BOTTOMPADDING", (0, 0), (-1, 0), 1),  # Header padding
                ("BACKGROUND", (0, 1), (-1, -1), colors.white),  # Body background white
                ("GRID", (0, 0), (-1, -1), 1, colors.black),  # Grid lines
            ]
        )

        table.setStyle(style)

        elements.append(table)
        # elements.append(Spacer(1, 12))
        elements.append(PageBreak())

    # elements.append(Spacer(1, 12))

    if os.path.exists(os.path.join(output, "CNV.npy")):
        elements.append(Paragraph("Copy Number Variation", underline_style))
        elements.append(
            Paragraph("These plots were generated with cnv_from_bam.", smaller_style)
        )

        cnv_summary = create_CNV_plot(CNVresult, cnv_dict)
        img = Image(cnv_summary, width=6 * inch, height=1.5 * inch)
        elements.append(img)
        elements.append(Spacer(1, 6))

        cnv_plots = create_CNV_plot_per_chromosome(CNVresult, cnv_dict)

        for contig, img_buf in cnv_plots:
            # elements.append(Paragraph(f"Copy Number Variation - {contig}", smaller_style))
            img = Image(img_buf, width=6 * inch, height=1.5 * inch)
            elements.append(img)
            # elements.append(Spacer(1, 6))

        if XYestimate != "Unknown":
            # if XYestimate == "XY":
            #    ui.icon("man").classes("text-4xl")
            # else:
            #    ui.icon("woman").classes("text-4xl")
            elements.append(
                Paragraph(f"Estimated Genetic Sex: {XYestimate}", smaller_style)
            )
        elements.append(
            Paragraph(f"Current Bin Width: {cnv_dict['bin_width']}", smaller_style)
        )
        elements.append(
            Paragraph(
                f"Current Variance: {round(cnv_dict['variance'], 3)}", smaller_style
            )
        )

        elements.append(Spacer(1, 12))
        elements.append(PageBreak())

    elements.append(Paragraph("Target Coverage", underline_style))
    if os.path.isfile(os.path.join(output, "coverage_main.csv")):
        elements.append(Paragraph("This plot was generated by ROBIN.", smaller_style))
        img_buf = coverage_plot(cov_df_main)
        width, height = A4
        img = Image(img_buf, width=width * 0.95, height=width, kind="proportional")
        elements.append(img)
        elements.append(Spacer(1, 12))

        elements.append(
            Paragraph(
                "Coverage over individual targets on each chromosome. Outliers are annotated by gene name.",
                smaller_style,
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
                smaller_style,
            )
        )
        # Get the outliers
        outliers = get_target_outliers(target_coverage_df)

        # Round the coverage values to 5 decimal places
        outliers["coverage"] = outliers["coverage"].apply(lambda x: round(x, 1))

        # Sort the dataframe by coverage in descending order
        outliers = outliers.sort_values(by="coverage", ascending=False)

        data = [outliers.columns.to_list()] + outliers.values.tolist()
        table = Table(data)
        # Add style to the table
        style = TableStyle(
            [
                (
                    "BACKGROUND",
                    (0, 0),
                    (-1, -1),
                    colors.white,
                ),  # Set background to white
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),  # Header text color
                ("TEXTCOLOR", (0, 1), (-1, -1), colors.black),  # Body text color
                ("ALIGN", (0, 0), (-1, -1), "CENTER"),  # Center align all cells
                ("FONTNAME", (0, 0), (-1, 0), "FiraSans-Bold"),  # Header font
                ("FONTNAME", (0, 1), (-1, -1), "FiraSans"),  # Body font
                ("FONTSIZE", (0, 0), (-1, 0), 6),  # Header font size
                ("FONTSIZE", (0, 1), (-1, -1), 5),  # Body font size
                ("BOTTOMPADDING", (0, 0), (-1, 0), 1),  # Header padding
                ("BACKGROUND", (0, 1), (-1, -1), colors.white),  # Body background white
                ("GRID", (0, 0), (-1, -1), 1, colors.black),  # Grid lines
            ]
        )

        table.setStyle(style)

        elements.append(table)
        elements.append(Spacer(1, 12))

        if os.path.exists(f"{output}/clair3/snpsift_output.vcf.csv"):
            elements.append(Paragraph("Pathogenic Variants", smaller_style))
            vcf = pd.read_csv(f"{output}/clair3/snpsift_output.vcf.csv")
            pathogenic_vcf = vcf[
                vcf["CLNSIG"].notna() & vcf["CLNSIG"].str.contains("pathogenic")
            ].loc[
                :,
                [
                    "CHROM",
                    "POS",
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                    "Allele",
                    "Annotation",
                    "Annotation_Impact",
                    "Gene_Name",
                    "Gene_ID",
                ],
            ]
            if len(pathogenic_vcf) > 0:
                data = [
                    pathogenic_vcf.columns.to_list()
                ] + pathogenic_vcf.values.tolist()
                table = Table(data)
                # Add style to the table
                style = TableStyle(
                    [
                        (
                            "BACKGROUND",
                            (0, 0),
                            (-1, -1),
                            colors.white,
                        ),  # Set background to white
                        (
                            "TEXTCOLOR",
                            (0, 0),
                            (-1, 0),
                            colors.black,
                        ),  # Header text color
                        (
                            "TEXTCOLOR",
                            (0, 1),
                            (-1, -1),
                            colors.black,
                        ),  # Body text color
                        ("ALIGN", (0, 0), (-1, -1), "CENTER"),  # Center align all cells
                        ("FONTNAME", (0, 0), (-1, 0), "FiraSans-Bold"),  # Header font
                        ("FONTNAME", (0, 1), (-1, -1), "FiraSans"),  # Body font
                        ("FONTSIZE", (0, 0), (-1, 0), 6),  # Header font size
                        ("FONTSIZE", (0, 1), (-1, -1), 5),  # Body font size
                        ("BOTTOMPADDING", (0, 0), (-1, 0), 1),  # Header padding
                        (
                            "BACKGROUND",
                            (0, 1),
                            (-1, -1),
                            colors.white,
                        ),  # Body background white
                        ("GRID", (0, 0), (-1, -1), 1, colors.black),  # Grid lines
                    ]
                )

                table.setStyle(style)
                elements.append(table)
        if os.path.exists(f"{output}/clair3/snpsift_indel_output.vcf.csv"):
            elements.append(Paragraph("Pathogenic Variants (in/dels)", smaller_style))
            vcf = pd.read_csv(f"{output}/clair3/snpsift_output.vcf.csv")
            pathogenic_vcf = vcf[
                vcf["CLNSIG"].notna() & vcf["CLNSIG"].str.contains("pathogenic")
            ].loc[
                :,
                [
                    "CHROM",
                    "POS",
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                    "Allele",
                    "Annotation",
                    "Annotation_Impact",
                    "Gene_Name",
                    "Gene_ID",
                ],
            ]
            if len(pathogenic_vcf) > 0:
                data = [
                    pathogenic_vcf.columns.to_list()
                ] + pathogenic_vcf.values.tolist()
                table = Table(data)
                # Add style to the table
                style = TableStyle(
                    [
                        (
                            "BACKGROUND",
                            (0, 0),
                            (-1, -1),
                            colors.white,
                        ),  # Set background to white
                        (
                            "TEXTCOLOR",
                            (0, 0),
                            (-1, 0),
                            colors.black,
                        ),  # Header text color
                        (
                            "TEXTCOLOR",
                            (0, 1),
                            (-1, -1),
                            colors.black,
                        ),  # Body text color
                        ("ALIGN", (0, 0), (-1, -1), "CENTER"),  # Center align all cells
                        ("FONTNAME", (0, 0), (-1, 0), "FiraSans-Bold"),  # Header font
                        ("FONTNAME", (0, 1), (-1, -1), "FiraSans"),  # Body font
                        ("FONTSIZE", (0, 0), (-1, 0), 6),  # Header font size
                        ("FONTSIZE", (0, 1), (-1, -1), 5),  # Body font size
                        ("BOTTOMPADDING", (0, 0), (-1, 0), 1),  # Header padding
                        (
                            "BACKGROUND",
                            (0, 1),
                            (-1, -1),
                            colors.white,
                        ),  # Body background white
                        ("GRID", (0, 0), (-1, -1), 1, colors.black),  # Grid lines
                    ]
                )

                table.setStyle(style)
                elements.append(table)
        elements.append(PageBreak())

    elements.append(Spacer(1, 12))
    """
    elements.append(Paragraph("Structural Variants", underline_style))

    if os.path.exists(os.path.join(output, "fusion_candidates.csv")):
        elements.append(
            Paragraph(
                "These are high confidence fusions. These fusions involve an interactions between genes in the candidate panel.",
                smaller_style,
            )
        )
        elements.append(
            Paragraph(
                f"There are {result[goodpairs].sort_values(by=7)['tag'].nunique()} high confidence fusions observed.",
                smaller_style,
            )
        )
        for gene_pair in result[goodpairs].sort_values(by=7)["tag"].unique():
            read_count = result[result["tag"].isin([gene_pair])][7].nunique()
            elements.append(
                Paragraph(f"{gene_pair}: {read_count} reads", smaller_style)
            )

        elements.append(Spacer(1, 12))
        for gene_pair in result[goodpairs].sort_values(by=7)["tag"].unique():
            elements.append(Paragraph(f"Fusion Candidate: {gene_pair}", smaller_style))
            for gene in result[goodpairs & result[goodpairs]["tag"].eq(gene_pair)][
                3
            ].unique():
                title = gene
                reads = result[goodpairs].sort_values(by=7)[
                    result[goodpairs].sort_values(by=7)[3].eq(gene)
                ]
                buf, buf2 = fusion_plot(title, reads)
                width, height = A4
                img = Image(buf, width=width * 0.6, height=width, kind="proportional")
                elements.append(img)
                width, height = A4
                img = Image(buf2, width=width * 0.6, height=width, kind="proportional")
                elements.append(img)
                elements.append(Spacer(1, 12))

    if os.path.exists(os.path.join(output, "fusion_candidates_master.csv")):
        elements.append(
            Paragraph(
                "These are low confidence fusions. These fusions involve an interactions between one target gene and any other gene in the genome.",
                smaller_style,
            )
        )
        elements.append(
            Paragraph(
                f"There are {result_all[goodpairs_all].sort_values(by=7)['tag'].nunique()} low confidence fusions observed.",
                smaller_style,
            )
        )
        for gene_pair in result_all[goodpairs_all].sort_values(by=7)["tag"].unique():
            read_count = result_all[result_all["tag"].isin([gene_pair])][7].nunique()
            elements.append(
                Paragraph(f"{gene_pair}: {read_count} reads", smaller_style)
            )

        elements.append(Spacer(1, 12))
        if len(result_all[goodpairs_all].sort_values(by=7)["tag"].unique()) > 0:
            for gene_pair in (
                result_all[goodpairs_all].sort_values(by=7)["tag"].unique()
            ):
                elements.append(
                    Paragraph(f"Fusion Candidate: {gene_pair}", smaller_style)
                )
                print (gene_pair)
                print(result_all)
                for gene in result_all[
                    goodpairs_all & result_all[goodpairs_all]["tag"].eq(gene_pair)
                ][3].unique():
                    title = gene
                    reads = result_all[goodpairs_all].sort_values(by=7)[
                        result_all[goodpairs_all].sort_values(by=7)[3].eq(gene)
                    ]
                    buf, buf2 = fusion_plot(title, reads)
                    width, height = A4
                    img = Image(
                        buf, width=width * 0.6, height=width, kind="proportional"
                    )
                    elements.append(img)
                    width, height = A4
                    img = Image(
                        buf2, width=width * 0.6, height=width, kind="proportional"
                    )
                    elements.append(img)
                    elements.append(Spacer(1, 12))
    """
    elements.append(Spacer(1, 12))

    last_seen = 0
    if not last_seen:
        for file in natsort.natsorted(os.listdir(output)):
            if file.endswith("_mgmt.csv"):
                count = int(file.split("_")[0])
                if count > last_seen:
                    results = pd.read_csv(os.path.join(output, file))
                    plot_out = os.path.join(output, file.replace(".csv", ".png"))
                    last_seen = count
        elements.append(PageBreak())

    if last_seen > 0:
        elements.append(Paragraph("MGMT Promoter Methylation", underline_style))
        image = Image(plot_out, 6 * inch, 4 * inch)  # Adjust the size as needed
        elements.append(image)

        data_list = [results.columns.values.tolist()] + results.values.tolist()

        table = Table(data_list)

        style = TableStyle(
            [
                (
                    "BACKGROUND",
                    (0, 0),
                    (-1, -1),
                    colors.white,
                ),  # Set background to white
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),  # Header text color
                ("TEXTCOLOR", (0, 1), (-1, -1), colors.black),  # Body text color
                ("ALIGN", (0, 0), (-1, -1), "CENTER"),  # Center align all cells
                ("FONTNAME", (0, 0), (-1, 0), "FiraSans-Bold"),  # Header font
                ("FONTNAME", (0, 1), (-1, -1), "FiraSans"),  # Body font
                ("FONTSIZE", (0, 0), (-1, 0), 6),  # Header font size
                ("FONTSIZE", (0, 1), (-1, -1), 5),  # Body font size
                ("BOTTOMPADDING", (0, 0), (-1, 0), 1),  # Header padding
                ("BACKGROUND", (0, 1), (-1, -1), colors.white),  # Body background white
                ("GRID", (0, 0), (-1, -1), 1, colors.black),  # Grid lines
            ]
        )
        table.setStyle(style)

        elements.append(table)

        elements.append(Spacer(1, 12))
        elements.append(PageBreak())

    # Add sections

    doc.multiBuild(elements, canvasmaker=header_footer_canvas_factory(sample_id))
    print(f"PDF created: {filename}")
    return filename


if __name__ == "__main__":
    # Create a sample plot
    # Register the custom fonts
    pdfmetrics.registerFont(
        TTFont(
            "FiraSans",
            os.path.join(
                os.path.dirname(os.path.abspath(fonts.__file__)),
                "fira-sans-v16-latin-regular.ttf",
            ),
        )
    )
    pdfmetrics.registerFont(
        TTFont(
            "FiraSans-Bold",
            os.path.join(
                os.path.dirname(os.path.abspath(fonts.__file__)),
                "fira-sans-v16-latin-700.ttf",
            ),
        )
    )

    # Update styles to use the custom font
    styles = getSampleStyleSheet()
    for style_name in styles.byName:
        styles[style_name].fontName = "FiraSans"

    # Define a smaller style for the header date
    smaller_style = ParagraphStyle(name="Smaller", parent=styles["Normal"], fontSize=8)

    # Define a bold style for the first header line
    bold_style = ParagraphStyle(
        name="Bold", parent=styles["Normal"], fontName="FiraSans-Bold"
    )

    # Define an underlined style for section headings
    underline_style = ParagraphStyle(
        name="Underline", parent=styles["Heading1"], underline=True
    )
    # Generate the PDF
    create_pdf(
        "sample_report.pdf",
        "/Users/mattloose/GIT/niceGUI/cnsmeth/fusion_output/ds1305_Intraop0047_b",
    )
