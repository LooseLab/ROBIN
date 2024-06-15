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
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import natsort
import seaborn as sns
import pandas as pd
import numpy as np
import io
import os
from datetime import datetime
import base64
import pickle

from robin.subpages.CNV_object import Result


def convert_to_space_separated_string(array):
    import ast

    # Convert array to list and extract the string
    string_repr = array.tolist()[0]

    # Evaluate the string to convert it to an actual list
    list_repr = ast.literal_eval(string_repr)

    # Join the elements of the list into a space-separated string
    return " ".join(list_repr)


async def generate_image(myplot):
    # Run the getDataURL method
    data_url = await myplot.run_chart_method(
        "getDataURL", {"type": "png", "pixelRatio": 2, "backgroundColor": "#fff"}
    )
    print(data_url)
    # Create an image element with the data URL
    # ui.image(data_url, width='600px', height='400px')
    # Decode and save the image to a file
    header, encoded = data_url.split(",", 1)
    data = base64.b64decode(encoded)
    with open("chart.png", "wb") as f:
        f.write(data)


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
        logo_path = "src/robin/images/MethBrain_small.png"  # Replace with the path to your logo
        max_logo_size = 40  # Maximum width and height in pixels
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
        page = "SampleID: %s - Page %s of %s" % ("SAMPLE", self._pageNumber, page_count)
        x = 190
        self.saveState()
        self.setStrokeColorRGB(0, 0, 0)
        self.setLineWidth(0.5)
        self.line(66, 78, A4[0] - 66, 78)
        self.setFont("FiraSans", 10)
        self.drawString(A4[0] - x, 65, page)
        self.restoreState()


def header_footer_canvas_factory(sample_id):
    def create_canvas(*args, **kwargs):
        return HeaderFooterCanvas(sample_id, *args, **kwargs)

    return create_canvas


def create_plot():
    plt.figure(figsize=(6, 4))
    x = np.linspace(0, 10, 100)
    y = np.sin(x)
    plt.plot(x, y, label="Sine wave")
    plt.title("Sample Plot")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.legend()
    plt.savefig("sample_plot.png")
    plt.close()


def target_distribution_plot(df):
    df["chrom"] = pd.Categorical(
        df["chrom"], categories=natsort.natsorted(df["chrom"].unique()), ordered=True
    )
    df = df.sort_values("chrom")

    # Generate the plot
    plt.figure(figsize=(14, 8))
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
    plt.savefig(buf, format="jpg", dpi=600)
    buf.seek(0)
    return buf


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


def coverage_plot(df):
    df = df[df["#rname"] != "chrM"].copy()

    # Sort chromosomes naturally
    df["#rname"] = pd.Categorical(
        df["#rname"], categories=natsort.natsorted(df["#rname"].unique()), ordered=True
    )
    df = df.sort_values("#rname")

    # Create subplots
    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1])

    # Plot number of reads per chromosome
    ax0 = plt.subplot(gs[0])
    sns.barplot(x="#rname", y="numreads", data=df, ax=ax0)
    ax0.set_title("Number of Reads per Chromosome")
    ax0.set_xlabel("")
    ax0.set_ylabel("Number of Reads")
    ax0.tick_params(axis="x", rotation=90)

    # Plot number of bases per chromosome
    ax1 = plt.subplot(gs[1])
    sns.barplot(x="#rname", y="covbases", data=df, ax=ax1)
    ax1.set_title("Number of Bases per Chromosome")
    ax1.set_xlabel("")
    ax1.set_ylabel("Number of Bases")
    ax1.tick_params(axis="x", rotation=90)

    # Plot mean depth per chromosome
    ax2 = plt.subplot(gs[2])
    sns.barplot(x="#rname", y="meandepth", data=df, ax=ax2)
    ax2.set_title("Mean Depth per Chromosome")
    ax2.set_xlabel("Chromosome")
    ax2.set_ylabel("Mean Depth")
    ax2.tick_params(axis="x", rotation=90)

    plt.tight_layout()
    # Save the plot as a JPG file
    buf = io.BytesIO()
    plt.savefig(buf, format="jpg", dpi=600)
    buf.seek(0)
    return buf


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
    plt.figure(figsize=(14, 6))
    sns.lineplot(data=df_filtered, x="timestamp", y="Value", hue="Condition")
    plt.title(f"{title} Classifications over Time")
    plt.xlabel("Timestamp")
    plt.ylabel("Value")
    plt.legend(title="Condition", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.xticks(rotation=45)

    # Format the x-axis with custom date format
    ax = plt.gca()
    ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter("%d %b %Y %H:%M"))

    plt.legend(title="Condition", loc="center right")

    plt.tight_layout()

    # Save the plot as a JPG file
    buf = io.BytesIO()
    plt.savefig(buf, format="jpg", dpi=600)
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
        s=4,
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
    plt.savefig(buf, format="jpg", dpi=600)
    buf.seek(0)
    return buf


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
    doc = SimpleDocTemplate(filename, pagesize=A4)
    elements = []

    if os.path.exists(os.path.join(output, "master.csv")):
        masterdf = pd.read_csv(
            os.path.join(output, "master.csv"), index_col=0, header=None
        )
    else:
        masterdf = None

    if masterdf is not None and isinstance(masterdf, pd.DataFrame):
        sample_id = convert_to_space_separated_string(
            masterdf.loc[(masterdf.index == "sample_ids")][1].values
        )

    else:
        sample_id = None

    import pprint

    pprint.pprint(masterdf)

    elements.append(Paragraph("Classification Summary", styles["Heading1"]))

    elements.append(Spacer(1, 12))

    elements.append(
        Paragraph("This sample has the following classifications:", styles["BodyText"])
    )
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
                styles["BodyText"],
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
                styles["BodyText"],
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
                styles["BodyText"],
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
            Paragraph(f"Estimated Genetic Sex: {XYestimate}", styles["BodyText"])
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
    else:
        elements.append(Paragraph("No Coverage Data Available", styles["BodyText"]))
    # Add summary section

    elements.append(
        Paragraph(
            "Fusion Candidates - using panel rCNS2 <br/> 0 high confidence fusions observed.<br/>0 low confidence fusions observed.",
            styles["BodyText"],
        )
    )

    if masterdf is not None and isinstance(masterdf, pd.DataFrame):
        elements.append(Paragraph("Run Data Summary", styles["Heading1"]))
        start_time = "Placeholder"
        elements.append(
            Paragraph(
                f"Sample ID: {sample_id}<br/>"
                f"Run Start: {convert_to_space_separated_string(masterdf.loc[(masterdf.index == 'run_time')][1].values)}<br/>"
                f"Run Folder: {masterdf.loc[(masterdf.index == 'watchfolder')][1].values}<br/>"
                f"Output Folder: {masterdf.loc[(masterdf.index == 'output')][1].values}<br/>"
                f"Target Panel: {masterdf.loc[(masterdf.index == 'target_panel')][1].values}<br/>"
                f"Reference: {masterdf.loc[(masterdf.index == 'reference')][1].values}<br/>"
                f"Sequencing Device: {convert_to_space_separated_string(masterdf.loc[(masterdf.index == 'devices')][1].values)}<br/>"
                f"Flowcell ID: {convert_to_space_separated_string(masterdf.loc[(masterdf.index == 'flowcell_ids')][1].values)}<br/>"
                f"Basecalling Model: {convert_to_space_separated_string(masterdf.loc[(masterdf.index == 'basecall_models')][1].values)}<br/>"
                f"Files Seen: {masterdf.loc[(masterdf.index == 'file_counters')][1].values}<br/>",
                styles["BodyText"],
            )
        )

    elements.append(Spacer(1, 12))

    elements.append(PageBreak())

    elements.append(Paragraph("Methylation Classifications", styles["Heading2"]))

    if os.path.exists(os.path.join(output, "sturgeon_scores.csv")):
        elements.append(Paragraph("Sturgeon Classification", styles["Heading3"]))
        elements.append(
            Paragraph("This plot was generated by sturgeon.", smaller_style)
        )
        img_buf = classification_plot(sturgeon_df_store, "Sturgeon", 0.05)
        width, height = A4
        img = Image(img_buf, width=width, height=width / 3, kind="proportional")
        elements.append(img)
        df = Sturgeonlastrow_plot.reset_index()
        df.columns = ["Classification", "Score"]
        df["Score"] = df["Score"].apply(lambda x: round(x, 5))
        data = [df.columns.to_list()] + df.values.tolist()
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
                ("FONTSIZE", (0, 0), (-1, 0), 8),  # Header font size
                ("FONTSIZE", (0, 1), (-1, -1), 8),  # Body font size
                ("BOTTOMPADDING", (0, 0), (-1, 0), 1),  # Header padding
                ("BACKGROUND", (0, 1), (-1, -1), colors.white),  # Body background white
                ("GRID", (0, 0), (-1, -1), 1, colors.black),  # Grid lines
            ]
        )

        table.setStyle(style)

        elements.append(table)
        elements.append(Spacer(1, 12))

    if os.path.exists(os.path.join(output, "nanoDX_scores.csv")):
        elements.append(Paragraph("NanoDX Classification", styles["Heading3"]))
        elements.append(Paragraph("This plot was generated by NanoDX.", smaller_style))
        img_buf = classification_plot(nanodx_df_store, "NanoDX", 0.05)
        width, height = A4
        img = Image(img_buf, width=width, height=width / 3, kind="proportional")
        elements.append(img)
        df = NanoDXlastrow_plot.reset_index()
        df.columns = ["Classification", "Score"]
        df["Score"] = df["Score"].apply(lambda x: round(x, 5))
        data = [df.columns.to_list()] + df.values.tolist()
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
                ("FONTSIZE", (0, 0), (-1, 0), 8),  # Header font size
                ("FONTSIZE", (0, 1), (-1, -1), 8),  # Body font size
                ("BOTTOMPADDING", (0, 0), (-1, 0), 1),  # Header padding
                ("BACKGROUND", (0, 1), (-1, -1), colors.white),  # Body background white
                ("GRID", (0, 0), (-1, -1), 1, colors.black),  # Grid lines
            ]
        )

        table.setStyle(style)

        elements.append(table)
        elements.append(Spacer(1, 12))

    if os.path.exists(os.path.join(output, "random_forest_scores.csv")):
        elements.append(Paragraph("Forest Classification", styles["Heading3"]))
        elements.append(
            Paragraph("This plot was generated by Random Forest.", smaller_style)
        )
        img_buf = classification_plot(rcns2_df_store, "Forest", 0.05)
        width, height = A4
        img = Image(img_buf, width=width, height=width / 3, kind="proportional")
        elements.append(img)
        df = Forestlastrow_plot.reset_index()
        df.columns = ["Classification", "Score"]
        df["Score"] = df["Score"].apply(lambda x: round(x, 5))
        data = [df.columns.to_list()] + df.values.tolist()
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
                ("FONTSIZE", (0, 0), (-1, 0), 8),  # Header font size
                ("FONTSIZE", (0, 1), (-1, -1), 8),  # Body font size
                ("BOTTOMPADDING", (0, 0), (-1, 0), 1),  # Header padding
                ("BACKGROUND", (0, 1), (-1, -1), colors.white),  # Body background white
                ("GRID", (0, 0), (-1, -1), 1, colors.black),  # Grid lines
            ]
        )

        table.setStyle(style)

        elements.append(table)
        elements.append(Spacer(1, 12))

    elements.append(Spacer(1, 12))

    if os.path.exists(os.path.join(output, "CNV.npy")):
        elements.append(Paragraph("Copy Number Variation", underline_style))
        elements.append(
            Paragraph("This plot was generated with cnv_from_bam.", smaller_style)
        )
        img_buf = create_CNV_plot(CNVresult, cnv_dict)
        width, height = A4
        img = Image(img_buf, width=width, height=width / 3, kind="proportional")
        elements.append(img)
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

    elements.append(Paragraph("Target Coverage", underline_style))
    if os.path.isfile(os.path.join(output, "coverage_main.csv")):
        elements.append(Paragraph("This plot was generated by ROBIN.", smaller_style))
        img_buf = coverage_plot(cov_df_main)
        width, height = A4
        img = Image(img_buf, width=width * 0.9, height=width, kind="proportional")
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
        outliers = get_target_outliers(target_coverage_df)
        outliers["coverage"] = outliers["coverage"].apply(lambda x: round(x, 5))
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
                ("FONTSIZE", (0, 0), (-1, 0), 8),  # Header font size
                ("FONTSIZE", (0, 1), (-1, -1), 8),  # Body font size
                ("BOTTOMPADDING", (0, 0), (-1, 0), 1),  # Header padding
                ("BACKGROUND", (0, 1), (-1, -1), colors.white),  # Body background white
                ("GRID", (0, 0), (-1, -1), 1, colors.black),  # Grid lines
            ]
        )

        table.setStyle(style)

        elements.append(table)

    elements.append(Spacer(1, 12))

    elements.append(Paragraph("Structural Variants", underline_style))

    elements.append(Spacer(1, 12))

    elements.append(Paragraph("MGMT Promoter Methylation", underline_style))

    elements.append(Spacer(1, 12))

    elements.append(Spacer(1, 12))

    # Add sections

    doc.multiBuild(elements, canvasmaker=header_footer_canvas_factory(sample_id))

    return filename


if __name__ == "__main__":
    # Create a sample plot
    # Register the custom fonts
    pdfmetrics.registerFont(
        TTFont("FiraSans", "src/robin/fonts/fira-sans-v16-latin-regular.ttf")
    )
    pdfmetrics.registerFont(
        TTFont("FiraSans-Bold", "src/robin/fonts/fira-sans-v16-latin-700.ttf")
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
    # Generate the PDF
    create_pdf("sample_report.pdf", "run2")
