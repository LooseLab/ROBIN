"""
Fusion Gene Identification Module

This module provides functionality for identifying gene fusion candidates from BAM files.
The primary components include running external tools to extract and process gene fusion
candidates, visualizing the results using a GUI, and managing the process asynchronously.

Dependencies:
    - pandas
    - os
    - sys
    - subprocess
    - gff3_parser
    - tempfile
    - random
    - nicegui
    - dna_features_viewer
    - matplotlib
    - click

Modules:
    - subpages.base_analysis: BaseAnalysis class from robin.subpages.base_analysis
    - theme: robin theme module
    - resources: robin resources module

Environment Variables:
    - CI: Set to "1"

Constants:
    - STRAND: Dictionary for mapping strand symbols to integers

Functions:
    - run_command(command: str): Executes a shell command and handles exceptions.
    - fusion_work(threads, bamfile, gene_bed, all_gene_bed, tempreadfile, tempbamfile, tempmappings, tempallmappings): Identifies fusion candidates from BAM files.
    - fusion_work_old(threads, bamfile, gene_bed, all_gene_bed, tempreadfile, tempbamfile, tempmappings, tempallmappings): Legacy function for identifying fusion candidates.

Classes:
    - FusionObject(BaseAnalysis): Manages the gene fusion identification process, including setting up the GUI and handling BAM file processing.

Command-line Interface:
    - main(port, threads, watchfolder, output, browse): CLI entry point for running the app, using Click for argument parsing.

Usage:
    The module can be run as a script to start the GUI for gene fusion identification, specifying
    options like the port, number of threads, watch folder, and output directory.
"""

import os
import sys
import subprocess
import gff3_parser
import tempfile
import random
import pysam
import logging
import networkx as nx
import pandas as pd
import click
from typing import Optional, Tuple
from nicegui import ui, run, app
from robin import theme, resources
from dna_features_viewer import GraphicFeature, GraphicRecord
from pathlib import Path
import matplotlib
from matplotlib.ticker import FuncFormatter
from robin.subpages.base_analysis import BaseAnalysis
from robin.utilities.decompress import decompress_gzip_file
from collections import Counter

matplotlib.use("agg")
from matplotlib import pyplot as plt

# Use the main logger configured in the main application
logger = logging.getLogger(__name__)

os.environ["CI"] = "1"
STRAND = {"+": 1, "-": -1}

decompress_gzip_file(
    os.path.join(
        os.path.dirname(os.path.abspath(resources.__file__)),
        "gencode.v45.basic.annotation.gff3.gz",
    )
)


def _get_reads(reads: pd.DataFrame) -> pd.DataFrame:
    """
    Get reads for a specific gene.

    Args:
        reads (pd.DataFrame): DataFrame with reads

    Returns:
        pd.DataFrame: DataFrame with reads
    """
    df = reads
    df.columns = [
        "chromosome",
        "start",
        "end",
        "gene",
        "chromosome2",
        "start2",
        "end2",
        "id",
        "quality",
        "strand",
        "read_start",
        "read_end",
        "secondary",
        "supplementary",
        "span",
        "tag",
        "color",
    ]
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)

    # Sort the DataFrame by chromosome, start, and end positions
    df = df.sort_values(by=["chromosome", "start", "end"])

    df = df.drop_duplicates(subset=["start2", "end2", "id"])

    # Group by chromosome and collapse ranges within each group
    result = (
        df.groupby("chromosome")
        .apply(lambda x: collapse_ranges(x, 10000))
        .reset_index(drop=True)
    )
    # return reads[reads[3].eq(gene)]
    return result


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
    colors = result.groupby(7).apply(lambda x: _generate_random_color())
    result = result.map(lambda x: x.strip() if isinstance(x, str) else x)
    result["Color"] = result[7].map(colors.get)
    goodpairs = result.groupby("tag")[7].transform("nunique") > 2
    # gene_pairs = result[goodpairs].sort_values(by=7)["tag"].unique().tolist()
    return result, goodpairs


def _generate_random_color() -> str:
    """
    Generates a random color for use in plotting.

    Returns:
        str: A random hex color code.
    """
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))


def get_gene_network(gene_pairs):
    G = nx.Graph()
    for pair in gene_pairs:
        G.add_edge(pair[0], pair[1])
    connected_components = list(nx.connected_components(G))
    return [list(component) for component in connected_components]


def extract_bam_info(bam_file):
    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Initialize lists to store information
    read_ids = []
    reference_ids = []
    reference_starts = []
    reference_ends = []
    read_starts = []
    read_ends = []
    strands = []
    mapping_qualities = []
    is_primary = []
    is_secondary = []
    is_supplementary = []

    # Iterate over each read in the BAM file
    for read in bam:
        # Append the read information to the lists
        read_ids.append(read.query_name)
        reference_ids.append(bam.get_reference_name(read.reference_id))
        reference_starts.append(read.reference_start)
        reference_ends.append(read.reference_end)
        read_starts.append(read.query_alignment_start)
        read_ends.append(read.query_alignment_end)
        strands.append("-" if read.is_reverse else "+")
        mapping_qualities.append(read.mapping_quality)
        # is_primary.append(read.is_primary)
        is_secondary.append(read.is_secondary)
        is_supplementary.append(read.is_supplementary)

    # Close the BAM file
    bam.close()

    # Create a DataFrame
    df = pd.DataFrame(
        {
            "read_id": read_ids,
            "reference_id": reference_ids,
            "reference_start": reference_starts,
            "reference_end": reference_ends,
            "read_start": read_starts,
            "read_end": read_ends,
            "strand": strands,
            "mapping_quality": mapping_qualities,
            # 'is_primary': is_primary,
            "is_secondary": is_secondary,
            "is_supplementary": is_supplementary,
        }
    )

    return df


# Function to collapse ranges within a fixed distance
def collapse_ranges(df, max_distance):
    collapsed = []
    current_range = None

    for _, row in df.iterrows():
        chromosome, start, end = row["chromosome"], row["start"], row["end"]

        if current_range is None:
            current_range = row
        else:
            if start <= current_range["end"] + max_distance:
                current_range["end"] = max(current_range["end"], end)
            else:
                collapsed.append(current_range)
                current_range = row

    if current_range is not None:
        collapsed.append(current_range)

    return pd.DataFrame(collapsed)


def run_command(command: str) -> None:
    """
    Run a shell command and handle exceptions.

    Args:
        command (str): The shell command to run.
    """
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        # logging.error(f"Command failed: {command}")
        raise e


def fusion_work(
    threads: int,
    bamfile: str,
    gene_bed: str,
    all_gene_bed: str,
    tempreadfile: str,
    tempbamfile: str,
    tempmappings: str,
    tempallmappings: str,
) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """
    Identify fusion candidates from BAM file based on gene BED files.

    Args:
        threads (int): Number of threads to use for parallel processing.
        bamfile (str): Path to the BAM file.
        gene_bed (str): Path to the gene BED file.
        all_gene_bed (str): Path to the BED file with all genes.
        tempreadfile (str): Path to the temporary file to store read names.
        tempbamfile (str): Path to the temporary BAM file.
        tempmappings (str): Path to the temporary file for gene mappings.
        tempallmappings (str): Path to the temporary file for all gene mappings.

    Returns:
        Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]: DataFrames with fusion candidates and all fusion candidates.
    """
    fusion_candidates: Optional[pd.DataFrame] = None
    fusion_candidates_all: Optional[pd.DataFrame] = None

    try:
        # Extract reads mapped with supplementary alignments
        run_command(
            f"samtools view -@{threads} -L {gene_bed} -d SA {bamfile} | cut -f1 > {tempreadfile}"
        )

        if os.path.getsize(tempreadfile) > 0:
            # Extract reads to a temporary BAM file
            run_command(
                f"samtools view -@{threads}  --write-index -N {tempreadfile} -o {tempbamfile} {bamfile}"
            )
            if os.path.getsize(tempbamfile) > 0:
                read_data = extract_bam_info(tempbamfile)
                # Intersect the gene BED file with the temporary BAM file
                run_command(
                    f"bedtools intersect -a {gene_bed} -b {tempbamfile} -wa -wb > {tempmappings}"
                )
                run_command(
                    f"bedtools intersect -a {all_gene_bed} -b {tempbamfile} -wa -wb > {tempallmappings}"
                )

                # Process gene mappings if available

                samfile = pysam.AlignmentFile(f"{tempbamfile}", "r")

                if os.path.getsize(tempmappings) > 0:
                    fusion_candidates = pd.read_csv(tempmappings, sep="\t", header=None)
                    # Merge on common columns
                    # Find columns in df1 that are not in df2
                    fusion_candidates.columns = [
                        "col1",
                        "col2",
                        "col3",
                        "col4",
                        "reference_id",
                        "reference_start",
                        "reference_end",
                        "read_id",
                        "mapping_quality",
                        "strand",
                    ]
                    additional_columns = [
                        col
                        for col in read_data.columns
                        if col not in fusion_candidates.columns
                    ]
                    fusion_candidates = pd.merge(
                        fusion_candidates,
                        read_data[
                            ["read_id", "reference_start", "reference_end"]
                            + additional_columns
                        ],
                        on=["read_id", "reference_start", "reference_end"],
                        how="left",
                    )
                    fusion_candidates = fusion_candidates[
                        fusion_candidates["mapping_quality"] > 40
                    ]
                    fusion_candidates["diff"] = (
                        fusion_candidates["reference_end"]
                        - fusion_candidates["reference_start"]
                    )
                    fusion_candidates = fusion_candidates[
                        fusion_candidates["diff"] > 100
                    ]

                # Process all gene mappings if available
                if os.path.getsize(tempallmappings) > 0:
                    fusion_candidates_all = pd.read_csv(
                        tempallmappings, sep="\t", header=None
                    )
                    fusion_candidates_all.columns = [
                        "col1",
                        "col2",
                        "col3",
                        "col4",
                        "reference_id",
                        "reference_start",
                        "reference_end",
                        "read_id",
                        "mapping_quality",
                        "strand",
                    ]

                    additional_columns = [
                        col
                        for col in read_data.columns
                        if col not in fusion_candidates_all.columns
                    ]
                    fusion_candidates_all = pd.merge(
                        fusion_candidates_all,
                        read_data[
                            ["read_id", "reference_start", "reference_end"]
                            + additional_columns
                        ],
                        on=["read_id", "reference_start", "reference_end"],
                        how="left",
                    )
                    fusion_candidates_all = fusion_candidates_all[
                        fusion_candidates_all["mapping_quality"] > 40
                    ]
                    fusion_candidates_all["diff"] = (
                        fusion_candidates_all["reference_end"]
                        - fusion_candidates_all["reference_start"]
                    )
                    fusion_candidates_all = fusion_candidates_all[
                        fusion_candidates_all["diff"] > 100
                    ]

                if os.path.isfile(f"{tempbamfile}.csi"):
                    logger.debug(f"removing {tempbamfile}.csi")
                    os.remove(f"{tempbamfile}.csi")

                samfile.close()

    except Exception:
        # logging.error(f"Error during fusion work: {e}")
        raise

    return fusion_candidates, fusion_candidates_all


class FusionObject(BaseAnalysis):
    """
    FusionObject handles the gene fusion identification process, including setting up the GUI
    and processing BAM files asynchronously.

    Attributes:
        target_panel (str): Name of the target panel.
        fusion_candidates (pd.DataFrame): DataFrame with fusion candidates.
        fusion_candidates_all (pd.DataFrame): DataFrame with all fusion candidates.
        fstable_all (ui.Table): UI table for displaying all fusion candidates.
        fstable (ui.Table): UI table for displaying fusion candidates.
        fstable_all_row_count (int): Row count for all fusion candidates table.
        all_candidates (int): Number of all fusion candidates.
        fstable_row_count (int): Row count for fusion candidates table.
        candidates (int): Number of fusion candidates.
        gene_gff3_2 (str): Path to the GFF3 file with gene annotations.
        gene_bed (str): Path to the gene BED file.
        all_gene_bed (str): Path to the BED file with all genes.
        gene_table (pd.DataFrame): DataFrame with gene annotations.
        gene_table_small (pd.DataFrame): Filtered DataFrame with gene annotations.
    """

    def __init__(self, *args, target_panel=None, **kwargs):
        self.target_panel = target_panel

        self.fusion_candidates = {}  # pd.DataFrame()
        self.fusion_candidates_all = {}  # pd.DataFrame()
        self.fstable_all = None
        self.fstable = None
        self.fstable_all_row_count = 0
        self.all_candidates = 0
        self.fstable_row_count = 0
        self.candidates = 0
        self.gene_gff3_2 = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)),
            "gencode.v45.basic.annotation.gff3",
        )

        if self.target_panel == "rCNS2":
            self.gene_bed = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)),
                "rCNS2_panel_name_uniq.bed",
            )
        elif self.target_panel == "AML":
            self.gene_bed = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)),
                "AML_panel_name_uniq.bed",
            )

        self.all_gene_bed = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)), "all_genes2.bed"
        )

        datafile = f"{self.target_panel}_data.csv.gz"

        if os.path.isfile(
            os.path.join(os.path.dirname(os.path.abspath(resources.__file__)), datafile)
        ):
            self.gene_table = pd.read_csv(
                os.path.join(
                    os.path.dirname(os.path.abspath(resources.__file__)), datafile
                )
            )
        else:
            # logging.info(
            #    f"This looks like the first time you have run the {self.target_panel} panel."
            # )
            # logging.info("Parsing GFF3")
            self.gene_table = gff3_parser.parse_gff3(
                self.gene_gff3_2, verbose=False, parse_attributes=True
            )

            self.gene_table_small = self.gene_table[
                self.gene_table["Type"].isin(["gene", "exon", "CDS"])
            ]
            self.gene_table_small = self.gene_table_small.drop(
                [
                    "Score",
                    "Phase",
                    "havana_gene",
                    "transcript_support_level",
                    "ont",
                    "transcript_id",
                    "hgnc_id",
                    "protein_id",
                    "havana_transcript",
                    "exon_number",
                    "artif_dupl",
                    "exon_id",
                    "gene_type",
                    "ID",
                    "gene_id",
                    "level",
                    "ccdsid",
                    "tag",
                    "transcript_name",
                    "Parent",
                ],
                axis=1,
            )
            self.gene_table_small.to_csv(
                os.path.join(
                    os.path.dirname(os.path.abspath(resources.__file__)), datafile
                ),
                index=False,
                compression="gzip",
            )
            self.gene_table = self.gene_table_small
        super().__init__(*args, **kwargs)

    def setup_ui(self) -> None:
        """
        Sets up the user interface for the Fusion Panel. This function creates the UI elements for the Fusion Panel.
        """
        if self.summary:
            with self.summary:
                with ui.card().classes('w-full p-4 mb-4'):
                    with ui.row().classes('w-full items-center justify-between'):
                        # Left side - Fusion Status
                        with ui.column().classes('gap-2'):
                            ui.label("Gene Fusion Analysis").classes('text-lg font-medium')
                            with ui.row().classes('items-center gap-2'):
                                ui.label(f"Panel: {self.target_panel}").classes('text-gray-600 font-medium')
                                ui.label("0").bind_text_from(
                                    self,
                                    "candidates",
                                    backward=lambda n: f"{n} high confidence"
                                ).classes('px-2 py-1 rounded bg-blue-100 text-blue-600')

                        # Right side - Additional metrics
                        with ui.column().classes('gap-2 text-right'):
                            ui.label("Analysis Details").classes('font-medium')
                            ui.label("0").bind_text_from(
                                self,
                                "all_candidates",
                                backward=lambda n: f"{n} low confidence fusions"
                            ).classes('text-gray-600')

                    # Bottom row - Information
                    with ui.row().classes('w-full mt-4 text-sm text-gray-500 justify-center'):
                        ui.label("Fusion candidates identified from reads with supplementary alignments")
        with ui.card().style("width: 100%"):
            ui.label("Gene Fusion Candidates").classes('text-sky-600 dark:text-white').style(
                "font-size: 150%; font-weight: 300"
            ).tailwind("drop-shadow", "font-bold")
            ui.label(
                "This panel identifies gene fusion candidates from the input bam files. "
                "The panel is split into two tabs, one for gene fusions between genes within the target panel and one for genome wide fusions. "
                "Fusions are identified on a streaming basis derived from reads with supplementary alignments. "
                "The plots are indicative of the presence of a fusion and should be interpreted with care. "
                "The tables show the reads that are indicative of a fusion."
            ).style("font-size: 125%; font-weight: 300")
            with ui.tabs().classes("w-full") as tabs:
                one = ui.tab("Within Target Fusions").style(
                    "font-size: 125%; font-weight: 300"
                )
                with one:
                    self.badge_one = (
                        ui.badge("0", color="red")
                        .bind_text_from(self, "candidates", backward=lambda n: f"{n}")
                        .props("floating rounded outline")
                    )
                two = ui.tab("Genome Wide Fusions").style(
                    "font-size: 125%; font-weight: 300"
                )
                with two:
                    self.badge_two = (
                        ui.badge("0", color="red")
                        .bind_text_from(
                            self, "all_candidates", backward=lambda n: f"{n}"
                        )
                        .props("floating rounded outline")
                    )
            with ui.tab_panels(tabs, value=one).classes("w-full"):
                with ui.tab_panel(one):
                    with ui.card().style("width: 100%"):
                        ui.label("Fusion Candidates (within targets)").style(
                            "font-size: 125%; font-weight: 300"
                        ).tailwind("drop-shadow", "font-bold")
                        ui.separator()
                        ui.label(
                            "Fusion Candidates are identified by looking for reads that map to two different genes from within the target panel. "
                            "Plots illustrate the alignments of reads to each region of the genome and individual reads are identified by colour. "
                            "These plots should be interpreted with care and are only indicative of the presence of a fusion."
                        ).style("font-size: 125%; font-weight: 300")
                        self.fusionplot = ui.row()
                        with self.fusionplot.classes("w-full"):
                            with ui.column().classes('gap-2'):
                                ui.label("Awaiting Fusion Plot").classes('text-lg font-medium')
                                ui.label("Plot will be displayed here when fusion data is available").classes('text-gray-600')
                        self.fusiontable = ui.row().classes("w-full")
                        with self.fusiontable:
                            with ui.column().classes('gap-2'):
                                ui.label("Awaiting Fusion Data").classes('text-lg font-medium')
                                ui.label("Table will be displayed here when fusion data is available").classes('text-gray-600')

                with ui.tab_panel(two):
                    with ui.card().style("width: 100%"):
                        ui.label("Fusion Candidates (genome wide)").style(
                            "font-size: 125%; font-weight: 300"
                        ).tailwind("drop-shadow", "font-bold")
                        ui.separator()
                        ui.label(
                            "Fusion Candidates are identified by looking for reads that map to at least one gene from within the target panel. "
                            "Plots illustrate the alignments of reads to each region of the genome and individual reads are identified by colour. "
                            "These plots should be interpreted with care and are only indicative of the presence of a fusion."
                        ).style("font-size: 125%; font-weight: 300")
                        self.fusionplot_all = ui.row()
                        with self.fusionplot_all.classes("w-full"):
                            with ui.column().classes('gap-2'):
                                ui.label("Awaiting Genome-wide Fusion Plot").classes('text-lg font-medium')
                                ui.label("Plot will be displayed here when genome-wide fusion data is available").classes('text-gray-600')
                        self.fusiontable_all = ui.row().classes("w-full")
                        with self.fusiontable_all:
                            with ui.column().classes('gap-2'):
                                ui.label("Awaiting Genome-wide Fusion Data").classes('text-lg font-medium')
                                ui.label("Table will be displayed here when genome-wide fusion data is available").classes('text-gray-600')

        if self.browse:
            self.show_previous_data()
        else:
            ui.timer(30, lambda: self.show_previous_data())

    def fusion_table_all(self) -> None:
        """
        Displays all fusion candidates in a table format.
        """
        # if not self.fusion_candidates_all.empty:
        if self.sampleID in self.fusion_candidates_all.keys():
            uniques_all = self.fusion_candidates_all[self.sampleID][
                "read_id"
            ].duplicated(keep=False)
            doubles_all = self.fusion_candidates_all[self.sampleID][uniques_all]
            counts_all = doubles_all.groupby("read_id")["col4"].transform("nunique")
            result_all = doubles_all[counts_all > 1]
            result_all.to_csv(
                os.path.join(
                    self.check_and_create_folder(self.output, self.sampleID),
                    "fusion_candidates_all.csv",
                ),
                index=False,
            )

    def update_fusion_table_all(self, result_all: pd.DataFrame) -> None:
        """
        Updates the UI table with all fusion candidates.

        Args:
            result_all (pd.DataFrame): DataFrame with all fusion candidates.
        """
        if result_all.shape[0] > self.fstable_all_row_count:
            self.fstable_all_row_count = result_all.shape[0]
            # self.fusiontable_all.clear()
            if not self.fstable_all:
                self.fusiontable_all.clear()
                with self.fusiontable_all:
                    self.fstable_all = (
                        ui.table.from_pandas(
                            result_all.sort_values(by=7).rename(
                                columns={
                                    0: "chromBED",
                                    1: "BS",
                                    2: "BE",
                                    3: "Gene",
                                    4: "chrom",
                                    5: "mS",
                                    6: "mE",
                                    7: "readID",
                                    8: "mapQ",
                                    9: "strand",
                                    10: "Read Map Start",
                                    11: "Read Map End",
                                    12: "Secondary",
                                    13: "Supplementary",
                                    14: "mapping span",
                                }
                            ),
                            pagination=25,
                        )
                        .props("dense")
                        .classes("w-full")
                        .style("height: 900px")
                        .style("font-size: 100%; font-weight: 300")
                    )
                    for col in self.fstable_all.columns:
                        col["sortable"] = True

                    with self.fstable_all.add_slot("top-right"):
                        with ui.input(placeholder="Search").props(
                            "type=search"
                        ).bind_value(self.fstable_all, "filter").add_slot("append"):
                            ui.icon("search")

                self.fusionplot_all.clear()
            else:
                self.fstable_all.update_rows(
                    result_all.sort_values(by=7)
                    .rename(
                        columns={
                            0: "chromBED",
                            1: "BS",
                            2: "BE",
                            3: "Gene",
                            4: "chrom",
                            5: "mS",
                            6: "mE",
                            7: "readID",
                            8: "mapQ",
                            9: "strand",
                            10: "Read Map Start",
                            11: "Read Map End",
                            12: "Secondary",
                            13: "Supplementary",
                            14: "mapping span",
                        }
                    )
                    .to_dict("records")
                )
                self.fstable_all.update()
                self.fusionplot_all.clear()

            result_all, goodpairs = _annotate_results(result_all)
            self.all_candidates = 0

            if not result_all.empty:
                with self.fusionplot_all.classes("w-full"):
                    gene_pairs = (
                        result_all[goodpairs].sort_values(by=7)["tag"].unique().tolist()
                    )
                    gene_pairs = [pair.split(", ") for pair in gene_pairs]
                    gene_groups_test = get_gene_network(gene_pairs)
                    gene_groups = []
                    for gene_group in gene_groups_test:
                        if (
                            len(
                                _get_reads(
                                    result_all[goodpairs][
                                        result_all[goodpairs][3].isin(gene_group)
                                    ]
                                )
                            )
                            > 1
                        ):
                            gene_groups.append(gene_group)
                    self.all_candidates = len(gene_groups)
                    with ui.row().classes("w-full"):
                        ui.select(
                            options=gene_groups,
                            with_input=True,
                            on_change=lambda e: show_gene_pair(
                                e.value, result_all, goodpairs
                            ),
                        ).classes("w-40")
                    with ui.row().classes("w-full"):
                        self.all_card = ui.card()
                        with self.all_card:
                            ui.label("Select gene pair to see results.").tailwind(
                                "drop-shadow", "font-bold"
                            )

                    def show_gene_pair(gene_group: str, result_all, goodpairs) -> None:
                        ui.notify(gene_group)
                        self.all_card.clear()
                        with self.all_card:
                            with ui.row():
                                ui.label(f"{gene_group}").tailwind(
                                    "drop-shadow", "font-bold"
                                )
                            with ui.row():
                                # for gene in gene_pair.split(", "):
                                #    reads = result_all[goodpairs].sort_values(by=7)[
                                #        result_all[goodpairs].sort_values(by=7)[3].str.contains(gene)
                                #    ]
                                reads = result_all[goodpairs][
                                    result_all[goodpairs][3].isin(gene_group)
                                ]
                                self.create_fusion_plot(gene_group, reads)

    def fusion_table(self) -> None:
        """
        Displays fusion candidates in a table format.
        """
        # if not self.fusion_candidates.empty:
        if self.sampleID in self.fusion_candidates.keys():
            uniques = self.fusion_candidates[self.sampleID]["read_id"].duplicated(
                keep=False
            )
            doubles = self.fusion_candidates[self.sampleID][uniques]
            counts = doubles.groupby("read_id")["col4"].transform("nunique")
            result = doubles[counts > 1]
            result.to_csv(
                os.path.join(
                    self.check_and_create_folder(self.output, self.sampleID),
                    "fusion_candidates_master.csv",
                ),
                index=False,
            )
        # self.update_fusion_table(result)

    def update_fusion_table(self, result: pd.DataFrame) -> None:
        """
        Updates the UI table with fusion candidates.

        Args:
            result (pd.DataFrame): DataFrame with fusion candidates.
        """
        if result.shape[0] > self.fstable_row_count:
            self.fstable_row_count = result.shape[0]
            if not self.fstable:
                self.fusiontable.clear()
                with self.fusiontable:
                    self.fstable = (
                        ui.table.from_pandas(
                            result.sort_values(by=7).rename(
                                columns={
                                    0: "chromBED",
                                    1: "BS",
                                    2: "BE",
                                    3: "Gene",
                                    4: "chrom",
                                    5: "mS",
                                    6: "mE",
                                    7: "readID",
                                    8: "mapQ",
                                    9: "strand",
                                    10: "Read Map Start",
                                    11: "Read Map End",
                                    12: "Secondary",
                                    13: "Supplementary",
                                    14: "mapping span",
                                }
                            ),
                            pagination=25,
                        )
                        .props("dense")
                        .classes("w-full")
                        .style("height: 900px")
                        .style("font-size: 100%; font-weight: 300")
                    )
                    for col in self.fstable.columns:
                        col["sortable"] = True

                    with self.fstable.add_slot("top-right"):
                        with ui.input(placeholder="Search").props(
                            "type=search"
                        ).bind_value(self.fstable, "filter").add_slot("append"):
                            ui.icon("search")

                self.fusionplot.clear()
            else:
                self.fstable.update_rows(
                    result.sort_values(by=7)
                    .rename(
                        columns={
                            0: "chromBED",
                            1: "BS",
                            2: "BE",
                            3: "Gene",
                            4: "chrom",
                            5: "mS",
                            6: "mE",
                            7: "readID",
                            8: "mapQ",
                            9: "strand",
                            10: "Read Map Start",
                            11: "Read Map End",
                            12: "Secondary",
                            13: "Supplementary",
                            14: "mapping span",
                        }
                    )
                    .to_dict("records")
                )

                self.fusionplot.clear()

            result, goodpairs = _annotate_results(result)
            # self.candidates = result[goodpairs].sort_values(by=7)["tag"].nunique()
            self.candidates = 0

            if not (result.empty):
                with self.fusionplot.classes("w-full"):
                    gene_pairs = (
                        result[goodpairs].sort_values(by=7)["tag"].unique().tolist()
                    )
                    gene_pairs = [pair.split(",") for pair in gene_pairs]
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
                    self.candidates = len(gene_groups)
                    with ui.row().classes("w-full"):
                        ui.select(
                            options=gene_groups,
                            with_input=True,
                            on_change=lambda e: show_gene_pair(
                                e.value, result, goodpairs
                            ),
                        ).classes("w-40")
                    with ui.row().classes("w-full"):
                        self.card = ui.card()
                        with self.card:
                            ui.label("Select gene pair to see results.").tailwind(
                                "drop-shadow", "font-bold"
                            )

                    def show_gene_pair(gene_group: str, result, goodpairs) -> None:
                        ui.notify(gene_group)
                        self.card.clear()
                        with self.card:
                            with ui.row():
                                ui.label(f"{gene_group}").tailwind(
                                    "drop-shadow", "font-bold"
                                )
                            with ui.row():
                                # for gene in gene_pair.split(", "):
                                #    reads = result_all[goodpairs].sort_values(by=7)[
                                #        result_all[goodpairs].sort_values(by=7)[3].str.contains(gene)
                                #    ]
                                reads = result[goodpairs][
                                    result[goodpairs][3].isin(gene_group)
                                ]
                                self.create_fusion_plot(gene_group, reads)

            """
            with self.fusionplot.classes("w-full"):
                for gene_pair in result[goodpairs].sort_values(by=7)["tag"].unique():
                    with ui.card():
                        with ui.row():
                            ui.label(f"{gene_pair}").tailwind(
                                "drop-shadow", "font-bold"
                            )
                        with ui.row():
                            for gene in result[
                                goodpairs & result[goodpairs]["tag"].eq(gene_pair)
                            ][3].unique():
                                self.create_fusion_plot(
                                    gene,
                                    result[goodpairs].sort_values(by=7)[
                                        result[goodpairs].sort_values(by=7)[3].eq(gene)
                                    ],
                                )
            """

    def create_fusion_plot(self, title: str, reads: pd.DataFrame) -> None:
        """
        Creates a fusion plot to illustrate the alignments of reads to each region of the genome.

        Args:
            title (str): Title of the plot.
            reads (pd.DataFrame): DataFrame with reads to be plotted.
        """
        with ui.card().classes("w-full no-shadow border-[2px]"):
            result = _get_reads(reads)
            df = reads

            # Function to rank overlapping ranges
            def rank_overlaps(df, start_col, end_col):
                # Sort by start and end columns
                df = df.sort_values(by=["gene", start_col, end_col]).reset_index(
                    drop=True
                )
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

            lines = df.sort_values(by=["id", "read_start"]).reset_index(drop=True)

            # Function to assign occurrence number
            def assign_occurrence(group):
                group["Occurrence"] = range(1, len(group) + 1)
                return group

            # Apply the function to each group
            lines = lines.groupby("id").apply(assign_occurrence).reset_index(drop=True)

            # Function to find join coordinates
            def find_joins(group):
                group = group.reset_index(drop=True)
                for i in range(len(group)):
                    if i < len(group) - 1:
                        if (
                            group.loc[i, "Occurrence"] == 1
                        ):  # The first read in the sequence of read mappings
                            group.loc[i, "Join_Gene"] = group.loc[i + 1, "gene"]
                            group.loc[i, "Join_Start"] = group.loc[
                                i, "start2"
                            ]  # Reference end of current read
                            group.loc[i, "Join_Chromosome"] = group.loc[
                                i + 1, "chromosome2"
                            ]
                            group.loc[i, "spanB"] = group.loc[i + 1, "span"]
                            group.loc[i, "start3"] = group.loc[i + 1, "start2"]
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
            lines = lines.groupby("id").apply(find_joins).reset_index(drop=True)
            # Remove rows containing NA values
            lines = lines.dropna()

            # Dict to store ax

            axdict = {}

            if len(result) > 1:
                gene_table = self.gene_table
                with ui.pyplot(figsize=(19, 5)).classes("w-full"):  # figsize=(20, 5)):
                    num_plots = 2 * len(result)
                    num_cols = len(result)  # 2  # Number of columns in the subplot grid
                    num_rows = (
                        num_plots + num_cols - 1
                    ) // num_cols  # Calculate the number of rows needed

                    for i, ax in enumerate(range(num_plots), start=1):
                        plt.subplot(num_rows, num_cols, i)
                        row, col = divmod(i - 1, num_cols)
                        # plt.title(f'Subplot {i} (Row {row + 1}, Col {col + 1})')
                        data = result.iloc[col]

                        chrom = data["chromosome"]
                        start = data["start"]
                        end = data["end"]

                        # start = df[df["chromosome"].eq(chrom)]["start2"].min()-100_000
                        # end = df[df["chromosome"].eq(chrom)]["end2"].max()+100_000

                        def human_readable_format(x, pos):
                            return f"{x / 1e6:.2f}"  # Mb'

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
                                                "family": "sans",
                                                "color": "black",
                                                "fontsize": 8,
                                            },
                                        )
                                    )

                            for index, row in (
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
                                        start=int(row["Start"]),
                                        end=int(row["End"]),
                                        strand=STRAND[row["Strand"]],
                                        thickness=4,
                                        color="#C0C0C0",
                                        # label=row["gene_name"],
                                        # fontdict = {
                                        #    'family': 'sans',
                                        #    'color':  'black',
                                        #    'fontsize': 8,
                                        # }
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

                            # plt.tight_layout()

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
                            )  # Set font size for x-axis labels in ax0
                            ax.set_xlabel(
                                f'Position (Mb) - {chrom} - {data["gene"]}', fontsize=10
                            )
                            ax.set_title(f'{data["gene"]}')

                            axdict[data["gene"]] = ax
                            # plt.tight_layout()
                    gene_counter = Counter()

                    for index, row in lines.iterrows():
                        # Coordinates in data space of each subplot

                        xyB = (row["Join_Start"], row["rank"])  # Point in subplot 1
                        xyA = (row["Join_End"], row["rankB"])

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
                            # Add the arc connection to the second subplot (ax2)
                            # axdict[row["Join_Gene"]].add_artist(con)
                            # plt.tight_layout()
                        """

    async def process_bam(self, bamfile: str, timestamp: str) -> None:
        """
        Processes a BAM file and identifies fusion candidates.

        Args:
            bamfile (str): Path to the BAM file to process.
            timestamp (str): Timestamp for the analysis.
        """
        tempreadfile = tempfile.NamedTemporaryFile(
            dir=self.check_and_create_folder(self.output, self.sampleID), suffix=".txt"
        )
        tempbamfile = tempfile.NamedTemporaryFile(
            dir=self.check_and_create_folder(self.output, self.sampleID), suffix=".bam"
        )
        tempmappings = tempfile.NamedTemporaryFile(
            dir=self.check_and_create_folder(self.output, self.sampleID), suffix=".txt"
        )
        tempallmappings = tempfile.NamedTemporaryFile(
            dir=self.check_and_create_folder(self.output, self.sampleID), suffix=".txt"
        )

        try:
            fusion_candidates, fusion_candidates_all = await run.cpu_bound(
                fusion_work,
                self.threads,
                bamfile,
                self.gene_bed,
                self.all_gene_bed,
                tempreadfile.name,
                tempbamfile.name,
                tempmappings.name,
                tempallmappings.name,
            )

            if fusion_candidates is not None:
                if self.sampleID not in self.fusion_candidates.keys():
                    self.fusion_candidates[self.sampleID] = fusion_candidates
                else:
                    self.fusion_candidates[self.sampleID] = pd.concat(
                        [self.fusion_candidates[self.sampleID], fusion_candidates]
                    ).reset_index(drop=True)

            if fusion_candidates_all is not None:
                if self.sampleID not in self.fusion_candidates_all.keys():
                    self.fusion_candidates_all[self.sampleID] = fusion_candidates_all
                else:
                    self.fusion_candidates_all[self.sampleID] = pd.concat(
                        [
                            self.fusion_candidates_all[self.sampleID],
                            fusion_candidates_all,
                        ]
                    ).reset_index(drop=True)

            self.fusion_table()
            self.fusion_table_all()

        except Exception:
            # logging.error(f"Error processing BAM file: {e}")
            raise
        finally:
            self.running = False

    def show_previous_data(self) -> None:
        """
        Displays previously analyzed data from the specified output folder.

        Args:
            output (str): Path to the output folder.
        """
        if not self.browse:
            for item in app.storage.general[self.mainuuid]:
                if item == "sample_ids":
                    for sample in app.storage.general[self.mainuuid][item]:
                        self.sampleID = sample
            output = self.output
        if self.browse:
            output = self.check_and_create_folder(self.output, self.sampleID)

        if self.check_file_time(os.path.join(output, "fusion_candidates_master.csv")):
            try:
                fusion_candidates = pd.read_csv(
                    os.path.join(output, "fusion_candidates_master.csv"),
                    dtype=str,
                    # names=["index", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, "diff"],
                    # names=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, "diff"],
                    header=None,
                    skiprows=1,
                )
                self.update_fusion_table(fusion_candidates)
            except pd.errors.EmptyDataError:
                pass

        if self.check_file_time(os.path.join(output, "fusion_candidates_all.csv")):
            try:
                fusion_candidates_all = pd.read_csv(
                    os.path.join(output, "fusion_candidates_all.csv"),
                    dtype=str,
                    # names=["index", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, "diff"],
                    # names=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, "diff"],
                    header=None,
                    skiprows=1,
                )
                self.update_fusion_table_all(fusion_candidates_all)
            except pd.errors.EmptyDataError:
                pass


def test_me(
    port: int,
    threads: int,
    watchfolder: str,
    output: str,
    reload: bool = False,
    browse: bool = False,
) -> None:
    """
    Sets up and runs the gene fusion identification application.

    Args:
        port (int): Port number for the server.
        threads (int): Number of threads to use for processing.
        watchfolder (str): Path to the folder to watch for new BAM files.
        output (str): Path to the output directory.
        reload (bool): Flag to reload the application on changes.
        browse (bool): Flag to enable browsing historic data.
    """
    my_connection = None
    with theme.frame("Fusion Gene Identification.", my_connection):
        TestObject = FusionObject(threads, output, progress=True)
    if not browse:
        path = watchfolder
        searchdirectory = os.fsencode(path)
        for root, d_names, f_names in os.walk(searchdirectory):
            directory = os.fsdecode(root)
            for f in f_names:
                filename = os.fsdecode(f)
                if filename.endswith(".bam"):
                    TestObject.add_bam(os.path.join(directory, filename))
    else:
        TestObject.progress_trackers.visible = False
        TestObject.show_previous_data(output)
    ui.run(port=port, reload=reload)


@click.command()
@click.option(
    "--port",
    default=12345,
    help="Port for GUI",
)
@click.option("--threads", default=4, help="Number of threads available.")
@click.argument(
    "watchfolder",
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, resolve_path=True, path_type=Path
    ),
    required=False,
)
@click.argument(
    "output",
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, resolve_path=True, path_type=Path
    ),
    required=False,
)
@click.option(
    "--browse",
    is_flag=True,
    show_default=True,
    default=False,
    help="Browse Historic Data.",
)
def main(port, threads, watchfolder, output, browse) -> None:
    """
    CLI entry point for running the gene fusion identification app.

    Args:
        port (int): The port to serve the app on.
        threads (int): Number of threads available for processing.
        watchfolder (str): Directory to watch for new BAM files.
        output (str): Directory to save output files.
        browse (bool): Enable browsing historic data.
    """
    if browse:
        click.echo("Browse mode is enabled. Only the output folder is required.")
        test_me(
            port=port,
            reload=False,
            threads=threads,
            watchfolder=None,
            output=watchfolder,
            browse=browse,
        )
    else:
        click.echo(f"Watchfolder: {watchfolder}, Output: {output}")
        if watchfolder is None or output is None:
            click.echo("Watchfolder and output are required when --browse is not set.")
            sys.exit(1)
        test_me(
            port=port,
            reload=False,
            threads=threads,
            watchfolder=watchfolder,
            output=output,
            browse=browse,
        )


if __name__ in {"__main__", "__mp_main__"}:
    print("GUI launched by auto-reload function.")
    main()
