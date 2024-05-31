from cnsmeth.subpages.base_analysis import BaseAnalysis
import os
import sys
import gff3_parser
import tempfile
import random
import asyncio
import pandas as pd
import click


from nicegui import ui, background_tasks, run
from cnsmeth import theme, resources
from dna_features_viewer import GraphicFeature, GraphicRecord
from pathlib import Path

import matplotlib

matplotlib.use("agg")
from matplotlib import pyplot as plt


os.environ["CI"] = "1"
STRAND = {"+": 1, "-": -1}


def fusion_work(
    threads,
    bamfile,
    gene_bed,
    all_gene_bed,
    tempreadfile,
    tempbamfile,
    tempmappings,
    tempallmappings,
):
    fusion_candidates = None
    fusion_candidates_all = None
    os.system(
        f"samtools view -@{threads} -L {gene_bed} -d SA {bamfile} | cut -f1 > {tempreadfile.name}"
    )
    if os.path.getsize(tempreadfile.name) > 0:
        os.system(
            f"samtools view -@{threads} -N {tempreadfile.name} -o {tempbamfile.name} {bamfile}"
        )
        if os.path.getsize(tempbamfile.name) > 0:
            os.system(
                f"bedtools intersect -a {gene_bed} -b {tempbamfile.name} -wa -wb > {tempmappings.name}"
            )
            os.system(
                f"bedtools intersect -a {all_gene_bed} -b {tempbamfile.name} -wa -wb > {tempallmappings.name}"
            )
            if os.path.getsize(tempmappings.name) > 0:
                fusion_candidates = pd.read_csv(
                    tempmappings.name, sep="\t", header=None
                )
                # Filter to only include good mappings
                fusion_candidates = fusion_candidates[fusion_candidates[8].gt(50)]
                # Filter to only keep reads that map to 1 kb or more of the reference.
                fusion_candidates["diff"] = fusion_candidates[6] - fusion_candidates[5]
                fusion_candidates = fusion_candidates[
                    fusion_candidates["diff"].gt(1000)
                ]
            if os.path.getsize(tempallmappings.name) > 0:
                fusion_candidates_all = pd.read_csv(
                    tempallmappings.name, sep="\t", header=None
                )

                fusion_candidates_all = fusion_candidates_all[
                    fusion_candidates_all[8].gt(59)
                ]

                fusion_candidates_all["diff"] = (
                    fusion_candidates_all[6] - fusion_candidates_all[5]
                )

                fusion_candidates_all = fusion_candidates_all[
                    fusion_candidates_all["diff"].gt(1000)
                ]
    return fusion_candidates, fusion_candidates_all


class Fusion_object(BaseAnalysis):
    def __init__(self, *args, target_panel=None, **kwargs):
        self.target_panel = target_panel

        self.fusion_candidates = pd.DataFrame()
        self.fusion_candidates_all = pd.DataFrame()
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
            print(
                f"This looks like the first time you have run the {self.target_panel} panel."
            )
            print("Parsing GFF3")
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

    def setup_ui(self):
        """
        Function to setup the UI for the Fusion Panel. This function creates the UI elements for the Fusion Panel.
        """
        if self.summary:
            with self.summary:
                ui.label(f"Fusion Candidates - using panel {self.target_panel}")
                with ui.row():
                    ui.label("0").bind_text_from(
                        self,
                        "candidates",
                        backward=lambda n: f"{n} high confidence fusions observed.",
                    )
                    ui.label("0").bind_text_from(
                        self,
                        "all_candidates",
                        backward=lambda n: f" {n} low confidence fusions observed.",
                    )
        with ui.card().style("width: 100%"):
            ui.label("Gene Fusion Candidates").style(
                "color: #6E93D6; font-size: 150%; font-weight: 300"
            ).tailwind("drop-shadow", "font-bold")
            ui.label(
                "This panel identifies gene fusion candidates from the input bam files. "
                "The panel is split into two tabs, one for gene fusions between genes within the target panel and one for genome wide fusions. "
                "Fusions are identified on a streaming basis derived from reads with supplementary alignments. "
                "The plots are indicative of the presence of a fusion and should be interpreted with care. "
                "The tables show the reads that are indicative of a fusion."
            ).style("color: #000000; font-size: 125%; font-weight: 300")
            with ui.tabs().classes("w-full") as tabs:
                one = ui.tab("Within Target Fusions").style(
                    "color: #000000; font-size: 125%; font-weight: 300"
                )
                with one:
                    self.badge_one = (
                        ui.badge("0", color="red")
                        .bind_text_from(self, "candidates", backward=lambda n: f"{n}")
                        .props("floating rounded outline")
                    )
                two = ui.tab("Genome Wide Fusions").style(
                    "color: #000000; font-size: 125%; font-weight: 300"
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
                            "color: #000000; font-size: 125%; font-weight: 300"
                        ).tailwind("drop-shadow", "font-bold")
                        ui.separator()
                        ui.label(
                            "Fusion Candidates are identified by looking for reads that map to two different genes from within the target panel. "
                            "Plots illustrate the alignments of reads to each region of the genome and individual reads are identified by colour. "
                            "These plots should be interpreted with care and are only indicative of the presence of a fusion."
                        ).style("color: #000000; font-size: 125%; font-weight: 300")
                        self.fusionplot = ui.row()
                        with self.fusionplot.classes("w-full"):
                            ui.label("Plot not yet available.").style(
                                "color: #000000; font-size: 125%; font-weight: 300"
                            )
                        self.fusiontable = ui.row().classes("w-full")
                        with self.fusiontable:
                            ui.label("Table not yet available.").style(
                                "color: #000000; font-size: 125%; font-weight: 300"
                            )

                with ui.tab_panel(two):
                    with ui.card().style("width: 100%"):
                        ui.label("Fusion Candidates (genome wide)").style(
                            "color: #000000; font-size: 125%; font-weight: 300"
                        ).tailwind("drop-shadow", "font-bold")
                        ui.separator()
                        ui.label(
                            "Fusion Candidates are identified by looking for reads that map to at least one gene from within the target panel. "
                            "Plots illustrate the alignments of reads to each region of the genome and individual reads are identified by colour. "
                            "These plots should be interpreted with care and are only indicative of the presence of a fusion."
                        ).style("color: #000000; font-size: 125%; font-weight: 300")
                        self.fusionplot_all = ui.row()
                        with self.fusionplot_all.classes("w-full"):
                            ui.label("Plot not yet available.").style(
                                "color: #000000; font-size: 125%; font-weight: 300"
                            )
                        self.fusiontable_all = ui.row().classes("w-full")
                        with self.fusiontable_all:
                            ui.label("Table not yet available.").style(
                                "color: #000000; font-size: 125%; font-weight: 300"
                            )
        if self.browse:
            self.show_previous_data(self.output)
        else:
            ui.timer(5, lambda: self.show_previous_data(self.output))

    def fusion_table_all(self):
        if not self.fusion_candidates_all.empty:
            uniques_all = self.fusion_candidates_all[7].duplicated(keep=False)
            doubles_all = self.fusion_candidates_all[uniques_all]
            counts_all = doubles_all.groupby(7)[3].transform("nunique")
            result_all = doubles_all[counts_all > 1]
            result_all.to_csv(
                os.path.join(self.output, "fusion_candidates_all.csv"), index=False
            )
            # self.update_fusion_table_all(result_all)

    def update_fusion_table_all(self, result_all):
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
                        }
                    )
                    .to_dict("records")
                )
                self.fstable_all.update()
                self.fusionplot_all.clear()

            result_all, goodpairs = self._annotate_results(result_all)
            self.all_candidates = (
                result_all[goodpairs].sort_values(by=7)["tag"].nunique()
            )

            if not result_all.empty:
                with self.fusionplot_all.classes("w-full"):
                    for gene_pair in (
                        result_all[goodpairs].sort_values(by=7)["tag"].unique()
                    ):
                        with ui.card():
                            with ui.row():
                                ui.label(f"{gene_pair}").tailwind(
                                    "drop-shadow", "font-bold"
                                )
                            with ui.row():
                                for gene in result_all[
                                    goodpairs
                                    & result_all[goodpairs]["tag"].eq(gene_pair)
                                ][3].unique():
                                    # self.create_fusion_plot(
                                    #    gene,
                                    #    result_all[goodpairs].sort_values(by=7)[
                                    #        result_all[goodpairs]
                                    #        .sort_values(by=7)[3]
                                    #        .eq(gene)
                                    #    ],
                                    # )
                                    title = gene
                                    reads = result_all[goodpairs].sort_values(by=7)[
                                        result_all[goodpairs]
                                        .sort_values(by=7)[3]
                                        .eq(gene)
                                    ]
                                    with ui.card().classes("no-shadow border-[2px]"):
                                        with ui.pyplot(figsize=(16, 2)):
                                            ax1 = plt.gca()
                                            features = []
                                            first_index = 0
                                            sequence_length = 100
                                            x_label = ""
                                            for index, row in self.gene_table[
                                                self.gene_table["gene_name"].eq(
                                                    title.strip()
                                                )
                                            ].iterrows():
                                                if row["Type"] == "gene":
                                                    x_label = row["Seqid"]
                                                    features.append(
                                                        GraphicFeature(
                                                            start=int(row["Start"]),
                                                            end=int(row["End"]),
                                                            strand=STRAND[
                                                                row["Strand"]
                                                            ],
                                                            thickness=4,
                                                            color="#ffd700",
                                                        )
                                                    )
                                                    first_index = (
                                                        int(row["Start"]) - 1000
                                                    )
                                                    sequence_length = (
                                                        int(row["End"])
                                                        - int(row["Start"])
                                                        + 2000
                                                    )
                                                if (
                                                    row["Type"] == "CDS"
                                                    and row["transcript_type"]
                                                    == "protein_coding"
                                                ):
                                                    features.append(
                                                        GraphicFeature(
                                                            start=int(row["Start"]),
                                                            end=int(row["End"]),
                                                            strand=STRAND[
                                                                row["Strand"]
                                                            ],
                                                            color="#ffcccc",
                                                        )
                                                    )
                                            record = GraphicRecord(
                                                sequence_length=sequence_length,
                                                first_index=first_index,
                                                features=features,
                                            )
                                            ax1.set_title(f"{title} - {x_label}")
                                            record.plot(ax=ax1)

                                        with ui.pyplot(figsize=(16, 1)):
                                            ax = plt.gca()
                                            features = []
                                            for index, row in reads.sort_values(
                                                by=7
                                            ).iterrows():
                                                features.append(
                                                    GraphicFeature(
                                                        start=int(row[5]),
                                                        end=int(row[6]),
                                                        strand=STRAND[row[9]],
                                                        color=row["Color"],
                                                    )
                                                )
                                            record = GraphicRecord(
                                                sequence_length=sequence_length,
                                                first_index=first_index,
                                                features=features,
                                            )
                                            record.plot(
                                                ax=ax, with_ruler=False, draw_line=True
                                            )

    def fusion_table(self):
        if not self.fusion_candidates.empty:
            uniques = self.fusion_candidates[7].duplicated(keep=False)
            doubles = self.fusion_candidates[uniques]
            counts = doubles.groupby(7)[3].transform("nunique")
            result = doubles[counts > 1]
            result.to_csv(
                os.path.join(self.output, "fusion_candidates_master.csv"), index=False
            )
        # self.update_fusion_table(result)

    def update_fusion_table(self, result):
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
                        }
                    )
                    .to_dict("records")
                )

                self.fusionplot.clear()

            result, goodpairs = self._annotate_results(result)
            self.candidates = result[goodpairs].sort_values(by=7)["tag"].nunique()
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

    def create_fusion_plot(self, title, reads):
        """
        Function to create a fusion plot. This function creates a plot to illustrate the alignments of reads to each region of the genome.
        :param title: The title of the plot.
        :param reads: The reads to be plotted.
        """
        with ui.card().classes("no-shadow border-[2px]"):
            with ui.pyplot(figsize=(16, 2)):
                ax1 = plt.gca()
                features = []
                first_index = 0
                sequence_length = 100
                x_label = ""
                for index, row in self.gene_table[
                    self.gene_table["gene_name"].eq(title.strip())
                ].iterrows():
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
                    if (
                        row["Type"] == "CDS"
                        and row["transcript_type"] == "protein_coding"
                    ):
                        features.append(
                            GraphicFeature(
                                start=int(row["Start"]),
                                end=int(row["End"]),
                                strand=STRAND[row["Strand"]],
                                color="#ffcccc",
                            )
                        )
                record = GraphicRecord(
                    sequence_length=sequence_length,
                    first_index=first_index,
                    features=features,
                )
                ax1.set_title(f"{title} - {x_label}")
                record.plot(ax=ax1)

            with ui.pyplot(figsize=(16, 1)):
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
                    sequence_length=sequence_length,
                    first_index=first_index,
                    features=features,
                )
                record.plot(ax=ax, with_ruler=False, draw_line=True)

    async def process_bam(self, bamfile, timestamp):
        """
        Function to process a bam file and identify fusion candidates.
        :param bamfile: The path to the bam file to process.
        """
        tempreadfile = tempfile.NamedTemporaryFile(dir=self.output, suffix=".txt")
        tempbamfile = tempfile.NamedTemporaryFile(dir=self.output, suffix=".bam")
        tempmappings = tempfile.NamedTemporaryFile(dir=self.output, suffix=".txt")
        tempallmappings = tempfile.NamedTemporaryFile(dir=self.output, suffix=".txt")


        #ToDo: Move this to another process

        fusion_candidates, fusion_candidates_all = fusion_work(self.threads, bamfile, self.gene_bed, self.all_gene_bed, tempreadfile, tempbamfile, tempmappings, tempallmappings)

        if fusion_candidates is not None:
            if self.fusion_candidates.empty:
                self.fusion_candidates = fusion_candidates
            else:
                self.fusion_candidates = pd.concat(
                    [self.fusion_candidates, fusion_candidates]
                ).reset_index(drop=True)

        if fusion_candidates_all is not None:
            if self.fusion_candidates_all.empty:
                self.fusion_candidates_all = fusion_candidates_all
            else:
                self.fusion_candidates_all = pd.concat(
                    [self.fusion_candidates_all, fusion_candidates_all]
                ).reset_index(drop=True)

        self.fusion_table()

        self.fusion_table_all()

        #await asyncio.sleep(0.1)
        self.running = False

    def _generate_random_color(self):
        """
        Function to generate a random color for use in plotting.
        """
        return "#{:06x}".format(random.randint(0, 0xFFFFFF))

    def _annotate_results(self, result):
        result_copy = result.copy()
        lookup = result_copy.groupby(7)[3].agg(lambda x: ",".join(set(x)))
        tags = result_copy[7].map(lookup.get)
        result_copy.loc[:, "tag"] = tags
        result = result_copy
        colors = result.groupby(7).apply(lambda x: self._generate_random_color())
        result["Color"] = result[7].map(colors.get)
        goodpairs = result.groupby("tag")[7].transform("nunique") > 1
        return result, goodpairs

    def count_shared_values(self, group):
        return group[7].nunique()

    def show_previous_data(self, output):
        if os.path.isfile(os.path.join(output, "fusion_candidates_master.csv")):
            fusion_candidates = pd.read_csv(
                os.path.join(output, "fusion_candidates_master.csv"),
                dtype=str,
                names=["index", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, "diff"],
                header=None,
                skiprows=1,
            )
            self.update_fusion_table(fusion_candidates)
        if os.path.isfile(os.path.join(output, "fusion_candidates_all.csv")):
            fusion_candidates_all = pd.read_csv(
                os.path.join(output, "fusion_candidates_all.csv"),
                dtype=str,
                names=["index", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, "diff"],
                header=None,
                skiprows=1,
            )
            self.update_fusion_table_all(fusion_candidates_all)


def test_me(
    port: int,
    threads: int,
    watchfolder: str,
    output: str,
    reload: bool = False,
    browse: bool = False,
):
    my_connection = None
    with theme.frame("Fusion Gene Identification.", my_connection):
        TestObject = Fusion_object(threads, output, progress=True)
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
        # print("Browse mode not implemented.")
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
def main(port, threads, watchfolder, output, browse):
    """
    Helper function to run the app.
    :param port: The port to serve the app on.
    :param reload: Should we reload the app on changes.
    :return:
    """
    if browse:
        # Handle the case when --browse is set
        click.echo("Browse mode is enabled. Only the output folder is required.")
        test_me(
            port=port,
            reload=False,
            threads=threads,
            # simtime=simtime,
            watchfolder=None,
            output=watchfolder,
            # sequencing_summary=sequencing_summary,
            # showerrors=showerrors,
            browse=browse,
            # exclude=exclude,
        )
        # Your logic for browse mode
    else:
        # Handle the case when --browse is not set
        click.echo(f"Watchfolder: {watchfolder}, Output: {output}")
        if watchfolder is None or output is None:
            click.echo("Watchfolder and output are required when --browse is not set.")
            sys.exit(1)
        test_me(
            port=port,
            reload=False,
            threads=threads,
            # simtime=simtime,
            watchfolder=watchfolder,
            output=output,
            # sequencing_summary=sequencing_summary,
            # showerrors=showerrors,
            browse=browse,
            # exclude=exclude,
        )


if __name__ in {"__main__", "__mp_main__"}:
    print("GUI launched by auto-reload function.")
    main()
