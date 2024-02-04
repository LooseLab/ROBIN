"""
This module contains the Fusion_Panel class which is used to create the Fusion Panel in the MethClass NiceGUI.
The Fusion Panel displays fusion candidates identified during a sequencing run.
Candidates are identified by filtering mapping files for reads with supplementary alignments indicative of a candidate fusion.
Fusions are filtered to those that occur between genes within the target panel and those that occur between genes outside of the target panel.
"""

import os
import gff3_parser
import tempfile
import random
import natsort
import pandas as pd


from nicegui import ui
from methnicegui import theme, resources
from dna_features_viewer import GraphicFeature, GraphicRecord

from matplotlib import pyplot as plt

os.environ["CI"] = "1"
STRAND = {"+": 1, "-": -1}


class Fusion_Panel:
    def __init__(self):
        """
        Constructor for the Fusion_Panel class. This class is used to create the Fusion Panel in the MethClass NiceGUI.
        We preload a bed file for the analysis and process the gff3 file to create a pandas dataframe of genes.
        This dataframe is stored for future use.
        """
        self.gene_bed = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)), "unique_genes.bed"
        )
        self.all_gene_bed = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)), "all_genes2.bed"
        )

        self.gene_gff3_2 = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)),
            "gencode.v45.basic.annotation.gff3",
        )

        if os.path.isfile(
            os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)), "data.csv.gz"
            )
        ):
            self.gene_table = pd.read_csv(
                os.path.join(
                    os.path.dirname(os.path.abspath(resources.__file__)), "data.csv.gz"
                )
            )
        else:
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
                    os.path.dirname(os.path.abspath(resources.__file__)), "data.csv.gz"
                ),
                index=False,
                compression="gzip",
            )
            self.gene_table = self.gene_table_small

    def setup_ui(self):
        """
        Function to setup the UI for the Fusion Panel. This function creates the UI elements for the Fusion Panel.
        """
        with ui.card().style("width: 100%"):
            with ui.tabs().classes("w-full") as tabs:
                one = ui.tab("Within Target Fusions")
                two = ui.tab("Genome Wide Fusions")
            with ui.tab_panels(tabs, value=one).classes("w-full"):
                with ui.tab_panel(one):
                    with ui.card().style("width: 100%"):
                        ui.label("Fusion Candidates (within targets)").tailwind(
                            "drop-shadow", "font-bold"
                        )
                        ui.separator()
                        ui.label(
                            "Fusion Candidates are identified by looking for reads that map to two different genes from within the target panel. "
                            "Plots illustrate the alignments of reads to each region of the genome and individual reads are identified by colour. "
                            "These plots should be interpreted with care and are only indicative of the presence of a fusion."
                        )
                        self.fusionplot = ui.row()
                        with self.fusionplot.classes("w-full"):
                            ui.label("Plot not yet available.")
                        self.fusiontable = ui.row()
                        with self.fusiontable:
                            ui.label("Table not yet available.")

                with ui.tab_panel(two):
                    with ui.card().style("width: 100%"):
                        ui.label("Fusion Candidates (genome wide)").tailwind(
                            "drop-shadow", "font-bold"
                        )
                        ui.separator()
                        ui.label(
                            "Fusion Candidates are identified by looking for reads that map to at least one gene from within the target panel. "
                            "Plots illustrate the alignments of reads to each region of the genome and individual reads are identified by colour. "
                            "These plots should be interpreted with care and are only indicative of the presence of a fusion."
                        )
                        self.fusionplot_all = ui.row()
                        with self.fusionplot_all.classes("w-full"):
                            ui.label("Plot not yet available.")
                        self.fusiontable_all = ui.row()
                        with self.fusiontable_all:
                            ui.label("Table not yet available.")

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
                plt.title(title)
                for index, row in self.gene_table[
                    self.gene_table["gene_name"].eq(title.strip())
                ].iterrows():
                    if row["Type"] == "gene":
                        x_label = row[0]
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
                ax1.set_xlabel(x_label)
                record.plot(ax=ax1)

            with ui.pyplot(figsize=(16, 1)):
                ax = plt.gca()
                features = []
                x_label = ""
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
                ax.set_xlabel(x_label)
                record.plot(ax=ax, with_ruler=False, draw_line=True)

    def parse_bams(self, bampath):
        """
        Function to parse bamfiles. This function parses bamfiles to identify fusion candidates.
        The function subsets the bamfiles with the bed file and only keeps reads with supplementary alignments.
        The function then merges the newly formed bamfiles and looks for fusions between the target gene panel and the merged bamfile.
        :param bampath: The path to the bamfiles.
        """
        bamfiles = natsort.natsorted(os.listdir(bampath))
        # Check if all bamfiles have already been subset - if not subset them with the bed file and only keep reads with supplementary alignments
        for file in bamfiles:
            if file.endswith(".bam"):
                if file[0].isdigit():
                    if not os.path.exists(os.path.join(bampath, f"subset_{file}")):
                        subset_file = f"subset_{file}"
                        tempreadfile = tempfile.NamedTemporaryFile()
                        os.system(
                            f"samtools view -L {self.gene_bed} -d SA {os.path.join(bampath, file)} | cut -f1 > {tempreadfile.name}"
                        )
                        os.system(
                            f"samtools view --write-index -N {tempreadfile.name} -o {os.path.join(bampath, subset_file)} {os.path.join(bampath, file)}"
                        )

        # Now we merge the newly formed bamfiles:
        os.system(
            f"samtools cat -o {os.path.join(bampath, 'merged.bam')} {os.path.join(bampath ,'subset_*.bam')}"
        )
        # This code will look for fusions between the target gene panel and the merged bamfile assuming that the fusion has occurred between these genes.
        os.system(
            f"bedtools intersect -a {self.gene_bed} -b {os.path.join(bampath, 'merged.bam')} -wa -wb > {os.path.join(bampath ,'mappings.txt')}"
        )
        # This code will look for fusions between all possible genes in the gtf file and the merged bamfile
        os.system(
            f"bedtools intersect -a {self.all_gene_bed} -b {os.path.join(bampath, 'merged.bam')} -wa -wb > {os.path.join(bampath ,'all_mappings.txt')}"
        )
        try:
            self.fusion_candidates = pd.read_csv(
                os.path.join(bampath, "mappings.txt"), sep="\t", header=None
            )
            self.fusion_candidates = self.fusion_candidates[
                self.fusion_candidates[8].gt(50)
            ]
            self.fusion_candidates["diff"] = (
                self.fusion_candidates[6] - self.fusion_candidates[5]
            )
            self.fusion_candidates = self.fusion_candidates[
                self.fusion_candidates["diff"].gt(1000)
            ]
            uniques = self.fusion_candidates[7].duplicated(keep=False)
            doubles = self.fusion_candidates[uniques]
            counts = doubles.groupby(7)[3].transform("nunique")
            result = doubles[counts > 1]
            self.fusiontable.clear()
            with self.fusiontable:
                f = ui.input("Filter")
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
                    )
                ).bind_filter_from(f, "value").classes("w-full")
            self.fusionplot.clear()

            result, goodpairs = self._annotate_results(result)

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
                                # for gene in result_all[goodpairs].sort_values(by=7)[3].unique():
                                # if len(result[goodpairs].sort_values(by=7)[
                                #           result[goodpairs].sort_values(by=7)[3].eq(gene)]) > 2:
                                self.create_fusion_plot(
                                    gene,
                                    result[goodpairs].sort_values(by=7)[
                                        result[goodpairs].sort_values(by=7)[3].eq(gene)
                                    ],
                                )

        except Exception as e:
            print(f"{e}")

        try:
            self.fusion_candidates_all = pd.read_csv(
                os.path.join(bampath, "all_mappings.txt"), sep="\t", header=None
            )

            self.fusion_candidates_all = self.fusion_candidates_all[
                self.fusion_candidates_all[8].gt(59)
            ]

            self.fusion_candidates_all["diff"] = (
                self.fusion_candidates_all[6] - self.fusion_candidates_all[5]
            )

            self.fusion_candidates_all = self.fusion_candidates_all[
                self.fusion_candidates_all["diff"].gt(1000)
            ]

            uniques_all = self.fusion_candidates_all[7].duplicated(keep=False)
            doubles_all = self.fusion_candidates_all[uniques_all]
            counts_all = doubles_all.groupby(7)[3].transform("nunique")
            result_all = doubles_all[counts_all > 1]

            self.fusiontable_all.clear()
            with self.fusiontable_all:
                f_all = ui.input("Filter")
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
                    )
                ).bind_filter_from(f_all, "value").classes("w-full")
            self.fusionplot_all.clear()

            result_all, goodpairs = self._annotate_results(result_all)

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
                                goodpairs & result_all[goodpairs]["tag"].eq(gene_pair)
                            ][3].unique():
                                self.create_fusion_plot(
                                    gene,
                                    result_all[goodpairs].sort_values(by=7)[
                                        result_all[goodpairs]
                                        .sort_values(by=7)[3]
                                        .eq(gene)
                                    ],
                                )
                    break
        except Exception as e:
            print(f"{e}")

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


def index_page() -> None:
    # from check_connection import ConnectionDialog
    # initial_ip = "127.0.0.1"
    # my_connection = ConnectionDialog(initial_ip)
    my_connection = None
    with theme.frame("MethClass Interactive", my_connection):
        # my_connection.connect_to_minknow()
        fusion = Fusion_Panel()
        fusion.setup_ui()
        fusion.parse_bams("/Users/mattloose/datasets/nicegui_ds1305/donebams/")
        # fusion.parse_bams("/Users/mattloose/006bams")
        # CNV_PLOT.create_cnv_scatter("CNV Scatter")
        # CNV_PLOT.cnv_plotting("/Users/mattloose/datasets/ds1305_sort.hg38.h2m.bam")
        # my_object = MinknowHistograms(my_connection.positions[0])


def run_class(port: int, reload: bool):
    """
    Helper function to run the app.
    :param port: The port to serve the app on.
    :param reload: Should we reload the app on changes.
    :return:
    """
    index_page()
    ui.run(port=port, reload=reload, title="MethClass NiceGUI")


# Entrypoint for when GUI is launched by the CLI.
# e.g.: python my_app/my_cli.py
if __name__ in {"__main__", "__mp_main__"}:
    """
    Entrypoint for when GUI is launched by the CLI
    :return: None
    """
    print("GUI launched by auto-reload")

    run_class(port=12398, reload=True)
