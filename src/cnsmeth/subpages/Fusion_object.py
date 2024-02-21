from cnsmeth.subpages.base_analysis import BaseAnalysis
import os
import gff3_parser
import tempfile
import random
import asyncio
import pandas as pd


from nicegui import ui
from cnsmeth import theme, resources
from dna_features_viewer import GraphicFeature, GraphicRecord

import matplotlib

matplotlib.use("agg")
from matplotlib import pyplot as plt


os.environ["CI"] = "1"
STRAND = {"+": 1, "-": -1}


class Fusion_object(BaseAnalysis):
    def __init__(self, *args, **kwargs):
        self.gene_bed = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)), "unique_genes.bed"
        )
        self.all_gene_bed = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)), "all_genes2.bed"
        )
        self.fusion_candidates = pd.DataFrame()
        self.fusion_candidates_all = pd.DataFrame()
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
        super().__init__(*args, **kwargs)

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
                        self.fusiontable = ui.row().classes("w-full")
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
                        self.fusiontable_all = ui.row().classes("w-full")
                        with self.fusiontable_all:
                            ui.label("Table not yet available.")

    def fusion_table_all(self):
        uniques_all = self.fusion_candidates_all[7].duplicated(keep=False)
        doubles_all = self.fusion_candidates_all[uniques_all]
        counts_all = doubles_all.groupby(7)[3].transform("nunique")
        result_all = doubles_all[counts_all > 1]
        result_all.to_csv(os.path.join(self.output, "fusion_candidates_all.csv"))
        if not result_all.empty:
            self.fusiontable_all.clear()
            with self.fusiontable_all:
                ui.aggrid.from_pandas(
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
                    theme="material",
                    options={
                        "defaultColDef": {
                            "flex": 1,
                            "minWidth": 150,
                            "sortable": True,
                            "resizable": True,
                        },
                        "columnDefs": [
                            {
                                "headerName": "Chromosome",
                                "field": "chromBED",
                                "filter": "agTextColumnFilter",
                                "floatingFilter": False,
                            },
                            {
                                "headerName": "BS",
                                "field": "BS",
                                "filter": "agNumberColumnFilter",
                                "floatingFilter": False,
                            },
                            {
                                "headerName": "BE",
                                "field": "BE",
                                "filter": "agNumberColumnFilter",
                                "floatingFilter": False,
                            },
                            {
                                "headerName": "Gene",
                                "field": "Gene",
                                "filter": "agTextColumnFilter",
                                "floatingFilter": False,
                            },
                            {
                                "headerName": "chrom",
                                "field": "chrom",
                                "filter": "agTextColumnFilter",
                                "floatingFilter": False,
                            },
                            {
                                "headerName": "mS",
                                "field": "mS",
                                "filter": "agNumberColumnFilter",
                                "floatingFilter": False,
                            },
                            {
                                "headerName": "mE",
                                "field": "mE",
                                "filter": "agNumberColumnFilter",
                                "floatingFilter": False,
                            },
                            {
                                "headerName": "readID",
                                "field": "readID",
                                "filter": "agTextColumnFilter",
                                "floatingFilter": False,
                            },
                            {
                                "headerName": "mapQ",
                                "field": "mapQ",
                                "filter": "agNumberColumnFilter",
                                "floatingFilter": False,
                            },
                            {
                                "headerName": "strand",
                                "field": "strand",
                                "filter": "agTextColumnFilter",
                                "floatingFilter": False,
                            },
                        ],
                        "pagination": True,
                        "paginationAutoPageSize": True,
                    },
                    auto_size_columns=True,
                ).classes("max-h-100 min-w-full").style("height: 900px")
            self.fusionplot_all.clear()

        result_all, goodpairs = self._annotate_results(result_all)

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
                                goodpairs & result_all[goodpairs]["tag"].eq(gene_pair)
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
                                    result_all[goodpairs].sort_values(by=7)[3].eq(gene)
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
                                                        strand=STRAND[row["Strand"]],
                                                        thickness=4,
                                                        color="#ffd700",
                                                    )
                                                )
                                                first_index = int(row["Start"]) - 1000
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
        uniques = self.fusion_candidates[7].duplicated(keep=False)
        doubles = self.fusion_candidates[uniques]
        counts = doubles.groupby(7)[3].transform("nunique")
        result = doubles[counts > 1]
        result.to_csv(os.path.join(self.output, "fusion_candidates_master.csv"))
        self.fusiontable.clear()
        with self.fusiontable:
            ui.aggrid.from_pandas(
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
                theme="material",
                options={
                    "defaultColDef": {
                        "sortable": True,
                        "resizable": True,
                    },
                    "columnDefs": [
                        {
                            "headerName": "Chromosome",
                            "field": "chromBED",
                            "filter": "agTextColumnFilter",
                            "floatingFilter": False,
                        },
                        {
                            "headerName": "BS",
                            "field": "BS",
                            "filter": "agNumberColumnFilter",
                            "floatingFilter": False,
                        },
                        {
                            "headerName": "BE",
                            "field": "BE",
                            "filter": "agNumberColumnFilter",
                            "floatingFilter": False,
                        },
                        {
                            "headerName": "Gene",
                            "field": "Gene",
                            "filter": "agTextColumnFilter",
                            "floatingFilter": False,
                        },
                        {
                            "headerName": "chrom",
                            "field": "chrom",
                            "filter": "agTextColumnFilter",
                            "floatingFilter": False,
                        },
                        {
                            "headerName": "mS",
                            "field": "mS",
                            "filter": "agNumberColumnFilter",
                            "floatingFilter": False,
                        },
                        {
                            "headerName": "mE",
                            "field": "mE",
                            "filter": "agNumberColumnFilter",
                            "floatingFilter": False,
                        },
                        {
                            "headerName": "readID",
                            "field": "readID",
                            "filter": "agTextColumnFilter",
                            "floatingFilter": False,
                        },
                        {
                            "headerName": "mapQ",
                            "field": "mapQ",
                            "filter": "agNumberColumnFilter",
                            "floatingFilter": False,
                        },
                        {
                            "headerName": "strand",
                            "field": "strand",
                            "filter": "agTextColumnFilter",
                            "floatingFilter": False,
                        },
                    ],
                    "pagination": True,
                    "paginationAutoPageSize": True,
                },
                auto_size_columns=False,
            ).classes("max-h-100 min-w-full").style("height: 900px")
        self.fusionplot.clear()

        result, goodpairs = self._annotate_results(result)

        with self.fusionplot.classes("w-full"):
            for gene_pair in result[goodpairs].sort_values(by=7)["tag"].unique():
                with ui.card():
                    with ui.row():
                        ui.label(f"{gene_pair}").tailwind("drop-shadow", "font-bold")
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

        tempreadfile = tempfile.NamedTemporaryFile(suffix=".txt")
        tempbamfile = tempfile.NamedTemporaryFile(suffix=".bam")
        tempmappings = tempfile.NamedTemporaryFile(suffix=".txt")
        tempallmappings = tempfile.NamedTemporaryFile(suffix=".txt")

        # Get the subset of reads that map to the target gene panel and have supplementary alignments
        os.system(
            f"samtools view -L {self.gene_bed} -d SA {bamfile} | cut -f1 > {tempreadfile.name}"
        )
        if os.path.getsize(tempreadfile.name) > 0:
            os.system(
                f"samtools view -N {tempreadfile.name} -o {tempbamfile.name} {bamfile}"
            )
            if os.path.getsize(tempbamfile.name) > 0:
                os.system(
                    f"bedtools intersect -a {self.gene_bed} -b {tempbamfile.name} -wa -wb > {tempmappings.name}"
                )
                os.system(
                    f"bedtools intersect -a {self.all_gene_bed} -b {tempbamfile.name} -wa -wb > {tempallmappings.name}"
                )
                if os.path.getsize(tempmappings.name) > 0:
                    fusion_candidates = pd.read_csv(
                        tempmappings.name, sep="\t", header=None
                    )
                    # Filter to only include good mappings
                    fusion_candidates = fusion_candidates[fusion_candidates[8].gt(50)]
                    # Filter to only keep reads that map to 1 kb or more of the reference.
                    fusion_candidates["diff"] = (
                        fusion_candidates[6] - fusion_candidates[5]
                    )
                    fusion_candidates = fusion_candidates[
                        fusion_candidates["diff"].gt(1000)
                    ]
                    if self.fusion_candidates.empty:
                        self.fusion_candidates = fusion_candidates
                    else:
                        self.fusion_candidates = pd.concat(
                            [self.fusion_candidates, fusion_candidates]
                        ).reset_index(drop=True)

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

                    if self.fusion_candidates_all.empty:
                        self.fusion_candidates_all = fusion_candidates_all
                    else:
                        self.fusion_candidates_all = pd.concat(
                            [self.fusion_candidates_all, fusion_candidates_all]
                        ).reset_index(drop=True)

                try:
                    self.fusion_table()
                except Exception as e:
                    print(f"Error 522: {e}")
                try:
                    self.fusion_table_all()
                except Exception as e:
                    print(f"Error 526: {e}")
        await asyncio.sleep(0.1)
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


def test_me():
    my_connection = None
    with theme.frame("Copy Number Variation Interactive", my_connection):
        TestObject = Fusion_object(progress=True)
        # path = "tests/static/bam"
        path = "/users/mattloose/datasets/ds1305_Intraop0006_A/20231123_1233_P2S-00770-A_PAS59057_b1e841e7/bam_pass"
        directory = os.fsencode(path)
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if filename.endswith(".bam"):
                TestObject.add_bam(os.path.join(path, filename))
    ui.run(port=12345)


# Entrypoint for when GUI is launched by the CLI.
# e.g.: python my_app/my_cli.py
if __name__ in {"__main__", "__mp_main__"}:
    """
    Entrypoint for when GUI is launched by the CLI
    :return: None
    """
    print("GUI launched by auto-reload")
    test_me()
