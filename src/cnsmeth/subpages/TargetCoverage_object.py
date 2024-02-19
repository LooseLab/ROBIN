from cnsmeth.subpages.base_analysis import BaseAnalysis
import natsort
from cnsmeth import theme, resources
import pandas as pd
import numpy as np
import os
import time
from nicegui import ui
from io import StringIO
import pysam
os.environ["CI"] = "1"


class TargetCoverage(BaseAnalysis):
    def __init__(self, *args, **kwargs):
        self.cov_df_main = pd.DataFrame()
        self.bedcov_df_main = pd.DataFrame()
        super().__init__(*args, **kwargs)

    def setup_ui(self):
        with ui.card().style("width: 100%"):
            ui.label("Coverage Data").tailwind("drop-shadow", "font-bold")
            with ui.grid(columns=2).classes("w-full h-auto"):
                with ui.column():
                    with ui.card().style("width: 100%"):
                        self.create_coverage_plot("Chromosome Coverage")
                with ui.column():
                    with ui.card().style("width: 100%"):
                        self.create_coverage_plot_targets("Target Coverage")
        with ui.card().style("width: 100%"):
            self.create_coverage_time_chart()
        with ui.card().style("width: 100%"):
            ui.label("Coverage over targets").tailwind("drop-shadow", "font-bold")
            self.targ_df = ui.row().classes("w-full").style("height: 900px")

    def create_coverage_plot(self, title):
        self.echart3 = (
            ui.echart(
                {
                    "grid": {"containLabel": True},
                    "title": {"text": title},
                    "toolbox": {"show": True, "feature": {"saveAsImage": {}}},
                    "yAxis": {"type": "value"},
                    "xAxis": {
                        "type": "category",
                        "data": [],
                        "axisTick": {"alignWithLabel": True},
                        "axisLabel": {
                            "interval": 0,
                            "rotate": 30,
                        },
                        "inverse": False,
                    },
                    #'legend': {},
                    "series": [],
                }
            )
            .style("height: 350px")
            .classes("border-double")
        )

    def create_coverage_plot_targets(self, title):
        self.echart4 = (
            ui.echart(
                {
                    "grid": {"containLabel": True},
                    "title": {"text": title},
                    "toolbox": {"show": True, "feature": {"saveAsImage": {}}},
                    "yAxis": {"type": "value"},
                    "xAxis": {
                        "type": "category",
                        "data": [],
                        "axisTick": {"alignWithLabel": True},
                        "axisLabel": {
                            "interval": 0,
                            "rotate": 30,
                        },
                        "inverse": False,
                    },
                    "legend": {"data": ["Off Target", "On Target"]},
                    "series": [],
                }
            )
            .style("height: 350px")
            .classes("border-double")
        )

    def create_coverage_time_chart(self):
        self.coverage_time_chart = (
            ui.echart(
                {
                    "grid": {"containLabel": True},
                    "title": {"text": "Coverage Over Time"},
                    "toolbox": {"show": True, "feature": {"saveAsImage": {}}},
                    "xAxis": {"type": "time"},
                    "yAxis": {"type": "value", "data": [], "inverse": False},
                    #'tooltip': {
                    #    'order': 'valueDesc',
                    #    'trigger': 'axis'
                    # },
                    "series": [
                        {
                            "type": "line",
                            "smooth": True,
                            "name": "Coverage",
                            "emphasis": {"focus": "series"},
                            "endLabel": {
                                "show": True,
                                "formatter": "{a}",
                                "distance": 20,
                            },
                            "lineStyle": {
                                "width": 2,
                            },
                            "data": [],
                        }
                    ],
                }
            )
            .style("height: 350px")
            .classes("border-double")
        )

    def update_coverage_plot(self, covdf):
        """
        Replaces the data in the RapidCNS2 plot.
        :param covdf: a pandas dataframe of
        :return:
        """
        sorteddf = covdf.sort_values(
            by="#rname",
            key=lambda x: np.argsort(natsort.index_natsorted(covdf["#rname"])),
        )
        sorteddf = sorteddf[sorteddf["#rname"] != "chrM"]
        self.echart3.options["title"]["text"] = "Per Chromosome Coverage"
        self.echart3.options["xAxis"]["data"] = sorteddf["#rname"].to_list()
        self.echart3.options["series"] = [
            {
                "type": "bar",
                "name": "Chromosome",
                "barWidth": "60%",
                "data": sorteddf["meandepth"].to_list(),
            }
        ]
        self.echart3.update()

    def update_coverage_plot_targets(self, covdf, bedcovdf):
        """
        Replaces the data in the RapidCNS2 plot.
        :param covdf: a pandas dataframe of
        :param bedcovdf: a pandas dataframe of
        :return:
        """
        bedcovdf["length"] = bedcovdf["endpos"] - bedcovdf["startpos"] + 1
        grouped = (
            bedcovdf.groupby("chrom")
            .agg({"bases": "sum", "length": "sum"})
            .reset_index()
        )
        groupeddf = grouped.sort_values(
            by="chrom",
            key=lambda x: np.argsort(natsort.index_natsorted(grouped["chrom"])),
        )
        groupeddf = groupeddf[groupeddf["chrom"] != "chrM"]
        groupeddf["meandepth"] = groupeddf["bases"] / groupeddf["length"]
        sorteddf = covdf.sort_values(
            by="#rname",
            key=lambda x: np.argsort(natsort.index_natsorted(covdf["#rname"])),
        )
        sorteddf = sorteddf[sorteddf["#rname"] != "chrM"]
        self.echart4.options["title"]["text"] = "Per Chromosome Target Coverage"
        self.echart4.options["xAxis"]["data"] = sorteddf["#rname"].to_list()
        self.echart4.options["series"] = [
            {
                "type": "scatter",
                "name": "Off Target",
                "symbolSize": 10,
                "data": sorteddf["meandepth"].to_list(),
            },
            {
                "type": "scatter",
                "name": "On Target",
                "symbolSize": 10,
                "data": groupeddf["meandepth"].to_list(),
            },
        ]
        self.echart4.update()

    def update_coverage_time_plot(self, covdf, timestamp):
        """
        Replaces the data in the RapidCNS2 plot.
        :param x: list of tumour types
        :param y: confidence scores for each tumour type
        :param count: the number of bams used to generate the plot
        :return:
        """
        bases = covdf["covbases"].sum()
        genome = covdf["endpos"].sum()
        coverage = bases / genome
        if timestamp:
            currenttime = timestamp * 1000
        else:
            currenttime = time.time() * 1000
        self.coverage_time_chart.options["series"][0]["data"].append(
            [currenttime, coverage]
        )
        self.coverage_time_chart.update()

    def update_target_coverage_table(self):
        with self.targ_df:
            self.targ_df.clear()
            ui.aggrid.from_pandas(
                self.target_coverage_df,
                theme="material",
                options={
                    "defaultColDef": {
                        "sortable": True,
                        "resizable": True,
                    },
                    "columnDefs": [
                        {
                            "headerName": "Chromosome",
                            "field": "chrom",
                            "filter": "agTextColumnFilter",
                            "floatingFilter": False,
                        },
                        {
                            "headerName": "Start",
                            "field": "startpos",
                            "filter": "agNumberColumnFilter",
                            "floatingFilter": False,
                        },
                        {
                            "headerName": "End",
                            "field": "endpos",
                            "filter": "agNumberColumnFilter",
                            "floatingFilter": False,
                        },
                        {
                            "headerName": "Gene/s",
                            "field": "name",
                            "filter": "agTextColumnFilter",
                            "floatingFilter": False,
                        },
                        {
                            "headerName": "Bases",
                            "field": "bases",
                            "filter": "agNumberColumnFilter",
                            "floatingFilter": False,
                        },
                        {
                            "headerName": "Length",
                            "field": "length",
                            "filter": "agNumberColumnFilter",
                            "floatingFilter": False,
                        },
                        {
                            "headerName": "Coverage",
                            "field": "coverage",
                            "filter": "agNumberColumnFilter",
                            "floatingFilter": False,
                        },
                    ],
                    "pagination": True,
                    "paginationAutoPageSize": True,
                },

            ).classes("w-full").style("height: 900px")


    def process_bam(self, bamfile, timestamp):
        pysam.index(f"{bamfile}")
        newcovdf = pd.read_csv(
            StringIO(pysam.coverage(f"{bamfile}")), sep="\t"
        )
        newcovdf.drop(
            columns=["coverage", "meanbaseq", "meanmapq"],
            inplace=True,
        )
        bedcovdf = pd.read_csv(
            StringIO(
                pysam.bedcov(
                    os.path.join(
                        os.path.dirname(
                            os.path.abspath(resources.__file__)
                        ),
                        "unique_genes.bed",
                    ),
                    f"{bamfile}",
                )
            ),
            names=["chrom", "startpos", "endpos", "name", "bases"],
            sep="\t",
        )

        if self.cov_df_main.empty:
            self.cov_df_main = newcovdf
            self.bedcov_df_main = bedcovdf
        else:
            merged_df = pd.merge(
                newcovdf,
                self.cov_df_main,
                on=["#rname", "startpos", "endpos"],
                suffixes=("_df1", "_df2"),
            )
            merged_df["numreads"] = (
                    merged_df["numreads_df1"] + merged_df["numreads_df2"]
            )
            merged_df["covbases"] = (
                    merged_df["covbases_df1"] + merged_df["covbases_df2"]
            )
            merged_df["meandepth"] = (
                    merged_df["meandepth_df1"] + merged_df["meandepth_df2"]
            )

            merged_df.drop(
                columns=[
                    "numreads_df1",
                    "numreads_df2",
                    "meandepth_df1",
                    "meandepth_df2",
                    "covbases_df1",
                    "covbases_df2",
                ],
                inplace=True,
            )
            self.cov_df_main = merged_df

            merged_bed_df = pd.merge(
                bedcovdf,
                self.bedcov_df_main,
                on=["chrom", "startpos", "endpos", "name"],
                suffixes=("_df1", "_df2"),
            )
            merged_bed_df["bases"] = merged_bed_df["bases_df1"] + merged_bed_df["bases_df2"]
            merged_bed_df.drop(columns=["bases_df1", "bases_df2"], inplace=True)
            self.bedcov_df_main = merged_bed_df

        if self.queue.empty() or self.bam_processed % 25 == 0:
            self.update_coverage_plot(self.cov_df_main)
            self.update_coverage_plot_targets(self.cov_df_main, self.bedcov_df_main)
            self.update_coverage_time_plot(self.cov_df_main, timestamp)

            self.target_coverage_df = self.bedcov_df_main
            self.target_coverage_df["coverage"] = (
                    self.target_coverage_df["bases"] / self.target_coverage_df["length"]
            )
            self.update_target_coverage_table()
        self.running = False



def test_me():
    my_connection = None
    with theme.frame("Copy Number Variation Interactive", my_connection):
        TestObject = TargetCoverage(progress=True)
        #path = "tests/static/bam"
        path = "/users/mattloose/datasets/ds1305_Intraop0006_A/20231123_1233_P2S-00770-A_PAS59057_b1e841e7/bam_pass"
        directory = os.fsencode(path)
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if filename.endswith(".bam"):
                TestObject.add_bam(os.path.join(path, filename))


# Entrypoint for when GUI is launched by the CLI.
# e.g.: python my_app/my_cli.py
if __name__ in {"__main__", "__mp_main__"}:
    """
    Entrypoint for when GUI is launched by the CLI
    :return: None
    """
    print("GUI launched by auto-reload")
    test_me()
    ui.run(port=12345)