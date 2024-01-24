# Python imports.
from __future__ import annotations


from nicegui import Tailwind, ui, app, run
import threading

import time
import datetime

from methnicegui import theme
from cnv_from_bam import iterate_bam_file
import os

os.environ["CI"] = "1"

import natsort
import numpy as np
import pandas as pd

from methnicegui import resources

import queue


class CNV_Plot:
    def __init__(self, **kwargs):
        print("CNV_Plot Initialised")
        self.gene_bed = pd.read_table(
            os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)), "unique_genes.bed"
            ),
            names=["chrom", "start_pos", "end_pos", "gene"],
            header=None,
            delim_whitespace=True,
        )
        self.cnv_dict = {}
        self.cnv_dict["bin_width"] = 0
        self.cnv_dict["variance"] = 0
        self.display_row = None
        self.cnv_queue = queue.Queue()
        self.worker = threading.Thread(target=self._cnv_plotting, args=())
        self.worker.daemon = True
        self.worker.start()

    def create_cnv_scatter(self, title):
        self.display_row = ui.row()
        with self.display_row:
            ui.label("Copy Number Variation").tailwind("drop-shadow", "font-bold")
        with ui.row():
            self.chrom_select = ui.select(
                options={"All": "All"},
                on_change=self._update_cnv_plot,
                label="Select Chromosome",
                value="All",
            ).style("width: 150px")
            self.gene_select = ui.select(
                options={"All": "All"},
                on_change=lambda e: self._update_cnv_plot()
                if e.value == "All"
                else self._update_cnv_plot(gene_target=e.value),
                label="Select Gene",
                value="All",
            ).style("width: 150px")
            ui.label().bind_text_from(
                self.cnv_dict, "bin_width", backward=lambda n: f"Bin Width: {n}"
            )
            ui.label().bind_text_from(
                self.cnv_dict, "variance", backward=lambda n: f"Variance: {round(n,3)}"
            )
        self.scatter_echart = (
            ui.echart(
                {
                    "grid": {"containLabel": True},
                    "title": {"text": title},
                    "toolbox": {"show": True, "feature": {"saveAsImage": {}}},
                    "xAxis": {
                        "type": "value",
                        "max": "dataMax",
                        "splitLine": {"show": False},
                    },
                    #'yAxis': {'axisLabel': {':formatter': 'value => "Ploidy" + value'}},
                    "yAxis": {
                        "type": "value",
                    },
                    #'legend': {'type':'scroll','formatter': '{name}','top': 'bottom'},
                    "dataZoom": [
                        {"type": "slider", "xAxisIndex": 0, "filterMode": "none"},
                        {
                            "type": "slider",
                            "yAxisIndex": 0,
                            "filterMode": "none",
                            "startValue": 0,
                            "endValue": 8,
                        },
                        {"type": "inside", "xAxisIndex": 0, "filterMode": "none"},
                        {"type": "inside", "yAxisIndex": 0, "filterMode": "none"},
                    ],
                    "series": [
                        {
                            "type": "scatter",
                            "symbolSize": 1,
                            "data": [],
                        }
                    ],
                }
            )
            .style("height: 450px")
            .classes("border-double")
        )

    def cnv_plotting(self, bam_path):
        self.cnv_queue.put(bam_path)

    def _cnv_plotting(self):
        while True:
            if self.display_row:
                if not self.cnv_queue.empty():
                    bam_path = self.cnv_queue.get()
                    self.display_row.clear()
                    with self.display_row:
                        ui.label("Copy Number Variation").tailwind(
                            "drop-shadow", "font-bold"
                        )
                        ui.label(
                            f"Last update: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
                        )
                        ui.label(f"{bam_path}")
                    self.result = iterate_bam_file(
                        bam_path,
                        _threads=8,
                        mapq_filter=60,  # , log_level=logging.getLevelName("WARN")
                    )
                    print(self.result)
                    self.cnv_dict["bin_width"] = self.result.bin_width
                    self.cnv_dict["variance"] = self.result.variance
                    # print(self.result)
                    self._update_cnv_plot()
            time.sleep(5)

    def _update_cnv_plot(self, gene_target=None):
        if self.result:
            total = 0
            valueslist = {"All": "All"}
            genevalueslist = {"All": "All"}
            self.chrom_filter = self.chrom_select.value

            min = 0
            max = "dataMax"

            if gene_target:
                print("Gene Target")
                start_pos = self.gene_bed.iloc[int(gene_target)].start_pos
                end_pos = self.gene_bed.iloc[int(gene_target)].end_pos
                chrom = self.gene_bed.iloc[int(gene_target)].chrom
                print(chrom)
                for counter, contig in enumerate(
                    natsort.natsorted(self.result.cnv), start=1
                ):
                    valueslist[counter] = contig
                    if contig == chrom:
                        break
                self.chrom_filter = counter
                print(self.chrom_filter)
                min = start_pos - 10 * self.cnv_dict["bin_width"]
                max = end_pos + 10 * self.cnv_dict["bin_width"]
            if self.chrom_filter == "All":
                counter = 0

                self.scatter_echart.options["title"][
                    "text"
                ] = "Copy Number Variation - All Chromosomes"
                self.scatter_echart.options["series"] = []
                for contig, cnv in natsort.natsorted(self.result.cnv.items()):
                    if contig == "chrM":
                        continue
                    counter += 1
                    valueslist[counter] = contig

                    data = list(
                        zip(
                            (np.arange(len(cnv)) + total) * self.cnv_dict["bin_width"],
                            cnv,
                        )
                    )

                    total += len(cnv)
                    self.scatter_echart.options["xAxis"]["max"] = max
                    self.scatter_echart.options["xAxis"]["min"] = min
                    self.scatter_echart.options["series"].append(
                        {
                            "type": "scatter",
                            "name": contig,
                            "data": data,
                            "symbolSize": 5,
                            "markLine": {
                                "lineStyle": {"width": 1},
                                "symbol": "none",
                                "label": {"formatter": contig},
                                "data": [
                                    {
                                        "name": contig,
                                        "xAxis": (total * self.cnv_dict["bin_width"]),
                                    }
                                ],
                            },
                        }
                    )
                    for index, gene in self.gene_bed[
                        self.gene_bed["chrom"] == contig
                    ].iterrows():
                        genevalueslist[index] = f"{gene.chrom} - {gene.gene}"

            else:
                self.scatter_echart.options["series"] = []

                for counter, contig in enumerate(
                    natsort.natsorted(self.result.cnv), start=1
                ):
                    valueslist[counter] = contig

                contig, cnv = natsort.natsorted(self.result.cnv.items())[
                    int(self.chrom_filter) - 1
                ]
                # self.log(self.gene_bed[self.gene_bed['chrom']==contig])
                data = list(
                    zip((np.arange(len(cnv)) + total) * self.cnv_dict["bin_width"], cnv)
                )
                if not gene_target:
                    min = 0
                    max = "dataMax"

                else:
                    if gene_target == "All":
                        min = 0
                        max = "dataMax"
                    else:
                        start_pos = self.gene_bed.iloc[int(gene_target)].start_pos
                        end_pos = self.gene_bed.iloc[int(gene_target)].end_pos

                        chrom = self.gene_bed.iloc[int(gene_target)].chrom
                        counter = 0
                        for counter, contig in enumerate(
                            natsort.natsorted(self.result.cnv), start=1
                        ):
                            if contig == chrom:
                                self.chrom_filter = counter
                                break
                        # self.update_plot(gene_target=[start_pos, end_pos])

                        min = start_pos - 10 * self.cnv_dict["bin_width"]
                        max = end_pos + 10 * self.cnv_dict["bin_width"]

                        if start_pos - min > 2_000_000:
                            min = start_pos - 2_000_000
                        if max - end_pos > 2_000_000:
                            max = end_pos + 2_000_000

                        if min < 0:
                            min = 0

                    # self.cnv_plot.update_plot(x=np.arange(len(cnv)) * self.result.bin_width, y=cnv,
                    #                          range=[gene_target[0] - 10*self.bin_width,
                    #                                 gene_target[1] + 10*self.bin_width])
                # self.cnv_plot.add_text(contig, total, y=y)
                # y = 7
                self.scatter_echart.options["title"][
                    "text"
                ] = f"Copy Number Variation - {contig}"
                self.scatter_echart.options["xAxis"]["max"] = max
                self.scatter_echart.options["xAxis"]["min"] = min
                self.scatter_echart.options["series"].append(
                    {
                        "type": "scatter",
                        "name": contig,
                        "data": data,
                        "symbolSize": 3,
                        "markArea": {
                            "itemStyle": {"color": "rgba(255, 173, 177, 0.4)"},
                            "data": [],
                        },
                    }
                )
                for index, gene in self.gene_bed[
                    self.gene_bed["chrom"] == contig
                ].iterrows():
                    genevalueslist[index] = f"{gene.chrom} - {gene.gene}"

                for _, row in self.gene_bed[
                    self.gene_bed["chrom"] == contig
                ].iterrows():
                    self.scatter_echart.options["series"][0]["markArea"]["data"].append(
                        [
                            {
                                "name": row["gene"],
                                "xAxis": row["start_pos"],
                            },
                            {
                                "xAxis": row["end_pos"],
                            },
                        ]
                    )

                    # self.cnv_plot.add_vline(row["start_pos"], color="green+")
                    # self.cnv_plot.add_text(row["gene"], row["start_pos"], y=y)
            # print (self.scatter_echart.options)
            self.chrom_select.set_options(valueslist)
            # self.chrom_select.update()
            self.gene_select.set_options(genevalueslist)
            # self.gene_select.update()
            self.scatter_echart.update()


def index_page() -> None:
    # from check_connection import ConnectionDialog
    # initial_ip = "127.0.0.1"
    # my_connection = ConnectionDialog(initial_ip)
    my_connection = None
    with theme.frame("MethClass Interactive", my_connection):
        # my_connection.connect_to_minknow()
        ui.label(f"Hello")
        CNV_PLOT = CNV_Plot()
        CNV_PLOT.create_cnv_scatter("CNV Scatter")
        CNV_PLOT.cnv_plotting("/Users/mattloose/datasets/ds1305_sort.hg38.h2m.bam")
        # my_object = MinknowHistograms(my_connection.positions[0])


def run_class(port: int, reload: bool):
    """
    Helper function to run the app.
    :param port: The port to serve the app on.
    :param reload: Should we reload the app on changes.
    :return:
    """
    index_page()
    ui.run(
        port=port, reload=reload, title="MethClass NiceGUI"
    )  # , native=True, fullscreen=False, window_size=(1200, 1000))


def main():  # , threads, simtime, watchfolder, output, sequencing_summary):
    from check_connection import ConnectionDialog

    """
    Entrypoint for when GUI is launched directly.
    :return: None
    """
    run_class(port=12398, reload=False)


# Entrypoint for when GUI is launched by the CLI.
# e.g.: python my_app/my_cli.py
if __name__ in {"__main__", "__mp_main__"}:
    """
    Entrypoint for when GUI is launched by the CLI
    :return: None
    """
    if __name__ == "__mp_main__":
        print("GUI launched by auto-reload")

    run_class(port=12398, reload=True)
