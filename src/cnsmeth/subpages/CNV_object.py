from cnv_from_bam import iterate_bam_file
from cnsmeth.subpages.base_analysis import BaseAnalysis
import natsort
from cnsmeth import theme, resources
import pandas as pd
import logging
import numpy as np
import os
import asyncio
from nicegui import ui
import click
from pathlib import Path

os.environ["CI"] = "1"


class Result:
    def __init__(self, cnv_dict):
        self.cnv = cnv_dict


def iterate_bam(bamfile, _threads=1, mapq_filter=60, copy_numbers=None):
    result = iterate_bam_file(
        bamfile,
        _threads=_threads,
        mapq_filter=mapq_filter,
        copy_numbers=copy_numbers,
        log_level=int(logging.ERROR),
    )
    return result, copy_numbers


def reduce_list(lst, max_length=500):
    while len(lst) > max_length:
        lst = lst[::2]
    return lst


class CNVAnalysis(BaseAnalysis):
    def __init__(self, *args, **kwargs):
        self.file_list = []
        self.cnv_dict = {}
        self.cnv_dict["bin_width"] = 0
        self.cnv_dict["variance"] = 0
        self.update_cnv_dict = {}
        self.result = None
        self.gene_bed = pd.read_table(
            os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)), "unique_genes.bed"
            ),
            names=["chrom", "start_pos", "end_pos", "gene"],
            header=None,
            delim_whitespace=True,
        )
        super().__init__(*args, **kwargs)

    async def process_bam(self, bamfile, timestamp):
        self.file_list.append(bamfile)
        # cnv_dict = self.update_cnv_dict.copy()
        # self.result, self.update_cnv_dict = await run.cpu_bound(iterate_bam, bamfile, _threads=self.threads, mapq_filter=60, copy_numbers=cnv_dict)
        #print (f"Processing {bamfile}, {timestamp}")
        self.result = iterate_bam_file(
            bamfile,
            _threads=self.threads,
            mapq_filter=60,
            copy_numbers=self.update_cnv_dict,
            log_level=int(logging.ERROR),
        )

        self.cnv_dict["bin_width"] = self.result.bin_width
        self.cnv_dict["variance"] = self.result.variance
        if self.summary:
            with self.summary:
                self.summary.clear()
                with ui.row():
                    ui.label(f"Current Bin Width: {self.result.bin_width}")
                    ui.label(f"Current Variance: {round(self.result.variance,3)}")
        np.save(os.path.join(self.output, 'CNV.npy'), self.result.cnv)
        np.save(os.path.join(self.output, 'CNV_dict.npy'), self.cnv_dict)

        # Only update the plot if the queue is empty?
        if self.bamqueue.empty() or self.bam_processed % 5 == 0:
            self._update_cnv_plot()
        else:
            await asyncio.sleep(0.05)
        self.running = False

    def setup_ui(self):
        self.display_row = ui.row()
        if self.summary:
            with self.summary:
                ui.label("No CNV data available.")
        with self.display_row:
            # self.progrock.visible = False
            ui.label("Copy Number Variation").style('color: #6E93D6; font-size: 150%; font-weight: 300').tailwind("drop-shadow", "font-bold")
        with ui.row():
            self.chrom_select = ui.select(
                options={"All": "All"},
                on_change=self._update_cnv_plot,
                label="Select Chromosome",
                value="All",
            ).style("width: 150px")
            self.gene_select = ui.select(
                options={"All": "All"},
                on_change=lambda e: (
                    self._update_cnv_plot()
                    if e.value == "All"
                    else self._update_cnv_plot(gene_target=e.value)
                ),
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
                    "animation": False,
                    "grid": {"containLabel": True},
                    "title": {"text": "CNV Scatter Plot"},
                    "toolbox": {"show": True, "feature": {"saveAsImage": {}}},
                    "xAxis": {
                        "type": "value",
                        "max": "dataMax",
                        "splitLine": {"show": False},
                    },
                    "yAxis": {
                        "type": "value",
                    },
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

    def _update_cnv_plot(self, gene_target=None):
        if self.result:
            total = 0
            valueslist = {"All": "All"}
            genevalueslist = {"All": "All"}
            self.chrom_filter = self.chrom_select.value

            min = 0
            max = "dataMax"

            if gene_target:
                start_pos = self.gene_bed.iloc[int(gene_target)].start_pos
                end_pos = self.gene_bed.iloc[int(gene_target)].end_pos
                chrom = self.gene_bed.iloc[int(gene_target)].chrom
                for counter, contig in enumerate(
                    natsort.natsorted(self.result.cnv), start=1
                ):
                    valueslist[counter] = contig
                    if contig == chrom:
                        break
                self.chrom_filter = counter
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

                    data = reduce_list(data)

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
                    if contig in ["chr7","chr10"]:
                        #print (self.scatter_echart.options["series"][-1])
                        pass
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

            self.chrom_select.set_options(valueslist)
            self.gene_select.set_options(genevalueslist)
            self.scatter_echart.update()

    def show_previous_data(self, output):
        result = np.load(os.path.join(output, 'CNV.npy'), allow_pickle='TRUE').item()
        self.result = Result(result)
        cnv_dict = np.load(os.path.join(output, 'CNV_dict.npy'), allow_pickle=True).item()
        self.cnv_dict["bin_width"] = cnv_dict["bin_width"]
        self.cnv_dict["variance"] = cnv_dict["variance"]
        self._update_cnv_plot()


def test_me(port: int, threads: int, watchfolder: str, output:str, reload: bool = False, browse: bool = False):
    my_connection = None
    with theme.frame("Copy Number Variation Testing.", my_connection):
        TestObject = CNVAnalysis(threads, output, progress=True)
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
        #print("Browse mode not implemented.")
        TestObject.progress_trackers.visible=False
        TestObject.show_previous_data(output)
    ui.run(port=port,reload=reload)

@click.command()
@click.option(
    "--port",
    default=12345,
    help="Port for GUI",
)
@click.option(
    "--threads",
    default=4,
    help="Number of threads available."
)
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
            #simtime=simtime,
            watchfolder=None,
            output=watchfolder,
            #sequencing_summary=sequencing_summary,
            #showerrors=showerrors,
            browse=browse,
            #exclude=exclude,
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
            #simtime=simtime,
            watchfolder=watchfolder,
            output=output,
            #sequencing_summary=sequencing_summary,
            #showerrors=showerrors,
            browse=browse,
            #exclude=exclude,
        )


if __name__ in {"__main__", "__mp_main__"}:
    print("GUI launched by auto-reload function.")
    main()
