from cnsmeth.subpages.base_analysis import BaseAnalysis
import natsort
from cnsmeth import theme, resources
import pandas as pd
import numpy as np
import os
import sys
import click
import time
from pathlib import Path
from nicegui import ui, run, background_tasks
from io import StringIO
import pysam
import asyncio
import tempfile
import shutil
import queue

from cnsmeth.subpages.SNPview import SNPview

os.environ["CI"] = "1"


def run_clair3(bamfile, bedfile, workdir, workdirout, threads, reference):
    # ToDo: handle any platform
    if sys.platform in ["darwin", "linux"]:
        runcommand = (
            f"docker run -it -v {workdir}:{workdir} "
            f"-v {workdirout}:{workdirout} "
            f"-v {reference}:{reference} "
            f"-v {reference}.fai:{reference}.fai "
            f"hkubal/clairs-to:latest "
            f"/opt/bin/run_clairs_to "
            f"--tumor_bam_fn {bamfile} "
            f"--ref_fn {reference} "
            f"--threads {threads} "
            f"--platform ont_r10_guppy_hac_5khz "
            f"--output_dir {workdirout} -b {bedfile}"
        )
        os.system(runcommand)
        shutil.copy2(f"{workdirout}/snv.vcf.gz", f"{workdirout}/output_done.vcf.gz")
        shutil.copy2(
            f"{workdirout}/indel.vcf.gz", f"{workdirout}/output_indel_done.vcf.gz"
        )

        # os.system(f"ls {workdirout}")
        command = f"snpEff hg38 {workdirout}/output_done.vcf.gz > {workdirout}/snpeff_output.vcf"
        os.system(command)
        # os.system(f"ls {workdirout}")
        command = f"SnpSift annotate {os.path.join(os.path.dirname(os.path.abspath(resources.__file__)),'clinvar.vcf')} {workdirout}/snpeff_output.vcf > {workdirout}/snpsift_output.vcf"
        os.system(command)
        command = f"snpEff hg38 {workdirout}/output_indel_done.vcf.gz > {workdirout}/snpeff_indel_output.vcf"
        os.system(command)
        # os.system(f"ls {workdirout}")
        command = f"SnpSift annotate {os.path.join(os.path.dirname(os.path.abspath(resources.__file__)), 'clinvar.vcf')} {workdirout}/snpeff_indel_output.vcf > {workdirout}/snpsift_indel_output.vcf"
        os.system(command)
        # os.system(f"ls {workdirout}")


def get_covdfs(bamfile):
    """
    This function runs modkit on a bam file and extracts the methylation data.
    """
    try:
        pysam.index(f"{bamfile}")
        newcovdf = pd.read_csv(StringIO(pysam.coverage(f"{bamfile}")), sep="\t")
        newcovdf.drop(
            columns=["coverage", "meanbaseq", "meanmapq"],
            inplace=True,
        )
        bedcovdf = pd.read_csv(
            StringIO(
                pysam.bedcov(
                    os.path.join(
                        os.path.dirname(os.path.abspath(resources.__file__)),
                        "unique_genes.bed",
                    ),
                    f"{bamfile}",
                )
            ),
            names=["chrom", "startpos", "endpos", "name", "bases"],
            sep="\t",
        )
        return newcovdf, bedcovdf
    except Exception as e:
        print(e)


def subset_bam(bamfile, targets, output):
    pysam.view("-L", f"{targets}", "-o", f"{output}", f"{bamfile}")


def sort_bam(bamfile, output, threads):
    pysam.sort(f"-@{threads}", "-o", output, bamfile)
    pysam.index(f"{output}")


def run_bedmerge(newcovdf, cov_df_main, bedcovdf, bedcov_df_main):
    merged_df = pd.merge(
        newcovdf,
        cov_df_main,
        on=["#rname", "startpos", "endpos"],
        suffixes=("_df1", "_df2"),
    )
    merged_df["numreads"] = merged_df["numreads_df1"] + merged_df["numreads_df2"]
    merged_df["covbases"] = merged_df["covbases_df1"] + merged_df["covbases_df2"]
    merged_df["meandepth"] = merged_df["meandepth_df1"] + merged_df["meandepth_df2"]

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

    merged_bed_df = pd.merge(
        bedcovdf,
        bedcov_df_main,
        on=["chrom", "startpos", "endpos", "name"],
        suffixes=("_df1", "_df2"),
    )
    merged_bed_df["bases"] = merged_bed_df["bases_df1"] + merged_bed_df["bases_df2"]
    merged_bed_df.drop(columns=["bases_df1", "bases_df2"], inplace=True)
    return merged_df, merged_bed_df


def run_bedtools(bamfile, bedfile, tempbamfile):
    """
    This function extracts the target sites from the bamfile.
    """
    try:
        os.system(f"bedtools intersect -a {bamfile} -b {bedfile} > {tempbamfile}")
        pysam.index(tempbamfile)
    except Exception as e:
        print(e)


class TargetCoverage(BaseAnalysis):
    def __init__(self, *args, target_panel=None, reference=None, **kwargs):
        self.callthreshold = 10
        self.clair3running = False
        self.targets_exceeding_threshold = 0
        self.targetbamfile = None
        self.covtable = None
        self.covtable_row_count = 0
        self.cov_df_main = pd.DataFrame()
        self.bedcov_df_main = pd.DataFrame()
        self.SNPqueue = queue.Queue()
        self.reference = reference
        if self.reference:
            self.snp_calling = True
        else:
            self.snp_calling = False
        if self.snp_calling:
            self.SNP_timer_run()
        self.target_panel = target_panel

        if self.target_panel == "rCNS2":
            self.bedfile = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)),
                "rCNS2_panel_name_uniq.bed",
            )
        elif self.target_panel == "AML":
            self.bedfile = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)),
                "AML_panel_name_uniq.bed",
            )
        super().__init__(*args, **kwargs)

    def SNP_timer_run(self):
        self.snp_timer = ui.timer(0.1, self._snp_worker)

    async def _snp_worker(self):
        """
        This function takes reads from the queue and adds them to the background thread for processing.
        """
        self.snp_timer.active = False
        if not self.SNPqueue.empty():
            while not self.SNPqueue.empty():
                gene_list, bamfile, bedfile = self.SNPqueue.get()

            workdirout = os.path.join(self.output, "clair3")
            if not os.path.exists(workdirout):
                os.mkdir(workdirout)
            # bamfile, bedfile, workdir, workdirout, threads
            self.clair3running = True
            shutil.copy2(bedfile, f"{bedfile}2")
            await background_tasks.create(
                run.cpu_bound(
                    run_clair3,
                    f"{bamfile}",
                    f"{bedfile}2",
                    self.output,
                    workdirout,
                    self.threads,
                    self.reference,
                )
            )
            if os.path.isfile(f"{workdirout}/snpsift_output.vcf"):
                self.SNPview.parse_vcf(f"{workdirout}/snpsift_output.vcf")
            if os.path.isfile(f"{workdirout}/snpsift_indel_output.vcf"):
                self.INDELview.parse_vcf(f"{workdirout}/snpsift_indel_output.vcf")

            self.clair3running = False

        else:
            await asyncio.sleep(1)
        self.snp_timer.active = True

    def setup_ui(self):
        if self.summary:
            with self.summary:
                ui.label("Current coverage estimates: Unknown")
        with ui.card().classes("w-full"):
            ui.label("Coverage Data").style(
                "color: #6E93D6; font-size: 150%; font-weight: 300"
            ).tailwind("drop-shadow", "font-bold")
            with ui.grid(columns=2).classes("w-full h-auto"):
                with ui.column():
                    with ui.card().classes("w-full"):
                        self.create_coverage_plot("Chromosome Coverage")
                with ui.column():
                    with ui.card().classes("w-full"):
                        self.create_coverage_plot_targets("Target Coverage")
        with ui.card().classes("w-full"):
            ui.label("Coverage over time").style(
                "color: #6E93D6; font-size: 150%; font-weight: 300"
            ).tailwind("drop-shadow", "font-bold")
            with ui.row().classes("w-full"):
                self.create_coverage_time_chart()
        with ui.card().classes("w-full"):
            ui.label("Coverage over targets").style(
                "color: #6E93D6; font-size: 150%; font-weight: 300"
            ).tailwind("drop-shadow", "font-bold")
            self.targ_df = ui.row().classes("w-full").style("height: 900px")
        with ui.card().classes("w-full"):
            ui.label("Candidate SNPs").style(
                "color: #6E93D6; font-size: 150%; font-weight: 300"
            ).tailwind("drop-shadow", "font-bold")
            with ui.row().classes("w-full"):
                self.SNPview = SNPview()
                self.SNPview.renderme()
            ui.label("Candidate IN/DELs").style(
                "color: #6E93D6; font-size: 150%; font-weight: 300"
            ).tailwind("drop-shadow", "font-bold")
            with ui.row().classes("w-full"):
                self.INDELview = SNPview()
                self.INDELview.renderme()
        if self.browse:
            self.show_previous_data(self.output)

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
        if not self.covtable:
            self.target_coverage_df["coverage"] = self.target_coverage_df[
                "coverage"
            ].round(2)
            with self.targ_df:
                self.targ_df.clear()
                self.covtable = (
                    ui.table.from_pandas(self.target_coverage_df, pagination=25)
                    .props("dense")
                    .classes("w-full")
                    .style("height: 900px")
                    .style("font-size: 100%; font-weight: 300")
                )
                for col in self.covtable.columns:
                    col["sortable"] = True
                    if col["name"] == "coverage":
                        col[":format"] = "value => value.toFixed(2)"
                self.covtable.add_slot(
                    "body-cell-coverage",
                    """
                    <q-td key="coverage" :props="props">
                        <q-badge :color="props.value < 10 ? 'red' : 'green'">
                            {{ props.value }}
                        </q-badge>
                    </q-td>
                """,
                )
                with self.covtable.add_slot("top-right"):
                    with ui.input(placeholder="Search").props("type=search").bind_value(
                        self.covtable, "filter"
                    ).add_slot("append"):
                        ui.icon("search")
        else:
            self.covtable.update_rows(self.target_coverage_df.to_dict(orient="records"))

    async def process_bam(self, bamfile, timestamp):
        newcovdf, bedcovdf = await background_tasks.create(
            run.cpu_bound(get_covdfs, bamfile)
        )

        tempbamfile = tempfile.NamedTemporaryFile(dir=self.output, suffix=".bam")

        await background_tasks.create(
            run.cpu_bound(run_bedtools, bamfile, self.bedfile, tempbamfile.name)
        )

        if pysam.AlignmentFile(tempbamfile.name, "rb").count(until_eof=True) > 0:
            if not self.targetbamfile:
                self.targetbamfile = os.path.join(self.output, "target.bam")
                shutil.copy2(tempbamfile.name, self.targetbamfile)
            else:
                tempbamholder = tempfile.NamedTemporaryFile(
                    dir=self.output, suffix=".bam"
                )
                pysam.cat(
                    "-o", tempbamholder.name, self.targetbamfile, tempbamfile.name
                )
                shutil.copy2(tempbamholder.name, self.targetbamfile)

        if self.cov_df_main.empty:
            self.cov_df_main = newcovdf
            self.bedcov_df_main = bedcovdf
        else:
            self.cov_df_main, self.bedcov_df_main = await background_tasks.create(
                run.cpu_bound(
                    run_bedmerge,
                    newcovdf,
                    self.cov_df_main,
                    bedcovdf,
                    self.bedcov_df_main,
                )
            )
            # self.cov_df_main, self.bedcov_df_main = run_bedmerge(newcovdf, self.cov_df_main, bedcovdf, self.bedcov_df_main)
        if self.bamqueue.empty() or self.bam_processed % 5 == 0:
            self.update_coverage_plot(self.cov_df_main)
            self.cov_df_main.to_csv(os.path.join(self.output, "coverage_main.csv"))
            # await asyncio.sleep(0.01)
            self.update_coverage_plot_targets(self.cov_df_main, self.bedcov_df_main)
            self.bedcov_df_main.to_csv(
                os.path.join(self.output, "bed_coverage_main.csv")
            )
            # await asyncio.sleep(0.01)
            self.update_coverage_time_plot(self.cov_df_main, timestamp)
            # await asyncio.sleep(0.01)
            self.target_coverage_df = self.bedcov_df_main
            self.target_coverage_df["coverage"] = (
                self.target_coverage_df["bases"] / self.target_coverage_df["length"]
            )
            self.target_coverage_df.to_csv(
                os.path.join(self.output, "target_coverage.csv")
            )
            if self.summary:
                with self.summary:
                    self.summary.clear()
                    with ui.row():
                        ui.label("Coverage Depths - ")
                        ui.label(
                            f"Global Estimated Coverage: {(self.cov_df_main['covbases'].sum()/self.cov_df_main['endpos'].sum()):.2f}x"
                        )
                        ui.label(
                            f"Targets Estimated Coverage: {(self.bedcov_df_main['bases'].sum()/self.bedcov_df_main['length'].sum()):.2f}x"
                        )
            run_list = self.target_coverage_df[
                self.target_coverage_df["coverage"].ge(self.callthreshold)
            ]
            if self.reference:
                if (
                    len(run_list) > self.targets_exceeding_threshold
                    and not self.clair3running
                ):
                    self.targets_exceeding_threshold = len(run_list)
                    run_list[["chrom", "startpos", "endpos"]].to_csv(
                        os.path.join(self.output, "targets_exceeding_threshold.bed"),
                        sep="\t",
                        header=None,
                        index=None,
                    )
                    clair3workdir = os.path.join(self.output, "clair3")
                    if not os.path.exists(clair3workdir):
                        os.mkdir(clair3workdir)
                    # await run.cpu_bound(subset_bam, self.targetbamfile, os.path.join(self.output, "targets_exceeding_threshold.bed"), os.path.join(self.output, "targets_exceeding.bam"))
                    await background_tasks.create(
                        run.cpu_bound(
                            sort_bam,
                            self.targetbamfile,
                            os.path.join(clair3workdir, "sorted_targets_exceeding.bam"),
                            self.threads,
                        )
                    )
                    if self.snp_calling:
                        self.SNPqueue.put(
                            [
                                run_list,
                                os.path.join(
                                    clair3workdir, "sorted_targets_exceeding.bam"
                                ),
                                os.path.join(
                                    self.output, "targets_exceeding_threshold.bed"
                                ),
                            ]
                        )

            self.update_target_coverage_table()

        await asyncio.sleep(0.05)
        self.running = False

    def show_previous_data(self, watchfolder):
        self.cov_df_main = pd.read_csv(os.path.join(watchfolder, "coverage_main.csv"))
        self.update_coverage_plot(self.cov_df_main)
        self.bedcov_df_main = pd.read_csv(
            os.path.join(watchfolder, "bed_coverage_main.csv")
        )
        self.update_coverage_plot_targets(self.cov_df_main, self.bedcov_df_main)
        self.target_coverage_df = pd.read_csv(
            os.path.join(watchfolder, "target_coverage.csv")
        )
        self.update_target_coverage_table()
        # if file.endswith("coverage_time_chart.csv"):
        #    self.coverage_time_chart = pd.read_csv(os.path.join(watchfolder, file))
        #    self.update_coverage_time_plot(self.cov_df_main, None)
        # self.running = False
        if os.path.isfile(f"{watchfolder}/clair3/snpsift_output.vcf"):
            self.SNPview.parse_vcf(f"{watchfolder}/clair3/snpsift_output.vcf")
        if os.path.isfile(f"{watchfolder}/clair3/snpsift_indel_output.vcf"):
            self.INDELview.parse_vcf(f"{watchfolder}/clair3/snpsift_indel_output.vcf")


def test_me(
    port: int,
    threads: int,
    watchfolder: str,
    output: str,
    reload: bool = False,
    browse: bool = False,
):
    my_connection = None
    with theme.frame("Target Coverage Data", my_connection):
        TestObject = TargetCoverage(threads, output, progress=True)
        # TestObject = MGMT_Object(threads, output, progress=True)
    if not browse:
        path = watchfolder
        searchdirectory = os.fsencode(path)
        for root, d_names, f_names in os.walk(searchdirectory):
            directory = os.fsdecode(root)
            for f in f_names:
                filename = os.fsdecode(f)
                if filename.endswith(".bam"):
                    TestObject.add_bam(os.path.join(directory, filename))
                    # break
    else:
        TestObject.progress_trackers.visible = False
        TestObject.show_previous_data(output)
    ui.run(port=port, reload=False)


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
def run_main(port, threads, watchfolder, output, browse):
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
    run_main()
