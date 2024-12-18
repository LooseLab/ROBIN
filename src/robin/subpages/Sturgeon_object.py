from robin.subpages.base_analysis import BaseAnalysis
import os
import sys
import tempfile
import time
import shutil
import pandas as pd
from nicegui import ui, app, run
from robin import theme
import pysam
from robin import models
from sturgeon.callmapping import (
    merge_probes_methyl_calls,
    probes_methyl_calls_to_bed,
)
import click
from pathlib import Path
from typing import List, Tuple
import logging

# Use the main logger configured in the main application
logger = logging.getLogger(__name__)


def run_probes_methyl_calls(merged_output_file, bed_output_file):
    probes_methyl_calls_to_bed(merged_output_file, bed_output_file)


def run_sturgeon_merge_probes(calls_per_probe_file, merged_output_file):
    merge_probes_methyl_calls(
        [calls_per_probe_file, merged_output_file],
        merged_output_file,
    )


def pysam_cat(tempbam, tomerge):
    pysam.cat("-o", tempbam, *tomerge)


def run_modkit(file, temp, threads):
    """
    This function runs modkit on a bam file and extracts the methylation data.
    Adjusts command based on modkit version.
    """
    try:
        # Get modkit version
        import subprocess
        version_output = subprocess.check_output(['modkit', '--version'], text=True).strip()
        version = version_output.split()[-1]  # Gets '0.4.1' from 'mod_kit 0.4.1'
        
        # Parse version number
        major, minor, *_ = version.split('.')
        version_num = float(f"{major}.{minor}")
        
        # Choose appropriate command based on version
        extract_cmd = "extract full" if version_num >= 0.4 else "extract"
        
        os.system(
            f"modkit {extract_cmd} --ignore h -t {threads} {file} {temp} "
            f"--force --suppress-progress >/dev/null 2>&1"
        )
    except Exception as e:
        print(e)
        pass


def run_sturgeon_predict(bedDir, dataDir, modelfile):
    os.system(
        f"sturgeon predict -i {bedDir} -o {dataDir} "
        f"--model-files {modelfile} >/dev/null 2>&1"
    )


def run_sturgeon_inputtobed(temp, temp2):
    try:
        os.system(
            f"sturgeon inputtobed -i {temp} -o {temp2} -s modkit "
            f"--reference-genome hg38 >/dev/null 2>&1"
        )
        # self.log(temp2)
    except Exception as e:
        print(e)
        # self.log(e)
        pass


class Sturgeon_object(BaseAnalysis):
    def __init__(self, *args, **kwargs):
        self.sturgeon_df_store = {}
        self.threshold = 0.05
        self.first_run = {}
        self.modelfile = os.path.join(
            os.path.dirname(os.path.abspath(models.__file__)), "general.zip"
        )
        self.dataDir = {}
        self.bedDir = {}
        self.st_num_probes = {}
        super().__init__(*args, **kwargs)

    def setup_ui(self):
        self.card = ui.card().classes('dark:bg-black').style("width: 100%")
        with self.card:
            with ui.grid(columns=8).classes("w-full h-auto"):
                with ui.card().classes(
                    f"min-[{self.MENU_BREAKPOINT+1}px]:col-span-3 max-[{self.MENU_BREAKPOINT}px]:col-span-8 dark:bg-black"
                ):
                    self.create_sturgeon_chart("Sturgeon")
                with ui.card().classes(
                    f"min-[{self.MENU_BREAKPOINT+1}px]:col-span-5 max-[{self.MENU_BREAKPOINT}px]:col-span-8 dark:bg-black"
                ):
                    self.create_sturgeon_time_chart("Sturgeon Time Series")
        if self.summary:
            with self.summary:
                ui.label("Sturgeon classification: Unknown")
        if self.browse:
            self.show_previous_data()
        else:
            ui.timer(5, lambda: self.show_previous_data())

    def show_previous_data(self):
        if not self.browse:
            for item in app.storage.general[self.mainuuid]:
                if item == "sample_ids":
                    for sample in app.storage.general[self.mainuuid][item]:
                        self.sampleID = sample
            output = self.output
        if self.browse:
            output = self.check_and_create_folder(self.output, self.sampleID)

        if self.check_file_time(os.path.join(output, "sturgeon_scores.csv")):
            self.sturgeon_df_store = pd.read_csv(
                os.path.join(os.path.join(output, "sturgeon_scores.csv")),
                index_col=0,
            )
            columns_greater_than_threshold = (
                self.sturgeon_df_store > self.threshold
            ).any()
            columns_not_greater_than_threshold = ~columns_greater_than_threshold
            result = self.sturgeon_df_store.columns[
                columns_not_greater_than_threshold
            ].tolist()
            self.update_sturgeon_time_chart(self.sturgeon_df_store.drop(columns=result))
            lastrow = self.sturgeon_df_store.iloc[-1].drop("number_probes")
            lastrow_plot = lastrow.sort_values(ascending=False).head(10)
            lastrow_plot_top = lastrow.sort_values(ascending=False).head(1)
            if self.summary:
                with self.summary:
                    self.summary.clear()
                    ui.label(
                        f"Sturgeon classification: {lastrow_plot_top.index[0]} - {lastrow_plot_top.values[0]:.2f}"
                    )
            self.update_sturgeon_plot(
                lastrow_plot.index.to_list(),
                list(lastrow_plot.values),
                "All",
                self.sturgeon_df_store.iloc[-1]["number_probes"],
            )

    async def process_bam(self, bamfile: List[Tuple[str, float]]) -> None:
        """
        Processes the BAM files and performs the analysis.

        Args:
            bamfile (List[Tuple[str, float]]): List of BAM files with their timestamps.
        """
        sampleID = self.sampleID
        # Initialize directories for each sampleID if not already present
        if sampleID not in self.dataDir.keys():
            self.dataDir[sampleID] = tempfile.TemporaryDirectory(
                dir=self.check_and_create_folder(self.output, sampleID)
            )
            self.bedDir[sampleID] = tempfile.TemporaryDirectory(
                dir=self.check_and_create_folder(self.output, sampleID)
            )
        tomerge = []
        latest_file = 0
        while len(bamfile) > 0:
            self.running = True
            (file, filetime) = bamfile.pop(0)
            if filetime > latest_file:
                latest_file = filetime
            tomerge.append(file)
            app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                "bams_in_processing"
            ] += 1
            if len(tomerge) > 100:
                break

        if latest_file:
            currenttime = latest_file * 1000
        else:
            currenttime = time.time() * 1000

        if len(tomerge) > 0:
            tempbam = tempfile.NamedTemporaryFile(
                dir=self.check_and_create_folder(self.output, sampleID)
            )

            await run.cpu_bound(pysam_cat, tempbam.name, tomerge)

            file = tempbam.name
            temp = tempfile.NamedTemporaryFile(
                dir=self.check_and_create_folder(self.output, sampleID)
            )
            with tempfile.TemporaryDirectory(
                dir=self.check_and_create_folder(self.output, sampleID)
            ) as temp2:
                await run.cpu_bound(run_modkit, file, temp.name, self.threads)

                await run.cpu_bound(run_sturgeon_inputtobed, temp.name, temp2)

                calls_per_probe_file = os.path.join(
                    temp2, "merged_probes_methyl_calls.txt"
                )
                merged_output_file = os.path.join(
                    self.dataDir[sampleID].name,
                    "_merged_probes_methyl_calls.txt",
                )

                if sampleID not in self.first_run.keys():
                    self.first_run[sampleID] = True
                    shutil.copyfile(calls_per_probe_file, merged_output_file)
                else:
                    await run.cpu_bound(
                        run_sturgeon_merge_probes,
                        calls_per_probe_file,
                        merged_output_file,
                    )

                bed_output_file = os.path.join(
                    self.bedDir[sampleID].name, "final_merged_probes_methyl_calls.bed"
                )

                await run.cpu_bound(
                    run_probes_methyl_calls, merged_output_file, bed_output_file
                )

                await run.cpu_bound(
                    run_sturgeon_predict,
                    self.bedDir[sampleID].name,
                    self.dataDir[sampleID].name,
                    self.modelfile,
                )
                if os.path.exists(
                    os.path.join(
                        self.dataDir[sampleID].name,
                        "final_merged_probes_methyl_calls_general.csv",
                    )
                ):
                    mydf = pd.read_csv(
                        os.path.join(
                            self.dataDir[sampleID].name,
                            "final_merged_probes_methyl_calls_general.csv",
                        )
                    )
                else:
                    self.running = False
                    return

                self.st_num_probes[sampleID] = mydf.iloc[-1]["number_probes"]
                # lastrow = mydf.iloc[-1].drop("number_probes")
                mydf_to_save = mydf
                mydf_to_save["timestamp"] = currenttime

                if sampleID not in self.sturgeon_df_store:
                    self.sturgeon_df_store[sampleID] = pd.DataFrame()

                # Exclude empty or all-NA entries before concatenation
                if not mydf_to_save.dropna(how="all").empty:
                    self.sturgeon_df_store[sampleID] = pd.concat(
                        [
                            self.sturgeon_df_store[sampleID],
                            mydf_to_save.set_index("timestamp"),
                        ]
                    )
                    self.sturgeon_df_store[sampleID].to_csv(
                        os.path.join(
                            self.check_and_create_folder(self.output, sampleID),
                            "sturgeon_scores.csv",
                        )
                    )
            app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                "bam_processed"
            ] += len(tomerge)

            app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                "bams_in_processing"
            ] -= len(tomerge)

        self.running = False

    def create_sturgeon_chart(self, title):
        self.echart2 = self.create_chart(title)

    def update_sturgeon_plot(self, x, y, count, st_num_probes):
        """
        Replaces the data in the RapidCNS2 plot.
        :param x: list of tumour types
        :param y: confidence scores for each tumour type
        :param count: the number of bams used to generate the plot
        :return:
        """
        self.echart2.options["title"][
            "text"
        ] = f"Sturgeon: processed {count} bams and found {int(st_num_probes)} probes"
        self.echart2.options["yAxis"]["data"] = x
        self.echart2.options["series"] = [
            {"type": "bar", "name": "Sturgeon", "data": y}
        ]
        self.echart2.update()

    def create_sturgeon_time_chart(self, title):
        self.sturgeon_time_chart = self.create_time_chart(title)

    def update_sturgeon_time_chart(self, datadf):
        """

        :param datadf: the data to plot
        :return:
        """
        self.sturgeon_time_chart.options["series"] = []
        for series, data in datadf.to_dict().items():
            # print(series)
            data_list = [[key, value] for key, value in data.items()]
            # print(data_list)
            if series != "number_probes":
                self.sturgeon_time_chart.options["series"].append(
                    {
                        "animation": False,
                        "type": "line",
                        "smooth": True,
                        "name": series,
                        "emphasis": {"focus": "series"},
                        "endLabel": {
                            "show": True,
                            "formatter": "{a}",
                            "distance": 20,
                        },
                        "lineStyle": {
                            "width": 2,
                        },
                        "data": data_list,
                    }
                )
        self.sturgeon_time_chart.update()


def test_me(
    port: int,
    threads: int,
    watchfolder: str,
    output: str,
    reload: bool = False,
    browse: bool = False,
):
    my_connection = None
    with theme.frame("Sturgeon Rapid CNS Diagnostic.", my_connection):
        TestObject = Sturgeon_object(threads, output, progress=True, batch=True)
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
        print("Browse mode not implemented.")
        TestObject.progress_trackers.visible = False
        # TestObject.show_previous_data(output)
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
def mainrun(port, threads, watchfolder, output, browse):
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
    mainrun()
