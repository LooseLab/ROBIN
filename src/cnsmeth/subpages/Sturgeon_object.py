from cnsmeth.subpages.base_analysis import BaseAnalysis
import os
import tempfile
import time
import shutil
import pandas as pd
from nicegui import ui, run
from cnsmeth import theme
import pysam
from cnsmeth import models
from sturgeon.callmapping import (
    merge_probes_methyl_calls,
    probes_methyl_calls_to_bed,
)
import asyncio


def run_probes_methyl_calls(merged_output_file, bed_output_file):
    probes_methyl_calls_to_bed(merged_output_file, bed_output_file)


def run_sturgeon_merge_probes(calls_per_probe_file, merged_output_file):
    merge_probes_methyl_calls(
        [calls_per_probe_file, merged_output_file],
        merged_output_file,
    )


def pysam_cat(tempbam, tomerge):
    pysam.cat("-o", tempbam, *tomerge)


def run_modkit(file, temp):
    """
    This function runs modkit on a bam file and extracts the methylation data.
    """
    try:
        os.system(
            f"modkit extract --ignore h -t {4} {file} {temp} "
            f"--force --suppress-progress >/dev/null 2>&1"
        )
        # self.log("Done processing bam file")
    except Exception as e:
        print(e)
        # self.log(e)
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


async def test_hello(file, temp):
    await run.cpu_bound(run_modkit, file, temp.name)
    ui.notify("hello")
    # await handle_click()


class Sturgeon_object(BaseAnalysis):
    def __init__(self, *args, **kwargs):
        self.dataDir = tempfile.TemporaryDirectory()
        self.bedDir = tempfile.TemporaryDirectory()
        self.sturgeon_df_store = pd.DataFrame()
        self.threshold = 0.05
        self.first_run = True
        self.loop = asyncio.get_event_loop()
        self.modelfile = os.path.join(
            os.path.dirname(os.path.abspath(models.__file__)), "general.zip"
        )
        super().__init__(*args, **kwargs)

    def setup_ui(self):
        self.card = ui.card().style("width: 100%")
        with self.card:
            with ui.grid(columns=8).classes("w-full h-auto"):
                with ui.card().classes("col-span-3"):
                    self.create_sturgeon_chart("Sturgeon")
                with ui.card().classes("col-span-5"):
                    self.create_sturgeon_time_chart()

    async def process_bam(self, bamfile):
        tomerge = []
        timestamp = None

        while len(bamfile) > 0:
            self.running = True
            (file, filetime) = bamfile.pop()
            tomerge.append(file)
            timestamp = filetime
            if len(tomerge) > 25:
                break

        if len(tomerge) > 0:
            tempbam = tempfile.NamedTemporaryFile()
            with self.card:
                ui.notify("Sturgeon: Merging bams")
            await run.cpu_bound(pysam_cat, tempbam.name, tomerge)

            file = tempbam.name

            temp = tempfile.NamedTemporaryFile()

            with tempfile.TemporaryDirectory() as temp2:
                await run.cpu_bound(run_modkit, file, temp.name)
                ui.notify("Sturgeon: Modkit Complete")
                await run.cpu_bound(run_sturgeon_inputtobed, temp.name, temp2)
                # ui.notify("Sturgeon: Inputtobed Complete")

                calls_per_probe_file = os.path.join(
                    temp2, "merged_probes_methyl_calls.txt"
                )
                merged_output_file = os.path.join(
                    self.dataDir.name,
                    "_merged_probes_methyl_calls.txt",
                )

                if self.first_run:
                    shutil.copyfile(calls_per_probe_file, merged_output_file)

                    self.first_run = False
                else:
                    await run.cpu_bound(
                        run_sturgeon_merge_probes,
                        calls_per_probe_file,
                        merged_output_file,
                    )
                    # merge_probes_methyl_calls(
                    #    [calls_per_probe_file, merged_output_file],
                    #    merged_output_file,
                    # )
                bed_output_file = os.path.join(
                    self.bedDir.name, "final_merged_probes_methyl_calls.bed"
                )

                await run.cpu_bound(
                    run_probes_methyl_calls, merged_output_file, bed_output_file
                )
                # probes_methyl_calls_to_bed(merged_output_file, bed_output_file)

                await run.cpu_bound(
                    run_sturgeon_predict,
                    self.bedDir.name,
                    self.dataDir.name,
                    self.modelfile,
                )

                ui.notify("Sturgeon: Prediction Complete")

                mydf = pd.read_csv(
                    os.path.join(
                        self.dataDir.name,
                        "final_merged_probes_methyl_calls_general.csv",
                    )
                )
                self.st_num_probes = mydf.iloc[-1]["number_probes"]
                lastrow = mydf.iloc[-1].drop("number_probes")
                lastrow_plot = lastrow.sort_values(ascending=False).head(10)

                mydf_to_save = mydf

                if timestamp:
                    mydf_to_save["timestamp"] = timestamp * 1000
                else:
                    mydf_to_save["timestamp"] = time.time() * 1000

                self.bam_processed += len(tomerge)
                self.bams_in_processing -= len(tomerge)

                self.sturgeon_df_store = pd.concat(
                    [self.sturgeon_df_store, mydf_to_save.set_index("timestamp")]
                )

                columns_greater_than_threshold = (
                    self.sturgeon_df_store > self.threshold
                ).any()
                columns_not_greater_than_threshold = ~columns_greater_than_threshold
                result = self.sturgeon_df_store.columns[
                    columns_not_greater_than_threshold
                ].tolist()

                self.update_sturgeon_time_chart(
                    self.sturgeon_df_store.drop(columns=result)
                )

                self.update_sturgeon_plot(
                    lastrow_plot.index.to_list(),
                    list(lastrow_plot.values),
                    self.bam_processed,
                    self.st_num_probes,
                )

        await asyncio.sleep(5)
        self.running = False

    def create_sturgeon_chart(self, title):
        self.echart2 = (
            ui.echart(
                {
                    "animation": False,
                    "grid": {"containLabel": True},
                    "title": {"text": title},
                    "toolbox": {"show": True, "feature": {"saveAsImage": {}}},
                    "xAxis": {"type": "value", "max": 1},
                    "yAxis": {"type": "category", "data": [], "inverse": True},
                    #'legend': {},
                    "series": [],
                }
            )
            .style("height: 320px")
            .classes("border-double")
        )

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

    def create_sturgeon_time_chart(self):
        self.sturgeon_time_chart = (
            ui.echart(
                {
                    "animation": False,
                    "grid": {"containLabel": True},
                    "title": {"text": "Sturgeon Over Time"},
                    "toolbox": {"show": True, "feature": {"saveAsImage": {}}},
                    "xAxis": {"type": "time"},
                    "yAxis": {"type": "value", "data": [], "inverse": False},
                    "series": [],
                }
            )
            .style("height: 320px")
            .classes("border-double")
        )

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


def test_ui():
    my_connection = None
    with theme.frame("Copy Number Variation Interactive", my_connection):
        TestObject = Sturgeon_object(progress=True, batch=True)
    path = "/users/mattloose/datasets/ds1305_Intraop0006_A/20231123_1233_P2S-00770-A_PAS59057_b1e841e7/bam_pass"
    # path = "tests/static/bam"
    directory = os.fsencode(path)
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".bam"):
            TestObject.add_bam(os.path.join(path, filename))
            time.sleep(0.001)
    ui.run(port=8082)


def start():
    test_ui()


if __name__ in ("__main__", "__mp_main__"):
    start()
