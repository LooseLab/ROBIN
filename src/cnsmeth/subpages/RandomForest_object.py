from cnsmeth.subpages.base_analysis import BaseAnalysis
import os
import tempfile
import time
import shutil
import pandas as pd
from nicegui import ui, run
from cnsmeth import theme, resources
import pysam
from cnsmeth import models
from sturgeon.callmapping import (
    merge_probes_methyl_calls,
    probes_methyl_calls_to_bed,
)

import yappi, tabulate

import asyncio
from cnsmeth import submodules

from cnsmeth.utilities.merge_bedmethyl import (
    merge_bedmethyl,
    save_bedmethyl,
    collapse_bedmethyl,
)

HVPATH = os.path.join(
    os.path.dirname(os.path.abspath(submodules.__file__)), "hv_rapidCNS2"
)

def run_probes_methyl_calls(merged_output_file, bed_output_file):
    probes_methyl_calls_to_bed(merged_output_file, bed_output_file)
def run_sturgeon_merge_probes(calls_per_probe_file,merged_output_file):
    merge_probes_methyl_calls(
        [calls_per_probe_file, merged_output_file],
        merged_output_file,
    )

def run_rcns2(rcns2folder, batch, bed, threads, showerrors):
    command = (
            f"Rscript {HVPATH}/bin/methylation_classification_nanodx_v0.1.R -s "
            + f"live_{batch} -o {rcns2folder} -i {bed} "
            + f"-p {HVPATH}/bin/top_probes_hm450.Rdata "
            + f"--training_data {HVPATH}/bin/capper_top_100k_betas_binarised.Rdata "
            + f"--array_file {HVPATH}/bin/HM450.hg38.manifest.gencode.v22.Rdata "
            + f"-t {threads} "
    )
    if not showerrors:
        command += ">/dev/null 2>&1"

    os.system(command)
    #returned_value = subprocess.call(
    #    command,
    #    shell=True,
    #)
    if showerrors:
        print(command)
        print(returned_value)

def run_samtools_sort(file, tomerge, sortfile):
    pysam.cat("-o", file, *tomerge)
    pysam.sort("--write-index", "-o", sortfile, file)

def run_modkit(bamfile, outbed):
    """
    This function runs modkit on a bam file and extracts the methylation data.
    """
    try:
        os.system(
            f"modkit pileup -t {8} --filter-threshold 0.73 --combine-mods {bamfile} "
            f"{outbed} #-suppress-progress  >/dev/null 2>&1 "
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
    #await handle_click()

class RandomForest_object(BaseAnalysis):
    def __init__(self, *args, **kwargs):
        self.dataDir = tempfile.TemporaryDirectory()
        self.bedDir = tempfile.TemporaryDirectory()
        self.rcns2_df_store = pd.DataFrame()
        self.threshold=0.05
        self.batch = 0
        self.offset = False
        self.first_run = True
        self.showerrors = False
        self.modelfile = os.path.join(
            os.path.dirname(os.path.abspath(models.__file__)), "general.zip"
        )
        super().__init__(*args, **kwargs)


    def setup_ui(self):
        self.card =  ui.card().style("width: 100%")
        with self.card:
            with ui.grid(columns=8).classes("w-full h-auto"):
                with ui.card().classes('col-span-3'):
                    self.create_rcns2_chart("Random Forest")
                with ui.card().classes('col-span-5'):
                    self.create_rcns2_time_chart()



    async def process_bam(self, bamfile):
        tomerge=[]
        timestamp=None

        while len(bamfile)>0:
            self.running=True
            (file,filetime) = bamfile.pop()
            tomerge.append(file)
            timestamp = filetime
            if len(tomerge)>100:
                break

        if len(tomerge)>0:
            tempbam = tempfile.NamedTemporaryFile(suffix=".bam")
            sortbam = tempfile.NamedTemporaryFile(suffix=".bam")
            tempbed = tempfile.NamedTemporaryFile(suffix=".bed")
            self.batch+=1
            await run.cpu_bound(run_samtools_sort, tempbam.name,tomerge,sortbam.name)

            await run.cpu_bound(run_modkit, sortbam.name, tempbed.name)

            if not self.first_run:
                bed_a = pd.read_table(
                    f"{tempbed.name}",
                    names=[
                        "chrom",
                        "start_pos",
                        "end_pos",
                        "mod",
                        "score",
                        "strand",
                        "start_pos2",
                        "end_pos2",
                        "colour",
                        "Nvalid",
                        "fraction",
                        "Nmod",
                        "Ncanon",
                        "Nother",
                        "Ndel",
                        "Nfail",
                        "Ndiff",
                        "Nnocall",
                    ],
                    dtype={
                        "chrom": "category",
                        "start_pos": "int32",
                        "end_pos": "int32",
                        "mod": "category",
                        "score": "int16",
                        "strand": "category",
                        "start_pos2": "int32",
                        "end_pos2": "int32",
                        "colour": "category",
                        "Nvalid": "int16",
                        "fraction": "float16",
                        "Nmod": "int16",
                        "Ncanon": "int16",
                        "Nother": "int16",
                        "Ndel": "int16",
                        "Nfail": "int16",
                        "Ndiff": "int16",
                        "Nnocall": "int16",
                    },
                    header=None,
                    delim_whitespace=True,
                )

                self.merged_bed_file = merge_bedmethyl(bed_a, self.merged_bed_file)
                save_bedmethyl(self.merged_bed_file, f"{tempbed.name}")

            else:
                self.merged_bed_file = pd.read_table(
                    f"{tempbed.name}",
                    names=[
                        "chrom",
                        "start_pos",
                        "end_pos",
                        "mod",
                        "score",
                        "strand",
                        "start_pos2",
                        "end_pos2",
                        "colour",
                        "Nvalid",
                        "fraction",
                        "Nmod",
                        "Ncanon",
                        "Nother",
                        "Ndel",
                        "Nfail",
                        "Ndiff",
                        "Nnocall",
                    ],
                    dtype={
                        "chrom": "category",
                        "start_pos": "int32",
                        "end_pos": "int32",
                        "mod": "category",
                        "score": "int16",
                        "strand": "category",
                        "start_pos2": "int32",
                        "end_pos2": "int32",
                        "colour": "category",
                        "Nvalid": "int16",
                        "fraction": "float16",
                        "Nmod": "int16",
                        "Ncanon": "int16",
                        "Nother": "int16",
                        "Ndel": "int16",
                        "Nfail": "int16",
                        "Ndiff": "int16",
                        "Nnocall": "int16",
                    },
                    header=None,
                    delim_whitespace=True,
                )
                self.first_run = False

            tempDir = tempfile.TemporaryDirectory()

            await run.cpu_bound(run_rcns2, tempDir.name, self.batch, tempbed.name, 8, self.showerrors)

            if os.path.isfile(f"{tempDir.name}/live_{self.batch}_votes.tsv"):
                scores = pd.read_table(
                    f"{tempDir.name}/live_{self.batch}_votes.tsv",
                    delim_whitespace=True,
                )
                scores_to_save = scores.drop(columns=["Freq"]).T

                if self.offset:
                    self.offset = time.time() - start_time
                    scores_to_save["timestamp"] = (time.time() + self.offset) * 1000
                else:
                    scores_to_save["timestamp"] = time.time() * 1000

                self.rcns2_df_store = pd.concat(
                    [self.rcns2_df_store, scores_to_save.set_index("timestamp")]
                )
                #self.rcns2_df_store.to_csv(
                #    os.path.join(self.resultfolder, "rcns2_scores.csv")
                #)

                columns_greater_than_threshold = (
                        self.rcns2_df_store > self.threshold * 100
                ).any()
                columns_not_greater_than_threshold = ~columns_greater_than_threshold
                result = self.rcns2_df_store.columns[
                    columns_not_greater_than_threshold
                ].tolist()

                self.update_rcns2_time_chart(
                    self.rcns2_df_store.drop(columns=result)
                )
                # with self.rcns2_container:
                #    self.rcns2_container.clear()
                #    ui.table.from_pandas(self.rcns2_df_store, pagination=3).classes('max-h-80')
                scores = scores.sort_values(by=["cal_Freq"], ascending=False).head(
                    10
                )
                self.bam_processed += len(tomerge)
                self.update_rcns2_plot(
                    scores.index.to_list(),
                    list(scores["cal_Freq"].values / 100),
                    self.bam_processed,
                )

        await asyncio.sleep(5)
        self.running=False


    def create_rcns2_chart(self, title):
        self.echart = (
            ui.echart(
                {
                    "grid": {"containLabel": True},
                    "title": {"text": title},
                    "toolbox": {"show": True, "feature": {"saveAsImage": {}}},
                    "xAxis": {"type": "value", "max": 1},
                    "yAxis": {"type": "category", "data": [], "inverse": True},
                    # 'legend': {},
                    "series": [],
                }
            )
            .style("height: 350px")
            .classes("border-double")
        )

    def update_rcns2_plot(self, x, y, count):
        """
        Replaces the data in the RapidCNS2 plot.
        :param x: list of tumour types
        :param y: confidence scores for each tumour type
        :param count: the number of bams used to generate the plot
        :return:
        """
        self.echart.options["title"]["text"] = f"RapidCNS2 - processed {count} Bams"
        self.echart.options["yAxis"]["data"] = x
        self.echart.options["series"] = [
            {"type": "bar", "name": "RapidCNS2", "data": y}
        ]
        self.echart.update()

    def create_rcns2_time_chart(self):
        self.rcns2_time_chart = (
            ui.echart(
                {
                    "animation": False,
                    "grid": {"containLabel": True},
                    "title": {"text": "RCNS2 Over Time"},
                    "toolbox": {"show": True, "feature": {"saveAsImage": {}}},
                    "xAxis": {"type": "time"},
                    "yAxis": {"type": "value", "data": [], "inverse": False},
                    #'tooltip': {
                    #    'order': 'valueDesc',
                    #    'trigger': 'axis'
                    # },
                    "series": [],
                }
            )
            .style("height: 350px")
            .classes("border-double")
        )

    def update_rcns2_time_chart(self, datadf):
        """

        :param datadf: the data to plot
        :return:
        """
        self.rcns2_time_chart.options["series"] = []
        for series, data in datadf.to_dict().items():
            data_list = [[key, value] for key, value in data.items()]
            # print(data_list)
            if series != "number_probes":
                self.rcns2_time_chart.options["series"].append(
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
        self.rcns2_time_chart.update()



def test_ui():
    my_connection = None
    with theme.frame("Copy Number Variation Interactive", my_connection):
        ui.button('start', on_click=start)
        ui.button('stop', on_click=stop)
        TestObject = RandomForest_object(progress=True, batch=True)
    path = "/users/mattloose/datasets/ds1305_Intraop0006_A/20231123_1233_P2S-00770-A_PAS59057_b1e841e7/bam_pass"
    #path = "tests/static/bam"
    directory = os.fsencode(path)
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".bam"):
            TestObject.add_bam(os.path.join(path, filename))
            time.sleep(0.001)
    ui.run(port=8082)
def start() -> None:
    yappi.clear_stats()
    yappi.start()

def stop() -> None:
    yappi.stop()
    table = [
        [str(v) for v in [stat.full_name, stat.ttot, stat.tsub, stat.tavg, stat.ncall]]
        for stat in yappi.get_func_stats()
        if 'python' not in stat.module
    ]
    print(tabulate.tabulate(table[:15], headers=['function', 'total', 'excl. sub', 'avg', 'ncall'], floatfmt='.4f'))
    yappi.get_thread_stats().print_all()



if __name__ in ("__main__", "__mp_main__"):
    test_ui()