from __future__ import annotations
from cnsmeth.subpages.base_analysis import BaseAnalysis

from nicegui import ui, run
import time
import os
import sys
import pysam
import pandas as pd
import shutil
import asyncio
import tempfile
from cnsmeth import models, theme, resources
from cnsmeth.submodules.nanoDX.workflow.scripts.NN_model import NN_classifier
from cnsmeth.utilities.merge_bedmethyl import (
    merge_bedmethyl,
    save_bedmethyl,
    collapse_bedmethyl,
)


def run_modkit(cpgs, sortfile, temp):
    """
    This function runs modkit on a bam file and extracts the methylation data.
    """
    try:
        os.system(
            f"modkit pileup --include-bed {cpgs} --filter-threshold 0.73 --combine-mods --only-tabs -t 4 {sortfile} {temp} --suppress-progress"
        )
    except Exception as e:
        print(e)


def run_samtools_sort(file, tomerge, sortfile):
    pysam.cat("-o", file, *tomerge)
    pysam.sort("--write-index", "-o", sortfile, file)




class NanoDX_object(BaseAnalysis):
    def __init__(self, *args, **kwargs):
        self.cpgs_file = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)),
            "hglft_genome_260e9_91a970_clean.bed",
        )
        self.cpgs = pd.read_csv(
            self.cpgs_file,
            sep="\t",
            header=None,
        )
        self.threshold = 0.05
        self.nanodx_bam_count = 0
        self.not_first_run = False
        self.nanodxfile = tempfile.NamedTemporaryFile()
        self.modelfile = os.path.join(
            os.path.dirname(os.path.abspath(models.__file__)), "Capper_et_al_NN.pkl"
        )
        self.nanodx_df_store = pd.DataFrame()
        self.NN = NN_classifier(self.modelfile)
        super().__init__(*args, **kwargs)

    def setup_ui(self):
        with ui.card().style("width: 100%"):
            with ui.grid(columns=8).classes("w-full h-auto"):
                with ui.card().classes("col-span-3"):
                    self.create_nanodx_chart("NanoDX")
                with ui.card().classes("col-span-5"):
                    self.create_nanodx_time_chart()

    async def process_bam(self, bamfile):
        tomerge = []
        timestamp = None
        while len(bamfile) > 0:
            self.running = True
            (file, filetime) = bamfile.pop()
            self.nanodx_bam_count += 1
            tomerge.append(file)
            timestamp = filetime
            if len(tomerge) > 50:
                break
        if len(tomerge) > 0:
            tempbam = tempfile.NamedTemporaryFile()
            sorttempbam = tempfile.NamedTemporaryFile()
            file = tempbam.name

            temp = tempfile.NamedTemporaryFile()

            sortfile = sorttempbam.name

            ui.notify("NanoDX: Merging bams")

            await run.cpu_bound(run_samtools_sort, file, tomerge, sortfile)

            ui.notify("NanoDX: Running modkit")

            await run.cpu_bound(run_modkit, self.cpgs_file, sortfile, temp.name)

            ui.notify("NanoDX: Merging bed files")

            if self.not_first_run:
                bed_a = pd.read_table(
                    f"{temp.name}",
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
                save_bedmethyl(self.merged_bed_file, self.nanodxfile.name)
            else:
                shutil.copy(f"{temp.name}", self.nanodxfile.name)
                self.merged_bed_file = pd.read_table(
                    self.nanodxfile.name,
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

                self.not_first_run = True

            self.merged_bed_file = collapse_bedmethyl(self.merged_bed_file)

            test_df = pd.merge(
                self.merged_bed_file,
                self.cpgs,
                left_on=["chrom", "start_pos"],
                right_on=[0, 1],
            )
            test_df.rename(
                columns={3: "probe_id", "fraction": "methylation_call"},
                inplace=True,
            )
            test_df.loc[test_df["methylation_call"] < 60, "methylation_call"] = -1
            test_df.loc[test_df["methylation_call"] >= 60, "methylation_call"] = 1

            try:
                predictions, class_labels, n_features = self.NN.predict(test_df)
            except Exception as e:
                print(e)
                test_df.to_csv("errordf.csv", sep=",", index=False, encoding="utf-8")
                # self.nanodx_status_txt["message"] = "Error generating predictions."
                sys.exit(1)

            nanoDX_df = pd.DataFrame({"class": class_labels, "score": predictions})

            nanoDX_save = nanoDX_df.set_index("class").T
            nanoDX_save["number_probes"] = n_features

            if timestamp:
                nanoDX_save["timestamp"] = timestamp * 1000
            else:
                nanoDX_save["timestamp"] = time.time() * 1000

            self.nanodx_df_store = pd.concat(
                [self.nanodx_df_store, nanoDX_save.set_index("timestamp")]
            )

            # self.nanodx_df_store.to_csv(
            #    os.path.join(self.resultfolder, "nanoDX_scores.csv")
            # )

            columns_greater_than_threshold = (
                self.nanodx_df_store > self.threshold
            ).any()
            columns_not_greater_than_threshold = ~columns_greater_than_threshold
            result = self.nanodx_df_store.columns[
                columns_not_greater_than_threshold
            ].tolist()

            self.update_nanodx_time_chart(self.nanodx_df_store.drop(columns=result))

            self.update_nanodx_plot(
                nanoDX_df["class"].head(10).values,
                nanoDX_df["score"].head(10).values,
                self.nanodx_bam_count,
                n_features,
            )

            self.bam_processed += len(tomerge)
            self.bams_in_processing -= len(tomerge)
        await asyncio.sleep(5)
        self.running = False

    def create_nanodx_chart(self, title):
        self.nanodxchart = (
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
            .style("height: 350px")
            .classes("border-double")
        )

    def update_nanodx_plot(self, x, y, count, n_features):
        """
        Replaces the data in the RapidCNS2 plot.
        :param x: list of tumour types
        :param y: confidence scores for each tumour type
        :param count: the number of bams used to generate the plot
        :param n_feature: the number of features detected during data analysis
        :return:
        """
        self.nanodxchart.options["title"][
            "text"
        ] = f"NanoDX: processed {count} bams and found {int(n_features)} features"
        self.nanodxchart.options["yAxis"]["data"] = x
        self.nanodxchart.options["series"] = [
            {"type": "bar", "name": "NanoDX", "data": y}
        ]
        self.nanodxchart.update()

    def create_nanodx_time_chart(self):
        self.nanodx_time_chart = (
            ui.echart(
                {
                    "animation": False,
                    "grid": {"containLabel": True},
                    "title": {"text": "NanoDX Over Time"},
                    "toolbox": {"show": True, "feature": {"saveAsImage": {}}},
                    "xAxis": {"type": "time"},
                    "yAxis": {"type": "value", "data": [], "inverse": False},
                    "series": [],
                }
            )
            .style("height: 350px")
            .classes("border-double")
        )

    def update_nanodx_time_chart(self, datadf):
        """

        :param datadf: the data to plot
        :return:
        """
        self.nanodx_time_chart.options["series"] = []
        for series, data in datadf.to_dict().items():
            # print(series)
            data_list = [[key, value] for key, value in data.items()]
            # print(data_list)
            if series != "number_probes":
                self.nanodx_time_chart.options["series"].append(
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
        self.nanodx_time_chart.update()


def test_ui():
    my_connection = None
    with theme.frame("Copy Number Variation Interactive", my_connection):
        TestObject = NanoDX_object(progress=True, batch=True)
    path = "/users/mattloose/datasets/ds1305_Intraop0006_A/20231123_1233_P2S-00770-A_PAS59057_b1e841e7/bam_pass"
    # path = "tests/static/bam"
    directory = os.fsencode(path)
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".bam"):
            TestObject.add_bam(os.path.join(path, filename))
            time.sleep(0.001)


def start():
    test_ui()
    ui.run(port=8082)


if __name__ in ("__main__", "__mp_main__"):
    start()
