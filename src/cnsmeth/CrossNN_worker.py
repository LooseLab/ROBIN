# Python imports.
from __future__ import annotations
from nicegui import ui
import threading
import time
import os, sys
import pysam
import pandas as pd
import shutil
import tempfile
from cnsmeth import models, theme, resources
from cnsmeth.submodules.nanoDX.workflow.scripts.NN_model import NN_classifier
from cnsmeth.utilities.merge_bedmethyl import (
    merge_bedmethyl,
    save_bedmethyl,
    collapse_bedmethyl,
)

from queue import Queue


class CrossNN_worker:
    def __init__(
        self,
        bamqueue,
        threads=4,
        output_folder=None,
        threshold=0.05,
        showerrors=False,
        browse=False,
    ):
        self.offset = None
        self.browse = browse
        self.bamqueue = bamqueue
        self.threads = threads
        self.threshold = threshold
        self.nanodx_bam_count = 0
        self.outputfolder = output_folder
        self.bedfoldercount = os.path.join(self.outputfolder, "bedscount")
        self.nanodxfolder = os.path.join(self.outputfolder, "nanodx")
        self.resultfolder = os.path.join(self.outputfolder, "results")
        self.running = False
        # ToDo: remove hard coded path
        self.cpgs_file = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)),
            "hglft_genome_260e9_91a970_clean.bed",
        )
        self.cpgs = pd.read_csv(
            self.cpgs_file,
            sep="\t",
            header=None,
        )

        self.modelfile = os.path.join(
            os.path.dirname(os.path.abspath(models.__file__)), "Capper_et_al_NN.pkl"
        )
        self.showerrors = showerrors
        self.result = None
        self.nanodx_df_store = pd.DataFrame()
        self.nanodx_status_txt = {"message": "Waiting for data."}

        self.NN = NN_classifier(self.modelfile)
        if not self.browse:
            if not os.path.exists(self.bedfoldercount):
                os.makedirs(self.bedfoldercount)

            if not os.path.exists(self.resultfolder):
                os.mkdir(self.resultfolder)

            if not os.path.exists(self.nanodxfolder):
                os.mkdir(self.nanodxfolder)

            self.nanodx_processing = threading.Thread(target=self.nanodx, args=())
            self.nanodx_processing.daemon = True
            self.nanodx_processing.start()

    def nanodx(self) -> None:
        """
        This function runs sturgeon on the bam files.
        It grabs all the bam files available in the bamforsturgeon list and cats them together.
        It then runs modkit extract and sturgeon inputtobed to generate sturgeon bed files
        (note these are not true bed files)
        It then merges the bed files together and runs sturgeon predict to
        generate the final sturgeon predictions.
        :param self:
        :return:
        """
        not_first_run = False
        run_count = 0
        self.nanodxfile = os.path.join(self.nanodxfolder, f"nanodx_{run_count}.bed")
        start_time = time.time()
        # worker = get_current_worker()
        while True:
            if self.bamqueue.qsize() > 0:
                self.running = True
                self.nanodxfinished = False
                bams = []

                self.nanodx_status_txt["message"] = (
                    f"Processing {self.bamqueue.qsize()} bam files for NanoDX."
                )

                # Pull all bam files from the queue and cat them together into a temporary file.
                while self.bamqueue.qsize() > 0:
                    bams.append(self.bamqueue.get())
                    self.nanodx_bam_count += 1
                run_count += 1
                tempbam = tempfile.NamedTemporaryFile()
                sorttempbam = tempfile.NamedTemporaryFile()
                file = tempbam.name
                sortfile = sorttempbam.name

                pysam.cat("-o", file, *bams)
                pysam.sort("--write-index", "-o", sortfile, file)

                self.nanodx_status_txt["message"] = (
                    "Bam files merged and sorted - starting extraction."
                )

                # Run modkit extract and sturgeon inputtobed to generate sturgeon bed files
                # Note: these are not the same as those required for NanoDX

                temp = tempfile.NamedTemporaryFile()

                try:
                    # os.system(
                    #    f"modkit extract --ignore h -t {int(self.threads/2)} {file} {temp.name} "
                    #    f"--force --suppress-progress >/dev/null 2>&1"
                    # )
                    os.system(
                        # f"modkit pileup --include-bed {self.cpgs_file} --filter-threshold 0.73 --combine-mods --cpg --reference /Users/mattloose/references/hg38_simple.fa --combine-strands --only-tabs -t 8 {sortfile} {temp.name} --suppress-progress"
                        f"modkit pileup --include-bed {self.cpgs_file} --filter-threshold 0.73 --combine-mods --only-tabs -t 8 {sortfile} {temp.name} --suppress-progress"
                    )
                    # self.log("Done processing bam file")
                except Exception as e:
                    print(e)
                    # self.log(e)
                    pass

                if not_first_run:
                    self.nanodx_status_txt["message"] = (
                        "Merging bed file with previous bed files."
                    )
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
                    save_bedmethyl(self.merged_bed_file, self.nanodxfile)
                else:
                    shutil.copy(f"{temp.name}", self.nanodxfile)
                    self.merged_bed_file = pd.read_table(
                        self.nanodxfile,
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

                    not_first_run = True

                self.nanodx_status_txt["message"] = "Prediction underway."

                self.merged_bed_file = collapse_bedmethyl(self.merged_bed_file)

                ### Adding in NanoDX predict
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
                # predictions, class_labels, n_features = self.NN.predict_from_bedMethyl("sorted.bam.CpG.450K.2.bed")

                try:
                    predictions, class_labels, n_features = self.NN.predict(test_df)
                except Exception as e:
                    print(e)
                    test_df.to_csv("errordf.csv", sep=',', index=False, encoding='utf-8')
                    self.nanodx_status_txt["message"] = "Error generating predictions."
                    sys.exit(1)

                nanoDX_df = pd.DataFrame({"class": class_labels, "score": predictions})

                nanoDX_save = nanoDX_df.set_index("class").T
                nanoDX_save["number_probes"] = n_features

                if self.offset:
                    self.offset = time.time() - start_time
                    nanoDX_save["timestamp"] = (time.time() + self.offset) * 1000
                else:
                    nanoDX_save["timestamp"] = time.time() * 1000

                self.nanodx_df_store = pd.concat(
                    [self.nanodx_df_store, nanoDX_save.set_index("timestamp")]
                )

                self.nanodx_df_store.to_csv(
                    os.path.join(self.resultfolder, "nanoDX_scores.csv")
                )

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

                self.nanodx_status_txt["message"] = (
                    "Predictions generated. Waiting for data."
                )
            self.running = False

            time.sleep(5)

            if self.bamqueue.qsize() == 0:
                self.nanodxfinished = True

    def status_nanodx(self):
        ui.label().bind_text_from(
            self.nanodx_status_txt,
            "message",
            backward=lambda n: f"NanoDX Status: {n}",
        )

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

    def playback_thread(self, data: pd.DataFrame):
        self.data = data
        nanodx_thread_processing = threading.Thread(
            target=self.playback_nanodx, args=()
        )
        nanodx_thread_processing.daemon = True
        nanodx_thread_processing.start()

    def playback_nanodx(self):
        latest_timestamps = self.data
        self.offset = 0
        for index, row in latest_timestamps.iterrows():
            current_time = time.time()
            print(f"Current time: {current_time}")
            print(f"Offset: {self.offset}")

            time_diff = row["file_produced"] - current_time - self.offset
            print(f"TimeDiff: {time_diff}")

            if self.offset == 0:
                self.offset = time_diff
                time_diff = 0

            while row["file_produced"] - current_time - self.offset > 0:
                time.sleep(0.1)
                if not self.running:
                    self.offset += 2
            self.bamqueue.put(row["full_path"])
            # self.bam_count["counter"] += 1
        # if "file" not in self.bam_count:
        #    self.bam_count["file"] = {}
        # self.bam_count["file"][row["full_path"]] = time.time()
        ## Expect a pandas dataframe containing files and times

    """
    def load_prior_data(self):
        self.sturgeon_df_store = pd.read_csv(
            os.path.join(self.resultfolder, "sturgeon_scores.csv")
        ).set_index("timestamp")
        columns_greater_than_threshold = (self.sturgeon_df_store > self.threshold).any()
        columns_not_greater_than_threshold = ~columns_greater_than_threshold
        result = self.sturgeon_df_store.columns[
            columns_not_greater_than_threshold
        ].tolist()

        self.update_sturgeon_time_chart(self.sturgeon_df_store.drop(columns=result))
        mydf = pd.read_csv(
            os.path.join(
                self.resultfolder,
                "final_merged_probes_methyl_calls_general.csv",
            )
        )

        self.st_num_probes = mydf.iloc[-1]["number_probes"]
        lastrow = mydf.iloc[-1].drop("number_probes")
        lastrow_plot = lastrow.sort_values(ascending=False).head(10)
        self.update_sturgeon_plot(
            lastrow_plot.index.to_list(),
            list(lastrow_plot.values),
            self.nanodx_bam_count,
        )

        self.sturgeon_status_txt["message"] = "Predictions Complete."

    def replay_prior_data(self):
        self.background_task = threading.Thread(target=self._replay_prior_data, args=())
        self.background_task.daemon = True
        self.background_task.start()

    def _replay_prior_data(self):
        '''
        Replay prior data from a file.
        :return:
        '''
        print("Replaying prior Sturgeon data")
        self.sturgeon_status_txt["message"] = "Replaying prior Sturgeon data."
        self.sturgeon_df_store = pd.read_csv(
            os.path.join(self.resultfolder, "sturgeon_scores.csv")
        )
        self.sturgeon_df_store["offset"] = (
            self.sturgeon_df_store["timestamp"] - self.sturgeon_df_store["timestamp"][0]
        )
        self.sturgeon_df_store.set_index("timestamp")

        self.sturgeon_status_txt["message"] = "Replaying prior data."

        scale_factor = 1200
        counter = 0
        current_time = time.time()
        for index, row in self.sturgeon_df_store.iterrows():
            counter += 1
            # print(index, row)
            if time.time() < current_time + (row["offset"] / 1000 / scale_factor):
                time.sleep(
                    current_time + (row["offset"] / 1000 / scale_factor) - time.time()
                )
            else:
                time.sleep(0.001)
            # time.sleep(row["offset"] / 1000 / scale_factor)

            temp_sturgeon_df_store = (
                self.sturgeon_df_store.head(counter)
                .drop(columns=["offset", "number_probes"])
                .set_index("timestamp")
            )
            columns_greater_than_threshold = (
                temp_sturgeon_df_store > self.threshold
            ).any()
            columns_not_greater_than_threshold = ~columns_greater_than_threshold
            result = temp_sturgeon_df_store.columns[
                columns_not_greater_than_threshold
            ].tolist()

            self.update_sturgeon_time_chart(temp_sturgeon_df_store.drop(columns=result))
            self.sturgeon_status_txt["message"] = (
                f"Replaying prior RCNS2 data - step {counter}."
            )
            lastrow = row.drop(["offset", "number_probes", "timestamp"])
            lastrow_plot = lastrow.sort_values(ascending=False).head(10)
            self.update_sturgeon_plot(
                lastrow_plot.index.to_list(),
                list(lastrow_plot.values),
                "replay",
            )
        self.sturgeon_status_txt["message"] = "Viewing historical Sturgeon data."
    """


def index_page() -> None:
    # from check_connection import ConnectionDialog
    # initial_ip = "127.0.0.1"
    # my_connection = ConnectionDialog(initial_ip)
    my_connection = None
    with theme.frame("MethClass Interactive", my_connection):
        # my_connection.connect_to_minknow()
        ui.label("Hello")

        bamfornanodx = Queue()
        nanodx_worker = CrossNN_worker(bamfornanodx, threads=4, output_folder="/tmp/")
        with ui.card().style("width: 100%"):
            nanodx_worker.status_nanodx()
            nanodx_worker.create_nanodx_chart("NanoDX")
        with ui.card().style("width: 100%"):
            nanodx_worker.create_nanodx_time_chart()

        bamfornanodx.put(
            "/Users/mattloose/GIT/niceGUI/cnsmeth/tests/static/bam/test_set.hg38.aa.sam.bam"
        )
        time.sleep(5)
        bamfornanodx.put(
            "/Users/mattloose/GIT/niceGUI/cnsmeth/tests/static/bam/test_set.hg38.ab.sam.bam"
        )


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
