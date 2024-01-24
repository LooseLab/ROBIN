from nicegui import Tailwind, ui, app
import threading
import time
import os, signal
import pysam
import pandas as pd
import shutil
import tempfile
from methnicegui import models, resources
from sturgeon.callmapping import (
    merge_probes_methyl_calls,
    probes_methyl_calls_to_bed,
)


class Sturgeon_worker:
    def __init__(
        self,
        bamqueue,
        threads=4,
        output_folder=None,
        threshold=0.05,
        showerrors=False,
        browse=False,
    ):
        self.browse = browse
        self.bamqueue = bamqueue
        self.threads = threads
        self.threshold = threshold
        self.st_bam_count = 0
        self.outputfolder = output_folder
        self.bedfoldercount = os.path.join(self.outputfolder, "bedscount")
        self.resultfolder = os.path.join(self.outputfolder, "results")
        self.modelfile = os.path.join(
            os.path.dirname(os.path.abspath(models.__file__)), "general.zip"
        )
        self.showerrors = showerrors
        self.result = None
        self.sturgeon_df_store = pd.DataFrame()
        self.sturgeon_status_txt = {"message": "Waiting for data."}
        if not self.browse:
            if not os.path.exists(self.bedfoldercount):
                os.makedirs(self.bedfoldercount)

            if not os.path.exists(self.resultfolder):
                os.mkdir(self.resultfolder)

            self.sturgeon_processing = threading.Thread(target=self.sturgeon, args=())
            self.sturgeon_processing.daemon = True
            self.sturgeon_processing.start()

    def sturgeon(self) -> None:
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
        first_run = True
        run_count = 0

        # worker = get_current_worker()
        while True:
            if self.bamqueue.qsize() > 0:
                self.stugeonfinished = False
                bams = []
                # self.log(f"there are {len(self.bamforsturgeon)} bams for sturgeon processing")
                self.sturgeon_status_txt[
                    "message"
                ] = f"Processing {self.bamqueue.qsize()} bam files for Sturgeon."
                while self.bamqueue.qsize() > 0:
                    bams.append(self.bamqueue.get())
                    self.st_bam_count += 1
                run_count += 1
                tempbam = tempfile.NamedTemporaryFile()

                pysam.cat("-o", tempbam.name, *bams)

                self.sturgeon_status_txt[
                    "message"
                ] = "Bam files merged - starting extraction."

                file = tempbam.name

                temp = tempfile.NamedTemporaryFile()

                with tempfile.TemporaryDirectory() as temp2:
                    try:
                        os.system(
                            f"modkit extract --ignore h -t {int(self.threads/2)} {file} {temp.name} "
                            f"--force --suppress-progress >/dev/null 2>&1"
                        )
                        # self.log("Done processing bam file")
                    except Exception as e:
                        print(e)
                        # self.log(e)
                        pass
                    try:
                        os.system(
                            f"sturgeon inputtobed -i {temp.name} -o {temp2} -s modkit "
                            f"--reference-genome hg38 >/dev/null 2>&1"
                        )
                        # self.log(temp2)
                    except Exception as e:
                        print(e)
                        # self.log(e)
                        pass

                    self.sturgeon_status_txt[
                        "message"
                    ] = "Sturgeon Bed files generated. Merging with previous."

                    calls_per_probe_file = os.path.join(
                        temp2, "merged_probes_methyl_calls.txt"
                    )

                    merged_output_file = os.path.join(
                        self.bedfoldercount,
                        f"{run_count}_merged_probes_methyl_calls.txt",
                    )
                    previous_merged_output_file = os.path.join(
                        self.bedfoldercount,
                        f"{run_count - 1}_merged_probes_methyl_calls.txt",
                    )

                    if first_run:
                        shutil.copyfile(calls_per_probe_file, merged_output_file)

                        first_run = False
                    else:
                        merge_probes_methyl_calls(
                            [calls_per_probe_file, previous_merged_output_file],
                            merged_output_file,
                        )

                    bed_output_file = os.path.join(
                        self.bedfoldercount, "final_merged_probes_methyl_calls.bed"
                    )
                    probes_methyl_calls_to_bed(merged_output_file, bed_output_file)
                    self.sturgeon_status_txt[
                        "message"
                    ] = "Bed files merged. Generating predictions."
                    os.system(
                        f"sturgeon predict -i {self.bedfoldercount} -o {self.resultfolder} "
                        f"--model-files {self.modelfile} >/dev/null 2>&1"
                    )
                    mydf = pd.read_csv(
                        os.path.join(
                            self.resultfolder,
                            "final_merged_probes_methyl_calls_general.csv",
                        )
                    )

                    self.st_num_probes = mydf.iloc[-1]["number_probes"]
                    lastrow = mydf.iloc[-1].drop("number_probes")
                    lastrow_plot = lastrow.sort_values(ascending=False).head(10)

                    mydf_to_save = mydf
                    mydf_to_save["timestamp"] = time.time() * 1000
                    self.sturgeon_df_store = pd.concat(
                        [self.sturgeon_df_store, mydf_to_save.set_index("timestamp")]
                    )
                    self.sturgeon_df_store.to_csv(
                        os.path.join(self.resultfolder, "sturgeon_scores.csv")
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
                        self.st_bam_count,
                    )

                    self.sturgeon_status_txt[
                        "message"
                    ] = "Predictions generated. Waiting for data."

            time.sleep(5)

            if self.bamqueue.qsize() == 0:
                self.stugeonfinished = True

    def status_sturgeon(self):
        ui.label().bind_text_from(
            self.sturgeon_status_txt,
            "message",
            backward=lambda n: f"Sturgeon Status: {n}",
        )

    def create_sturgeon_chart(self, title):
        self.echart2 = (
            ui.echart(
                {
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

    def update_sturgeon_plot(self, x, y, count):
        """
        Replaces the data in the RapidCNS2 plot.
        :param x: list of tumour types
        :param y: confidence scores for each tumour type
        :param count: the number of bams used to generate the plot
        :return:
        """
        self.echart2.options["title"][
            "text"
        ] = f"Sturgeon: processed {count} bams and found {int(self.st_num_probes)} probes"
        self.echart2.options["yAxis"]["data"] = x
        self.echart2.options["series"] = [
            {"type": "bar", "name": "Sturgeon", "data": y}
        ]
        self.echart2.update()

    def create_sturgeon_time_chart(self):
        self.sturgeon_time_chart = (
            ui.echart(
                {
                    "grid": {"containLabel": True},
                    "title": {"text": "Sturgeon Over Time"},
                    "toolbox": {"show": True, "feature": {"saveAsImage": {}}},
                    "xAxis": {"type": "time"},
                    "yAxis": {"type": "value", "data": [], "inverse": False},
                    "series": [],
                }
            )
            .style("height: 350px")
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
            self.st_bam_count,
        )

        self.sturgeon_status_txt["message"] = "Predictions Complete."

    def replay_prior_data(self):
        self.background_tast = threading.Thread(target=self._replay_prior_data, args=())
        self.background_tast.daemon = True
        self.background_tast.start()

    def _replay_prior_data(self):
        """
        Replay prior data from a file.
        :return:
        """
        print("Replaying prior Sturgeon data")
        self.sturgeon_status_txt["message"] = f"Replaying prior Sturgeon data."
        self.sturgeon_df_store = pd.read_csv(
            os.path.join(self.resultfolder, "sturgeon_scores.csv")
        ).set_index("timestamp")
        min_index_value = self.sturgeon_df_store.index.min()
        max_index_value = self.sturgeon_df_store.index.max()
        start_time = time.time()
        elapsed_time = 0
        scale_factor = 300
        scaled_elapsed_time = elapsed_time * scale_factor

        print(max_index_value - min_index_value)

        df_len = 0

        while (1000 * scaled_elapsed_time) + min_index_value < max_index_value:
            print("replaying_data")
            elapsed_time = time.time() - start_time
            scaled_elapsed_time = elapsed_time * scale_factor
            print(scaled_elapsed_time)

            temp_sturgeon_df_store = self.sturgeon_df_store[
                self.sturgeon_df_store.index
                < (min_index_value + (1000 * scaled_elapsed_time))
            ]

            columns_greater_than_threshold = (
                temp_sturgeon_df_store > self.threshold
            ).any()
            columns_not_greater_than_threshold = ~columns_greater_than_threshold
            result = temp_sturgeon_df_store.columns[
                columns_not_greater_than_threshold
            ].tolist()
            if temp_sturgeon_df_store.drop(columns=result).shape[0] > df_len:
                self.update_sturgeon_time_chart(
                    temp_sturgeon_df_store.drop(columns=result)
                )
                df_len = temp_sturgeon_df_store.drop(columns=result).shape[0]
                self.sturgeon_status_txt[
                    "message"
                ] = f"Replaying prior Sturgeon data - step {df_len}."
                lastrow = temp_sturgeon_df_store.iloc[-1].drop("number_probes")
                lastrow_plot = lastrow.sort_values(ascending=False).head(10)
                self.update_sturgeon_plot(
                    lastrow_plot.index.to_list(),
                    list(lastrow_plot.values),
                    "replay",
                )

            time.sleep(0.5)
        self.sturgeon_status_txt["message"] = f"Viewing historical Sturgeon data."
