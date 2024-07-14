from robin.subpages.base_analysis import BaseAnalysis
import os
import tempfile
import time
import pandas as pd
from nicegui import ui, app, run
from robin import theme, resources
import pysam
from robin import models
from sturgeon.callmapping import (
    merge_probes_methyl_calls,
    probes_methyl_calls_to_bed,
)

import yappi
import tabulate


from robin import submodules

from robin.utilities.merge_bedmethyl import (
    merge_bedmethyl,
    save_bedmethyl,
)
from typing import List, Tuple


HVPATH = os.path.join(
    os.path.dirname(os.path.abspath(submodules.__file__)), "hv_rapidCNS2"
)


def run_probes_methyl_calls(merged_output_file, bed_output_file):
    probes_methyl_calls_to_bed(merged_output_file, bed_output_file)


def run_sturgeon_merge_probes(calls_per_probe_file, merged_output_file):
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
    # print(command)

    os.system(command)


def run_samtools_sort(file, tomerge, sortfile, threads, regions):
    pysam.cat("-o", file, *tomerge)
    # pysam.sort("-@", f"{threads}", "--write-index", "-o", sortfile, file)
    intermediate_bam = tempfile.NamedTemporaryFile(suffix=".bam")
    command = (
        f"samtools sort -@{threads} --write-index -o {intermediate_bam.name} {file}"
    )
    os.system(command)
    command2 = f"samtools view -b -L {regions} -o {sortfile} --write-index {intermediate_bam.name} "
    os.system(command2)


def run_modkit(bamfile, outbed, cpgs, threads, showerrors):
    """
    This function runs modkit on a bam file and extracts the methylation data.
    """
    try:
        command = (
            f"modkit pileup -t {threads} --include-bed {cpgs} --filter-threshold 0.73 --combine-mods {bamfile} "
            f"{outbed} "
        )
        if not showerrors:
            command += "--suppress-progress  >/dev/null 2>&1"
        os.system(command)
        # self.log("Done processing bam file")
    except Exception as e:
        print(e)
        # self.log(e)
        pass


class RandomForest_object(BaseAnalysis):
    def __init__(self, *args, showerrors=False, **kwargs):
        self.rcns2_df_store = {}
        self.threshold = 0.05
        self.bambatch = {}
        self.cpgs_file = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)),
            "hglft_genome_260e9_91a970_clean.bed",
        )
        self.offset = False
        self.first_run = {}
        self.showerrors = showerrors
        self.modelfile = os.path.join(
            os.path.dirname(os.path.abspath(models.__file__)), "general.zip"
        )
        self.dataDir = {}
        self.bedDir = {}
        self.merged_bed_file = {}
        super().__init__(*args, **kwargs)

    def setup_ui(self):
        self.card = ui.card().style("width: 100%")
        with self.card:
            with ui.grid(columns=8).classes("w-full h-auto"):
                with ui.card().classes(
                    f"min-[{self.MENU_BREAKPOINT+1}px]:col-span-3 max-[{self.MENU_BREAKPOINT}px]:col-span-8"
                ):
                    self.create_rcns2_chart("Random Forest")
                with ui.card().classes(
                    f"min-[{self.MENU_BREAKPOINT+1}px]:col-span-5 max-[{self.MENU_BREAKPOINT}px]:col-span-8"
                ):
                    self.create_rcns2_time_chart("Random Forest Time Series")
        if self.summary:
            with self.summary:
                ui.label("Forest classification: Unknown")
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
        try:
            if self.check_file_time(os.path.join(output, "random_forest_scores.csv")):
                self.rcns2_df_store = pd.read_csv(
                    os.path.join(output, "random_forest_scores.csv"),
                    index_col=0,
                )
                columns_greater_than_threshold = (
                    self.rcns2_df_store > self.threshold
                ).any()
                columns_not_greater_than_threshold = ~columns_greater_than_threshold
                result = self.rcns2_df_store.columns[
                    columns_not_greater_than_threshold
                ].tolist()
                self.update_rcns2_time_chart(self.rcns2_df_store.drop(columns=result))
                lastrow = self.rcns2_df_store.iloc[-1]  # .drop("number_probes")
                lastrow_plot = lastrow.sort_values(ascending=False).head(10)
                lastrow_plot_top = lastrow.sort_values(ascending=False).head(1)
                if self.summary:
                    with self.summary:
                        self.summary.clear()
                        ui.label(
                            f"Forest classification: {lastrow_plot_top.index[0]} - {lastrow_plot_top.values[0]:.2f}"
                        )
                self.update_rcns2_plot(
                    lastrow_plot.index.to_list(),
                    list(lastrow_plot.values / 100),
                    "All",
                )
        except FileNotFoundError:
            pass

    async def process_bam(self, bamfile: List[Tuple[str, float]]) -> None:
        """
        Processes the BAM files and performs the NanoDX analysis.

        Args:
            bamfile (List[Tuple[str, float]]): List of BAM files with their timestamps.
        """
        sampleID = self.sampleID
        if sampleID not in self.dataDir.keys():
            self.dataDir[sampleID] = tempfile.TemporaryDirectory(
                dir=self.check_and_create_folder(self.output, sampleID)
            )
        if sampleID not in self.bedDir.keys():
            self.bedDir[sampleID] = tempfile.TemporaryDirectory(
                dir=self.check_and_create_folder(self.output, sampleID)
            )
        tomerge = []
        latest_file = 0
        while len(bamfile) > 0:
            self.running = True
            (file, filetime) = bamfile.pop()
            if filetime > latest_file:
                latest_file = filetime
            tomerge.append(file)
            app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                "bams_in_processing"
            ] += 1
            if len(tomerge) > 50:
                break
        if latest_file:
            currenttime = latest_file * 1000
        else:
            currenttime = time.time() * 1000
        if len(tomerge) > 0:
            tempbam = tempfile.NamedTemporaryFile(
                dir=self.check_and_create_folder(self.output, sampleID),
                suffix=".bam",
            )
            sortbam = tempfile.NamedTemporaryFile(
                dir=self.check_and_create_folder(self.output, sampleID),
                suffix=".bam",
            )
            tempbed = tempfile.NamedTemporaryFile(
                dir=self.check_and_create_folder(self.output, sampleID),
                suffix=".bed",
            )
            self.bambatch[sampleID] = self.bambatch.get(sampleID, 0) + 1

            await run.cpu_bound(
                run_samtools_sort,
                tempbam.name,
                tomerge,
                sortbam.name,
                self.threads,
                self.cpgs_file,
            )

            await run.cpu_bound(
                run_modkit,
                sortbam.name,
                tempbed.name,
                self.cpgs_file,
                self.threads,
                self.showerrors,
            )

            try:
                os.remove(f"{sortbam.name}.csi")
            except FileNotFoundError:
                pass

            if sampleID in self.first_run.keys() and not self.first_run[sampleID]:
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
                    sep="\s+",
                )

                self.merged_bed_file[sampleID] = await run.cpu_bound(
                    merge_bedmethyl, bed_a, self.merged_bed_file.get(sampleID)
                )
                save_bedmethyl(self.merged_bed_file[sampleID], f"{tempbed.name}")
            else:
                self.merged_bed_file[sampleID] = pd.read_table(
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
                    sep="\s+",
                )
                self.first_run[sampleID] = False

            tempDir = tempfile.TemporaryDirectory(
                dir=self.check_and_create_folder(self.output, sampleID)
            )

            await run.cpu_bound(
                run_rcns2,
                tempDir.name,
                self.bambatch[sampleID],
                tempbed.name,
                self.threads,
                self.showerrors,
            )

            if os.path.isfile(
                f"{tempDir.name}/live_{self.bambatch[sampleID]}_votes.tsv"
            ):
                scores = pd.read_table(
                    f"{tempDir.name}/live_{self.bambatch[sampleID]}_votes.tsv",
                    sep="\s+",
                )
                scores_to_save = scores.drop(columns=["Freq"]).T
                scores_to_save["timestamp"] = currenttime

                if sampleID not in self.rcns2_df_store:
                    self.rcns2_df_store[sampleID] = pd.DataFrame()

                self.rcns2_df_store[sampleID] = pd.concat(
                    [
                        self.rcns2_df_store[sampleID],
                        scores_to_save.set_index("timestamp"),
                    ]
                )

                self.rcns2_df_store[sampleID].to_csv(
                    os.path.join(
                        self.check_and_create_folder(self.output, sampleID),
                        "random_forest_scores.csv",
                    )
                )

                app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                    "bam_processed"
                ] += len(tomerge)
                app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                    "bams_in_processing"
                ] -= len(tomerge)

            self.running = False

    def create_rcns2_chart(self, title):
        self.echart = self.create_chart(title)

    def update_rcns2_plot(self, x, y, count):
        """
        Replaces the data in the RapidCNS2 plot.
        :param x: list of tumour types
        :param y: confidence scores for each tumour type
        :param count: the number of bams used to generate the plot
        :return:
        """
        self.echart.options["title"]["text"] = f"Random Forest - processed {count} Bams"
        self.echart.options["yAxis"]["data"] = x
        self.echart.options["series"] = [
            {"type": "bar", "name": "Random Forest", "data": y}
        ]
        self.echart.update()

    def create_rcns2_time_chart(self, title):
        self.rcns2_time_chart = self.create_time_chart(title)

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
        ui.button("start", on_click=start)
        ui.button("stop", on_click=stop)
        TestObject = RandomForest_object(progress=True, batch=True)
    path = "/users/mattloose/datasets/ds1305_Intraop0006_A/20231123_1233_P2S-00770-A_PAS59057_b1e841e7/bam_pass"
    # path = "tests/static/bam"
    directory = os.fsencode(path)
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".bam"):
            TestObject.add_bam(os.path.join(path, filename))
            time.sleep(0.001)
    ui.run(port=8082, reload=False)


def start() -> None:
    yappi.clear_stats()
    yappi.start()


def stop() -> None:
    yappi.stop()
    table = [
        [str(v) for v in [stat.full_name, stat.ttot, stat.tsub, stat.tavg, stat.ncall]]
        for stat in yappi.get_func_stats()
        if "python" not in stat.module
    ]
    print(
        tabulate.tabulate(
            table[:15],
            headers=["function", "total", "excl. sub", "avg", "ncall"],
            floatfmt=".4f",
        )
    )
    yappi.get_thread_stats().print_all()


if __name__ in ("__main__", "__mp_main__"):
    test_ui()
