from __future__ import annotations
from cnsmeth.subpages.base_analysis import BaseAnalysis

from nicegui import ui, app, run
import time
import os
import sys
import click
from pathlib import Path
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


def run_modkit(cpgs, sortfile, temp, threads):
    """
    This function runs modkit on a bam file and extracts the methylation data.
    """
    try:
        os.system(
            f"modkit pileup --include-bed {cpgs} --filter-threshold 0.73 --combine-mods --only-tabs -t {threads} {sortfile} {temp} --suppress-progress >/dev/null 2>&1"
        )
    except Exception as e:
        print(e)


def run_samtools_sort(file, tomerge, sortfile, threads):
    pysam.cat("-o", file, *tomerge)
    pysam.sort("-@", f"{threads}", "--write-index", "-o", sortfile, file)


def classification(modelfile, test_df):
    NN = NN_classifier(modelfile)
    try:
        predictions, class_labels, n_features = NN.predict(test_df)
    except Exception as e:
        print(e)
        test_df.to_csv("errordf.csv", sep=",", index=False, encoding="utf-8")
        #sys.exit(1)
    return predictions, class_labels, n_features


class NanoDX_object(BaseAnalysis):
    def __init__(self, *args, model="Capper_et_al_NN.pkl", **kwargs):
        self.cpgs_file = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)),
            "hglft_genome_260e9_91a970_clean.bed",
        )
        self.cpgs = pd.read_csv(
            self.cpgs_file,
            sep="\t",
            header=None,
        )
        self.model = model
        self.threshold = 0.05
        self.nanodx_bam_count = 0
        self.not_first_run = False
        self.modelfile = os.path.join(
            os.path.dirname(os.path.abspath(models.__file__)), self.model
        )
        self.nanodx_df_store = pd.DataFrame()
        #self.NN = NN_classifier(self.modelfile)
        super().__init__(*args, **kwargs)
        self.nanodxfile = tempfile.NamedTemporaryFile(dir=self.output, suffix=".nanodx")

    def setup_ui(self):
        with ui.card().style("width: 100%"):
            with ui.grid(columns=8).classes("w-full h-auto"):
                with ui.card().classes("col-span-3"):
                    self.create_nanodx_chart("NanoDX")
                with ui.card().classes("col-span-5"):
                    self.create_nanodx_time_chart("NanoDX Time Series")
        if self.summary:
            with self.summary:
                ui.label("NanoDX classification: Unknown")
        if self.browse:
            self.show_previous_data(self.output)
        else:
            ui.timer(5, lambda: self.show_previous_data(self.output))

    def show_previous_data(self, output):
        if self.check_file_time(os.path.join(output, "nanoDX_scores.csv")):
            self.nanodx_df_store = pd.read_csv(
                os.path.join(os.path.join(self.output, "nanoDX_scores.csv")),
                index_col=0,
            )
            columns_greater_than_threshold = (
                self.nanodx_df_store > self.threshold
            ).any()
            columns_not_greater_than_threshold = ~columns_greater_than_threshold
            result = self.nanodx_df_store.columns[
                columns_not_greater_than_threshold
            ].tolist()
            self.update_nanodx_time_chart(self.nanodx_df_store.drop(columns=result))
            lastrow = self.nanodx_df_store.iloc[-1].drop("number_probes")
            lastrow_plot = lastrow.sort_values(ascending=False).head(10)
            lastrow_plot_top = lastrow.sort_values(ascending=False).head(1)
            if self.summary:
                with self.summary:
                    self.summary.clear()
                    ui.label(
                        f"NanoDX classification: {lastrow_plot_top.index[0]} - {lastrow_plot_top.values[0]:.2f}"
                    )
            self.update_nanodx_plot(
                lastrow_plot.index.to_list(),
                list(lastrow_plot.values),
                "All",
                self.nanodx_df_store.iloc[-1]["number_probes"],
            )

    async def process_bam(self, bamfile):
        tomerge = []
        # timestamp = None
        while len(bamfile) > 0:
            self.running = True
            (file, filetime) = bamfile.pop()
            self.nanodx_bam_count += 1
            tomerge.append(file)
            # timestamp = filetime

            if len(tomerge) > 5:
                break
        app.storage.general[self.mainuuid][self.name]["counters"][
            "bams_in_processing"
        ] += len(tomerge)

        if len(tomerge) > 0:
            tempbam = tempfile.NamedTemporaryFile(dir=self.output, suffix=".bam")
            sorttempbam = tempfile.NamedTemporaryFile(dir=self.output, suffix=".bam")
            file = tempbam.name

            temp = tempfile.NamedTemporaryFile(dir=self.output)

            sortfile = sorttempbam.name

            await run.cpu_bound(
                run_samtools_sort, file, tomerge, sortfile, self.threads
            )

            await run.cpu_bound(
                run_modkit, self.cpgs_file, sortfile, temp.name, self.threads
            )

            try:
                os.remove(f"{sortfile}.csi")
            except FileNotFoundError:
                pass

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
                    sep="\s+",
                    # delim_whitespace=True,
                )
                self.merged_bed_file = await run.cpu_bound(
                    merge_bedmethyl, bed_a, self.merged_bed_file
                )
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
                    sep="\s+",
                    # delim_whitespace=True,
                )

                self.not_first_run = True
            self.merged_bed_file = await run.cpu_bound(
                collapse_bedmethyl, self.merged_bed_file
            )
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

            predictions, class_labels, n_features = await run.cpu_bound(
                classification, self.modelfile, test_df
            )

            # try:
            #    predictions, class_labels, n_features = self.NN.predict(test_df)
            # except Exception as e:
            #    print(e)
            #    test_df.to_csv("errordf.csv", sep=",", index=False, encoding="utf-8")
            #    # self.nanodx_status_txt["message"] = "Error generating predictions."
            #    sys.exit(1)
            nanoDX_df = pd.DataFrame({"class": class_labels, "score": predictions})

            nanoDX_save = nanoDX_df.set_index("class").T
            nanoDX_save["number_probes"] = n_features
            # if timestamp:
            #    nanoDX_save["timestamp"] = timestamp * 1000
            # else:
            nanoDX_save["timestamp"] = time.time() * 1000

            self.nanodx_df_store = pd.concat(
                [self.nanodx_df_store, nanoDX_save.set_index("timestamp")]
            )

            self.nanodx_df_store.to_csv(os.path.join(self.output, "nanoDX_scores.csv"))

            app.storage.general[self.mainuuid][self.name]["counters"][
                "bam_processed"
            ] += len(tomerge)
            app.storage.general[self.mainuuid][self.name]["counters"][
                "bams_in_processing"
            ] -= len(tomerge)
        self.running = False

    def create_nanodx_chart(self, title):
        self.nanodxchart = self.create_chart(title)

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

    def create_nanodx_time_chart(self, title):
        self.nanodx_time_chart = self.create_time_chart(title)

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


def test_me(
    port: int,
    threads: int,
    watchfolder: str,
    output: str,
    reload: bool = False,
    browse: bool = False,
):
    # if __name__ == '__mp_main__':
    my_connection = None
    with theme.frame("Target Coverage Data", my_connection):
        TestObject = NanoDX_object(threads, output, progress=True, batch=True)
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
