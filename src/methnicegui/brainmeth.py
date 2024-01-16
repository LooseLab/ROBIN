import click
from nicegui import Tailwind, ui, app
from nicegui.events import ValueChangeEventArguments
from io import StringIO

from random import random
from pathlib import Path

from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer
import logging
from methnicegui import images
import threading
import time
import os, signal
import pysam
import pandas as pd
import shutil
import subprocess
import natsort
import numpy as np
import tempfile
from cnv_from_bam import iterate_bam_file

from methnicegui import models, resources

from methnicegui.merge_bedmethyl import merge_bedmethyl, save_bedmethyl

from sturgeon.callmapping import (
    merge_probes_methyl_calls,
    probes_methyl_calls_to_bed,
)

os.environ["CI"] = "1"


# Very basic example GUI in a class.
class BrainMeth:
    def __init__(
        self,
        threads=4,
        simtime=False,
        watchfolder=None,
        output=None,
        sequencing_summary=None,
        showerrors=False,
    ):
        self.threshold = 0.001
        self.runfinished = False
        self.rcns2finished = False
        self.stugeonfinished = False
        self.nofiles = False
        self.showerrors=showerrors
        self.bamfolder = watchfolder
        self.threads = threads
        self.imagefile = os.path.join(
            os.path.dirname(os.path.abspath(images.__file__)), "MethBrain_small.png"
        )
        self.bam_count = {"counter": 0}
        self.rapidcns_status_txt = {"message": "Waiting for data."}
        self.sturgeon_status_txt = {"message": "Waiting for data."}
        self.rcns2_bam_count = 0
        self.st_bam_count = 0
        self.st_num_probes = 0
        self.simtime = simtime
        self.cnv_dict = {}
        self.cnv_dict["bin_width"] = 0
        self.cnv_dict["variance"] = 0
        self.bamforcns = []
        self.bamforsturgeon = []
        self.outputfolder = output
        self.bedfoldercount = os.path.join(self.outputfolder, "bedscount")
        if not os.path.exists(self.bedfoldercount):
            os.makedirs(self.bedfoldercount)
        self.resultfolder = os.path.join(self.outputfolder, "results")
        if not os.path.exists(self.resultfolder):
            os.mkdir(self.resultfolder)
        self.subsetbamfolder = os.path.join(self.outputfolder, "subsetbams")
        if not os.path.exists(self.subsetbamfolder):
            os.mkdir(self.subsetbamfolder)
        self.donebamfolder = os.path.join(self.outputfolder, "donebams")
        if not os.path.exists(self.donebamfolder):
            os.mkdir(self.donebamfolder)
        self.targetsbamfolder = os.path.join(self.outputfolder, "targetsbams")
        if not os.path.exists(self.targetsbamfolder):
            os.mkdir(self.targetsbamfolder)
        self.sortedbamfile = os.path.join(self.subsetbamfolder, "sorted.bam")
        self.merged_bam_file = None
        self.modelfile = os.path.join(
            os.path.dirname(os.path.abspath(models.__file__)), "general.zip"
        )
        self.rcns2folder = os.path.join(self.outputfolder, "rcns2")
        if not os.path.exists(self.rcns2folder):
            os.mkdir(self.rcns2folder)
        self.result = None
        self.sturgeon_df_store = pd.DataFrame()
        self.rcns2_df_store = pd.DataFrame()
        self.target_coverage_df = pd.DataFrame()
        self.gene_bed = pd.read_table(
            os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)), "unique_genes.bed"
            ),
            names=["chrom", "start_pos", "end_pos", "gene"],
            header=None,
            delim_whitespace=True,
        )

        self.event_handler = BamEventHandler(self.bam_count, logger=None)
        self.observer = Observer()
        self.observer.schedule(self.event_handler, self.bamfolder, recursive=True)
        self.observer.start()
        self.bam_processing = threading.Thread(target=self.process_bams, args=())
        self.bam_processing.daemon = True
        self.bam_processing.start()
        self.rapidcns2_processing = threading.Thread(target=self.rapid_cns2, args=())
        self.rapidcns2_processing.daemon = True
        self.rapidcns2_processing.start()
        self.sturgeon_processing = threading.Thread(target=self.sturgeon, args=())
        self.sturgeon_processing.daemon = True
        self.sturgeon_processing.start()
        self.sequence_summary = sequencing_summary
        self.check_for_existing_bams = threading.Thread(
            target=self.check_existing_bams,
            args=(),
            kwargs={"sequencing_summary": self.sequence_summary},
        )
        self.check_for_existing_bams.daemon = True
        self.check_for_existing_bams.start()

        with ui.dialog() as dialog, ui.card():
            ui.label("Quitting the app will stop monitoring the code. Are you sure?")
            ui.button("Cancel", on_click=dialog.close).props("outline").classes(
                "shadow-lg"
            )
            ui.button(
                "Really Quit", icon="logout", on_click=self.cleanup_and_exit
            ).props("outline").classes("shadow-lg")

        with ui.header(fixed=True).classes(replace="row items-center p-2") as header:
            with ui.grid(columns=2).style("width: 100%"):
                ui.label("Real Time Brain Tumour Nanopore Classifier").tailwind(
                    "text-2xl font-bold font-italic drop-shadow"
                )
                with ui.row().classes("ml-auto"):
                    self.dark = ui.dark_mode()
                    ui.switch(
                        "Dark Mode", on_change=self.dark_mode
                    ).classes('ml-4 bg-transparent').props('color="black"')  # .classes('ml-4')#.props('outline')
                    ui.button(
                        "Quit", icon="logout", on_click=dialog.open
                    )  # .classes('ml-4')#.props('outline') #.classes('shadow-lg')
                    ui.image(self.imagefile).style("width: 40px")

        with ui.card().style("width: 100%"):
            ui.label(
                "This tool enables classification of brain tumours in real time from Oxford Nanopore Data."
            ).tailwind("drop-shadow", "font-bold")
            with ui.row():
                ui.label(f"Monitoring the path:{watchfolder}").tailwind("drop-shadow")
                ui.label(f"Outputting to:{output}").tailwind("drop-shadow")
                ui.label().bind_text_from(
                    self.bam_count, "counter", backward=lambda n: f"BAM files seen: {n}"
                ).tailwind("drop-shadow")

        with ui.tabs().classes('w-full') as tabs:
            methylation = ui.tab('Methylation Classification')
            copy_numer = ui.tab('Copy Number Variation')
            coverage = ui.tab('Target Coverage')
            mgmt = ui.tab('MGMT')
        with ui.tab_panels(tabs, value=methylation).classes('w-full'):
            with ui.tab_panel(methylation):
                with ui.card().style("width: 100%"):
                    with ui.grid(columns=2).classes("w-full h-auto"):
                        with ui.column():
                            with ui.card().style("width: 100%"):
                                ui.label().bind_text_from(
                                    self.rapidcns_status_txt,
                                    "message",
                                    backward=lambda n: f"RapidCNS2 Status: {n}",
                                )
                                self.create_chart("RapidCNS2")
                        with ui.column():
                            with ui.card().style("width: 100%"):
                                ui.label().bind_text_from(
                                    self.sturgeon_status_txt,
                                    "message",
                                    backward=lambda n: f"Sturgeon Status: {n}",
                                )
                                self.create_chart2("Sturgeon")

                with ui.card().style("width: 100%"):
                    self.create_sturgeon_time_chart()

                with ui.card().style("width: 100%"):
                    self.create_rcns2_time_chart()

            with ui.tab_panel(copy_numer):
                with ui.card().style("width: 100%"):
                    with ui.row():
                        self.chrom_select = ui.select(
                            options={"All": "All"},
                            on_change=self.update_cnv_plot,
                            label="Select Chromosome",
                            value="All",
                        ).style("width: 150px")
                        self.gene_select = ui.select(
                            options={"All": "All"},
                            on_change=lambda e: self.update_cnv_plot()
                            if e.value == "All"
                            else self.update_cnv_plot(gene_target=e.value),
                            label="Select Gene",
                            value="All",
                        ).style("width: 150px")
                        ui.label().bind_text_from(
                            self.cnv_dict, "bin_width", backward=lambda n: f"Bin Width: {n}"
                        )
                        ui.label().bind_text_from(
                            self.cnv_dict, "variance", backward=lambda n: f"Variance: {n}"
                        )
                    self.create_scatter("Copy Number Variation")

            with ui.tab_panel(coverage):
                with ui.card().style("width: 100%"):
                    ui.label("Coverage Data").tailwind("drop-shadow", "font-bold")
                    with ui.grid(columns=2).classes("w-full h-auto"):
                        with ui.column():
                            with ui.card().style("width: 100%"):
                                self.create_coverage_plot("Chromosome Coverage")
                        with ui.column():
                            with ui.card().style("width: 100%"):
                                self.create_coverage_plot_targets("Target Coverage")
                with ui.card().style("width: 100%"):
                    self.create_coverage_time_chart()
                with ui.card().style("width: 100%"):
                    self.f = ui.input('Filter')
                    self.targ_df = ui.row()

            self.mgmt = ui.tab_panel(mgmt)
            with self.mgmt:
                with ui.card().style("width: 100%"):
                    ui.label("MGMT Methylation").tailwind("drop-shadow", "font-bold")
                    #self.mgmtcontent=ui.column()
                    #with self.mgmtcontent:
                    self.mgmtplot = ui.row()
                    with self.mgmtplot.classes('w-full'):
                        ui.label("Plot not yet available.")
                        #ui.image("/tmp/run2/targetsbams/25_sorted.png").props("fit=scale-down")
                    self.mgmtable = ui.row()
                    with self.mgmtable:
                        ui.label("Table not yet available.")






        with ui.footer():
            ui.label(
                "Some aspects of this application are ©Looselab - all analyses provided for research use only."
            ).tailwind("text-sm font-bold font-italic drop-shadow")

    def cleanup_and_exit(self):
        self.observer.stop()
        self.observer.join()
        app.shutdown()

    def process_bams(self) -> None:
        """
        This function processes the bam files and adds them to the bamforcns and bamforsturgeon lists.
        These lists are then processed by the rapid_cns2 and sturgeon functions.
        #ToDo: switch to using threadsafe queues.
        :param self:
        :return:
        """

        # First we look for existing files:
        # self.check_existing_bams()

        while True:
            # self.log("Processing bam files")
            # self.log(self.bam_count)
            # self.bam_count = self.bam_count
            if "file" in self.bam_count:
                while len(self.bam_count["file"]) > 0:
                    self.nofiles = False
                    file = self.bam_count["file"].popitem()
                    if file[1] > time.time() - 5:
                        time.sleep(5)
                    self.bamforcns.append(file[0])
                    self.bamforsturgeon.append(file[0])
                time.sleep(1)
                self.nofiles = True
            time.sleep(1)

    def check_existing_bams(self, sequencing_summary=None):
        file_endings = {".bam"}
        if sequencing_summary:
            # get the current time
            now = time.time()
            print("Checking for existing bams against sequencing summary")
            df = pd.read_csv(
                sequencing_summary,
                delimiter="\t",
                usecols=["filename_bam", "template_start", "template_duration"],
            )
            df["template_end"] = df["template_start"] + df["template_duration"]

            df.drop(columns=["template_start", "template_duration"], inplace=True)
            latest_timestamps = (
                df.groupby("filename_bam")["template_end"]
                .max()
                .reset_index()
                .sort_values(by="template_end")
                .reset_index(drop=True)
            )
            latest_timestamps["full_path"] = ""
            latest_timestamps["file_produced"] = latest_timestamps["template_end"] + now
            print(latest_timestamps)
            for path, dirs, files in os.walk(self.bamfolder):
                for f in files:
                    if "".join(Path(f).suffixes) in file_endings:
                        latest_timestamps.loc[
                            latest_timestamps["filename_bam"] == f, "full_path"
                        ] = os.path.join(path, f)
                        # latest_timestamps[latest_timestamps['filename_bam'] == f]['full_path'] = os.path.join(path, f)
            print(latest_timestamps)
            for index, row in latest_timestamps.iterrows():
                current_time = time.time()
                time_diff = row["file_produced"] - current_time
                print(f"Sleeping for {time_diff} seconds.")
                if time_diff > 0:
                    time.sleep(time_diff)
                self.bam_count["counter"] += 1
                if "file" not in self.bam_count:
                    self.bam_count["file"] = {}
                self.bam_count["file"][row["full_path"]] = time.time()
                # print (index,row)
            self.runfinished = True
            #os.kill(os.getpid(), signal.SIGINT)
        else:
            # print ("Checking Bams")
            for path, dirs, files in os.walk(self.bamfolder):
                for f in files:
                    if "".join(Path(f).suffixes) in file_endings:
                        self.bam_count["counter"] += 1
                        if "file" not in self.bam_count:
                            self.bam_count["file"] = {}
                        self.bam_count["file"][os.path.join(path, f)] = time.time()
                        if self.simtime:
                            time.sleep(1)

    def rapid_cns2(self):
        """
        This function runs RapidCNS2 on the bam files.
        It runs as a worker to not block the main thread.
        It grabs batches of bam files from the bamforcns list once the list contains 5 or more BAM files.
        :return:
        """
        cov_df = None  # This will be a pandas dataframe to store coverage information.
        not_first_run = False
        batch = 0
        wait_for_batch_size = 2
        while True:
            if len(self.bamforcns) > wait_for_batch_size:
                self.rcns2finished = False
                wait_for_batch_size = 0
                batch += 1
                print("Running RapidCNS")
                bams = []
                print(f"there are {len(self.bamforcns)} bams for processing")
                self.rapidcns_status_txt[
                    "message"
                ] = f"Processing {len(self.bamforcns)} bam files for RapidCNS2."
                while len(self.bamforcns) > 0:
                    bams.append(self.bamforcns.pop())
                    self.rcns2_bam_count += 1

                # self.log.info(bams)

                self.rapidcnsbamfile = os.path.join(
                    self.subsetbamfolder, f"{batch}_rapidcns.bam"
                )

                self.previous_rapidcnsbamfile = os.path.join(
                    self.subsetbamfolder, f"{batch-1}_rapidcns.bam"
                )

                tempbam = tempfile.NamedTemporaryFile()

                self.rapidcns_status_txt["message"] = "Merging bam files."

                #The bam file created here consists of all the reads from the bam files in the most recent batch.

                pysam.cat("-o", tempbam.name, *bams)

                self.rapidcns_status_txt["message"] = "Sorting bam files."

                os.system(
                    f"samtools sort --write-index -@{self.threads} -o {self.sortedbamfile} {tempbam.name} "
                    f">/dev/null 2>&1"
                )

                self.rapidcns_status_txt["message"] = "Running modkit pileup."

                os.system(
                    f"modkit pileup -t {self.threads} --filter-threshold 0.73 --combine-mods {self.sortedbamfile} "
                    f"{self.rapidcnsbamfile}.bed --suppress-progress  >/dev/null 2>&1 "
                )  #

                if not_first_run:
                    self.rapidcns_status_txt[
                        "message"
                    ] = "Merging bed file with previous bed files."
                    bed_a = pd.read_table(
                        f"{self.rapidcnsbamfile}.bed",
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
                    # Doing this is slow - we should just hold the file in memory surely?
                    # bedB = pd.read_table(f"{self.previous_rapidcnsbamfile}.bed",
                    # names = ['chrom', 'start_pos', 'end_pos', 'mod', 'score', 'strand',
                    # 'start_pos2', 'end_pos2', 'colour','Nvalid', 'fraction', 'Nmod', 'Ncanon',
                    # 'Nother', 'Ndel', 'Nfail', 'Ndiff', 'Nnocall'], header=None, delim_whitespace=True)
                    self.merged_bed_file = merge_bedmethyl(bed_a, self.merged_bed_file)
                    save_bedmethyl(self.merged_bed_file, f"{self.rapidcnsbamfile}.bed")
                    newcovdf = pd.read_csv(
                        StringIO(pysam.coverage(f"{self.sortedbamfile}")), sep="\t"
                    )
                    newcovdf.drop(
                        columns=["coverage", "meanbaseq", "meanmapq"],
                        inplace=True,
                    )
                    merged_df = pd.merge(
                        covdf,
                        newcovdf,
                        on=["#rname", "startpos", "endpos"],
                        suffixes=("_df1", "_df2"),
                    )
                    merged_df["numreads"] = (
                        merged_df["numreads_df1"] + merged_df["numreads_df2"]
                    )
                    merged_df["covbases"] = (
                            merged_df["covbases_df1"] + merged_df["covbases_df2"]
                    )
                    merged_df["meandepth"] = (
                        merged_df["meandepth_df1"] + merged_df["meandepth_df2"]
                    )

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
                    covdf = merged_df

                    newbedcovdf = pd.read_csv(
                        StringIO(
                            pysam.bedcov(
                                os.path.join(
                                    os.path.dirname(
                                        os.path.abspath(resources.__file__)
                                    ),
                                    "unique_genes.bed",
                                ),
                                f"{self.sortedbamfile}",
                            )
                        ),
                        names=["chrom", "startpos", "endpos", "name", "bases"],
                        sep="\t",
                    )
                    merged_df = pd.merge(
                        bedcovdf,
                        newbedcovdf,
                        on=["chrom", "startpos", "endpos", "name"],
                        suffixes=("_df1", "_df2"),
                    )
                    merged_df["bases"] = merged_df["bases_df1"] + merged_df["bases_df2"]
                    merged_df.drop(columns=["bases_df1", "bases_df2"], inplace=True)
                    bedcovdf = merged_df
                else:
                    self.merged_bed_file = pd.read_table(
                        f"{self.rapidcnsbamfile}.bed",
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
                    covdf = pd.read_csv(
                        StringIO(pysam.coverage(f"{self.sortedbamfile}")), sep="\t"
                    )
                    covdf.drop(
                        columns=["coverage", "meanbaseq", "meanmapq"],
                        inplace=True,
                    )
                    bedcovdf = pd.read_csv(
                        StringIO(
                            pysam.bedcov(
                                os.path.join(
                                    os.path.dirname(
                                        os.path.abspath(resources.__file__)
                                    ),
                                    "unique_genes.bed",
                                ),
                                f"{self.sortedbamfile}",
                            )
                        ),
                        names=["chrom", "startpos", "endpos", "name", "bases"],
                        sep="\t",
                    )
                    not_first_run = True
                self.update_coverage_plot(covdf)
                self.update_coverage_plot_targets(covdf, bedcovdf)
                self.update_coverage_time_plot(covdf)
                self.target_coverage_df=bedcovdf
                self.target_coverage_df['coverage']=self.target_coverage_df['bases']/self.target_coverage_df['length']
                with self.targ_df:
                    self.targ_df.clear()
                    ui.table.from_pandas(self.target_coverage_df, pagination=10, row_key='name').bind_filter_from(self.f, 'value').classes('w-full')
                # self.log("Merged Bed File Info:")
                # self.log(self.merged_bed_file.info())

                self.rapidcns_status_txt[
                    "message"
                ] = "Running RapidCNS2 methylation classification."

                command = f"Rscript ../hv_rapidCNS2/bin/methylation_classification_nanodx_v0.1.R -s " +\
                    f"live_{batch} -o {self.rcns2folder} -i {self.rapidcnsbamfile}.bed " +\
                    f"-p ../hv_rapidCNS2/bin/top_probes_hm450.Rdata " +\
                    f"--training_data ../hv_rapidCNS2/bin/capper_top_100k_betas_binarised.Rdata " +\
                    f"--array_file ../hv_rapidCNS2/bin/HM450.hg38.manifest.gencode.v22.Rdata " +\
                    f"-t {self.threads} "

                if not self.showerrors:
                    command += ">/dev/null 2>&1"

                returned_value = subprocess.call(command,
                    shell=True,
                )
                if self.showerrors:
                    print(command)
                    print(returned_value)
                # self.log(returned_value)

                scores = pd.read_table(
                    f"{self.rcns2folder}/live_{batch}_votes.tsv", delim_whitespace=True
                )
                scores_to_save = scores.drop(columns=["Freq"]).T
                scores_to_save["timestamp"] = time.time() * 1000
                self.rcns2_df_store = pd.concat(
                    [self.rcns2_df_store, scores_to_save.set_index("timestamp")]
                )
                self.rcns2_df_store.to_csv(
                    os.path.join(self.resultfolder, "rcns2_scores.csv")
                )

                columns_greater_than_threshold = (self.rcns2_df_store > self.threshold*100).any()
                columns_not_greater_than_threshold = ~columns_greater_than_threshold
                result = self.rcns2_df_store.columns[columns_not_greater_than_threshold].tolist()


                self.update_rcns2_time_chart(self.rcns2_df_store.drop(columns=result))
                # with self.rcns2_container:
                #    self.rcns2_container.clear()
                #    ui.table.from_pandas(self.rcns2_df_store, pagination=3).classes('max-h-80')
                scores = scores.sort_values(by=["cal_Freq"], ascending=False).head(10)
                # print(scores.index.to_list())
                # print(list(scores["cal_Freq"].values / 100))
                self.update_rcns2_plot(
                    scores.index.to_list(),
                    list(scores["cal_Freq"].values / 100),
                    self.rcns2_bam_count,
                )

                self.rapidcns_status_txt[
                    "message"
                ] = "RapidCNS2 methylation classification done. Waiting for data."


                os.rename(
                    f"{self.sortedbamfile}",
                    os.path.join(self.donebamfolder, f"{batch}_sorted.bam"),
                )
                os.rename(
                    f"{self.sortedbamfile}.csi",
                    os.path.join(self.donebamfolder, f"{batch}_sorted.bam.csi"),
                )
                self.cnv_plotting()
                self.keep_regions(os.path.join(self.donebamfolder, f"{batch}_sorted.bam"), batch)
                self.mgmtmethylpredict(self.rapidcnsbamfile)

                pass
            time.sleep(5)

            if len(self.bamforcns)==0:
                self.rcns2finished = True
                if self.nofiles and self.stugeonfinished and self.runfinished and self.rcns2finished:
                    print ("All done")
                    os.kill(os.getpid(), signal.SIGINT)

    def mgmtmethylpredict(self, bamfile):
        print("Running MGMT predictor")
        self.rapidcns_status_txt["message"] = "Running MGMT predictor."
        MGMT_BED = "../hv_rapidCNS2/bin/mgmt_hg38.bed"
        os.system(
            f"bedtools intersect -a {bamfile}.bed -b {MGMT_BED} > {self.resultfolder}/mgmt_result.bed"
        )
        print(f"bedtools intersect -a {bamfile}.bed -b {MGMT_BED} > {self.resultfolder}/mgmt_result.bed")

        if os.path.getsize(f"{self.resultfolder}/mgmt_result.bed") > 0:

            cmd = f"Rscript ../hv_rapidCNS2/bin/mgmt_pred_v0.3.R --input={self.resultfolder}/mgmt_result.bed --out_dir={self.resultfolder} --probes=../hv_rapidCNS2/bin/mgmt_probes.Rdata --model=../hv_rapidCNS2/bin/mgmt_137sites_mean_model.Rdata --sample=live_analysis"
            os.system(
                cmd
            )
            results = pd.read_csv(os.path.join(self.resultfolder,"live_analysis_mgmt_status.csv"))
            self.mgmtable.clear()
            with self.mgmtable:
                ui.table.from_pandas(results)
            print(cmd)
            print ("MGMT predictor done")
            self.rapidcns_status_txt["message"] = "MGMT predictor done."

        else:
            print ("No MGMT sites yet found.")
            self.rapidcns_status_txt["message"] = "No MGMT sites yet found."




    def keep_regions(self, bamtoextract, batch):
        print("Keeping regions")
        bam_out = os.path.join(self.targetsbamfolder, f"{batch}_sorted.bam")
        plot_out = os.path.join(self.targetsbamfolder, f"{batch}_sorted.png")
        bedfile = os.path.join(os.path.dirname(os.path.abspath(resources.__file__)),'unique_genes.bed')
        os.system(
            f"samtools view --write-index -L {bedfile} -@{self.threads} -o {bam_out} {bamtoextract} "
            #f">/dev/null 2>&1"
        )
        merged_bam_out = os.path.join(self.targetsbamfolder, f"{batch}_merged.bam")
        to_be_merged = os.path.join(self.targetsbamfolder, f"*_sorted.bam")
        os.system(
            f"samtools merge --write-index -@{self.threads} -f {merged_bam_out} {to_be_merged}"
            #f">/dev/null 2>&1"
        )
        os.system(
            f"rm {to_be_merged}"
        )
        os.system(
            f"mv {merged_bam_out} {bam_out}"
        )
        os.system(
            f"mv {merged_bam_out}.csi {bam_out}.csi"
        )
        print(f"methylartist locus -i chr10:129466536-129467536 -b {bam_out} -o {plot_out}  --motif CG --mods m")
        os.system(
            f"methylartist locus -i chr10:129466536-129467536 -b {bam_out} -o {plot_out}  --motif CG --mods m"
        )
        if os.path.exists(plot_out):
            self.mgmtplot.clear()
            with self.mgmtplot.classes('w-full'):
                ui.image(plot_out).props('fit=scale-down')


    def cnv_plotting(self):
        bam_path = Path(self.donebamfolder)
        self.result = iterate_bam_file(
            bam_path, _threads=4, mapq_filter=60, log_level=logging.getLevelName("WARN")
        )
        self.cnv_dict["bin_width"] = self.result.bin_width
        self.cnv_dict["variance"] = self.result.variance
        # print(self.result)
        self.update_cnv_plot()

    def update_cnv_plot(self, gene_target=None):
        if self.result:
            # print (self.scatter_echart.options)
            total = 0
            y = 5
            # self.cnv_plot.clear_data()
            valueslist = {"All": "All"}
            genevalueslist = {"All": "All"}
            self.chrom_filter = self.chrom_select.value

            min = 0
            max = "dataMax"

            if gene_target:
                print ("Gene Target")
                start_pos = self.gene_bed.iloc[int(gene_target)].start_pos
                end_pos = self.gene_bed.iloc[int(gene_target)].end_pos
                chrom = self.gene_bed.iloc[int(gene_target)].chrom
                print(chrom)
                for counter, contig in enumerate(
                    natsort.natsorted(self.result.cnv), start=1
                ):
                    valueslist[counter] = contig
                self.chrom_filter = counter
                print(self.chrom_filter)
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

                    # self.cnv_plot.update_plot(x=np.arange(len(cnv)) * self.result.bin_width, y=cnv,
                    #                          range=[gene_target[0] - 10*self.bin_width,
                    #                                 gene_target[1] + 10*self.bin_width])
                # self.cnv_plot.add_text(contig, total, y=y)
                # y = 7
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

                    # self.cnv_plot.add_vline(row["start_pos"], color="green+")
                    # self.cnv_plot.add_text(row["gene"], row["start_pos"], y=y)
            # print (self.scatter_echart.options)
            self.chrom_select.set_options(valueslist)
            # self.chrom_select.update()
            self.gene_select.set_options(genevalueslist)
            # self.gene_select.update()
            self.scatter_echart.update()

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
            if len(self.bamforsturgeon) > 0:
                self.stugeonfinished = False
                bams = []
                # self.log(f"there are {len(self.bamforsturgeon)} bams for sturgeon processing")
                self.sturgeon_status_txt[
                    "message"
                ] = f"Processing {len(self.bamforsturgeon)} bam files for Sturgeon."
                while len(self.bamforsturgeon) > 0:
                    bams.append(self.bamforsturgeon.pop())
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

                    columns_greater_than_threshold = (self.sturgeon_df_store > self.threshold).any()
                    columns_not_greater_than_threshold = ~columns_greater_than_threshold
                    result = self.sturgeon_df_store.columns[columns_not_greater_than_threshold].tolist()

                    self.update_sturgeon_time_chart(self.sturgeon_df_store.drop(columns=result))

                    self.update_sturgeon_plot(
                        lastrow_plot.index.to_list(),
                        list(lastrow_plot.values),
                        self.st_bam_count,
                    )

                    self.sturgeon_status_txt[
                        "message"
                    ] = "Predictions generated. Waiting for data."

            time.sleep(5)

            if len(self.bamforsturgeon) == 0:
                self.stugeonfinished = True

    def dark_mode(self, event: ValueChangeEventArguments):
        if event.value:
            self.dark.enable()
        else:
            self.dark.disable()

    def create_sturgeon_time_chart(self):
        self.sturgeon_time_chart = (
            ui.echart(
                {
                    "grid": {"containLabel": True},
                    "title": {"text": "Sturgeon Over Time"},
                    'toolbox': {
                        'show': True,
                        'feature': {
                            'saveAsImage': {}
                        }
                    },
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

    def create_rcns2_time_chart(self):
        self.rcns2_time_chart = (
            ui.echart(
                {
                    "grid": {"containLabel": True},
                    "title": {"text": "RCNS2 Over Time"},
                    'toolbox': {
                        'show': True,
                        'feature': {
                            'saveAsImage': {}
                        }
                    },
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
            # print(series)
            data_list = [[key, value] for key, value in data.items()]
            # print(data_list)
            if series != "number_probes":
                self.rcns2_time_chart.options["series"].append(
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
        self.rcns2_time_chart.update()

    def create_scatter(self, title):
        self.scatter_echart = (
            ui.echart(
                {
                    "grid": {"containLabel": True},
                    "title": {"text": title},
                    'toolbox': {
                        'show': True,
                        'feature': {
                            'saveAsImage': {}
                        }
                    },
                    "xAxis": {"type": "value", "max": "dataMax"},
                    #'yAxis': {'axisLabel': {':formatter': 'value => "Ploidy" + value'}},
                    "yAxis": {"type": "value"},
                    #'legend': {'type':'scroll','formatter': '{name}','top': 'bottom'},
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
                            "data": [],
                        }
                    ],
                }
            )
            .style("height: 350px")
            .classes("border-double")
        )

    def create_chart(self, title):
        self.echart = (
            ui.echart(
                {
                    "grid": {"containLabel": True},
                    "title": {"text": title},
                    'toolbox': {
                        'show': True,
                        'feature': {
                            'saveAsImage': {}
                        }
                    },
                    "xAxis": {"type": "value", "max": 1},
                    "yAxis": {"type": "category", "data": [], "inverse": True},
                    #'legend': {},
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

    def create_coverage_plot_targets(self, title):
        self.echart4 = (
            ui.echart(
                {
                    "grid": {"containLabel": True},
                    "title": {"text": title},
                    'toolbox': {
                        'show': True,
                        'feature': {
                            'saveAsImage': {}
                        }
                    },
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
        self.echart4.options["title"]["text"] = f"Per Chromosome Target Coverage"
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

    def create_coverage_plot(self, title):
        self.echart3 = (
            ui.echart(
                {
                    "grid": {"containLabel": True},
                    "title": {"text": title},
                    'toolbox': {
                        'show': True,
                        'feature': {
                            'saveAsImage': {}
                        }
                    },
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
        self.echart3.options["title"]["text"] = f"Per Chromosome Coverage"
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

    def create_chart2(self, title):
        self.echart2 = (
            ui.echart(
                {
                    "grid": {"containLabel": True},
                    "title": {"text": title},
                    'toolbox': {
                        'show': True,
                        'feature': {
                            'saveAsImage': {}
                        }
                    },
                    "xAxis": {"type": "value", "max": 1},
                    "yAxis": {"type": "category", "data": [], "inverse": True},
                    #'legend': {},
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
                    'toolbox': {
                        'show': True,
                        'feature': {
                            'saveAsImage': {}
                        }
                    },
                    "xAxis": {"type": "time"},
                    "yAxis": {"type": "value", "data": [], "inverse": False},
                    #'tooltip': {
                    #    'order': 'valueDesc',
                    #    'trigger': 'axis'
                    # },
                    "series": [{
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
                        }],
                    }
                )
            .style("height: 350px")
            .classes("border-double")
        )

    def update_coverage_time_plot(self, covdf):
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
        currenttime = time.time() * 1000
        self.coverage_time_chart.options["series"][0]["data"].append([currenttime, coverage])
        self.coverage_time_chart.update()

    def update(self):
        for series in self.echart.options["series"]:
            for i in range(len(series["data"])):
                series["data"][i] = random()
        self.echart.update()

    def add(self):
        self.echart.options["series"].append(
            {"type": "bar", "name": "Gamma", "data": [0.5, 0.6]}
        )
        self.echart.update()

    async def choose_file(self):
        files = await app.native.main_window.create_file_dialog(allow_multiple=True)
        for file in files:
            ui.notify(file)


class BamEventHandler(FileSystemEventHandler):
    def __init__(self, bam_count, logger=None):
        """This class handles events from the file system watcher."""
        self.bam_count = bam_count
        # self.logger = logger
        # self.logger.info(type(self.bam_count))

    def on_created(self, event):
        """This will add a file which is added to the watchfolder to the creates and the info file."""
        # self.logger.info(f"Processing created file {event.src_path}")
        if event.src_path.endswith(".bam"):
            # self.logger.info(f"Created file {event.src_path}")
            self.bam_count["counter"] += 1
            # self.logger.info(f"Bam count is {self.bam_count}")
            if "file" not in self.bam_count:
                self.bam_count["file"] = {}
            self.bam_count["file"][event.src_path] = time.time()
