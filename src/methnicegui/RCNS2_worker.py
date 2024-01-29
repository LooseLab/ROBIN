from nicegui import Tailwind, ui, app

from io import StringIO

import threading
import time
import os, signal
import pysam
import pandas as pd
import subprocess
import tempfile

from methnicegui import models, resources

from methnicegui.merge_bedmethyl import merge_bedmethyl, save_bedmethyl
from methnicegui import submodules

os.environ["CI"] = "1"

HVPATH = os.path.join(
    os.path.dirname(os.path.abspath(submodules.__file__)), "hv_rapidCNS2"
)


class RCNS2_worker:
    def __init__(
        self,
        bamqueue,
        cnv,
        target_coverage,
        mgmt_panel,
        fusion_panel,
        threads=4,
        output_folder=None,
        threshold=0.05,
        showerrors=False,
        browse=False,
    ):
        self.bamqueue = bamqueue
        self.cnv = cnv
        self.threads = threads
        self.outputfolder = output_folder
        self.target_coverage = target_coverage
        self.threshold = threshold
        self.mgmt_panel = mgmt_panel
        self.fusion_panel = fusion_panel
        self.nofiles = False
        self.browse = browse
        self.rcns2_bam_count = 0
        self.showerrors = showerrors
        self.rapidcns_status_txt = {"message": "Waiting for data."}
        self.st_bam_count = 0
        self.resultfolder = os.path.join(self.outputfolder, "results")
        self.subsetbamfolder = os.path.join(self.outputfolder, "subsetbams")
        self.donebamfolder = os.path.join(self.outputfolder, "donebams")
        self.targetsbamfolder = os.path.join(self.outputfolder, "targetsbams")
        self.sortedbamfile = os.path.join(self.subsetbamfolder, "sorted.bam")
        self.merged_bam_file = None
        self.rcns2folder = os.path.join(self.outputfolder, "rcns2")
        self.rcns2_df_store = pd.DataFrame()
        if not self.browse:
            if not os.path.exists(self.resultfolder):
                os.mkdir(self.resultfolder)
            if not os.path.exists(self.subsetbamfolder):
                os.mkdir(self.subsetbamfolder)
            if not os.path.exists(self.donebamfolder):
                os.mkdir(self.donebamfolder)
            if not os.path.exists(self.targetsbamfolder):
                os.mkdir(self.targetsbamfolder)
            if not os.path.exists(self.rcns2folder):
                os.mkdir(self.rcns2folder)
            self.rapidcns2_processing = threading.Thread(
                target=self.rapid_cns2, args=()
            )
            self.rapidcns2_processing.daemon = True
            self.rapidcns2_processing.start()
        else:
            print("Browse mode - not running RapidCNS2 live")

    def rapid_cns2(self):
        """
        This function runs RapidCNS2 on the bam files.
        It runs as a worker to not block the main thread.
        It grabs batches of bam files from the bamforcns list once the list contains 2 or more BAM files.
        :return:
        """
        cov_df = None  # This will be a pandas dataframe to store coverage information.
        not_first_run = False
        batch = 0
        wait_for_batch_size = 2
        while True:
            if self.bamqueue.qsize() > wait_for_batch_size:
                self.rcns2finished = False
                wait_for_batch_size = 0
                batch += 1
                print("Running RapidCNS")
                bams = []
                print(f"there are {self.bamqueue.qsize()} bams for processing")
                self.rapidcns_status_txt[
                    "message"
                ] = f"Processing {self.bamqueue.qsize()} bam files for RapidCNS2."
                while self.bamqueue.qsize() > 0:
                    bams.append(self.bamqueue.get())
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

                # The bam file created here consists of all the reads from the bam files in the most recent batch.

                pysam.cat("-o", tempbam.name, *bams)

                self.rapidcns_status_txt["message"] = "Sorting bam files."

                os.system(
                    f"samtools sort --write-index -@{self.threads} -o {self.sortedbamfile} {tempbam.name} "
                    f">/dev/null 2>&1"
                )

                self.rapidcns_status_txt["message"] = "Running modkit pileup."
                if self.showerrors:
                    print(
                        f"modkit pileup -t {self.threads} --filter-threshold 0.73 --combine-mods {self.sortedbamfile} "
                        f"{self.rapidcnsbamfile}.bed --suppress-progress  >/dev/null 2>&1 "
                    )
                os.system(
                    f"modkit pileup -t {self.threads} --filter-threshold 0.73 --combine-mods {self.sortedbamfile} "
                    f"{self.rapidcnsbamfile}.bed --suppress-progress  >/dev/null 2>&1 "
                )  #

                # ToDo: Implement saving of target coverage data.

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
                self.target_coverage.update_coverage_plot(covdf)
                self.target_coverage.update_coverage_plot_targets(covdf, bedcovdf)
                self.target_coverage.update_coverage_time_plot(covdf)
                self.target_coverage_df = bedcovdf
                self.target_coverage_df["coverage"] = (
                    self.target_coverage_df["bases"] / self.target_coverage_df["length"]
                )
                with self.target_coverage.targ_df:
                    self.target_coverage.targ_df.clear()
                    ui.table.from_pandas(
                        self.target_coverage_df, pagination=10, row_key="name"
                    ).bind_filter_from(self.target_coverage.f, "value").classes(
                        "w-full"
                    )
                # self.log("Merged Bed File Info:")
                # self.log(self.merged_bed_file.info())

                self.rapidcns_status_txt[
                    "message"
                ] = "Running RapidCNS2 methylation classification."

                command = (
                    f"Rscript {HVPATH}/bin/methylation_classification_nanodx_v0.1.R -s "
                    + f"live_{batch} -o {self.rcns2folder} -i {self.rapidcnsbamfile}.bed "
                    + f"-p {HVPATH}/bin/top_probes_hm450.Rdata "
                    + f"--training_data {HVPATH}/bin/capper_top_100k_betas_binarised.Rdata "
                    + f"--array_file {HVPATH}/bin/HM450.hg38.manifest.gencode.v22.Rdata "
                    + f"-t {self.threads} "
                )

                if not self.showerrors:
                    command += ">/dev/null 2>&1"

                returned_value = subprocess.call(
                    command,
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

                columns_greater_than_threshold = (
                    self.rcns2_df_store > self.threshold * 100
                ).any()
                columns_not_greater_than_threshold = ~columns_greater_than_threshold
                result = self.rcns2_df_store.columns[
                    columns_not_greater_than_threshold
                ].tolist()

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
                self.cnv.cnv_plotting(self.donebamfolder, folder=True)

                self.keep_regions(
                    os.path.join(self.donebamfolder, f"{batch}_sorted.bam"), batch
                )

                self.fusion_panel.parse_bams(self.donebamfolder)
                self.mgmtmethylpredict(self.rapidcnsbamfile)

                pass
            time.sleep(5)

            if self.bamqueue.qsize() == 0:
                self.rcns2finished = True
                if (
                    self.nofiles
                    and self.stugeonfinished
                    and self.runfinished
                    and self.rcns2finished
                ):
                    print("All done")
                    os.kill(os.getpid(), signal.SIGINT)

    def mgmtmethylpredict(self, bamfile):
        print("Running MGMT predictor")
        self.rapidcns_status_txt["message"] = "Running MGMT predictor."
        MGMT_BED = f"{HVPATH}/bin/mgmt_hg38.bed"
        os.system(
            f"bedtools intersect -a {bamfile}.bed -b {MGMT_BED} > {self.resultfolder}/mgmt_result.bed"
        )
        print(
            f"bedtools intersect -a {bamfile}.bed -b {MGMT_BED} > {self.resultfolder}/mgmt_result.bed"
        )

        if os.path.getsize(f"{self.resultfolder}/mgmt_result.bed") > 0:
            cmd = f"Rscript {HVPATH}/bin/mgmt_pred_v0.3.R --input={self.resultfolder}/mgmt_result.bed --out_dir={self.resultfolder} --probes={HVPATH}/bin/mgmt_probes.Rdata --model={HVPATH}/bin/mgmt_137sites_mean_model.Rdata --sample=live_analysis"
            os.system(cmd)
            results = pd.read_csv(
                os.path.join(self.resultfolder, "live_analysis_mgmt_status.csv")
            )
            self.mgmt_panel.mgmtable.clear()
            with self.mgmt_panel.mgmtable:
                ui.table.from_pandas(results)
            print(cmd)
            print("MGMT predictor done")
            self.rapidcns_status_txt["message"] = "MGMT predictor done."

        else:
            print("No MGMT sites yet found.")
            self.rapidcns_status_txt["message"] = "No MGMT sites yet found."

    def keep_regions(self, bamtoextract, batch):
        print("Keeping regions")
        bam_out = os.path.join(self.targetsbamfolder, f"{batch}_sorted.bam")
        plot_out = os.path.join(self.targetsbamfolder, f"{batch}_sorted.png")
        bedfile = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)), "unique_genes.bed"
        )
        os.system(
            f"samtools view --write-index -L {bedfile} -@{self.threads} -o {bam_out} {bamtoextract} "
            # f">/dev/null 2>&1"
        )
        merged_bam_out = os.path.join(self.targetsbamfolder, f"{batch}_merged.bam")
        to_be_merged = os.path.join(self.targetsbamfolder, f"*_sorted.bam")
        os.system(
            f"samtools merge --write-index -@{self.threads} -f {merged_bam_out} {to_be_merged}"
            # f">/dev/null 2>&1"
        )
        os.system(f"rm {to_be_merged}")
        os.system(f"mv {merged_bam_out} {bam_out}")
        os.system(f"mv {merged_bam_out}.csi {bam_out}.csi")
        print(
            f"methylartist locus -i chr10:129466536-129467536 -b {bam_out} -o {plot_out}  --motif CG --mods m"
        )
        os.system(
            f"methylartist locus -i chr10:129466536-129467536 -b {bam_out} -o {plot_out}  --motif CG --mods m"
        )
        if os.path.exists(plot_out):
            self.mgmt_panel.mgmtplot.clear()
            with self.mgmt_panel.mgmtplot.classes("w-full"):
                ui.image(plot_out).props("fit=scale-down")

    def status_rcns2(self):
        ui.label().bind_text_from(
            self.rapidcns_status_txt,
            "message",
            backward=lambda n: f"RapidCNS2 Status: {n}",
        )

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


    def load_prior_data(self):
        """
        Load prior data from a file.
        :return:
        """
        print("Loading prior data")
        self.rcns2_df_store = pd.read_csv(
            os.path.join(self.resultfolder, "rcns2_scores.csv")
        ).set_index("timestamp")
        print(self.rcns2_df_store)
        columns_greater_than_threshold = (
            self.rcns2_df_store > self.threshold * 100
        ).any()
        columns_not_greater_than_threshold = ~columns_greater_than_threshold
        result = self.rcns2_df_store.columns[
            columns_not_greater_than_threshold
        ].tolist()

        self.update_rcns2_time_chart(self.rcns2_df_store.drop(columns=result))

        batch = find_largest_integer(self.rcns2folder, "live_", "_votes.tsv")

        if batch:
            scores = pd.read_table(
                f"{self.rcns2folder}/live_{batch}_votes.tsv", delim_whitespace=True
            )
            scores = scores.sort_values(by=["cal_Freq"], ascending=False).head(10)
            self.update_rcns2_plot(
                scores.index.to_list(),
                list(scores["cal_Freq"].values / 100),
                "ALL",
            )

        self.rapidcns_status_txt["message"] = "Predictions Complete."
        ##Here we have to pass a folder of files rather than an individual file.
        self.cnv.update_cnv_dict = {}
        self.cnv.cnv_plotting(self.donebamfolder, folder=True)
        plot_out = os.path.join(self.targetsbamfolder, f"{batch}_sorted.png")
        if os.path.exists(plot_out):
            self.mgmt_panel.mgmtplot.clear()
            with self.mgmt_panel.mgmtplot.classes("w-full"):
                ui.image(plot_out).props("fit=scale-down")
            results = pd.read_csv(
                os.path.join(self.resultfolder, "live_analysis_mgmt_status.csv")
            )
            self.mgmt_panel.mgmtable.clear()
            with self.mgmt_panel.mgmtable:
                ui.table.from_pandas(results)
        print(self.donebamfolder)
        self.fusion_panel.parse_bams(self.donebamfolder)

    def replay_prior_data(self):
        self.background_task = threading.Thread(target=self._replay_prior_data, args=())
        self.background_task.daemon = True
        self.background_task.start()

    def _replay_prior_data(self):
        """
        Replay prior data from a file.
        :return:
        """
        print("Replaying prior data")
        self.rcns2_df_store = pd.read_csv(
            os.path.join(self.resultfolder, "rcns2_scores.csv")
        )
        self.rcns2_df_store["offset"] = (
            self.rcns2_df_store["timestamp"] - self.rcns2_df_store["timestamp"][0]
        )
        self.rcns2_df_store.set_index("timestamp")
        # bamsubsets = self.cnv.replay(self.donebamfolder)
        self.rapidcns_status_txt["message"] = f"Replaying prior data."

        scale_factor = 1200
        counter = 0
        current_time = time.time()
        #self.cnv.update_cnv_dict = {}
        for index, row in self.rcns2_df_store.iterrows():
            counter += 1
            # print (index, row)

            if time.time() < current_time + (row["offset"] / 1000 / scale_factor):
                time.sleep(
                    current_time + (row["offset"] / 1000 / scale_factor) - time.time()
                )
            else:
                time.sleep(0.1)

            #self.cnv.cnv_plotting(
            #    os.path.join(self.donebamfolder, f"{counter}_sorted.bam")
            #)
            temp_rcns2_df_store = (
                self.rcns2_df_store.head(counter)
                .drop(columns=["offset"])
                .set_index("timestamp")
            )
            columns_greater_than_threshold = (
                temp_rcns2_df_store > self.threshold * 100
            ).any()
            columns_not_greater_than_threshold = ~columns_greater_than_threshold
            result = temp_rcns2_df_store.columns[
                columns_not_greater_than_threshold
            ].tolist()

            self.update_rcns2_time_chart(temp_rcns2_df_store.drop(columns=result))
            self.rapidcns_status_txt[
                "message"
            ] = f"Replaying prior RCNS2 data - step {counter}."
            lastrow = row.drop(["offset", "timestamp"])
            lastrow_plot = lastrow.sort_values(ascending=False).head(10)
            self.update_rcns2_plot(
                lastrow_plot.index.to_list(),
                list(lastrow_plot.values / 100),
                "replay",
            )
            plot_out = os.path.join(self.targetsbamfolder, f"{counter}_sorted.png")
            if os.path.exists(plot_out):
                self.mgmt_panel.mgmtplot.clear()
                with self.mgmt_panel.mgmtplot.classes("w-full"):
                    ui.image(plot_out).props("fit=scale-down")
                results = pd.read_csv(
                    os.path.join(self.resultfolder, "live_analysis_mgmt_status.csv")
                )
                self.mgmt_panel.mgmtable.clear()
                with self.mgmt_panel.mgmtable:
                    ui.table.from_pandas(results)

        self.rapidcns_status_txt["message"] = f"Viewing historical RCNS2 data."
        #replaycontrol.visible = True


def find_largest_integer(directory_path, prefix, suffix):
    files = os.listdir(directory_path)

    # Filter files that match the pattern live_N_sample.txt
    matching_files = [
        file for file in files if file.startswith(prefix) and file.endswith(suffix)
    ]

    if not matching_files:
        print("No matching files found.")
        return None

    # Extract integers from the file names
    integers = [int(file.split("_")[1]) for file in matching_files]

    # Find the largest integer
    largest_integer = max(integers)

    return largest_integer
