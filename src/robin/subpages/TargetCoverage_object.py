"""
This code has the following known issues:
1. The code lacks documentation.
2. The code will load earlier iterations of the data on first load. This is not ideal!

"""

from robin.subpages.base_analysis import BaseAnalysis
import natsort
from robin import theme, resources
import pandas as pd
import numpy as np
import os
import sys
import click
import time
from pathlib import Path
from nicegui import ui, run, app  # , background_tasks
from io import StringIO
import pysam
import tempfile
import shutil
import queue


os.environ["CI"] = "1"


def process_annotations(record: dict) -> dict:
    """
    This function takes a dictionary record from a vcf file and explodes the record into multiple records based on the contents of the INFO field.
    We expect some unit of 16 entries in the INFO field. Where there are multiples of 16 entries, we split them into a new record entry for that specific mutation.
    """
    # print(record["INFO"])
    if "INFO" not in record.keys():
        return {}, {}
    annotations = record["INFO"]
    # This dictionary holds the information for a single record
    rec_dict = {}
    # This dictionary holds one or more records derived from the annotation field.
    ann_dict = {}
    for ann in annotations.split(";"):
        if "=" in ann:
            mykey = ann.split("=")[0]
            myvalue = ann.split("=")[1]
            if "|" in myvalue:
                if mykey == "ANN":
                    if len(myvalue.split("|")) == 16:
                        try:
                            myvalues = myvalue.split("|")
                            count = 0
                            ann_dict[count] = dict()
                            ann_dict[count]["Allele"] = myvalues[0]
                            ann_dict[count]["Annotation"] = myvalues[1]
                            ann_dict[count]["Annotation_Impact"] = myvalues[2]
                            ann_dict[count]["Gene_Name"] = myvalues[3]
                            ann_dict[count]["Gene_ID"] = myvalues[4]
                            ann_dict[count]["Feature_Type"] = myvalues[5]
                            ann_dict[count]["Feature_ID"] = myvalues[6]
                            ann_dict[count]["Transcript_BioType"] = myvalues[7]
                            ann_dict[count]["Rank"] = myvalues[8]
                            ann_dict[count]["HGVS.c"] = myvalues[9]
                            ann_dict[count]["HGVS.p"] = myvalues[10]
                            ann_dict[count]["cDNA.pos / cDNA.length"] = myvalues[11]
                            ann_dict[count]["CDS.pos / CDS.length"] = myvalues[12]
                            ann_dict[count]["AA.pos / AA.length"] = myvalues[13]
                            ann_dict[count]["Distance"] = myvalues[14]
                            ann_dict[count]["ERRORS / WARNINGS / INFO"] = myvalues[15]
                        except Exception as e:
                            print(f"Error: {e}")
                            sys.exit()
                    elif len(myvalue.split("|")) > 16:
                        count = 0
                        for chunk in myvalue.split(","):
                            try:
                                myvalues = chunk.split("|")
                                ann_dict[count] = dict()
                                ann_dict[count]["Allele"] = myvalues[0]
                                ann_dict[count]["Annotation"] = myvalues[1]
                                ann_dict[count]["Annotation_Impact"] = myvalues[2]
                                ann_dict[count]["Gene_Name"] = myvalues[3]
                                ann_dict[count]["Gene_ID"] = myvalues[4]
                                ann_dict[count]["Feature_Type"] = myvalues[5]
                                ann_dict[count]["Feature_ID"] = myvalues[6]
                                ann_dict[count]["Transcript_BioType"] = myvalues[7]
                                ann_dict[count]["Rank"] = myvalues[8]
                                ann_dict[count]["HGVS.c"] = myvalues[9]
                                ann_dict[count]["HGVS.p"] = myvalues[10]
                                ann_dict[count]["cDNA.pos / cDNA.length"] = myvalues[11]
                                ann_dict[count]["CDS.pos / CDS.length"] = myvalues[12]
                                ann_dict[count]["AA.pos / AA.length"] = myvalues[13]
                                ann_dict[count]["Distance"] = myvalues[14]
                                ann_dict[count]["ERRORS / WARNINGS / INFO"] = myvalues[
                                    15
                                ]
                                count += 1
                            except Exception as e:
                                print(f"Error: {e}")
                                sys.exit()
            else:
                rec_dict[mykey] = myvalue
    return ann_dict, rec_dict


def parse_vcf(vcf_file):
    header = "CHROM POS ID REF ALT QUAL FILTER INFO FORMAT GT".split()
    vcf = pd.read_csv(vcf_file, delimiter="\t", comment="#", names=header)
    result = None
    if len(vcf) > 0:
        explodedvcf = []
        for record in vcf.to_dict("records"):
            result, result2 = process_annotations(record)
            if len(result) > 1:
                for res in result:
                    dat = {**record, **result[res], **result2}
                    explodedvcf.append(dat)

        vcf = pd.DataFrame.from_records(explodedvcf)
        if "INFO" in vcf.columns:
            vcf = vcf.drop(columns=["INFO"]).drop_duplicates()
        else:
            vcf = vcf.drop_duplicates()

        def set_unique_values(series):
            set_series = set(series)
            if len(set_series) > 0:
                return "{}".format(", ".join(map(str, set_series)))
            return None

        if len(vcf) > 0:
            shared_columns = [
                "CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "FORMAT",
                "GT",
                "Allele",
            ]

            # Define columns to be aggregated
            non_shared_columns = [
                col for col in vcf.columns if col not in shared_columns
            ]

            vcf = vcf.replace({np.nan: None})
            result = (
                vcf.groupby(shared_columns)[non_shared_columns]
                .agg(set_unique_values)
                .reset_index()
            )

            vcf = result

            try:
                vcf.to_csv(f"{vcf_file}.csv", index=False)
                # print(f"VCF file saved as {vcf_file}.csv")
            except Exception as e:
                print(e)
                sys.exit(1)


def run_clair3(bamfile, bedfile, workdir, workdirout, threads, reference):
    # ToDo: handle any platform
    if sys.platform in ["darwin", "linux"]:
        runcommand = (
            f"docker run -it -v {workdir}:{workdir} "
            f"-v {workdirout}:{workdirout} "
            f"-v {reference}:{reference} "
            f"-v {reference}.fai:{reference}.fai "
            f"hkubal/clairs-to:latest "
            f"/opt/bin/run_clairs_to "
            f"--tumor_bam_fn {bamfile} "
            f"--ref_fn {reference} "
            f"--threads {threads} "
            f"--platform ont_r10_guppy_hac_5khz "
            f"--output_dir {workdirout} -b {bedfile}"
            # f" >/dev/null 2>&1"
        )
        os.system(runcommand)
        shutil.copy2(f"{workdirout}/snv.vcf.gz", f"{workdirout}/output_done.vcf.gz")
        shutil.copy2(
            f"{workdirout}/indel.vcf.gz", f"{workdirout}/output_indel_done.vcf.gz"
        )

        command = f"snpEff -q hg38 {workdirout}/output_done.vcf.gz > {workdirout}/snpeff_output.vcf"
        os.system(command)
        command = f"SnpSift annotate {os.path.join(os.path.dirname(os.path.abspath(resources.__file__)),'clinvar.vcf')} {workdirout}/snpeff_output.vcf > {workdirout}/snpsift_output.vcf"
        os.system(command)
        command = f"snpEff -q hg38 {workdirout}/output_indel_done.vcf.gz > {workdirout}/snpeff_indel_output.vcf"
        os.system(command)
        command = f"SnpSift annotate {os.path.join(os.path.dirname(os.path.abspath(resources.__file__)), 'clinvar.vcf')} {workdirout}/snpeff_indel_output.vcf > {workdirout}/snpsift_indel_output.vcf"
        os.system(command)
        parse_vcf(f"{workdirout}/snpsift_output.vcf")
        parse_vcf(f"{workdirout}/snpsift_indel_output.vcf")


def get_covdfs(bamfile):
    """
    This function runs modkit on a bam file and extracts the methylation data.
    """
    try:
        # no_secondary_bam = tempfile.NamedTemporaryFile(dir=output, suffix=".bam").name
        # command = f"samtools view -@2 -h -F 0x100 --write-index {bamfile} -o {no_secondary_bam}"
        # os.system(command)
        newcovdf = pd.read_csv(StringIO(pysam.coverage(f"{bamfile}")), sep="\t")
        newcovdf.drop(
            columns=["coverage", "meanbaseq", "meanmapq"],
            inplace=True,
        )

        bedcovdf = pd.read_csv(
            StringIO(
                pysam.bedcov(
                    os.path.join(
                        os.path.dirname(os.path.abspath(resources.__file__)),
                        "unique_genes.bed",
                    ),
                    f"{bamfile}",
                )
            ),
            names=["chrom", "startpos", "endpos", "name", "bases"],
            sep="\t",
        )
    except Exception as e:
        print(e)
        return None
    return newcovdf, bedcovdf


def subset_bam(bamfile, targets, output):
    pysam.view("-L", f"{targets}", "-o", f"{output}", f"{bamfile}")


def sort_bam(bamfile, output, threads):
    pysam.sort(f"-@{threads}", "-o", output, bamfile)
    pysam.index(f"{output}", f"{output}.bai")
    print(f"Sorted bam file saved as {output}")


def run_bedmerge(newcovdf, cov_df_main, bedcovdf, bedcov_df_main):
    merged_df = pd.merge(
        newcovdf,
        cov_df_main,
        on=["#rname", "startpos", "endpos"],
        suffixes=("_df1", "_df2"),
    )
    merged_df["numreads"] = merged_df["numreads_df1"] + merged_df["numreads_df2"]
    merged_df["covbases"] = merged_df["covbases_df1"] + merged_df["covbases_df2"]
    merged_df["meandepth"] = merged_df["meandepth_df1"] + merged_df["meandepth_df2"]

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

    merged_bed_df = pd.merge(
        bedcovdf,
        bedcov_df_main,
        on=["chrom", "startpos", "endpos", "name"],
        suffixes=("_df1", "_df2"),
    )
    merged_bed_df["bases"] = merged_bed_df["bases_df1"] + merged_bed_df["bases_df2"]
    merged_bed_df.drop(columns=["bases_df1", "bases_df2"], inplace=True)
    return merged_df, merged_bed_df


def run_bedtools(bamfile, bedfile, tempbamfile):
    """
    This function extracts the target sites from the bamfile.
    """
    try:
        os.system(f"bedtools intersect -a {bamfile} -b {bedfile} > {tempbamfile}")
        pysam.index(tempbamfile)
    except Exception as e:
        print(e)


class TargetCoverage(BaseAnalysis):
    def __init__(self, *args, target_panel=None, reference=None, **kwargs):
        self.callthreshold = 10
        self.clair3running = False
        self.targets_exceeding_threshold = {}
        self.targetbamfile = {}
        self.covtable = None
        self.covtable_row_count = 0
        self.coverage_over_time = {} #np.empty((0, 2))
        self.cov_df_main = {}
        self.bedcov_df_main = {}
        self.SNPqueue = queue.Queue()
        self.reference = reference
        if self.reference:
            self.snp_calling = True
        else:
            self.snp_calling = False
        if self.snp_calling:
            self.SNP_timer_run()
        self.target_panel = target_panel

        if self.target_panel == "rCNS2":
            self.bedfile = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)),
                "rCNS2_panel_name_uniq.bed",
            )
        elif self.target_panel == "AML":
            self.bedfile = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)),
                "AML_panel_name_uniq.bed",
            )
        super().__init__(*args, **kwargs)

    def SNP_timer_run(self):
        self.snp_timer = ui.timer(0.1, self._snp_worker)

    async def _snp_worker(self):
        """
        This function takes reads from the queue and adds them to the background thread for processing.
        """
        self.snp_timer.active = False
        if not self.SNPqueue.empty():
            while not self.SNPqueue.empty():
                gene_list, bamfile, bedfile = self.SNPqueue.get()

            workdirout = os.path.join(
                self.check_and_create_folder(self.output, self.sampleID), "clair3"
            )
            if not os.path.exists(workdirout):
                os.mkdir(workdirout)
            # bamfile, bedfile, workdir, workdirout, threads
            self.clair3running = True
            shutil.copy2(bedfile, f"{bedfile}2")
            await run.cpu_bound(
                run_clair3,
                f"{bamfile}",
                f"{bedfile}2",
                self.check_and_create_folder(self.output, self.sampleID),
                workdirout,
                self.threads,
                self.reference,
            )
            self.clair3running = False
        # else:
        #    await asyncio.sleep(1)
        self.snp_timer.active = True

    def setup_ui(self):
        if self.summary:
            with self.summary:
                ui.label("Current coverage estimates: Unknown")
        with ui.card().classes("w-full"):
            ui.label("Coverage Data").style(
                "color: #6E93D6; font-size: 150%; font-weight: 300"
            ).tailwind("drop-shadow", "font-bold")
            with ui.grid(columns=2).classes("w-full h-auto"):
                with ui.card().classes(
                    f"min-[{self.MENU_BREAKPOINT+1}px]:col-span-1 max-[{self.MENU_BREAKPOINT}px]:col-span-2"
                ):
                    self.create_coverage_plot("Chromosome Coverage")

                with ui.card().classes(
                    f"min-[{self.MENU_BREAKPOINT+1}px]:col-span-1 max-[{self.MENU_BREAKPOINT}px]:col-span-2"
                ):
                    self.create_coverage_plot_targets("Target Coverage")
        with ui.card().classes("w-full"):
            ui.label("Target Outliers").style(
                "color: #6E93D6; font-size: 150%; font-weight: 300"
            ).tailwind("drop-shadow", "font-bold")
            with ui.column().classes("w-full"):
                with ui.card().classes("w-full"):
                    self.create_target_boxplot()
        with ui.card().classes("w-full"):
            ui.label("Coverage over time").style(
                "color: #6E93D6; font-size: 150%; font-weight: 300"
            ).tailwind("drop-shadow", "font-bold")
            with ui.column().classes("w-full"):
                with ui.card().classes("w-full"):
                    self.create_coverage_time_chart()
        with ui.card().classes("w-full"):
            ui.label("Coverage over targets").style(
                "color: #6E93D6; font-size: 150%; font-weight: 300"
            ).tailwind("drop-shadow", "font-bold")
            self.targ_df = ui.row().classes("w-full").style("height: 900px")
        with ui.card().classes("w-full"):
            ui.label("Candidate SNPs").style(
                "color: #6E93D6; font-size: 150%; font-weight: 300"
            ).tailwind("drop-shadow", "font-bold")
            with ui.row().classes("w-full"):
                self.SNPplaceholder = (
                    ui.card().tight().classes("w-full overflow-x-auto")
                )
                with self.SNPplaceholder:
                    ui.label(
                        "Candidate SNPs will be displayed here. SNPs are called based on available data at that time."
                    )
                # self.SNPview = SNPview(self.SNPplaceholder)
                # ui.timer(0.1,lambda: self.SNPview.renderme(), once=True)
            ui.label("Candidate IN/DELs").style(
                "color: #6E93D6; font-size: 150%; font-weight: 300"
            ).tailwind("drop-shadow", "font-bold")
            with ui.row().classes("w-full"):
                self.INDELplaceholder = (
                    ui.card().tight().classes("w-full overflow-x-auto")
                )

                with self.INDELplaceholder:
                    ui.label(
                        "Candidate IN/DELs will be displayed here. IN/DELs are called based on available data at that time."
                    )
                # self.INDELview = SNPview(self.INDELplaceholder)
                # ui.timer(0.1,lambda: self.INDELview.renderme(), once=True)
        if self.browse:
            ui.timer(0.1, callback=self.show_previous_data, once=True)
        else:
            ui.timer(30, lambda: self.show_previous_data())

    def create_coverage_plot(self, title):
        self.echart3 = (
            ui.echart(
                {
                    "textStyle": {"fontFamily": "Fira Sans, Fira Mono"},
                    "grid": {"containLabel": True},
                    "title": {"text": title},
                    "toolbox": {"show": True, "feature": {"saveAsImage": {}}},
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

    def create_coverage_plot_targets(self, title):
        self.echart4 = (
            ui.echart(
                {
                    "textStyle": {"fontFamily": "Fira Sans, Fira Mono"},
                    "grid": {"containLabel": True},
                    "title": {"text": title},
                    "toolbox": {"show": True, "feature": {"saveAsImage": {}}},
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

    def create_coverage_time_chart(self):
        self.coverage_time_chart = (
            ui.echart(
                {
                    "textStyle": {"fontFamily": "Fira Sans, Fira Mono"},
                    "grid": {"containLabel": True},
                    "title": {"text": "Coverage Over Time"},
                    "toolbox": {"show": True, "feature": {"saveAsImage": {}}},
                    "xAxis": {"type": "time"},
                    "yAxis": {"type": "value", "data": [], "inverse": False},
                    #'tooltip': {
                    #    'order': 'valueDesc',
                    #    'trigger': 'axis'
                    # },
                    "series": [
                        {
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
                        }
                    ],
                }
            )
            .style("height: 350px")
            .classes("border-double")
        )

    def create_target_boxplot(self):
        self.target_boxplot = (
            ui.echart(
                {
                    "textStyle": {"fontFamily": "Fira Sans, Fira Mono"},
                    "title": [
                        {
                            "text": "Target Coverage",
                        },
                    ],
                    "dataset": [
                        {
                            "id": "raw",
                            "dimensions": [
                                "chrom",
                                "min",
                                "Q1",
                                "median",
                                "Q3",
                                "max",
                                "chrom_index",
                            ],
                            "source": [],
                        },
                        {
                            "id": "rawdata",
                            "dimensions": ["chrom", "coverage", "name"],
                            "source": [],
                        },
                        {
                            "id": "outliers",
                            "dimensions": ["chrom", "coverage", "name"],
                            "source": [],
                        },
                        {
                            "id": "globaloutliers",
                            "dimensions": ["chrom", "coverage", "name"],
                            "source": [],
                        },
                    ],
                    "tooltip": {"trigger": "item", "axisPointer": {"type": "shadow"}},
                    "dataZoom": [
                        {
                            "type": "slider",
                            "yAxisIndex": 0,
                        }
                    ],
                    "grid": {"left": "10%", "right": "10%", "bottom": "15%"},
                    "xAxis": {
                        "type": "category",
                        "name": "Chromosome",
                        "boundaryGap": True,
                        "nameGap": 30,
                        "splitArea": {"show": False},
                        "splitLine": {"show": False},
                    },
                    "yAxis": {
                        # "type": 'value',
                        "name": "Coverage",
                        # "splitArea": {
                        #    "show": True
                        # }
                    },
                    "legend": {
                        "data": ["boxplot", "raw data", "outliers", "global outliers"],
                        "selected": {
                            "boxplot": True,
                            "global outliers": True,
                            "raw data": False,
                            "outliers": False,
                        },
                    },
                    "series": [
                        {
                            "name": "boxplot",
                            "type": "boxplot",
                            "datasetId": "raw",
                            "encode": {
                                "y": ["min", "Q1", "median", "Q3", "max"],
                                "x": "chrom",
                                "itemName": ["chrom"],
                                "tooltip": ["min", "Q1", "median", "Q3", "max"],
                            },
                        },
                        {
                            "name": "raw data",
                            "type": "scatter",
                            "datasetId": "rawdata",
                            "encode": {
                                "y": "coverage",
                                "x": "chrom",
                                "label": "name",
                                "itemName": "name",
                                "tooltip": ["chrom", "name", "coverage"],
                            },
                            "label": {
                                "show": True,
                                "position": "right",
                                "itemName": "name",
                                "color": "black",
                                "fontSize": 16,
                            },
                        },
                        {
                            "name": "outliers",
                            "type": "scatter",
                            "datasetId": "outliers",
                            "encode": {
                                "y": "coverage",
                                "x": "chrom",
                                "label": "name",
                                "itemName": "name",
                                "tooltip": ["chrom", "name", "coverage"],
                            },
                            "label": {
                                "show": True,
                                "position": "right",
                                "itemName": "name",
                                "color": "black",
                                "fontSize": 16,
                            },
                        },
                        {
                            "name": "global outliers",
                            "type": "scatter",
                            "datasetId": "globaloutliers",
                            "encode": {
                                "y": "coverage",
                                "x": "chrom",
                                "label": "name",
                                "itemName": "name",
                                "tooltip": ["chrom", "name", "coverage"],
                            },
                            "label": {
                                "show": True,
                                "position": "right",
                                "itemName": "name",
                                "color": "black",
                                "fontSize": 16,
                            },
                        },
                    ],
                }
            )
            .style("height: 500px")
            .classes("border-double")
        )

    def update_target_boxplot(self, dataframe):
        """
        Replaces the data in the RapidCNS2 plot.
        :param dataframe: a pandas dataframe of
        :return:
        """
        if "coverage" in dataframe.columns:
            # Naturally sort the unique chromosome values
            sorted_chroms = natsort.natsorted(dataframe["chrom"].unique())

            # Create a lookup index for chromosomes
            chrom_lookup = {chrom: idx for idx, chrom in enumerate(sorted_chroms)}

            # Add the index column to the original DataFrame
            dataframe["index"] = dataframe["chrom"].map(chrom_lookup)

            # Aggregate the data by 'chrom'
            aggregated_data = (
                dataframe.groupby("chrom")
                .agg(
                    min_coverage=("coverage", "min"),
                    Q1_coverage=("coverage", lambda x: np.percentile(x, 25)),
                    median_coverage=("coverage", "median"),
                    Q3_coverage=("coverage", lambda x: np.percentile(x, 75)),
                    max_coverage=("coverage", "max"),
                    index=(
                        "index",
                        "first",
                    ),  # Get the corresponding index for the chromosome
                )
                .reset_index()
            )

            # Sort the aggregated data by 'chrom' naturally using natsort
            aggregated_data["chrom"] = pd.Categorical(
                aggregated_data["chrom"],
                categories=natsort.natsorted(aggregated_data["chrom"].unique()),
                ordered=True,
            )
            aggregated_data = aggregated_data.sort_values("chrom").reset_index(
                drop=True
            )

            # Format the result as required
            result = [
                ["chrom", "min", "Q1", "median", "Q3", "max", "index"]
            ] + aggregated_data.values.tolist()

            df = dataframe

            # Function to identify outliers within each chromosome group
            def identify_outliers_per_chromosome(df):
                outliers = pd.DataFrame()
                for chrom in df["chrom"].unique():
                    chrom_data = df[df["chrom"] == chrom]
                    Q1 = chrom_data["coverage"].quantile(0.25)
                    Q3 = chrom_data["coverage"].quantile(0.75)
                    IQR = Q3 - Q1
                    lower_bound = Q1 - 1.5 * IQR
                    upper_bound = Q3 + 1.5 * IQR
                    chrom_outliers = chrom_data[
                        (chrom_data["coverage"] < lower_bound)
                        | (chrom_data["coverage"] > upper_bound)
                    ]
                    outliers = pd.concat([outliers, chrom_outliers])
                return outliers

            def identify_outliers(df):
                Q1 = df["coverage"].quantile(0.25)
                Q3 = df["coverage"].quantile(0.75)
                IQR = Q3 - Q1
                lower_bound = Q1 - 1.5 * IQR
                upper_bound = Q3 + 1.5 * IQR
                return df[
                    (df["coverage"] < lower_bound) | (df["coverage"] > upper_bound)
                ]

            outliers = identify_outliers_per_chromosome(df)
            globaloutliers = identify_outliers(df)

            self.target_boxplot.options["dataset"][0]["source"] = result
            self.target_boxplot.options["dataset"][1]["source"] = dataframe[
                ["chrom", "coverage", "name"]
            ].values.tolist()
            self.target_boxplot.options["dataset"][2]["source"] = outliers[
                ["chrom", "coverage", "name"]
            ].values.tolist()
            self.target_boxplot.options["dataset"][3]["source"] = globaloutliers[
                ["chrom", "coverage", "name"]
            ].values.tolist()
            self.target_boxplot.update()

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
        self.echart3.options["title"]["text"] = "Per Chromosome Coverage"
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
        self.echart4.options["title"]["text"] = "Per Chromosome Target Coverage"
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

    def update_coverage_time_plot(self):
        if self.browse:
            filepath = os.path.join(
                self.check_and_create_folder(self.output, self.sampleID),
                "coverage_time_chart.npy",
            )
        else:
            filepath = os.path.join(self.output, "coverage_time_chart.npy")
        if os.path.isfile(filepath):
            self.coverage_over_time = np.load(filepath)
            self.coverage_time_chart.options["series"][0][
                "data"
            ] = self.coverage_over_time.tolist()
            self.coverage_time_chart.update()

    def update_target_coverage_table(self):
        if not self.covtable:
            self.target_coverage_df["coverage"] = self.target_coverage_df[
                "coverage"
            ].round(2)
            with self.targ_df:
                self.targ_df.clear()
                self.covtable = (
                    ui.table.from_pandas(self.target_coverage_df, pagination=25)
                    .props("dense")
                    .classes("w-full")
                    .style("height: 900px")
                    .style("font-size: 100%; font-weight: 300")
                )
                for col in self.covtable.columns:
                    col["sortable"] = True
                    if col["name"] == "coverage":
                        col[":format"] = "value => value.toFixed(2)"
                self.covtable.add_slot(
                    "body-cell-coverage",
                    """
                    <q-td key="coverage" :props="props">
                        <q-badge :color="props.value < 10 ? 'red' : 'green'">
                            {{ props.value }}
                        </q-badge>
                    </q-td>
                """,
                )
                with self.covtable.add_slot("top-right"):
                    with ui.input(placeholder="Search").props("type=search").bind_value(
                        self.covtable, "filter"
                    ).add_slot("append"):
                        ui.icon("search")
        else:
            self.covtable.update_rows(self.target_coverage_df.to_dict(orient="records"))

    async def process_bam(self, bamfile, timestamp):
        # loop = asyncio.get_event_loop()
        # newcovdf, bedcovdf = await loop.run_in_executor(None, get_covdfs, bamfile)

        newcovdf, bedcovdf = await run.cpu_bound(get_covdfs, bamfile)

        tempbamfile = tempfile.NamedTemporaryFile(
            dir=self.check_and_create_folder(self.output, self.sampleID), suffix=".bam"
        )
        # await loop.run_in_executor(
        #    None, run_bedtools, bamfile, self.bedfile, tempbamfile.name
        # )
        # run_bedtools(bamfile, self.bedfile, tempbamfile.name)
        await run.cpu_bound(run_bedtools, bamfile, self.bedfile, tempbamfile.name)
        # )

        if pysam.AlignmentFile(tempbamfile.name, "rb").count(until_eof=True) > 0:
            if self.sampleID not in self.targetbamfile.keys():
                self.targetbamfile[self.sampleID] = os.path.join(
                    self.check_and_create_folder(self.output, self.sampleID),
                    "target.bam",
                )
                shutil.copy2(tempbamfile.name, self.targetbamfile[self.sampleID])
                os.remove(f"{tempbamfile.name}.bai")
            else:
                tempbamholder = tempfile.NamedTemporaryFile(
                    dir=self.check_and_create_folder(self.output, self.sampleID),
                    suffix=".bam",
                )
                pysam.cat(
                    "-o", tempbamholder.name, self.targetbamfile[self.sampleID], tempbamfile.name
                )
                shutil.copy2(tempbamholder.name, self.targetbamfile[self.sampleID])
                try:
                    os.remove(f"{tempbamholder.name}.bai")
                except FileNotFoundError:
                    pass
        else:
            os.remove(f"{tempbamfile.name}.bai")

        try:
            os.remove(
                f"{tempbamfile.name}.bai"
            )  # Ensure removal of .bai file if exists
        except FileNotFoundError:
            pass

        if self.sampleID not in self.cov_df_main.keys():
            self.cov_df_main[self.sampleID] = newcovdf
            self.bedcov_df_main[self.sampleID] = bedcovdf
        else:
            self.cov_df_main[self.sampleID], self.bedcov_df_main[self.sampleID] = await run.cpu_bound(
                run_bedmerge,
                newcovdf,
                self.cov_df_main[self.sampleID],
                bedcovdf,
                self.bedcov_df_main[self.sampleID],
            )


        bases = self.cov_df_main[self.sampleID]["covbases"].sum()
        genome = self.cov_df_main[self.sampleID]["endpos"].sum()
        coverage = bases / genome
        # if timestamp:
        #    currenttime = timestamp * 1000
        # else:
        currenttime = time.time() * 1000
        if self.sampleID not in self.coverage_over_time.keys():
            self.coverage_over_time[self.sampleID] = np.empty((0, 2))
        self.coverage_over_time[self.sampleID] = np.vstack(
            [self.coverage_over_time[self.sampleID], [(currenttime, coverage)]]
        )

        np.save(
            os.path.join(
                self.check_and_create_folder(self.output, self.sampleID),
                "coverage_time_chart.npy",
            ),
            self.coverage_over_time[self.sampleID],
        )

        self.cov_df_main[self.sampleID].to_csv(
            os.path.join(
                self.check_and_create_folder(self.output, self.sampleID),
                "coverage_main.csv",
            ),
            index=False,
        )
        # await asyncio.sleep(0.01)

        # self.update_coverage_plot_targets(self.cov_df_main, self.bedcov_df_main)
        self.bedcov_df_main[self.sampleID].to_csv(
            os.path.join(
                self.check_and_create_folder(self.output, self.sampleID),
                "bed_coverage_main.csv",
            ),
            index=False,
        )
        # await asyncio.sleep(0.01)
        # self.update_coverage_time_plot(self.cov_df_main, timestamp)
        # await asyncio.sleep(0.01)
        self.target_coverage_df = self.bedcov_df_main[self.sampleID]
        self.target_coverage_df["length"] = (
            self.target_coverage_df["endpos"] - self.target_coverage_df["startpos"] + 1
        )

        self.target_coverage_df["coverage"] = (
            self.target_coverage_df["bases"] / self.target_coverage_df["length"]
        )
        self.target_coverage_df.to_csv(
            os.path.join(
                self.check_and_create_folder(self.output, self.sampleID),
                "target_coverage.csv",
            ),
            index=False,
        )
        #if self.summary:
        #    with self.summary:
        #        self.summary.clear()
        #        with ui.row():
        #            ui.label("Coverage Depths - ")
        #            ui.label(
        #                f"Global Estimated Coverage: {(self.cov_df_main['covbases'].sum()/self.cov_df_main['endpos'].sum()):.2f}x"
        #            )
        #            ui.label(
        #                f"Targets Estimated Coverage: {(self.bedcov_df_main['bases'].sum()/self.bedcov_df_main['length'].sum()):.2f}x"
        #            )
        run_list = self.target_coverage_df[
            self.target_coverage_df["coverage"].ge(self.callthreshold)
        ]
        if self.sampleID not in self.targets_exceeding_threshold.keys():
            self.targets_exceeding_threshold[self.sampleID] = 0
        if self.reference:
            if (
                len(run_list) > self.targets_exceeding_threshold[self.sampleID]
                and not self.clair3running
            ):
                self.targets_exceeding_threshold[self.sampleID] = len(run_list)
                run_list[["chrom", "startpos", "endpos"]].to_csv(
                    os.path.join(
                        self.check_and_create_folder(self.output, self.sampleID),
                        "targets_exceeding_threshold.bed",
                    ),
                    sep="\t",
                    header=None,
                    index=None,
                )
                clair3workdir = os.path.join(
                    self.check_and_create_folder(self.output, self.sampleID), "clair3"
                )
                if not os.path.exists(clair3workdir):
                    os.mkdir(clair3workdir)

                await run.cpu_bound(
                    sort_bam,
                    self.targetbamfile[self.sampleID],
                    os.path.join(clair3workdir, "sorted_targets_exceeding.bam"),
                    self.threads,
                )

                if self.snp_calling:
                    self.SNPqueue.put(
                        [
                            run_list,
                            os.path.join(clair3workdir, "sorted_targets_exceeding.bam"),
                            os.path.join(
                                self.check_and_create_folder(
                                    self.output, self.sampleID
                                ),
                                "targets_exceeding_threshold.bed",
                            ),
                        ]
                    )

        # ToDo: Reinstate this line later.
        # self.update_target_coverage_table()

        # await asyncio.sleep(0.5)
        self.running = False

    def show_previous_data(self):
        if not self.browse:
            for item in app.storage.general[self.mainuuid]:
                if item == "sample_ids":
                    for sample in app.storage.general[self.mainuuid][item]:
                        self.sampleID = sample
            output = self.output
        if self.browse:
            output = self.check_and_create_folder(self.output, self.sampleID)
        # ToDo: This function needs to run in background threads.
        if self.check_file_time(os.path.join(output, "coverage_main.csv")):
            self.cov_df_main = pd.read_csv(os.path.join(output, "coverage_main.csv"))
            self.update_coverage_plot(self.cov_df_main)

        if self.check_file_time(os.path.join(output, "bed_coverage_main.csv")):
            self.bedcov_df_main = pd.read_csv(
                os.path.join(output, "bed_coverage_main.csv")
            )
            self.update_coverage_plot_targets(self.cov_df_main, self.bedcov_df_main)
            self.update_target_boxplot(self.bedcov_df_main)
        if self.check_file_time(os.path.join(output, "target_coverage.csv")):
            self.target_coverage_df = pd.read_csv(
                os.path.join(output, "target_coverage.csv")
            )
            self.update_target_coverage_table()
            self.update_coverage_time_plot()
            if self.summary:
                with self.summary:
                    self.summary.clear()
                    with ui.row():
                        ui.label("Coverage Depths - ")
                        if len(self.cov_df_main) > 0:
                            ui.label(
                                f"Global Estimated Coverage: {(self.cov_df_main['covbases'].sum()/self.cov_df_main['endpos'].sum()):.2f}x"
                            )
                            if len(self.bedcov_df_main) > 0:
                                ui.label(
                                    f"Targets Estimated Coverage: {(self.bedcov_df_main['bases'].sum()/self.bedcov_df_main['length'].sum()):.2f}x"
                                )
                                ui.label(
                                    f"Estimated enrichment: {(self.bedcov_df_main['bases'].sum()/self.bedcov_df_main['length'].sum())/(self.cov_df_main['covbases'].sum()/self.cov_df_main['endpos'].sum()):.2f}x"
                                )
                            else:
                                ui.label("Targets Estimated Coverage: Calculating....")
                        else:
                            ui.label("No data available")


        if self.check_file_time(f"{output}/clair3/snpsift_output.vcf.csv"):
            df = pd.read_csv(f"{output}/clair3/snpsift_output.vcf.csv", low_memory=False)

            # Define a function to check for 'pathogenic'
            def contains_pathogenic(text):
                if pd.isna(text):
                    return False
                return 'pathogenic' in str(text).lower()

            # Apply the function to each cell and filter rows
            #vcf = df[df.map(contains_pathogenic).any(axis=1)]
            vcf = df

            if len(vcf) > 0:
                self.SNPplaceholder.clear()
                with self.SNPplaceholder:
                    self.snptable = (
                        ui.table.from_pandas(vcf, pagination=25)
                        .props("dense")
                        .style("height: 900px")
                        .style("font-size: 100%; font-weight: 300")
                    )
                    for col in self.snptable.columns:
                        col["sortable"] = True

                    def toggle(column: dict, visible: bool) -> None:
                        column["classes"] = "" if visible else "hidden"
                        column["headerClasses"] = "" if visible else "hidden"
                        self.snptable.update()

                    def set_pathogenic(value: bool) -> None:
                        self.snptable.filter = "pathogenic" if value else None

                    with self.snptable.add_slot("top-left"):

                        def toggle_fs() -> None:
                            self.snptable.toggle_fullscreen()
                            button.props(
                                "icon=fullscreen_exit"
                                if self.snptable.is_fullscreen
                                else "icon=fullscreen"
                            )

                        button = ui.button(
                            "Toggle fullscreen",
                            icon="fullscreen",
                            on_click=toggle_fs,
                        ).props("flat")

                        with ui.button(icon="menu"):
                            with ui.menu(), ui.column().classes("gap-0 p-2"):
                                for column in self.snptable.columns:
                                    ui.switch(
                                        column["label"],
                                        value=True,
                                        on_change=lambda e, column=column: toggle(
                                            column, e.value
                                        ),
                                    )

                        ui.switch(
                            "Show potentially pathogenic SNPs only",
                            value=False,
                            on_change=lambda e: set_pathogenic(e.value),
                        )

                    with self.snptable.add_slot("top-right"):
                        with ui.input(placeholder="Search").props("type=search").bind_value(
                            self.snptable, "filter"
                        ).add_slot("append"):
                            ui.icon("search")

        if self.check_file_time(f"{output}/clair3/snpsift_indel_output.vcf.csv"):
            df = pd.read_csv(f"{output}/clair3/snpsift_indel_output.vcf.csv", low_memory=False)

            # Define a function to check for 'pathogenic'
            def contains_pathogenic(text):
                if pd.isna(text):
                    return False
                return 'pathogenic' in str(text).lower()

            # Apply the function to each cell and filter rows
            #vcfindel = df[df.map(contains_pathogenic).any(axis=1)]

            vcfindel = df

            if len(vcfindel) > 0:
                self.INDELplaceholder.clear()
                with self.INDELplaceholder:
                    self.indeltable = (
                        ui.table.from_pandas(vcfindel, pagination=25)
                        .props("dense")
                        .style("height: 900px")
                        .style("font-size: 100%; font-weight: 300")
                    )
                    for col in self.indeltable.columns:
                        col["sortable"] = True

                    def idtoggle(column: dict, visible: bool) -> None:
                        column["classes"] = "" if visible else "hidden"
                        column["headerClasses"] = "" if visible else "hidden"
                        self.indeltable.update()

                    def set_idpathogenic(value: bool) -> None:
                        self.indeltable.filter = "pathogenic" if value else None

                    with self.indeltable.add_slot("top-left"):

                        def toggle_idfs() -> None:
                            self.indeltable.toggle_fullscreen()
                            button.props(
                                "icon=fullscreen_exit"
                                if self.indeltable.is_fullscreen
                                else "icon=fullscreen"
                            )

                        button = ui.button(
                            "Toggle fullscreen",
                            icon="fullscreen",
                            on_click=toggle_idfs,
                        ).props("flat")

                        with ui.button(icon="menu"):
                            with ui.menu(), ui.column().classes("gap-0 p-2"):
                                for column in self.indeltable.columns:
                                    ui.switch(
                                        column["label"],
                                        value=True,
                                        on_change=lambda e, column=column: idtoggle(
                                            column, e.value
                                        ),
                                    )

                        ui.switch(
                            "Show potentially pathogenic IN/DELs only",
                            value=False,
                            on_change=lambda e: set_idpathogenic(e.value),
                        )

                    with self.indeltable.add_slot("top-right"):
                        with ui.input(placeholder="Search").props("type=search").bind_value(
                            self.indeltable, "filter"
                        ).add_slot("append"):
                            ui.icon("search")

def test_me(
    port: int,
    threads: int,
    watchfolder: str,
    output: str,
    reload: bool = False,
    browse: bool = False,
):
    my_connection = None
    with theme.frame("Target Coverage Data", my_connection):
        TestObject = TargetCoverage(threads, output, progress=True)
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
