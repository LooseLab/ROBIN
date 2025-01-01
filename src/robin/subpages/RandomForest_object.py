"""
A module for analyzing and visualizing Random Forest methylation classification results.

This module provides functionality for processing BAM files to analyze methylation patterns
using a Random Forest classifier. It includes capabilities for:

- Real-time BAM file processing
- Methylation data extraction using modkit
- Random Forest classification
- Time series analysis of classification confidence

The module integrates with the NiceGUI framework for interactive visualization
and uses temporary files for efficient data processing.

Dependencies
-----------
- pysam
- pandas
- nicegui
- R (with required packages)
- modkit

Notes
-----
The module requires proper configuration of input/output directories and
assumes the presence of necessary R scripts and model files.
"""

from robin.subpages.base_analysis import BaseAnalysis
import os
import tempfile
import time
import pandas as pd
from nicegui import ui, app, run
from robin import theme, resources
import pysam
import logging
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

# Use the main logger configured in the main application
logger = logging.getLogger(__name__)

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
    logger.debug(command)
    # print(command)

    os.system(command)


def run_samtools_sort(file, tomerge, sortfile, threads, regions):
    pysam.cat("-o", file, *tomerge)
    # pysam.sort("-@", f"{threads}", "--write-index", "-o", sortfile, file)
    intermediate_bam = tempfile.NamedTemporaryFile(suffix=".bam")
    command = (
        f"samtools sort -@{threads} --write-index -o {intermediate_bam.name} {file}"
    )
    logger.debug(command)
    os.system(command)
    command2 = f"samtools view -b -L {regions} -o {sortfile} --write-index {intermediate_bam.name} "
    logger.debug(command2)
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
        logger.debug(command)
        os.system(command)
        # self.log("Done processing bam file")
    except Exception as e:
        print(e)
        # self.log(e)
        pass


class RandomForest_object(BaseAnalysis):
    """
    A class for processing and visualizing Random Forest methylation classification results.

    This class extends BaseAnalysis to provide specialized functionality for
    Random Forest methylation analysis. It handles real-time processing of BAM files,
    methylation data extraction, and visualization of classification results.

    Parameters
    ----------
    *args
        Variable length argument list passed to BaseAnalysis
    showerrors : bool, optional
        Whether to show error messages (default: False)
    **kwargs
        Arbitrary keyword arguments passed to BaseAnalysis

    Attributes
    ----------
    rcns2_df_store : dict
        Store for Random Forest dataframes
    threshold : float
        Confidence threshold for predictions (default: 0.5)
    bambatch : dict
        Dictionary tracking BAM file batches
    cpgs_file : str
        Path to the CpG sites BED file
    offset : bool
        Flag for offset calculations
    first_run : dict
        Tracks first run status for each sample
    showerrors : bool
        Whether to show error messages
    modelfile : str
        Path to the model file
    dataDir : dict
        Dictionary of temporary directories for data
    bedDir : dict
        Dictionary of temporary directories for BED files
    merged_bed_file : dict
        Dictionary of merged BED files
    """
    def __init__(self, *args, showerrors=False, **kwargs):
        self.rcns2_df_store = {}
        self.threshold = 0.5
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
        """
        Load and display previously generated Random Forest analysis results.
        """
        try:
            logger.info("Starting show_previous_data")
            
            # Get output path
            if not self.browse:
                for item in app.storage.general[self.mainuuid]:
                    if item == "sample_ids":
                        for sample in app.storage.general[self.mainuuid][item]:
                            self.sampleID = sample
                output = self.output
                logger.debug(f"Using output path: {output}")
            if self.browse:
                output = self.check_and_create_folder(self.output, self.sampleID)
                logger.debug(f"Using browse output path: {output}")

            scores_file = os.path.join(output, "random_forest_scores.csv")
            logger.info(f"Looking for scores file: {scores_file}")

            if os.path.exists(scores_file):
                logger.info(f"Loading scores from: {scores_file}")
                try:
                    self.rcns2_df_store = pd.read_csv(scores_file, index_col=0)
                    logger.info(f"DataFrame loaded with shape: {self.rcns2_df_store.shape}")
                    
                    if not self.rcns2_df_store.empty:
                        columns_greater_than_threshold = (self.rcns2_df_store > self.threshold).any()
                        columns_not_greater_than_threshold = ~columns_greater_than_threshold
                        result = self.rcns2_df_store.columns[columns_not_greater_than_threshold].tolist()
                        logger.debug(f"Filtered columns: {result}")
                        
                        # Update time series chart
                        filtered_df = self.rcns2_df_store.drop(columns=result)
                        logger.info(f"Updating time chart with filtered data shape: {filtered_df.shape}")
                        self.update_rcns2_time_chart(filtered_df)
                        
                        # Get last row data
                        lastrow = self.rcns2_df_store.iloc[-1]
                        logger.debug(f"Last row before processing: {lastrow.to_dict()}")
                        
                        n_features = lastrow.get("number_probes", 0)
                        if "number_probes" in lastrow.index:
                            lastrow = lastrow.drop("number_probes")
                            logger.info(f"Dropped number_probes column, features found: {n_features}")
                        
                        lastrow_plot = lastrow.sort_values(ascending=False).head(10)
                        lastrow_plot_top = lastrow.sort_values(ascending=False).head(1)
                        logger.info(f"Top classification: {lastrow_plot_top.index[0]} ({lastrow_plot_top.values[0]:.1f}%)")
                        
                        # Update summary card
                        if self.summary:
                            with self.summary:
                                self.summary.clear()
                                classification_text = f"Forest classification: {lastrow_plot_top.index[0]}"
                                logger.debug(f"Creating summary card with text: {classification_text}")
                                self.create_summary_card(
                                    classification_text=classification_text,
                                    confidence_value=lastrow_plot_top.values[0] / 100,
                                    features_found=int(n_features) if n_features else None
                                )
                        
                        # Update bar plot
                        logger.info("Updating bar plot with top 10 classifications")
                        self.update_rcns2_plot(
                            lastrow_plot.index.to_list(),
                            list(lastrow_plot.values),
                            "All",
                            n_features if n_features else None
                        )
                        logger.info("Completed show_previous_data successfully")
                except Exception as e:
                    logger.error(f"Error reading scores file: {str(e)}", exc_info=True)
            else:
                logger.debug(f"No scores file found at: {scores_file}")
                
        except Exception as e:
            logger.error(f"Error in show_previous_data: {str(e)}", exc_info=True)
            raise

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
            (file, filetime) = bamfile.pop(0)
            if filetime > latest_file:
                latest_file = filetime
            tomerge.append(file)
            app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                "bams_in_processing"
            ] += 1
            if len(tomerge) > 100:
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
        """
        Create a bar chart for displaying Random Forest classification results.
        """
        self.echart = self.create_chart(title)
        self.echart.options.update({
            "backgroundColor": "transparent",
            "title": {
                "text": title,
                "left": "center",
                "top": 20,  # Increased for consistency
                "textStyle": {
                    "fontSize": 16,
                    "fontWeight": "normal",
                    "color": "#000000"  # Explicit color for better contrast
                }
            },
            "tooltip": {
                "trigger": "axis",
                "axisPointer": {"type": "shadow"},
                "formatter": "{b}: {c}%",
                "textStyle": {"fontSize": 14}
            },
            "grid": {
                "left": "15%",  # Standardized margin
                "right": "10%",
                "bottom": "10%",
                "top": "25%",  # Increased for title and legend
                "containLabel": True
            },
            "xAxis": {
                "type": "value",
                "min": 0,
                "max": 100,
                "interval": 20,
                "axisLabel": {
                    "fontSize": 12,
                    "formatter": "{value}%",
                    "color": "#666666"  # Subtle color for axis labels
                },
                "splitLine": {
                    "show": True,
                    "lineStyle": {
                        "color": "#E0E0E0",
                        "type": "dashed"
                    }
                }
            },
            "yAxis": {
                "type": "category",
                "inverse": True,
                "data": [],
                "axisLabel": {
                    "fontSize": 12,
                    "width": 250,
                    "overflow": "break",
                    "interval": 0,
                    "align": "right",
                    "color": "#666666"
                }
            },
            "series": [{
                "type": "bar",
                "name": "Confidence",
                "barMaxWidth": "60%",  # Standardized bar width
                "itemStyle": {
                    "color": "#007AFF",
                    "borderRadius": [0, 4, 4, 0]
                },
                "label": {
                    "show": True,
                    "position": "right",
                    "formatter": "{c}%",
                    "fontSize": 12,
                    "color": "#666666"
                },
                "data": []
            }],
            "aria": {
                "enabled": True,
                "decal": {
                    "show": True
                }
            }
        })

    def create_rcns2_time_chart(self, title):
        """
        Create a time series chart for Random Forest results.
        """
        self.rcns2_time_chart = self.create_time_chart(title)
        self.rcns2_time_chart.options.update({
            "backgroundColor": "transparent",
            "title": {
                "text": title,
                "left": "center",
                "top": 20,
                "textStyle": {
                    "fontSize": 16,
                    "fontWeight": "normal",
                    "color": "#000000"
                }
            },
            "tooltip": {
                "trigger": "axis",
                "axisPointer": {"type": "line"},
                "textStyle": {"fontSize": 14}
            },
            "grid": {
                "left": "15%",
                "right": "15%",
                "bottom": "10%",
                "top": "25%",
                "containLabel": True
            },
            "legend": {
                "type": "scroll",
                "orient": "horizontal",
                "top": 50,
                "textStyle": {
                    "fontSize": 12,
                    "color": "#666666"
                }
            },
            "xAxis": {
                "type": "time",
                "axisLabel": {
                    "fontSize": 12,
                    "formatter": "{HH}:{mm}",
                    "color": "#666666"
                },
                "splitLine": {
                    "show": True,
                    "lineStyle": {
                        "color": "#E0E0E0",
                        "type": "dashed"
                    }
                }
            },
            "yAxis": {
                "type": "value",
                "min": 0,
                "max": 100,
                "interval": 20,
                "axisLabel": {
                    "fontSize": 12,
                    "formatter": "{value}%",
                    "color": "#666666"
                },
                "splitLine": {
                    "show": True,
                    "lineStyle": {
                        "color": "#E0E0E0",
                        "type": "dashed"
                    }
                }
            },
            "aria": {
                "enabled": True,
                "decal": {
                    "show": True
                }
            }
        })

    def update_rcns2_plot(self, x, y, count, n_features=None):
        """
        Update the Random Forest visualization plot with new data.

        Parameters
        ----------
        x : List[str]
            List of tumor types
        y : List[float]
            Confidence scores for each tumor type (already in percentages)
        count : str
            Number of BAM files processed
        n_features : int, optional
            Number of features detected
        """
        try:
            logger.debug(f"Updating plot with {len(x)} categories")
            logger.debug(f"Input values - x: {x}, y: {y}")
            
            # Values are already percentages, just format them
            formatted_values = [float(f"{val:.1f}") for val in y]
            logger.debug(f"Formatted values: {formatted_values}")
            
            # Sort the data in descending order
            sorted_indices = sorted(range(len(formatted_values)), key=lambda k: formatted_values[k], reverse=True)
            sorted_values = [formatted_values[i] for i in sorted_indices]
            sorted_labels = [x[i] for i in sorted_indices]
            logger.debug(f"Sorted values: {sorted_values}")
            logger.debug(f"Sorted labels: {sorted_labels}")
            
            # Create descriptive title with key information
            title_text = f"Random Forest Analysis Results\n{count} samples processed"
            if n_features:
                title_text += f" • {int(n_features)} features found"
            logger.debug(f"Setting title: {title_text}")
            
            self.echart.options["title"]["text"] = title_text
            self.echart.options["yAxis"]["data"] = sorted_labels
            self.echart.options["series"][0].update({
                "data": sorted_values,
                "itemStyle": {
                    "color": "#007AFF",  # iOS blue
                    "borderRadius": [0, 4, 4, 0]
                }
            })
            self.echart.update()
            logger.debug("Plot updated successfully")
            
        except Exception as e:
            logger.error(f"Error in update_rcns2_plot: {str(e)}", exc_info=True)
            raise

    def update_rcns2_time_chart(self, datadf):
        """
        Update the time series chart with new data.

        Parameters
        ----------
        datadf : pd.DataFrame
            DataFrame containing time series data for visualization (values already in percentages)
        """
        try:
            logger.debug(f"Updating time chart with DataFrame shape: {datadf.shape}")
            logger.debug(f"DataFrame columns: {datadf.columns.tolist()}")
            
            self.rcns2_time_chart.options["series"] = []
            
            # iOS color palette for multiple series
            colors = [
                "#007AFF",  # Blue
                "#34C759",  # Green
                "#FF9500",  # Orange
                "#FF2D55",  # Red
                "#5856D6",  # Purple
                "#FF3B30",  # Red-Orange
                "#5AC8FA",  # Light Blue
                "#4CD964",  # Light Green
            ]
            
            for idx, (series, data) in enumerate(datadf.to_dict().items()):
                if series != "number_probes":
                    logger.debug(f"Processing series: {series}")
                    # Values are already percentages, just format them
                    data_list = [[key, float(f"{value:.1f}")] for key, value in data.items()]
                    logger.debug(f"First few data points for {series}: {data_list[:3]}")
                    
                    self.rcns2_time_chart.options["series"].append({
                        "name": series,
                        "type": "line",
                        "smooth": True,
                        "animation": False,
                        "symbolSize": 6,
                        "emphasis": {
                            "focus": "series",
                            "itemStyle": {
                                "borderWidth": 2
                            }
                        },
                        "endLabel": {
                            "show": True,
                            "formatter": "{a}: {c}%",
                            "distance": 10,
                            "fontSize": 12
                        },
                        "lineStyle": {
                            "width": 2,
                            "color": colors[idx % len(colors)]
                        },
                        "itemStyle": {
                            "color": colors[idx % len(colors)]
                        },
                        "data": data_list
                    })
            
            # Update chart title with summary
            latest_data = datadf.iloc[-1]  # Get latest data
            logger.debug(f"Latest data: {latest_data.to_dict()}")
            
            if "number_probes" in latest_data:
                latest_data = latest_data.drop("number_probes")
            
            max_confidence = latest_data.max()  # Get maximum confidence (already percentage)
            max_type = latest_data.idxmax()  # Get type with maximum confidence
            logger.debug(f"Max confidence: {max_confidence:.1f}% for type: {max_type}")
            
            self.rcns2_time_chart.options["title"]["text"] = (
                f"Classification Confidence Over Time\n"
                f"Current highest confidence: {max_type} ({max_confidence:.1f}%)"
            )
            
            self.rcns2_time_chart.update()
            logger.debug("Time chart updated successfully")
            
        except Exception as e:
            logger.error(f"Error in update_rcns2_time_chart: {str(e)}", exc_info=True)
            raise


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
