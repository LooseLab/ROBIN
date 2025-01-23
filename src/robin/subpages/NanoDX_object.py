"""
NanoDX Analysis Module

This module provides functionality for analyzing NanoDX methylation data. The primary
components include running external tools to extract and process methylation sites from BAM files,
visualizing the results using a GUI, and managing the process asynchronously.

Dependencies:
    - pandas: Data manipulation and analysis.
    - os: Interaction with the operating system.
    - sys: System-specific parameters and functions.
    - time: Time-related functions.
    - nicegui: GUI creation.
    - pysam: BAM file processing.
    - shutil: File operations.
    - tempfile: Temporary file creation.
    - click: Command-line interface creation.
    - pathlib: File system paths.
    - typing: Type hinting.
    - robin: Custom modules for the specific workflow.

Modules:
    - subpages.base_analysis: BaseAnalysis class from robin.subpages.base_analysis.
    - theme: robin theme module.
    - resources: robin resources module.
    - NN_model: NN_classifier class from robin.submodules.nanoDX.workflow.scripts.NN_model.
    - merge_bedmethyl: Functions for merging, saving, and collapsing bed methylation data.

Environment Variables:
    - CI: Set to "1".

Functions:
    - run_modkit: Executes modkit to extract methylation data from a BAM file.
    - run_samtools_sort: Sorts BAM files using Samtools.
    - classification: Runs classification on the extracted data using a neural network model.

Classes:
    - NanoDX_object: Manages the NanoDX analysis process, including setting up the GUI and handling BAM file processing.

Command-line Interface:
    - run_main: CLI entry point for running the app, using Click for argument parsing.

Usage:
    The module can be run as a script to start the GUI for NanoDX analysis, specifying
    options like the port, number of threads, watch folder, and output directory.
"""

from __future__ import annotations
from robin.subpages.base_analysis import BaseAnalysis
from nicegui import ui, app, run
import time
import os
import sys
import click
from pathlib import Path
import numpy as np
import pysam
import pandas as pd
import shutil
import tempfile
from robin import models, theme, resources
import asyncio
import logging

# Use the main logger configured in the main application
logger = logging.getLogger(__name__)

# Code to edit the NN_model to use cpu and not gpu
edit_file_path = os.path.join(
    os.path.dirname(os.path.abspath(theme.__file__)),
    "submodules",
    "nanoDX",
    "workflow",
    "scripts",
    "NN_model.py",
)
if not os.path.exists(edit_file_path):
    raise FileNotFoundError(f"The file {edit_file_path} does not exist.")

with open(edit_file_path, "r") as file:
    lines = file.readlines()

# Check if the line matches the target content
line_number = 28
target_content = (
    'self.device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")'
)
if lines[line_number - 1].strip() == target_content.strip():
    target_substring = '("cuda:0" if torch.cuda.is_available() else "cpu")'
    new_substring = '("cpu")'
    lines[line_number - 1] = lines[line_number - 1].replace(
        target_substring, new_substring
    )
    with open(edit_file_path, "w") as file:
        file.writelines(lines)
    logger.info(f"Line {line_number} was modified.")
else:
    logger.info(f"Line {line_number} was not modified.")

from robin.submodules.nanoDX.workflow.scripts.NN_model import NN_classifier
from robin.utilities.merge_bedmethyl import (
    merge_bedmethyl,
    save_bedmethyl,
    collapse_bedmethyl,
)
from typing import List, Tuple, Optional, Dict, Any


def load_modkit_data(parquet_path):
    for attempt in range(5):  # Retry up to 5 times
        try:
            merged_modkit_df = pd.read_parquet(parquet_path)
            logger.debug("Successfully read the Parquet file.")
            break
        except Exception as e:
            logger.debug(f"Attempt {attempt+1}: File not ready ({e}). Retrying...")
            time.sleep(10)
    else:
        logger.debug("Failed to read Parquet file after multiple attempts.")
        return None

    column_names = [
        "chrom", "chromStart", "chromEnd", "mod_code", "score_bed", "strand",
        "thickStart", "thickEnd", "color", "valid_cov", "percent_modified",
        "n_mod", "n_canonical", "n_othermod", "n_delete", "n_fail",
        "n_diff", "n_nocall"
    ]

    # Keep only the original 18 columns
    return merged_modkit_df[column_names]

def run_modkit(cpgs: str, sortfile: str, temp: str, threads: int) -> None:
    """
    Executes modkit on a bam file and extracts the methylation data.

    Args:
        cpgs (str): Path to the CpG BED file.
        sortfile (str): Path to the sorted BAM file.
        temp (str): Path to the temporary output file.
        threads (int): Number of threads to use.
    """
    try:
        logger.debug(
            f"Running modkit with the following parameters: cpgs={cpgs}, sortfile={sortfile}, temp={temp}, threads={threads}"
        )

        if logger.isEnabledFor(logging.DEBUG):
            command = f"modkit pileup --include-bed {cpgs} --filter-threshold 0.73 --combine-mods --mixed-delim -t {threads} {sortfile} {temp}"
        else:
            command = f"modkit pileup --include-bed {cpgs} --filter-threshold 0.73 --combine-mods --mixed-delim -t {threads} {sortfile} {temp} --suppress-progress >/dev/null 2>&1"

        logger.debug(f"Executing command: {command}")
        os.system(command)
        logger.debug("modkit command executed successfully.")
    except Exception:
        logger.error("An error occurred while running modkit", exc_info=True)


def run_samtools_sort(
    file: str, tomerge: List[str], sortfile: str, threads: int
) -> None:
    """
    Sorts BAM files using Samtools.

    Args:
        file (str): Path to the output BAM file.
        tomerge (List[str]): List of BAM files to merge.
        sortfile (str): Path to the sorted BAM file.
        threads (int): Number of threads to use.
    """
    logger.debug(
        f"Running samtools sort with the following parameters: file={file}, tomerge={tomerge}, sortfile={sortfile}, threads={threads}"
    )
    pysam.cat("-o", file, *tomerge)
    pysam.sort("-@", f"{threads}", "--write-index", "-o", sortfile, file)
    logger.debug("samtools sort command executed successfully.")


def classification(
    modelfile: str, test_df: pd.DataFrame
) -> Tuple[np.ndarray, np.ndarray, int]:
    """
    Runs classification on the extracted data using a neural network model.

    Args:
        modelfile (str): Path to the neural network model file.
        test_df (pd.DataFrame): DataFrame containing the test data.

    Returns:
        Tuple[np.ndarray, np.ndarray, int]: Predictions, class labels, and number of features.
    """
    logger.debug(f"Running classification with model file: {modelfile}")
    NN = NN_classifier(modelfile)
    try:
        predictions, class_labels, n_features = NN.predict(test_df)
        logger.debug("Classification executed successfully.")
    except Exception:
        logger.error("An error occurred during classification", exc_info=True)
        test_df.to_csv("errordf.csv", sep=",", index=False, encoding="utf-8")
        # sys.exit(1)
    return predictions, class_labels, n_features


class NanoDX_object(BaseAnalysis):
    """
    NanoDX_object handles the NanoDX analysis process, including setting up the GUI
    and processing BAM files asynchronously.

    Attributes:
        cpgs_file (str): Path to the CpG BED file.
        cpgs (pd.DataFrame): DataFrame containing CpG data.
        model (str): Name of the neural network model file.
        threshold (float): Threshold for classification confidence.
        nanodx_bam_count (int): Counter for the number of BAM files processed.
        not_first_run (bool): Flag indicating whether it's the first run.
        modelfile (str): Path to the neural network model file.
        nanodx_df_store (pd.DataFrame): DataFrame for storing NanoDX results.
        nanodxfile (tempfile.NamedTemporaryFile): Temporary file for NanoDX results.
    """

    def __init__(self, *args, model: str = "Capper_et_al_NN.pkl", **kwargs) -> None:
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
        self.dataDir = {}
        self.bedDir = {}
        self.threshold = 0.05
        self.nanodx_bam_count = {}
        self.not_first_run = {}  # False
        self.modelfile = os.path.join(
            os.path.dirname(os.path.abspath(models.__file__)), self.model
        )
        logger.info(f"model file: {self.modelfile}")
        self.nanodx_df_store = {}  # pd.DataFrame()
        self.nanodxfile = {}
        self.merged_bed_file = {}
        super().__init__(*args, **kwargs)
        if self.model != "Capper_et_al_NN.pkl":
            self.storefile = "PanNanoDX_scores.csv"
        else:
            self.storefile = "NanoDX_scores.csv"
        #    self.output = f"{self.output}_PanCan"

    def __del__(self):
        if self.nanodxfile:
            pass
            #self.nanodxfile.close()

    def setup_ui(self) -> None:
        """
        Sets up the user interface for the NanoDX analysis.
        """
        with ui.card().classes("w-full p-2"):
            with ui.grid(columns=8).classes("w-full h-auto gap-2"):
                with ui.card().classes(
                    f"min-[{self.MENU_BREAKPOINT+1}px]:col-span-3 max-[{self.MENU_BREAKPOINT}px]:col-span-8 shadow-lg rounded-lg p-2"
                ):
                    self.create_nanodx_chart("NanoDX")
                with ui.card().classes(
                    f"min-[{self.MENU_BREAKPOINT+1}px]:col-span-5 max-[{self.MENU_BREAKPOINT}px]:col-span-8 shadow-lg rounded-lg p-2"
                ):
                    self.create_nanodx_time_chart("NanoDX Time Series")
        if self.summary:
            with self.summary:
                ui.label(f"NanoDX classification {self.model}: Unknown")
        if self.browse:
            self.show_previous_data()
        else:
            ui.timer(5, lambda: self.show_previous_data())

    def show_previous_data(self):
        """
        Load and display previously generated NanoDX analysis results.
        """
        if not self.browse:
            for item in app.storage.general[self.mainuuid]:
                if item == "sample_ids":
                    for sample in app.storage.general[self.mainuuid][item]:
                        self.sampleID = sample
            output = self.output
        if self.browse:
            output = self.check_and_create_folder(self.output, self.sampleID)

        if self.check_file_time(os.path.join(output, self.storefile)):
            self.nanodx_df_store = pd.read_csv(
                os.path.join(output, self.storefile),
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
            
            # Update summary with new card
            if self.summary:
                with self.summary:
                    self.summary.clear()
                    classification_text = f"NanoDX classification: {lastrow_plot_top.index[0]}"
                    self.create_summary_card(
                        classification_text=classification_text,
                        confidence_value=lastrow_plot_top.values[0],
                        model_name=self.model,
                        features_found=int(self.nanodx_df_store.iloc[-1]["number_probes"])
                    )
            
            self.update_nanodx_plot(
                lastrow_plot.index.to_list(),
                list(lastrow_plot.values),
                "All",
                self.nanodx_df_store.iloc[-1]["number_probes"],
            )

    async def process_bam(self, bamfile: List[Tuple[str, float]]) -> None:
        """
        Process BAM files and perform NanoDX analysis.

        This method handles the complete workflow of processing BAM files:
        1. Merging BAM files
        2. Extracting methylation data using modkit
        3. Converting to BED format
        4. Running Sturgeon predictions
        5. Updating visualizations

        Parameters
        ----------
        bamfile : List[Tuple[str, float]]
            List of tuples containing BAM file paths and their timestamps

        Notes
        -----
        The method processes files in batches and updates progress trackers.
        Results are stored in temporary directories and visualized in real-time.
        """
        sampleID = self.sampleID
        # Initialize directories for each sampleID if not already present
        if sampleID not in self.dataDir.keys():
            self.dataDir[sampleID] = tempfile.TemporaryDirectory(
                dir=self.check_and_create_folder(self.output, sampleID)
            )
            self.bedDir[sampleID] = tempfile.TemporaryDirectory(
                dir=self.check_and_create_folder(self.output, sampleID)
            )
        tomerge = []
        latest_file = 0
        if app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                    "bam_count"
                ] > 0:
            if latest_file:
                currenttime = latest_file * 1000
            else:
                currenttime = time.time() * 1000
                
            parquet_path = os.path.join(self.check_and_create_folder(self.output, sampleID), f"{sampleID}.parquet")
            
            try:
                if self.check_file_time(parquet_path):
                    tomerge_length_file = os.path.join(self.check_and_create_folder(self.output, sampleID), "tomerge_length.txt")
                    with open(tomerge_length_file, "r") as f:
                        tomerge_length = int(f.readline().strip().split(": ")[1])
                    
                    merged_modkit_df = await run.cpu_bound(load_modkit_data, parquet_path)
                    
                    merged_modkit_df.rename(columns={"chromStart":"start_pos", "chromEnd":"end_pos", "mod_code":"mod", "thickStart":"start_pos2", "thickEnd":"end_pos2", "color":"colour"}, inplace=True)
                    merged_modkit_df.rename(columns={'n_canonical':'Ncanon', 'n_delete':'Ndel', 'n_diff':'Ndiff', 'n_fail':'Nfail', 'n_mod':'Nmod', 'n_nocall':'Nnocall', 'n_othermod':'Nother', 'valid_cov':'Nvalid', 'percent_modified':'score'}, inplace=True)
                    
                    nanodx_df= await run.cpu_bound(collapse_bedmethyl, merged_modkit_df)
                    test_df = pd.merge(
                        nanodx_df,
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
                    
                    nanoDX_df = pd.DataFrame({"class": class_labels, "score": predictions})
                    nanoDX_save = nanoDX_df.set_index("class").T
                    nanoDX_save["number_probes"] = n_features
                    nanoDX_save["timestamp"] = currenttime

                    if sampleID not in self.nanodx_df_store.keys():
                        self.nanodx_df_store[sampleID] = pd.DataFrame()
                    self.nanodx_df_store[sampleID] = pd.concat(
                        [self.nanodx_df_store[sampleID], nanoDX_save.set_index("timestamp")]
                    )

                    self.nanodx_df_store[sampleID].to_csv(
                        os.path.join(
                            self.check_and_create_folder(self.output, sampleID),
                            self.storefile,
                        )
                    )
                    
                    app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                        "bam_processed"
                    ] = tomerge_length

            except Exception as e:
                print(e)

        self.running = False


    def create_nanodx_chart(self, title):
        """
        Create a bar chart for displaying NanoDX classification results.
        """
        self.nanodxchart = self.create_chart(title)
        self.nanodxchart.options.update({
            "backgroundColor": "transparent",
            "title": {
                "text": title,
                "left": "center",
                "top": 10,
                "textStyle": {
                    "fontSize": 16,
                    "fontWeight": "normal",
                    "color": "#000000"
                },
                "subtextStyle": {
                    "fontSize": 12,
                    "color": "#666666",
                    "align": "center"
                }
            },
            "tooltip": {
                "trigger": "axis",
                "axisPointer": {"type": "shadow"},
                "formatter": "{b}: {c}%",
                "textStyle": {"fontSize": 14}
            },
            "grid": {
                "left": "5%",
                "right": "5%",
                "bottom": "5%",
                "top": "25%",
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
                "barMaxWidth": "60%",
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
            }]
        })

    def create_nanodx_time_chart(self, title):
        """
        Create a time series chart for NanoDX results.
        """
        self.nanodx_time_chart = self.create_time_chart(title)
        self.nanodx_time_chart.options.update({
            "backgroundColor": "transparent",
            "title": {
                "text": title,
                "left": "center",
                "top": 5,
                "textStyle": {
                    "fontSize": 16,
                    "fontWeight": "normal",
                    "color": "#000000"
                },
                "padding": [0, 0, 20, 0]
            },
            "tooltip": {
                "trigger": "axis",
                "axisPointer": {"type": "line"},
                "textStyle": {"fontSize": 14}
            },
            "grid": {
                "left": "5%",
                "right": "5%",
                "bottom": "5%",
                "top": "25%",
                "containLabel": True
            },
            "legend": {
                "type": "scroll",
                "orient": "horizontal",
                "top": 45,
                "width": "90%",
                "left": "center",
                "textStyle": {
                    "fontSize": 12,
                    "color": "#666666"
                },
                "pageButtonPosition": "end",
                "pageButtonGap": 5,
                "pageButtonItemGap": 5,
                "pageIconColor": "#666666",
                "pageIconInactiveColor": "#aaa",
                "pageIconSize": 12,
                "pageTextStyle": {
                    "color": "#666666"
                },
                "itemGap": 25,
                "itemWidth": 14,
                "itemHeight": 14,
                "selectedMode": True
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
            }
        })

    def update_nanodx_plot(self, x, y, count, n_features):
        """
        Update the NanoDX visualization plot with new data.

        Parameters
        ----------
        x : List[str]
            List of tumor types
        y : List[float]
            Confidence scores for each tumor type
        count : str
            Number of BAM files processed
        n_features : int
            Number of features detected
        """
        # Convert values to percentages and format
        formatted_values = [float(f"{val * 100:.1f}") for val in y]
        
        # Sort the data in descending order
        sorted_indices = sorted(range(len(formatted_values)), key=lambda k: formatted_values[k], reverse=True)
        sorted_values = [formatted_values[i] for i in sorted_indices]
        sorted_labels = [x[i] for i in sorted_indices]
        
        # Create descriptive title with key information
        title_text = (
            f"NanoDX Analysis Results\n"
            f"{count} samples processed • {int(n_features)} features found"
        )
        
        self.nanodxchart.options["title"].update({
            "text": title_text,
            "subtext": f"Model: {self.model}"
        })
        self.nanodxchart.options["yAxis"]["data"] = sorted_labels
        self.nanodxchart.options["series"][0].update({
            "data": sorted_values,
            "itemStyle": {
                "color": "#007AFF",
                "borderRadius": [0, 4, 4, 0]
            }
        })
        self.nanodxchart.update()

    def update_nanodx_time_chart(self, datadf):
        """
        Update the time series chart with new data.

        Parameters
        ----------
        datadf : pd.DataFrame
            DataFrame containing time series data for visualization
        """
        try:
            logger.debug(f"Updating NanoDX time chart with DataFrame shape: {datadf.shape}")
            logger.debug(f"DataFrame columns: {datadf.columns.tolist()}")

            if datadf.empty:
                logger.warning("No data available for NanoDX time series chart")
                self.nanodx_time_chart.options["title"]["text"] = "No Data Available"
                self.nanodx_time_chart.options["series"] = []
                self.nanodx_time_chart.update()
                return

            self.nanodx_time_chart.options["series"] = []
            
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
                    # Convert values to percentages
                    data_list = [[key, float(f"{value * 100:.1f}")] for key, value in data.items()]
                    logger.debug(f"First few data points for {series}: {data_list[:3]}")
                    
                    self.nanodx_time_chart.options["series"].append({
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
            if not datadf.empty:
                latest_data = datadf.iloc[-1]
                if "number_probes" in latest_data:
                    latest_data = latest_data.drop("number_probes")  # Remove number_probes from latest data
                
                if not latest_data.empty and not latest_data.isna().all():
                    max_confidence = latest_data.max() * 100  # Convert to percentage
                    max_type = latest_data.idxmax()
                    self.nanodx_time_chart.options["title"]["text"] = (
                        f"Classification Confidence Over Time\n"
                        f"Current highest confidence: {max_type} ({max_confidence:.1f}%)"
                    )
                else:
                    self.nanodx_time_chart.options["title"]["text"] = "Classification Confidence Over Time"
            
            self.nanodx_time_chart.update()
            
        except Exception as e:
            logger.error(f"Error updating NanoDX time chart: {str(e)}", exc_info=True)
            self.nanodx_time_chart.options["title"]["text"] = "Error Updating Chart"
            self.nanodx_time_chart.options["series"] = []
            self.nanodx_time_chart.update()


def test_me(
    port: int,
    threads: int,
    watchfolder: str,
    output: str,
    reload: bool = False,
    browse: bool = False,
) -> None:
    """
    Sets up and runs the NanoDX analysis application.

    Args:
        port (int): Port number for the server.
        threads (int): Number of threads to use for processing.
        watchfolder (str): Path to the folder to watch for new BAM files.
        output (str): Path to the output directory.
        reload (bool): Flag to reload the application on changes.
        browse (bool): Flag to enable browsing historic data.
    """
    my_connection: Optional[Dict[str, Any]] = None
    with theme.frame("Target Coverage Data", my_connection):
        TestObject = NanoDX_object(threads, output, progress=True, batch=True)
    if not browse:
        path = watchfolder
        searchdirectory = os.fsencode(path)
        for root, d_names, f_names in os.walk(searchdirectory):
            directory = os.fsdecode(root)
            for f in f_names:
                filename = os.fsdecode(f)
                if filename.endswith(".bam"):
                    TestObject.add_bam(os.path.join(directory, filename))
    else:
        TestObject.progress_trackers.visible = False
        TestObject.show_previous_data(output)
    ui.run(port=port, reload=reload)


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
def run_main(
    port: int, threads: int, watchfolder: str, output: str, browse: bool
) -> None:
    """
    CLI entry point for running the NanoDX analysis app.

    Args:
        port (int): The port to serve the app on.
        threads (int): Number of threads available for processing.
        watchfolder (str): Directory to watch for new BAM files.
        output (str): Directory to save output files.
        browse (bool): Enable browsing historic data.
    """
    if browse:
        click.echo("Browse mode is enabled. Only the output folder is required.")
        test_me(
            port=port,
            reload=False,
            threads=threads,
            watchfolder=None,
            output=watchfolder,
            browse=browse,
        )
    else:
        click.echo(f"Watchfolder: {watchfolder}, Output: {output}")
        if watchfolder is None or output is None:
            click.echo("Watchfolder and output are required when --browse is not set.")
            sys.exit(1)
        test_me(
            port=port,
            reload=False,
            threads=threads,
            watchfolder=watchfolder,
            output=output,
            browse=browse,
        )


if __name__ in {"__main__", "__mp_main__"}:
    print("GUI launched by auto-reload function.")
    run_main()
