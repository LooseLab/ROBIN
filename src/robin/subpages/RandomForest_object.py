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

from robin.subpages.base_analysis import BaseAnalysis, BaseVis
import os
import tempfile
import time
import pandas as pd
from nicegui import ui, app, run, background_tasks
from robin import resources
import logging
from robin import models

from robin import submodules

from robin.utilities.merge_bedmethyl import (
    collapse_bedmethyl,
    load_minimal_modkit_data,
    reconstruct_full_bedmethyl_data,
)

from typing import List, Tuple
from robin.core.state import state, ProcessState

# Use the main logger configured in the main application
logger = logging.getLogger(__name__)

HVPATH = os.path.join(
    os.path.dirname(os.path.abspath(submodules.__file__)), "hv_rapidCNS2"
)


def run_rcns2(rcns2folder, batch, bed, threads, showerrors):
    """
    Run the Random Forest R script on the methylation data.

    Parameters
    ----------
    rcns2folder : str
        Directory for R script output
    batch : int
        Batch number for output file naming
    bed : str
        Path to input BED file
    threads : int
        Number of threads to use
    showerrors : bool
        Whether to show error messages
    """
    try:
        logger.info(f"Starting run_rcns2 with bed file: {bed}")
        logger.info(f"Output directory: {rcns2folder}")
        logger.info(f"Batch number: {batch}")

        # Check if R script exists
        r_script_path = f"{HVPATH}/bin/methylation_classification_nanodx_v0.2.R"
        if not os.path.exists(r_script_path):
            logger.error(f"R script not found at: {r_script_path}")
            raise FileNotFoundError(f"R script not found: {r_script_path}")

        # Check if other required files exist
        required_files = {
            "probes": f"{HVPATH}/bin/top_probes_hm450.Rdata",
            "training_data": f"{HVPATH}/bin/capper_top_100k_betas_binarised.Rdata",
            "array_file": f"{HVPATH}/bin/HM450.hg38.manifest.gencode.v22.Rdata",
        }

        for file_type, file_path in required_files.items():
            if not os.path.exists(file_path):
                logger.error(f"{file_type} file not found at: {file_path}")
                raise FileNotFoundError(f"{file_type} file not found: {file_path}")

        # Run the R script
        command = (
            f"Rscript {r_script_path} -s "
            + f"live_{batch} -o {rcns2folder} -i {bed} "
            + f"-p {HVPATH}/bin/top_probes_hm450.Rdata "
            + f"--training_data {HVPATH}/bin/capper_top_100k_betas_binarised.Rdata "
            + f"--array_file {HVPATH}/bin/HM450.hg38.manifest.gencode.v22.Rdata "
            + f"-t {threads} "
        )

        logger.info(f"Executing R command: {command}")

        # Execute command and capture output
        import subprocess

        try:
            result = subprocess.run(
                command.split(), capture_output=True, text=True, check=True
            )
            logger.info("R script executed successfully")
            logger.debug(f"R script stdout: {result.stdout}")
            # print(f"R script stdout: {result.stdout}")
            if result.stderr:
                pass
                # logger.warning(f"R script stderr: {result.stderr}")
        except subprocess.CalledProcessError as e:
            logger.error(f"R script failed with return code {e.returncode}")
            logger.error(f"stdout: {e.stdout}")
            logger.error(f"stderr: {e.stderr}")
            raise
            pass

        # Check if output file was created
        expected_output = f"{rcns2folder}/live_{batch}_votes.tsv"
        if os.path.exists(expected_output):
            logger.info(f"Output file created successfully: {expected_output}")
        else:
            logger.error(f"Expected output file not found: {expected_output}")
            raise FileNotFoundError(
                f"R script did not create expected output file: {expected_output}"
            )

    except Exception as e:
        logger.error(f"Error in run_rcns2: {str(e)}", exc_info=True)
        raise


def load_modkit_data(parquet_path):
    """
    Load minimal bedmethyl data for RandomForest analysis.

    This function loads only the essential columns needed for RandomForest classification:
    - chrom: Chromosome name
    - chromStart: Start position
    - percent_modified: Primary methylation data
    - mod_code: Modification code
    - strand: Strand information

    Args:
        parquet_path (str): Path to the parquet file containing bedmethyl data

    Returns:
        pd.DataFrame: DataFrame with minimal columns, sorted by chrom and chromStart
    """
    return load_minimal_modkit_data(parquet_path)


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
        super().__init__(*args, **kwargs)
        # Remove state tracking for Random Forest Analysis
        # state.start_process("Random Forest Analysis", ProcessType.BATCH)
        state.set_process_state("Random Forest Analysis", ProcessState.WAITING_FOR_DATA)
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

        # Set RandomForest-specific confidence thresholds
        # RandomForest confidence is reported as percentage (0-100) but converted to 0-1 in create_summary_card
        kwargs["high_confidence_threshold"] = (
            0.85  # RandomForest-specific high confidence threshold
        )
        kwargs["medium_confidence_threshold"] = (
            0.65  # RandomForest-specific medium confidence threshold
        )

    async def process_bam(
        self, bamfile: List[Tuple[str, float]], timestamp: float = None
    ) -> None:
        """
        Process BAM files and perform Random Forest analysis.

        Parameters
        ----------
        bamfile : List[Tuple[str, float]]
            List of tuples containing BAM file paths and their timestamps
            Each tuple contains (file_path, timestamp)
        timestamp : float, optional
            Optional timestamp override for the current processing batch
        """
        state.set_process_state("Random Forest Analysis", ProcessState.RUNNING)
        if not self.parquetqueue.empty():
            num_bam_files_seen = 0
            while not self.parquetqueue.empty():
                parquet_path, sampleID, file_count = self.parquetqueue.get_nowait()
                num_bam_files_seen += file_count
                
            if timestamp:
                currenttime = timestamp * 1000
            else:
                currenttime = timestamp * 1000 if timestamp else time.time() * 1000

            parquet_path = os.path.join(
                self.check_and_create_folder(self.output, sampleID),
                f"{sampleID}.parquet",
            )
            
            async def forest_bam_background_work(sampleID, parquet_path):
                try:
                    if self.check_file_time(parquet_path):
                        logger.debug("Parquet file exists and is ready for processing")
                        tomerge_length_file = os.path.join(
                            self.check_and_create_folder(self.output, sampleID),
                            "tomerge_length.txt",
                        )
                        try:
                            with open(tomerge_length_file, "r") as f:
                                tomerge_length = int(
                                    f.readline().strip().split(": ")[1]
                                )
                            logger.info(f"Number of files to merge: {tomerge_length}")
                        except FileNotFoundError:
                            logger.warning(
                                "tomerge_length.txt not found yet, waiting for file creation"
                            )
                            return
                        except Exception as e:
                            logger.error(f"Error reading tomerge_length.txt: {str(e)}")
                            return

                        logger.debug("Loading modkit data from parquet file...")
                        merged_modkit_df = await run.cpu_bound(
                            load_modkit_data, parquet_path
                        )

                        # For RandomForest, we need to reconstruct the full data structure
                        # since the R script expects all columns. We'll use the minimal data
                        # but need to add the missing columns for compatibility
                        full_modkit_df = await run.cpu_bound(
                            reconstruct_full_bedmethyl_data, merged_modkit_df
                        )

                        forest_dx = await run.cpu_bound(
                            collapse_bedmethyl, full_modkit_df
                        )

                        merged_modkit_df = forest_dx
                        if merged_modkit_df is None:
                            logger.error("Failed to load modkit data from parquet file")
                            return
                        logger.info(
                            f"Loaded modkit data with shape: {merged_modkit_df.shape}"
                        )

                        # Create a temporary directory for R script output
                        rcns2folder = tempfile.mkdtemp(
                            dir=self.check_and_create_folder(self.output, sampleID)
                        )
                        logger.info(
                            f"Created temporary directory for R output: {rcns2folder}"
                        )

                        # Write the BED file to the output folder as "RandomForestBed.bed"
                        randomforest_bed_output = os.path.join(
                            self.check_and_create_folder(self.output, sampleID),
                            "RandomForestBed.bed",
                        )
                        logger.info(f"Writing BED file to: {randomforest_bed_output}")

                        merged_modkit_df.to_csv(
                            randomforest_bed_output, sep="\t", index=False, header=False
                        )

                        # Initialize batch number if not exists
                        if sampleID not in self.bambatch:
                            self.bambatch[sampleID] = 1
                        logger.info(f"Using batch number: {self.bambatch[sampleID]}")

                        # Run the R script
                        logger.info("Attempting to run R script...")
                        try:
                            await run.cpu_bound(
                                run_rcns2,
                                rcns2folder,
                                self.bambatch[sampleID],
                                randomforest_bed_output,
                                self.threads,
                                self.showerrors,
                            )
                            logger.info("R script completed successfully")
                        except Exception as e:
                            logger.error(
                                f"Error running R script: {str(e)}", exc_info=True
                            )
                            raise

                        votes_file = (
                            f"{rcns2folder}/live_{self.bambatch[sampleID]}_votes.tsv"
                        )
                        if os.path.isfile(votes_file):
                            logger.info(f"Found votes file: {votes_file}")
                            scores = pd.read_table(votes_file, sep="\s+")
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

                            # Save results
                            output_file = os.path.join(
                                self.check_and_create_folder(self.output, sampleID),
                                "random_forest_scores.csv",
                            )
                            logger.info(f"Saving results to: {output_file}")
                            self.rcns2_df_store[sampleID].to_csv(output_file)
                        else:
                            logger.error(f"Votes file not found: {votes_file}")

                            # Counter updated automatically by BaseAnalysis._batch_worker()

                except Exception as e:
                    logger.error(f"Error in process_bam: {str(e)}", exc_info=True)

                    
            await background_tasks.create(
                forest_bam_background_work(sampleID, parquet_path)
            )
            self.running = False
            app.storage.general[self.mainuuid][sampleID][self.name]["counters"]["bams_in_processing"] -= num_bam_files_seen
            app.storage.general[self.mainuuid][sampleID][self.name]["counters"]["bam_processed"] += num_bam_files_seen

    async def stop_analysis(self):
        """Stop the Random Forest analysis."""
        state.set_process_state("Random Forest Analysis", ProcessState.STOPPING)
        state.stop_process("Random Forest Analysis")
        await super().stop_analysis()


class RandomForestVis(BaseVis):
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
        super().__init__(*args, **kwargs)
        # Remove state tracking for Random Forest Analysis
        # state.start_process("Random Forest Analysis", ProcessType.BATCH)
        state.set_process_state("Random Forest Analysis", ProcessState.WAITING_FOR_DATA)
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

        # Set RandomForest-specific confidence thresholds
        # RandomForest confidence is reported as percentage (0-100) but converted to 0-1 in create_summary_card
        kwargs["high_confidence_threshold"] = (
            0.85  # RandomForest-specific high confidence threshold
        )
        kwargs["medium_confidence_threshold"] = (
            0.65  # RandomForest-specific medium confidence threshold
        )

    async def setup_ui(self):
        self.card = ui.card().classes("w-full p-2")
        with self.card:
            with ui.grid(columns=8).classes("w-full h-auto gap-2"):
                with ui.card().classes(
                    f"min-[{self.MENU_BREAKPOINT+1}px]:col-span-3 max-[{self.MENU_BREAKPOINT}px]:col-span-8 p-2"
                ):
                    self.create_rcns2_chart("Random Forest")
                with ui.card().classes(
                    f"min-[{self.MENU_BREAKPOINT+1}px]:col-span-5 max-[{self.MENU_BREAKPOINT}px]:col-span-8 p-2"
                ):
                    self.create_rcns2_time_chart("Random Forest Time Series")
        if self.summary:
            with self.summary:
                ui.label("Forest classification: Unknown")
        # await ui.context.client.connected()
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
                    logger.info(
                        f"DataFrame loaded with shape: {self.rcns2_df_store.shape}"
                    )

                    if not self.rcns2_df_store.empty:
                        columns_greater_than_threshold = (
                            self.rcns2_df_store > self.threshold
                        ).any()
                        columns_not_greater_than_threshold = (
                            ~columns_greater_than_threshold
                        )
                        result = self.rcns2_df_store.columns[
                            columns_not_greater_than_threshold
                        ].tolist()
                        logger.debug(f"Filtered columns: {result}")

                        # Update time series chart
                        filtered_df = self.rcns2_df_store.drop(columns=result)
                        logger.info(
                            f"Updating time chart with filtered data shape: {filtered_df.shape}"
                        )
                        self.update_rcns2_time_chart(filtered_df)

                        # Get last row data
                        lastrow = self.rcns2_df_store.iloc[-1]
                        logger.debug(f"Last row before processing: {lastrow.to_dict()}")

                        n_features = lastrow.get("number_probes", 0)
                        if "number_probes" in lastrow.index:
                            lastrow = lastrow.drop("number_probes")
                            logger.info(
                                f"Dropped number_probes column, features found: {n_features}"
                            )

                        lastrow_plot = lastrow.sort_values(ascending=False).head(10)
                        lastrow_plot_top = lastrow.sort_values(ascending=False).head(1)
                        logger.info(
                            f"Top classification: {lastrow_plot_top.index[0]} ({lastrow_plot_top.values[0]:.1f}%)"
                        )

                        # Update summary card
                        if self.summary:
                            with self.summary:
                                self.summary.clear()
                                classification_text = f"Forest classification: {lastrow_plot_top.index[0]}"
                                logger.debug(
                                    f"Creating summary card with text: {classification_text}"
                                )
                                self.create_summary_card(
                                    classification_text=classification_text,
                                    confidence_value=lastrow_plot_top.values[0] / 100,
                                    features_found=(
                                        int(n_features) if n_features else None
                                    ),
                                )

                        # Update bar plot
                        logger.info("Updating bar plot with top 10 classifications")
                        self.update_rcns2_plot(
                            lastrow_plot.index.to_list(),
                            list(lastrow_plot.values),
                            "All",
                            n_features if n_features else None,
                        )
                        logger.info("Completed show_previous_data successfully")
                except Exception as e:
                    logger.error(f"Error reading scores file: {str(e)}", exc_info=True)
            else:
                logger.debug(f"No scores file found at: {scores_file}")

        except Exception as e:
            logger.error(f"Error in show_previous_data: {str(e)}", exc_info=True)
            raise

    def create_rcns2_chart(self, title):
        """
        Create a bar chart for displaying Random Forest classification results.
        """
        self.echart = self.create_chart(title)
        self.echart.options.update(
            {
                "backgroundColor": "transparent",
                "title": {
                    "text": title,
                    "left": "center",
                    "top": 10,
                    "textStyle": {
                        "fontSize": 16,
                        "fontWeight": "normal",
                        "color": "#000000",
                    },
                },
                "tooltip": {
                    "trigger": "axis",
                    "axisPointer": {"type": "shadow"},
                    "formatter": "{b}: {c}%",
                    "textStyle": {"fontSize": 14},
                },
                "grid": {
                    "left": "5%",
                    "right": "5%",
                    "bottom": "5%",
                    "top": "25%",
                    "containLabel": True,
                },
                "xAxis": {
                    "type": "value",
                    "min": 0,
                    "max": 100,
                    "interval": 20,
                    "axisLabel": {
                        "fontSize": 12,
                        "formatter": "{value}%",
                        "color": "#666666",
                    },
                    "splitLine": {
                        "show": True,
                        "lineStyle": {"color": "#E0E0E0", "type": "dashed"},
                    },
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
                        "color": "#666666",
                    },
                },
                "series": [
                    {
                        "type": "bar",
                        "name": "Confidence",
                        "barMaxWidth": "60%",
                        "itemStyle": {"color": "#007AFF", "borderRadius": [0, 4, 4, 0]},
                        "label": {
                            "show": True,
                            "position": "right",
                            "formatter": "{c}%",
                            "fontSize": 12,
                            "color": "#666666",
                        },
                        "data": [],
                    }
                ],
            }
        )

    def create_rcns2_time_chart(self, title):
        """
        Create a time series chart for Random Forest results.
        """
        self.rcns2_time_chart = self.create_time_chart(title)
        self.rcns2_time_chart.options.update(
            {
                "backgroundColor": "transparent",
                "title": {
                    "text": title,
                    "left": "center",
                    "top": 5,
                    "textStyle": {
                        "fontSize": 16,
                        "fontWeight": "normal",
                        "color": "#000000",
                    },
                    "padding": [0, 0, 20, 0],  # Add padding below title
                },
                "tooltip": {
                    "trigger": "axis",
                    "axisPointer": {"type": "line"},
                    "textStyle": {"fontSize": 14},
                },
                "grid": {
                    "left": "5%",
                    "right": "5%",
                    "bottom": "5%",
                    "top": "25%",
                    "containLabel": True,
                },
                "legend": {
                    "type": "scroll",
                    "orient": "horizontal",
                    "top": 45,  # Increased from 35
                    "width": "90%",
                    "left": "center",
                    "textStyle": {"fontSize": 12, "color": "#666666"},
                    "pageButtonPosition": "end",
                    "pageButtonGap": 5,
                    "pageButtonItemGap": 5,
                    "pageIconColor": "#666666",
                    "pageIconInactiveColor": "#aaa",
                    "pageIconSize": 12,
                    "pageTextStyle": {"color": "#666666"},
                    "itemGap": 25,
                    "itemWidth": 14,
                    "itemHeight": 14,
                    "selectedMode": True,
                },
                "xAxis": {
                    "type": "time",
                    "axisLabel": {
                        "fontSize": 12,
                        "formatter": "{HH}:{mm}",
                        "color": "#666666",
                    },
                    "splitLine": {
                        "show": True,
                        "lineStyle": {"color": "#E0E0E0", "type": "dashed"},
                    },
                },
                "yAxis": {
                    "type": "value",
                    "min": 0,
                    "max": 100,
                    "interval": 20,
                    "axisLabel": {
                        "fontSize": 12,
                        "formatter": "{value}%",
                        "color": "#666666",
                    },
                    "splitLine": {
                        "show": True,
                        "lineStyle": {"color": "#E0E0E0", "type": "dashed"},
                    },
                },
            }
        )

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

            # Sort the data in descending order and take top 10
            sorted_indices = sorted(
                range(len(formatted_values)),
                key=lambda k: formatted_values[k],
                reverse=True,
            )[:10]
            sorted_values = [formatted_values[i] for i in sorted_indices]
            sorted_labels = [x[i] for i in sorted_indices]
            logger.debug(f"Top 10 sorted values: {sorted_values}")
            logger.debug(f"Top 10 sorted labels: {sorted_labels}")

            # Create descriptive title with key information
            title_text = f"Random Forest Analysis Results\n{count} samples processed"
            if n_features:
                title_text += f" • {int(n_features)} features found"
            logger.debug(f"Setting title: {title_text}")

            self.echart.options["title"]["text"] = title_text
            self.echart.options["yAxis"]["data"] = sorted_labels
            self.echart.options["series"][0].update(
                {
                    "data": sorted_values,
                    "itemStyle": {
                        "color": "#007AFF",  # iOS blue
                        "borderRadius": [0, 4, 4, 0],
                    },
                }
            )
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

            # Ensure DataFrame is sorted by index (timestamp)
            datadf = datadf.sort_index()

            # Get the top 10 diagnoses based on the latest data point
            latest_data = datadf.iloc[-1]
            if "number_probes" in latest_data:
                latest_data = latest_data.drop("number_probes")
            top_10_diagnoses = latest_data.nlargest(10).index.tolist()

            # Filter dataframe to only include top 10 diagnoses
            filtered_df = datadf[top_10_diagnoses]

            for idx, (series, data) in enumerate(filtered_df.to_dict().items()):
                logger.debug(f"Processing series: {series}")
                # Values are already percentages, just format them and ensure sorted by timestamp
                data_list = [
                    [key, float(f"{value:.1f}")]
                    for key, value in sorted(data.items())  # Sort by timestamp
                ]
                logger.debug(f"First few data points for {series}: {data_list[:3]}")

                self.rcns2_time_chart.options["series"].append(
                    {
                        "name": series,
                        "type": "line",
                        "smooth": True,
                        "animation": False,
                        "symbolSize": 6,
                        "emphasis": {
                            "focus": "series",
                            "itemStyle": {"borderWidth": 2},
                        },
                        "endLabel": {
                            "show": True,
                            "formatter": "{a}: {c}%",
                            "distance": 10,
                            "fontSize": 12,
                        },
                        "lineStyle": {"width": 2, "color": colors[idx % len(colors)]},
                        "itemStyle": {"color": colors[idx % len(colors)]},
                        "data": data_list,
                    }
                )

            # Update chart title with summary
            max_confidence = (
                latest_data.max()
            )  # Get maximum confidence (already percentage)
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
