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
import shutil

from robin import submodules

from robin.utilities.merge_bedmethyl import (
    merge_bedmethyl,
    save_bedmethyl,
    collapse_bedmethyl,
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
        
        """
        # Check if input file exists
        if not os.path.exists(bed):
            logger.error(f"Input BED file does not exist: {bed}")
            raise FileNotFoundError(f"BED file not found: {bed}")
            
        # Read the bed file
        logger.debug("Reading BED file...")
        bed_df = pd.read_csv(bed, sep="\t")
        logger.info(f"Read BED file with shape: {bed_df.shape}")
        logger.debug(f"BED file columns: {bed_df.columns.tolist()}")
        
        # Ensure numeric types for position columns
        logger.debug("Converting position columns to numeric...")
        bed_df["chromStart"] = pd.to_numeric(bed_df["chromStart"], errors="coerce")
        bed_df["chromEnd"] = pd.to_numeric(bed_df["chromEnd"], errors="coerce")
        
        # Drop any rows with NaN values after conversion
        bed_df = bed_df.dropna(subset=["chromStart", "chromEnd"])
        
        # Rename columns to match expected format
        logger.debug("Renaming columns...")
        bed_df = bed_df.rename(columns={
            "chromStart": "start",  # Changed from start_pos to start
            "chromEnd": "end",      # Changed from end_pos to end
            "percent_modified": "methylation_call"
        })
        
        
        # Convert methylation calls to binary format (0/1)
        logger.debug("Converting methylation calls to binary format...")
        bed_df["methylation_call"] = (bed_df["methylation_call"] >= 60).astype(int)
        
        # Ensure all required columns are present
        required_columns = ["chrom", "start", "end", "methylation_call"]
        missing_columns = [col for col in required_columns if col not in bed_df.columns]
        if missing_columns:
            logger.error(f"Missing required columns: {missing_columns}")
            raise ValueError(f"Missing required columns: {missing_columns}")
        
        # Save in format expected by R script
        converted_bed = os.path.join(os.path.dirname(bed), "converted_for_r.bed")
        logger.info(f"Saving converted BED file to: {converted_bed}")
        # Save only required columns in correct order
        bed_df.to_csv(converted_bed, sep="\t", index=False, header=False)
        """
        
        # Check if R script exists
        r_script_path = f"{HVPATH}/bin/methylation_classification_nanodx_v0.1.R"
        if not os.path.exists(r_script_path):
            logger.error(f"R script not found at: {r_script_path}")
            raise FileNotFoundError(f"R script not found: {r_script_path}")
            
        # Check if other required files exist
        required_files = {
            "probes": f"{HVPATH}/bin/top_probes_hm450.Rdata",
            "training_data": f"{HVPATH}/bin/capper_top_100k_betas_binarised.Rdata",
            "array_file": f"{HVPATH}/bin/HM450.hg38.manifest.gencode.v22.Rdata"
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
                command.split(),
                capture_output=True,
                text=True,
                check=True
            )
            logger.info("R script executed successfully")
            logger.debug(f"R script stdout: {result.stdout}")
            if result.stderr:
                logger.warning(f"R script stderr: {result.stderr}")
        except subprocess.CalledProcessError as e:
            logger.error(f"R script failed with return code {e.returncode}")
            logger.error(f"stdout: {e.stdout}")
            logger.error(f"stderr: {e.stderr}")
            raise
        
        # Check if output file was created
        expected_output = f"{rcns2folder}/live_{batch}_votes.tsv"
        if os.path.exists(expected_output):
            logger.info(f"Output file created successfully: {expected_output}")
        else:
            logger.error(f"Expected output file not found: {expected_output}")
            raise FileNotFoundError(f"R script did not create expected output file: {expected_output}")
        
    except Exception as e:
        logger.error(f"Error in run_rcns2: {str(e)}", exc_info=True)
        raise


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
        logger.error(e)
        # self.log(e)
        pass


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
        Process BAM files and perform Random Forest analysis.

        Parameters
        ----------
        bamfile : List[Tuple[str, float]]
            List of tuples containing BAM file paths and their timestamps
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
            logger.info(f"Processing parquet file: {parquet_path}")
            
            try:
                if self.check_file_time(parquet_path):
                    logger.debug("Parquet file exists and is ready for processing")
                    tomerge_length_file = os.path.join(self.check_and_create_folder(self.output, sampleID), "tomerge_length.txt")
                    try:
                        with open(tomerge_length_file, "r") as f:
                            tomerge_length = int(f.readline().strip().split(": ")[1])
                        logger.info(f"Number of files to merge: {tomerge_length}")
                    except FileNotFoundError:
                        logger.warning("tomerge_length.txt not found yet, waiting for file creation")
                        return
                    except Exception as e:
                        logger.error(f"Error reading tomerge_length.txt: {str(e)}")
                        return
                    
                    logger.debug("Loading modkit data from parquet file...")
                    merged_modkit_df = await run.cpu_bound(load_modkit_data, parquet_path)
                    print (merged_modkit_df.columns)
                    merged_modkit_df.rename(columns={"chromStart":"start_pos", "chromEnd":"end_pos", "mod_code":"mod", "thickStart":"start_pos2", "thickEnd":"end_pos2", "color":"colour"}, inplace=True)
                    merged_modkit_df.rename(columns={'n_canonical':'Ncanon', 'n_delete':'Ndel', 'n_diff':'Ndiff', 'n_fail':'Nfail', 'n_mod':'Nmod', 'n_nocall':'Nnocall', 'n_othermod':'Nother', 'valid_cov':'Nvalid', 'percent_modified':'score'}, inplace=True)
                   
                    forest_dx = await run.cpu_bound(collapse_bedmethyl, merged_modkit_df)
                    merged_modkit_df = forest_dx
                    if merged_modkit_df is None:
                        logger.error("Failed to load modkit data from parquet file")
                        return
                    logger.info(f"Loaded modkit data with shape: {merged_modkit_df.shape}")
                    
                    # Create a temporary directory for R script output
                    rcns2folder = tempfile.mkdtemp(dir=self.check_and_create_folder(self.output, sampleID))
                    logger.info(f"Created temporary directory for R output: {rcns2folder}")
                    
                    # Write the BED file to the output folder as "RandomForestBed.bed"
                    randomforest_bed_output = os.path.join(self.check_and_create_folder(self.output, sampleID), "RandomForestBed.bed")
                    logger.info(f"Writing BED file to: {randomforest_bed_output}")
                    
                    
                    
                    merged_modkit_df.to_csv(randomforest_bed_output, sep="\t", index=False, header=False)

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
                        logger.error(f"Error running R script: {str(e)}", exc_info=True)
                        raise

                    votes_file = f"{rcns2folder}/live_{self.bambatch[sampleID]}_votes.tsv"
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

                    app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                        "bam_processed"
                    ] = tomerge_length

            except Exception as e:
                logger.error(f"Error in process_bam: {str(e)}", exc_info=True)
                

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
                "top": 10,
                "textStyle": {
                    "fontSize": 16,
                    "fontWeight": "normal",
                    "color": "#000000"
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
                "top": 5,
                "textStyle": {
                    "fontSize": 16,
                    "fontWeight": "normal",
                    "color": "#000000"
                },
                "padding": [0, 0, 20, 0]  # Add padding below title
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
                "top": 45,  # Increased from 35
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
            sorted_indices = sorted(range(len(formatted_values)), key=lambda k: formatted_values[k], reverse=True)[:10]
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
            
            # Get the top 10 diagnoses based on the latest data point
            latest_data = datadf.iloc[-1]
            if "number_probes" in latest_data:
                latest_data = latest_data.drop("number_probes")
            top_10_diagnoses = latest_data.nlargest(10).index.tolist()
            
            # Filter dataframe to only include top 10 diagnoses
            filtered_df = datadf[top_10_diagnoses]
            
            for idx, (series, data) in enumerate(filtered_df.to_dict().items()):
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
