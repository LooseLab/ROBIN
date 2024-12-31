"""
A module for analyzing and visualizing Sturgeon methylation data in real-time.

This module provides functionality for processing BAM files to analyze methylation patterns
using the Sturgeon tool suite. It includes capabilities for:

- Real-time BAM file processing
- Methylation data extraction using modkit
- Sturgeon prediction and visualization
- Time series analysis of methylation patterns

The module integrates with the NiceGUI framework for interactive visualization
and uses temporary files for efficient data processing.

Dependencies
-----------
- pysam
- pandas
- nicegui
- sturgeon
- modkit

Notes
-----
The module requires proper configuration of input/output directories and
assumes the presence of necessary model files for Sturgeon predictions.
"""

from robin.subpages.base_analysis import BaseAnalysis
import os
import sys
import tempfile
import time
import shutil
import pandas as pd
from nicegui import ui, app, run
from robin import theme
import pysam
from robin import models
from sturgeon.callmapping import (
    merge_probes_methyl_calls,
    probes_methyl_calls_to_bed,
)
import click
from pathlib import Path
from typing import List, Tuple
import logging

# Use the main logger configured in the main application
logger = logging.getLogger(__name__)


def run_probes_methyl_calls(merged_output_file, bed_output_file):
    """
    Convert merged methylation calls to BED format.

    Parameters
    ----------
    merged_output_file : str
        Path to the input file containing merged methylation calls
    bed_output_file : str
        Path where the output BED file should be written

    Notes
    -----
    This function is a wrapper around sturgeon's probes_methyl_calls_to_bed function.
    """
    probes_methyl_calls_to_bed(merged_output_file, bed_output_file)


def run_sturgeon_merge_probes(calls_per_probe_file, merged_output_file):
    """
    Merge multiple probe methylation call files.

    Parameters
    ----------
    calls_per_probe_file : str
        Path to the input file containing methylation calls per probe
    merged_output_file : str
        Path where the merged output should be written

    Notes
    -----
    This function merges multiple probe methylation call files into a single output file
    using sturgeon's merge_probes_methyl_calls function.
    """
    merge_probes_methyl_calls(
        [calls_per_probe_file, merged_output_file],
        merged_output_file,
    )


def pysam_cat(tempbam, tomerge):
    """
    Concatenate multiple BAM files using pysam.

    Parameters
    ----------
    tempbam : str
        Path where the concatenated BAM file should be written
    tomerge : list
        List of BAM file paths to merge

    Notes
    -----
    This function uses pysam's cat functionality to merge multiple BAM files
    into a single output file.
    """
    pysam.cat("-o", tempbam, *tomerge)


def run_modkit(file, temp, threads):
    """
    Run modkit to extract methylation data from a BAM file.

    Parameters
    ----------
    file : str
        Path to input BAM file
    temp : str
        Path for temporary output
    threads : int
        Number of threads to use for processing

    Notes
    -----
    This function adapts its command based on the installed modkit version.
    For versions >= 0.4, it uses 'extract full', otherwise just 'extract'.
    """
    try:
        # Get modkit version
        import subprocess
        version_output = subprocess.check_output(['modkit', '--version'], text=True).strip()
        version = version_output.split()[-1]  # Gets '0.4.1' from 'mod_kit 0.4.1'
        
        # Parse version number
        major, minor, *_ = version.split('.')
        version_num = float(f"{major}.{minor}")
        
        # Choose appropriate command based on version
        extract_cmd = "extract full" if version_num >= 0.4 else "extract"
        
        os.system(
            f"modkit {extract_cmd} --ignore h -t {threads} {file} {temp} "
            f"--force --suppress-progress >/dev/null 2>&1"
        )
    except Exception as e:
        print(e)
        pass


def run_sturgeon_predict(bedDir, dataDir, modelfile):
    """
    Run Sturgeon prediction on methylation data.

    Parameters
    ----------
    bedDir : str
        Directory containing BED format methylation data
    dataDir : str
        Directory where prediction output will be saved
    modelfile : str
        Path to the Sturgeon model file to use for predictions

    Notes
    -----
    This function executes the Sturgeon predict command with the specified
    model file and suppresses terminal output.
    """
    os.system(
        f"sturgeon predict -i {bedDir} -o {dataDir} "
        f"--model-files {modelfile} >/dev/null 2>&1"
    )


def run_sturgeon_inputtobed(temp, temp2):
    """
    Convert Sturgeon input format to BED format.

    Parameters
    ----------
    temp : str
        Path to input file in modkit format
    temp2 : str
        Directory where BED format output will be saved

    Notes
    -----
    This function converts modkit format methylation data to BED format
    using hg38 as the reference genome. Errors are caught and printed.
    """
    try:
        os.system(
            f"sturgeon inputtobed -i {temp} -o {temp2} -s modkit "
            f"--reference-genome hg38 >/dev/null 2>&1"
        )
    except Exception as e:
        print(e)
        pass


class Sturgeon_object(BaseAnalysis):
    """
    A class for processing and visualizing Sturgeon methylation analysis results.

    This class extends BaseAnalysis to provide specialized functionality for
    Sturgeon methylation analysis. It handles real-time processing of BAM files,
    methylation data extraction, and visualization of results.

    Parameters
    ----------
    *args
        Variable length argument list passed to BaseAnalysis
    **kwargs
        Arbitrary keyword arguments passed to BaseAnalysis

    Attributes
    ----------
    sturgeon_df_store : dict
        Store for Sturgeon dataframes
    threshold : float
        Confidence threshold for predictions (default: 0.05)
    first_run : dict
        Tracks first run status for each sample
    modelfile : str
        Path to the Sturgeon model file
    dataDir : dict
        Dictionary of temporary directories for data
    bedDir : dict
        Dictionary of temporary directories for BED files
    st_num_probes : dict
        Dictionary tracking number of probes per sample
    """
    def __init__(self, *args, **kwargs):
        self.sturgeon_df_store = {}
        self.threshold = 0.05
        self.first_run = {}
        self.modelfile = os.path.join(
            os.path.dirname(os.path.abspath(models.__file__)), "general.zip"
        )
        self.dataDir = {}
        self.bedDir = {}
        self.st_num_probes = {}
        super().__init__(*args, **kwargs)

    def setup_ui(self):
        """
        Set up the user interface components for Sturgeon visualization.

        This method creates and arranges the UI cards and grids for displaying
        Sturgeon analysis results. It sets up:
        - Main card with dark theme
        - Grid layout for charts
        - Sturgeon classification chart
        - Time series chart
        - Initial classification label

        Notes
        -----
        The layout is responsive and adjusts based on the MENU_BREAKPOINT value.
        """
        self.card = ui.card().classes('dark:bg-black').style("width: 100%")
        with self.card:
            with ui.grid(columns=8).classes("w-full h-auto gap-4"):
                with ui.card().classes(
                    f"min-[{self.MENU_BREAKPOINT+1}px]:col-span-3 max-[{self.MENU_BREAKPOINT}px]:col-span-8 dark:bg-black shadow-lg rounded-lg"
                ).style('background-color: #FFFFFF; padding: 16px;'):
                    self.create_sturgeon_chart("Sturgeon")
                with ui.card().classes(
                    f"min-[{self.MENU_BREAKPOINT+1}px]:col-span-5 max-[{self.MENU_BREAKPOINT}px]:col-span-8 dark:bg-black shadow-lg rounded-lg"
                ).style('background-color: #FFFFFF; padding: 16px;'):
                    self.create_sturgeon_time_chart("Sturgeon Time Series")
        if self.summary:
            with self.summary:
                ui.label(f"Sturgeon classification: Unknown")
        if self.browse:
            self.show_previous_data()
        else:
            ui.timer(5, lambda: self.show_previous_data())

    def show_previous_data(self):
        """
        Load and display previously generated Sturgeon analysis results.
        """
        if not self.browse:
            for item in app.storage.general[self.mainuuid]:
                if item == "sample_ids":
                    for sample in app.storage.general[self.mainuuid][item]:
                        self.sampleID = sample
            output = self.output
        if self.browse:
            output = self.check_and_create_folder(self.output, self.sampleID)

        if self.check_file_time(os.path.join(output, "sturgeon_scores.csv")):
            self.sturgeon_df_store = pd.read_csv(
                os.path.join(os.path.join(output, "sturgeon_scores.csv")),
                index_col=0,
            )
            columns_greater_than_threshold = (
                self.sturgeon_df_store > self.threshold
            ).any()
            columns_not_greater_than_threshold = ~columns_greater_than_threshold
            result = self.sturgeon_df_store.columns[
                columns_not_greater_than_threshold
            ].tolist()
            self.update_sturgeon_time_chart(self.sturgeon_df_store.drop(columns=result))
            lastrow = self.sturgeon_df_store.iloc[-1].drop("number_probes")
            lastrow_plot = lastrow.sort_values(ascending=False).head(10)
            lastrow_plot_top = lastrow.sort_values(ascending=False).head(1)
            
            # Update summary with new card
            if self.summary:
                with self.summary:
                    self.summary.clear()
                    classification_text = f"Sturgeon classification: {lastrow_plot_top.index[0]}"
                    self.create_summary_card(
                        classification_text=classification_text,
                        confidence_value=lastrow_plot_top.values[0],
                        features_found=int(self.sturgeon_df_store.iloc[-1]["number_probes"])
                    )
            
            self.update_sturgeon_plot(
                lastrow_plot.index.to_list(),
                list(lastrow_plot.values),
                "All",
                self.sturgeon_df_store.iloc[-1]["number_probes"],
            )

    async def process_bam(self, bamfile: List[Tuple[str, float]]) -> None:
        """
        Process BAM files and perform Sturgeon analysis.

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
                dir=self.check_and_create_folder(self.output, sampleID)
            )

            await run.cpu_bound(pysam_cat, tempbam.name, tomerge)

            file = tempbam.name
            temp = tempfile.NamedTemporaryFile(
                dir=self.check_and_create_folder(self.output, sampleID)
            )
            with tempfile.TemporaryDirectory(
                dir=self.check_and_create_folder(self.output, sampleID)
            ) as temp2:
                await run.cpu_bound(run_modkit, file, temp.name, self.threads)

                await run.cpu_bound(run_sturgeon_inputtobed, temp.name, temp2)

                calls_per_probe_file = os.path.join(
                    temp2, "merged_probes_methyl_calls.txt"
                )
                merged_output_file = os.path.join(
                    self.dataDir[sampleID].name,
                    "_merged_probes_methyl_calls.txt",
                )

                if sampleID not in self.first_run.keys():
                    self.first_run[sampleID] = True
                    shutil.copyfile(calls_per_probe_file, merged_output_file)
                else:
                    await run.cpu_bound(
                        run_sturgeon_merge_probes,
                        calls_per_probe_file,
                        merged_output_file,
                    )

                bed_output_file = os.path.join(
                    self.bedDir[sampleID].name, "final_merged_probes_methyl_calls.bed"
                )

                await run.cpu_bound(
                    run_probes_methyl_calls, merged_output_file, bed_output_file
                )

                await run.cpu_bound(
                    run_sturgeon_predict,
                    self.bedDir[sampleID].name,
                    self.dataDir[sampleID].name,
                    self.modelfile,
                )
                if os.path.exists(
                    os.path.join(
                        self.dataDir[sampleID].name,
                        "final_merged_probes_methyl_calls_general.csv",
                    )
                ):
                    mydf = pd.read_csv(
                        os.path.join(
                            self.dataDir[sampleID].name,
                            "final_merged_probes_methyl_calls_general.csv",
                        )
                    )
                else:
                    self.running = False
                    return

                self.st_num_probes[sampleID] = mydf.iloc[-1]["number_probes"]
                # lastrow = mydf.iloc[-1].drop("number_probes")
                mydf_to_save = mydf
                mydf_to_save["timestamp"] = currenttime

                if sampleID not in self.sturgeon_df_store:
                    self.sturgeon_df_store[sampleID] = pd.DataFrame()

                # Exclude empty or all-NA entries before concatenation
                if not mydf_to_save.dropna(how="all").empty:
                    self.sturgeon_df_store[sampleID] = pd.concat(
                        [
                            self.sturgeon_df_store[sampleID],
                            mydf_to_save.set_index("timestamp"),
                        ]
                    )
                    self.sturgeon_df_store[sampleID].to_csv(
                        os.path.join(
                            self.check_and_create_folder(self.output, sampleID),
                            "sturgeon_scores.csv",
                        )
                    )
            app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                "bam_processed"
            ] += len(tomerge)

            app.storage.general[self.mainuuid][sampleID][self.name]["counters"][
                "bams_in_processing"
            ] -= len(tomerge)

        self.running = False

    def create_sturgeon_chart(self, title):
        """
        Create a bar chart for displaying Sturgeon classification results.
        """
        self.echart2 = self.create_chart(title)
        self.echart2.options.update({
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
                "axisPointer": {"type": "shadow"},
                "formatter": "{b}: {c}%",
                "textStyle": {"fontSize": 14}
            },
            "grid": {
                "left": "15%",
                "right": "10%",
                "bottom": "10%",
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
            }],
            "aria": {
                "enabled": True,
                "decal": {
                    "show": True
                }
            }
        })

    def create_sturgeon_time_chart(self, title):
        """
        Create a time series chart for Sturgeon results.
        """
        self.sturgeon_time_chart = self.create_time_chart(title)
        self.sturgeon_time_chart.options.update({
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

    def update_sturgeon_plot(self, x, y, count, st_num_probes):
        """
        Update the Sturgeon visualization plot with new data.

        Parameters
        ----------
        x : List[str]
            List of tumor types
        y : List[float]
            Confidence scores for each tumor type
        count : str
            Number of BAM files processed
        st_num_probes : int
            Number of probes found
        """
        # Convert values to percentages and format
        formatted_values = [float(f"{val * 100:.1f}") for val in y]
        
        # Sort the data in descending order
        sorted_indices = sorted(range(len(formatted_values)), key=lambda k: formatted_values[k], reverse=True)
        sorted_values = [formatted_values[i] for i in sorted_indices]
        sorted_labels = [x[i] for i in sorted_indices]
        
        # Create descriptive title with key information
        title_text = (
            f"Sturgeon Analysis Results\n"
            f"{count} samples processed • {int(st_num_probes)} probes found"
        )
        
        self.echart2.options["title"]["text"] = title_text
        self.echart2.options["yAxis"]["data"] = sorted_labels
        self.echart2.options["series"][0].update({
            "data": sorted_values,
            "itemStyle": {
                "color": "#007AFF",
                "borderRadius": [0, 4, 4, 0]
            }
        })
        self.echart2.update()

    def update_sturgeon_time_chart(self, datadf):
        """
        Update the time series chart with new data.

        Parameters
        ----------
        datadf : pd.DataFrame
            DataFrame containing time series data for visualization
        """
        self.sturgeon_time_chart.options["series"] = []
        
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
                # Convert values to percentages
                data_list = [[key, float(f"{value * 100:.1f}")] for key, value in data.items()]
                self.sturgeon_time_chart.options["series"].append({
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
        latest_data = datadf.iloc[-1].drop("number_probes")  # Remove number_probes from latest data
        max_confidence = latest_data.max() * 100  # Convert to percentage
        max_type = latest_data.idxmax()
        self.sturgeon_time_chart.options["title"]["text"] = (
            f"Classification Confidence Over Time\n"
            f"Current highest confidence: {max_type} ({max_confidence:.1f}%)"
        )
        
        self.sturgeon_time_chart.update()


def test_me(
    port: int,
    threads: int,
    watchfolder: str,
    output: str,
    reload: bool = False,
    browse: bool = False,
):
    """
    Test function to run the Sturgeon analysis interface.

    Parameters
    ----------
    port : int
        Port number for the web interface
    threads : int
        Number of processing threads
    watchfolder : str
        Directory to watch for new BAM files
    output : str
        Directory for output files
    reload : bool, optional
        Whether to enable auto-reload (default: False)
    browse : bool, optional
        Whether to run in browse mode (default: False)

    Notes
    -----
    This function sets up the web interface and starts processing BAM files
    either from a watch folder or in browse mode.
    """
    my_connection = None
    with theme.frame("Sturgeon Rapid CNS Diagnostic.", my_connection):
        TestObject = Sturgeon_object(threads, output, progress=True, batch=True)
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
        print("Browse mode not implemented.")
        TestObject.progress_trackers.visible = False
        # TestObject.show_previous_data(output)
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
def mainrun(port, threads, watchfolder, output, browse):
    """
    Main entry point for the Sturgeon analysis application.

    This function handles command line arguments and launches the application
    either in watch mode or browse mode.

    Parameters
    ----------
    port : int
        Port number for the web interface
    threads : int
        Number of processing threads
    watchfolder : Path
        Directory to watch for new BAM files
    output : Path
        Directory for output files
    browse : bool
        Whether to run in browse mode

    Notes
    -----
    In watch mode, both watchfolder and output are required.
    In browse mode, only output is required.
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
    mainrun()
