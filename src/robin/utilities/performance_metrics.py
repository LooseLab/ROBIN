"""
Performance Metrics Module
=========================

This module provides a reusable performance metrics tracking system for ROBIN pages.
It handles data collection, persistence, and UI components for performance monitoring.

Features
--------
* Real-time performance tracking
* Data persistence in both numpy and CSV formats
* Interactive UI components
* Statistical analysis
* Time series visualization
"""

import os
import time
import numpy as np
import pandas as pd
import logging
from nicegui import ui
from pathlib import Path
from typing import Dict, List, Optional, Union

logger = logging.getLogger(__name__)


class PerformanceMetrics:
    """
    A reusable class for tracking and displaying performance metrics.

    This class provides functionality for:
    * Collecting performance data
    * Saving and loading metrics
    * Creating and updating UI components
    * Statistical analysis
    * Time series visualization
    """

    def __init__(
        self, output_dir: Union[str, Path], metrics: Optional[Dict[str, List]] = None
    ):
        """
        Initialize the performance metrics tracker.

        Parameters
        ----------
        output_dir : Union[str, Path]
            Directory where performance data will be stored
        metrics : Optional[Dict[str, List]]
            Initial metrics dictionary (optional)
        """
        self.output_dir = Path(output_dir)
        self.perf_dir = self.output_dir / "performance"
        self.perf_dir.mkdir(exist_ok=True)

        # Initialize metrics dictionary with empty lists
        self.metrics = {
            "bam_merge_times": [],
            "df_merge_times": [],
            "total_processing_times": [],
            "timestamps": [],
            "bam_file_sizes": [],
            "bam_file_names": [],
        }

        # Update with any provided metrics
        if metrics:
            for key, value in metrics.items():
                if key in self.metrics:
                    self.metrics[key] = value

        # UI components
        self.performance_container = None
        self.performance_chart = None
        self.stats_cards = {}
        self.ui_initialized = False

        logger.info("Performance metrics initialized")

    def start_tracking(self) -> float:
        """
        Start tracking performance for a new operation.

        Returns
        -------
        float
            Start time in seconds
        """
        return time.time()

    def record_metrics(
        self,
        processing_time: float,
        file_path: Optional[str] = None,
        additional_metrics: Optional[Dict] = None,
    ) -> None:
        """
        Record performance metrics for an operation.

        Parameters
        ----------
        processing_time : float
            Time taken for the operation in seconds
        file_path : Optional[str]
            Path to the processed file (optional)
        additional_metrics : Optional[Dict]
            Additional metrics to record (optional)
        """
        current_time = time.time() * 1000  # Convert to milliseconds

        # Get current length of arrays
        current_length = len(self.metrics["timestamps"])

        # Append new values
        self.metrics["total_processing_times"].append(processing_time)
        self.metrics["timestamps"].append(current_time)

        if file_path:
            try:
                file_size = os.path.getsize(file_path)
                self.metrics["bam_file_sizes"].append(file_size)
                self.metrics["bam_file_names"].append(os.path.basename(file_path))
            except Exception as e:
                logger.warning(f"Could not get file size for {file_path}: {e}")
                self.metrics["bam_file_sizes"].append(None)
                self.metrics["bam_file_names"].append(None)
        else:
            self.metrics["bam_file_sizes"].append(None)
            self.metrics["bam_file_names"].append(None)

        # Initialize additional metrics if provided
        if additional_metrics:
            for key, value in additional_metrics.items():
                if key not in self.metrics:
                    self.metrics[key] = [None] * current_length
                self.metrics[key].append(value)

        # Ensure all arrays are the same length
        new_length = len(self.metrics["timestamps"])
        for key in self.metrics:
            if len(self.metrics[key]) < new_length:
                self.metrics[key].extend([None] * (new_length - len(self.metrics[key])))

        logger.info(f"Recorded metrics for operation at {current_time}")
        self.save_metrics()

    def save_metrics(self) -> None:
        """Save performance metrics to disk."""
        try:
            # Verify all arrays are the same length
            lengths = {
                key: len(values)
                for key, values in self.metrics.items()
                if isinstance(values, list)
            }
            if not all(length == lengths["timestamps"] for length in lengths.values()):
                logger.error("Arrays are not the same length before saving")
                return

            # Save numpy arrays
            for key, values in self.metrics.items():
                if isinstance(values, list) and values:
                    np.save(self.perf_dir / f"{key}.npy", np.array(values))

            # Save CSV
            perf_df = pd.DataFrame(
                {
                    "timestamp": self.metrics["timestamps"],
                    "file_name": self.metrics["bam_file_names"],
                    "file_size_mb": [
                        size / (1024 * 1024) if size is not None else None
                        for size in self.metrics["bam_file_sizes"]
                    ],
                    "bam_merge_time_s": self.metrics["bam_merge_times"],
                    "df_merge_time_s": self.metrics["df_merge_times"],
                    "total_processing_time_s": self.metrics["total_processing_times"],
                }
            )
            perf_df.to_csv(self.perf_dir / "performance_metrics.csv", index=False)

            logger.info("Performance metrics saved successfully")

        except Exception as e:
            logger.error(f"Error saving performance metrics: {e}")

    def load_metrics(self) -> bool:
        """
        Load performance metrics from disk.

        Returns
        -------
        bool
            True if metrics were loaded successfully, False otherwise
        """
        try:
            if not self.perf_dir.exists():
                return False

            # Load numpy arrays with allow_pickle=True
            for key in self.metrics.keys():
                file_path = self.perf_dir / f"{key}.npy"
                if file_path.exists():
                    self.metrics[key] = np.load(file_path, allow_pickle=True).tolist()

            # Load additional data from CSV
            csv_path = self.perf_dir / "performance_metrics.csv"
            if csv_path.exists():
                perf_df = pd.read_csv(csv_path)
                for col in perf_df.columns:
                    if col not in self.metrics:
                        self.metrics[col] = perf_df[col].tolist()

            logger.info("Performance metrics loaded successfully")
            return True

        except Exception as e:
            logger.error(f"Error loading performance metrics: {e}")
            return False

    def create_ui(self) -> None:
        """Create the performance metrics UI components."""
        with ui.card().classes("w-full p-4 mb-4") as self.performance_container:
            # Header with title and legend
            with ui.row().classes("w-full items-center justify-between mb-4"):
                ui.label("Performance Metrics").classes("text-lg font-medium")
                with ui.button(icon="help_outline", color="gray").props("flat round"):
                    with ui.menu():
                        with ui.card().classes("p-4 text-sm"):
                            ui.label("Performance Metrics Legend").classes(
                                "font-medium mb-2"
                            )
                            ui.label(
                                "• BAM Merge Time: Time taken to merge BAM files and extract target regions"
                            ).classes("mb-1")
                            ui.label(
                                "• DataFrame Merge Time: Time taken to merge and process coverage data"
                            ).classes("mb-1")
                            ui.label(
                                "• Total Processing Time: Total time for complete BAM processing pipeline"
                            ).classes("mb-1")
                            ui.separator().classes("my-2")
                            ui.label("Statistics").classes("font-medium mb-2")
                            ui.label("• Min: Minimum recorded time").classes("mb-1")
                            ui.label("• Max: Maximum recorded time").classes("mb-1")
                            ui.label("• Mean: Average processing time").classes("mb-1")
                            ui.label("• Median: Middle value of all times").classes(
                                "mb-1"
                            )

            # Stats cards
            with ui.grid(columns=3).classes("w-full gap-4 mb-4"):
                self._create_stats_card(
                    "BAM Merge Times",
                    "bam_merge_times",
                    "s",
                    "Time taken to merge BAM files and extract target regions",
                )
                self._create_stats_card(
                    "DataFrame Merge Times",
                    "df_merge_times",
                    "s",
                    "Time taken to merge and process coverage data",
                )
                self._create_stats_card(
                    "Total Processing Times",
                    "total_processing_times",
                    "s",
                    "Total time for complete BAM processing pipeline",
                )

            # Time series chart
            with ui.card().classes("w-full p-4 border rounded-lg shadow-sm"):
                self._create_time_series_chart()

        self.ui_initialized = True
        logger.info("Performance metrics UI created")

    def _create_stats_card(
        self, title: str, metric_key: str, unit: str, description: str
    ) -> None:
        """Create a stats card for a specific metric."""
        with ui.card().classes("p-4 border rounded-lg shadow-sm") as card:
            # Store the title and description in the card's metadata
            card.metadata = {"title": title, "description": description, "unit": unit}
            with ui.row().classes("items-center gap-2 mb-2"):
                ui.label(title).classes("font-medium")
                with ui.tooltip(description):
                    ui.icon("info_outline").classes("text-gray-500")

            # Create a container for stats with pre-defined labels
            with ui.column().classes("w-full") as stats_container:
                # Create a grid for the statistics with pre-defined structure
                with ui.grid(columns=2).classes("w-full gap-1"):
                    # Headers (static)
                    ui.label("Statistic").classes("text-sm font-medium text-gray-600")
                    ui.label(f"Time ({unit})").classes(
                        "text-sm font-medium text-gray-600 text-right"
                    )

                    # Create value labels that will be updated
                    ui.label("Min").classes("text-sm")
                    min_value = ui.label("--").classes("text-sm text-right")
                    ui.label("Max").classes("text-sm")
                    max_value = ui.label("--").classes("text-sm text-right")
                    ui.label("Mean").classes("text-sm")
                    mean_value = ui.label("--").classes("text-sm text-right")
                    ui.label("Median").classes("text-sm")
                    median_value = ui.label("--").classes("text-sm text-right")

            # Store the value labels for updates
            self.stats_cards[metric_key] = {
                "card": card,
                "stats_container": stats_container,
                "values": {
                    "min": min_value,
                    "max": max_value,
                    "mean": mean_value,
                    "median": median_value,
                },
            }

    def _create_time_series_chart(self) -> None:
        """Create the time series chart."""
        self.performance_chart = ui.echart(
            {
                "backgroundColor": "transparent",
                "textStyle": {
                    "fontFamily": "SF Pro Text, -apple-system, BlinkMacSystemFont, Helvetica, Arial, sans-serif",
                    "fontSize": 12,
                },
                "grid": {
                    "left": "5%",
                    "right": "5%",
                    "bottom": "5%",
                    "top": "15%",
                    "containLabel": True,
                },
                "title": {
                    "text": "Processing Times Over Time",
                    "left": "center",
                    "top": 5,
                    "textStyle": {
                        "fontSize": 16,
                        "fontWeight": "500",
                        "color": "#1D1D1F",
                    },
                },
                "legend": {
                    "data": [
                        "BAM Merge Time",
                        "DataFrame Merge Time",
                        "Total Processing Time",
                    ],
                    "top": 35,
                    "textStyle": {"color": "#86868B"},
                    ":formatter": """
                    function(name) {
                        const descriptions = {
                            'BAM Merge Time': 'Time taken to merge BAM files and extract target regions',
                            'DataFrame Merge Time': 'Time taken to merge and process coverage data',
                            'Total Processing Time': 'Total time for complete BAM processing pipeline'
                        };
                        return descriptions[name] || name;
                    }
                """,
                },
                "tooltip": {
                    "trigger": "axis",
                    "backgroundColor": "rgba(255, 255, 255, 0.9)",
                    "borderColor": "#E5E5EA",
                    "textStyle": {"color": "#1D1D1F"},
                    ":formatter": """
                    function(params) {
                        if (!params || !params.length || !params[0].value) return '';
                        
                        var date = new Date(params[0].value[0]);
                        var timeStr = date.toLocaleTimeString();
                        var result = timeStr + '<br/>';
                        
                        params.forEach(function(param) {
                            if (param && param.value) {
                                var color = param.color;
                                var marker = '<span style="display:inline-block;margin-right:4px;border-radius:10px;width:10px;height:10px;background-color:' + color + '"></span>';
                                var value = param.value[1];
                                var valueStr = value !== null && value !== undefined ? value.toFixed(2) + 's' : 'N/A';
                                result += marker + param.seriesName + ': ' + valueStr + '<br/>';
                            }
                        });
                        
                        return result;
                    }
                """,
                },
                "xAxis": {
                    "type": "time",
                    "axisLabel": {
                        "color": "#86868B",
                        "fontSize": 12,
                        "formatter": "{HH}:{mm}:{ss}",
                    },
                },
                "yAxis": {
                    "type": "value",
                    "name": "Time (seconds)",
                    "nameTextStyle": {
                        "color": "#86868B",
                        "fontSize": 12,
                        "padding": [0, 30, 0, 0],
                    },
                    "axisLabel": {"color": "#86868B"},
                    "splitLine": {
                        "show": True,
                        "lineStyle": {"type": "dashed", "color": "#E5E5EA"},
                    },
                },
                "series": [
                    {
                        "name": "BAM Merge Time",
                        "type": "line",
                        "smooth": True,
                        "showSymbol": True,
                        "connectNulls": True,
                        "data": [],
                        "lineStyle": {"width": 2, "color": "#007AFF"},
                        "itemStyle": {"color": "#007AFF"},
                    },
                    {
                        "name": "DataFrame Merge Time",
                        "type": "line",
                        "smooth": True,
                        "showSymbol": True,
                        "connectNulls": True,
                        "data": [],
                        "lineStyle": {"width": 2, "color": "#34C759"},
                        "itemStyle": {"color": "#34C759"},
                    },
                    {
                        "name": "Total Processing Time",
                        "type": "line",
                        "smooth": True,
                        "showSymbol": True,
                        "connectNulls": True,
                        "data": [],
                        "lineStyle": {"width": 2, "color": "#FF9500"},
                        "itemStyle": {"color": "#FF9500"},
                    },
                ],
            }
        ).style("height: 300px")

    def update_ui(self) -> None:
        """Update the performance metrics UI with current data."""
        if not self.ui_initialized:
            return

        try:
            # Update stats cards
            for metric_key, card_data in self.stats_cards.items():
                if metric_key in self.metrics and self.metrics[metric_key]:
                    # Filter out None values before calculating statistics
                    values = np.array(
                        [v for v in self.metrics[metric_key] if v is not None]
                    )
                    if len(values) > 0:
                        # Update only the value labels
                        card_data["values"]["min"].text = f"{np.min(values):.2f}"
                        card_data["values"]["max"].text = f"{np.max(values):.2f}"
                        card_data["values"]["mean"].text = f"{np.mean(values):.2f}"
                        card_data["values"]["median"].text = f"{np.median(values):.2f}"
                    else:
                        # If no valid values, show dashes
                        card_data["values"]["min"].text = "--"
                        card_data["values"]["max"].text = "--"
                        card_data["values"]["mean"].text = "--"
                        card_data["values"]["median"].text = "--"

            # Update time series chart
            if "timestamps" in self.metrics and self.metrics["timestamps"]:
                # Filter out None values and create valid data points
                valid_indices = [
                    i for i, t in enumerate(self.metrics["timestamps"]) if t is not None
                ]

                if valid_indices:
                    timestamps = [self.metrics["timestamps"][i] for i in valid_indices]

                    # Get valid data points for each series
                    series_data = []
                    for series_name in [
                        "bam_merge_times",
                        "df_merge_times",
                        "total_processing_times",
                    ]:
                        if series_name in self.metrics:
                            series_values = self.metrics[series_name]
                            # Create data points, using null for missing values
                            data = []
                            for i, ts in enumerate(timestamps):
                                if i < len(series_values):
                                    data.append([ts, series_values[i]])
                                else:
                                    data.append([ts, None])
                            series_data.append(data)
                        else:
                            series_data.append([])

                    # Update each series data in the chart options
                    for i, data in enumerate(series_data):
                        if i < len(self.performance_chart.options["series"]):
                            self.performance_chart.options["series"][i]["data"] = data

                    # Update the chart with the modified options
                    self.performance_chart.update()

        except Exception as e:
            logger.error(f"Error updating performance metrics UI: {e}")

    def get_stats(self) -> Dict[str, Dict[str, float]]:
        """
        Calculate statistics for all metrics.

        Returns
        -------
        Dict[str, Dict[str, float]]
            Dictionary containing statistics for each metric
        """
        stats = {}
        for key, values in self.metrics.items():
            if isinstance(values, list) and values:
                values_array = np.array(values)
                stats[key] = {
                    "min": float(np.min(values_array)),
                    "max": float(np.max(values_array)),
                    "mean": float(np.mean(values_array)),
                    "median": float(np.median(values_array)),
                }
        return stats
