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
from typing import Dict, List, Optional, Union, Tuple
import psutil
import subprocess

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
            "python_ram_usage": [],  # New metric for Python RAM usage
            "r_ram_usage": [],      # New metric for R RAM usage
            "ram_timestamps": [],    # New timestamps for RAM measurements
        }

        # Update with any provided metrics
        if metrics:
            for key, value in metrics.items():
                if key in self.metrics:
                    self.metrics[key] = value

        # UI components
        self.performance_container = None
        self.performance_chart = None
        self.ram_chart = None  # New chart for RAM usage
        self.stats_cards = {}
        self.ui_initialized = False
        self.ram_tracking_timer = None  # Timer for RAM tracking

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

    def start_ram_tracking(self, interval_seconds: int = 10) -> None:
        """
        Start tracking RAM usage at regular intervals.

        Parameters
        ----------
        interval_seconds : int
            Interval in seconds between RAM usage measurements
        """
        if self.ram_tracking_timer is None:
            self.ram_tracking_timer = ui.timer(interval_seconds, self._track_ram_usage)
            logger.info(f"Started RAM tracking with {interval_seconds}s interval")

    def stop_ram_tracking(self) -> None:
        """Stop tracking RAM usage."""
        if self.ram_tracking_timer is not None:
            self.ram_tracking_timer.active = False
            self.ram_tracking_timer = None
            logger.info("Stopped RAM tracking")

    def _get_r_memory_usage(self) -> Optional[float]:
        """
        Get memory usage of R processes in MB.

        Returns
        -------
        Optional[float]
            Total memory usage of R processes in MB, or None if no R processes found
        """
        try:
            # Use ps to find R processes and their memory usage
            result = subprocess.run(
                ["ps", "-eo", "pid,comm,rss"],
                capture_output=True,
                text=True,
                check=True
            )
            
            # Parse output to find R processes
            total_r_memory = 0
            for line in result.stdout.splitlines()[1:]:  # Skip header
                parts = line.split()
                if len(parts) >= 3 and "R" in parts[1]:
                    # RSS is in KB, convert to MB
                    total_r_memory += float(parts[2]) / 1024
            
            return total_r_memory if total_r_memory > 0 else None
        except Exception as e:
            logger.error(f"Error getting R memory usage: {e}")
            return None

    def _track_ram_usage(self) -> None:
        """Track RAM usage of Python and R processes."""
        try:
            current_time = time.time() * 1000  # Convert to milliseconds
            
            # Get Python process memory usage
            python_memory = psutil.Process().memory_info().rss / (1024 * 1024)  # Convert to MB
            
            # Get R process memory usage
            r_memory = self._get_r_memory_usage()
            
            # Record metrics
            self.metrics["python_ram_usage"].append(python_memory)
            self.metrics["r_ram_usage"].append(r_memory)
            self.metrics["ram_timestamps"].append(current_time)
            
            # Save metrics
            self.save_metrics()
            
            # Update UI if initialized
            if self.ui_initialized and self.ram_chart is not None:
                self._update_ram_chart()
                
        except Exception as e:
            logger.error(f"Error tracking RAM usage: {e}")

    def create_ui(self, parent=None) -> None:
        """Create the performance metrics UI components. If parent is provided, use it as the context."""
        context = parent if parent is not None else ui.card().classes("w-full p-4 mb-4")
        with context as self.performance_container:
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

            # Time series charts
            with ui.grid(columns=2).classes("w-full gap-4"):
                # Processing times chart
                with ui.card().classes("w-full p-4 border rounded-lg shadow-sm"):
                    self._create_time_series_chart()
                
                # RAM usage chart
                with ui.card().classes("w-full p-4 border rounded-lg shadow-sm"):
                    self._create_ram_chart()

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

    def _create_ram_chart(self) -> None:
        """Create the RAM usage time series chart."""
        self.ram_chart = ui.echart(
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
                    "text": "RAM Usage Over Time",
                    "left": "center",
                    "top": 5,
                    "textStyle": {
                        "fontSize": 16,
                        "fontWeight": "500",
                        "color": "#1D1D1F",
                    },
                },
                "legend": {
                    "data": ["Python RAM Usage", "R RAM Usage"],
                    "top": 35,
                    "textStyle": {"color": "#86868B"},
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
                                var valueStr = value !== null && value !== undefined ? value.toFixed(2) + ' MB' : 'N/A';
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
                    "name": "Memory (MB)",
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
                        "name": "Python RAM Usage",
                        "type": "line",
                        "smooth": True,
                        "showSymbol": True,
                        "connectNulls": True,
                        "data": [],
                        "lineStyle": {"width": 2, "color": "#007AFF"},
                        "itemStyle": {"color": "#007AFF"},
                    },
                    {
                        "name": "R RAM Usage",
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

    def _update_ram_chart(self) -> None:
        """Update the RAM usage chart with current data."""
        if not self.ui_initialized or self.ram_chart is None:
            return

        try:
            if "ram_timestamps" in self.metrics and self.metrics["ram_timestamps"]:
                # Filter out None values and create valid data points
                valid_indices = [
                    i for i, t in enumerate(self.metrics["ram_timestamps"]) if t is not None
                ]

                if valid_indices:
                    timestamps = [self.metrics["ram_timestamps"][i] for i in valid_indices]

                    # Get valid data points for each series
                    python_data = []
                    r_data = []
                    
                    for i, ts in enumerate(timestamps):
                        if i < len(self.metrics["python_ram_usage"]):
                            python_data.append([ts, self.metrics["python_ram_usage"][i]])
                        else:
                            python_data.append([ts, None])
                            
                        if i < len(self.metrics["r_ram_usage"]):
                            r_data.append([ts, self.metrics["r_ram_usage"][i]])
                        else:
                            r_data.append([ts, None])

                    # Update chart series data
                    self.ram_chart.options["series"][0]["data"] = python_data
                    self.ram_chart.options["series"][1]["data"] = r_data

                    # Update the chart
                    self.ram_chart.update()

        except Exception as e:
            logger.error(f"Error updating RAM chart: {e}")

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

            # Update RAM chart
            self._update_ram_chart()

        except Exception as e:
            logger.error(f"Error updating performance metrics UI: {e}")

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
                    "python_ram_mb": self.metrics["python_ram_usage"],
                    "r_ram_mb": self.metrics["r_ram_usage"],
                    "ram_timestamp": self.metrics["ram_timestamps"],
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
