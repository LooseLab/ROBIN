"""
MNPFlex Analysis Object for ROBIN

This module provides the MNPFlex analysis object for ROBIN, which handles integration with
the MNP-FLEX platform (https://mnp-flex.org/) for methylation analysis.

Key Components:

1. **MNPFlex_Object Class**:
   - Handles authentication and communication with MNP-FLEX API
   - Processes methylation data and sends to MNP-FLEX
   - Manages report downloads and data extraction
   - Provides visualization of results

Dependencies:
- `nicegui` (ui, app)
- `robin.utilities.mnp_flex.APIClient`
- `robin.reporting.pdf_extractor.PDFExtractor`
"""

import os
import logging
import json
from datetime import datetime
from typing import Optional, Dict, Any

from nicegui import ui

from robin.subpages.base_analysis import BaseAnalysis
from robin.utilities.mnp_flex import APIClient as MnpFlexClient
from robin.reporting.pdf_extractor import PDFExtractor
from robin.core.state import state, ProcessState

# Use the main logger configured in the main application
logger = logging.getLogger(__name__)


class MNPFlex_Object(BaseAnalysis):
    """
    Class for handling MNP-FLEX methylation analysis integration.

    This class extends BaseAnalysis to provide specialized functionality for
    MNP-FLEX methylation analysis. It handles authentication, data submission,
    and result visualization.

    Parameters
    ----------
    *args
        Variable length argument list passed to BaseAnalysis
    mnpflex_config : dict, optional
        Configuration dictionary containing MNP-FLEX credentials
    **kwargs
        Arbitrary keyword arguments passed to BaseAnalysis

    Attributes
    ----------
    mnpflex_config : dict
        Configuration for MNP-FLEX API access
    mnpflex_client : MnpFlexClient
        Client for interacting with MNP-FLEX API
    classification_charts : dict
        Dictionary of charts for different classification levels
    mgmt_chart : ui.echart
        Chart for MGMT methylation visualization
    """

    def __init__(
        self, *args, mnpflex_config: Optional[Dict[str, Any]] = None, **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.mnpflex_config = mnpflex_config or {}
        self.mnpflex_client = None
        self.classification_charts = {}
        self.mgmt_chart = None

        # Initialize MNP-FLEX client if config is provided
        if self.mnpflex_config and all(self.mnpflex_config.values()):
            self.initialize_mnpflex_client()

    def initialize_mnpflex_client(self) -> None:
        """Initialize the MNP-FLEX client with provided credentials."""
        try:
            self.mnpflex_client = MnpFlexClient(
                base_url="https://mnp-flex.org", verify_ssl=True
            )
            self.mnpflex_client.authenticate(
                username=self.mnpflex_config["mnpuser"],
                password=self.mnpflex_config["mnppass"],
                client_id="ROBIN",
                client_secret="SECRET",
            )
            logger.info("MNP-FLEX client initialized successfully")
        except Exception as e:
            logger.error(f"Failed to initialize MNP-FLEX client: {str(e)}")
            self.mnpflex_client = None

    def create_dynamic_classification_chart(
        self, key: str, y_max: float = 1
    ) -> ui.echart:
        """Create a chart for a specific classification key with dynamic y-axis max."""
        return ui.echart(
            {
                "backgroundColor": "transparent",
                "title": {
                    "text": f"{key.title()} Classification Over Time",
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
                    "top": 45,
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
                        "formatter": "{yyyy}-{MM}-{dd} {HH}:{mm}",
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
                    "max": y_max,
                    "interval": y_max / 5 if y_max > 0 else 1,
                    "axisLabel": {
                        "fontSize": 12,
                        "formatter": "{value}",
                        "color": "#666666",
                    },
                    "splitLine": {
                        "show": True,
                        "lineStyle": {"color": "#E0E0E0", "type": "dashed"},
                    },
                },
                "series": [],
            }
        ).classes("w-full h-64")

    def create_mgmt_chart(self) -> ui.echart:
        """Create a chart for displaying MGMT methylation trends."""
        return ui.echart(
            {
                "backgroundColor": "transparent",
                "title": {
                    "text": "MGMT Methylation Over Time",
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
                "xAxis": {
                    "type": "time",
                    "axisLabel": {
                        "fontSize": 12,
                        "formatter": "{yyyy}-{MM}-{dd} {HH}:{mm}",
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
                "series": [
                    {
                        "name": "Methylation",
                        "type": "line",
                        "smooth": True,
                        "animation": False,
                        "symbolSize": 6,
                        "emphasis": {
                            "focus": "series",
                            "itemStyle": {"borderWidth": 2},
                        },
                        "lineStyle": {"width": 2, "color": "#007AFF"},
                        "itemStyle": {"color": "#007AFF"},
                        "data": [],
                    }
                ],
            }
        ).classes("w-full h-64")

    def update_summary_card(self):
        """Update the summary card with the latest report and classification information."""
        if not (
            self._summary_status_label
            and self._summary_reports_label
            and self._summary_last_label
        ):
            return
        sample_dir = self.output
        try:
            mnpflex_reports = [
                f
                for f in os.listdir(sample_dir)
                if f.endswith(".pdf") and "mnpflex" in f.lower()
            ]
            if mnpflex_reports:
                self._summary_status_label.text = "Status: Results Available"
                self._summary_reports_label.text = (
                    f"Reports Generated: {len(mnpflex_reports)}"
                )
                # Find latest report
                latest_report = max(
                    mnpflex_reports,
                    key=lambda f: os.path.getmtime(os.path.join(sample_dir, f)),
                )
                latest_time = os.path.getmtime(os.path.join(sample_dir, latest_report))
                dt = datetime.fromtimestamp(latest_time)
                self._summary_last_label.text = (
                    f"Last Analysis: {dt.strftime('%Y-%m-%d %H:%M:%S')}"
                )
            else:
                self._summary_status_label.text = "Status: Awaiting Data"
                self._summary_reports_label.text = "Reports Generated: 0"
                self._summary_last_label.text = "Last Analysis: --"

            # Now update latest classification and mgmt status
            json_path = os.path.join(
                self.output, "extracted_data", "extracted_data.json"
            )
            if os.path.exists(json_path):
                with open(json_path, "r") as f:
                    data = json.load(f)
                if data:
                    # Find the latest entry by timestamp
                    latest_entry = max(data, key=lambda x: x.get("timestamp", ""))
                    # Update classification labels
                    for key in ["superfamily", "family", "class", "subclass"]:
                        label = "--"
                        score = ""
                        if (
                            "classification" in latest_entry
                            and key in latest_entry["classification"]
                        ):
                            label = latest_entry["classification"][key].get(
                                "label", "--"
                            )
                            score_val = latest_entry["classification"][key].get("score")
                            if score_val is not None:
                                score = f" ({score_val:.2f})"
                        if key in self._summary_classification_labels:
                            self._summary_classification_labels[key].text = (
                                f"{key.title()}: {label}{score}"
                            )
                    # Update MGMT label
                    mgmt = latest_entry.get("mgmt_info", {})
                    mgmt_status = mgmt.get("status", "--")
                    avg = mgmt.get("average_methylation")
                    if avg is not None:
                        self._summary_mgmt_label.text = (
                            f"MGMT: {mgmt_status} (Avg: {avg:.1f}%)"
                        )
                    else:
                        self._summary_mgmt_label.text = f"MGMT: {mgmt_status}"
                else:
                    for key in self._summary_classification_labels:
                        self._summary_classification_labels[key].text = (
                            f"{key.title()}: --"
                        )
                    if self._summary_mgmt_label:
                        self._summary_mgmt_label.text = "MGMT: --"
            else:
                for key in self._summary_classification_labels:
                    self._summary_classification_labels[key].text = f"{key.title()}: --"
                if self._summary_mgmt_label:
                    self._summary_mgmt_label.text = "MGMT: --"
        except Exception as e:
            self._summary_status_label.text = f"Status: Error - {e}"
            self._summary_reports_label.text = "Reports Generated: --"
            self._summary_last_label.text = "Last Analysis: --"
            for key in self._summary_classification_labels:
                self._summary_classification_labels[key].text = f"{key.title()}: --"
            if self._summary_mgmt_label:
                self._summary_mgmt_label.text = "MGMT: --"

    def update_extracted_data_charts(self) -> None:
        """Update dynamic charts with data from extracted_data.json."""
        self.update_summary_card()
        try:
            logger.info("Starting chart update")
            if not hasattr(self, "classification_charts") or not hasattr(
                self, "mgmt_chart"
            ):
                logger.warning("Charts not initialized")
                return

            json_path = os.path.join(
                self.output, "extracted_data", "extracted_data.json"
            )
            logger.info(f"Looking for data file at: {json_path}")

            if not os.path.exists(json_path):
                logger.info(f"No data file found at {json_path}")
                return

            with open(json_path, "r") as f:
                data = json.load(f)

            if not data:
                logger.info("No data found in JSON file")
                return

            logger.info(f"Processing {len(data)} data points")

            # Build a dict: key -> label -> list of (timestamp, score)
            key_label_series = {}
            mgmt_data = []
            colors = [
                "#007AFF",
                "#34C759",
                "#FF9500",
                "#FF2D55",
                "#5856D6",
                "#FF3B30",
                "#5AC8FA",
                "#4CD964",
            ]

            for entry in data:
                timestamp = entry.get("timestamp")
                # Handle classification data
                for key, info in entry.get("classification", {}).items():
                    label = info.get("label")
                    score = info.get("score")
                    if label is None or score is None:
                        continue
                    if key not in key_label_series:
                        key_label_series[key] = {}
                    if label not in key_label_series[key]:
                        key_label_series[key][label] = []
                    key_label_series[key][label].append([timestamp, score])

                # Handle MGMT data
                mgmt = entry.get("mgmt_info", {})
                avg = mgmt.get("average_methylation")
                if avg is not None and timestamp is not None:
                    mgmt_data.append([timestamp, avg])

            # Sort all series by timestamp
            for key in key_label_series:
                for label in key_label_series[key]:
                    key_label_series[key][label].sort()
            mgmt_data.sort()

            # Update classification charts
            for key, label_dict in key_label_series.items():
                all_scores = [
                    score for points in label_dict.values() for _, score in points
                ]
                y_max = max(all_scores) if all_scores else 1
                y_max = max(1, int(y_max + 0.5))  # Round up for clarity

                if key not in self.classification_charts:
                    self.classification_charts[key] = (
                        self.create_dynamic_classification_chart(key, y_max=y_max)
                    )

                chart = self.classification_charts[key]
                chart.options["yAxis"]["max"] = y_max
                chart.options["yAxis"]["interval"] = y_max / 5 if y_max > 0 else 1

                chart_series = []
                for idx, (label, data_points) in enumerate(label_dict.items()):
                    chart_series.append(
                        {
                            "name": label,
                            "type": "line",
                            "smooth": True,
                            "animation": False,
                            "symbolSize": 6,
                            "emphasis": {
                                "focus": "series",
                                "itemStyle": {"borderWidth": 2},
                            },
                            "lineStyle": {
                                "width": 2,
                                "color": colors[idx % len(colors)],
                            },
                            "itemStyle": {"color": colors[idx % len(colors)]},
                            "data": data_points,
                        }
                    )
                chart.options["series"] = chart_series
                chart.update()

            # Update MGMT chart
            if mgmt_data:
                logger.info("Updating MGMT chart with %d points", len(mgmt_data))
                self.mgmt_chart.options["series"][0]["data"] = mgmt_data
                self.mgmt_chart.update()

        except Exception as e:
            logger.error(
                f"Error updating extracted data charts: {str(e)}", exc_info=True
            )

    async def process_bam(self, bamfile: str, timestamp: float) -> None:
        """Process a BAM file and submit to MNP-FLEX."""
        try:
            if not self.mnpflex_client:
                logger.warning("MNP-FLEX client not initialized, skipping processing")
                return

            # Get sample ID
            sample_id = self.force_sampleid if self.force_sampleid else self.sampleID
            if not sample_id:
                logger.error("No sample ID available")
                return

            # Prepare the file for MNP-FLEX
            clean_file = await self.prepare_file_for_mnpflex(bamfile)
            if not clean_file:
                logger.error("Failed to prepare file for MNP-FLEX")
                return

            # Upload to MNP-FLEX
            response = self.mnpflex_client.upload_sample(
                file_path=clean_file,
                sample_name=f"{sample_id}.mnpFlex",
                disclaimer_confirmed=True,
            )

            if not response or "id" not in response:
                logger.error("Failed to upload sample to MNP-FLEX")
                return

            sample_id = response["id"]
            sample = self.mnpflex_client.get_sample(sample_id)
            result_status = sample["bed_file_sample"]["analysis_status"]

            # Wait for analysis to complete
            while result_status == "initialized":
                await ui.run_javascript(
                    "await new Promise(resolve => setTimeout(resolve, 1000))"
                )
                sample = self.mnpflex_client.get_sample(sample_id)
                result_status = sample["bed_file_sample"]["analysis_status"]

            if result_status == "Analysis error":
                logger.error(f"Analysis error for {sample_id}")
                self.mnpflex_client.delete_sample(sample_id)
                return

            # Get and save report
            report_content = self.mnpflex_client.get_sample_report(sample_id)
            if isinstance(report_content, bytes):
                report_path = os.path.join(
                    self.output,
                    f'{sample["sample_name"]}_{datetime.now().strftime("%Y%m%d_%H%M%S")}.pdf',
                )
                with open(report_path, "wb") as file:
                    file.write(report_content)
                logger.info(f"Report saved as {report_path}")

                # Extract and store data from the report
                extractor = PDFExtractor(os.path.join(self.output, "extracted_data"))
                data = extractor.extract_from_pdf(report_path)
                if data:
                    extractor.save_data(data)
                    logger.info(f"Extracted and stored data from report: {report_path}")

            # Clean up
            self.mnpflex_client.delete_sample(sample_id)

        except Exception as e:
            logger.error(f"Error processing BAM file: {str(e)}")

    async def prepare_file_for_mnpflex(self, bamfile: str) -> Optional[str]:
        """Prepare a file for MNP-FLEX submission."""
        # Implementation of file preparation logic
        # This would include any necessary file format conversion or cleaning
        # Return the path to the prepared file
        return bamfile

    def setup_ui(self) -> None:
        """Set up the user interface for MNP-FLEX analysis."""
        # Add summary card
        self._summary_status_label = None
        self._summary_reports_label = None
        self._summary_last_label = None
        self._summary_classification_labels = {}
        self._summary_mgmt_label = None
        if self.summary:
            with self.summary:
                with ui.card().classes("w-full p-4 mb-4"):
                    with ui.row().classes("w-full items-center justify-between"):
                        # Left side - MNP-FLEX Status and Results
                        with ui.column().classes("gap-2"):
                            ui.label("MNP-FLEX Methylation Analysis").classes(
                                "text-lg font-medium"
                            )
                            # RESEARCH USE ONLY indicator
                            with ui.row().classes("items-center gap-2"):
                                ui.icon("warning").classes("text-red-600 text-base")
                                ui.label("RESEARCH USE ONLY").classes(
                                    "text-red-700 font-semibold text-sm"
                                )
                            with ui.row().classes("items-center gap-2"):
                                self._summary_status_label = ui.label(
                                    "Status: Awaiting Data"
                                ).classes("text-gray-600")
                                ui.label("--").classes(
                                    "px-2 py-1 rounded bg-gray-100 text-gray-600"
                                )
                            # Latest classification results
                            for key in ["superfamily", "family", "class", "subclass"]:
                                self._summary_classification_labels[key] = ui.label(
                                    f"{key.title()}: --"
                                ).classes("text-gray-700 text-sm")
                            self._summary_mgmt_label = ui.label("MGMT: --").classes(
                                "text-gray-700 text-sm"
                            )

                        # Right side - Analysis metrics
                        with ui.column().classes("gap-2 text-right"):
                            ui.label("Analysis Details").classes("font-medium")
                            self._summary_reports_label = ui.label(
                                "Reports Generated: --"
                            ).classes("text-gray-600")
                            self._summary_last_label = ui.label(
                                "Last Analysis: --"
                            ).classes("text-gray-600")

                    # Bottom row - Information
                    with ui.row().classes(
                        "w-full mt-4 text-sm text-gray-500 justify-center"
                    ):
                        ui.label(
                            "Methylation-based classification using MNP-FLEX platform"
                        )

        # Add warning banner
        with ui.card().classes("w-full mb-4 bg-red-50 border-2 border-red-200"):
            with ui.row().classes("items-center gap-2"):
                ui.icon("warning").classes("text-red-600 text-2xl")
                with ui.column().classes("gap-1"):
                    ui.label("⚠️ IMPORTANT WARNING").classes(
                        "text-red-800 font-bold text-lg"
                    )
                    ui.label("MNP-FLEX is for RESEARCH USE ONLY").classes(
                        "text-red-700 font-semibold"
                    )
                    ui.label(
                        "Users must have their own account with Epignostix to use this service"
                    ).classes("text-red-700")
                    ui.link(
                        "https://mnp-flex.org/", "Visit MNP-FLEX website", new_tab=True
                    ).classes("text-red-600 hover:text-red-800")

        with ui.card().classes("w-full mb-4"):
            ui.label("About MNP-FLEX").classes("text-lg font-semibold mb-2")
            ui.label(
                "MNP-FLEX is a methylation analysis platform that provides detailed insights into DNA methylation patterns. Visit "
            ).classes("text-gray-700")
            ui.link(
                "https://mnp-flex.org/", "https://mnp-flex.org/", new_tab=True
            ).classes("text-blue-600 hover:text-blue-800")
            ui.label(" for more information.").classes("text-gray-700")

            with ui.column().classes("mt-4 space-y-2"):
                ui.label("Important Notes:").classes("font-semibold text-gray-700")
                ui.label(
                    "• Users must have their own MNP-FLEX account to use this service"
                ).classes("text-gray-600")
                ui.label("• All results are for research use only").classes(
                    "text-gray-600"
                )
                ui.label(
                    "• Results are automatically downloaded as PDF reports"
                ).classes("text-gray-600")

            with ui.column().classes("mt-4 space-y-2"):
                ui.label("Reference:").classes("font-semibold text-gray-700")
                ui.label("MNP-FLEX is based on the following publication:").classes(
                    "text-gray-600"
                )
                with ui.row().classes("items-start gap-2"):
                    ui.icon("article").classes("text-blue-600 mt-1")
                    with ui.column():
                        ui.link(
                            "https://www.nature.com/articles/s41591-025-03562-5",
                            "MNP-FLEX: A methylation-based platform for rapid and accurate classification of brain tumors",
                            new_tab=True,
                        ).classes("text-blue-600 hover:text-blue-800")
                        ui.label("Nature Medicine (2025)").classes(
                            "text-gray-600 text-sm"
                        )
                        ui.label("DOI: 10.1038/s41591-025-03562-5").classes(
                            "text-gray-600 text-sm"
                        )

        # Create refreshable container for reports
        @ui.refreshable
        def show_mnpflex_reports():
            sample_dir = self.output
            try:
                with ui.card().classes("w-full mb-4"):
                    with ui.row().classes("items-center gap-2"):
                        ui.icon("folder").classes("text-blue-600")
                        ui.label("Searching for reports in:").classes("text-gray-600")
                    ui.label(sample_dir).classes(
                        "text-sm font-mono text-gray-700 ml-8 break-all"
                    )

                mnpflex_reports = [
                    f
                    for f in os.listdir(sample_dir)
                    if f.endswith(".pdf") and "mnpflex" in f.lower()
                ]
                if mnpflex_reports:
                    # Parse and sort reports
                    report_info = []
                    for report in mnpflex_reports:
                        try:
                            parts = report.split("_")
                            sample_id = parts[0]

                            if (
                                len(parts) >= 5
                                and len(parts[3]) == 8
                                and len(parts[4].replace(".pdf", "")) == 6
                            ):
                                sites = parts[2]
                                date = parts[3]
                                time = parts[4].replace(".pdf", "")
                                dt = datetime.strptime(
                                    f"{date}_{time}", "%Y%m%d_%H%M%S"
                                )
                                formatted_date = dt.strftime("%Y-%m-%d %H:%M:%S")
                            else:
                                file_path = os.path.join(sample_dir, report)
                                dt = datetime.fromtimestamp(os.path.getmtime(file_path))
                                formatted_date = dt.strftime("%Y-%m-%d %H:%M:%S")
                                try:
                                    sites = parts[1].split(".")[-1]
                                except Exception as e:
                                    print(
                                        f"Error parsing report filename {report}: {e}"
                                    )
                                    sites = "N/A"

                            report_info.append(
                                {
                                    "filename": report,
                                    "sample_id": sample_id,
                                    "sites": sites,
                                    "datetime": dt,
                                    "formatted_date": formatted_date,
                                }
                            )
                        except Exception as e:
                            logger.error(
                                f"Error parsing report filename {report}: {str(e)}"
                            )
                            continue

                    report_info.sort(key=lambda x: x["datetime"], reverse=True)

                    # Display latest report
                    if report_info:
                        latest = report_info[0]
                        with ui.card().classes("w-full mb-4 bg-blue-50"):
                            ui.label("Latest Report").classes(
                                "text-lg font-semibold text-blue-800"
                            )
                            with ui.row().classes("items-center gap-2"):
                                ui.icon("description").classes("text-blue-600")
                                ui.label(f"Sample: {latest['sample_id']}").classes(
                                    "text-blue-800"
                                )
                            with ui.row().classes("items-center gap-2"):
                                ui.icon("schedule").classes("text-blue-600")
                                ui.label(
                                    f"Generated: {latest['formatted_date']}"
                                ).classes("text-blue-800")
                            with ui.row().classes("items-center gap-2"):
                                ui.icon("analytics").classes("text-blue-600")
                                ui.label(f"Sites analyzed: {latest['sites']}").classes(
                                    "text-blue-800"
                                )
                            ui.button(
                                "Download Latest Report",
                                on_click=lambda r=os.path.join(
                                    sample_dir, latest["filename"]
                                ): self.download_mnpflex_report(r),
                            ).classes("bg-blue-600 hover:bg-blue-700 mt-2")

                    # Display older reports
                    if len(report_info) > 1:
                        with ui.expansion("Previous Reports", icon="history").classes(
                            "w-full"
                        ):
                            with ui.column().classes("w-full gap-2"):
                                for report in report_info[1:]:
                                    with ui.card().classes("w-full"):
                                        with ui.row().classes(
                                            "items-center justify-between w-full"
                                        ):
                                            with ui.column().classes("gap-1"):
                                                ui.label(
                                                    f"Sample: {report['sample_id']}"
                                                ).classes("text-sm font-medium")
                                                ui.label(
                                                    f"Generated: {report['formatted_date']}"
                                                ).classes("text-sm text-gray-600")
                                                ui.label(
                                                    f"Sites: {report['sites']}"
                                                ).classes("text-sm text-gray-600")
                                            ui.button(
                                                "Download",
                                                on_click=lambda r=os.path.join(
                                                    sample_dir, report["filename"]
                                                ): self.download_mnpflex_report(r),
                                            ).classes("bg-gray-100 hover:bg-gray-200")
                else:
                    ui.label("No MNPFlex reports available yet.").classes(
                        "text-gray-500"
                    )
            except FileNotFoundError:
                ui.label("Sample directory not found.").classes("text-gray-500")

        # Initial display of reports
        show_mnpflex_reports()

        # Set up periodic refresh
        refresh_timer = ui.timer(30.0, show_mnpflex_reports.refresh)

        # Add extracted data visualization section
        with ui.card().classes("w-full p-4 bg-white rounded-lg shadow-sm mt-4"):
            with ui.expansion("Extracted Data Visualization", icon="analytics").classes(
                "w-full"
            ):
                with ui.column().classes("w-full gap-4"):
                    # Create charts
                    for key in ["superfamily", "family", "class", "subclass"]:
                        self.classification_charts[key] = (
                            self.create_dynamic_classification_chart(key)
                        )
                    self.mgmt_chart = self.create_mgmt_chart()

                    # Set up auto-refresh timer
                    ui.timer(10.0, self.update_extracted_data_charts)

        # Clean up timer when tab is closed
        def cleanup():
            refresh_timer.active = False

        ui.on("disconnect", cleanup)

    async def download_mnpflex_report(self, report_path: str) -> None:
        """Download an MNP-FLEX report."""
        if os.path.exists(report_path):
            ui.download(report_path)
        else:
            ui.notify(f"Report not found at {report_path}", type="negative")

    async def stop_analysis(self) -> None:
        """Stop the MNP-FLEX analysis."""
        state.set_process_state("MNPFlex Analysis", ProcessState.STOPPING)
        state.stop_process("MNPFlex Analysis")
        await super().stop_analysis()
