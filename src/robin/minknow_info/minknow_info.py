"""
MinKNOW Information Module

This module provides a comprehensive interface for interacting with Oxford Nanopore's MinKNOW software.
It handles real-time monitoring of sequencing runs, data visualization, and device management through
a modern web interface built with NiceGUI.

The module includes classes for managing MinKNOW connections, monitoring device status,
tracking sequencing progress, and displaying real-time statistics and plots.

Key Features:
- Real-time monitoring of MinKNOW sequencing runs
- Interactive data visualization with plots and statistics
- Device and flowcell status tracking
- Run configuration and management
- Automated data collection and analysis

Dependencies:
- nicegui: For building the web interface
- minknow_api: For communicating with MinKNOW
- numpy: For numerical computations
- logging: For application logging
"""

# Python imports.
from __future__ import annotations
from nicegui import ui, app
from robin import images
from pathlib import Path
import threading
from datetime import datetime
from robin import theme
from minknow_api.manager import Manager
import os
import asyncio
import re
import logging
from contextlib import contextmanager
from collections import deque
import time

# MinKNOW API Imports
from robin.utilities.camera import Camera
from typing import Sequence
import uuid

import cv2
import zxingcpp

import numpy as np


import base64

# We need `find_protocol` to search for the required protocol given a kit + product code.
from minknow_api.tools import protocols

UNIQUE_ID: str = str(uuid.uuid4())


class ExperimentSpec(object):
    """
    A class to hold experiment specifications for a MinKNOW sequencing position.
    
    This class maintains information about a sequencing position and its associated
    protocol ID for experiment execution.
    
    Attributes:
        position: The sequencing position object from MinKNOW
        protocol_id (str): The identifier for the protocol to be run
    """
    
    def __init__(self, position):
        """
        Initialize an ExperimentSpec instance.
        
        Args:
            position: A MinKNOW position object representing a sequencing position
        """
        self.position = position
        self.protocol_id = ""


ExperimentSpecs = Sequence[ExperimentSpec]


# Determine which protocol to run for each experiment, and add its ID to experiment_specs
def add_protocol_ids(experiment_specs, kit, basecall_config, expected_flowcell_id):
    """
    Add protocol IDs to experiment specifications based on flowcell and kit information.
    
    This function validates flowcell presence and compatibility, then finds and assigns
    the appropriate protocol ID for each experiment specification.
    
    Args:
        experiment_specs (list): List of ExperimentSpec objects
        kit (str): The sequencing kit identifier
        basecall_config (str): Basecalling configuration name
        expected_flowcell_id (str): The expected flowcell ID to validate against
        
    Returns:
        bool: True if protocol IDs were successfully added, False otherwise
        
    Raises:
        None: Errors are handled internally and reported via UI notifications
    """
    
    for spec in experiment_specs:
        # Connect to the sequencing position:
        position_connection = spec.position.connect()
        flow_cell_info = position_connection.device.get_flow_cell_info()
        if flow_cell_info.flow_cell_id != expected_flowcell_id:
            ui.notify(
                f"Flowcell {expected_flowcell_id} is not found in position {spec.position}. Please check.",
                type="negative",
            )
            return
        if not flow_cell_info.has_flow_cell:
            ui.notify(
                "No flow cell present in position {}".format(spec.position),
                type="negative",
            )
            return

        product_code = flow_cell_info.user_specified_product_code
        if not product_code:
            product_code = flow_cell_info.product_code

        # Find the protocol identifier for the required protocol:
        protocol_info = protocols.find_protocol(
            position_connection,
            product_code=product_code,
            kit=kit,
            basecalling=True,
            basecall_config=basecall_config,
            barcoding=False,
            barcoding_kits=[],
            force_reload=True,
            experiment_type="sequencing",
        )

        if not protocol_info:
            ui.notify(
                "Failed to find protocol for position %s" % (spec.position),
                type="negative",
            )
            return

        # Store the identifier for later:
        spec.protocol_id = protocol_info.identifier

    return True


@contextmanager
def disable(button: ui.button):
    """
    Context manager to temporarily disable a UI button.
    
    This ensures the button is re-enabled even if an exception occurs.
    
    Args:
        button (ui.button): The button to disable
        
    Yields:
        None
        
    Example:
        with disable(my_button):
            # Perform some operation
            # Button will be re-enabled after the block
    """
    button.disable()
    try:
        yield
    finally:
        button.enable()


# Define the path to the image file used in the header and footer
IMAGEFILE = os.path.join(
    os.path.dirname(os.path.abspath(images.__file__)), "ROBIN_logo_small.png"
)


DEVICEDICT = {
    "FLONGLE": "Flongle",
    "GRIDION": "GridION",
    "MINION": "MinION",
    "PROMETHION": "PromethION",
    "P2_INTEGRATED": "P2Integrated",
    "P2_SOLO": "P2Solo",
}


class ErrorChecker:
    def __init__(self, *elements) -> None:
        self.elements = elements

    @property
    def no_errors(self) -> bool:
        return all(
            validation(element.value)
            for element in self.elements
            for validation in element.validation.values()
        )


def identify_device(device_type, device_name):
    if device_type == "MINION":
        return f"{DEVICEDICT[device_type]}", "MinION-flow-cell-lights-on-whiteBG.jpg"
    elif device_type == "GRIDION":
        return (
            f"{DEVICEDICT[device_type]}",
            "GridION_Front_Square_Elevated_Closed_Flow Cells 2_White.jpg",
        )
    elif device_type == "P2_INTEGRATED":
        return (
            f"{DEVICEDICT[device_type]}",
            "p2_-left-45_screen-up-1_splash_transparent.png",
        )
    elif device_type == "P2_SOLO":
        if device_name.startswith("P2S"):
            return "P2 Solo", "P2_Solo_Left-45_Open_Full.png"
        else:
            return (f"{device_type}", None)
    elif device_type == "PROMETHION":
        return (
            f"{DEVICEDICT[device_type]}",
            "PromethION_24_Right 45 Elevated_Closed_Full_white.jpg",
        )
    else:
        return (f"{device_type}", None)


class Minknow_Info:
    def __init__(
        self,
        position,
        centreID,
        kit,
        reference,
        basecall_config,
        bed_file,
        experiment_duration,
        dev=False,
    ):
        self.dev = dev
        self.color = "text-blue"
        self.position = position
        self.connection = self.position.connect()
        self.device = self.position.device_type
        self.show = False

        # Initialize yield tracking with lists instead of deques to maintain all history
        self.timestamps = []
        self.read_counts = []
        self.pass_reads = []
        self.pass_bases = []
        self.last_update = time.time()
        self.update_interval = 10  # Update plot every 10 seconds
        self.acquisition_run_id = None

        # Initialize current values
        self.Read_Count = 0
        self.Pass_Read_Count = 0
        self.Pass_Bases = 0
        self.Mean_Basecall_Speed = 0
        self.N50 = 0
        self.Estimated_N50 = 0
        self.Basecall_Speed = 0
        self.Channel_Count = 0
        self.Experiment_Group = ""
        self.Sample_ID = ""
        self.Flowcell_Type = ""
        self.running_kit = ""
        self.output_folder = ""
        self.basecalling_config_filename = ""
        self.start_time = None

        self.basecall_config = basecall_config
        self.centreID = centreID
        self.kit = kit
        self.reference = reference
        self.bed_file = bed_file
        self.experiment_duration = experiment_duration

        self.name, self.image = identify_device(
            self.position.device_type, self.position.name
        )
        if self.image:
            self.deviceicon = os.path.join(
                os.path.dirname(os.path.abspath(images.__file__)),
                "ONTimages",
                self.image,
            )
        else:
            self.deviceicon = os.path.join(
                os.path.dirname(os.path.abspath(images.__file__)),
                "ONTimages",
                "unknown_sequencer.png",
            )
        self.check_instance = threading.Thread(
            target=self.stream_instance_activity, args=()
        )
        self.check_instance.daemon = True
        self.check_instance.start()
        self.render_me()
        self._status_updates = []
        ui.timer(1.0, self._update_all_status)

    def render_me(self):
        # Initialize camera handlers
        sample_camera = Camera(
            icon="photo_camera",
            icon_color="primary",
            background_color="rgba(66, 165, 245, 0.3)",
            for_id="samplefileupload",
            canvas_id="samplecanvas",
            on_change=lambda: generate_sampleID(),
        )

        flowcell_camera = Camera(
            icon="photo_camera",
            icon_color="primary",
            background_color="rgba(66, 165, 245, 0.3)",
            for_id="flowcellfileupload",
            canvas_id="flowcellcanvas",
            on_change=lambda: generate_flowcellID(),
        )

        # Create containers for both panels with full width
        with ui.column().classes('w-full h-full'):
            # Setup Panel - only visible when not running
            setup_container = ui.element('div').classes('w-full h-full').bind_visibility_from(self, 'show', value=False)
            monitor_container = ui.element('div').classes('w-full h-full').bind_visibility_from(self, 'show')

            with monitor_container:
                with ui.card().classes('w-full h-full'):
                    # Device Header
                    with ui.row().classes('items-center q-pa-md'):
                        with ui.avatar(size='xl').classes('q-mr-md'):
                            ui.image(self.deviceicon)
                        with ui.column():
                            ui.label(f'{self.name} - {self.position.name}').classes('text-h6')
                            ui.label('MinKNOW Monitoring').classes(f'text-caption {self.color}')
                            ui.label(f'Device Type: {self.position.device_type}').classes('text-body2')

                    # Status Indicators
                    with ui.row().classes('q-pa-md q-gutter-md justify-between'):
                        self._render_status_chip('Running', self, 'show')
                        self._render_status_chip('Basecalling', self, 'Basecall_Speed', 
                                              lambda v: 'positive' if v > 0 else 'negative')

                    # Monitoring Grid
                    with ui.grid(columns=3).classes('q-pa-md q-gutter-md').bind_visibility_from(self, 'show'):
                        # Basic Info
                        self._render_monitoring_tile('Experiment Group', 'Experiment_Group')
                        self._render_monitoring_tile('Sample ID', 'Sample_ID')
                        self._render_monitoring_tile('Flowcell Type', 'Flowcell_Type')
                    
                    # Yield Plot Container
                    with ui.card().classes('w-full q-ma-md').style('min-height: 500px'):
                        with ui.card_section():
                            ui.label('Sequencing Yield').classes('text-h6')
                        
                        # Initialize the yield plot
                        logging.debug("=== Initializing Yield Plot ===")
                        logging.debug("Initial timestamps: %d", len(self.timestamps))
                        logging.debug("Initial read counts: %d", len(self.read_counts))
                        logging.debug("Initial pass reads: %d", len(self.pass_reads))
                        logging.debug("Initial pass bases: %d", len(self.pass_bases))
                        initial_options = self.create_yield_plot()
                        self.yield_plot = ui.echart(options=initial_options).classes('w-full').style('height: 400px')
                        logging.debug("Plot initialized with ID: %s", self.yield_plot.id)
                        logging.debug("=== Plot Initialization Complete ===")

            with setup_container:
                with ui.card().classes('w-full h-full'):
                    with ui.card_section():
                        ui.label("Sample Setup").classes('text-h5 text-primary q-mb-md')
                    
                    with ui.card_section().classes('q-pa-md'):
                        with ui.stepper().props('vertical').classes('w-full') as stepper:
                            # Sample ID Step
                            with ui.step('Sample ID').classes('text-body1 text-weight-medium'):
                                ui.label('Sample ID Guidelines').classes('text-subtitle1 text-weight-medium q-mb-sm')
                                ui.label('IDs should not include human identifiable information').classes('text-caption q-mb-md')
                                
                                with ui.row().classes('items-center w-full q-gutter-md'):
                                    sampleid = ui.input(
                                        placeholder='Enter sample ID',
                                        validation={
                                            "Too short": lambda value: len(value) >= 5,
                                            "No Flowcell IDs": lambda value: re.match(r"^[A-Za-z]{3}\d{5}$", value) is None,
                                            "Alphanumeric only": lambda value: value.isalnum() and not any(char.isspace() for char in value),
                                        }
                                    ).props('outlined dense').classes('w-full')
                                    
                                    sample_camera.show_camera()
                                    
                                with ui.stepper_navigation():
                                    checker = ErrorChecker(sampleid)
                                    ui.button('Next', on_click=stepper.next).bind_enabled_from(checker, 'no_errors')

                            # Flowcell ID Step
                            with ui.step('Flowcell').classes('text-body1 text-weight-medium'):
                                ui.label('Flowcell ID Entry').classes('text-subtitle1 text-weight-medium q-mb-sm')
                                
                                with ui.row().classes('items-center w-full q-gutter-md'):
                                    flowcellid = ui.input(
                                        placeholder='Enter flowcell ID',
                                        validation={
                                            "Not a valid flowcell ID": lambda value: re.match(
                                                r"^[A-Za-z]{3}\d{5}$", value
                                            ) is not None
                                        },
                                    ).props('outlined dense').classes('w-full')
                                    
                                    flowcell_camera.show_camera()
                                
                                with ui.stepper_navigation():
                                    ui.button('Back', on_click=stepper.previous).props('flat')
                                    checkerflowcell = ErrorChecker(flowcellid)
                                    ui.button('Next', on_click=stepper.next).bind_enabled_from(
                                        checkerflowcell, 'no_errors'
                                    )

                            # Device Position Step
                            with ui.step('Device Position').classes('text-body1 text-weight-medium'):
                                ui.label('Position Selection').classes('text-subtitle1 text-weight-medium q-mb-sm')
                                ui.label('Select the correct device position').classes('text-caption q-mb-md')
                                
                                run_button = ui.button(
                                    'Start Run',
                                    on_click=lambda: self.start_run(
                                        position=self.position.name,
                                        reference=self.reference,
                                        sample_id=sampleid.value,
                                        flowcell_id=flowcellid.value,
                                        kit=self.kit,
                                        basecall_config=self.basecall_config,
                                        centreID=self.centreID,
                                        experiment_duration=self.experiment_duration,
                                        bed_file=self.bed_file,
                                    )
                                ).props('color=primary')

                                ui.radio(
                                    [self.position.name], 
                                    value=self.position.name
                                ).on('update:model-value', lambda: run_button.enable()
                                ).classes('q-mt-md')

                                with ui.stepper_navigation():
                                    ui.button('Back', on_click=stepper.previous).props('flat')
                                    run_button.disable()

                        # Settings Summary
                        with ui.expansion('Fixed Settings', icon='settings').classes('w-full q-mt-md'):
                            with ui.card().classes('q-pa-md'):
                                self._render_setting_item('Device', self.device)
                                current_date = datetime.now()
                                self._render_setting_item(
                                    'Experiment Group ID', 
                                    f"{self.centreID}_{current_date.strftime('%B')}_{current_date.year}"
                                )
                                self._render_setting_item('Kit', self.kit)
                                self._render_setting_item('Reference', self.reference)

    def _render_setting_item(self, label: str, value: str):
        """Helper method to render consistent setting items"""
        with ui.row().classes('items-center justify-between w-full q-py-sm'):
            ui.label(label).classes('text-body2')
            ui.label(str(value)).classes('text-body2 text-weight-medium')

    def _render_monitoring_tile(self, label: str, bind_property: str, format_func=str):
        """Helper method to render monitoring information tiles"""
        with ui.card().classes('q-pa-sm'):
            ui.label(label).classes(f'text-caption {self.color}')
            ui.label('--').bind_text_from(
                self, bind_property,
                lambda n: format_func(n) if n is not None else '--'
            ).classes('text-body1 text-weight-medium')

    def _render_status_chip(self, label: str, obj, bind_property: str, 
                           color_func=lambda v: 'positive' if v else 'negative'):
        """Helper method to render status indicator chips"""
        status_button = ui.button(
            label,
            icon='circle'
        ).props('flat dense').classes('q-ma-xs')
        
        # Store the update function and button in a list for batch updates
        if not hasattr(self, '_status_updates'):
            self._status_updates = []
        
        def update_status():
            if hasattr(obj, bind_property):
                value = getattr(obj, bind_property)
                color = color_func(value)
                status_button.props(f'color={color}')
        
        self._status_updates.append(update_status)
        return status_button

    def _update_all_status(self):
        """Update all status indicators in a single timer tick"""
        for update_func in self._status_updates:
            try:
                update_func()
            except Exception as e:
                logger.error(f"Error updating status: {e}", exc_info=True)

    def create_yield_plot(self):
        """Create an interactive yield plot using echarts"""
        # Format timestamps for x-axis
        timestamps = [t.strftime('%H:%M:%S') for t in self.timestamps]
        
        # Determine appropriate scale for bases
        max_bases = max(self.pass_bases) if self.pass_bases else 0
        if max_bases >= 1:  # Greater than 1 Gb
            base_scale = 1
            base_unit = 'Gb'
            formatted_bases = self.pass_bases
        else:
            base_scale = 1000  # Convert to Mb
            base_unit = 'Mb'
            formatted_bases = [b * 1000 for b in self.pass_bases]  # Convert Gb to Mb
        
        options = {
            'animation': False,  # Disable animations for smoother updates
            'title': {
                'text': 'Sequencing Yield Over Time',
                'left': 'center'
            },
            'tooltip': {
                'trigger': 'axis',
                'axisPointer': {
                    'type': 'cross'
                },
                ':formatter': f'''
                    function(params) {{
                        let result = params[0].axisValueLabel + '<br/>';
                        for (let i = 0; i < params.length; i++) {{
                            let value = params[i].value;
                            let marker = params[i].marker;
                            let name = params[i].seriesName;
                            if (name.includes('Bases')) {{
                                value = value.toFixed(1) + ' {base_unit}';
                            }} else {{
                                value = value.toLocaleString();
                            }}
                            result += marker + ' ' + name + ': ' + value + '<br/>';
                        }}
                        return result;
                    }}
                '''
            },
            'legend': {
                'data': ['Total Reads', 'Pass Reads', f'Pass Bases ({base_unit})'],
                'top': '25px'
            },
            'grid': {
                'left': '3%',
                'right': '4%',
                'bottom': '3%',
                'top': '15%',
                'containLabel': True
            },
            'xAxis': {
                'type': 'category',
                'boundaryGap': False,
                'data': timestamps,
                'axisLabel': {
                    'rotate': 45
                }
            },
            'yAxis': [
                {
                    'type': 'value',
                    'name': 'Read Count',
                    'position': 'left',
                    'axisLine': {
                        'show': True,
                        'lineStyle': {
                            'color': '#5470C6'
                        }
                    },
                    'axisLabel': {
                        ':formatter': 'value => value.toLocaleString()'
                    }
                },
                {
                    'type': 'value',
                    'name': f'Bases ({base_unit})',
                    'position': 'right',
                    'axisLine': {
                        'show': True,
                        'lineStyle': {
                            'color': '#91CC75'
                        }
                    },
                    'axisLabel': {
                        ':formatter': 'value => value.toFixed(1)'
                    }
                }
            ],
            'series': [
                {
                    'name': 'Total Reads',
                    'type': 'line',
                    'smooth': True,
                    'data': list(self.read_counts),
                    'itemStyle': {'color': '#5470C6'},
                    'showSymbol': False  # Hide symbols for smoother line
                },
                {
                    'name': 'Pass Reads',
                    'type': 'line',
                    'smooth': True,
                    'data': list(self.pass_reads),
                    'itemStyle': {'color': '#91CC75'},
                    'showSymbol': False  # Hide symbols for smoother line
                },
                {
                    'name': f'Pass Bases ({base_unit})',
                    'type': 'line',
                    'smooth': True,
                    'yAxisIndex': 1,
                    'data': list(formatted_bases),
                    'itemStyle': {'color': '#FAC858'},
                    'showSymbol': False  # Hide symbols for smoother line
                }
            ]
        }
        
        return options

    def update_yield_data(self, read_count, pass_read_count, pass_bases):
        """
        Update the yield data with new values from the sequencing run.
        
        This method updates the internal data structures with new sequencing
        statistics and triggers UI updates if sufficient time has passed since
        the last update.
        
        Args:
            read_count (int): Total number of reads
            pass_read_count (int): Number of reads passing quality filters
            pass_bases (int): Total number of bases in passing reads
            
        Note:
            Updates are rate-limited by self.update_interval to prevent
            overwhelming the UI.
        """
        current_time = datetime.now()
        
        # Only update if enough time has passed and we're not getting historical data
        if time.time() - self.last_update >= self.update_interval:
            # Check if this timestamp would be newer than our last one
            if not self.timestamps or current_time > self.timestamps[-1]:
                # Update data arrays
                self.timestamps.append(current_time)
                self.read_counts.append(read_count)
                self.pass_reads.append(pass_read_count)
                self.pass_bases.append(pass_bases / 1e9)  # Convert to Gb
                self.last_update = time.time()
                
                # Update current values
                self.Read_Count = read_count
                self.Pass_Read_Count = pass_read_count
                self.Pass_Bases = pass_bases
                
                # Update the plot data
                if hasattr(self, 'yield_plot'):
                    # Determine appropriate scale
                    max_bases = max(self.pass_bases)
                    if max_bases >= 1:  # Greater than 1 Gb
                        base_scale = 1
                        base_unit = 'Gb'
                        formatted_bases = self.pass_bases
                    else:
                        base_scale = 1000  # Convert to Mb
                        base_unit = 'Mb'
                        formatted_bases = [b * 1000 for b in self.pass_bases]  # Convert Gb to Mb
                    
                    # Update x-axis data
                    new_timestamps = [t.strftime('%H:%M:%S') for t in self.timestamps]
                    self.yield_plot.options['xAxis']['data'] = new_timestamps
                    
                    # Update y-axis label
                    self.yield_plot.options['yAxis'][1]['name'] = f'Bases ({base_unit})'
                    
                    # Update legend
                    self.yield_plot.options['legend']['data'] = ['Total Reads', 'Pass Reads', f'Pass Bases ({base_unit})']
                    
                    # Update series data
                    self.yield_plot.options['series'][0]['data'] = list(self.read_counts)
                    self.yield_plot.options['series'][1]['data'] = list(self.pass_reads)
                    self.yield_plot.options['series'][2]['data'] = list(formatted_bases)
                    self.yield_plot.options['series'][2]['name'] = f'Pass Bases ({base_unit})'
                    
                    # Trigger the update
                    self.yield_plot.update()

    def start_acquisition_stream(self):
        """Start streaming acquisition output data"""
        logging.info("=== Starting Acquisition Stream ===")
        logging.debug(f"Acquisition Run ID: {self.acquisition_run_id}")
        
        if not self.acquisition_run_id:
            logging.warning("No acquisition run ID available, cannot start stream")
            return

        try:
            # Create a proper DataSelection object for historical data
            from minknow_api.statistics_pb2 import DataSelection
            data_selection = DataSelection()
            data_selection.start = 0  # From the beginning
            data_selection.step = 60  # Get data every minute
            logging.debug(f"Created DataSelection: start={data_selection.start}, step={data_selection.step}")

            logging.debug("Setting up acquisition stream...")
            # Set up the acquisition output stream
            self._acquisition_stream = self.connection.statistics.stream_acquisition_output(
                acquisition_run_id=self.acquisition_run_id,
                data_selection=data_selection
            )
            logging.debug("Acquisition stream created successfully")

            # Start a thread to process the stream
            logging.debug("Starting acquisition processing thread...")
            self._acquisition_thread = threading.Thread(
                target=self._process_acquisition_stream,
                daemon=True
            )
            self._acquisition_thread.start()
            logging.debug("Acquisition processing thread started")
            logging.info("=== Acquisition Stream Setup Complete ===")
        except Exception as e:
            logging.error(f"Error starting acquisition stream: {str(e)}")

    def _process_acquisition_stream(self):
        """Process the acquisition output stream data"""
        try:
            for data in self._acquisition_stream:
                try:
                    if hasattr(data, 'snapshots') and data.snapshots:
                        # The data structure has snapshots nested inside snapshots
                        for outer_snapshot in data.snapshots:
                            if hasattr(outer_snapshot, 'snapshots'):
                                for snapshot in outer_snapshot.snapshots:
                                    try:
                                        # Get the time in seconds from the acquisition start
                                        acquisition_time = getattr(snapshot, 'seconds', None)
                                        if acquisition_time is None:
                                            logging.debug("No valid timestamp found in snapshot, using current time")
                                            current_time = datetime.now()
                                        else:
                                            # Calculate the actual timestamp
                                            current_time = datetime.fromtimestamp(self.start_time.timestamp() + float(acquisition_time))
                                        
                                        # Check if yield_summary exists
                                        if not hasattr(snapshot, 'yield_summary'):
                                            logging.debug("No yield_summary in snapshot, available fields: %s", 
                                                      [field.name for field in snapshot.DESCRIPTOR.fields])
                                            continue
                                        
                                        yield_summary = snapshot.yield_summary
                                        
                                        # Check if yield_summary has the required fields
                                        if not all(hasattr(yield_summary, attr) for attr in ['read_count', 'basecalled_pass_read_count', 'basecalled_pass_bases']):
                                            logging.debug("Missing required fields in yield_summary, available fields: %s",
                                                      [field.name for field in yield_summary.DESCRIPTOR.fields])
                                            continue
                                        
                                        # Update current values
                                        self.Read_Count = yield_summary.read_count
                                        self.Pass_Read_Count = yield_summary.basecalled_pass_read_count
                                        self.Pass_Bases = yield_summary.basecalled_pass_bases
                                        
                                        # Only append if it's a new timestamp
                                        if not self.timestamps or current_time > self.timestamps[-1]:
                                            self.timestamps.append(current_time)
                                            self.read_counts.append(yield_summary.read_count)
                                            self.pass_reads.append(yield_summary.basecalled_pass_read_count)
                                            self.pass_bases.append(yield_summary.basecalled_pass_bases / 1e9)  # Convert to Gb
                                            
                                            logging.debug("Added data point - Time: %s, Reads: %s, Pass: %s, Bases: %.1fGb",
                                                      current_time.strftime('%H:%M:%S'),
                                                      f"{yield_summary.read_count:,}",
                                                      f"{yield_summary.basecalled_pass_read_count:,}",
                                                      yield_summary.basecalled_pass_bases/1e9)
                                            
                                            # Update the plot data
                                            if hasattr(self, 'yield_plot'):
                                                # Update x-axis data
                                                new_timestamps = [t.strftime('%H:%M:%S') for t in self.timestamps]
                                                self.yield_plot.options['xAxis']['data'] = new_timestamps
                                                
                                                # Update series data
                                                self.yield_plot.options['series'][0]['data'] = list(self.read_counts)
                                                self.yield_plot.options['series'][1]['data'] = list(self.pass_reads)
                                                self.yield_plot.options['series'][2]['data'] = list(self.pass_bases)
                                                
                                                # Trigger the update
                                                self.yield_plot.update()
                                    
                                    except Exception as e:
                                        logging.error("Error processing individual snapshot: %s", str(e))
                                        if hasattr(snapshot, 'DESCRIPTOR'):
                                            logging.debug("Available snapshot fields: %s", 
                                                      [field.name for field in snapshot.DESCRIPTOR.fields])
                            else:
                                logging.debug("No nested snapshots found in outer snapshot")
                            
                except Exception as e:
                    logging.error("Error processing acquisition data: %s", str(e))
        except Exception as e:
            logging.error("Error in acquisition stream: %s", str(e))

    def stream_instance_activity(self) -> None:
        """
        This function will stream instance activity from the minknow api.
        We configure a connection and give it a handler.
        This allows us to call cancel on the handler on exit to escape the loop.
        """
        try:
            logging.info("=== Starting Instance Activity Stream ===")
            self._stream_instance_activity = (
                self.connection.instance.stream_instance_activity()
            )
            for info in self._stream_instance_activity:
                try:
                    if info.HasField("device_info"):
                        self.Channel_Count = info.device_info.device_info.max_channel_count
                    if info.HasField("acquisition_run_info"):
                        logging.info("Received acquisition_run_info:")
                        self.basecalling_config_filename = (
                            info.acquisition_run_info.config_summary.basecalling_config_filename
                        )
                        self.start_time = datetime.fromtimestamp(
                            info.acquisition_run_info.start_time.seconds
                            + info.acquisition_run_info.start_time.nanos / 1e9
                        )
                        self.output_folder = (
                            info.acquisition_run_info.config_summary.reads_directory
                        )
                        logging.debug("  Config: %s", self.basecalling_config_filename)
                        logging.debug("  Start Time: %s", self.start_time)
                        logging.debug("  Output Folder: %s", self.output_folder)
                        
                        # Start acquisition stream when we get run info
                        if hasattr(info.acquisition_run_info, 'run_id'):
                            logging.debug("  Found Run ID: %s", info.acquisition_run_info.run_id)
                            self.acquisition_run_id = info.acquisition_run_info.run_id
                            self.start_acquisition_stream()
                        else:
                            logging.warning("  No Run ID found in acquisition_run_info")
                    if info.HasField("basecall_speed"):
                        self.Mean_Basecall_Speed = info.basecall_speed.mean_basecall_speed
                    if info.HasField("n50"):
                        self.N50 = info.n50.n50
                        self.Estimated_N50 = info.n50.estimated_n50
                    if info.HasField("protocol_run_info"):
                        self.running_kit = info.protocol_run_info.meta_info.tags[
                            "kit"
                        ].string_value
                        self.Flowcell_Type = info.protocol_run_info.meta_info.tags[
                            "flow cell"
                        ].string_value
                        if info.protocol_run_info.phase != 0:
                            self.show = True
                        else:
                            self.show = False

                        self.Experiment_Group = (
                            info.protocol_run_info.user_info.protocol_group_id.value
                        )
                        self.Sample_ID = info.protocol_run_info.user_info.sample_id.value
                    if info.HasField("yield_summary"):
                        try:
                            # Update yield data for plotting
                            self.Read_Count = info.yield_summary.read_count
                            self.Pass_Read_Count = info.yield_summary.basecalled_pass_read_count
                            self.Pass_Bases = info.yield_summary.basecalled_pass_bases
                            
                            # Update the plot data
                            self.update_yield_data(
                                self.Read_Count,
                                self.Pass_Read_Count,
                                self.Pass_Bases
                            )
                        except Exception as e:
                            logging.error("Error processing yield summary: %s", str(e))
                except Exception as e:
                    logging.error("Error processing stream info: %s", str(e))
        except Exception as e:
            logging.error("Error in stream activity: %s", str(e))

    async def start_run(
        self,
        position=None,
        reference=None,
        sample_id=None,
        centreID=None,
        flowcell_id=None,
        kit=None,
        basecall_config=None,
        experiment_duration=None,
        bed_file=None,
    ):
        ui.notify(f"Starting Run {sample_id} on {flowcell_id}!", type="positive")
        # ToDo: At every stage we need to confirm that the correct values have been entered.
        # position = "1B"
        kit = kit
        ###Memo to self... basecall config must not include cfg
        basecall_config = basecall_config
        alignment_reference = reference
        bed_file = bed_file
        experiment_duration = experiment_duration
        current_date = datetime.now()
        centreID = centreID
        experiment_group_id = (
            f"{centreID}_{current_date.strftime('%B')}_{current_date.year}"
        )
        # sample_id = "SAMPLE_ID"
        experiment_specs = []
        self.connection_ip = "127.0.0.1"
        self.manager = Manager(host=self.connection_ip)
        # Add all the positions to the list:
        for pos in self.manager.flow_cell_positions():
            if pos.name == position:
                experiment_specs.append(ExperimentSpec(position=pos))
        # Check if the flowcell ID is correct

        if add_protocol_ids(experiment_specs, kit, basecall_config, flowcell_id):

            # Build arguments for starting protocol:
            alignment_args = protocols.AlignmentArgs(
                reference_files=[alignment_reference],
                bed_file=bed_file,
            )

            basecalling_args = protocols.BasecallingArgs(
                config=basecall_config,
                barcoding=None,
                alignment=alignment_args,
            )

            read_until_args = protocols.ReadUntilArgs(
                filter_type="enrich",
                reference_files=[alignment_reference],
                bed_file=bed_file,
                first_channel=None,
                last_channel=None,
            )

            bam_arguments = protocols.OutputArgs(
                reads_per_file=4000,
                batch_duration="1",
            )
            pod5_arguments = protocols.OutputArgs(
                reads_per_file=4000, batch_duration="1"
            )

            # Now start the protocol(s):
            for spec in experiment_specs:
                position_connection = spec.position.connect()

                # Generate stop criteria for use by Run Until
                # The `runtime` is in seconds, while the `experiment_duration` is in hours
                stop_criteria = protocols.CriteriaValues(
                    runtime=int(experiment_duration * 60 * 60)
                )

                run_id = protocols.start_protocol(
                    position_connection,
                    identifier=spec.protocol_id,
                    sample_id=sample_id,
                    experiment_group=experiment_group_id,
                    barcode_info=None,
                    basecalling=basecalling_args,
                    read_until=read_until_args,
                    fastq_arguments=None,
                    fast5_arguments=None,
                    pod5_arguments=pod5_arguments,
                    bam_arguments=bam_arguments,
                    disable_active_channel_selection=False,
                    mux_scan_period=1.5,
                    stop_criteria=stop_criteria,
                    args=[],  # Any extra args passed.
                )

                flow_cell_info = position_connection.device.get_flow_cell_info()

                ui.notify(
                    f"Started protocol:\n    run_id={run_id}\n    position={spec.position.name}\n    flow_cell_id={flow_cell_info.flow_cell_id}\n",
                    multi_line=True,
                    type="positive",
                )
        else:
            ui.notify("Run Start Failed", type="negative")


@ui.page("/", response_timeout=30)
def index_page() -> None:
    initial_ip = "127.0.0.1"
    my_connection = Manager(host=initial_ip)
    with theme.frame(
        "<strong><font color='#000000'>R</font></strong>apid nanop<strong><font color='#000000'>O</font></strong>re <strong><font color='#000000'>B</font></strong>rain intraoperat<strong><font color='#000000'>I</font></strong>ve classificatio<strong><font color='#000000'>N</font></strong>",
        smalltitle="<strong><font color='#000000'>R.O.B.I.N</font></strong>",
    ):
        # my_connection.connect_to_minknow()
        positions = list(my_connection.flow_cell_positions())
        ui.label(f"{positions[0]}")
        display_object = Minknow_Info(positions[0])


def run_class(port: int, reload: bool):
    """
    Helper function to run the app.
    :param port: The port to serve the app on.
    :param reload: Should we reload the app on changes.
    :return:
    """
    ui.add_css("""
        .monitoring-tile {
            transition: all 0.3s ease;
        }
        .monitoring-tile:hover {
            transform: translateY(-2px);
            box-shadow: 0 4px 8px rgba(0,0,0,0.1);
        }
        .status-chip {
            transition: background-color 0.3s ease;
        }
        .device-header {
            border-bottom: 1px solid #e0e0e0;
        }
        @media (max-width: 600px) {
            .monitoring-grid {
                grid-template-columns: repeat(2, 1fr) !important;
            }
        }
    """)
    
    # Register some fonts that we might need later on.
    app.add_static_files("/fonts", str(Path(__file__).parent.parent / "fonts"))

    ui.run(
        port=port,
        reload=reload,
        title="Readfish NiceGUI",
        storage_secret="waynesworld",
    )  # , native=True, fullscreen=False, window_size=(1200, 1000))


def main():
    """
    Main entry point for the MinKNOW monitoring application.
    
    This function initializes and runs the web interface for monitoring
    MinKNOW sequencing runs. It sets up the necessary routes and starts
    the web server.
    
    Returns:
        None
    """
    run_class(port=12398, reload=False)


# Entrypoint for when GUI is launched by the CLI.
# e.g.: python my_app/my_cli.py
if __name__ in {"__main__", "__mp_main__"}:
    """
    Entrypoint for when GUI is launched by the CLI
    :return: None
    """
    if __name__ == "__mp_main__":
        logging.info("GUI launched by auto-reload")

    run_class(port=12398, reload=True)
