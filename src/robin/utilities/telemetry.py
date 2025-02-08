"""
Telemetry module for collecting anonymous usage statistics.

This module collects non-personally identifiable information about how ROBIN is used,
including when it is run, where it was run, and what arguments were used.
No personal data, results, or file contents are ever collected.
"""

import json
import platform
import uuid
import logging
import time
from datetime import datetime
from typing import Dict, Any, Optional
import requests
from pathlib import Path
import os
import geocoder
from nicegui import ui
from packaging import version
from robin import __about__
import psutil

class Telemetry:
    def __init__(self, endpoint_url: Optional[str] = None):
        """
        Initialize the telemetry system.
        
        Args:
            endpoint_url: Optional URL for the telemetry endpoint. If not provided,
                        telemetry will be disabled.
        """
        self.endpoint_url = endpoint_url
        self.session_id = str(uuid.uuid4())
        self.start_time = datetime.utcnow()
        self._location_info = self._get_location_info()
        self._version_info = self._get_version_info()

    def _get_location_info(self) -> Dict[str, str]:
        """
        Get approximate location information using geocoder.
        Only collects country and city information, never stores the IP address.
        
        Returns:
            Dict containing location information (country, city).
        """
        try:
            #print("Getting location information...")
            # Use ipapi.co provider which is free and doesn't require an API key
            g = geocoder.ip('me')
            
            if g.ok:
                # Get the raw response which contains all available fields
                raw = g.raw
                location_info = {
                    'country': raw.get('country', 'Unknown'),
                    'city': raw.get('city', 'Unknown'),
                    'region': raw.get('region', 'Unknown'),
                    'latlng': g.latlng if g.latlng else None
                }
                #print(f"Successfully retrieved location info: {location_info}")
                return location_info
            else:
                #print("Could not determine location")
                raise Exception("Geolocation failed")
                
        except Exception as e:
            #print(f"Error getting location info: {str(e)}")
            return {
                'country': 'Unknown',
                'city': 'Unknown',
                'region': 'Unknown',
                'latlng': None
            }

    def _get_version_info(self) -> Dict[str, Any]:
        """
        Get version information and status compared to the latest release.
        
        Returns:
            Dict containing version information and status.
        """
        try:
            # Get the remote version from GitHub
            response = requests.get('https://raw.githubusercontent.com/LooseLab/ROBIN/main/src/robin/__about__.py')
            response.raise_for_status()
            
            # Extract version from the response text
            remote_version_str = None
            for line in response.text.split('\n'):
                if line.startswith('__version__'):
                    remote_version_str = line.split('=')[1].strip().strip('"').strip("'")
                    break
            
            if not remote_version_str:
                return {
                    'local_version': __about__.__version__,
                    'remote_version': 'unknown',
                    'status': 'unknown'
                }

            local_version = version.parse(__about__.__version__)
            remote_version = version.parse(remote_version_str)

            if local_version == remote_version:
                status = 'current'
            elif local_version < remote_version:
                status = 'outdated'
            else:
                status = 'development'

            return {
                'local_version': str(local_version),
                'remote_version': str(remote_version),
                'status': status
            }
            
        except Exception as e:
            logging.warning(f"Failed to check version: {str(e)}")
            return {
                'local_version': __about__.__version__,
                'remote_version': 'unknown',
                'status': 'unknown'
            }

    def collect_system_info(self) -> Dict[str, Any]:
        """
        Collect anonymous system information.
        
        Returns:
            Dict containing system information.
        """
        return {
            "os": platform.system(),
            "os_version": platform.release(),
            "python_version": platform.python_version(),
            "machine_arch": platform.machine(),
            "processor": platform.processor(),
            "session_id": self.session_id,
            "timestamp_utc": self.start_time.isoformat(),
            "location": self._location_info,
            "version_info": self._version_info
        }

    def collect_run_info(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """
        Collect information about how the program was run.
        
        Args:
            args: Dictionary of command line arguments
            
        Returns:
            Dict containing run information.
        """
        # Filter out any potentially sensitive information
        safe_args = {
            "threads": args.get("threads"),
            "target_panel": args.get("target_panel"),
            "simtime": args.get("simtime"),
            "showerrors": args.get("showerrors"),
            "browse": args.get("browse"),
            "exclude": args.get("exclude"),
            "experiment_duration": args.get("experiment_duration")
        }
        
        return {
            "arguments": safe_args,
            "version": args.get("version", "unknown")
        }

    def send_telemetry(self, run_args: Dict[str, Any]) -> bool:
        """
        Send telemetry data to the collection endpoint.
        
        Args:
            run_args: Dictionary of run arguments
            
        Returns:
            bool: True if telemetry was sent successfully, False otherwise
        """
        logging.info("Sending telemetry data...")
        logging.debug(f"Endpoint URL: {self.endpoint_url}")
        
        # Create the telemetry data in the exact format expected by the API
        inner_data = {
            "id": str(uuid.uuid4()),  # Generate a new UUID for each record
            "message_type": "robin_init",
            "data": {
                "system_info": {
                    "os": platform.system(),
                    "os_version": platform.release(),
                    "python_version": platform.python_version(),
                    "machine_arch": platform.machine(),
                    "processor": platform.processor(),
                    "session_id": self.session_id,
                    "timestamp_utc": self.start_time.isoformat(),
                    "location": {
                        "country": self._location_info.get('country', 'Unknown'),
                        "city": self._location_info.get('city', 'Unknown'),
                        "region": self._location_info.get('region', 'Unknown'),
                        "latlng": [str(coord) for coord in self._location_info.get('latlng', [])] if self._location_info.get('latlng') else None
                    },
                    "version_info": self._version_info
                },
                "run_info": {
                    "arguments": {
                        "threads": run_args.get("threads"),
                        "target_panel": run_args.get("target_panel"),
                        "simtime": run_args.get("simtime"),
                        "showerrors": run_args.get("showerrors"),
                        "browse": run_args.get("browse"),
                        "exclude": run_args.get("exclude", []),
                        "experiment_duration": run_args.get("experiment_duration")
                    },
                    "version": run_args.get("version", "unknown")
                }
            }
        }

        # Wrap in body field for Lambda, converting inner_data to a JSON string
        telemetry_data = {
            "body": json.dumps(inner_data)
        }

        # Log the telemetry data
        logging.debug("\nTelemetry Data Being Collected:")
        logging.debug("--------------------------------")
        logging.info("Telemetry data is collected anonymously to help us understand how ROBIN is used.")
        logging.info("It is not used to identify you or your data.")
        logging.info("It is only used to help us improve ROBIN.")
        logging.info("You can opt out by running with --no-telemetry.")
        logging.debug("The information collected is as follows:")
        logging.debug("--------------------------------")
        logging.debug(json.dumps(inner_data, indent=2))

        if not self.endpoint_url:
            return False

        try:
            # Prepare headers
            headers = {
                'Content-Type': 'application/json'
            }
            
            # Send data to endpoint using POST
            response = requests.post(
                self.endpoint_url,
                headers=headers,
                json=telemetry_data,
                timeout=5  # 5 second timeout
            )
            
            if response.status_code != 200:
                logging.warning(f"Failed to send telemetry. Status code: {response.status_code}, Response: {response.text}")
                return False

            # Parse the response
            response_data = response.json()
            if response_data.get('statusCode') == 400:
                logging.warning(f"Failed to process telemetry. Response: {response.text}")
                return False

            logging.info(f"Telemetry sent successfully. Status code: {response.status_code}")
            logging.debug(f"Response: {response.text}")
            logging.debug(f"Record ID: {inner_data['id']}")
            return True
            
        except Exception as e:
            logging.warning(f"Failed to send telemetry: {str(e)}")
            return False

    def save_local_telemetry(self, run_args: Dict[str, Any], output_dir: Path) -> bool:
        """
        Save telemetry data locally if endpoint is not available.
        
        Args:
            run_args: Dictionary of run arguments
            output_dir: Directory to save telemetry data
            
        Returns:
            bool: True if telemetry was saved successfully, False otherwise
        """
        try:
            telemetry_data = {
                "message_type": "robin_init",
                "system_info": self.collect_system_info(),
                "run_info": self.collect_run_info(run_args)
            }
            
            # Create telemetry directory if it doesn't exist
            telemetry_dir = output_dir / "telemetry"
            telemetry_dir.mkdir(exist_ok=True)
            
            # Save telemetry data with timestamp and session ID
            filename = f"telemetry_{self.start_time.strftime('%Y%m%d_%H%M%S')}_{self.session_id[:8]}.json"
            with open(telemetry_dir / filename, 'w') as f:
                json.dump(telemetry_data, f, indent=2)
                
            return True
            
        except Exception as e:
            logging.warning(f"Failed to save local telemetry: {str(e)}")
            return False

    def get_map_coordinates(self) -> Optional[tuple[float, float]]:
        """
        Get the coordinates for the map display.
        
        Returns:
            Optional tuple of (latitude, longitude) or None if location is unknown
        """
        return self._location_info.get('latlng')

    def get_location_description(self) -> str:
        """
        Get a human-readable description of the location.
        
        Returns:
            String describing the location
        """
        parts = []
        if self._location_info['city'] != 'Unknown':
            parts.append(self._location_info['city'])
        if self._location_info['region'] != 'Unknown':
            parts.append(self._location_info['region'])
        if self._location_info['country'] != 'Unknown':
            parts.append(self._location_info['country'])
        
        return ', '.join(parts) if parts else 'Unknown Location'

    def create_map_element(self) -> None:
        """
        Create a NiceGUI leaflet map element showing the user's location.
        Must be called within a NiceGUI context.
        """
        try:
            coords = self.get_map_coordinates()
            location_desc = self.get_location_description()
            
            logging.info(f"Creating map with coordinates: {coords}")
            logging.info(f"Location description: {location_desc}")
            logging.info(f"Location info: {self._location_info}")
            
            with ui.card().classes('w-192'):  # Double the width from w-96 to w-192
                ui.label('Your Location').classes(
                                    "text-sky-600 dark:text-white"
                                ).style("font-size: 150%; font-weight: 300").tailwind(
                                    "drop-shadow", "font-bold"
                                )#.classes('text-xl font-bold')
                ui.label(location_desc)
                
                if coords and len(coords) == 2 and all(isinstance(x, (int, float)) for x in coords):
                    logging.info("Creating leaflet map with valid coordinates")
                    m = ui.leaflet().classes('w-full h-64')
                    m.set_center(coords)  # Set the center coordinates
                    m.set_zoom(10)  # Set the zoom level
                    # Just create a basic marker at the coordinates
                    m.marker(coords)
                    # Add the location description below the map
                    ui.label(f'ROBIN Instance at {location_desc}').classes('text-sm text-gray-600')
                else:
                    logging.warning(f"Invalid coordinates received: {coords}")
                    ui.label('Location map not available - invalid coordinates')
                    
        except Exception as e:
            logging.error(f"Failed to create map element: {str(e)}", exc_info=True)
            with ui.card().classes('w-192'):  # Also update the error card width
                ui.label('Location').classes('text-xl font-bold')
                ui.label('Error creating location map')

    def send_run_telemetry(self, methnice_instance) -> bool:
        """
        Send run telemetry data to the collection endpoint.
        
        Args:
            methnice_instance: Instance of Methnice class to get current metrics from
            
        Returns:
            bool: True if telemetry was sent successfully, False otherwise
        """
        if not self.endpoint_url:
            return False

        try:
            logging.debug("\n=== Sending Periodic Telemetry Update ===")
            logging.debug(f"Session ID: {self.session_id}")
            logging.debug(f"Time since start: {round((datetime.utcnow() - self.start_time).total_seconds() / 60, 1)} minutes")
            
            # Get hardware metrics
            memory_mb = round(psutil.Process().memory_info().rss / (1024 * 1024), 2)
            cpu_pct = round(psutil.cpu_percent(), 2)
            disk_read = round(psutil.disk_io_counters().read_bytes / (1024 * 1024), 2)
            disk_write = round(psutil.disk_io_counters().write_bytes / (1024 * 1024), 2)
            
            logging.debug("\nSystem Metrics:")
            logging.debug(f"Memory Usage: {memory_mb} MB")
            logging.debug(f"CPU Usage: {cpu_pct}%")
            logging.debug(f"Disk Read: {disk_read} MB")
            logging.debug(f"Disk Write: {disk_write} MB")
            logging.debug(f"Threads: {len(psutil.Process().threads())}")
            
            # Get analysis metrics
            reads_processed = len(methnice_instance.robin.processed_files) if hasattr(methnice_instance.robin, 'processed_files') else 0
            sample_count = len(methnice_instance.robin.samples) if hasattr(methnice_instance.robin, 'samples') else 0
            queue_size = len(methnice_instance.robin.file_queue) if hasattr(methnice_instance.robin, 'file_queue') else 0
            
            logging.debug("\nAnalysis Status:")
            logging.debug(f"Reads Processed: {reads_processed}")
            logging.debug(f"Samples: {sample_count}")
            logging.debug(f"Queue Size: {queue_size}")
            logging.debug(f"Active Modules: {[m for m in ['sturgeon', 'forest', 'nanodx', 'pannanodx', 'cnv', 'fusion', 'coverage', 'mgmt'] if m not in methnice_instance.exclude]}")

            # Create the telemetry data in the format expected by the API
            inner_data = {
                "id": str(uuid.uuid4()),  # Generate a new UUID for each record
                "message_type": "robin_run",
                "data": {
                    "system_info": {
                        "os": platform.system(),
                        "os_version": platform.release(),
                        "python_version": platform.python_version(),
                        "machine_arch": platform.machine(),
                        "processor": platform.processor(),
                        "session_id": self.session_id,
                        "timestamp_utc": datetime.utcnow().isoformat(),
                        "location": {
                            "country": self._location_info.get('country', 'Unknown'),
                            "city": self._location_info.get('city', 'Unknown'),
                            "region": self._location_info.get('region', 'Unknown'),
                            "latlng": [str(coord) for coord in self._location_info.get('latlng', [])] if self._location_info.get('latlng') else None
                        },
                        "version_info": self._version_info
                    },
                    "run_info": {
                        "arguments": {
                            "threads": methnice_instance.threads,
                            "target_panel": methnice_instance.target_panel,
                            "simtime": methnice_instance.simtime,
                            "showerrors": methnice_instance.showerrors,
                            "browse": methnice_instance.browse,
                            "exclude": methnice_instance.exclude,
                            "experiment_duration": methnice_instance.experiment_duration
                        },
                        "version": __about__.__version__,
                        "run_id": str(uuid.uuid4()),
                        "session_id": self.session_id,
                        "timestamp_utc": datetime.utcnow().isoformat()
                    },
                    "hardware_metrics": {
                        "peak_memory_mb": memory_mb,
                        "cpu_percent": cpu_pct,
                        "disk_read_mb": disk_read,
                        "disk_write_mb": disk_write,
                        "thread_count": len(psutil.Process().threads())
                    },
                    "analysis_metrics": {
                        "target_panel": methnice_instance.target_panel,
                        "reads_processed": reads_processed,
                        "duration_seconds": round((datetime.utcnow() - self.start_time).total_seconds(), 2),
                        "active_modules": [m for m in ["sturgeon", "forest", "nanodx", "pannanodx", "cnv", "fusion", "coverage", "mgmt"] 
                                        if m not in methnice_instance.exclude],
                        "success": True,
                        "error_type": None
                    },
                    "configuration": {
                        "is_live": not methnice_instance.browse,
                        "adaptive_sampling": methnice_instance.readfish_toml is not None,
                        "basecalling_config": {
                            "model": methnice_instance.basecall_config
                        } if methnice_instance.basecall_config else None,
                        "experiment_duration": methnice_instance.experiment_duration
                    },
                    "performance": {
                        "processing_rate": 0.0,  # TODO: Calculate from methnice instance
                        "sample_count": sample_count,
                        "queue_size": queue_size
                    }
                }
            }

            logging.debug("\nComplete Telemetry Data:")
            logging.debug("=" * 40)
            logging.debug(json.dumps(inner_data, indent=2))
            logging.debug("=" * 40)

            # Wrap in body field for Lambda
            telemetry_data = {
                "body": json.dumps(inner_data)
            }

            # Send data to endpoint using POST
            response = requests.post(
                self.endpoint_url,
                headers={'Content-Type': 'application/json'},
                json=telemetry_data,
                timeout=5  # 5 second timeout
            )
            
            if response.status_code != 200:
                logging.warning("\nFailed to send telemetry!")
                logging.warning(f"Status code: {response.status_code}")
                logging.warning(f"Response: {response.text}")
                return False

            # Parse the response
            response_data = response.json()
            if response_data.get('statusCode') == 400:
                logging.warning("\nServer rejected telemetry data!")
                logging.warning(f"Error: {response_data.get('body')}")
                return False

            logging.info("Telemetry update sent successfully!")
            logging.debug("=" * 40 + "\n")
            return True
            
        except Exception as e:
            logging.warning(f"Error sending telemetry: {str(e)}")
            return False