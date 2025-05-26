"""
pdf_extractor.py

This module handles extraction of data from PDF reports and stores it for analysis.
"""

import os
import logging
import re
import json
from datetime import datetime
import PyPDF2
from typing import Dict, Optional

logger = logging.getLogger(__name__)


class PDFExtractor:
    """Class to extract and store data from PDF reports."""

    def __init__(self, storage_dir: str):
        """Initialize the PDF extractor.

        Args:
            storage_dir: Directory to store extracted data
        """
        self.storage_dir = storage_dir
        self.data_file = os.path.join(storage_dir, "extracted_data.json")
        self._ensure_storage_dir()

    def _ensure_storage_dir(self):
        """Ensure the storage directory exists."""
        os.makedirs(self.storage_dir, exist_ok=True)

    def extract_from_pdf(self, pdf_path: str) -> Optional[Dict]:
        """Extract data from a PDF report.

        Args:
            pdf_path: Path to the PDF file

        Returns:
            Dictionary containing extracted data or None if extraction fails
        """
        try:
            with open(pdf_path, "rb") as file:
                reader = PyPDF2.PdfReader(file)
                full_text = ""
                for page in reader.pages:
                    full_text += page.extract_text()

                # Extract sample ID
                sample_id = re.search(r"Sample ID:\s+([a-f0-9\-]+)", full_text)
                sample_id = sample_id.group(1) if sample_id else "Not found"

                # Extract classification data using original pattern
                classification = re.findall(
                    r"Methylation (superfamily|family|class|subclass)\s+([A-Za-z \-]+)\s+([0-9.]+)",
                    full_text,
                )
                classification_dict = {
                    label_type: {"label": label.strip(), "score": float(score)}
                    for label_type, label, score in classification
                }

                # Extract MGMT data
                mgmt_sites = re.search(r"Number of sites used = (\d+)", full_text)
                mgmt_avg = re.search(r"Average methylation = ([0-9.]+)%", full_text)
                mgmt_status = re.search(
                    r"Predicted MGMT promoter status= (\w+)", full_text
                )

                mgmt_info = {
                    "sites_used": int(mgmt_sites.group(1)) if mgmt_sites else None,
                    "average_methylation": (
                        float(mgmt_avg.group(1)) if mgmt_avg else None
                    ),
                    "status": mgmt_status.group(1) if mgmt_status else "Not found",
                }

                # Get timestamp from filename or current time
                timestamp = datetime.now()
                filename_match = re.search(r"(\d{8}_\d{6})", os.path.basename(pdf_path))
                if filename_match:
                    try:
                        timestamp = datetime.strptime(
                            filename_match.group(1), "%Y%m%d_%H%M%S"
                        )
                    except ValueError:
                        pass

                # Log the extracted data for debugging
                logging.info(f"Extracted data from {pdf_path}:")
                logging.info(f"Full text length: {len(full_text)} characters")
                logging.info(f"Classification matches: {classification}")
                logging.info(f"Classification dict: {classification_dict}")
                logging.info(f"MGMT info: {mgmt_info}")

                return {
                    "sample_id": sample_id,
                    "timestamp": timestamp.isoformat(),
                    "classification": classification_dict,
                    "mgmt_info": mgmt_info,
                    "source_file": pdf_path,
                }

        except Exception as e:
            logging.error(f"Error extracting data from PDF {pdf_path}: {str(e)}")
            return None

    def save_data(self, data: Dict) -> bool:
        """Save extracted data to storage.

        Args:
            data: Dictionary containing extracted data

        Returns:
            bool: True if save was successful, False otherwise
        """
        try:
            # Load existing data
            existing_data = []
            if os.path.exists(self.data_file):
                with open(self.data_file, "r") as f:
                    existing_data = json.load(f)

            # Check if this report has already been processed
            source_file = data.get("source_file")
            if source_file and any(
                d.get("source_file") == source_file for d in existing_data
            ):
                logger.info(
                    f"Report {source_file} has already been processed, skipping"
                )
                return True

            # Append new data
            existing_data.append(data)

            # Save updated data
            with open(self.data_file, "w") as f:
                json.dump(existing_data, f, indent=2)

            return True

        except Exception as e:
            logger.error(f"Error saving extracted data: {str(e)}")
            return False

    def get_data(self) -> list:
        """Get all stored data.

        Returns:
            list: List of all extracted data dictionaries
        """
        try:
            if not os.path.exists(self.data_file):
                return []

            with open(self.data_file, "r") as f:
                return json.load(f)

        except Exception as e:
            logger.error(f"Error reading data file: {str(e)}")
            return []

    def process_existing_reports(self, reports_dir: str) -> int:
        """Process all existing PDF reports in a directory.

        Args:
            reports_dir: Directory containing PDF reports

        Returns:
            int: Number of reports successfully processed
        """
        try:
            # Get list of all PDF files in the directory and its subdirectories
            pdf_files = []
            for root, _, files in os.walk(reports_dir):
                for file in files:
                    if file.lower().endswith(".pdf"):
                        pdf_files.append(os.path.join(root, file))

            if not pdf_files:
                logger.info(f"No PDF reports found in {reports_dir}")
                return 0

            # Get list of already processed files
            processed_files = set()
            if os.path.exists(self.data_file):
                with open(self.data_file, "r") as f:
                    data = json.load(f)
                    processed_files = {
                        d.get("source_file") for d in data if d.get("source_file")
                    }

            # Process new files
            processed_count = 0
            for pdf_file in pdf_files:
                if pdf_file not in processed_files:
                    logger.info(f"Processing report: {pdf_file}")
                    data = self.extract_from_pdf(pdf_file)
                    if data and self.save_data(data):
                        processed_count += 1
                        logger.info(f"Successfully processed {pdf_file}")
                    else:
                        logger.warning(f"Failed to process {pdf_file}")

            logger.info(f"Processed {processed_count} new reports")
            return processed_count

        except Exception as e:
            logger.error(f"Error processing existing reports: {str(e)}")
            return 0
