#!/usr/bin/env python3
"""
NanoDX and PanNanoDX analysis module for robin.

This module provides neural network-based classification for methylation data
using NanoDX and PanNanoDX models from the robin package.
"""

import gc
import json
import logging
import os
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Tuple

import numpy as np
import pandas as pd

# Import robin utilities
try:
    # ToDo: Resolve models into robin
    from robin import models, resources
    from robin.analysis.utilities.merge_bedmethyl import (
        collapse_minimal_bedmethyl,
        load_modkit_data,
    )
    from robin.submodules.nanoDX.workflow.scripts.NN_model import NN_classifier
except ImportError as e:
    logging.warning(f"Some robin dependencies not available: {e}")
    collapse_minimal_bedmethyl = None
    NN_classifier = None
    resources = None
    models = None
    load_modkit_data = None

from robin.logging_config import get_job_logger


def _is_fail_only_expected(job) -> bool:
    """True when errors are expected for fail-only BAM submissions."""
    try:
        if not bool(job.context.metadata.get("fail_only_bam_submission", False)):
            return False
        fp = getattr(job.context, "filepath", "") or ""
        base = os.path.basename(fp).lower()
        return ("fail" in base) and ("pass" not in base)
    except Exception:
        return False


@dataclass
class NanodxMetadata:
    """Container for NanoDX analysis metadata and results"""

    sample_id: str
    parquet_path: str
    analysis_timestamp: float
    model_used: str
    predictions: Optional[np.ndarray] = None
    class_labels: Optional[np.ndarray] = None
    n_features: int = 0
    processing_steps: list = field(default_factory=list)
    error_message: Optional[str] = None
    results: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        if self.processing_steps is None:
            self.processing_steps = []
        if self.results is None:
            self.results = {}


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
    logger = logging.getLogger("robin.nanodx")
    logger.debug(f"Running classification with model file: {modelfile}")

    # Save DataFrame as a temporary CSV file
    with tempfile.NamedTemporaryFile(delete=True, suffix=".csv") as tmp_input:
        input_file = tmp_input.name
        test_df.to_csv(input_file, sep=",", index=False, encoding="utf-8")

        # Create a temporary output file for the prediction results
        with tempfile.NamedTemporaryFile(delete=True, suffix=".json") as tmp_output:
            output_file = tmp_output.name

            # Run classification in separate process
            try:
                # Get the path to this script
                current_script = os.path.abspath(__file__)

                # Run the classification in a separate process
                result = subprocess.run(
                    [
                        sys.executable,  # Use the same Python interpreter
                        current_script,
                        "classify",  # Command to run classification
                        "--model",
                        modelfile,
                        "--input",
                        input_file,
                        "--output",
                        output_file,
                    ],
                    capture_output=True,
                    text=True,
                    timeout=300,  # 5 minute timeout
                )

                if result.returncode != 0:
                    logger.error(f"Classification failed: {result.stderr}")
                    # Fallback to direct classification if subprocess fails
                    logger.info("Falling back to direct classification")
                    return _direct_classification(modelfile, test_df)

                # Read the classification results
                with open(output_file, "r") as f:
                    result_data = json.load(f)

                # Extract results
                predictions = np.array(result_data["predictions"])
                class_labels = np.array(result_data["class_labels"])
                n_features = result_data["n_features"]

                logger.debug("Classification executed successfully via subprocess.")

            except subprocess.TimeoutExpired:
                logger.error("Classification timed out after 5 minutes")
                # Fallback to direct classification
                logger.info("Falling back to direct classification")
                return _direct_classification(modelfile, test_df)
            except Exception as e:
                logger.error(f"Subprocess classification failed: {str(e)}")
                # Fallback to direct classification
                logger.info("Falling back to direct classification")
                return _direct_classification(modelfile, test_df)

    return predictions, class_labels, n_features


def _direct_classification(
    modelfile: str, test_df: pd.DataFrame
) -> Tuple[np.ndarray, np.ndarray, int]:
    """
    Direct classification function as fallback when subprocess fails.

    Args:
        modelfile (str): Path to the neural network model file.
        test_df (pd.DataFrame): DataFrame containing the test data.

    Returns:
        Tuple[np.ndarray, np.ndarray, int]: Predictions, class labels, and number of features.
    """
    logger = logging.getLogger("robin.nanodx")
    logger.debug(f"Running direct classification with model file: {modelfile}")

    if NN_classifier is None:
        logger.error("NN_classifier not available")
        return np.array([]), np.array([]), 0

    NN = NN_classifier(modelfile)
    try:
        predictions, class_labels, n_features = NN.predict(test_df)
        logger.debug("Direct classification executed successfully.")
    except Exception:
        logger.error("An error occurred during direct classification", exc_info=True)
        test_df.to_csv("errordf.csv", sep=",", index=False, encoding="utf-8")
        # Return empty results on error
        return np.array([]), np.array([]), 0
    finally:
        NN = None
        del NN
        gc.collect()

    # Convert lists to numpy arrays to match type hints
    if isinstance(predictions, list):
        predictions = np.array(predictions)
    if isinstance(class_labels, list):
        class_labels = np.array(class_labels)

    return predictions, class_labels, n_features


class NanodxAnalysis:
    """NanoDX analysis worker"""

    def __init__(self, work_dir=None, model: str = "Capper_et_al_NN_v2.pkl"):
        self.work_dir = work_dir or os.getcwd()
        self.model = model

        # Find required files and paths
        self.cpgs_file = self._find_cpgs_file()
        self.modelfile = self._find_model_file()

        # Load CPGs data
        self.cpgs = self._load_cpgs_data()

        # Determine store file name
        if self.model != "Capper_et_al_NN_v2.pkl":
            self.storefile = "PanNanoDX_scores.csv"
        else:
            self.storefile = "NanoDX_scores.csv"

        # Use logger for initialization
        logger = logging.getLogger("robin.nanodx")
        logger.info(f"NanoDX Analysis initialized with model: {self.model}")
        logger.debug(f"Work directory: {self.work_dir}")
        logger.debug(f"CPGs file: {self.cpgs_file}")
        logger.debug(f"Model file: {self.modelfile}")
        logger.debug(f"Store file: {self.storefile}")

    def _find_cpgs_file(self) -> str:
        """Find the CPGs file from robin resources"""
        logger = logging.getLogger("robin.nanodx")

        if resources is not None:
            try:
                cpgs_path = os.path.join(
                    os.path.dirname(os.path.abspath(resources.__file__)),
                    "hglft_genome_260e9_91a970_clean.bed",
                )
                if os.path.exists(cpgs_path):
                    return cpgs_path
            except Exception:
                pass

        # Fallback paths
        possible_paths = [
            "hglft_genome_260e9_91a970_clean.bed",
            "data/hglft_genome_260e9_91a970_clean.bed",
            "/usr/local/share/hglft_genome_260e9_91a970_clean.bed",
        ]

        for path in possible_paths:
            if os.path.exists(path):
                return path

        # If not found, create a placeholder (this will cause an error later)
        logger.warning("CPGs file not found, will use placeholder")
        return "hglft_genome_260e9_91a970_clean.bed"

    def _find_model_file(self) -> str:
        """Find the model file from robin models"""
        logger = logging.getLogger("robin.nanodx")

        if models is not None:
            try:
                model_path = str(models.DIR / self.model)
                if os.path.exists(model_path):
                    return model_path
            except Exception:
                pass

        # Fallback paths
        possible_paths = [
            self.model,
            f"models/{self.model}",
            f"data/{self.model}",
            f"/usr/local/share/{self.model}",
        ]

        for path in possible_paths:
            if os.path.exists(path):
                return path

        # If not found, create a placeholder (this will cause an error later)
        logger.warning(f"Model file {self.model} not found, will use placeholder")
        return self.model

    def _load_cpgs_data(self) -> pd.DataFrame:
        """Load CPGs data from file"""
        logger = logging.getLogger("robin.nanodx")

        try:
            cpgs = pd.read_csv(
                self.cpgs_file,
                sep="\t",
                header=None,
            )
            logger.debug(f"Loaded CPGs data with shape: {cpgs.shape}")
            return cpgs
        except Exception as e:
            logger.error(f"Failed to load CPGs data: {e}")
            # Return empty DataFrame as fallback
            return pd.DataFrame()

    def process_parquet_file(self, parquet_path: str, sample_id: str) -> NanodxMetadata:
        """Process a single parquet file for NanoDX analysis"""
        logger = logging.getLogger("robin.nanodx")

        logger.info(f"Processing parquet file: {parquet_path} for NanoDX analysis")
        start_time = time.time()

        nanodx_result = NanodxMetadata(
            sample_id=sample_id,
            parquet_path=parquet_path,
            analysis_timestamp=start_time,
            model_used=self.model,
        )

        logger.info(f"Starting NanoDX analysis for sample: {sample_id}")
        logger.debug(f"Parquet file: {os.path.basename(parquet_path)}")
        logger.debug(f"Model: {self.model}")

        try:
            # Step 1: Check if parquet file exists and is readable
            if not os.path.exists(parquet_path):
                raise FileNotFoundError(f"Parquet file not found: {parquet_path}")

            nanodx_result.processing_steps.append("file_validation")

            # Step 2: Load modkit data from parquet
            if load_modkit_data is None:
                raise ImportError("load_modkit_data function not available")

            logger.info("Loading modkit data from parquet...")
            merged_modkit_df = load_modkit_data(parquet_path)
            nanodx_result.processing_steps.append("modkit_data_loaded")

            # Step 3: Collapse minimal bedmethyl data
            if collapse_minimal_bedmethyl is None:
                raise ImportError("collapse_minimal_bedmethyl function not available")

            logger.info("Collapsing minimal bedmethyl data...")
            nanodx_df = collapse_minimal_bedmethyl(merged_modkit_df)
            nanodx_result.processing_steps.append("bedmethyl_collapsed")

            # Step 4: Merge with CPGs data
            logger.info("Merging with CPGs data...")
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
            nanodx_result.processing_steps.append("cpgs_merged")

            # Step 5: Apply methylation call thresholds
            logger.info("Applying methylation call thresholds...")
            test_df.loc[test_df["methylation_call"] < 60, "methylation_call"] = -1
            test_df.loc[test_df["methylation_call"] >= 60, "methylation_call"] = 1
            nanodx_result.processing_steps.append("thresholds_applied")

            # Step 6: Run classification
            logger.info("Running neural network classification...")
            predictions, class_labels, n_features = _direct_classification(
                self.modelfile, test_df
            )

            nanodx_result.predictions = predictions
            nanodx_result.class_labels = class_labels
            nanodx_result.n_features = n_features
            nanodx_result.processing_steps.append("classification_complete")

            # Step 7: Create results DataFrame
            logger.info("Creating results DataFrame...")
            nanoDX_df = pd.DataFrame({"class": class_labels, "score": predictions})
            nanoDX_save = nanoDX_df.set_index("class").T
            nanoDX_save["number_probes"] = n_features
            nanoDX_save["timestamp"] = start_time * 1000  # Convert to milliseconds

            # Step 8: Save results
            logger.info("Saving results...")
            sample_dir = os.path.join(self.work_dir, sample_id)
            os.makedirs(sample_dir, exist_ok=True)

            store_path = os.path.join(sample_dir, self.storefile)

            if os.path.exists(store_path):
                nanodx_df_store = pd.read_csv(store_path, index_col=0)
            else:
                nanodx_df_store = pd.DataFrame()

            nanodx_df_store = pd.concat(
                [
                    nanodx_df_store,
                    nanoDX_save.set_index("timestamp"),
                ]
            )

            nanodx_df_store.to_csv(store_path)
            nanodx_result.processing_steps.append("results_saved")

            # Store results in metadata
            nanodx_result.results = {
                "predictions": predictions.tolist() if len(predictions) > 0 else [],
                "class_labels": class_labels.tolist() if len(class_labels) > 0 else [],
                "n_features": n_features,
                "store_path": store_path,
                "model_used": self.model,
            }

            logger.info(f"NanoDX analysis complete for {sample_id}")
            logger.info(f"   Results saved to: {store_path}")
            logger.info(f"   Number of features: {n_features}")
            logger.info(f"   Model used: {self.model}")

            return nanodx_result

        except Exception as e:
            logger.error(f"Error in NanoDX analysis for {sample_id}: {e}")
            nanodx_result.error_message = str(e)
            nanodx_result.processing_steps.append("analysis_failed")
            return nanodx_result


class PanNanodxAnalysis(NanodxAnalysis):
    """PanNanoDX analysis worker"""

    def __init__(self, work_dir=None, model: str = "pancan_devel_v5i_NN_v2.pkl"):
        super().__init__(work_dir=work_dir, model=model)

        # Override store file name for PanNanoDX
        self.storefile = "PanNanoDX_scores.csv"

        # Use logger for initialization
        logger = logging.getLogger("robin.pannanodx")
        logger.info(f"PanNanoDX Analysis initialized with model: {self.model}")
        logger.debug(f"Work directory: {self.work_dir}")
        logger.debug(f"CPGs file: {self.cpgs_file}")
        logger.debug(f"Model file: {self.modelfile}")
        logger.debug(f"Store file: {self.storefile}")


def process_multiple_files(
    parquet_paths, metadata_list, work_dir, logger, model="Capper_et_al_NN_v2.pkl"
):
    """
    Process multiple parquet files for NanoDX analysis.

    This function processes multiple parquet files for the same sample.
    Each file is processed individually and results are accumulated in
    the same NanoDX_scores.csv or PanNanoDX_scores.csv file.

    Args:
        parquet_paths: List of paths to parquet files
        metadata_list: List of metadata dictionaries (one per parquet file)
        work_dir: Working directory
        logger: Logger instance
        model: Model to use ("Capper_et_al_NN_v2.pkl" for NanoDX or "pancan_devel_v5i_NN_v2.pkl" for PanNanoDX)

    Returns:
        Dictionary with aggregated NanoDX analysis results
    """
    if not parquet_paths or not metadata_list:
        raise ValueError("parquet_paths and metadata_list must not be empty")

    if len(parquet_paths) != len(metadata_list):
        raise ValueError("parquet_paths and metadata_list must have the same length")

    # Get sample ID from first metadata (assuming all parquet files are from same sample)
    sample_id = metadata_list[0].get("sample_id", "unknown")

    # Determine analysis type based on model
    is_pannanodx = model != "Capper_et_al_NN_v2.pkl"
    analysis_type = "PanNanoDX" if is_pannanodx else "NanoDX"

    logger.info(
        f"🧠 Starting multi-file {analysis_type} analysis for sample: {sample_id}"
    )
    logger.info(f"Processing {len(parquet_paths)} parquet files for sample {sample_id}")
    logger.info(f"Using model: {model}")

    # Log essential metadata only
    for i, (parquet_path, metadata) in enumerate(zip(parquet_paths, metadata_list)):
        logger.debug(f"Parquet file {i+1}: {os.path.basename(parquet_path)}")

    analysis_result = {
        "sample_id": sample_id,
        "parquet_paths": parquet_paths,
        "analysis_timestamp": time.time(),
        "processing_steps": [],
        "error_message": None,
        "files_processed": 0,
        "total_files": len(parquet_paths),
        "model_used": model,
        "analysis_type": analysis_type,
        "store_path": None,
        "total_features": 0,
    }

    try:
        # Create sample-specific output directory
        sample_dir = os.path.join(work_dir, sample_id)
        os.makedirs(sample_dir, exist_ok=True)

        # Determine store file name
        if is_pannanodx:
            storefile = "PanNanoDX_scores.csv"
        else:
            storefile = "NanoDX_scores.csv"

        store_path = os.path.join(sample_dir, storefile)
        analysis_result["store_path"] = store_path

        logger.info(f"Created output directory: {sample_dir}")
        logger.info(f"Store file: {storefile}")
        analysis_result["processing_steps"].append("directory_created")

        # Initialize NanoDX analysis
        if is_pannanodx:
            nanodx_analyzer = PanNanodxAnalysis(work_dir=work_dir, model=model)
        else:
            nanodx_analyzer = NanodxAnalysis(work_dir=work_dir, model=model)

        logger.info(f"Initialized {analysis_type} analyzer")
        analysis_result["processing_steps"].append("analyzer_initialized")

        # Process each parquet file individually
        logger.info("Processing parquet files individually")
        processed_files = 0
        total_features = 0

        for i, (parquet_path, metadata) in enumerate(zip(parquet_paths, metadata_list)):
            logger.info(
                f"Processing parquet file {i+1}/{len(parquet_paths)}: {os.path.basename(parquet_path)}"
            )

            try:
                # Check if parquet file exists
                if not os.path.exists(parquet_path):
                    logger.warning(f"Parquet file not found: {parquet_path}")
                    continue

                # Process the parquet file
                nanodx_result = nanodx_analyzer.process_parquet_file(
                    parquet_path, sample_id
                )

                if nanodx_result.error_message:
                    logger.warning(
                        f"Error processing {os.path.basename(parquet_path)}: {nanodx_result.error_message}"
                    )
                    continue

                processed_files += 1
                total_features += nanodx_result.n_features
                logger.debug(
                    f"Successfully processed file {i+1}: {os.path.basename(parquet_path)}"
                )
                logger.debug(f"Features from this file: {nanodx_result.n_features}")

            except Exception as e:
                logger.warning(
                    f"Error processing {os.path.basename(parquet_path)}: {e}"
                )
                continue

        if processed_files == 0:
            analysis_result["error_message"] = (
                "No files could be processed successfully"
            )
            analysis_result["processing_steps"].append("no_files_processed")
            return analysis_result

        analysis_result["files_processed"] = processed_files
        analysis_result["total_features"] = total_features
        analysis_result["processing_steps"].append("files_processed")

        # Check if scores file was created
        if os.path.exists(store_path):
            analysis_result["processing_steps"].append("scores_file_created")
            logger.info(f"{analysis_type} scores accumulated in: {store_path}")

            # Load and log summary of results
            try:
                scores_df = pd.read_csv(store_path, index_col=0)
                logger.info(f"Scores file contains {len(scores_df)} entries")
                logger.debug(f"Scores file shape: {scores_df.shape}")
            except Exception as e:
                logger.warning(f"Could not read scores file for summary: {e}")
        else:
            logger.warning(f"No {storefile} file was created")

        analysis_result["processing_steps"].append("analysis_complete")
        logger.info(f"Multi-file {analysis_type} analysis completed for {sample_id}")
        logger.info(
            f"Files successfully processed: {analysis_result['files_processed']}/{analysis_result['total_files']}"
        )
        logger.info(f"Total features processed: {analysis_result['total_features']}")
        logger.info(f"Output directory: {sample_dir}")

        return analysis_result

    except Exception as e:
        logger.error(
            f"Error in multi-file {analysis_type} analysis for {sample_id}: {e}"
        )
        analysis_result["error_message"] = str(e)
        analysis_result["processing_steps"].append("analysis_failed")
        return analysis_result


def nanodx_handler(job, work_dir=None):
    """
    Handler function for NanoDX analysis jobs.
    This function processes parquet files for NanoDX classification.

    Args:
        job: The workflow job containing file and metadata
        work_dir: Optional base directory for output
    """
    # Get job-specific logger first, before any potential exceptions
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
    suppress_expected = _is_fail_only_expected(job)

    try:
        # Check if this is a batched job
        batched_job = job.context.metadata.get("_batched_job")
        if batched_job:
            batch_size = batched_job.get_file_count()
            sample_id = batched_job.get_sample_id()
            batch_id = batched_job.batch_id
            logger.info(
                f"Processing NanoDX analysis batch: {batch_size} files for sample '{sample_id}' (batch_id: {batch_id})"
            )

            # Get all filepaths in the batch
            filepaths = batched_job.get_filepaths()

            # Log individual files in the batch
            for i, filepath in enumerate(filepaths):
                logger.info(
                    f"  Batch file {i+1}/{batch_size}: {os.path.basename(filepath)}"
                )

            # Prepare metadata list and extract parquet paths from bed conversion results
            metadata_list = []
            parquet_paths = []

            for i, bam_path in enumerate(filepaths):
                # Get metadata from preprocessing for this specific file
                file_metadata = batched_job.contexts[i].metadata.get("bam_metadata", {})

                # Get sample ID from preprocessing results for this specific file
                file_context = batched_job.contexts[i]
                file_sample_id = file_context.get_sample_id()

                # Use the sample ID from the file's context (which should have preprocessing results)
                if file_sample_id != "unknown":
                    file_metadata["sample_id"] = file_sample_id
                else:
                    file_metadata["sample_id"] = sample_id

                # Get parquet path from bed conversion results for this file
                bed_conversion_result = batched_job.contexts[i].results.get(
                    "bed_conversion", {}
                )
                parquet_path = bed_conversion_result.get("parquet_path")

                if parquet_path:
                    parquet_paths.append(parquet_path)
                    metadata_list.append(file_metadata)
                    logger.debug(
                        f"Found parquet path for file {i+1}: {os.path.basename(parquet_path)}"
                    )
                else:
                    logger.warning(
                        f"No parquet path found for file {i+1}: {os.path.basename(bam_path)}"
                    )

            if not parquet_paths:
                error_msg = (
                    "No parquet paths found from bed conversion results in batch"
                )
                if suppress_expected:
                    logger.warning(
                        f"{error_msg} (expected for fail-only BAM submission)"
                    )
                    job.context.add_result(
                        "nanodx_analysis",
                        {"status": "expected_failure", "reason": error_msg},
                    )
                else:
                    logger.error(error_msg)
                    job.context.add_error("nanodx_analysis", error_msg)
                return

            # Determine work directory for the batch
            if work_dir is None:
                # Default to first parquet file directory
                batch_work_dir = os.path.dirname(parquet_paths[0])
            else:
                # Use specified work directory, create if it doesn't exist
                os.makedirs(work_dir, exist_ok=True)
                batch_work_dir = work_dir
                logger.debug(f"Using specified work directory: {batch_work_dir}")

            # Process all parquet files in the batch using the new aggregated function
            logger.info(
                f"Processing {len(parquet_paths)} parquet files as aggregated batch for sample '{sample_id}'"
            )
            batch_result = process_multiple_files(
                parquet_paths=parquet_paths,
                metadata_list=metadata_list,
                work_dir=batch_work_dir,
                logger=logger,
                model="Capper_et_al_NN_v2.pkl",  # NanoDX model
            )

            # Store batch results in job context (maintain compatibility with existing structure)
            job.context.add_metadata(
                "nanodx_analysis",
                {
                    "batch_result": batch_result,  # Single aggregated result
                    "batch_size": batch_size,
                    "sample_id": sample_id,
                    "batch_id": batch_id,
                    "files_processed": batch_result.get(
                        "files_processed", len(parquet_paths)
                    ),
                    "total_files": batch_result.get("total_files", len(parquet_paths)),
                },
            )

            logger.info(
                f"Completed NanoDX analysis batch processing: {batch_size} files for sample '{sample_id}'"
            )
            logger.info(
                f"Files successfully processed: {batch_result.get('files_processed', len(parquet_paths))}/{batch_result.get('total_files', len(parquet_paths))}"
            )
            logger.info(
                f"Total features processed: {batch_result.get('total_features', 0)}"
            )

            if batch_result.get("error_message"):
                if suppress_expected:
                    logger.warning(
                        f"Batch processing completed with expected errors (fail-only BAM submission): {batch_result['error_message']}"
                    )
                    job.context.add_result(
                        "nanodx_analysis",
                        {
                            "status": "expected_failure",
                            "error_message": batch_result["error_message"],
                            "sample_id": sample_id,
                        },
                    )
                else:
                    logger.error(
                        f"Batch processing completed with errors: {batch_result['error_message']}"
                    )
                    job.context.add_error(
                        "nanodx_analysis", batch_result["error_message"]
                    )
            else:
                logger.info(
                    "Batch processing completed successfully with aggregated NanoDX analysis"
                )
                job.context.add_result(
                    "nanodx_analysis",
                    {
                        "status": "success",
                        "sample_id": sample_id,
                        "analysis_time": batch_result.get("analysis_timestamp", 0),
                        "model_used": batch_result.get(
                            "model_used", "Capper_et_al_NN_v2.pkl"
                        ),
                        "n_features": batch_result.get("total_features", 0),
                        "store_path": batch_result.get("store_path", ""),
                        "processing_steps": batch_result.get("processing_steps", []),
                        "files_processed": batch_result.get(
                            "files_processed", len(parquet_paths)
                        ),
                        "total_files": batch_result.get(
                            "total_files", len(parquet_paths)
                        ),
                    },
                )

            return

        else:
            # Single file processing (backward compatibility)
            # Get metadata from previous steps
            bed_conversion_result = job.context.results.get("bed_conversion", {})
            bam_metadata = job.context.metadata.get("bam_metadata", {})

            # Extract parquet path from bed conversion results
            parquet_path = bed_conversion_result.get("parquet_path")
            if not parquet_path:
                raise ValueError("No parquet path found from bed conversion step")

            logger.info(
                f"Starting NanoDX analysis for: {os.path.basename(parquet_path)}"
            )

            # Get sample ID
            sample_id = bam_metadata.get("sample_id", "unknown")

            # Determine work directory
            if work_dir is None:
                work_dir = os.path.dirname(parquet_path)
            else:
                os.makedirs(work_dir, exist_ok=True)
                logger.debug(f"Using specified work directory: {work_dir}")

            # Create NanoDX analysis instance
            nanodx_analyzer = NanodxAnalysis(work_dir=work_dir)

            # Process the parquet file
            nanodx_result = nanodx_analyzer.process_parquet_file(
                parquet_path, sample_id
            )

            # Store results in job context
            job.context.add_metadata("nanodx_analysis", nanodx_result.results)
            job.context.add_metadata(
                "nanodx_processing_steps", nanodx_result.processing_steps
            )

            if nanodx_result.error_message:
                job.context.add_error("nanodx_analysis", nanodx_result.error_message)
                logger.error(f"NanoDX analysis failed: {nanodx_result.error_message}")
            else:
                job.context.add_result(
                    "nanodx_analysis",
                    {
                        "status": "success",
                        "sample_id": nanodx_result.sample_id,
                        "analysis_time": nanodx_result.analysis_timestamp,
                        "model_used": nanodx_result.model_used,
                        "n_features": nanodx_result.n_features,
                        "store_path": nanodx_result.results.get("store_path"),
                        "processing_steps": nanodx_result.processing_steps,
                    },
                )
                logger.info(
                    f"NanoDX analysis complete for {os.path.basename(parquet_path)}"
                )
                logger.info(f"Sample ID: {nanodx_result.sample_id}")
                logger.info(f"Model used: {nanodx_result.model_used}")
                logger.info(f"Number of features: {nanodx_result.n_features}")
                logger.debug(
                    f"Processing steps: {', '.join(nanodx_result.processing_steps)}"
                )

    except Exception as e:
        if suppress_expected:
            logger.warning(f"Expected NanoDX failure for fail-only BAM submission: {e}")
            job.context.add_result(
                "nanodx_analysis",
                {"status": "expected_failure", "error_message": str(e)},
            )
            return
        job.context.add_error("nanodx_analysis", str(e))
        logger.error(f"Error in NanoDX analysis for {job.context.filepath}: {e}")


def pannanodx_handler(job, work_dir=None):
    """
    Handler function for PanNanoDX analysis jobs.
    This function processes parquet files for PanNanoDX classification.

    Args:
        job: The workflow job containing file and metadata
        work_dir: Optional base directory for output
    """
    # Get job-specific logger first, before any potential exceptions
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
    suppress_expected = _is_fail_only_expected(job)

    try:
        # Check if this is a batched job
        batched_job = job.context.metadata.get("_batched_job")
        if batched_job:
            batch_size = batched_job.get_file_count()
            sample_id = batched_job.get_sample_id()
            batch_id = batched_job.batch_id
            logger.info(
                f"Processing PanNanoDX analysis batch: {batch_size} files for sample '{sample_id}' (batch_id: {batch_id})"
            )

            # Get all filepaths in the batch
            filepaths = batched_job.get_filepaths()

            # Log individual files in the batch
            for i, filepath in enumerate(filepaths):
                logger.info(
                    f"  Batch file {i+1}/{batch_size}: {os.path.basename(filepath)}"
                )

            # Prepare metadata list and extract parquet paths from bed conversion results
            metadata_list = []
            parquet_paths = []

            for i, bam_path in enumerate(filepaths):
                # Get metadata from preprocessing for this specific file
                file_metadata = batched_job.contexts[i].metadata.get("bam_metadata", {})

                # Get sample ID from preprocessing results for this specific file
                file_context = batched_job.contexts[i]
                file_sample_id = file_context.get_sample_id()

                # Use the sample ID from the file's context (which should have preprocessing results)
                if file_sample_id != "unknown":
                    file_metadata["sample_id"] = file_sample_id
                else:
                    file_metadata["sample_id"] = sample_id

                # Get parquet path from bed conversion results for this file
                bed_conversion_result = batched_job.contexts[i].results.get(
                    "bed_conversion", {}
                )
                parquet_path = bed_conversion_result.get("parquet_path")

                if parquet_path:
                    parquet_paths.append(parquet_path)
                    metadata_list.append(file_metadata)
                    logger.debug(
                        f"Found parquet path for file {i+1}: {os.path.basename(parquet_path)}"
                    )
                else:
                    logger.warning(
                        f"No parquet path found for file {i+1}: {os.path.basename(bam_path)}"
                    )

            if not parquet_paths:
                error_msg = (
                    "No parquet paths found from bed conversion results in batch"
                )
                if suppress_expected:
                    logger.warning(
                        f"{error_msg} (expected for fail-only BAM submission)"
                    )
                    job.context.add_result(
                        "pannanodx_analysis",
                        {"status": "expected_failure", "reason": error_msg},
                    )
                else:
                    logger.error(error_msg)
                    job.context.add_error("pannanodx_analysis", error_msg)
                return

            # Determine work directory for the batch
            if work_dir is None:
                # Default to first parquet file directory
                batch_work_dir = os.path.dirname(parquet_paths[0])
            else:
                # Use specified work directory, create if it doesn't exist
                os.makedirs(work_dir, exist_ok=True)
                batch_work_dir = work_dir
                logger.debug(f"Using specified work directory: {batch_work_dir}")

            # Process all parquet files in the batch using the new aggregated function
            logger.info(
                f"Processing {len(parquet_paths)} parquet files as aggregated batch for sample '{sample_id}'"
            )
            batch_result = process_multiple_files(
                parquet_paths=parquet_paths,
                metadata_list=metadata_list,
                work_dir=batch_work_dir,
                logger=logger,
                model="pancan_devel_v5i_NN_v2.pkl",  # PanNanoDX model
            )

            # Store batch results in job context (maintain compatibility with existing structure)
            job.context.add_metadata(
                "pannanodx_analysis",
                {
                    "batch_result": batch_result,  # Single aggregated result
                    "batch_size": batch_size,
                    "sample_id": sample_id,
                    "batch_id": batch_id,
                    "files_processed": batch_result.get(
                        "files_processed", len(parquet_paths)
                    ),
                    "total_files": batch_result.get("total_files", len(parquet_paths)),
                },
            )

            logger.info(
                f"Completed PanNanoDX analysis batch processing: {batch_size} files for sample '{sample_id}'"
            )
            logger.info(
                f"Files successfully processed: {batch_result.get('files_processed', len(parquet_paths))}/{batch_result.get('total_files', len(parquet_paths))}"
            )
            logger.info(
                f"Total features processed: {batch_result.get('total_features', 0)}"
            )

            if batch_result.get("error_message"):
                if suppress_expected:
                    logger.warning(
                        f"Batch processing completed with expected errors (fail-only BAM submission): {batch_result['error_message']}"
                    )
                    job.context.add_result(
                        "pannanodx_analysis",
                        {
                            "status": "expected_failure",
                            "error_message": batch_result["error_message"],
                            "sample_id": sample_id,
                        },
                    )
                else:
                    logger.error(
                        f"Batch processing completed with errors: {batch_result['error_message']}"
                    )
                    job.context.add_error(
                        "pannanodx_analysis", batch_result["error_message"]
                    )
            else:
                logger.info(
                    "Batch processing completed successfully with aggregated PanNanoDX analysis"
                )
                job.context.add_result(
                    "pannanodx_analysis",
                    {
                        "status": "success",
                        "sample_id": sample_id,
                        "analysis_time": batch_result.get("analysis_timestamp", 0),
                        "model_used": batch_result.get(
                            "model_used", "pancan_devel_v5i_NN_v2.pkl"
                        ),
                        "n_features": batch_result.get("total_features", 0),
                        "store_path": batch_result.get("store_path", ""),
                        "processing_steps": batch_result.get("processing_steps", []),
                        "files_processed": batch_result.get(
                            "files_processed", len(parquet_paths)
                        ),
                        "total_files": batch_result.get(
                            "total_files", len(parquet_paths)
                        ),
                    },
                )

            return

        else:
            # Single file processing (backward compatibility)
            # Get metadata from previous steps
            bed_conversion_result = job.context.results.get("bed_conversion", {})
            bam_metadata = job.context.metadata.get("bam_metadata", {})

            # Extract parquet path from bed conversion results
            parquet_path = bed_conversion_result.get("parquet_path")
            if not parquet_path:
                raise ValueError("No parquet path found from bed conversion step")

            logger.info(
                f"Starting PanNanoDX analysis for: {os.path.basename(parquet_path)}"
            )

            # Get sample ID
            sample_id = bam_metadata.get("sample_id", "unknown")

            # Determine work directory
            if work_dir is None:
                work_dir = os.path.dirname(parquet_path)
            else:
                os.makedirs(work_dir, exist_ok=True)
                logger.debug(f"Using specified work directory: {work_dir}")

            # Create PanNanoDX analysis instance
            pannanodx_analyzer = PanNanodxAnalysis(work_dir=work_dir)

            # Process the parquet file
            nanodx_result = pannanodx_analyzer.process_parquet_file(
                parquet_path, sample_id
            )

            # Store results in job context
            job.context.add_metadata("pannanodx_analysis", nanodx_result.results)
            job.context.add_metadata(
                "pannanodx_processing_steps", nanodx_result.processing_steps
            )

            if nanodx_result.error_message:
                job.context.add_error("pannanodx_analysis", nanodx_result.error_message)
                logger.error(
                    f"PanNanoDX analysis failed: {nanodx_result.error_message}"
                )
            else:
                job.context.add_result(
                    "pannanodx_analysis",
                    {
                        "status": "success",
                        "sample_id": nanodx_result.sample_id,
                        "analysis_time": nanodx_result.analysis_timestamp,
                        "model_used": nanodx_result.model_used,
                        "n_features": nanodx_result.n_features,
                        "store_path": nanodx_result.results.get("store_path"),
                        "processing_steps": nanodx_result.processing_steps,
                    },
                )
                logger.info(
                    f"PanNanoDX analysis complete for {os.path.basename(parquet_path)}"
                )
                logger.info(f"Sample ID: {nanodx_result.sample_id}")
                logger.info(f"Model used: {nanodx_result.model_used}")
                logger.info(f"Number of features: {nanodx_result.n_features}")
                logger.debug(
                    f"Processing steps: {', '.join(nanodx_result.processing_steps)}"
                )

    except Exception as e:
        if suppress_expected:
            logger.warning(
                f"Expected PanNanoDX failure for fail-only BAM submission: {e}"
            )
            job.context.add_result(
                "pannanodx_analysis",
                {"status": "expected_failure", "error_message": str(e)},
            )
            return
        job.context.add_error("pannanodx_analysis", str(e))
        logger.error(f"Error in PanNanoDX analysis for {job.context.filepath}: {e}")


# Command line interface for classification subprocess
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="NanoDX classification")
    parser.add_argument("command", choices=["classify"])
    parser.add_argument("--model", required=True, help="Model file path")
    parser.add_argument("--input", required=True, help="Input CSV file")
    parser.add_argument("--output", required=True, help="Output JSON file")

    args = parser.parse_args()

    if args.command == "classify":
        # Read input data
        test_df = pd.read_csv(args.input)

        # Run classification
        predictions, class_labels, n_features = _direct_classification(
            args.model, test_df
        )

        # Save results
        result_data = {
            "predictions": predictions.tolist(),
            "class_labels": class_labels.tolist(),
            "n_features": n_features,
        }

        with open(args.output, "w") as f:
            json.dump(result_data, f)
