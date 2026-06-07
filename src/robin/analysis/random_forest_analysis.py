#!/usr/bin/env python3
"""
Random Forest methylation classification analysis module for robin.

This module provides functionality for processing BAM files to analyze methylation patterns
using a Random Forest classifier. It includes capabilities for:

- Real-time BAM file processing
- Methylation data extraction using modkit
- Random Forest classification
- Time series analysis of classification confidence

The module integrates with the robin workflow system and uses temporary files for efficient data processing.

Dependencies
-----------
- pysam
- pandas
- R (with required packages)
- modkit
- robin package

Notes
-----
The module requires proper configuration of input/output directories and
assumes the presence of necessary R scripts and model files.
"""

import os
import tempfile
import time
import subprocess
import logging
import shutil
from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List

import pandas as pd

from robin.logging_config import get_job_logger

# Import robin utilities
try:
    from robin.analysis.utilities.merge_bedmethyl import (
        collapse_bedmethyl,
        load_minimal_modkit_data,
        reconstruct_full_bedmethyl_data,
    )
except ImportError as e:
    logging.warning(f"Some robin dependencies not available: {e}")
    collapse_bedmethyl = None
    load_minimal_modkit_data = None
    reconstruct_full_bedmethyl_data = None


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


def _compute_hvpath() -> Optional[str]:
    """Resolve the hv_rapidCNS2 path from this package's submodules.

    Looks up: robin/submodules/hv_rapidCNS2
    Returns the path if it exists, else None.
    """
    base_dir = os.path.dirname(os.path.dirname(__file__))
    candidate = os.path.join(base_dir, "submodules", "hv_rapidCNS2")
    if os.path.exists(candidate):
        return candidate
    return None


HVPATH = _compute_hvpath()




@dataclass
class RandomForestMetadata:
    """Container for Random Forest analysis metadata and results"""

    sample_id: str
    parquet_path: str
    analysis_timestamp: float
    batch_number: int
    votes_file_path: Optional[str] = None
    scores_file_path: Optional[str] = None
    bed_file_path: Optional[str] = None
    processing_steps: List[str] = field(default_factory=list)
    error_message: Optional[str] = None
    results: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        if self.processing_steps is None:
            self.processing_steps = []
        if self.results is None:
            self.results = {}


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
    logger = logging.getLogger("robin.random_forest")

    try:
        logger.info(f"Starting run_rcns2 with bed file: {bed}")
        logger.info(f"Output directory: {rcns2folder}")
        logger.info(f"Batch number: {batch}")

        # Check if R script exists
        r_script_path = f"{HVPATH}/bin/methylation_classification_nanodx_v0.3.R"
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
        try:
            result = subprocess.run(
                command.split(), capture_output=True, text=True, check=True
            )
            logger.info("R script executed successfully")
            logger.debug(f"R script stdout: {result.stdout}")
            if result.stderr and showerrors:
                logger.warning(f"R script stderr: {result.stderr}")
        except subprocess.CalledProcessError as e:
            logger.error(f"R script failed with return code {e.returncode}")
            logger.error(f"stdout: {e.stdout}")
            logger.error(f"stderr: {e.stderr}")
            # Always log stderr on failure, regardless of showerrors setting
            if e.stderr:
                logger.error(f"R script error output: {e.stderr}")
            raise

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
    if load_minimal_modkit_data is None:
        raise ImportError(
            "load_minimal_modkit_data function not available from robin package"
        )
    return load_minimal_modkit_data(parquet_path)


class RandomForestAnalysis:
    """
    Random Forest methylation classification analysis worker.

    This class provides specialized functionality for Random Forest methylation analysis.
    It handles real-time processing of BAM files, methylation data extraction, and
    classification results.
    """

    def __init__(self, work_dir=None, threads=4, showerrors=False):
        self.work_dir = work_dir or os.getcwd()
        self.threads = threads
        self.showerrors = showerrors

        # Initialize batch tracking
        self.bambatch = {}

        # Use logger for initialization
        logger = logging.getLogger("robin.random_forest")
        logger.info("Random Forest Analysis initialized")
        logger.debug(f"Work directory: {self.work_dir}")
        logger.debug(f"Threads: {self.threads}")
        logger.debug(f"Show errors: {self.showerrors}")

        if HVPATH is None:
            logger.warning("HVPATH not available - R scripts may not work")
        else:
            logger.debug(f"HVPATH: {HVPATH}")

    def process_parquet_file(
        self, parquet_path: str, sample_id: str
    ) -> RandomForestMetadata:
        """Process a single parquet file for Random Forest analysis"""
        logger = logging.getLogger("robin.random_forest")

        logger.info(
            f"Processing parquet file: {parquet_path} for Random Forest analysis"
        )
        start_time = time.time()

        # Initialize batch number if not exists
        if sample_id not in self.bambatch:
            self.bambatch[sample_id] = 1

        random_forest_result = RandomForestMetadata(
            sample_id=sample_id,
            parquet_path=parquet_path,
            analysis_timestamp=start_time,
            batch_number=self.bambatch[sample_id],
        )

        logger.info(f"Starting Random Forest analysis for sample: {sample_id}")
        logger.debug(f"Parquet file: {os.path.basename(parquet_path)}")
        logger.debug(f"Batch number: {self.bambatch[sample_id]}")

        try:
            # Step 1: Check if parquet file exists and is readable
            if not os.path.exists(parquet_path):
                raise FileNotFoundError(f"Parquet file not found: {parquet_path}")

            random_forest_result.processing_steps.append("file_validation")

            # Step 2: Load modkit data from parquet
            logger.info("Loading modkit data from parquet...")
            merged_modkit_df = load_modkit_data(parquet_path)
            random_forest_result.processing_steps.append("modkit_data_loaded")

            # Step 3: Reconstruct full bedmethyl data for RandomForest
            logger.info("Reconstructing full bedmethyl data...")
            full_modkit_df = reconstruct_full_bedmethyl_data(merged_modkit_df)
            random_forest_result.processing_steps.append("bedmethyl_reconstructed")

            # Step 4: Collapse bedmethyl data
            logger.info("Collapsing bedmethyl data...")
            forest_dx = collapse_bedmethyl(full_modkit_df)
            merged_modkit_df = forest_dx
            random_forest_result.processing_steps.append("bedmethyl_collapsed")

            if merged_modkit_df is None:
                raise ValueError("Failed to load modkit data from parquet file")

            logger.info(f"Loaded modkit data with shape: {merged_modkit_df.shape}")

            # Step 5: Create sample output directory
            sample_dir = os.path.join(self.work_dir, sample_id)
            os.makedirs(sample_dir, exist_ok=True)
            logger.info(f"Created sample output directory: {sample_dir}")

            # Step 6: Create temporary directory for R script output
            rcns2folder = tempfile.mkdtemp(dir=sample_dir)
            logger.info(f"Created temporary directory for R output: {rcns2folder}")

            try:
                # Step 7: Write the BED file
                randomforest_bed_output = os.path.join(
                    sample_dir, "RandomForestBed.bed"
                )
                logger.info(f"Writing BED file to: {randomforest_bed_output}")

                # Create BED file format expected by the R script
                # The R script expects columns 1:3, 6, 10:11 to be:
                # Column 1: chr, Column 2: start, Column 3: end, Column 6: strand, Column 10: cov, Column 11: methylation_percent
                bed_df = merged_modkit_df.copy()

                # Create a BED file with the correct number of columns
                bed_df = bed_df.rename(
                    columns={
                        "chrom": "chr",
                        "start_pos": "start",
                        "end_pos": "end",
                        "Nvalid": "cov",
                        "score": "methylation_percent",  # The reconstructed data has 'score' not 'fraction'
                    }
                )

                # Create a BED file with at least 11 columns
                # The R script reads: columns 1:3 (chr, start, end), column 6 (strand), columns 10:11 (cov, methylation_percent)
                final_bed_df = pd.DataFrame()
                final_bed_df["chr"] = bed_df["chr"]  # Column 1
                final_bed_df["start"] = bed_df["start"]  # Column 2
                final_bed_df["end"] = bed_df["end"]  # Column 3
                final_bed_df["dummy1"] = "."  # Column 4 (dummy)
                final_bed_df["dummy2"] = 0  # Column 5 (dummy)
                final_bed_df["strand"] = bed_df["strand"]  # Column 6
                final_bed_df["dummy3"] = "."  # Column 7 (dummy)
                final_bed_df["dummy4"] = 0  # Column 8 (dummy)
                final_bed_df["dummy5"] = "."  # Column 9 (dummy)
                final_bed_df["cov"] = bed_df["cov"]  # Column 10
                final_bed_df["methylation_percent"] = bed_df[
                    "methylation_percent"
                ]  # Column 11

                # Write with header as expected by the R script
                final_bed_df.to_csv(
                    randomforest_bed_output, sep="\t", index=False, header=True
                )

                # Validate BED file format
                logger.info(f"BED file created with shape: {final_bed_df.shape}")
                logger.info(f"BED file columns: {final_bed_df.columns.tolist()}")

                # Check if we have the required columns for R script
                required_columns = [
                    "chr",
                    "start",
                    "end",
                    "strand",
                    "cov",
                    "methylation_percent",
                ]
                missing_required = [
                    col for col in required_columns if col not in final_bed_df.columns
                ]
                if missing_required:
                    raise ValueError(
                        f"Missing required columns for R script: {missing_required}"
                    )

                # Check if BED file is not empty
                if final_bed_df.empty:
                    raise ValueError(
                        "BED file is empty - no methylation data to process"
                    )

                # Verify the file was written correctly
                if os.path.exists(randomforest_bed_output):
                    file_size = os.path.getsize(randomforest_bed_output)
                    logger.info(
                        f"BED file written successfully, size: {file_size} bytes"
                    )
                else:
                    raise FileNotFoundError("BED file was not created")

                random_forest_result.bed_file_path = randomforest_bed_output
                random_forest_result.processing_steps.append("bed_file_created")

                # Step 8: Run the R script
                logger.info("Attempting to run R script...")
                try:
                    run_rcns2(
                        rcns2folder,
                        self.bambatch[sample_id],
                        randomforest_bed_output,
                        self.threads,
                        self.showerrors,
                    )
                    logger.info("R script completed successfully")
                    random_forest_result.processing_steps.append("r_script_completed")
                except Exception as e:
                    logger.error(f"Error running R script: {str(e)}", exc_info=True)
                    raise

                # Step 9: Process results
                votes_file = f"{rcns2folder}/live_{self.bambatch[sample_id]}_votes.tsv"
                if os.path.isfile(votes_file):
                    logger.info(f"Found votes file: {votes_file}")
                    scores = pd.read_table(votes_file, sep=r"\s+")
                    logger.debug(f"Original scores shape: {scores.shape}")
                    logger.debug(f"Original scores columns: {list(scores.columns)}")
                    logger.debug(f"Original scores head:\n{scores.head()}")

                    # Process scores correctly - drop Freq column and transpose
                    scores_to_save = scores.drop(columns=["Freq"]).T
                    logger.debug(f"After transpose shape: {scores_to_save.shape}")
                    logger.debug(
                        f"After transpose columns: {list(scores_to_save.columns)}"
                    )
                    logger.debug(f"After transpose index: {list(scores_to_save.index)}")

                    # Add timestamp column
                    scores_to_save["timestamp"] = (
                        start_time * 1000
                    )  # Convert to milliseconds

                    # Reorder columns to put timestamp first
                    cols = scores_to_save.columns.tolist()
                    cols.insert(0, cols.pop(cols.index("timestamp")))
                    scores_to_save = scores_to_save[cols]

                    # Load existing results if available and accumulate
                    output_file = os.path.join(sample_dir, "random_forest_scores.csv")
                    if os.path.exists(output_file):
                        try:
                            logger.info(f"Loading existing results from: {output_file}")
                            existing_scores = pd.read_csv(output_file)
                            logger.info(
                                f"Found {len(existing_scores)} existing results"
                            )

                            # Append new results to existing data
                            combined_scores = pd.concat(
                                [existing_scores, scores_to_save], ignore_index=True
                            )
                            logger.info(
                                f"Combined results: {len(combined_scores)} total entries"
                            )

                            # Save accumulated results
                            combined_scores.to_csv(output_file, index=False)
                            logger.info(f"Saved accumulated results to: {output_file}")

                        except Exception as e:
                            logger.warning(
                                f"Error loading existing results, saving new results only: {e}"
                            )
                            scores_to_save.to_csv(output_file, index=False)
                            logger.info(f"Saved new results to: {output_file}")
                    else:
                        # First run - save new results
                        scores_to_save.to_csv(output_file, index=False)
                        logger.info(f"Saved initial results to: {output_file}")

                    random_forest_result.votes_file_path = votes_file
                    random_forest_result.scores_file_path = output_file
                    random_forest_result.processing_steps.append("results_saved")

                    # Store results in metadata
                    random_forest_result.results = {
                        "analysis_time": time.time() - start_time,
                        "batch_number": self.bambatch[sample_id],
                        "votes_file": votes_file,
                        "scores_file": output_file,
                        "bed_file": randomforest_bed_output,
                        "scores_shape": scores_to_save.shape,
                        "processing_steps": random_forest_result.processing_steps.copy(),
                    }

                    # Increment batch number for next run
                    self.bambatch[sample_id] += 1

                    logger.info(
                        f"Random Forest analysis completed for {sample_id} in {time.time() - start_time:.2f}s"
                    )

                else:
                    raise FileNotFoundError(f"Votes file not found: {votes_file}")

            finally:
                # Clean up temporary directory
                try:
                    if os.path.exists(rcns2folder):
                        shutil.rmtree(rcns2folder, ignore_errors=True)
                        logger.debug(f"Cleaned up temporary directory: {rcns2folder}")
                except Exception as e:
                    logger.warning(
                        f"Failed to clean up temporary directory {rcns2folder}: {e}"
                    )

        except Exception as e:
            error_msg = f"Random Forest analysis failed for {sample_id}: {str(e)}"
            logger.error(error_msg, exc_info=True)
            random_forest_result.error_message = error_msg
            random_forest_result.processing_steps.append("analysis_failed")

        return random_forest_result


def process_multiple_files(
    parquet_paths,
    metadata_list,
    work_dir,
    logger,
    threads=4,
    showerrors=False,
    job_id=None,
):
    """
    Process multiple parquet files for Random Forest analysis.
    
    This function processes multiple parquet files for the same sample.
    Each file is processed individually and results are accumulated in
    the same random_forest_scores.csv file.

    Args:
        parquet_paths: List of paths to parquet files
        metadata_list: List of metadata dictionaries (one per parquet file)
        work_dir: Working directory
        logger: Logger instance
        threads: Number of threads for R script execution
        showerrors: Whether to show R script error messages

    Returns:
        Dictionary with aggregated Random Forest analysis results
    """
    if not parquet_paths or not metadata_list:
        raise ValueError("parquet_paths and metadata_list must not be empty")
    
    if len(parquet_paths) != len(metadata_list):
        raise ValueError("parquet_paths and metadata_list must have the same length")
    
    # Get sample ID from first metadata (assuming all parquet files are from same sample)
    sample_id = metadata_list[0].get("sample_id", "unknown")
    
    logger.info(f"🌲 Starting multi-file Random Forest analysis for sample: {sample_id}")
    logger.info(f"Processing {len(parquet_paths)} parquet files for sample {sample_id}")
    logger.info(f"Using {threads} threads for R script execution")
    
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
        "total_batches": 0,
        "scores_file_path": None,
        "bed_file_path": None,
    }

    try:
        # Create sample-specific output directory
        sample_dir = os.path.join(work_dir, sample_id)
        os.makedirs(sample_dir, exist_ok=True)
        
        scores_file_path = os.path.join(sample_dir, "random_forest_scores.csv")
        analysis_result["scores_file_path"] = scores_file_path
        
        logger.info(f"Created output directory: {sample_dir}")
        logger.info(f"Scores file: {scores_file_path}")
        analysis_result["processing_steps"].append("directory_created")

        # Initialize Random Forest analysis
        random_forest_analyzer = RandomForestAnalysis(
            work_dir=work_dir, 
            threads=threads, 
            showerrors=showerrors
        )
        
        logger.info("Initialized Random Forest analyzer")
        analysis_result["processing_steps"].append("analyzer_initialized")

        # Process each parquet file individually
        logger.info("Processing parquet files individually")
        processed_files = 0
        total_batches = 0
        
        for i, (parquet_path, metadata) in enumerate(zip(parquet_paths, metadata_list)):
            logger.info(f"Processing parquet file {i+1}/{len(parquet_paths)}: {os.path.basename(parquet_path)}")
            
            try:
                # Check if parquet file exists
                if not os.path.exists(parquet_path):
                    logger.warning(f"Parquet file not found: {parquet_path}")
                    continue
                
                # Process the parquet file
                rf_result = random_forest_analyzer.process_parquet_file(parquet_path, sample_id)
                
                if rf_result.error_message:
                    logger.warning(f"Error processing {os.path.basename(parquet_path)}: {rf_result.error_message}")
                    continue
                
                processed_files += 1
                total_batches += 1
                logger.debug(f"Successfully processed file {i+1}: {os.path.basename(parquet_path)}")
                logger.debug(f"Batch number: {rf_result.batch_number}")
                if job_id is not None:
                    try:
                        from robin.workflow_ray import notify_coordinator_files_completed

                        notify_coordinator_files_completed(
                            "random_forest", 1, job_id=job_id
                        )
                    except Exception:
                        pass
                
            except Exception as e:
                logger.warning(f"Error processing {os.path.basename(parquet_path)}: {e}")
                continue

        if processed_files == 0:
            analysis_result["error_message"] = "No files could be processed successfully"
            analysis_result["processing_steps"].append("no_files_processed")
            return analysis_result

        analysis_result["files_processed"] = processed_files
        analysis_result["total_batches"] = total_batches
        analysis_result["processing_steps"].append("files_processed")

        # Check if scores file was created
        if os.path.exists(scores_file_path):
            analysis_result["processing_steps"].append("scores_file_created")
            logger.info(f"Random Forest scores accumulated in: {scores_file_path}")
            
            # Load and log summary of results
            try:
                scores_df = pd.read_csv(scores_file_path)
                logger.info(f"Scores file contains {len(scores_df)} entries")
                logger.debug(f"Scores file shape: {scores_df.shape}")
                logger.debug(f"Scores file columns: {scores_df.columns.tolist()}")
            except Exception as e:
                logger.warning(f"Could not read scores file for summary: {e}")
        else:
            logger.warning("No random_forest_scores.csv file was created")

        # Set bed file path (from last processed file)
        bed_file_path = os.path.join(sample_dir, "RandomForestBed.bed")
        if os.path.exists(bed_file_path):
            analysis_result["bed_file_path"] = bed_file_path
            logger.debug(f"BED file available: {bed_file_path}")

        analysis_result["processing_steps"].append("analysis_complete")
        logger.info(f"Multi-file Random Forest analysis completed for {sample_id}")
        logger.info(f"Files successfully processed: {analysis_result['files_processed']}/{analysis_result['total_files']}")
        logger.info(f"Total batches processed: {analysis_result['total_batches']}")
        logger.info(f"Output directory: {sample_dir}")
        
        return analysis_result

    except Exception as e:
        logger.error(f"Error in multi-file Random Forest analysis for {sample_id}: {e}")
        analysis_result["error_message"] = str(e)
        analysis_result["processing_steps"].append("analysis_failed")
        return analysis_result


def random_forest_handler(job, work_dir=None):
    """
    Handler function for Random Forest analysis jobs.
    This function processes parquet files for Random Forest methylation analysis.

    Args:
        job: The workflow job containing file and metadata
        work_dir: Optional base directory for output (defaults to parquet file directory)
    """
    # Get job-specific logger
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
    suppress_expected = _is_fail_only_expected(job)

    try:
        # Check if this is a batched job
        batched_job = job.context.metadata.get("_batched_job")
        if batched_job:
            batch_size = batched_job.get_file_count()
            sample_id = batched_job.get_sample_id()
            batch_id = batched_job.batch_id
            logger.info(f"Processing Random Forest analysis batch: {batch_size} files for sample '{sample_id}' (batch_id: {batch_id})")
            
            # Get all filepaths in the batch
            filepaths = batched_job.get_filepaths()
            
            # Log individual files in the batch
            for i, filepath in enumerate(filepaths):
                logger.info(f"  Batch file {i+1}/{batch_size}: {os.path.basename(filepath)}")
            
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
                bed_conversion_result = batched_job.contexts[i].results.get("bed_conversion", {})
                parquet_path = bed_conversion_result.get("parquet_path")
                
                if parquet_path:
                    parquet_paths.append(parquet_path)
                    metadata_list.append(file_metadata)
                    logger.debug(f"Found parquet path for file {i+1}: {os.path.basename(parquet_path)}")
                else:
                    logger.warning(f"No parquet path found for file {i+1}: {os.path.basename(bam_path)}")
            
            if not parquet_paths:
                error_msg = "No parquet paths found from bed conversion results in batch"
                if suppress_expected:
                    logger.warning(f"{error_msg} (expected for fail-only BAM submission)")
                    job.context.add_result(
                        "random_forest_analysis",
                        {"status": "expected_failure", "reason": error_msg},
                    )
                else:
                    logger.error(error_msg)
                    job.context.add_error("random_forest_analysis", error_msg)
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
            logger.info(f"Processing {len(parquet_paths)} parquet files as aggregated batch for sample '{sample_id}'")
            batch_result = process_multiple_files(
                parquet_paths=parquet_paths,
                metadata_list=metadata_list,
                work_dir=batch_work_dir,
                logger=logger,
                threads=4,  # Default thread count
                showerrors=False,  # Default error display setting
                job_id=job.job_id,
            )
            
            # Store batch results in job context (maintain compatibility with existing structure)
            job.context.add_metadata("random_forest_analysis", {
                "batch_result": batch_result,  # Single aggregated result
                "batch_size": batch_size,
                "sample_id": sample_id,
                "batch_id": batch_id,
                "files_processed": batch_result.get("files_processed", len(parquet_paths)),
                "total_files": batch_result.get("total_files", len(parquet_paths))
            })
            
            logger.info(f"Completed Random Forest analysis batch processing: {batch_size} files for sample '{sample_id}'")
            logger.info(f"Files successfully processed: {batch_result.get('files_processed', len(parquet_paths))}/{batch_result.get('total_files', len(parquet_paths))}")
            logger.info(f"Total batches processed: {batch_result.get('total_batches', 0)}")
            
            if batch_result.get("error_message"):
                if suppress_expected:
                    logger.warning(
                        f"Batch processing completed with expected errors (fail-only BAM submission): {batch_result['error_message']}"
                    )
                    job.context.add_result(
                        "random_forest_analysis",
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
                        "random_forest_analysis", batch_result["error_message"]
                    )
            else:
                logger.info("Batch processing completed successfully with aggregated Random Forest analysis")
                job.context.add_result(
                    "random_forest_analysis",
                    {
                        "status": "success",
                        "sample_id": sample_id,
                        "analysis_time": batch_result.get("analysis_timestamp", 0),
                        "total_batches": batch_result.get("total_batches", 0),
                        "processing_steps": batch_result.get("processing_steps", []),
                        "scores_file": batch_result.get("scores_file_path", ""),
                        "bed_file": batch_result.get("bed_file_path", ""),
                        "files_processed": batch_result.get("files_processed", len(parquet_paths)),
                        "total_files": batch_result.get("total_files", len(parquet_paths)),
                    },
                )
            
            return
            
        else:
            # Single file processing (backward compatibility)
            # Get the parquet file path from bed_conversion results
            bed_conversion_results = job.context.results.get("bed_conversion", {})
            parquet_path = bed_conversion_results.get("parquet_path")

            if not parquet_path:
                raise ValueError("No parquet file path found in bed_conversion results")

            logger.info(
                f"Starting Random Forest analysis for: {os.path.basename(parquet_path)}"
            )

            # Get metadata from preprocessing
            bam_metadata = job.context.metadata.get("bam_metadata", {})
            sample_id = bam_metadata.get("sample_id", "unknown")

            # Determine work directory
            if work_dir is None:
                # Default to parquet file directory
                work_dir = os.path.dirname(parquet_path)
            else:
                # Use specified work directory, create if it doesn't exist
                os.makedirs(work_dir, exist_ok=True)
                logger.debug(f"Using specified work directory: {work_dir}")

            # Create Random Forest analysis instance
            random_forest_analyzer = RandomForestAnalysis(work_dir=work_dir)

            # Process the parquet file
            result = random_forest_analyzer.process_parquet_file(parquet_path, sample_id)

            # Store results in job context
            job.context.add_metadata("random_forest_analysis", result.results)
            job.context.add_metadata(
                "random_forest_processing_steps", result.processing_steps
            )

            if result.error_message:
                job.context.add_error("random_forest_analysis", result.error_message)
                logger.error(f"Random Forest analysis failed: {result.error_message}")
            else:
                job.context.add_result(
                    "random_forest_analysis",
                    {
                        "status": "success",
                        "sample_id": result.sample_id,
                        "analysis_time": result.results.get("analysis_time", 0),
                        "batch_number": result.batch_number,
                        "processing_steps": result.processing_steps,
                        "votes_file": result.votes_file_path,
                        "scores_file": result.scores_file_path,
                        "bed_file": result.bed_file_path,
                    },
                )
                logger.info(
                    f"Random Forest analysis complete for {os.path.basename(parquet_path)}"
                )
                try:
                    from robin.workflow_ray import notify_coordinator_files_completed

                    notify_coordinator_files_completed(
                        "random_forest", 1, job_id=job.job_id
                    )
                except Exception:
                    pass
                logger.info(f"Sample ID: {result.sample_id}")
                logger.info(f"Batch Number: {result.batch_number}")
                logger.debug(f"Processing steps: {', '.join(result.processing_steps)}")
                logger.debug(
                    f"Output directory: {os.path.dirname(result.scores_file_path or '')}"
                )

    except Exception as e:
        if suppress_expected:
            logger.warning(
                f"Expected Random Forest failure for fail-only BAM submission: {e}"
            )
            job.context.add_result(
                "random_forest_analysis",
                {"status": "expected_failure", "error_message": str(e)},
            )
            return
        job.context.add_error("random_forest_analysis", str(e))
        logger.error(f"Error in Random Forest analysis for {job.context.filepath}: {e}")
