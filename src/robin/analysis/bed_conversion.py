#!/usr/bin/env python3
"""
BAM to parquet file conversion module for robin.

Requires Python 3.12+. Automated conversion of BAM files to parquet using
matkit; integrates with robin's preprocessing pipeline and CPG master file.
"""
from __future__ import annotations

import sys
if sys.version_info < (3, 12):
    raise RuntimeError("robin bed_conversion requires Python 3.12 or newer")

import os
import time
import tempfile
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional

# Import robin utilities and resources
try:
    from robin.analysis.utilities.matkit import run_matkit
    from robin import resources
except ImportError:
    run_matkit = None
    resources = None

# Import local utilities for merge_modkit_files
try:
    from robin.analysis.temp_utilities import merge_modkit_files
except ImportError:
    merge_modkit_files = None

from robin.logging_config import get_job_logger

# Module-level logger to avoid repeated lookups
_LOGGER = logging.getLogger("robin.analysis.bed_conversion")

# Cache CPGs master file path across instances
_CPGS_MASTER_FILE_CACHE: Optional[str] = None


def _is_fail_only_expected(job) -> bool:
    """True when conversion 'failures' are expected for fail-only BAM submissions."""
    try:
        if not bool(job.context.metadata.get("fail_only_bam_submission", False)):
            return False
        fp = getattr(job.context, "filepath", "") or ""
        base = os.path.basename(fp).lower()
        return ("fail" in base) and ("pass" not in base)
    except Exception:
        return False


@dataclass(slots=True)
class BedConversionMetadata:
    """Container for BAM to parquet conversion metadata and results."""

    sample_id: str
    bam_path: str
    analysis_timestamp: float
    parquet_path: Optional[str] = None
    processing_steps: List[str] = field(default_factory=list)
    error_message: Optional[str] = None
    results: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        if self.processing_steps is None:
            self.processing_steps = []
        if self.results is None:
            self.results = {}


class BedConversionAnalysis:
    """BAM to parquet conversion analysis worker"""

    def __init__(self, work_dir=None, threads=1):
        self.work_dir = work_dir or os.getcwd()
        self.threads = threads
        
        # Find required files and paths
        self.cpgs_master_file = self._find_cpgs_master_file()

        # Initialize counter for incremental file naming
        self.file_counter = 1

        # Use logger for initialization
        self.logger = _LOGGER
        self.logger.info("BAM to parquet conversion analysis initialized")
        self.logger.debug(f"Work directory: {self.work_dir}")
        self.logger.debug(f"CPGs master file: {self.cpgs_master_file}")
        self.logger.debug(f"Threads: {self.threads}")

    def _find_cpgs_master_file(self) -> str:
        """Find the parquet filter file (CPGs to retain) from robin resources.
        Only parquet_filter.txt is used; no fallback to other files."""
        global _CPGS_MASTER_FILE_CACHE
        if _CPGS_MASTER_FILE_CACHE is not None:
            return _CPGS_MASTER_FILE_CACHE

        if resources is None:
            raise FileNotFoundError(
                "robin.resources not available; cannot locate parquet_filter.txt"
            )

        resources_dir = os.path.dirname(os.path.abspath(resources.__file__))
        parquet_filter = os.path.join(resources_dir, "parquet_filter.txt")
        if not os.path.exists(parquet_filter):
            raise FileNotFoundError(
                f"parquet_filter.txt not found at {parquet_filter}. "
                "Bed conversion requires parquet_filter.txt in robin resources."
            )

        _CPGS_MASTER_FILE_CACHE = parquet_filter
        return _CPGS_MASTER_FILE_CACHE

    def _get_next_file_number(self) -> int:
        """Get the next file number for incremental naming"""
        current = self.file_counter
        self.file_counter += 1
        return current

    def process_bam_file(
        self, bam_path: str, metadata: Dict[str, Any]
    ) -> BedConversionMetadata:
        """Process a single BAM file for parquet conversion"""
        logger = self.logger
        logger.info(f"Processing BAM file: {bam_path} for parquet conversion")
        sample_id = metadata.get("sample_id", "unknown")
        start_time = time.time()

        # Get next file number for incremental naming
        file_number = self._get_next_file_number()

        bed_result = BedConversionMetadata(
            sample_id=sample_id, bam_path=bam_path, analysis_timestamp=start_time
        )

        logger.info(f"Starting BAM to parquet conversion for sample: {sample_id}")
        logger.debug(f"BAM file: {os.path.basename(bam_path)}")
        logger.debug(f"Work directory: {self.work_dir}")

        try:
            # Step 1: Check if BAM file exists and is readable
            if not os.path.exists(bam_path):
                raise FileNotFoundError(f"BAM file not found: {bam_path}")

            bed_result.processing_steps.append("file_validation")

            # Step 2: Create sample-specific directory
            sample_dir = os.path.join(self.work_dir, sample_id)
            os.makedirs(sample_dir, exist_ok=True)
            logger.info(f"Created sample directory: {sample_dir}")

            bed_result.processing_steps.append("directory_created")

            # Step 3: Process BAM file with matkit
            logger.info("Converting BAM to parquet using matkit...")

            # Create output path specific to this sample
            parquet_path = os.path.join(
                self.work_dir,
                sample_id,
                f"{sample_id}.parquet",  # Use sampleID for the output filename
            )

            # Process the BAM file - this follows the exact same pattern as the working code
            data = self._process_bams([bam_path], sample_dir)

            # Update state with new data - this creates the parquet file using the same approach
            state = self._update_state(parquet_path, data, sample_id, file_number)

            bed_result.parquet_path = state
            bed_result.processing_steps.append("parquet_conversion_complete")

            logger.info(f"BAM to parquet conversion complete for {sample_id}")
            logger.info(f"   Parquet file: {parquet_path}")
            logger.info(f"   Processed files: {len(data)}")

            return bed_result

        except Exception as e:
            logger.error(f"Error in BAM to parquet conversion for {sample_id}: {e}")
            bed_result.error_message = str(e)
            bed_result.processing_steps.append("conversion_failed")
            return bed_result

    def _process_bams(self, bams: List[str], work_dir: str) -> List[str]:
        """Process BAM files using matkit and return list of processed file paths"""
        logger = self.logger
        processed_files: List[str] = []
        run_matkit_callable = run_matkit  # local binding
        threads = max(1, int(self.threads or 1))

        def cleanup_temp_files(files: List[str]) -> None:
            for file_path in files:
                try:
                    os.remove(file_path)
                    logger.debug(f"Deleted temporary file: {file_path}")
                except OSError as e:
                    logger.warning(f"Failed to delete temporary file {file_path}: {e}")

        def process_single_bam(bam: str) -> str:
            logger.debug(f"Processing BAM file: {bam}")
            # Temporary parquet path: matkit writes 8-column parquet so merge_modkit_files can read directly (no CSV parse).
            temp_file = tempfile.NamedTemporaryFile(
                suffix=".parquet", prefix="modkit_", dir=work_dir, delete=False
            )
            temp_file.close()

            try:
                # Run matkit on the BAM file (writes parquet when path ends in .parquet)
                if run_matkit_callable is not None:
                    logger.debug("Running matkit callable")
                    run_matkit_callable(bam, temp_file.name)
                else:
                    # Fallback if robin is not available (write empty parquet not used in real runs)
                    import pyarrow as pa
                    import pyarrow.parquet as pq
                    empty = pa.table(
                        {
                            "chrom": pa.array([], type=pa.binary()),
                            "chromStart": pa.array([], type=pa.int64()),
                            "mod_code": pa.array([], type=pa.binary()),
                            "strand": pa.array([], type=pa.binary()),
                            "valid_cov": pa.array([], type=pa.uint32()),
                            "percent_modified": pa.array([], type=pa.float32()),
                            "n_mod": pa.array([], type=pa.uint32()),
                            "n_canonical": pa.array([], type=pa.uint32()),
                        }
                    )
                    pq.write_table(empty, temp_file.name)

                logger.debug(f"Successfully processed: {temp_file.name}")
                return temp_file.name

            except Exception as e:
                logger.error(f"Failed to process {bam}: {e}")
                # Clean up the temporary file
                try:
                    os.remove(temp_file.name)
                except Exception as e:
                    logger.error(f"Failed to delete temporary file {temp_file.name}: {e}")
                    pass
                raise

        if threads == 1 or len(bams) <= 1:
            return [process_single_bam(bam) for bam in bams]

        # Parallel path: tolerate per-file failures so one bad BAM does not fail the batch
        try:
            with ThreadPoolExecutor(max_workers=threads) as executor:
                futures = {executor.submit(process_single_bam, bam): bam for bam in bams}
                for future in as_completed(futures):
                    try:
                        processed_files.append(future.result())
                    except Exception as e:
                        bam_path = futures[future]
                        logger.warning(f"Failed to process BAM {os.path.basename(bam_path)}: {e}")
            return processed_files
        except Exception:
            cleanup_temp_files(processed_files)
            raise

    def _update_state(
        self, state: str, data: List[str], sample_id: str, file_number: int
    ) -> str:
        """Update the state with new data from BAM processing - creates parquet file directly.

        ``file_number`` is retained for callers (incremental naming); merge uses ``len(data)``
        as the number of BAM-derived matkit inputs.
        """
        logger = _LOGGER

        if merge_modkit_files is None:
            logger.error(
                "robin.analysis.temp_utilities.merge_modkit_files not available"
            )
            raise ImportError(
                "robin.analysis.temp_utilities.merge_modkit_files is required for parquet creation"
            )

        # Configuration for mnpflex - allow env var overrides
        mnpflex_config = {
            "mnpuser": os.getenv("MNPFLEX_USER"),
            "mnppass": os.getenv("MNPFLEX_PASS"),
        }

        logger.debug(f"Creating parquet file: {state}")
        logger.debug(
            "Processing %s matkit files (file_number=%s)",
            len(data),
            file_number,
        )

        # One matkit output file per BAM; merge needs the true count for metadata/cumulative stats
        num_bam_files = len(data)
        # Use merge_modkit_files exactly like the working code
        # This function is designed to create binary parquet files, not text files
        merge_modkit_files(
            data,
            state,  # Output parquet file path
            state,  # State parquet file path (same as output)
            self.cpgs_master_file,
            sample_id,
            self.work_dir,
            mnpflex_config,
            num_bam_files,
        )

        # When no reads passed the QS filter, all per-BAM parquets are empty and merge
        # returns without writing; treat as non-fatal and carry on so next batches can run.
        if not os.path.exists(state):
            logger.info(
                "No parquet written for this batch (no data after filtering); "
                "continuing so next batches can be processed."
            )
            for file_path in data:
                try:
                    os.remove(file_path)
                    logger.debug(f"Deleted temporary file: {file_path}")
                except OSError as e:
                    logger.warning(f"Failed to delete temporary file {file_path}: {e}")
            return state

        parquet_size = os.path.getsize(state)
        logger.debug(
            f"Parquet file created successfully: {state} ({parquet_size} bytes)"
        )

        # Clean up temporary files - exactly like the working code
        for file_path in data:
            try:
                os.remove(file_path)
                logger.debug(f"Deleted temporary file: {file_path}")
            except OSError as e:
                logger.warning(f"Failed to delete temporary file {file_path}: {e}")

        return state

def process_multiple_files(
    bam_paths, metadata_list, work_dir, logger, threads=4, expected_fail_only: bool = False
):
    """
    Process multiple BAM files for bed conversion analysis.
    
    This function processes multiple BAM files for the same sample.
    Each file is processed individually and results are accumulated in
    the same parquet file.

    Args:
        bam_paths: List of paths to BAM files
        metadata_list: List of metadata dictionaries (one per BAM file)
        work_dir: Working directory
        logger: Logger instance
        threads: Number of threads for processing

    Returns:
        Dictionary with aggregated bed conversion results
    """
    if not bam_paths or not metadata_list:
        raise ValueError("bam_paths and metadata_list must not be empty")
    
    if len(bam_paths) != len(metadata_list):
        raise ValueError("bam_paths and metadata_list must have the same length")
    
    # Get sample ID from first metadata (assuming all BAM files are from same sample)
    sample_id = metadata_list[0].get("sample_id", "unknown")
    
    logger.info(f"🔄 Starting multi-file bed conversion for sample: {sample_id}")
    logger.info(f"Processing {len(bam_paths)} BAM files for sample {sample_id}")
    logger.info(f"Using {threads} threads for processing")
    
    analysis_result = {
        "sample_id": sample_id,
        "bam_paths": bam_paths,
        "analysis_timestamp": time.time(),
        "processing_steps": [],
        "error_message": None,
        "files_processed": 0,
        "total_files": len(bam_paths),
        "parquet_path": None,
    }

    try:
        # Create sample-specific output directory
        sample_dir = os.path.join(work_dir, sample_id)
        os.makedirs(sample_dir, exist_ok=True)
        
        parquet_path = os.path.join(sample_dir, f"{sample_id}.parquet")
        analysis_result["parquet_path"] = parquet_path
        
        logger.info(f"Created output directory: {sample_dir}")
        logger.info(f"Parquet file: {parquet_path}")
        analysis_result["processing_steps"].append("directory_created")

        # Initialize bed conversion analysis
        bed_analyzer = BedConversionAnalysis(work_dir=work_dir, threads=threads)
        
        logger.info("Initialized bed conversion analyzer")
        analysis_result["processing_steps"].append("analyzer_initialized")

        # Single pass: collect existing paths and warn on misses
        valid_bam_paths: List[str] = []
        for p in bam_paths:
            if os.path.exists(p):
                valid_bam_paths.append(p)
            else:
                logger.warning(f"BAM file not found: {p}")
        if not valid_bam_paths:
            analysis_result["error_message"] = "No BAM files found or all paths missing"
            analysis_result["processing_steps"].append("no_files_processed")
            return analysis_result

        logger.info(f"Processing {len(valid_bam_paths)} BAM files (up to {bed_analyzer.threads} in parallel)")
        all_processed_data = bed_analyzer._process_bams(valid_bam_paths, sample_dir)
        processed_files = len(all_processed_data)

        if processed_files == 0:
            analysis_result["error_message"] = "No files could be processed successfully"
            analysis_result["processing_steps"].append("no_files_processed")
            return analysis_result

        analysis_result["files_processed"] = processed_files
        analysis_result["processing_steps"].append("files_processed")

        # Create parquet file from all processed data
        if all_processed_data:
            logger.info(f"Creating parquet file from {len(all_processed_data)} processed files")
            try:
                # Use merge_modkit_files to create the final parquet file
                bed_analyzer._update_state(
                    parquet_path, 
                    all_processed_data, 
                    sample_id, 
                    processed_files  # Use number of processed files as file_number
                )
                
                analysis_result["processing_steps"].append("parquet_created")
                logger.info(f"Parquet file created successfully: {parquet_path}")
                
                # Verify parquet file was created
                if os.path.exists(parquet_path):
                    parquet_size = os.path.getsize(parquet_path)
                    logger.info(f"Parquet file size: {parquet_size} bytes")
                else:
                    msg = "Parquet file was not created"
                    if expected_fail_only:
                        logger.warning(f"{msg} (expected for fail-only BAM submission)")
                    else:
                        logger.error(msg)
                    analysis_result["error_message"] = "Parquet file creation failed"
                    return analysis_result
                    
            except Exception as e:
                logger.error(f"Error creating parquet file: {e}")
                analysis_result["error_message"] = f"Parquet creation failed: {str(e)}"
                return analysis_result
        else:
            analysis_result["error_message"] = "No processed data available for parquet creation"
            analysis_result["processing_steps"].append("no_data_for_parquet")
            return analysis_result

        analysis_result["processing_steps"].append("analysis_complete")
        logger.info(f"Multi-file bed conversion completed for {sample_id}")
        logger.info(f"Files successfully processed: {analysis_result['files_processed']}/{analysis_result['total_files']}")
        logger.info(f"Output parquet file: {analysis_result['parquet_path']}")
        
        return analysis_result

    except Exception as e:
        logger.error(f"Error in multi-file bed conversion for {sample_id}: {e}")
        analysis_result["error_message"] = str(e)
        analysis_result["processing_steps"].append("analysis_failed")
        return analysis_result


def bed_conversion_handler(job, work_dir=None):
    """
    Handler function for BAM to parquet conversion jobs.
    This function processes BAM files for parquet conversion analysis.

    Args:
        job: The workflow job containing file and metadata
        work_dir: Optional base directory for output (defaults to BAM file directory)
    """
    try:
        # Get job-specific logger
        logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
        suppress_expected = _is_fail_only_expected(job)
        
        # Check if this is a batched job
        batched_job = job.context.metadata.get("_batched_job")
        if batched_job:
            batch_size = batched_job.get_file_count()
            sample_id = batched_job.get_sample_id()
            batch_id = batched_job.batch_id
            logger.info(f"Processing bed conversion batch: {batch_size} files for sample '{sample_id}' (batch_id: {batch_id})")
            
            # Get all filepaths in the batch
            filepaths = batched_job.get_filepaths()
            
            # Log individual files in the batch
            for i, filepath in enumerate(filepaths):
                logger.info(f"  Batch file {i+1}/{batch_size}: {os.path.basename(filepath)}")
            
            # Prepare metadata list for all BAM files in the batch (list comp inlined in 3.12)
            def _batch_metadata(i: int) -> dict:
                ctx = batched_job.contexts[i]
                sid = ctx.get_sample_id()
                return {**ctx.metadata.get("bam_metadata", {}), "sample_id": sid if sid != "unknown" else sample_id}
            metadata_list = [_batch_metadata(i) for i in range(len(filepaths))]

            # Determine work directory for the batch
            if work_dir is None:
                # Default to first BAM file directory
                batch_work_dir = os.path.dirname(filepaths[0])
            else:
                # Use specified work directory, create if it doesn't exist
                os.makedirs(work_dir, exist_ok=True)
                batch_work_dir = work_dir
                logger.debug(f"Using specified work directory: {batch_work_dir}")
            
            # Process all BAM files in the batch using the new aggregated function
            logger.info(f"Processing {batch_size} BAM files as aggregated batch for sample '{sample_id}'")
            batch_result = process_multiple_files(
                bam_paths=filepaths,
                metadata_list=metadata_list,
                work_dir=batch_work_dir,
                logger=logger,
                threads=1,  # Default thread count
                expected_fail_only=bool(suppress_expected),
            )
            
            # Store batch results in job context (maintain compatibility with existing structure)
            job.context.add_metadata("bed_conversion", {
                "batch_result": batch_result,  # Single aggregated result
                "batch_size": batch_size,
                "sample_id": sample_id,
                "batch_id": batch_id,
                "files_processed": batch_result.get("files_processed", batch_size),
                "total_files": batch_result.get("total_files", batch_size)
            })
            
            logger.info(f"Completed bed conversion batch processing: {batch_size} files for sample '{sample_id}'")
            logger.info(f"Files successfully processed: {batch_result.get('files_processed', batch_size)}/{batch_result.get('total_files', batch_size)}")
            
            if batch_result.get("error_message"):
                if suppress_expected:
                    logger.warning(
                        f"Batch processing completed with expected errors (fail-only BAM submission): {batch_result['error_message']}"
                    )
                    job.context.add_result(
                        "bed_conversion",
                        {
                            "status": "expected_failure",
                            "error_message": batch_result["error_message"],
                            "sample_id": sample_id,
                            "parquet_path": batch_result.get("parquet_path", ""),
                        },
                    )
                else:
                    logger.error(
                        f"Batch processing completed with errors: {batch_result['error_message']}"
                    )
                    job.context.add_error("bed_conversion", batch_result["error_message"])
            else:
                logger.info("Batch processing completed successfully with aggregated bed conversion")
                job.context.add_result(
                    "bed_conversion",
                    {
                        "status": "success",
                        "sample_id": sample_id,
                        "analysis_time": batch_result.get("analysis_timestamp", 0),
                        "parquet_path": batch_result.get("parquet_path", ""),
                        "processing_steps": batch_result.get("processing_steps", []),
                        "files_processed": batch_result.get("files_processed", batch_size),
                        "total_files": batch_result.get("total_files", batch_size),
                    },
                )
            
            return
            
        else:
            # Single file processing (backward compatibility)
            bam_path = job.context.filepath

            logger.info(
                f"Starting BAM to parquet conversion for: {os.path.basename(bam_path)}"
            )

            # Get metadata from preprocessing
            bam_metadata = job.context.metadata.get("bam_metadata", {})

            # Determine work directory
            if work_dir is None:
                # Default to BAM file directory
                work_dir = os.path.dirname(bam_path)
            else:
                # Use specified work directory, create if it doesn't exist
                os.makedirs(work_dir, exist_ok=True)
                logger.debug(f"Using specified work directory: {work_dir}")

            # Create BAM to parquet conversion analysis instance
            bed_analyzer = BedConversionAnalysis(work_dir=work_dir)

            # Process the BAM file
            bed_result = bed_analyzer.process_bam_file(bam_path, bam_metadata)

            # Store results in job context
            job.context.add_metadata("bed_conversion", bed_result.results)
            job.context.add_metadata("bed_processing_steps", bed_result.processing_steps)

            if bed_result.error_message:
                if suppress_expected:
                    logger.warning(
                        f"Expected BAM->parquet issue for fail-only BAM submission: {bed_result.error_message}"
                    )
                    job.context.add_result(
                        "bed_conversion",
                        {
                            "status": "expected_failure",
                            "error_message": bed_result.error_message,
                            "sample_id": bed_result.sample_id,
                            "parquet_path": bed_result.parquet_path,
                        },
                    )
                else:
                    job.context.add_error("bed_conversion", bed_result.error_message)
                    logger.error(
                        f"BAM to parquet conversion failed: {bed_result.error_message}"
                    )
            else:
                job.context.add_result(
                    "bed_conversion",
                    {
                        "status": "success",
                        "sample_id": bed_result.sample_id,
                        "analysis_time": bed_result.analysis_timestamp,
                        "parquet_path": bed_result.parquet_path,
                        "processing_steps": bed_result.processing_steps,
                    },
                )
                logger.info(
                    f"BAM to parquet conversion complete for {os.path.basename(bam_path)}"
                )
                logger.info(f"Sample ID: {bed_result.sample_id}")
                logger.info(f"Parquet file: {bed_result.parquet_path}")
                logger.debug(f"Processing steps: {', '.join(bed_result.processing_steps)}")

    except Exception as e:
        suppress_expected = _is_fail_only_expected(job)
        if suppress_expected:
            logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
            logger.warning(
                f"Expected bed conversion failure for fail-only BAM submission: {e}"
            )
            job.context.add_result(
                "bed_conversion",
                {"status": "expected_failure", "error_message": str(e)},
            )
            return
        job.context.add_error("bed_conversion", str(e))
        logger = logging.getLogger("robin.analysis.bed_conversion")
        logger.error(
            f"Error in BAM to parquet conversion for {job.context.filepath}: {e}"
        )
