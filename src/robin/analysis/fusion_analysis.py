#!/usr/bin/env python3
"""
Simplified Fusion Analysis Module for robin.

Requires Python 3.12+ (PEP 709 inlined comprehensions, match/case type patterns,
faster interpreter). Workflow integration and standalone fusion detection;
core logic in fusion_work.py.

Key functions:
- Workflow integration (fusion_handler)
- Standalone file processing (process_single_file)
- Metadata management (FusionMetadata)
"""
from __future__ import annotations

import os
import sys

if sys.version_info < (3, 12):
    raise RuntimeError("robin fusion analysis requires Python 3.12 or newer")
import tempfile
import logging
import time
import json
from typing import Dict, Any, Optional, List, Tuple, Set
from dataclasses import dataclass
import numpy as np
import pandas as pd
import pysam
from robin.logging_config import get_job_logger

# Import core fusion detection logic from fusion_work.py
from robin.analysis.fusion_work import (
    process_bam_file,
    _generate_output_files,
    FusionMetadata,
    GeneRegion,
    _setup_file_paths,
    _load_bed_regions,
    _ensure_gene_regions_loaded,
    has_supplementary_alignments,
    find_reads_with_supplementary,
    _find_gene_intersections,
    _process_reads_for_fusions,
    _optimize_fusion_dataframe,
    _filter_fusion_candidates,
    process_bam_for_fusions_work,
    _merge_fusion_metadata_objects,
    _load_fusion_metadata,
    preprocess_fusion_data_standalone,
    process_bam_with_staging,
    accumulate_fusion_candidates,
    finalize_fusion_accumulation_for_sample,
)

logger = logging.getLogger(__name__)


def json_serializable(obj):
    """Convert numpy/pandas types to JSON-serializable (3.12: type patterns, no guards)."""
    match obj:
        case int():
            return int(obj)
        case float():
            return float(obj)
        case _ if isinstance(obj, np.ndarray):
            return obj.tolist()
        case _ if isinstance(obj, pd.Series):
            return obj.tolist()
        case _ if isinstance(obj, pd.DataFrame):
            return obj.to_dict("records")
        case _:
            return obj


def _get_fusion_batch_size(default: int = 10) -> int:
    """
    Get the workflow batch size for fusion jobs.
    Falls back to a reasonable default if workflow config is unavailable.
    """
    try:
        from robin.workflow_ray import BATCH_CONFIG

        return int(BATCH_CONFIG.get("fusion", {}).get("max_batch_size", default))
    except Exception:
        return default


def _merge_fusion_metadata(
    new_metadata: FusionMetadata, existing_metadata: Optional[Dict]
) -> FusionMetadata:
    """
    Merge new fusion metadata with existing metadata.

    Args:
        new_metadata: New fusion metadata object
        existing_metadata: Existing metadata dictionary

    Returns:
        Merged FusionMetadata object
    """
    if existing_metadata is None:
        return new_metadata

    # Merge processing steps
    existing_steps = existing_metadata.get("processing_steps", [])
    new_steps = new_metadata.processing_steps
    merged_steps = list(set(existing_steps + new_steps))  # Remove duplicates

    # Merge fusion data
    existing_fusion_data = existing_metadata.get("fusion_data", {})
    new_fusion_data = new_metadata.fusion_data or {}

    # Combine target candidates
    existing_target = existing_fusion_data.get("target_candidates", [])
    new_target = new_fusion_data.get("target_candidates", [])
    merged_target = existing_target + new_target

    # Combine genome-wide candidates
    existing_genome = existing_fusion_data.get("genome_wide_candidates", [])
    new_genome = new_fusion_data.get("genome_wide_candidates", [])
    merged_genome = existing_genome + new_genome

    # Update the new metadata object
    new_metadata.processing_steps = merged_steps
    new_metadata.fusion_data = {
        "target_candidates": merged_target,
        "genome_wide_candidates": merged_genome,
    }

    # Update analysis results counts
    if new_metadata.analysis_results:
        new_metadata.analysis_results["target_candidates_count"] = len(merged_target)
        new_metadata.analysis_results["genome_wide_candidates_count"] = len(
            merged_genome
        )

    logging.info(
        f"Merged fusion metadata: {len(existing_steps)} existing steps + {len(new_steps)} new steps"
    )

    return new_metadata


def _perform_fusion_analysis(
    file_path: str,
    temp_dir: str,
    metadata: Dict[str, Any],
    fusion_metadata: FusionMetadata,
    target_panel: str,
    has_supplementary: bool = False,
    supplementary_read_ids: List[str] = [],
    work_dir: str = None,
) -> Dict[str, Any]:
    """
    Perform the actual fusion analysis using fusion_work.py.

    Args:
        file_path: Path to the file to process
        temp_dir: Temporary directory for intermediate files
        metadata: File metadata
        fusion_metadata: Fusion metadata object
        target_panel: Target panel to use
        has_supplementary: Whether file has supplementary reads
        supplementary_read_ids: List of supplementary read IDs
        work_dir: Working directory for saving metadata

    Returns:
        Dictionary with analysis results
    """
    fusion_metadata.processing_steps.append("analysis_started")
    results = process_bam_file(
        file_path,
        temp_dir,
        metadata,
        fusion_metadata,
        target_panel,
        has_supplementary,
        supplementary_read_ids,
        work_dir,
    )
    return results


def process_single_file(
    file_path: str,
    metadata: Dict[str, Any],
    work_dir: str,
    target_panel: str,
    logger=None,
    has_supplementary: bool = False,
    supplementary_read_ids: List[str] = [],
) -> Dict[str, Any]:
    """
    Process a single file for fusion analysis (standalone function-based approach).

    Args:
        file_path: Path to the file to process
        metadata: Metadata about the file
        work_dir: Working directory for output
        logger: Logger instance
        target_panel: Target panel to use for analysis
        has_supplementary: Whether file has supplementary reads
        supplementary_read_ids: List of supplementary read IDs

    Returns:
        Dictionary with processing results
    """
    try:
        if logger:
            logger.info(f"Starting fusion analysis for {file_path}")
            logger.info(f"Using work directory: {work_dir}")
            logger.info(f"Target panel: {target_panel}")

        # Extract sample ID from metadata
        sample_id = metadata.get("sample_id", "unknown")

        # Create metadata object
        fusion_metadata = FusionMetadata(
            sample_id=sample_id,
            file_path=file_path,
            analysis_timestamp=time.time(),
            target_panel=target_panel,
        )

        fusion_metadata.processing_steps.append("started")

        # Check if file exists
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")

        # Check if file is readable
        if not os.access(file_path, os.R_OK):
            raise PermissionError(f"File not readable: {file_path}")

        # Check if file is a BAM file
        if not file_path.lower().endswith(".bam"):
            if logger:
                logger.warning(f"File does not have .bam extension: {file_path}")

        # Create temporary directory
        with tempfile.TemporaryDirectory() as temp_dir:
            fusion_metadata.processing_steps.append("temp_dir_created")

            # Perform fusion analysis using fusion_work.py
            analysis_results = _perform_fusion_analysis(
                file_path,
                temp_dir,
                metadata,
                fusion_metadata,
                target_panel,
                has_supplementary,
                supplementary_read_ids,
                work_dir,
            )

            # Try to get reference from environment (fallback - handlers should pass it)
            reference = None
            try:
                reference = os.environ.get("robin_REFERENCE")
                if reference:
                    reference = os.path.expanduser(reference)
            except Exception:
                pass
            
            # Generate output files using fusion_work.py
            # For single file processing, don't generate master BED (it should be done at batch end)
            output_files = _generate_output_files(
                sample_id, analysis_results, fusion_metadata, work_dir, reference=reference, generate_master_bed=False
            )

            # Update metadata with results (now includes merged data)
            fusion_metadata.analysis_results = analysis_results
            fusion_metadata.target_fusion_path = output_files.get("target_fusion")
            fusion_metadata.genome_wide_fusion_path = output_files.get(
                "genome_wide_fusion"
            )

            # Update analysis results with final merged counts
            if fusion_metadata.fusion_data:
                fusion_metadata.analysis_results["target_candidates_count"] = len(
                    fusion_metadata.fusion_data.get("target_candidates", [])
                )
                fusion_metadata.analysis_results["genome_wide_candidates_count"] = len(
                    fusion_metadata.fusion_data.get("genome_wide_candidates", [])
                )

            fusion_metadata.processing_steps.append("completed")

            if logger:
                logger.info(f"Fusion analysis completed for {file_path}")

            # Convert to dictionary
            result = {
                "success": fusion_metadata.error_message is None,
                "sample_id": fusion_metadata.sample_id,
                "file_path": fusion_metadata.file_path,
                "analysis_timestamp": fusion_metadata.analysis_timestamp,
                "target_fusion_path": fusion_metadata.target_fusion_path,
                "genome_wide_fusion_path": fusion_metadata.genome_wide_fusion_path,
                "processing_steps": fusion_metadata.processing_steps,
                "error_message": fusion_metadata.error_message,
                "analysis_results": fusion_metadata.analysis_results,
            }

            if fusion_metadata.error_message:
                if logger:
                    logger.error(
                        f"Fusion analysis failed: {fusion_metadata.error_message}"
                    )
                    logger.error(
                        f"Processing steps completed: {fusion_metadata.processing_steps}"
                    )
            else:
                if logger:
                    logger.info(
                        f"Fusion analysis completed successfully for {file_path}"
                    )
                    logger.info(f"Processing steps: {fusion_metadata.processing_steps}")

            return result

    except Exception as e:
        import traceback

        error_msg = f"Error in fusion analysis: {str(e)}"
        if logger:
            logger.error(error_msg)
            logger.error(f"Traceback: {traceback.format_exc()}")
        return {
            "success": False,
            "error_message": error_msg,
            "file_path": file_path,
        }


def process_multiple_files(bam_paths, metadata_list, work_dir, logger, target_panel=None, reference=None):
    """
    Process multiple BAM files for fusion analysis using staged processing.
    
    This function processes multiple BAM files for the same sample using the
    existing staging infrastructure. Each file is processed individually and
    staged, then all staged files are accumulated in a single batch operation.
    Only files with supplementary reads are processed.

    Args:
        bam_paths: List of paths to BAM files
        metadata_list: List of metadata dictionaries (one per BAM file)
        work_dir: Working directory
        logger: Logger instance
        target_panel: Target panel type (rCNS2, AML, or custom panel)
        reference: Optional path to reference genome

    Returns:
        Dictionary with aggregated fusion analysis results
    """
    if not bam_paths or not metadata_list:
        raise ValueError("bam_paths and metadata_list must not be empty")
    
    if len(bam_paths) != len(metadata_list):
        raise ValueError("bam_paths and metadata_list must have the same length")
    
    # Get sample ID from first metadata (assuming all BAMs are from same sample)
    sample_id = metadata_list[0].get("sample_id", "unknown")
    
    logger.info(f"🔗 Starting multi-file fusion analysis for sample: {sample_id}")
    logger.info(f"Processing {len(bam_paths)} BAM files for sample {sample_id}")
    
    # Log essential metadata only (if debug logging is enabled)
    # JobLogger wraps a standard logger, access it via .logger attribute
    if hasattr(logger, 'logger') and logger.logger.isEnabledFor(logging.DEBUG):
        for i, (bam_path, metadata) in enumerate(zip(bam_paths, metadata_list)):
            logger.debug(f"BAM file {i+1}: {os.path.basename(bam_path)}")

    analysis_result = {
        "sample_id": sample_id,
        "bam_paths": bam_paths,
        "analysis_timestamp": time.time(),
        "processing_steps": [],
        "error_message": None,
        "files_processed": 0,
        "files_with_supplementary": 0,
        "total_files": len(bam_paths),
        "fusion_data": {},
        "target_fusion_path": None,
        "genome_wide_fusion_path": None,
    }

    try:
        # Filter files that have supplementary reads
        valid_bam_paths = []
        valid_metadata_list = []
        
        for bam_path, metadata in zip(bam_paths, metadata_list):
            has_supplementary = metadata.get("has_supplementary_reads", False)
            if has_supplementary:
                valid_bam_paths.append(bam_path)
                valid_metadata_list.append(metadata)
                analysis_result["files_with_supplementary"] += 1
                logger.debug(f"BAM {os.path.basename(bam_path)}: has supplementary reads")
            else:
                logger.debug(f"BAM {os.path.basename(bam_path)}: no supplementary reads - skipping")
        
        if not valid_bam_paths:
            logger.info(f"No BAM files with supplementary reads found for {sample_id}")
            analysis_result["error_message"] = "No supplementary reads found in any BAM files"
            analysis_result["processing_steps"].append("no_supplementary_reads")
            return analysis_result

        logger.info(f"Processing {len(valid_bam_paths)} valid BAM files out of {len(bam_paths)} total")
        analysis_result["files_processed"] = len(valid_bam_paths)
        analysis_result["processing_steps"].append("supplementary_reads_found")

        # Pre-initialize batch-level resources (done once per batch, not per file)
        # 1. Ensure gene regions are loaded (cached, but ensure it's done once)
        from robin.analysis.fusion_work import _ensure_gene_regions_loaded
        _ensure_gene_regions_loaded(target_panel)
        
        # 2. Pre-create staging directory (avoids repeated os.makedirs calls)
        from robin.analysis.fusion_work import _get_staging_dir
        staging_dir = _get_staging_dir(work_dir, sample_id)
        logger.debug(f"Pre-created staging directory: {staging_dir}")
        
        # Process each BAM file individually using staging
        logger.info("Processing files with fusion staging (fast path)")
        processed_files = 0
        # Collect all supplementary read IDs across all BAM files for consolidated output
        all_supplementary_read_ids = set()
        staging_batch_size = _get_fusion_batch_size()
        
        for i, (bam_path, metadata) in enumerate(zip(valid_bam_paths, valid_metadata_list)):
            logger.info(f"Processing BAM file {i+1}/{len(valid_bam_paths)}: {os.path.basename(bam_path)}")
            
            try:
                # Get supplementary read information
                has_supplementary = metadata.get("has_supplementary_reads", False)
                supplementary_read_ids = metadata.get("supplementary_read_ids", [])
                supp_ids_path = metadata.get("supplementary_read_ids_path")
                
                # Load supplementary read IDs from file if available
                if (not supplementary_read_ids) and supp_ids_path and os.path.exists(supp_ids_path):
                    try:
                        with open(supp_ids_path, "r") as f:
                            supplementary_read_ids = [
                                line.strip() for line in f if line.strip()
                            ]
                    except Exception as e:
                        logger.warning(f"Could not read supplementary_read_ids from {supp_ids_path}: {e}")
                
                # Collect supplementary read IDs for consolidated output
                if supplementary_read_ids:
                    all_supplementary_read_ids.update(supplementary_read_ids)
                
                # Create fusion metadata
                fusion_metadata = FusionMetadata(
                    sample_id=sample_id,
                    file_path=bam_path,
                    analysis_timestamp=time.time(),
                    target_panel=target_panel,
                )
                
                # Create temporary directory for processing
                with tempfile.TemporaryDirectory() as temp_dir:
                    # Process with staging
                    # Note: should_accumulate is ignored here since we force accumulation at batch end
                    analysis_results, _ = process_bam_with_staging(
                        bam_path,
                        temp_dir,
                        metadata,
                        fusion_metadata,
                        target_panel,
                        has_supplementary=has_supplementary,
                        supplementary_read_ids=supplementary_read_ids,
                        work_dir=work_dir,
                        batch_size=staging_batch_size,
                    )
                    
                    if analysis_results.get("error_message"):
                        logger.warning(f"Error processing {os.path.basename(bam_path)}: {analysis_results['error_message']}")
                        continue
                    
                    processed_files += 1
                    logger.debug(f"Successfully staged file {i+1}: {os.path.basename(bam_path)}")
                    
                    # Keep supplementary IDs file for later searching/debugging
                    # File is preserved at: {work_dir}/{sample_id}/supplementary_read_ids.txt
                
            except Exception as e:
                logger.warning(f"Error processing {os.path.basename(bam_path)}: {e}")
                continue

        if processed_files == 0:
            analysis_result["error_message"] = "No files could be processed successfully"
            analysis_result["processing_steps"].append("no_files_processed")
            return analysis_result

        analysis_result["files_processed"] = processed_files
        analysis_result["processing_steps"].append("files_staged")

        # Write consolidated supplementary read IDs file (aggregates all BAM files)
        if all_supplementary_read_ids:
            try:
                sample_dir = os.path.join(work_dir, sample_id)
                os.makedirs(sample_dir, exist_ok=True)
                supp_output_path = os.path.join(sample_dir, "supplementary_read_ids.txt")
                # Write one ID per line (sorted for easier searching)
                tmp_path = supp_output_path + ".tmp"
                with open(tmp_path, "w") as f:
                    for rid in sorted(all_supplementary_read_ids):
                        f.write(f"{rid}\n")
                os.replace(tmp_path, supp_output_path)
                logger.info(f"Saved {len(all_supplementary_read_ids)} unique supplementary read IDs to {supp_output_path}")
            except Exception as e:
                logger.warning(f"Could not write consolidated supplementary_read_ids file: {e}")

        # Force accumulation of all staged files
        # Expand reference path if provided
        if reference:
            reference = os.path.expanduser(reference)
        
        logger.info(f"Accumulating {processed_files} staged fusion files for sample {sample_id}")
        accumulation_result = accumulate_fusion_candidates(
            work_dir, sample_id, target_panel, force=True, batch_size=1, reference=reference
        )
        
        if accumulation_result.get("status") != "success":
            analysis_result["error_message"] = f"Accumulation failed: {accumulation_result.get('error', 'Unknown error')}"
            analysis_result["processing_steps"].append("accumulation_failed")
            return analysis_result

        analysis_result["processing_steps"].append("accumulation_complete")
        logger.info(f"Fusion accumulation completed: {accumulation_result}")
        
        # Generate final master BED file now that all files are processed
        try:
            from robin.analysis.master_bed_generator import generate_master_bed
            from robin.analysis.fusion_work import _load_analysis_counter
            
            # Get analysis counter
            analysis_counter = _load_analysis_counter(sample_id, work_dir)
            
            logger.info(f"Generating final master BED file for sample {sample_id} (counter: {analysis_counter})")
            master_bed_path = generate_master_bed(
                sample_id=sample_id,
                work_dir=work_dir,
                analysis_counter=analysis_counter,
                target_panel=target_panel,
                logger_instance=logger,
                reference=reference,
            )
            if master_bed_path:
                logger.info(f"Final master BED file generated: {master_bed_path}")
                analysis_result["processing_steps"].append("master_bed_generated")
        except Exception as e:
            logger.warning(f"Could not generate final master BED file: {e}")

        # Load final accumulated data for result metadata
        sample_output_dir = os.path.join(work_dir, sample_id)
        
        # Set output file paths
        analysis_result["target_fusion_path"] = os.path.join(sample_output_dir, "target_fusion.csv")
        analysis_result["genome_wide_fusion_path"] = os.path.join(sample_output_dir, "genome_wide_fusion.csv")
        
        # Store final results
        analysis_result["fusion_data"] = {
            "target_candidates_count": accumulation_result.get("target_candidates_count", 0),
            "genome_wide_candidates_count": accumulation_result.get("genome_wide_candidates_count", 0),
            "files_with_supplementary": analysis_result["files_with_supplementary"],
            "files_processed": processed_files,
        }
        
        analysis_result["processing_steps"].append("analysis_complete")
        logger.info(f"Multi-file fusion analysis completed for {sample_id}")
        logger.info(f"Files successfully processed: {analysis_result['files_processed']}/{analysis_result['total_files']}")
        logger.info(f"Files with supplementary reads: {analysis_result['files_with_supplementary']}")
        logger.info(f"Target fusion candidates: {analysis_result['fusion_data']['target_candidates_count']}")
        logger.info(f"Genome-wide fusion candidates: {analysis_result['fusion_data']['genome_wide_candidates_count']}")
        
        return analysis_result

    except Exception as e:
        logger.error(f"Error in multi-file fusion analysis for {sample_id}: {e}")
        analysis_result["error_message"] = str(e)
        analysis_result["processing_steps"].append("analysis_failed")
        return analysis_result


def fusion_handler(job, work_dir=None, target_panel=None):
    """
    Handler function for fusion analysis jobs in the workflow system.

    Args:
        job: Job object from the workflow system
        work_dir: Working directory for output
        target_panel: Target gene panel for fusion analysis (required)
    """
    # Validate required parameters
    if not target_panel:
        raise ValueError("target_panel is required for fusion analysis")
    
    try:
        # Get logger with proper parameters
        logger = get_job_logger(str(job.job_id), "fusion", job.context.filepath)
        
        # Check if this is a batched job
        batched_job = job.context.metadata.get("_batched_job")
        if batched_job:
            batch_size = batched_job.get_file_count()
            sample_id = batched_job.get_sample_id()
            batch_id = batched_job.batch_id
            logger.info(f"Processing fusion analysis batch: {batch_size} files for sample '{sample_id}' (batch_id: {batch_id})")
            
            # Get all filepaths in the batch
            filepaths = batched_job.get_filepaths()
            
            # Log individual files in the batch
            for i, filepath in enumerate(filepaths):
                logger.info(f"  Batch file {i+1}/{batch_size}: {os.path.basename(filepath)}")
            
            # Prepare metadata list for all BAM files in the batch
            metadata_list = []
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
                
                metadata_list.append(file_metadata)
            
            # Determine work directory for the batch
            if work_dir is None:
                # Default to first BAM file directory
                batch_work_dir = os.path.dirname(filepaths[0])
            else:
                # Use specified work directory, create if it doesn't exist
                os.makedirs(work_dir, exist_ok=True)
                batch_work_dir = work_dir
                logger.debug(f"Using specified work directory: {batch_work_dir}")
            
            # Log and validate target panel
            job_panel = batched_job.contexts[0].metadata.get("target_panel", target_panel)
            if job_panel != target_panel:
                logger.warning(f"Panel mismatch: job metadata has '{job_panel}' but handler received '{target_panel}'. Using '{job_panel}' from metadata.")
                target_panel = job_panel
            logger.info(f"Using target panel: {target_panel}")
            
            # Get reference from job metadata
            reference = job.context.metadata.get("reference")
            if reference:
                reference = os.path.expanduser(reference)
                logger.debug(f"Using reference genome from job metadata: {reference}")
            
            # Process all BAM files in the batch using the new aggregated function
            logger.info(f"Processing {batch_size} BAM files as aggregated batch for sample '{sample_id}'")
            batch_result = process_multiple_files(
                bam_paths=filepaths,
                metadata_list=metadata_list,
                work_dir=batch_work_dir,
                logger=logger,
                target_panel=target_panel,
                reference=reference
            )
            
            # Store batch results in job context (maintain compatibility with existing structure)
            job.context.add_metadata("fusion_analysis", {
                "batch_result": batch_result,  # Single aggregated result
                "batch_size": batch_size,
                "sample_id": sample_id,
                "batch_id": batch_id,
                "files_processed": batch_result.get("files_processed", batch_size),
                "total_files": batch_result.get("total_files", batch_size)
            })
            
            logger.info(f"Completed fusion analysis batch processing: {batch_size} files for sample '{sample_id}'")
            logger.info(f"Files successfully processed: {batch_result.get('files_processed', batch_size)}/{batch_result.get('total_files', batch_size)}")
            logger.info(f"Files with supplementary reads: {batch_result.get('files_with_supplementary', 0)}")
            
            if batch_result.get("error_message"):
                logger.error(f"Batch processing completed with errors: {batch_result['error_message']}")
                job.context.add_error("fusion_analysis", batch_result["error_message"])
            else:
                logger.info("Batch processing completed successfully with aggregated fusion analysis")
                job.context.add_result(
                    "fusion_analysis",
                    {
                        "success": True,
                        "sample_id": sample_id,
                        "analysis_time": batch_result.get("analysis_timestamp", 0),
                        "target_candidates_count": batch_result.get("fusion_data", {}).get("target_candidates_count", 0),
                        "genome_wide_candidates_count": batch_result.get("fusion_data", {}).get("genome_wide_candidates_count", 0),
                        "processing_steps": batch_result.get("processing_steps", []),
                        "target_fusion_path": batch_result.get("target_fusion_path", ""),
                        "genome_wide_fusion_path": batch_result.get("genome_wide_fusion_path", ""),
                        "files_processed": batch_result.get("files_processed", batch_size),
                        "total_files": batch_result.get("total_files", batch_size),
                        "files_with_supplementary": batch_result.get("files_with_supplementary", 0),
                    },
                )
            
            return
            
        else:
            # Single file processing (backward compatibility)
            # Extract file path and metadata from job
            file_path = job.context.filepath
            metadata = job.context.metadata.get("bam_metadata", {})

            # Log and validate target panel
            job_panel = job.context.metadata.get("target_panel", target_panel)
            if job_panel != target_panel:
                logger.warning(f"Panel mismatch: job metadata has '{job_panel}' but handler received '{target_panel}'. Using '{job_panel}' from metadata.")
                target_panel = job_panel
            logger.info(f"Using target panel: {target_panel}")
            logger.info(f"DEBUG: Job metadata: {job.context.metadata}")

            # Access supplementary read information
            has_supplementary = metadata.get("has_supplementary_reads", False)
            if has_supplementary:
                logger.info(f"Starting fusion analysis for {file_path}")
                logger.info(f"Metadata: {metadata}")

                # Prefer disk-based supplementary read IDs if available to avoid large in-memory lists
                supplementary_read_ids = metadata.get("supplementary_read_ids", [])
                supp_ids_path = metadata.get("supplementary_read_ids_path")
                if (
                    (not supplementary_read_ids)
                    and supp_ids_path
                    and os.path.exists(supp_ids_path)
                ):
                    try:
                        with open(supp_ids_path, "r") as f:
                            supplementary_read_ids = [
                                line.strip() for line in f if line.strip()
                            ]
                    except Exception as e:
                        logger.warning(
                            f"Could not read supplementary_read_ids from {supp_ids_path}: {e}"
                        )

                # Set default work directory if not provided
                if work_dir is None:
                    work_dir = "fusion_output"

                logger.info(f"Using work directory: {work_dir}")
                logger.info(f"Using target panel: {target_panel}")

                # Use staging-based processing for performance
                logger.info("Using fusion staging-based processing (fast path)")
                staging_batch_size = _get_fusion_batch_size()
                
                # Create fusion metadata
                sample_id = metadata.get("sample_id", "unknown")
                fusion_metadata = FusionMetadata(
                    sample_id=sample_id,
                    file_path=file_path,
                    analysis_timestamp=time.time(),
                    target_panel=target_panel,
                )
                
                # Create temporary directory for processing
                import tempfile
                with tempfile.TemporaryDirectory() as temp_dir:
                    # Process with staging
                    analysis_results, should_accumulate = process_bam_with_staging(
                        file_path,
                        temp_dir,
                        metadata,
                        fusion_metadata,
                        target_panel,
                        has_supplementary=has_supplementary,
                        supplementary_read_ids=supplementary_read_ids,
                        work_dir=work_dir,
                        batch_size=staging_batch_size,
                    )
                    
                    # Convert to result format
                    result = {
                        "success": True,
                        "sample_id": sample_id,
                        "file_path": file_path,
                        "analysis_timestamp": time.time(),
                        "processing_steps": ["staging_complete"],
                        "error_message": None,
                        "analysis_results": analysis_results,
                    }
                    
                    # Add result to job context
                    job.context.add_result("fusion_analysis", result)
                    
                    # Trigger accumulation if threshold reached
                    # The accumulation function re-checks the staging threshold before proceeding.
                    if should_accumulate:
                        logger.info("Fusion accumulation threshold reached - attempting batch accumulation")
                        # Get reference from job metadata if available
                        reference = job.context.metadata.get("reference")
                        if reference:
                            reference = os.path.expanduser(reference)
                            logger.debug(f"Using reference genome from job metadata: {reference}")
                        
                        accumulation_result = accumulate_fusion_candidates(
                            work_dir, sample_id, target_panel, force=False, batch_size=staging_batch_size, reference=reference
                        )
                        if accumulation_result.get("status") == "below_threshold":
                            logger.debug("Accumulation skipped - below threshold")
                        else:
                            logger.info(f"Fusion accumulation result: {accumulation_result}")
                        job.context.add_metadata("fusion_accumulation_result", accumulation_result)
                    
                    # Store flag for potential end-of-queue accumulation
                    job.context.add_metadata("needs_final_fusion_accumulation", True)

                if result["success"]:
                    logger.info(f"Fusion analysis completed successfully for {file_path}")
                    logger.info(f"Results: {result}")
                    # Keep supplementary IDs file for later searching/debugging
                    # File is preserved at: {work_dir}/{sample_id}/supplementary_read_ids.txt
                else:
                    error_msg = result.get("error_message", "Unknown error")
                    logger.error(f"Fusion analysis failed for {file_path}: {error_msg}")
                    logger.error(f"Full result: {result}")
                    job.context.add_error("fusion_analysis", error_msg)
            else:
                logger.info(f"No supplementary reads found for {file_path}")
                # This is not an error - just a normal skip condition
                job.context.add_result(
                    "fusion_analysis",
                    {
                        "success": True,
                        "skipped": True,
                        "reason": "No supplementary reads found",
                        "file_path": file_path,
                    },
                )

    except Exception as e:
        import traceback

        error_msg = f"Error in fusion handler: {str(e)}"

        # Try to get a logger, but don't fail if we can't
        try:
            logger = get_job_logger(str(job.job_id), "fusion", job.context.filepath)
            logger.error(error_msg)
            logger.error(f"Traceback: {traceback.format_exc()}")
            logger.error(
                f"Job details - ID: {job.job_id}, Type: {job.job_type}, File: {job.context.filepath}"
            )
        except Exception as logger_error:
            # If we can't get a logger, print to stderr
            import sys

            print(f"FUSION HANDLER ERROR: {error_msg}", file=sys.stderr)
            print(f"LOGGER ERROR: {logger_error}", file=sys.stderr)
            print(f"TRACEBACK: {traceback.format_exc()}", file=sys.stderr)

        job.context.add_error("fusion_analysis", error_msg)


def _get_available_panels() -> List[str]:
    """Get list of available panels from resources directory."""
    panels = ["rCNS2", "AML"]  # Built-in panels
    
    try:
        from pathlib import Path
        # Look for the resources directory relative to this file
        current_file = Path(__file__)
        resources_dir = current_file.parent.parent.parent / "robin" / "resources"
        
        if resources_dir.exists():
            # Look for custom panels (files ending with _panel_name_uniq.bed)
            for bed_file in resources_dir.glob("*_panel_name_uniq.bed"):
                panel_name = bed_file.stem.replace("_panel_name_uniq", "")
                if panel_name not in panels:
                    panels.append(panel_name)
            
            panels.sort()
        
    except Exception:
        # Fallback to built-in panels only
        pass
    
    return panels


if __name__ == "__main__":
    """
    Standalone CLI for fusion analysis testing.
    
    Usage:
        # Analyze a single BAM file
        python -m robin.analysis.fusion_analysis file.bam --target-panel rCNS2 --work-dir output/
        
        # Analyze all BAM files in a folder
        python -m robin.analysis.fusion_analysis /path/to/bam/folder --target-panel rCNS2 --work-dir output/
        
        # Analyze specific BAM files
        python -m robin.analysis.fusion_analysis file1.bam file2.bam --target-panel rCNS2 --work-dir output/
    """
    import argparse
    import glob
    from pathlib import Path
    
    parser = argparse.ArgumentParser(
        description="Standalone fusion analysis for BAM files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze a single BAM file
  python -m robin.analysis.fusion_analysis sample.bam --target-panel rCNS2 --work-dir output/
  
  # Analyze all BAM files in a folder
  python -m robin.analysis.fusion_analysis /path/to/bam/folder --target-panel rCNS2 --work-dir output/
  
  # Analyze multiple specific BAM files
  python -m robin.analysis.fusion_analysis file1.bam file2.bam file3.bam --target-panel AML --work-dir output/
  
  # With reference genome
  python -m robin.analysis.fusion_analysis folder/ --target-panel rCNS2 --work-dir output/ --reference /path/to/reference.fa
        """
    )
    
    parser.add_argument(
        "input",
        nargs="+",
        help="BAM file(s) or folder containing BAM files to analyze"
    )
    
    parser.add_argument(
        "--target-panel",
        type=str,
        choices=_get_available_panels(),
        required=True,
        help=f"Target gene panel for fusion analysis. Available: {', '.join(_get_available_panels())}"
    )
    
    parser.add_argument(
        "--work-dir",
        type=str,
        default="fusion_output",
        help="Working directory for output files (default: fusion_output)"
    )
    
    parser.add_argument(
        "--reference",
        type=str,
        help="Path to reference genome FASTA file (optional, for master BED generation)"
    )
    
    parser.add_argument(
        "--sample-id",
        type=str,
        help="Sample ID to use (default: auto-detect from first BAM file name)"
    )
    
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    
    args = parser.parse_args()
    
    # Set up logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logger = logging.getLogger("fusion_analysis_standalone")
    
    # Collect BAM files
    bam_files = []
    for input_path in args.input:
        path = Path(input_path)
        if path.is_file():
            if path.suffix.lower() == ".bam":
                bam_files.append(str(path))
            else:
                logger.warning(f"Skipping non-BAM file: {input_path}")
        elif path.is_dir():
            # Find all BAM files in directory
            bam_files.extend(glob.glob(str(path / "*.bam")))
        else:
            logger.warning(f"Path does not exist: {input_path}")
    
    if not bam_files:
        logger.error("No BAM files found to process!")
        sys.exit(1)
    
    logger.info(f"Found {len(bam_files)} BAM file(s) to process")
    
    # Expand reference path if provided
    reference = None
    if args.reference:
        reference = os.path.expanduser(args.reference)
        if not os.path.exists(reference):
            logger.error(f"Reference genome file not found: {reference}")
            sys.exit(1)
    
    # Create work directory
    work_dir = os.path.expanduser(args.work_dir)
    os.makedirs(work_dir, exist_ok=True)
    logger.info(f"Using work directory: {work_dir}")
    
    # Prepare metadata for each BAM file
    metadata_list = []
    sample_id = args.sample_id
    
    for bam_path in bam_files:
        # Check for supplementary reads
        logger.info(f"Checking {os.path.basename(bam_path)} for supplementary reads...")
        has_supplementary = has_supplementary_alignments(bam_path)
        
        if not has_supplementary:
            logger.warning(f"No supplementary reads found in {os.path.basename(bam_path)} - will be skipped")
        
        # Get supplementary read IDs
        supplementary_read_ids = []
        if has_supplementary:
            try:
                supplementary_read_ids = list(find_reads_with_supplementary(bam_path))
                logger.info(f"Found {len(supplementary_read_ids)} reads with supplementary alignments")
            except Exception as e:
                logger.warning(f"Error finding supplementary reads: {e}")
                has_supplementary = False
        
        # Auto-detect sample ID from first file if not provided
        if sample_id is None:
            # Extract sample ID from filename (remove .bam extension and common prefixes)
            filename = os.path.basename(bam_path)
            sample_id = filename.replace(".bam", "").replace(".BAM", "")
            # Remove common prefixes
            for prefix in ("sample_", "Sample_", "SAMPLE_"):
                if sample_id.startswith(prefix):
                    sample_id = sample_id.removeprefix(prefix)
                    break
            logger.info(f"Auto-detected sample ID: {sample_id}")
        
        metadata = {
            "sample_id": sample_id,
            "has_supplementary_reads": has_supplementary,
            "supplementary_read_ids": supplementary_read_ids,
            "file_path": bam_path,
        }
        metadata_list.append(metadata)
    
    # Process all files
    logger.info(f"Starting fusion analysis for sample: {sample_id}")
    logger.info(f"Target panel: {args.target_panel}")
    
    try:
        result = process_multiple_files(
            bam_paths=bam_files,
            metadata_list=metadata_list,
            work_dir=work_dir,
            logger=logger,
            target_panel=args.target_panel,
            reference=reference
        )
        
        if result.get("error_message"):
            logger.error(f"Analysis failed: {result['error_message']}")
            sys.exit(1)
        
        # Print summary
        print("\n" + "="*60)
        print("FUSION ANALYSIS COMPLETE")
        print("="*60)
        print(f"Sample ID: {result.get('sample_id', 'unknown')}")
        print(f"Files processed: {result.get('files_processed', 0)}/{result.get('total_files', 0)}")
        print(f"Files with supplementary reads: {result.get('files_with_supplementary', 0)}")
        print(f"Target fusion candidates: {result.get('fusion_data', {}).get('target_candidates_count', 0)}")
        print(f"Genome-wide fusion candidates: {result.get('fusion_data', {}).get('genome_wide_candidates_count', 0)}")
        print(f"\nOutput files:")
        if result.get('target_fusion_path'):
            print(f"  - Target fusions: {result['target_fusion_path']}")
        if result.get('genome_wide_fusion_path'):
            print(f"  - Genome-wide fusions: {result['genome_wide_fusion_path']}")
        print("="*60)
        
    except KeyboardInterrupt:
        logger.info("Analysis interrupted by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Error during analysis: {e}", exc_info=True)
        sys.exit(1)