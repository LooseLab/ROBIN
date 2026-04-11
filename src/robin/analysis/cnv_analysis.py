#!/usr/bin/env python3
"""
Copy Number Variation (CNV) Analysis Module for robin

This module provides comprehensive CNV analysis following the exact architecture
specified in the CNV BAM Processing Pipeline Documentation.

Key Features
-----------
- Two-pass CNV analysis (sample vs reference)
- Dynamic bin width calculation
- Breakpoint detection using Kernel Change Point Detection
- State persistence and incremental processing
- Comprehensive output file generation
- Configurable execution mode: direct or subprocess (USE_CNV_SUBPROCESS flag)

Classes
-------
CNVMetadata
    Container for CNV analysis metadata and results.

CNVAnalysis
    Main analysis class that processes BAM files for CNV analysis.

Dependencies
-----------
- pandas: Data manipulation and analysis
- numpy: Numerical computations
- pysam: BAM file processing
- ruptures: Change point detection
- scipy: Scientific computing
- cnv_from_bam: CNV extraction from BAM files (required)
- logging: Logging for debugging and monitoring
- typing: Type hints
- tempfile: Temporary file creation
- pathlib: File system paths
- os: Operating system interface
- time: Time utilities
- json: JSON serialization
- pickle: Python object serialization
- gc: Garbage collection

Example Usage
-----------
.. code-block:: python

    from robin.analysis.cnv_analysis import CNVAnalysis

    # Initialize analysis
    cnv_analysis = CNVAnalysis(
        work_dir="output/",
        ref_cnv_dict_path="HG01280_control_new.pkl"
    )

    # Process BAM files
    cnv_analysis.process_bam("sample.bam")

Notes
-----
The module follows the robin framework patterns for:
- Integration with workflow system
- Worker process management
- State tracking and persistence
- Error handling and logging
- Output generation and file management

Requirements
-----------
- cnv_from_bam module must be available for CNV extraction
- Reference CNV data must be available in robin module resources

Authors
-------
Matt Loose
"""

import os
import json
import logging
import time
import pickle
import gc
import subprocess
import sys
from typing import Dict, Optional, List, Tuple
from dataclasses import dataclass

import numpy as np
import pysam
from scipy.ndimage import uniform_filter1d

import ruptures as rpt
from robin.logging_config import get_job_logger
import robin.resources as resources

os.environ["CI"] = "1"


# Global configuration for performance tuning
CHUNK_SIZE = 1000  # For processing large arrays in chunks
# Removed MAX_WORKERS - no longer using parallel processing

# Configuration flag for CNV analysis execution mode
# Set to True to use subprocess execution (original approach)
# Set to False to use direct execution (new approach, default)
# To revert to subprocess mode, change this to: USE_CNV_SUBPROCESS = True
USE_CNV_SUBPROCESS = False


# Global cache for reference CNV dict to avoid reloading for every sample
_ref_cnv_dict_cache = None
_ref_cnv_dict_path_cache = None

# Global cache for per-sample data to avoid reloading across BAM files for the same sample
_sample_cache = {}
_current_sample_id = None


def get_cached_ref_cnv_dict(ref_cnv_dict_path: str, logger) -> dict:
    """
    Get reference CNV dictionary with caching to avoid reloading for every sample.
    
    Args:
        ref_cnv_dict_path: Path to reference CNV pickle file
        logger: Logger instance
        
    Returns:
        Reference CNV dictionary
    """
    global _ref_cnv_dict_cache, _ref_cnv_dict_path_cache
    
    # Check if we need to load the reference dict
    if _ref_cnv_dict_cache is None or _ref_cnv_dict_path_cache != ref_cnv_dict_path:
        logger.info(f"Loading reference CNV dict from {ref_cnv_dict_path} (first time or path changed)")
        load_start = time.time()
        
        with open(ref_cnv_dict_path, "rb") as f:
            _ref_cnv_dict_cache = pickle.load(f)
        _ref_cnv_dict_path_cache = ref_cnv_dict_path
        
        load_time = time.time() - load_start
        logger.info(f"Reference CNV dict loaded in {load_time:.3f}s (cached for subsequent samples)")
    else:
        logger.debug("Using cached reference CNV dict")
    
    return _ref_cnv_dict_cache

def clear_ref_cnv_dict_cache():
    """Clear the reference CNV dict cache (useful when switching reference files)"""
    global _ref_cnv_dict_cache, _ref_cnv_dict_path_cache
    _ref_cnv_dict_cache = None
    _ref_cnv_dict_path_cache = None

def get_sample_cache(sample_id: str) -> dict:
    """
    Get the cache for a specific sample, creating it if it doesn't exist.
    
    Args:
        sample_id: Sample identifier
        
    Returns:
        Dictionary containing cached data for the sample
    """
    global _sample_cache
    if sample_id not in _sample_cache:
        _sample_cache[sample_id] = {
            'copy_numbers': None,
            'copy_numbers_path': None,
            'analysis_counter': None,
            'last_accessed': time.time(),
            'bam_count': 0
        }
    _sample_cache[sample_id]['last_accessed'] = time.time()
    return _sample_cache[sample_id]

def update_sample_cache(sample_id: str, **kwargs):
    """
    Update the cache for a specific sample.
    
    Args:
        sample_id: Sample identifier
        **kwargs: Key-value pairs to update in the cache
    """
    cache = get_sample_cache(sample_id)
    cache.update(kwargs)

def clear_sample_cache(sample_id: str = None):
    """
    Clear the cache for a specific sample or all samples.
    
    Args:
        sample_id: Sample identifier to clear, or None to clear all
    """
    global _sample_cache, _current_sample_id
    if sample_id is None:
        _sample_cache.clear()
        _current_sample_id = None
    else:
        _sample_cache.pop(sample_id, None)
        if _current_sample_id == sample_id:
            _current_sample_id = None

def set_current_sample(sample_id: str, logger) -> bool:
    """
    Set the current sample being processed and manage cache transitions.
    
    Args:
        sample_id: Sample identifier
        logger: Logger instance
        
    Returns:
        True if this is a new sample (cache was cleared), False if same sample
    """
    global _current_sample_id
    
    if _current_sample_id != sample_id:
        logger.info(f"Switching from sample '{_current_sample_id}' to '{sample_id}' - managing cache")
        _current_sample_id = sample_id
        
        # Update BAM count for the new sample
        cache = get_sample_cache(sample_id)
        cache['bam_count'] += 1
        
        return True  # New sample
    else:
        # Same sample, just increment BAM count
        cache = get_sample_cache(sample_id)
        cache['bam_count'] += 1
        logger.debug(f"Processing BAM file {cache['bam_count']} for sample '{sample_id}'")
        return False  # Same sample

def cleanup_sample_cache_on_completion(sample_id: str, logger) -> None:
    """
    Clean up sample cache when all BAM files for a sample are complete.
    This should be called by the workflow system when a sample is fully processed.
    
    Args:
        sample_id: Sample identifier
        logger: Logger instance
    """
    cache = get_sample_cache(sample_id)
    bam_count = cache.get('bam_count', 0)
    logger.info(f"Sample '{sample_id}' completed processing {bam_count} BAM files - cache can be cleaned up")
    
    # Note: We don't actually clear the cache here as it might be needed for other operations
    # The workflow system should call clear_sample_cache() when appropriate

def get_sample_cache_stats() -> dict:
    """
    Get statistics about the current sample cache for monitoring purposes.
    
    Returns:
        Dictionary with cache statistics
    """
    global _sample_cache, _current_sample_id
    return {
        'current_sample': _current_sample_id,
        'cached_samples': list(_sample_cache.keys()),
        'cache_size': len(_sample_cache),
        'sample_details': {
            sid: {
                'bam_count': cache.get('bam_count', 0),
                'last_accessed': cache.get('last_accessed', 0),
                'has_copy_numbers': cache.get('copy_numbers') is not None,
                'analysis_counter': cache.get('analysis_counter')
            }
            for sid, cache in _sample_cache.items()
        }
    }

def run_cnv_analysis_direct(
    bam_path,
    copy_numbers,
    ref_cnv_dict,
    logger,
    threads=1,
    mapq_filter=60,
    sample_id: str = None,
):
    """
    Run CNV analysis directly using cnv_from_bam without subprocess isolation.
    OPTIMIZED: Accepts data directly to eliminate unnecessary file I/O.

    Args:
        bam_path: Path to BAM file
        copy_numbers: Copy numbers dictionary (passed directly, no file I/O)
        ref_cnv_dict: Reference CNV dictionary (passed directly, no file I/O)
        logger: Logger instance
        threads: Number of threads to use
        mapq_filter: Mapping quality filter
        sample_id: Sample identifier

    Returns:
        Dictionary with analysis results or None if failed
    """
    import cnv_from_bam
    
    try:
        # Use provided copy_numbers dict directly (no file I/O)
        if copy_numbers is None:
            copy_numbers = {}
            logger.debug(f"No existing copy_numbers, starting fresh")
        else:
            logger.debug(f"Using provided copy_numbers dict directly")

        # Use provided reference CNV dict directly (no file I/O)
        ref_cnv_dict_loaded = ref_cnv_dict
        logger.debug(f"Using provided reference CNV dict directly")

        # First pass: process sample with accumulated copy numbers
        logger.debug(f"Starting Pass 1: Sample CNV extraction with {threads} threads")
        pass1_start = time.time()
        result = cnv_from_bam.iterate_bam_file(
            bam_path,
            _threads=threads,
            mapq_filter=mapq_filter,
            copy_numbers=copy_numbers,
            log_level=int(logging.ERROR),
        )
        pass1_time = time.time() - pass1_start
        logger.info(f"Pass 1 completed in {pass1_time:.2f}s (bin_width: {result.bin_width}, variance: {result.variance:.6f})")

        # Second pass: process against reference using the same bin width
        logger.debug(f"Starting Pass 2: Reference CNV extraction with {threads} threads")
        pass2_start = time.time()
        result2 = cnv_from_bam.iterate_bam_file(
            bam_path,
            _threads=threads,
            mapq_filter=mapq_filter,
            copy_numbers=ref_cnv_dict_loaded,
            log_level=int(logging.ERROR),
            bin_width=result.bin_width,  # Use the same bin width as the sample
        )
        pass2_time = time.time() - pass2_start
        logger.info(f"Pass 2 completed in {pass2_time:.2f}s")
        logger.info(f"Total CNV extraction time: {pass1_time + pass2_time:.2f}s")

        # Prepare results
        analysis_results = {
            "success": True,
            "r_cnv": result.cnv,
            "r_bin": result.bin_width,
            "r_var": result.variance,
            "genome_length": result.genome_length,
            "r2_cnv": result2.cnv,
            "updated_copy_numbers": copy_numbers,  # The mutated copy_numbers
            "timing": {
                "pass1_time": pass1_time,
                "pass2_time": pass2_time,
                "total_time": pass1_time + pass2_time,
            }
        }

        logger.debug("CNV analysis completed successfully (direct mode)")
        return analysis_results

    except Exception as e:
        logger.error(f"Error in direct CNV analysis: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return {
            "success": False,
            "error": str(e),
            "error_type": type(e).__name__,
        }

def run_cnv_analysis_subprocess(
    bam_path,
    copy_numbers,
    ref_cnv_dict,
    temp_dir,
    logger,
    threads=1,
    mapq_filter=60,
    sample_id: str = None,
    copy_numbers_path: str = None,
    timeout: int = 3600,
):
    """
    Run CNV analysis using cnv_from_bam in a subprocess to avoid signal handling issues.

    Args:
        bam_path: Path to BAM file
        copy_numbers: Copy numbers dictionary (used if copy_numbers_path not provided)
        ref_cnv_dict: Reference CNV dictionary or path to reference pickle
        temp_dir: Temporary directory for intermediate files
        logger: Logger instance
        threads: Number of threads to use
        mapq_filter: Mapping quality filter
        sample_id: Sample identifier
        copy_numbers_path: Path to per-sample copy_numbers file (OPTIMIZED APPROACH)
        timeout: Subprocess timeout in seconds (default: 3600 = 1 hour)

    Returns:
        Dictionary with analysis results or None if failed
    """
    # Ensure temp directory exists
    os.makedirs(temp_dir, exist_ok=True)

    # Use cached reference CNV dict to avoid reloading for every sample
    if isinstance(ref_cnv_dict, str) and os.path.exists(ref_cnv_dict):
        # Load and cache the reference dict
        ref_cnv_dict_loaded = get_cached_ref_cnv_dict(ref_cnv_dict, logger)
        ref_cnv_dict_path = os.path.join(temp_dir, "ref_cnv_dict.pkl")
        with open(ref_cnv_dict_path, "wb") as f:
            pickle.dump(ref_cnv_dict_loaded, f, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        # Already a dict, serialize it
        ref_cnv_dict_path = os.path.join(temp_dir, "ref_cnv_dict.pkl")
        with open(ref_cnv_dict_path, "wb") as f:
            pickle.dump(ref_cnv_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

    # Get path to the subprocess script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    subprocess_script = os.path.join(script_dir, "cnv_subprocess.py")

    # Build command using per-sample copy_numbers file (OPTIMIZED)
    cmd = [
        sys.executable,  # Use the same Python interpreter
        subprocess_script,
        "--bam-path",
        bam_path,
        "--ref-cnv-dict-path",
        ref_cnv_dict_path,
        "--output-dir",
        temp_dir,
        "--threads",
        str(threads),
        "--mapq-filter",
        str(mapq_filter),
    ]
    
    # Use per-sample copy_numbers file if provided (OPTIMIZED APPROACH)
    if copy_numbers_path:
        cmd.extend(["--copy-numbers-path", copy_numbers_path])
        logger.debug(f"Using per-sample copy_numbers file: {copy_numbers_path}")
    else:
        # Fallback: serialize copy_numbers to temp file
        temp_copy_numbers_path = os.path.join(temp_dir, "temp_copy_numbers.pkl")
        with open(temp_copy_numbers_path, "wb") as f:
            pickle.dump(copy_numbers, f, protocol=pickle.HIGHEST_PROTOCOL)
        cmd.extend(["--copy-numbers-path", temp_copy_numbers_path])
        logger.debug("Using temporary copy_numbers file (fallback)")

    logger.debug(f"Running CNV analysis in subprocess: {' '.join(cmd)}")
    
    # Log file size for performance monitoring
    try:
        bam_size_mb = os.path.getsize(bam_path) / (1024 * 1024)
        logger.info(f"Processing BAM file: {bam_size_mb:.1f} MB")
        
        # Warn about very large files that might take a long time
        if bam_size_mb > 500:
            logger.warning(
                f"Large BAM file detected ({bam_size_mb:.1f} MB) - CNV analysis may take several minutes. "
                f"Timeout set to {timeout}s."
            )
    except Exception:
        pass

    # Run subprocess with stdout/stderr redirected to files to avoid capture overhead
    stdout_log_path = os.path.join(temp_dir, "cnv_subprocess.stdout.log")
    stderr_log_path = os.path.join(temp_dir, "cnv_subprocess.stderr.log")
    
    try:
        with open(stdout_log_path, "w") as _out, open(stderr_log_path, "w") as _err:
            result = subprocess.run(
                cmd,
                stdout=_out,
                stderr=_err,
                timeout=timeout,  # Configurable timeout
            )

        if result.returncode != 0:
            logger.error(f"Subprocess failed with return code {result.returncode}")
            try:
                with open(stdout_log_path, "r") as f:
                    stdout_content = f.read()
                if stdout_content.strip():
                    logger.error(f"subprocess stdout:\n{stdout_content}")
            except Exception:
                pass
            try:
                with open(stderr_log_path, "r") as f:
                    stderr_content = f.read()
                if stderr_content.strip():
                    logger.error(f"subprocess stderr:\n{stderr_content}")
            except Exception:
                pass
            return None

        # Load results
        results_path = os.path.join(temp_dir, "cnv_analysis_results.pkl")
        if os.path.exists(results_path):
            with open(results_path, "rb") as f:
                result = pickle.load(f)
            return result
        else:
            logger.error("Results file not found")
            return None
    
    except subprocess.TimeoutExpired:
        logger.error(
            f"CNV analysis subprocess timed out after {timeout}s for {os.path.basename(bam_path)}. "
            f"This usually indicates a very large file or system resource constraints. "
            f"Consider increasing timeout or checking system resources."
        )
        # Try to read logs for diagnostic info
        try:
            with open(stderr_log_path, "r") as f:
                stderr_content = f.read()
            if stderr_content.strip():
                logger.error(f"subprocess stderr before timeout:\n{stderr_content[-1000:]}")  # Last 1000 chars
        except Exception:
            pass
        return None
    except Exception as e:
        logger.error(f"Error running CNV analysis subprocess: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return None


class Result:
    """
    A class to store CNV results.
    """

    def __init__(self, cnv_dict: dict) -> None:
        self.cnv = cnv_dict


@dataclass
class CNVMetadata:
    """Container for CNV analysis metadata and results"""

    sample_id: str
    bam_path: str
    analysis_timestamp: float
    cnv_data_path: Optional[str] = None
    cnv_plot_path: Optional[str] = None
    sex_estimate: Optional[str] = None
    breakpoints: Optional[List[Dict]] = None
    chromosome_stats: Optional[Dict] = None
    processing_steps: List[str] = None
    error_message: Optional[str] = None

    def __post_init__(self):
        if self.processing_steps is None:
            self.processing_steps = []


def run_ruptures(
    r_cnv: list, penalty_value: int, bin_width: int
) -> List[Tuple[int, int]]:
    """
    Detect change points in a time series using kernel CPD with a linear (L2) kernel (mean-shift detection).

    Args:
        r_cnv: A list containing the time series data.
        penalty_value (int): The penalty value for the change point detection algorithm.
        bin_width (int): The width of the bins used to scale the detected change points.

    Returns:
        List[Tuple[float, float]]: A list of tuples where each tuple represents a detected change point range
                                    as (start, end) with respect to the bin width.
    """
    # Convert to numpy array once for better performance
    r_cnv_array = np.array(r_cnv)

    # Linear kernel (L2 cost): detects mean shifts; faster than RBF for long signals
    algo_c = rpt.KernelCPD(kernel="linear").fit(r_cnv_array)

    # Predict the change points using the provided penalty value
    ruptures_result = algo_c.predict(pen=penalty_value)

    # Compute the ranges around each change point using vectorized operations
    return [
        (cp * bin_width - bin_width, cp * bin_width + bin_width)
        for cp in ruptures_result
    ]


class CNV_Difference:
    """
    A class to store CNV difference data.
    """

    def __init__(self, *args, **kwargs) -> None:
        self.cnv = {}


def moving_average(data: np.ndarray, n: int = 3) -> np.ndarray:
    """
    Calculate the moving average of a given data array using scipy's uniform_filter1d.

    Args:
        data (np.ndarray): The data array.
        n (int): The window size for the moving average.

    Returns:
        np.ndarray: The array of moving averages.
    """
    return uniform_filter1d(data, size=n, mode="nearest")


def pad_arrays(
    arr1: np.ndarray, arr2: np.ndarray, pad_value: int = 0
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Pad two arrays to the same length with a specified value.

    Args:
        arr1 (np.ndarray): The first array.
        arr2 (np.ndarray): The second array.
        pad_value (int): The value to pad with.

    Returns:
        Tuple[np.ndarray, np.ndarray]: The padded arrays.
    """
    len1, len2 = len(arr1), len(arr2)
    max_len = max(len1, len2)

    # Use numpy's pad function more efficiently
    if len1 < max_len:
        arr1 = np.pad(
            arr1, (0, max_len - len1), mode="constant", constant_values=pad_value
        )
    if len2 < max_len:
        arr2 = np.pad(
            arr2, (0, max_len - len2), mode="constant", constant_values=pad_value
        )

    return arr1, arr2


def has_reads(bam_file: str) -> bool:
    """Quickly checks if a BAM file has any reads."""
    try:
        with pysam.AlignmentFile(bam_file, "rb") as bam:
            for _ in bam.fetch(until_eof=True):
                return True
        return False
    except Exception as e:
        # Replace print with logging
        logging.error(f"Error checking reads in {bam_file}: {e}")
        return False


def estimate_sex_from_cnv(cnv_data: Dict, logger) -> str:
    """
    Estimate genetic sex (XX or XY) based on CNV data.

    Args:
        cnv_data: Dictionary containing normalized CNV data
        logger: Logger instance

    Returns:
        String indicating estimated sex
    """
    try:
        if "chrX" not in cnv_data or "chrY" not in cnv_data:
            return "Unknown"

        # Calculate average CNV values for X and Y chromosomes more efficiently
        x_data = [i for i in cnv_data["chrX"] if i != 0]
        y_data = [i for i in cnv_data["chrY"] if i != 0]

        if not x_data or not y_data:
            return "Unknown"

        x_avg = np.mean(x_data)
        y_avg = np.mean(y_data)

        logger.debug(f"X chromosome average: {x_avg:.3f}")
        logger.debug(f"Y chromosome average: {y_avg:.3f}")

        # Sex estimation logic from original code
        if x_avg >= 0.1 and y_avg <= -0.1:
            return "Female"
        elif x_avg <= 0.1 and y_avg >= -0.2:
            return "Male"
        elif x_avg >= 0.1 and y_avg >= -0.1:
            return "Male (query X/Y copy number changes)"
        elif x_avg > 0.1 and y_avg > 0.1:
            return "Male (query X/Y copy number changes)"
        elif x_avg < 0.1 and y_avg < -0.2:
            return "Unknown (Query XY copy number changes)"
        else:
            return "Unknown"

    except Exception as e:
        logger.error(f"Error estimating sex: {e}")
        return "Unknown"


def detect_breakpoints_from_cnv(cnv_data: Dict, bin_width: int, logger) -> List[Dict]:
    """
    Detect breakpoints in CNV data using change point detection.

    Args:
        cnv_data: Dictionary containing CNV data
        bin_width: Width of bins
        logger: Logger instance

    Returns:
        List of detected breakpoints
    """
    breakpoints = []

    try:
        logger.debug("Detecting breakpoints in CNV data")

        # Simple sequential processing for all datasets
        for key, data in cnv_data.items():
            if key != "chrM" and len(data) > 3:
                try:
                    paired_changepoints = run_ruptures(
                        data, penalty_value=10, bin_width=bin_width
                    )

                    for start, end in paired_changepoints:
                        if start >= 0 and end > start:
                            breakpoints.append(
                                {
                                    "chromosome": key,
                                    "start": start,
                                    "end": end,
                                    "length": end - start,
                                }
                            )
                except Exception as e:
                    logger.debug(f"Error processing chromosome {key}: {e}")

        logger.debug(f"Detected {len(breakpoints)} breakpoints")
        return breakpoints

    except Exception as e:
        logger.error(f"Error in breakpoint detection: {e}")
        return []


def calculate_chromosome_stats_from_cnv(cnv_data: Dict, logger) -> Dict:
    """
    Calculate chromosome-wide statistics from normalized CNV data.

    Args:
        cnv_data: Dictionary containing normalized CNV data
        logger: Logger instance

    Returns:
        Dictionary of chromosome statistics
    """
    stats = {}

    try:
        logger.debug("Calculating chromosome statistics")

        for chrom, data in cnv_data.items():
            if len(data) == 0:
                continue

            # Calculate basic statistics more efficiently using numpy
            data_array = np.array(data)
            stats[chrom] = {
                "mean": float(np.mean(data_array)),
                "std": float(np.std(data_array)),
                "min": float(np.min(data_array)),
                "max": float(np.max(data_array)),
                "length": len(data_array),
            }

        logger.debug(f"Calculated statistics for {len(stats)} chromosomes")
        return stats

    except Exception as e:
        logger.error(f"Error calculating chromosome stats: {e}")
        return {}


def load_analysis_counter(sample_id: str, work_dir: str, logger) -> int:
    """Load the analysis counter for a sample from disk with caching"""
    # Check cache first
    cache = get_sample_cache(sample_id)
    if cache['analysis_counter'] is not None:
        logger.debug(f"Using cached analysis counter for {sample_id}: {cache['analysis_counter']}")
        return cache['analysis_counter']
    
    # Load from disk
    counter_file = os.path.join(work_dir, sample_id, "cnv_analysis_counter.txt")
    result = 0
    if os.path.exists(counter_file):
        try:
            with open(counter_file, "r") as f:
                result = int(f.read().strip())
        except (ValueError, IOError) as e:
            logger.warning(f"Error loading counter for {sample_id}: {e}")
    
    # Cache the result
    cache['analysis_counter'] = result
    logger.debug(f"Loaded and cached analysis counter for {sample_id}: {result}")
    return result


def save_analysis_counter(sample_id: str, counter: int, work_dir: str, logger) -> None:
    """Save the analysis counter for a sample to disk and update cache"""
    counter_file = os.path.join(work_dir, sample_id, "cnv_analysis_counter.txt")
    try:
        # Create sample directory if it doesn't exist
        os.makedirs(os.path.dirname(counter_file), exist_ok=True)
        with open(counter_file, "w") as f:
            f.write(str(counter))
        
        # Update cache
        update_sample_cache(sample_id, analysis_counter=counter)
        logger.debug(f"Saved and cached analysis counter for {sample_id}: {counter}")
    except IOError as e:
        logger.error(f"Error saving counter for {sample_id}: {e}")


def find_significant_regions(values, chrom, bin_width):
    """
    Find regions with significant CNV changes and merge adjacent regions.

    Args:
        values: CNV values for a chromosome
        chrom: Chromosome name
        bin_width: Width of bins in base pairs

    Returns:
        List of significant regions (merged adjacent regions)
    """
    regions = []
    threshold = 0.5  # CNV change threshold

    # Use numpy for better performance
    values_array = np.array(values)
    significant_indices = np.where(np.abs(values_array) > threshold)[0]

    # If no significant regions, return empty list
    if len(significant_indices) == 0:
        return regions

    # Merge adjacent indices into continuous regions
    current_start_idx = significant_indices[0]
    current_type = "gain" if values_array[significant_indices[0]] > 0 else "loss"

    for i in range(1, len(significant_indices)):
        idx = significant_indices[i]
        prev_idx = significant_indices[i - 1]
        idx_type = "gain" if values_array[idx] > 0 else "loss"

        # If adjacent and same type, continue the region
        if idx == prev_idx + 1 and idx_type == current_type:
            continue
        else:
            # Save the previous region
            regions.append(
                {
                    "chromosome": chrom,
                    "start": current_start_idx * bin_width,
                    "end": (significant_indices[i - 1] + 1) * bin_width,
                    "type": current_type,
                }
            )
            # Start new region
            current_start_idx = idx
            current_type = idx_type

    # Don't forget the last region
    regions.append(
        {
            "chromosome": chrom,
            "start": current_start_idx * bin_width,
            "end": (significant_indices[-1] + 1) * bin_width,
            "type": current_type,
        }
    )

    return regions


def generate_bed_files(bed_dir, analysis_counter, breakpoints, cnv_data, bin_width, logger):
    """
    Generate BED files for CNV regions and breakpoints.

    Args:
        bed_dir: BED files directory
        analysis_counter: Current analysis counter
        breakpoints: Detected breakpoints
        cnv_data: CNV data
        bin_width: Width of bins in base pairs
        logger: Logger instance
    """
    try:
        # Generate BED file for CNV regions
        bed_file = os.path.join(bed_dir, f"new_file_{analysis_counter:03d}.bed")
        with open(bed_file, "w") as f:
            for chrom, values in cnv_data.items():
                if len(values) > 0:
                    # Find regions with significant CNV changes
                    significant_regions = find_significant_regions(values, chrom, bin_width)
                    for region in significant_regions:
                        f.write(
                            f"{chrom}\t{region['start']}\t{region['end']}\t{region['type']}\n"
                        )

        # Generate BED file for breakpoints
        if breakpoints:
            breakpoint_file = os.path.join(
                bed_dir, f"breakpoints_{analysis_counter:03d}.bed"
            )
            with open(breakpoint_file, "w") as f:
                for bp in breakpoints:
                    f.write(
                        f"{bp['chromosome']}\t{bp['start']}\t{bp['end']}\tbreakpoint\n"
                    )

        logger.debug(f"Generated BED files in {bed_dir}")

    except Exception as e:
        logger.error(f"Error generating BED files: {e}")


def persist_ichor_py_outputs(
    sample_dir,
    result3_cnv: Dict,
    bin_width: int,
    logger,
) -> bool:
    """
    Run the ichor_py grid fit on normalized CNV (CNV3) and write
    ``ichor_py_params.json`` and ``ichor_py_segments.tsv`` under ``sample_dir``.

    Use this from the GUI to backfill ichor outputs without re-running the full
    CNV pipeline, as long as CNV3-like per-chromosome arrays and ``bin_width``
    are available.
    """
    try:
        from robin.analysis.ichor_py import run_ichor_like_from_cnv_dict

        ichor_out = run_ichor_like_from_cnv_dict(
            result3_cnv,
            bin_width=int(bin_width),
            normal_grid=(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
            ploidy_grid=(1.8, 2.0, 2.2, 2.5, 3.0, 3.5),
            min_cn=0,
            max_cn=6,
            self_transition=0.995,
        )
        ffit = ichor_out.fit
        params_payload = {
            "method": "ichor_py_grid_v0",
            "normal_fraction": float(ffit.normal_fraction),
            "tumor_fraction": float(ffit.tumor_fraction),
            "tumor_ploidy": float(ffit.tumor_ploidy),
            "log_likelihood": float(ffit.log_likelihood),
            "altered_genome_fraction": float(
                ichor_out.segments.altered_genome_fraction
            ),
            "n_segments": int(ichor_out.segments.n_segments),
            "n_altered_segments": int(ichor_out.segments.n_altered_segments),
            "state_copy_numbers": [
                int(v) for v in ffit.state_space.copy_numbers.tolist()
            ],
            "state_labels": [str(v) for v in ffit.state_space.labels],
        }
        if getattr(ffit, "tumor_ploidy_hmm_only", None) is not None:
            params_payload["tumor_ploidy_hmm_only"] = float(ffit.tumor_ploidy_hmm_only)
            params_payload["normal_fraction_hmm_only"] = float(
                ffit.normal_fraction_hmm_only
            )
            params_payload["tumor_fraction_hmm_only"] = float(
                1.0 - float(ffit.normal_fraction_hmm_only)
            )
            params_payload["log_likelihood_hmm_only"] = float(
                ffit.log_likelihood_hmm_only
            )
        if getattr(ffit, "xy_sex_used", None) is not None:
            params_payload["xy_sex_used"] = ffit.xy_sex_used
            params_payload["xy_sigma"] = ffit.xy_sigma
            params_payload["xy_weight"] = ffit.xy_weight
            if ffit.combined_score is not None:
                params_payload["combined_score"] = float(ffit.combined_score)
            if ffit.xy_penalty_at_best is not None:
                params_payload["xy_penalty_at_best"] = float(ffit.xy_penalty_at_best)
        params_path = os.path.join(sample_dir, "ichor_py_params.json")
        with open(params_path, "w", encoding="utf-8") as f:
            json.dump(params_payload, f, indent=2)

        seg_df = ichor_out.segments.segments.copy()
        if not seg_df.empty:
            seg_df = seg_df.sort_values(
                ["chrom", "start"], ascending=[True, True]
            ).reset_index(drop=True)
        seg_path = os.path.join(sample_dir, "ichor_py_segments.tsv")
        seg_df.to_csv(seg_path, sep="\t", index=False)
        logger.info(
            "[ichor_py] Saved params to %s and segments to %s",
            params_path,
            seg_path,
        )
        return True
    except Exception as e:
        logger.warning("[ichor_py] fit/save failed: %s", e)
        return False


def save_cnv_files(
    sample_dir,
    analysis_counter,
    r_cnv,
    r2_cnv,
    result3_cnv,
    bin_width,
    variance,
    breakpoints,
    sex_estimate,
    copy_numbers,
    logger,
    reference: Optional[str] = None,
    generate_master_bed: bool = False,
):
    """
    Save CNV data files in the specified format from the documentation.

    Args:
        sample_dir: Sample output directory
        analysis_counter: Current analysis counter
        r_cnv: Sample CNV data
        r2_cnv: Reference CNV data
        result3_cnv: Normalized CNV data
        bin_width: Bin width
        variance: Variance
        breakpoints: Detected breakpoints
        sex_estimate: Sex estimation result
        copy_numbers: Updated copy numbers for this sample (saved separately by caller)
        logger: Logger instance
    """
    try:
        # Save CNV data files as specified in documentation
        np.save(os.path.join(sample_dir, "CNV.npy"), r_cnv)
        np.save(os.path.join(sample_dir, "CNV2.npy"), r2_cnv)
        np.save(os.path.join(sample_dir, "CNV3.npy"), result3_cnv)

        # Save processing metadata
        cnv_dict = {
            "bin_width": bin_width,
            "variance": variance,
            "analysis_counter": analysis_counter,
        }
        np.save(os.path.join(sample_dir, "CNV_dict.npy"), cnv_dict)

        # Save breakpoint data as structured array
        if breakpoints:
            breakpoint_array = np.array(
                [(bp["chromosome"], bp["start"], bp["end"]) for bp in breakpoints],
                dtype=[("name", "U10"), ("start", "i8"), ("end", "i8")],
            )
            np.save(os.path.join(sample_dir, "cnv_data_array.npy"), breakpoint_array)

        # Helper for atomic pickle writes to avoid readers seeing partial files
        def _atomic_pickle_dump(obj, target_path: str) -> None:
            tmp_path = target_path + ".tmp"
            with open(tmp_path, "wb") as tf:
                pickle.dump(obj, tf, protocol=pickle.HIGHEST_PROTOCOL)
                tf.flush()
                os.fsync(tf.fileno())
            os.replace(tmp_path, target_path)

        # Save sex estimation (atomic)
        _atomic_pickle_dump(sex_estimate, os.path.join(sample_dir, "XYestimate.pkl"))

        # Note: copy_numbers is now saved separately by the caller as a per-sample file
        # This eliminates the multi-sample dictionary bottleneck

        # Create state directory for CNV detector
        state_dir = os.path.join(sample_dir, "cnv_detector_state")
        os.makedirs(state_dir, exist_ok=True)

        # Save state metadata
        state_metadata = {
            "analysis_counter": analysis_counter,
            "timestamp": time.time(),
            "bin_width": bin_width,
            "variance": variance,
            "breakpoints_count": len(breakpoints),
            "sex_estimate": sex_estimate,
        }
        _atomic_pickle_dump(
            state_metadata, os.path.join(state_dir, "tracker_metadata.pkl")
        )

        # Create BED files directory
        bed_dir = os.path.join(sample_dir, "bed_files")
        os.makedirs(bed_dir, exist_ok=True)

        # Generate BED files for CNV regions and breakpoints
        generate_bed_files(bed_dir, analysis_counter, breakpoints, result3_cnv, bin_width, logger)

        # Run python-native ichor-like purity/ploidy fit and persist compact outputs.
        # This uses normalized CNV (CNV3) and autosomes by default.
        persist_ichor_py_outputs(sample_dir, result3_cnv, int(bin_width), logger)

        # Generate master BED file only if requested (should only be done once per batch at the end)
        # Use async (non-blocking) generation to avoid blocking the analysis pipeline
        if generate_master_bed:
            try:
                from robin.analysis.master_bed_generator import (
                    generate_master_bed_async,
                    _try_get_target_panel_from_fusion_metadata,
                )
                
                # Extract sample_id from sample_dir
                sample_id = os.path.basename(sample_dir)
                # work_dir is the parent of sample_dir
                work_dir = os.path.dirname(sample_dir)
                
                # Try to get target_panel from fusion metadata if available
                target_panel = _try_get_target_panel_from_fusion_metadata(sample_id, work_dir)
                
                # Generate asynchronously (non-blocking)
                generate_master_bed_async(
                    sample_id=sample_id,
                    work_dir=work_dir,
                    analysis_counter=analysis_counter,
                    target_panel=target_panel,
                    logger_instance=logger,
                    reference=reference,
                )
            except Exception as e:
                logger.warning(f"Could not start async master BED generation: {e}")

        logger.debug(f"Saved all CNV files to {sample_dir}")

    except Exception as e:
        logger.error(f"Error saving CNV files: {e}")




def process_single_bam(bam_path, metadata, work_dir, logger, threads=2, reference: Optional[str] = None):
    """
    Process a single BAM file for CNV analysis using the complete pipeline.

    Args:
        bam_path: Path to BAM file
        metadata: BamMetadata object
        work_dir: Working directory
        logger: Logger instance
        threads: Number of threads to use for CNV analysis (default: 2, configurable)

    Returns:
        Dictionary with CNV analysis results
    """
    sample_id = metadata.extracted_data.get("sample_id", "unknown")
    
    # Set current sample and check if this is a new sample
    is_new_sample = set_current_sample(sample_id, logger)
    
    if is_new_sample:
        logger.info(f"🧬 Starting CNV analysis for NEW sample: {sample_id}")
    else:
        logger.info(f"🧬 Continuing CNV analysis for sample: {sample_id}")

    # Log essential metadata only
    logger.debug(f"BAM file: {metadata.file_path} ({metadata.file_size:,} bytes)")
    logger.debug(f"Sample ID: {sample_id}")
    logger.debug(f"Threads: {threads}")

    # Create sample-specific output directory
    sample_output_dir = os.path.join(work_dir, sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)

    # Load analysis counter from cache or disk
    analysis_counter = load_analysis_counter(sample_id, work_dir, logger)
    logger.debug(f"Analysis counter: {analysis_counter}")

    analysis_result = {
        "sample_id": sample_id,
        "bam_path": bam_path,
        "analysis_timestamp": time.time(),
        "analysis_counter": analysis_counter,
        "cnv_data_path": None,
        "sex_estimate": None,
        "breakpoints": [],
        "chromosome_stats": {},
        "processing_steps": [],
        "error_message": None,
    }

    try:
        # Use the sample output directory directly for intermediate files
        temp_dir = sample_output_dir

        # Check if BAM has reads
        if not has_reads(bam_path):
            logger.warning(f"No reads found in BAM file for {sample_id}")
            analysis_result["error_message"] = "No reads found in BAM file"
            analysis_result["processing_steps"].append("no_reads_found")
            return analysis_result

        analysis_result["processing_steps"].append("reads_found")

        # Check file size and adjust timeout accordingly
        bam_size_mb = os.path.getsize(bam_path) / (1024 * 1024)
        logger.debug(f"BAM file size: {bam_size_mb:.1f} MB")
        
        # Adaptive timeout based on file size
        # Base: 1 hour for files up to 100 MB
        # Add 1 minute per additional 10 MB for large files
        base_timeout = 3600  # 1 hour
        if bam_size_mb > 100:
            extra_timeout = int((bam_size_mb - 100) / 10) * 60  # 1 min per 10 MB
            adaptive_timeout = base_timeout + extra_timeout
            logger.info(
                f"Large file ({bam_size_mb:.1f} MB) - using extended timeout: {adaptive_timeout}s ({adaptive_timeout/60:.1f} min)"
            )
        else:
            adaptive_timeout = base_timeout

        # Use per-sample copy_numbers file (OPTIMIZED APPROACH)
        # This eliminates the multi-sample dictionary bottleneck
        copy_numbers_path = os.path.join(sample_output_dir, f"{sample_id}_copy_numbers.pkl")
        
        # Check cache first for copy_numbers
        cache = get_sample_cache(sample_id)
        if cache['copy_numbers'] is not None and cache['copy_numbers_path'] == copy_numbers_path:
            logger.debug(f"Using cached copy_numbers for {sample_id}")
            copy_numbers = cache['copy_numbers']
        else:
            # Load from disk or start fresh
            copy_numbers = None
            
            # Backward compatibility: Migrate from old multi-sample dict to per-sample file
            legacy_dict_path = os.path.join(sample_output_dir, "update_cnv_dict.pkl")
            if not os.path.exists(copy_numbers_path) and os.path.exists(legacy_dict_path):
                logger.info("Migrating from legacy multi-sample dict to per-sample file...")
                try:
                    with open(legacy_dict_path, "rb") as f:
                        update_cnv_dict = pickle.load(f)
                    if sample_id in update_cnv_dict:
                        # Extract this sample's data and save to per-sample file
                        with open(copy_numbers_path, "wb") as f:
                            pickle.dump(update_cnv_dict[sample_id], f, protocol=pickle.HIGHEST_PROTOCOL)
                        logger.info(f"Migrated copy_numbers for {sample_id} to per-sample file")
                        
                        # Optionally remove legacy file after migration (commented out for safety)
                        # os.remove(legacy_dict_path)
                except Exception as e:
                    logger.warning(f"Could not migrate legacy copy_numbers: {e}")
            
            if os.path.exists(copy_numbers_path):
                logger.debug(f"Loading copy_numbers from disk: {copy_numbers_path}")
                try:
                    with open(copy_numbers_path, "rb") as f:
                        copy_numbers = pickle.load(f)
                    logger.debug(f"Loaded copy_numbers for {sample_id}")
                except Exception as e:
                    logger.warning(f"Error loading copy_numbers for {sample_id}: {e}")
                    copy_numbers = {}
            else:
                logger.debug("No previous copy_numbers found, starting fresh")
                copy_numbers = {}
            
            # Cache the loaded copy_numbers
            update_sample_cache(sample_id, copy_numbers=copy_numbers, copy_numbers_path=copy_numbers_path)

        # Load reference CNV dict once (cached for subsequent samples)
        ref_cnv_path = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)),
            "HG01280_control_new.pkl",
        )
        ref_cnv_dict_loaded = get_cached_ref_cnv_dict(ref_cnv_path, logger)

        # Process BAM file with cnv_from_bam using configurable execution mode
        execution_mode = "subprocess" if USE_CNV_SUBPROCESS else "direct"
        logger.debug(f"Processing BAM file with cnv_from_bam ({execution_mode}, {threads} threads)")
        analysis_start = time.time()
        try:
            # Run CNV analysis using configurable execution mode
            if USE_CNV_SUBPROCESS:
                # Use subprocess execution (original approach)
                subprocess_result = run_cnv_analysis_subprocess(
                    bam_path,
                    {},  # Empty dict (not used when copy_numbers_path is provided)
                    ref_cnv_path,
                    temp_dir,
                    logger,
                    threads=threads,  # Use configurable threads parameter
                    mapq_filter=60,
                    sample_id=sample_id,
                    copy_numbers_path=copy_numbers_path,  # PER-SAMPLE FILE (OPTIMIZED)
                    timeout=adaptive_timeout,  # Adaptive timeout based on file size
                )
            else:
                # Use direct execution (OPTIMIZED: pass data directly, no file I/O)
                subprocess_result = run_cnv_analysis_direct(
                    bam_path,
                    copy_numbers,  # Pass copy_numbers dict directly
                    ref_cnv_dict_loaded,  # Pass reference dict directly
                    logger,
                    threads=threads,  # Use configurable threads parameter
                    mapq_filter=60,
                    sample_id=sample_id,
                )

            if subprocess_result is None or not subprocess_result.get("success", False):
                error_msg = (
                    subprocess_result.get("error", "Unknown error")
                    if subprocess_result
                    else "Subprocess failed"
                )
                raise RuntimeError(f"CNV analysis subprocess failed: {error_msg}")

            # Extract results from subprocess
            r_cnv = subprocess_result["r_cnv"]
            r_bin = subprocess_result["r_bin"]
            r_var = subprocess_result["r_var"]
            genome_length = subprocess_result["genome_length"]
            r2_cnv = subprocess_result["r2_cnv"]
            updated_copy_numbers = subprocess_result["updated_copy_numbers"]
            
            # Log timing information from analysis
            if "timing" in subprocess_result:
                timing = subprocess_result["timing"]
                logger.info(f"CNV analysis timing: Pass1={timing['pass1_time']:.2f}s, Pass2={timing['pass2_time']:.2f}s, Total={timing['total_time']:.2f}s")
            
            analysis_elapsed = time.time() - analysis_start
            logger.info(f"CNV analysis ({execution_mode}) completed in {analysis_elapsed:.2f}s (total with overhead)")

            # Save updated copy_numbers back to per-sample file (OPTIMIZED)
            save_start = time.time()
            with open(copy_numbers_path, "wb") as f:
                pickle.dump(updated_copy_numbers, f, protocol=pickle.HIGHEST_PROTOCOL)
            save_time = time.time() - save_start
            logger.debug(f"Saved updated copy_numbers to {copy_numbers_path} in {save_time:.3f}s")
            
            # Update cache with the saved copy_numbers
            update_sample_cache(sample_id, copy_numbers=updated_copy_numbers)
            
            # Update the local copy_numbers variable for subsequent processing
            copy_numbers = updated_copy_numbers

            analysis_result["processing_steps"].append("cnv_data_extracted")
            logger.debug("CNV data extracted successfully")
            logger.debug(
                f"Bin width: {r_bin:,}, Variance: {r_var:.6f}, Genome length: {genome_length:,}"
            )

        except Exception as e:
            logger.error(f"Error extracting CNV data: {e}")
            analysis_result["error_message"] = f"CNV extraction failed: {str(e)}"
            analysis_result["processing_steps"].append("cnv_extraction_failed")
            return analysis_result


        # Calculate normalized CNV data (difference between sample and reference)
        logger.debug("Calculating normalized CNV data")
        
        result3_cnv = {}
        for key in r_cnv.keys():
            if key != "chrM" and key in r2_cnv:
                moving_avg_data1 = moving_average(r_cnv[key])
                moving_avg_data2 = moving_average(r2_cnv[key])
                moving_avg_data1, moving_avg_data2 = pad_arrays(
                    moving_avg_data1, moving_avg_data2
                )
                result3_cnv[key] = moving_avg_data1 - moving_avg_data2

        analysis_result["processing_steps"].append("normalized_cnv_calculated")

        # Estimate sex from CNV data
        sex_estimate = estimate_sex_from_cnv(result3_cnv, logger)
        analysis_result["sex_estimate"] = sex_estimate
        analysis_result["processing_steps"].append("sex_estimated")

        # Detect breakpoints
        breakpoints = detect_breakpoints_from_cnv(r_cnv, r_bin, logger)
        analysis_result["breakpoints"] = breakpoints
        analysis_result["processing_steps"].append("breakpoints_detected")

        # Calculate chromosome statistics
        chromosome_stats = calculate_chromosome_stats_from_cnv(result3_cnv, logger)
        analysis_result["chromosome_stats"] = chromosome_stats
        analysis_result["processing_steps"].append("chromosome_stats_calculated")

        # Save CNV data in the specified format
        analysis_counter += 1

        # Save CNV data files as specified in the documentation
        # Don't generate master BED here - it should only be generated once per batch
        save_cnv_files(
            sample_output_dir,
            analysis_counter,
            r_cnv,
            r2_cnv,
            result3_cnv,
            r_bin,
            r_var,
            breakpoints,
            sex_estimate,
            updated_copy_numbers,  # Pass updated copy_numbers (already saved to per-sample file above)
            logger,
            reference=reference,
            generate_master_bed=False,  # Single BAM - master BED should be generated at batch end
        )

        analysis_result["cnv_data_path"] = os.path.join(
            sample_output_dir, f"{analysis_counter}_cnv_data.json"
        )
        analysis_result["processing_steps"].append("cnv_data_saved")

        # Save updated analysis counter to disk
        save_analysis_counter(sample_id, analysis_counter, work_dir, logger)
        analysis_result["analysis_counter"] = analysis_counter
        logger.debug(f"Saved analysis counter: {analysis_counter}")

        # Force garbage collection
        gc.collect()

        analysis_result["processing_steps"].append("cnv_analysis_complete")

        # Replace print with logging
        logger.debug(f"Analysis result: {analysis_result}")
        logger.info(f"CNV analysis complete for {sample_id}")
        logger.info(f"Sex Estimate: {analysis_result['sex_estimate']}")
        logger.info(f"Breakpoints: {len(analysis_result['breakpoints'])} detected")
        
        return analysis_result

    except Exception as e:
        logger.error(f"Error in CNV analysis for {sample_id}: {e}")
        analysis_result["error_message"] = str(e)
        analysis_result["processing_steps"].append("analysis_failed")
        return analysis_result


def process_multiple_bams(bam_paths, metadata_list, work_dir, logger, threads=2, reference: Optional[str] = None):
    """
    Process multiple BAM files for CNV analysis using aggregated CNV data.
    
    This function processes multiple BAM files for the same sample, accumulating
    CNV data across all files before performing downstream analysis. This is
    more efficient than processing each BAM file individually and then trying
    to combine results.

    Args:
        bam_paths: List of paths to BAM files
        metadata_list: List of BamMetadata objects (one per BAM file)
        work_dir: Working directory
        logger: Logger instance
        threads: Number of threads to use for CNV analysis (default: 2)

    Returns:
        Dictionary with aggregated CNV analysis results
    """
    if not bam_paths or not metadata_list:
        raise ValueError("bam_paths and metadata_list must not be empty")
    
    if len(bam_paths) != len(metadata_list):
        raise ValueError("bam_paths and metadata_list must have the same length")
    
    # Get sample ID from first metadata (assuming all BAMs are from same sample)
    sample_id = metadata_list[0].extracted_data.get("sample_id", "unknown")
    
    # Set current sample and check if this is a new sample
    is_new_sample = set_current_sample(sample_id, logger)
    
    if is_new_sample:
        logger.info(f"🧬 Starting multi-BAM CNV analysis for NEW sample: {sample_id}")
    else:
        logger.info(f"🧬 Continuing multi-BAM CNV analysis for sample: {sample_id}")

    logger.info(f"Processing {len(bam_paths)} BAM files for sample {sample_id}")
    
    # Log essential metadata only
    for i, (bam_path, metadata) in enumerate(zip(bam_paths, metadata_list)):
        logger.debug(f"BAM file {i+1}: {metadata.file_path} ({metadata.file_size:,} bytes)")
    
    logger.debug(f"Sample ID: {sample_id}")
    logger.debug(f"Threads: {threads}")

    # Create sample-specific output directory
    sample_output_dir = os.path.join(work_dir, sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)

    # Load analysis counter from cache or disk
    analysis_counter = load_analysis_counter(sample_id, work_dir, logger)
    logger.debug(f"Analysis counter: {analysis_counter}")

    analysis_result = {
        "sample_id": sample_id,
        "bam_paths": bam_paths,
        "analysis_timestamp": time.time(),
        "analysis_counter": analysis_counter,
        "cnv_data_path": None,
        "sex_estimate": None,
        "breakpoints": [],
        "chromosome_stats": {},
        "processing_steps": [],
        "error_message": None,
        "files_processed": 0,
        "total_files": len(bam_paths),
    }

    try:
        # Use the sample output directory directly for intermediate files
        temp_dir = sample_output_dir

        # Check if any BAM has reads
        valid_bam_paths = []
        valid_metadata_list = []
        
        for bam_path, metadata in zip(bam_paths, metadata_list):
            if has_reads(bam_path):
                valid_bam_paths.append(bam_path)
                valid_metadata_list.append(metadata)
            else:
                logger.warning(f"No reads found in BAM file: {bam_path}")
        
        if not valid_bam_paths:
            logger.warning(f"No valid BAM files found for {sample_id}")
            analysis_result["error_message"] = "No reads found in any BAM files"
            analysis_result["processing_steps"].append("no_reads_found")
            return analysis_result

        logger.info(f"Processing {len(valid_bam_paths)} valid BAM files out of {len(bam_paths)} total")
        analysis_result["files_processed"] = len(valid_bam_paths)
        analysis_result["processing_steps"].append("reads_found")

        # Use per-sample copy_numbers file (OPTIMIZED APPROACH)
        copy_numbers_path = os.path.join(sample_output_dir, f"{sample_id}_copy_numbers.pkl")
        
        # Check cache first for copy_numbers
        cache = get_sample_cache(sample_id)
        if cache['copy_numbers'] is not None and cache['copy_numbers_path'] == copy_numbers_path:
            logger.debug(f"Using cached copy_numbers for {sample_id}")
            copy_numbers = cache['copy_numbers']
        else:
            # Load from disk or start fresh
            copy_numbers = None
            
            # Backward compatibility: Migrate from old multi-sample dict to per-sample file
            legacy_dict_path = os.path.join(sample_output_dir, "update_cnv_dict.pkl")
            if not os.path.exists(copy_numbers_path) and os.path.exists(legacy_dict_path):
                logger.info("Migrating from legacy multi-sample dict to per-sample file...")
                try:
                    with open(legacy_dict_path, "rb") as f:
                        update_cnv_dict = pickle.load(f)
                    if sample_id in update_cnv_dict:
                        # Extract this sample's data and save to per-sample file
                        with open(copy_numbers_path, "wb") as f:
                            pickle.dump(update_cnv_dict[sample_id], f, protocol=pickle.HIGHEST_PROTOCOL)
                        logger.info(f"Migrated copy_numbers for {sample_id} to per-sample file")
                except Exception as e:
                    logger.warning(f"Could not migrate legacy copy_numbers: {e}")
            
            if os.path.exists(copy_numbers_path):
                logger.debug(f"Loading copy_numbers from disk: {copy_numbers_path}")
                try:
                    with open(copy_numbers_path, "rb") as f:
                        copy_numbers = pickle.load(f)
                    logger.debug(f"Loaded copy_numbers for {sample_id}")
                except Exception as e:
                    logger.warning(f"Error loading copy_numbers for {sample_id}: {e}")
                    copy_numbers = {}
            else:
                logger.debug("No previous copy_numbers found, starting fresh")
                copy_numbers = {}
            
            # Cache the loaded copy_numbers
            update_sample_cache(sample_id, copy_numbers=copy_numbers, copy_numbers_path=copy_numbers_path)

        # Load reference CNV dict once (cached for subsequent samples)
        ref_cnv_path = os.path.join(
            os.path.dirname(os.path.abspath(resources.__file__)),
            "HG01280_control_new.pkl",
        )
        ref_cnv_dict_loaded = get_cached_ref_cnv_dict(ref_cnv_path, logger)

        # Process all BAM files with cnv_from_bam using configurable execution mode
        execution_mode = "subprocess" if USE_CNV_SUBPROCESS else "direct"
        logger.debug(f"Processing {len(valid_bam_paths)} BAM files with cnv_from_bam ({execution_mode}, {threads} threads)")
        
        # Calculate adaptive timeout based on total file sizes
        total_bam_size_mb = sum(os.path.getsize(bam_path) for bam_path in valid_bam_paths) / (1024 * 1024)
        logger.debug(f"Total BAM file size: {total_bam_size_mb:.1f} MB")
        
        # Adaptive timeout based on total file size
        base_timeout = 3600  # 1 hour
        if total_bam_size_mb > 100:
            extra_timeout = int((total_bam_size_mb - 100) / 10) * 60  # 1 min per 10 MB
            adaptive_timeout = base_timeout + extra_timeout
            logger.info(
                f"Large total file size ({total_bam_size_mb:.1f} MB) - using extended timeout: {adaptive_timeout}s ({adaptive_timeout/60:.1f} min)"
            )
        else:
            adaptive_timeout = base_timeout

        analysis_start = time.time()
        try:
            # Process each BAM file, accumulating CNV data in copy_numbers
            for i, (bam_path, metadata) in enumerate(zip(valid_bam_paths, valid_metadata_list)):
                logger.info(f"Processing BAM file {i+1}/{len(valid_bam_paths)}: {os.path.basename(bam_path)}")
                
                # Run CNV analysis using configurable execution mode
                if USE_CNV_SUBPROCESS:
                    # Use subprocess execution (original approach)
                    subprocess_result = run_cnv_analysis_subprocess(
                        bam_path,
                        {},  # Empty dict (not used when copy_numbers_path is provided)
                        ref_cnv_path,
                        temp_dir,
                        logger,
                        threads=threads,
                        mapq_filter=60,
                        sample_id=sample_id,
                        copy_numbers_path=copy_numbers_path,  # PER-SAMPLE FILE (OPTIMIZED)
                        timeout=adaptive_timeout,  # Adaptive timeout based on file size
                    )
                else:
                    # Use direct execution (OPTIMIZED: pass data directly, no file I/O)
                    subprocess_result = run_cnv_analysis_direct(
                        bam_path,
                        copy_numbers,  # Pass copy_numbers dict directly - this accumulates data
                        ref_cnv_dict_loaded,  # Pass reference dict directly
                        logger,
                        threads=threads,
                        mapq_filter=60,
                        sample_id=sample_id,
                    )

                if subprocess_result is None or not subprocess_result.get("success", False):
                    error_msg = (
                        subprocess_result.get("error", "Unknown error")
                        if subprocess_result
                        else "Subprocess failed"
                    )
                    raise RuntimeError(f"CNV analysis subprocess failed for {os.path.basename(bam_path)}: {error_msg}")

                # Extract results from subprocess (only need the final accumulated copy_numbers)
                updated_copy_numbers = subprocess_result["updated_copy_numbers"]
                
                # Log timing information from analysis
                if "timing" in subprocess_result:
                    timing = subprocess_result["timing"]
                    logger.debug(f"BAM {i+1} CNV analysis timing: Pass1={timing['pass1_time']:.2f}s, Pass2={timing['pass2_time']:.2f}s, Total={timing['total_time']:.2f}s")
                
                # Update copy_numbers for next iteration
                copy_numbers = updated_copy_numbers
                
                logger.info(f"Completed BAM file {i+1}/{len(valid_bam_paths)}: {os.path.basename(bam_path)}")

            # After processing all BAM files, get the final aggregated results
            # We need to run one final analysis to get the aggregated CNV data
            logger.info("Generating final aggregated CNV results...")
            
            # Use the last BAM file to generate final results (copy_numbers now contains aggregated data)
            final_bam_path = valid_bam_paths[-1]
            
            if USE_CNV_SUBPROCESS:
                # For subprocess mode, we need to run one more analysis to get final results
                final_result = run_cnv_analysis_subprocess(
                    final_bam_path,
                    {},
                    ref_cnv_path,
                    temp_dir,
                    logger,
                    threads=threads,
                    mapq_filter=60,
                    sample_id=sample_id,
                    copy_numbers_path=copy_numbers_path,
                    timeout=adaptive_timeout,
                )
            else:
                # For direct mode, run one final analysis to get aggregated results
                final_result = run_cnv_analysis_direct(
                    final_bam_path,
                    copy_numbers,  # This now contains aggregated data from all BAM files
                    ref_cnv_dict_loaded,
                    logger,
                    threads=threads,
                    mapq_filter=60,
                    sample_id=sample_id,
                )

            if final_result is None or not final_result.get("success", False):
                error_msg = (
                    final_result.get("error", "Unknown error")
                    if final_result
                    else "Final analysis failed"
                )
                raise RuntimeError(f"Final CNV analysis failed: {error_msg}")

            # Extract final aggregated results
            r_cnv = final_result["r_cnv"]
            r_bin = final_result["r_bin"]
            r_var = final_result["r_var"]
            genome_length = final_result["genome_length"]
            r2_cnv = final_result["r2_cnv"]
            final_copy_numbers = final_result["updated_copy_numbers"]
            
            analysis_elapsed = time.time() - analysis_start
            logger.info(f"Multi-BAM CNV analysis ({execution_mode}) completed in {analysis_elapsed:.2f}s (total with overhead)")

            # Save updated copy_numbers back to per-sample file (OPTIMIZED)
            save_start = time.time()
            with open(copy_numbers_path, "wb") as f:
                pickle.dump(final_copy_numbers, f, protocol=pickle.HIGHEST_PROTOCOL)
            update_sample_cache(sample_id, copy_numbers=final_copy_numbers)
            copy_numbers = final_copy_numbers
            logger.info(
                f"[cnv] Save copy_numbers + cache update completed in {time.time() - save_start:.2f}s"
            )

            analysis_result["processing_steps"].append("cnv_data_extracted")
            logger.debug("Aggregated CNV data extracted successfully")
            logger.debug(
                f"Bin width: {r_bin:,}, Variance: {r_var:.6f}, Genome length: {genome_length:,}"
            )

        except Exception as e:
            logger.error(f"Error extracting aggregated CNV data: {e}")
            analysis_result["error_message"] = f"CNV extraction failed: {str(e)}"
            analysis_result["processing_steps"].append("cnv_extraction_failed")
            return analysis_result

        # Calculate normalized CNV data (difference between sample and reference)
        logger.debug("Calculating normalized CNV data from aggregated results")
        t0 = time.time()
        result3_cnv = {}
        for key in r_cnv.keys():
            if key != "chrM" and key in r2_cnv:
                moving_avg_data1 = moving_average(r_cnv[key])
                moving_avg_data2 = moving_average(r2_cnv[key])
                moving_avg_data1, moving_avg_data2 = pad_arrays(
                    moving_avg_data1, moving_avg_data2
                )
                result3_cnv[key] = moving_avg_data1 - moving_avg_data2
        logger.info(f"[cnv] Normalized CNV (moving avg + diff) completed in {time.time() - t0:.2f}s")

        analysis_result["processing_steps"].append("normalized_cnv_calculated")

        # Estimate sex from aggregated CNV data
        t0 = time.time()
        sex_estimate = estimate_sex_from_cnv(result3_cnv, logger)
        logger.info(f"[cnv] Sex estimate completed in {time.time() - t0:.2f}s")
        analysis_result["sex_estimate"] = sex_estimate
        analysis_result["processing_steps"].append("sex_estimated")

        # Detect breakpoints from aggregated CNV data
        t0 = time.time()
        breakpoints = detect_breakpoints_from_cnv(r_cnv, r_bin, logger)
        logger.info(f"[cnv] Breakpoint detection completed in {time.time() - t0:.2f}s")
        analysis_result["breakpoints"] = breakpoints
        analysis_result["processing_steps"].append("breakpoints_detected")

        # Calculate chromosome statistics from aggregated CNV data
        t0 = time.time()
        chromosome_stats = calculate_chromosome_stats_from_cnv(result3_cnv, logger)
        logger.info(f"[cnv] Chromosome stats completed in {time.time() - t0:.2f}s")
        analysis_result["chromosome_stats"] = chromosome_stats
        analysis_result["processing_steps"].append("chromosome_stats_calculated")

        # Save CNV data in the specified format
        analysis_counter += 1

        # Save CNV data files as specified in the documentation
        # Generate master BED here since this is batch processing (all BAMs processed)
        t0 = time.time()
        save_cnv_files(
            sample_output_dir,
            analysis_counter,
            r_cnv,
            r2_cnv,
            result3_cnv,
            r_bin,
            r_var,
            breakpoints,
            sex_estimate,
            final_copy_numbers,  # Pass final aggregated copy_numbers
            logger,
            reference=reference,
            generate_master_bed=True,  # Batch processing - generate master BED once at the end
        )
        logger.info(f"[cnv] save_cnv_files (incl. master BED) completed in {time.time() - t0:.2f}s")

        analysis_result["cnv_data_path"] = os.path.join(
            sample_output_dir, f"{analysis_counter}_cnv_data.json"
        )
        analysis_result["processing_steps"].append("cnv_data_saved")

        # Save updated analysis counter to disk
        t0 = time.time()
        save_analysis_counter(sample_id, analysis_counter, work_dir, logger)
        analysis_result["analysis_counter"] = analysis_counter
        gc.collect()
        logger.info(f"[cnv] Save analysis counter + gc completed in {time.time() - t0:.2f}s")

        analysis_result["processing_steps"].append("cnv_analysis_complete")

        # Replace print with logging
        logger.debug(f"Analysis result: {analysis_result}")
        logger.info(f"Multi-BAM CNV analysis complete for {sample_id}")
        logger.info(f"Files processed: {analysis_result['files_processed']}/{analysis_result['total_files']}")
        logger.info(f"Sex Estimate: {analysis_result['sex_estimate']}")
        logger.info(f"Breakpoints: {len(analysis_result['breakpoints'])} detected")
        
        return analysis_result

    except Exception as e:
        logger.error(f"Error in multi-BAM CNV analysis for {sample_id}: {e}")
        analysis_result["error_message"] = str(e)
        analysis_result["processing_steps"].append("analysis_failed")
        return analysis_result


def cnv_handler(job, work_dir=None, target_panel=None, threads=2):
    """
    Handler function for CNV analysis jobs.
    This function processes BAM files for copy number variation analysis.

    Args:
        job: The workflow job containing file and metadata
        work_dir: Optional base directory for output (defaults to BAM file directory)
        target_panel: Target panel type (required)
        threads: Number of threads to use for CNV analysis (default: 2, configurable via job metadata)
    """
    # Validate required parameters
    if not target_panel:
        raise ValueError("target_panel is required for CNV analysis")
    
    # Get job-specific logger
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
    
    # Check if this is a batched job
    batched_job = job.context.metadata.get("_batched_job")
    if batched_job:
        batch_size = batched_job.get_file_count()
        sample_id = batched_job.get_sample_id()
        batch_id = batched_job.batch_id
        logger.info(f"Processing CNV batch: {batch_size} files for sample '{sample_id}' (batch_id: {batch_id})")
        
        # Get all filepaths in the batch
        filepaths = batched_job.get_filepaths()
        
        # Log individual files in the batch
        for i, filepath in enumerate(filepaths):
            logger.info(f"  Batch file {i+1}/{batch_size}: {os.path.basename(filepath)}")
        
        # Prepare metadata list for all BAM files in the batch
        metadata_list = []
        for i, bam_path in enumerate(filepaths):
            # Get metadata from preprocessing for this specific file
            # Note: Each file in the batch should have its own metadata
            file_metadata = batched_job.contexts[i].metadata.get("bam_metadata", {})
            
            # Get sample ID from preprocessing results for this specific file
            file_context = batched_job.contexts[i]
            file_sample_id = file_context.get_sample_id()
            
            # Use the sample ID from the file's context (which should have preprocessing results)
            if file_sample_id != "unknown":
                file_metadata["sample_id"] = file_sample_id
            else:
                file_metadata["sample_id"] = sample_id
            
            # Create BamMetadata object for compatibility
            from robin.analysis.bam_preprocessor import BamMetadata
            
            metadata = BamMetadata(
                file_path=bam_path,
                file_size=batched_job.contexts[i].metadata.get("file_size", 0),
                creation_time=batched_job.contexts[i].metadata.get("creation_time", time.time()),
                extracted_data=file_metadata,
            )
            metadata_list.append(metadata)
        
        # Determine work directory for the batch
        if work_dir is None:
            # Default to first BAM file directory
            batch_work_dir = os.path.dirname(filepaths[0])
        else:
            # Use specified work directory, create if it doesn't exist
            os.makedirs(work_dir, exist_ok=True)
            batch_work_dir = work_dir
            logger.debug(f"Using specified work directory: {batch_work_dir}")
        
        # Get reference from job metadata
        reference = job.context.metadata.get("reference")
        if reference:
            # Expand user home directory if present
            reference = os.path.expanduser(reference)
            logger.debug(f"Using reference genome from job metadata: {reference}")
        
        # Process all BAM files in the batch using the new aggregated function
        logger.info(f"Processing {batch_size} BAM files as aggregated batch for sample '{sample_id}'")
        
        batch_result = process_multiple_bams(
            bam_paths=filepaths,
            metadata_list=metadata_list,
            work_dir=batch_work_dir,
            logger=logger,
            threads=threads,
            reference=reference,
        )
        
        # Store batch results in job context (maintain compatibility with existing structure)
        job.context.add_metadata("cnv_analysis", {
            "batch_result": batch_result,  # Single aggregated result
            "batch_size": batch_size,
            "sample_id": sample_id,
            "batch_id": batch_id,
            "files_processed": batch_result.get("files_processed", batch_size),
            "total_files": batch_result.get("total_files", batch_size)
        })
        
        logger.info(f"Completed CNV batch processing: {batch_size} files for sample '{sample_id}'")
        logger.info(f"Files successfully processed: {batch_result.get('files_processed', batch_size)}/{batch_result.get('total_files', batch_size)}")
        
        if batch_result.get("error_message"):
            logger.error(f"Batch processing completed with errors: {batch_result['error_message']}")
        else:
            logger.info("Batch processing completed successfully with aggregated CNV analysis")
        
        return
        
    else:
        # CNV jobs should always be batched - raise an error if not
        error_msg = f"CNV job received without batching metadata. Expected batched job but got single file: {os.path.basename(job.context.filepath)}"
        logger.error(error_msg)
        raise ValueError(error_msg)
