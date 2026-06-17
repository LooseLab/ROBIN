"""BAM file preprocessing and metadata extraction for robin.

Requires Python 3.12+ (slots=True, str.removeprefix, PEP 709 comprehensions).
"""
from __future__ import annotations

import os
import sys
import time
import re
import hashlib
from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List, Tuple
from pathlib import Path

if sys.version_info < (3, 12):
    raise RuntimeError("robin bam_preprocessor requires Python 3.12 or newer")

import pysam
from dateutil import parser
from robin.analysis.master_csv_manager import MasterCSVManager
from robin.logging_config import get_job_logger


# ============================================================================
# CONSTANTS AND CONFIGURATION
# ============================================================================

# Compile regex pattern once for barcode detection (module-level constant)
_BARCODE_PATTERN = re.compile(r"_barcode(\d{1,2})$")

# Pre-compile common string operations
_PASS_STR = "pass"
_FAIL_STR = "fail"
_RG_TAG = "RG"
_ST_TAG = "st"
_BASECALL_MODEL_PREFIX = "basecall_model="
_RUNID_PREFIX = "runid="
_MODBASE_MODELS_KEY = "modbase_models"
_CPG_MODBASE_MARKER = "5mCG_5hmCG"
_ALL_CONTEXT_MODBASE_MARKER = "5mC_5hmC"
_UNRESOLVED_MODBASE_MODEL_VALUES = frozenset({"modbase_model_version_id"})

# Optional: configure BAM read threads via environment variable (pysam/htslib BGZF threads)
# Set LJ_BAM_THREADS=4 (or higher) to enable multi-threaded decompression when reading BAMs.
# Defaults to single-threaded when unset or invalid to maintain behavior.
try:
    _LJ_BAM_THREADS = int(os.getenv("LJ_BAM_THREADS", "0"))
    if _LJ_BAM_THREADS < 0:
        _LJ_BAM_THREADS = 0
except ValueError:
    _LJ_BAM_THREADS = 0

# When True, BAMs with >50k reads are processed (downstream batching uses batch size 1 per file).
# Set ROBIN_PROCESS_LARGE_BAMS=1 to enable.
_PROCESS_LARGE_BAMS_INDIVIDUALLY = os.getenv("ROBIN_PROCESS_LARGE_BAMS", "0").strip().lower() in ("1", "true", "yes", "on")

# Constants for MGMT locus (chr10:129,466,536-129,467,536)
_MGMT_CHR = "chr10"
_MGMT_START = 129466536
_MGMT_END = 129467536

# ============================================================================
# DATA STRUCTURES
# ============================================================================


@dataclass(slots=True)
class BamMetadata:
    """Container for BAM file metadata and extracted data."""

    file_path: str
    file_size: int
    creation_time: float
    extracted_data: Dict[str, Any] = field(default_factory=dict)
    processing_steps: List[str] = field(default_factory=list)

    def __post_init__(self):
        if self.extracted_data is None:
            self.extracted_data = {}
        if self.processing_steps is None:
            self.processing_steps = []


def _persist_supplementary_read_ids(
    metadata: BamMetadata,
    bam_path: str,
    work_dir: str,
) -> None:
    """Persist a complete supplementary-read ID set under a BAM-specific path."""
    supp_ids = metadata.extracted_data.get("supplementary_read_ids", [])
    metadata.extracted_data["supplementary_read_ids_complete"] = True
    metadata.extracted_data["supplementary_read_ids_count"] = len(supp_ids)
    if not supp_ids:
        metadata.extracted_data.pop("supplementary_read_ids", None)
        metadata.extracted_data.pop("supplementary_read_ids_path", None)
        return

    sample_id = metadata.extracted_data.get("sample_id", "unknown")
    supp_dir = os.path.join(work_dir, sample_id, "_supplementary_read_ids")
    os.makedirs(supp_dir, exist_ok=True)
    path_hash = hashlib.sha256(os.path.abspath(bam_path).encode("utf-8")).hexdigest()[:16]
    supp_path = os.path.join(
        supp_dir,
        f"{os.path.basename(bam_path)}.{path_hash}.txt",
    )
    tmp_path = f"{supp_path}.tmp"
    with open(tmp_path, "w") as f:
        for read_id in sorted(supp_ids):
            f.write(f"{read_id}\n")
    os.replace(tmp_path, supp_path)
    metadata.extracted_data["supplementary_read_ids_path"] = supp_path
    metadata.extracted_data.pop("supplementary_read_ids", None)


# ============================================================================
# CORE EXTRACTION FUNCTIONS
# ============================================================================


def get_rg_tags_from_bam(sam_file) -> Optional[Tuple[Optional[str], ...]]:
    """
    Extracts RG (Read Group) tags from the BAM file header.
    This is the first step in metadata extraction.

    Args:
        sam_file: pysam AlignmentFile object

    Returns:
        Optional[Tuple[Optional[str], ...]]: A tuple containing RG tag information
                                             or None if no RG tags are found.
    """
    if not sam_file:
        return None

    rg_tags = sam_file.header.get(_RG_TAG, [])
    if not rg_tags:
        return None

    # Process only the first RG tag for efficiency
    rg_tag = rg_tags[0]
    id_tag = rg_tag.get("ID")
    dt_tag = rg_tag.get("DT")
    ds_tag = rg_tag.get("DS", "")

    ds_fields = {}
    for token in ds_tag.split():
        key, separator, value = token.partition("=")
        if separator:
            ds_fields[key] = value
    basecall_model_tag = ds_fields.get(_BASECALL_MODEL_PREFIX.removesuffix("="))
    runid_tag = ds_fields.get(_RUNID_PREFIX.removesuffix("="))
    modbase_models_tag = ds_fields.get(_MODBASE_MODELS_KEY)

    lb_tag = rg_tag.get("LB")
    pl_tag = rg_tag.get("PL")
    pm_tag = rg_tag.get("PM")
    pu_tag = rg_tag.get("PU")
    al_tag = rg_tag.get("al")

    return (
        id_tag,
        dt_tag,
        basecall_model_tag,
        runid_tag,
        lb_tag,
        pl_tag,
        pm_tag,
        pu_tag,
        al_tag,
        modbase_models_tag,
    )


def _is_unresolved_modbase_model(modbase_models: Optional[str]) -> bool:
    """Return True when the BAM only records an unresolved modbase model placeholder."""
    if not modbase_models:
        return False
    return modbase_models.strip().lower() in _UNRESOLVED_MODBASE_MODEL_VALUES


def _get_modbase_model_warning_level(modbase_models: Optional[str]) -> str:
    """Return ``info`` for unknown model placeholders, else ``warning``."""
    if modbase_models and _is_unresolved_modbase_model(modbase_models):
        return "info"
    return "warning"


def _get_modbase_model_warning(modbase_models: Optional[str]) -> Optional[str]:
    """Return a user-facing warning for unsupported methylation model settings."""
    if modbase_models and _CPG_MODBASE_MARKER in modbase_models:
        return None
    if modbase_models and _ALL_CONTEXT_MODBASE_MARKER in modbase_models:
        return (
            f"BAM uses all-context methylation calling ({modbase_models}). "
            "Methylation classifications may be incorrect and slower than expected. "
            "Use 5mCG_5hmCG modbase calling, which restricts methylation calling "
            "to CpG contexts."
        )
    if modbase_models and _is_unresolved_modbase_model(modbase_models):
        return (
            "The BAM header does not record which modbase model was used. "
            "Compatibility with the recommended 5mCG_5hmCG configuration could not "
            "be verified."
        )
    if modbase_models:
        return (
            f"BAM reports modbase models without 5mCG_5hmCG ({modbase_models}). "
            "Methylation classifications may be incorrect. Use 5mCG_5hmCG "
            "modbase calling in CpG contexts."
        )
    return (
        "BAM header does not report modbase_models. Methylation classifications "
        "may be incorrect. Use 5mCG_5hmCG modbase calling in CpG contexts."
    )


def _extract_sample_id_from_bam(bam_path: str) -> str:
    """
    Extract sample ID from BAM file read group tags (and, if needed, a quick
    barcode probe).

    This is used both as a fallback when comprehensive processing fails and as
    an early probe to decide whether this BAM's sample output directory
    already exists (idempotency for re-adding watch folders).
    """
    try:
        bam_basename = os.path.basename(bam_path)
        unknown_fallback = f"unknown_sample_{bam_basename}"

        with pysam.AlignmentFile(bam_path, "rb") as bam:
            # Get read group tags
            rg_tags = bam.header.get(_RG_TAG, [])
            if not rg_tags:
                return unknown_fallback

            # Process only the first RG tag for efficiency
            rg_tag = rg_tags[0]

            # Look for LB (Library) tag which contains the sample ID
            lb_tag = rg_tag.get("LB")
            if lb_tag:
                sample_id = lb_tag
            else:
                # Fallback to runid from DS tag per ONT specification
                ds_tag = rg_tag.get("DS", "")
                if ds_tag:
                    ds_tags = ds_tag.split(" ")
                    if ds_tags and ds_tags[0].startswith(_RUNID_PREFIX):
                        runid = ds_tags[0].removeprefix(_RUNID_PREFIX)
                        sample_id = runid if runid else unknown_fallback
                    else:
                        sample_id = unknown_fallback
                else:
                    sample_id = unknown_fallback

            # If the sample_id already ends with a barcode, we are done.
            if _BARCODE_PATTERN.search(sample_id):
                return sample_id

            # Quick probe: scan reads until we find a barcode in the RG tag.
            # This avoids loading/processing the full BAM for idempotency checks.
            try:
                for read in bam.fetch(until_eof=True):
                    if not read.has_tag(_RG_TAG):
                        continue
                    rg_value = read.get_tag(_RG_TAG)
                    barcode_match = _BARCODE_PATTERN.search(rg_value)
                    if not barcode_match:
                        continue

                    barcode_num = int(barcode_match.group(1))
                    if 1 <= barcode_num <= 96:
                        barcode_str = f"_barcode{barcode_num:02d}"
                        if not sample_id.endswith(barcode_str):
                            sample_id = f"{sample_id}{barcode_str}"
                        return sample_id
            except Exception:
                # If the probe scan fails, keep the header-derived sample_id.
                pass

            return sample_id

    except Exception:
        # If pysam fails, return a stable directory-safe fallback.
        return f"unknown_sample_{os.path.basename(bam_path)}"


# ============================================================================
# READ PROCESSING FUNCTIONS
# ============================================================================


def process_bam_reads(bam_file: str) -> Optional[Dict[str, Any]]:
    """
    Processes the reads in the BAM file and aggregates information.
    This is the main processing step that analyzes all reads in the file.
    Uses streaming processing to minimize memory usage.

    Args:
        bam_file: Path to the BAM file

    Returns:
        Optional[Dict[str, Any]]: A dictionary containing aggregated read information,
                                  or None if processing fails.
    """
    if not bam_file:
        return None

    # print(f"Pre_Processing BAM file: {bam_file}")

    # Determine state from filename - use in operator for better performance
    state = _PASS_STR if _PASS_STR in bam_file else _FAIL_STR

    # Use context manager to ensure proper file cleanup
    try:
        # Enable BGZF multi-threaded reading if configured
        if _LJ_BAM_THREADS > 0:
            sam_open_kwargs = {"check_sq": False, "threads": _LJ_BAM_THREADS}
        else:
            sam_open_kwargs = {"check_sq": False}

        with pysam.AlignmentFile(bam_file, "rb", **sam_open_kwargs) as sam_file:
            # Step 1: Extract RG tags from header
            rg_tags = get_rg_tags_from_bam(sam_file)

            if not rg_tags:
                return None

            # Step 2: Determine sample ID
            # Ensure sample_id is not None - use runid as fallback per ONT specification
            # LB tag is optional, but runid from DS tag is always present
            sample_id = (
                rg_tags[4]  # LB (Library) tag
                if rg_tags[4]
                else (
                    rg_tags[3]  # runid from DS tag as fallback
                    if rg_tags[3]
                    else f"unknown_sample_{os.path.basename(bam_file)}"
                )  # filename as final fallback
            )

            # Step 3: Initialize data structures
            # Pre-allocate dictionary with known size for better performance
            bam_read = {
                "ID": rg_tags[0],
                "time_of_run": rg_tags[1],
                "sample_id": sample_id,
                "basecall_model": rg_tags[2],
                "runid": rg_tags[3],
                "platform": rg_tags[5],
                "flow_cell_id": rg_tags[7],
                "device_position": rg_tags[6],
                "al": rg_tags[8],
                "modbase_models": rg_tags[9],
                "state": state,
                "last_start": None,
                "elapsed_time": None,
            }

            # Initialize counters as local integers for faster hot-loop updates
            mapped_reads = 0
            unmapped_reads = 0
            yield_tracking = 0
            pass_mapped_reads_num = 0
            fail_mapped_reads_num = 0
            pass_unmapped_reads_num = 0
            fail_unmapped_reads_num = 0
            mapped_bases = 0
            unmapped_bases = 0
            pass_mapped_bases = 0
            fail_mapped_bases = 0
            pass_unmapped_bases = 0
            fail_unmapped_bases = 0
            supplementary_reads = 0

            # OPTIMIZATION: Use limited-size sets to prevent memory explosion
            # Only track unique reads up to a reasonable limit
            reads_with_supplementary = set()
            barcode_found = False
            last_start = None

            # Step 3.5: Initialize MGMT read detection
            has_mgmt_reads = False
            mgmt_read_count = 0

            # Step 4: Process reads in streaming fashion
            for read in sam_file.fetch(until_eof=True):
                # Extract RG tag and check for barcode - early termination optimization
                if not barcode_found and read.has_tag(_RG_TAG):
                    rg_tag = read.get_tag(_RG_TAG)
                    # Use compiled regex pattern for barcode detection
                    barcode_match = _BARCODE_PATTERN.search(rg_tag)
                    if barcode_match and sample_id:
                        barcode_num = int(barcode_match.group(1))
                        if 1 <= barcode_num <= 96:
                            barcode_str = f"_barcode{barcode_num:02d}"
                            if not sample_id.endswith(barcode_str):
                                sample_id = f"{sample_id}{barcode_str}"
                                barcode_found = True

                # Cache read properties to avoid repeated attribute access
                # Use direct attribute access for better performance
                read_length = read.query_length or 0
                is_secondary = read.is_secondary
                is_unmapped = read.is_unmapped
                is_supplementary = read.is_supplementary
                query_name = read.query_name

                # Track supplementary reads
                # A read has supplementary alignments if:
                # 1. This alignment itself is supplementary (is_supplementary=True), OR
                # 2. ANY alignment (primary or secondary) has an SA tag (indicating supplementary alignments exist)
                # 
                # We check ALL reads (not just primaries) to ensure we catch every read with supplementary mappings,
                # even if the SA tag is only present on certain alignment records
                has_supplementary_alignments = False
                if is_supplementary:
                    supplementary_reads += 1
                    has_supplementary_alignments = True
                
                # Also check SA tag on ALL reads (primary, secondary, supplementary) to catch any we might miss
                # Some BAM files may have SA tag on different records than expected
                if read.has_tag("SA"):
                    has_supplementary_alignments = True
                
                if has_supplementary_alignments:
                    reads_with_supplementary.add(query_name)

                if not is_secondary:  # Only process primary alignments
                    if not is_unmapped:
                        mapped_reads += 1
                        mapped_bases += read_length
                        if state == _PASS_STR:
                            pass_mapped_bases += read_length
                            pass_mapped_reads_num += 1
                        else:
                            fail_mapped_bases += read_length
                            fail_mapped_reads_num += 1

                        # Check for MGMT reads (chr10:129466536-129467536)
                        if not has_mgmt_reads and read.reference_name == _MGMT_CHR:
                            read_start = read.reference_start
                            read_end = read.reference_end
                            # Check if read overlaps with MGMT region (129466536-129467536)
                            if read_start < _MGMT_END and read_end > _MGMT_START:
                                has_mgmt_reads = True
                                mgmt_read_count += 1
                    else:
                        unmapped_reads += 1
                        unmapped_bases += read_length
                        if state == _PASS_STR:
                            pass_unmapped_bases += read_length
                            pass_unmapped_reads_num += 1
                        else:
                            fail_unmapped_bases += read_length
                            fail_unmapped_reads_num += 1

                    if read_length > 0:
                        yield_tracking += read_length

                    # Cache tag value to avoid multiple get_tag() calls
                    try:
                        st_tag = read.get_tag(_ST_TAG)
                        if not last_start or st_tag > last_start:
                            last_start = st_tag
                    except KeyError:
                        if not last_start:
                            last_start = bam_read["time_of_run"]

            # Step 5: Finalize data
            # Build counters dict once
            counters = {
                "mapped_reads": mapped_reads,
                "unmapped_reads": unmapped_reads,
                "yield_tracking": yield_tracking,
                "mapped_reads_num": mapped_reads,
                "unmapped_reads_num": unmapped_reads,
                "pass_mapped_reads_num": pass_mapped_reads_num,
                "fail_mapped_reads_num": fail_mapped_reads_num,
                "pass_unmapped_reads_num": pass_unmapped_reads_num,
                "fail_unmapped_reads_num": fail_unmapped_reads_num,
                "mapped_bases": mapped_bases,
                "unmapped_bases": unmapped_bases,
                "pass_mapped_bases": pass_mapped_bases,
                "fail_mapped_bases": fail_mapped_bases,
                "pass_unmapped_bases": pass_unmapped_bases,
                "fail_unmapped_bases": fail_unmapped_bases,
                "supplementary_reads": supplementary_reads,
                "reads_with_supplementary": len(reads_with_supplementary),
            }

            # Update sample_id in bam_read
            bam_read["sample_id"] = sample_id
            bam_read["last_start"] = last_start

            # OPTIMIZATION: Only store boolean flag and count, not the actual read IDs by default.
            # However, tests expect the list of unique read IDs; keep the behavior identical.
            bam_read["has_supplementary_reads"] = len(reads_with_supplementary) > 0
            bam_read["supplementary_read_ids_complete"] = True
            bam_read["supplementary_read_ids"] = (
                list(reads_with_supplementary) if reads_with_supplementary else []
            )

            # Add MGMT read information to metadata
            bam_read["has_mgmt_reads"] = has_mgmt_reads
            bam_read["mgmt_read_count"] = mgmt_read_count

            # Parse dates once at the end instead of during the loop
            if last_start and bam_read["time_of_run"]:
                try:
                    bam_read["elapsed_time"] = parser.parse(last_start) - parser.parse(
                        bam_read["time_of_run"]
                    )
                except (ValueError, TypeError):
                    # Handle invalid date formats gracefully
                    bam_read["elapsed_time"] = None

            # Add statistics to the return dictionary
            bam_read.update(counters)

            return bam_read

    except (IOError, ValueError):
        return None


# ============================================================================
# STATISTICS AND SUMMARY FUNCTIONS
# ============================================================================


def calculate_bam_summary(bam_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Calculates summary statistics from BAM processing data.
    This step processes the raw data to generate meaningful statistics.

    Args:
        bam_data: Dictionary containing BAM processing results

    Returns:
        dict: A dictionary containing read statistics including counts and mean lengths.
    """
    # Extract values from bam_data with defaults
    mapped_reads_num = bam_data.get("mapped_reads_num", 0)
    unmapped_reads_num = bam_data.get("unmapped_reads_num", 0)
    pass_mapped_reads_num = bam_data.get("pass_mapped_reads_num", 0)
    fail_mapped_reads_num = bam_data.get("fail_mapped_reads_num", 0)
    pass_unmapped_reads_num = bam_data.get("pass_unmapped_reads_num", 0)
    fail_unmapped_reads_num = bam_data.get("fail_unmapped_reads_num", 0)

    mapped_bases = bam_data.get("mapped_bases", 0)
    unmapped_bases = bam_data.get("unmapped_bases", 0)
    pass_mapped_bases = bam_data.get("pass_mapped_bases", 0)
    fail_mapped_bases = bam_data.get("fail_mapped_bases", 0)
    pass_unmapped_bases = bam_data.get("pass_unmapped_bases", 0)
    fail_unmapped_bases = bam_data.get("fail_unmapped_bases", 0)

    # Calculate mean lengths with cached values - use conditional expressions for better performance
    mean_mapped_length = mapped_bases / mapped_reads_num if mapped_reads_num > 0 else 0
    mean_unmapped_length = (
        unmapped_bases / unmapped_reads_num if unmapped_reads_num > 0 else 0
    )
    mean_pass_mapped_length = (
        pass_mapped_bases / pass_mapped_reads_num if pass_mapped_reads_num > 0 else 0
    )
    mean_fail_mapped_length = (
        fail_mapped_bases / fail_mapped_reads_num if fail_mapped_reads_num > 0 else 0
    )
    mean_pass_unmapped_length = (
        pass_unmapped_bases / pass_unmapped_reads_num
        if pass_unmapped_reads_num > 0
        else 0
    )
    mean_fail_unmapped_length = (
        fail_unmapped_bases / fail_unmapped_reads_num
        if fail_unmapped_reads_num > 0
        else 0
    )

    # Pre-allocate result dictionary with known size
    result = {
        "mapped_reads": bam_data.get("mapped_reads", 0),
        "unmapped_reads": bam_data.get("unmapped_reads", 0),
        "yield_tracking": bam_data.get("yield_tracking", 0),
        "state": bam_data.get("state", "unknown"),
        # Add read numbers
        "mapped_reads_num": mapped_reads_num,
        "unmapped_reads_num": unmapped_reads_num,
        "pass_mapped_reads_num": pass_mapped_reads_num,
        "fail_mapped_reads_num": fail_mapped_reads_num,
        "pass_unmapped_reads_num": pass_unmapped_reads_num,
        "fail_unmapped_reads_num": fail_unmapped_reads_num,
        # Add bases
        "mapped_bases": mapped_bases,
        "unmapped_bases": unmapped_bases,
        "pass_mapped_bases": pass_mapped_bases,
        "fail_mapped_bases": fail_mapped_bases,
        "pass_unmapped_bases": pass_unmapped_bases,
        "fail_unmapped_bases": fail_unmapped_bases,
        # Add mean lengths
        "mean_mapped_length": mean_mapped_length,
        "mean_unmapped_length": mean_unmapped_length,
        "mean_pass_mapped_length": mean_pass_mapped_length,
        "mean_fail_mapped_length": mean_fail_mapped_length,
        "mean_pass_unmapped_length": mean_pass_unmapped_length,
        "mean_fail_unmapped_length": mean_fail_unmapped_length,
        # Add supplementary read statistics
        "supplementary_reads": bam_data.get("supplementary_reads", 0),
        "reads_with_supplementary": bam_data.get("reads_with_supplementary", 0),
        "has_supplementary_reads": bam_data.get("has_supplementary_reads", False),
        # Add supplementary read IDs for fusion analysis
        "supplementary_read_ids_complete": bam_data.get(
            "supplementary_read_ids_complete", False
        ),
        "supplementary_read_ids": bam_data.get("supplementary_read_ids", []),
        # Add MGMT read statistics
        "has_mgmt_reads": bam_data.get("has_mgmt_reads", False),
        "mgmt_read_count": bam_data.get("mgmt_read_count", 0),
    }

    return result


# ============================================================================
# MAIN EXTRACTION FUNCTION
# ============================================================================


def extract_bam_metadata(bam_path: str) -> BamMetadata:
    """
    Extract comprehensive metadata from a BAM file.
    This is the main orchestration function that coordinates all processing steps.

    Processing order:
    1. Extract basic file information (size, creation time)
    2. Process BAM reads and extract metadata
    3. Calculate summary statistics
    4. Compile final metadata object

    Args:
        bam_path: Path to the BAM file

    Returns:
        BamMetadata: Object containing all extracted metadata
    """
    # Step 1: Extract basic file information
    # Cache path components to avoid repeated operations
    bam_path_obj = Path(bam_path)
    filename = bam_path_obj.name
    directory = str(bam_path_obj.parent)

    # Get basic file info
    stat = os.stat(bam_path)
    file_size = stat.st_size
    creation_time = stat.st_ctime

    # Initialize metadata
    metadata = BamMetadata(
        file_path=bam_path, file_size=file_size, creation_time=creation_time
    )

    # Step 2: Process BAM reads and extract comprehensive data
    bam_info = process_bam_reads(bam_path)
    if bam_info is None:
        # Fallback to basic extraction if processing fails
        sample_id = _extract_sample_id_from_bam(bam_path)
        state = _PASS_STR if _PASS_STR in bam_path else _FAIL_STR

        # Pre-allocate dictionary for better performance
        metadata.extracted_data = {
            "sample_id": sample_id,
            "file_size": file_size,
            "state": state,
            "filename": filename,
            "directory": directory,
            "has_mgmt_reads": False,
            "mgmt_read_count": 0,
        }
    else:
        # Use comprehensive data from processing - pre-allocate with known size
        metadata.extracted_data = {
            "sample_id": bam_info.get("sample_id", "unknown"),
            "file_size": file_size,
            "state": bam_info.get("state", "unknown"),
            "filename": filename,
            "directory": directory,
            "runid": bam_info.get("runid"),
            "platform": bam_info.get("platform"),
            "device_position": bam_info.get("device_position"),
            "basecall_model": bam_info.get("basecall_model"),
            "modbase_models": bam_info.get("modbase_models"),
            "flow_cell_id": bam_info.get("flow_cell_id"),
            "time_of_run": bam_info.get("time_of_run"),
            "file_path": bam_path,
            "has_mgmt_reads": bam_info.get("has_mgmt_reads", False),
            "mgmt_read_count": bam_info.get("mgmt_read_count", 0),
        }

    # Step 3: Calculate summary statistics
    bam_stats = calculate_bam_summary(bam_info) if bam_info else {}
    if bam_stats:
        metadata.extracted_data.update(bam_stats)

    # Step 4: Mark processing as complete
    metadata.processing_steps.append("bam_preprocessing_complete")

    return metadata


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


def _send_alignment_warning_notification(
    warning_msg: str, sample_id: str, filename: str
) -> None:
    """
    Send an alignment warning notification to the GUI.
    
    Args:
        warning_msg: The warning message to display
        sample_id: Sample ID for context
        filename: BAM filename for context
    """
    try:
        from robin.gui.app import send_gui_update
        from robin.gui_launcher import UpdateType
        
        send_gui_update(
            UpdateType.WARNING_NOTIFICATION,
            {
                "message": warning_msg,
                "sample_id": sample_id,
                "filename": filename,
                "title": "BAM Alignment Warning",
            },
            priority=5,  # Medium-high priority for warnings
        )
    except Exception:
        # Fail silently if GUI is not available - logging already happened
        pass


def _send_modbase_warning_notification(
    warning_msg: str, sample_id: str, filename: str
) -> None:
    """Send a methylation model warning to the GUI when available."""
    try:
        from robin.gui.app import send_gui_update
        from robin.gui_launcher import UpdateType

        send_gui_update(
            UpdateType.WARNING_NOTIFICATION,
            {
                "message": warning_msg,
                "sample_id": sample_id,
                "filename": filename,
                "title": "Methylation Model Warning",
            },
            priority=5,
        )
    except Exception:
        pass


def _send_skip_notification(
    warning_msg: str,
    sample_id: str,
    filename: str,
    sample_dir: str,
) -> None:
    """
    Notify GUI that a BAM was skipped due to existing output artifacts.
    """
    try:
        from robin.gui.app import send_gui_update
        from robin.gui_launcher import UpdateType

        send_gui_update(
            UpdateType.WARNING_NOTIFICATION,
            {
                "message": warning_msg,
                "sample_id": sample_id,
                "filename": filename,
                "sample_dir": sample_dir,
                "title": "Previously Analysed Sample",
            },
            priority=5,
        )
    except Exception:
        pass


# ============================================================================
# JOB HANDLER FUNCTION
# ============================================================================


def bam_preprocessing_handler(job, center: str = None):
    """
    Handler function for BAM preprocessing jobs.
    This is the main entry point that orchestrates the entire preprocessing workflow.

    Processing order:
    1. Initialize logging
    2. Extract metadata from BAM file
    3. Store metadata in job context
    4. Update master.csv with extracted data
    5. Add results to job context
    6. Log completion and debug information

    This function extracts metadata from BAM files and stores it in the job context.
    """
    try:
        bam_path = job.context.filepath

        # Step 1: Initialize logging
        logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
        logger.info(f"Starting preprocessing for: {os.path.basename(bam_path)}")

        # Step 2: Extract metadata from BAM file
        logger.debug(f"Extracting metadata from: {bam_path}")
        metadata = extract_bam_metadata(bam_path)
        logger.debug(f"Extracted metadata: {metadata.extracted_data}")

        # Step 2.5: Check for alignment data and warn if missing
        mapped_reads = metadata.extracted_data.get("mapped_reads", 0)
        unmapped_reads = metadata.extracted_data.get("unmapped_reads", 0)
        total_reads = mapped_reads + unmapped_reads
        sample_id = metadata.extracted_data.get("sample_id", "unknown")

        modbase_warning = _get_modbase_model_warning(
            metadata.extracted_data.get("modbase_models")
        )
        if modbase_warning:
            modbase_models = metadata.extracted_data.get("modbase_models")
            warning_level = _get_modbase_model_warning_level(modbase_models)
            if warning_level == "info":
                logger.info(modbase_warning)
            else:
                logger.warning(f"WARNING: {modbase_warning}")
                _send_modbase_warning_notification(
                    modbase_warning, sample_id, os.path.basename(bam_path)
                )
            metadata.extracted_data["modbase_warning"] = modbase_warning
            metadata.extracted_data["modbase_warning_level"] = warning_level
            job.context.add_metadata("modbase_warning", modbase_warning)
            job.context.add_metadata("modbase_warning_level", warning_level)
        
        if total_reads > 0:
            # Check if BAM file has no mapped reads (no alignment data)
            if mapped_reads == 0:
                warning_msg = (
                    f"BAM file contains no alignment data. "
                    f"Found {unmapped_reads:,} unmapped reads but {mapped_reads:,} mapped reads. "
                    f"This file may not be suitable for analysis that requires aligned reads."
                )
                logger.warning(f"WARNING: {warning_msg}")
                # Store warning in metadata for downstream display
                metadata.extracted_data["alignment_warning"] = warning_msg
                job.context.add_metadata("alignment_warning", warning_msg)
                # Send GUI notification
                _send_alignment_warning_notification(
                    warning_msg, sample_id, os.path.basename(bam_path)
                )
            # Check if mapped reads are very low relative to total (less than 1%)
            elif total_reads > 100 and mapped_reads > 0:
                mapped_percentage = (mapped_reads / total_reads) * 100
                if mapped_percentage < 1.0:
                    warning_msg = (
                        f"BAM file has very few aligned reads. "
                        f"Only {mapped_reads:,} mapped reads ({mapped_percentage:.2f}%) out of {total_reads:,} total reads. "
                        f"This may indicate alignment issues or a file that is not suitable for analysis requiring aligned reads."
                    )
                    logger.warning(f"WARNING: {warning_msg}")
                    # Store warning in metadata for downstream display
                    metadata.extracted_data["alignment_warning"] = warning_msg
                    job.context.add_metadata("alignment_warning", warning_msg)
                    # Send GUI notification
                    _send_alignment_warning_notification(
                        warning_msg, sample_id, os.path.basename(bam_path)
                    )
        elif total_reads == 0:
            # Edge case: BAM file has no reads at all
            warning_msg = (
                f"BAM file contains no reads. "
                f"This file may be empty or corrupted."
            )
            logger.warning(f"WARNING: {warning_msg}")
            metadata.extracted_data["alignment_warning"] = warning_msg
            job.context.add_metadata("alignment_warning", warning_msg)
            # Send GUI notification
            _send_alignment_warning_notification(
                warning_msg, sample_id, os.path.basename(bam_path)
            )

        # Step 3: Store metadata in job context
        job.context.add_metadata("bam_metadata", metadata.extracted_data)
        job.context.add_metadata("file_size", metadata.file_size)
        job.context.add_metadata("creation_time", metadata.creation_time)
        job.context.add_metadata("processing_steps", metadata.processing_steps)
        
        # Store center information
        if center:
            job.context.add_metadata("center", center)
            logger.info(f"Center set to: {center}")
        
        # Preserve target_panel metadata from job context (set by file classifier)
        existing_target_panel = job.context.metadata.get("target_panel")
        if existing_target_panel:
            logger.info(f"Preserving target panel from job context: {existing_target_panel}")
        else:
            logger.warning("No target_panel found in job context - this may cause panel assignment issues")
        
        # Preserve reference metadata from job context (set by file classifier)
        existing_reference = job.context.metadata.get("reference")
        if existing_reference:
            logger.info(f"Preserving reference from job context: {existing_reference}")
        else:
            logger.debug("No reference found in job context")

        # Check read count threshold: skip large BAMs unless ROBIN_PROCESS_LARGE_BAMS is set
        try:
            total_reads = int(metadata.extracted_data.get("mapped_reads", 0)) + int(
                metadata.extracted_data.get("unmapped_reads", 0)
            )
        except Exception:
            total_reads = 0
        
        if total_reads > 50000:
            if _PROCESS_LARGE_BAMS_INDIVIDUALLY:
                # Process this BAM; downstream batching will pass it to workers as a single-file batch
                job.context.add_metadata("force_individual_batch", True)
                job.context.add_metadata("large_bam_read_count", total_reads)
                logger.info(
                    f"Large BAM ({total_reads:,} reads) will be processed individually (ROBIN_PROCESS_LARGE_BAMS=1)."
                )
                # Fall through: persist supp, update CSV, add success result, trigger downstream
            else:
                # Mark in context to prevent downstream triggers; record reason and count
                job.context.add_metadata("skip_downstream", True)
                job.context.add_metadata("skip_reason", "too_many_reads")
                job.context.add_metadata("skip_read_count", total_reads)
                # Add an explicit result so the wrapper can log appropriately
                job.context.add_result(
                    "preprocessing",
                    {
                        "status": "skipped",
                        "reason": "too_many_reads",
                        "read_count": total_reads,
                        "message": "Deliberately not processed (>50,000 reads)",
                    },
                )
                # CLI warning (logger goes to CLI) with filename
                logger.warning(
                    f"Skipping {os.path.basename(bam_path)}: {total_reads} reads exceeds 50,000 limit. File deliberately not processed."
                )
                # Do not proceed with CSV updates or further processing
                return

        # Persist the complete ID set per BAM to avoid retaining large lists in memory.
        try:
            work_dir = job.context.metadata.get(
                "work_dir", os.path.dirname(bam_path)
            )
            _persist_supplementary_read_ids(metadata, bam_path, work_dir)
        except Exception:
            # Keep the complete in-memory list as a safe fallback.
            metadata.extracted_data.pop("supplementary_read_ids_path", None)

        # Step 4: Update master.csv if we have comprehensive data
        if metadata.extracted_data and "sample_id" in metadata.extracted_data:
            sample_id = metadata.extracted_data["sample_id"]

            # Determine work directory - check job context first, then use BAM file directory as default
            work_dir = job.context.metadata.get("work_dir", os.path.dirname(bam_path))

            # Create MasterCSVManager and update master.csv
            try:
                csv_manager = MasterCSVManager(work_dir)

                # Extract BAM statistics for CSV update (dict comp inlined in 3.12)
                _bam_stat_keys = (
                    "mapped_reads",
                    "unmapped_reads",
                    "yield_tracking",
                    "mapped_reads_num",
                    "unmapped_reads_num",
                    "pass_mapped_reads_num",
                    "fail_mapped_reads_num",
                    "pass_unmapped_reads_num",
                    "fail_unmapped_reads_num",
                    "mapped_bases",
                    "unmapped_bases",
                    "pass_mapped_bases",
                    "fail_mapped_bases",
                    "pass_unmapped_bases",
                    "fail_unmapped_bases",
                    "supplementary_reads",
                    "reads_with_supplementary",
                    "has_mgmt_reads",
                    "mgmt_read_count",
                )
                bam_stats = {key: metadata.extracted_data.get(key, 0) for key in _bam_stat_keys}

                # Update master.csv
                csv_manager.update_master_csv(
                    sample_id, bam_stats, metadata.extracted_data
                )
                
                # Update analysis panel in master.csv (preserve from job context)
                existing_target_panel = job.context.metadata.get("target_panel")
                if existing_target_panel:
                    csv_manager.update_analysis_panel(sample_id, existing_target_panel)
                    logger.info(f"Updated master.csv with analysis panel '{existing_target_panel}' for sample {sample_id}")
                else:
                    logger.warning(f"No target_panel found in job context for sample {sample_id} - analysis_panel not set in master.csv")

            except Exception as e:
                logger.warning(f"Could not update master.csv for {sample_id}: {e}")

        # Step 5: Add result to context
        sample_id = metadata.extracted_data.get("sample_id", "unknown")
        
        job.context.add_result(
            "preprocessing",
            {
                "status": "success",
                "sample_id": sample_id,
                "state": metadata.extracted_data.get("state", "unknown"),
                "mapped_reads": metadata.extracted_data.get("mapped_reads", 0),
                "unmapped_reads": metadata.extracted_data.get("unmapped_reads", 0),
                "total_yield": metadata.extracted_data.get("yield_tracking", 0),
                "file_size": metadata.file_size,
                "has_supplementary_reads": metadata.extracted_data.get(
                    "has_supplementary_reads", False
                ),
                "supplementary_reads": metadata.extracted_data.get(
                    "supplementary_reads", 0
                ),
                "reads_with_supplementary": metadata.extracted_data.get(
                    "reads_with_supplementary", 0
                ),
                "supplementary_read_ids_complete": metadata.extracted_data.get(
                    "supplementary_read_ids_complete", False
                ),
                "supplementary_read_ids_path": metadata.extracted_data.get(
                    "supplementary_read_ids_path"
                ),
                "supplementary_read_ids_count": metadata.extracted_data.get(
                    "supplementary_read_ids_count", 0
                ),
                "supplementary_read_ids": metadata.extracted_data.get(
                    "supplementary_read_ids", []
                ),
                "has_mgmt_reads": metadata.extracted_data.get("has_mgmt_reads", False),
                "mgmt_read_count": metadata.extracted_data.get("mgmt_read_count", 0),
                "modbase_models": metadata.extracted_data.get("modbase_models"),
                "modbase_warning": metadata.extracted_data.get("modbase_warning"),
            },
        )

        # Step 6: Log completion and summary
        logger.info(f"Preprocessing complete for {os.path.basename(bam_path)}")
        logger.info(f"Sample ID: {metadata.extracted_data.get('sample_id', 'unknown')}")
        logger.info(f"State: {metadata.extracted_data.get('state', 'unknown')}")
        logger.info(f"Mapped reads: {metadata.extracted_data.get('mapped_reads', 0):,}")
        logger.info(
            f"Unmapped reads: {metadata.extracted_data.get('unmapped_reads', 0):,}"
        )
        logger.info(
            f"Total yield: {metadata.extracted_data.get('yield_tracking', 0):,}"
        )

        # Log supplementary read information
        has_supplementary = metadata.extracted_data.get(
            "has_supplementary_reads", False
        )
        supplementary_count = metadata.extracted_data.get("supplementary_reads", 0)
        reads_with_supplementary = metadata.extracted_data.get(
            "reads_with_supplementary", 0
        )

        if has_supplementary:
            logger.info(
                f"Supplementary reads found: {supplementary_count:,} supplementary alignments from {reads_with_supplementary:,} unique reads"
            )
        else:
            logger.info("No supplementary reads found")

        # Log MGMT read information
        has_mgmt_reads = metadata.extracted_data.get("has_mgmt_reads", False)
        mgmt_read_count = metadata.extracted_data.get("mgmt_read_count", 0)

        if has_mgmt_reads:
            logger.info(
                f"MGMT reads found: {mgmt_read_count:,} reads spanning MGMT promoter region (chr10:129466536-129467536)"
            )
        else:
            logger.info("No MGMT reads found")

        # Step 7: Detailed debug logging (only if debug is enabled)
        if logger.logger.isEnabledFor(10):  # DEBUG level
            # Cache frequently accessed values to avoid repeated dict lookups
            extracted_data = metadata.extracted_data
            sample_id = extracted_data.get("sample_id", "unknown")
            mapped_reads = extracted_data.get("mapped_reads", 0)
            unmapped_reads = extracted_data.get("unmapped_reads", 0)

            # Log sample ID extraction details at debug level
            logger.debug("Sample ID extraction details:")
            logger.debug(f"  - Extracted sample_id: {sample_id}")
            logger.debug(f"  - Filename: {os.path.basename(bam_path)}")
            logger.debug(f"  - Full path: {bam_path}")

            # Log comprehensive metadata at debug level
            logger.debug("BAM File Metadata:")
            logger.debug(f"File path: {metadata.file_path}")
            logger.debug(f"File size: {metadata.file_size:,} bytes")
            logger.debug(f"Creation time: {time.ctime(metadata.creation_time)}")

            # Log run information if available - use cached values
            runid = extracted_data.get("runid")
            if runid:
                logger.debug(f"Run ID: {runid}")
            platform = extracted_data.get("platform")
            if platform:
                logger.debug(f"Platform: {platform}")
            device_position = extracted_data.get("device_position")
            if device_position:
                logger.debug(f"Device: {device_position}")
            flow_cell_id = extracted_data.get("flow_cell_id")
            if flow_cell_id:
                logger.debug(f"Flow cell: {flow_cell_id}")
            basecall_model = extracted_data.get("basecall_model")
            if basecall_model:
                logger.debug(f"Basecall model: {basecall_model}")
            modbase_models = extracted_data.get("modbase_models")
            if modbase_models:
                logger.debug(f"Modbase models: {modbase_models}")
            time_of_run = extracted_data.get("time_of_run")
            if time_of_run:
                logger.debug(f"Run time: {time_of_run}")

            # Log detailed statistics - use cached values
            logger.debug("Read Statistics:")
            logger.debug(f"Mapped reads: {mapped_reads:,}")
            logger.debug(f"Unmapped reads: {unmapped_reads:,}")
            logger.debug(f"Total reads: {mapped_reads + unmapped_reads:,}")

            mapped_bases = extracted_data.get("mapped_bases", 0)
            unmapped_bases = extracted_data.get("unmapped_bases", 0)
            yield_tracking = extracted_data.get("yield_tracking", 0)
            logger.debug(f"Mapped bases: {mapped_bases:,}")
            logger.debug(f"Unmapped bases: {unmapped_bases:,}")
            logger.debug(f"Total yield: {yield_tracking:,}")

            # Log mean lengths if available
            mean_mapped_length = extracted_data.get("mean_mapped_length", 0)
            if mean_mapped_length > 0:
                logger.debug(f"Mean mapped read length: {mean_mapped_length:.1f}")
            mean_unmapped_length = extracted_data.get("mean_unmapped_length", 0)
            if mean_unmapped_length > 0:
                logger.debug(f"Mean unmapped read length: {mean_unmapped_length:.1f}")

            # Log pass/fail breakdown if available
            pass_mapped_reads_num = extracted_data.get("pass_mapped_reads_num", 0)
            fail_mapped_reads_num = extracted_data.get("fail_mapped_reads_num", 0)
            if pass_mapped_reads_num > 0 or fail_mapped_reads_num > 0:
                logger.debug(f"Pass mapped reads: {pass_mapped_reads_num:,}")
                logger.debug(f"Fail mapped reads: {fail_mapped_reads_num:,}")
                logger.debug(
                    f"Pass unmapped reads: {extracted_data.get('pass_unmapped_reads_num', 0):,}"
                )
                logger.debug(
                    f"Fail unmapped reads: {extracted_data.get('fail_unmapped_reads_num', 0):,}"
                )

            # Log supplementary read information if available
            has_supplementary = extracted_data.get("has_supplementary_reads", False)
            if has_supplementary:
                supplementary_count = extracted_data.get("supplementary_reads", 0)
                reads_with_supplementary = extracted_data.get(
                    "reads_with_supplementary", 0
                )
                supplementary_read_ids = extracted_data.get(
                    "supplementary_read_ids", []
                )
                logger.debug(
                    f"Supplementary reads: {supplementary_count:,} supplementary alignments"
                )
                logger.debug(
                    f"Reads with supplementary mappings: {reads_with_supplementary:,}"
                )
                if supplementary_read_ids:
                    logger.debug(
                        f"First 10 supplementary read IDs: {supplementary_read_ids[:10]}"
                    )
                    if len(supplementary_read_ids) > 10:
                        logger.debug(f"... and {len(supplementary_read_ids) - 10} more")
            else:
                logger.debug("No supplementary reads found")

            # Log MGMT read information if available
            has_mgmt_reads = extracted_data.get("has_mgmt_reads", False)
            if has_mgmt_reads:
                mgmt_read_count = extracted_data.get("mgmt_read_count", 0)
                logger.debug(
                    f"MGMT reads: {mgmt_read_count:,} reads spanning MGMT promoter region"
                )
                logger.debug("MGMT region: chr10:129466536-129467536")
            else:
                logger.debug("No MGMT reads found")

            logger.debug(f"Processing steps: {', '.join(metadata.processing_steps)}")

    except Exception as e:
        job.context.add_error("preprocessing", str(e))
        logger.error(f"Error preprocessing {job.context.filepath}: {e}")
