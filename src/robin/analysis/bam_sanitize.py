"""Skip BAM records that htslib rejects (CIGAR vs query length)."""

from __future__ import annotations

import logging
import os
from typing import Iterator

import pysam

logger = logging.getLogger("robin.bam")

_CORRUPT_ALIGNMENT_MARKERS = (
    "cigar and query sequence lengths differ",
    "error -4 while reading file",
    "error while reading file",
    "error reading from input file",
    "failed to create index",
)


def is_corrupt_alignment_error(exc: BaseException) -> bool:
    """True when htslib/samtools aborted on a malformed alignment record."""
    text = " ".join(str(exc).lower().split())
    if any(marker in text for marker in _CORRUPT_ALIGNMENT_MARKERS):
        return True
    # IteratorRowRegion: "error while reading file\\nb'/path.bam': -4"
    return "-4" in text and "reading file" in text


def is_samtools_tool_error(exc: BaseException) -> bool:
    """True for pysam SamtoolsError (coverage/bedcov/index) failures."""
    if type(exc).__name__ == "SamtoolsError":
        return True
    return "samtools returned with error" in str(exc).lower()


def is_missing_index_error(exc: BaseException) -> bool:
    """True when pysam fetch was used on an unindexed BAM."""
    return "without index" in str(exc).lower()


def bam_index_exists(bam_path: str) -> bool:
    """Return True if a sibling .bai or .csi index is present."""
    return os.path.exists(f"{bam_path}.bai") or os.path.exists(f"{bam_path}.csi")


def iter_alignments(
    alignment_file: pysam.AlignmentFile,
    *fetch_args,
    **fetch_kwargs,
) -> Iterator[pysam.AlignedSegment]:
    """
    Yield alignments, skipping records htslib rejects.

    ``bam_read1`` consumes a corrupt record before returning -4, so the next
    call continues at the following alignment. This does not rewrite the BAM.
    """
    if fetch_args or fetch_kwargs:
        iterator = alignment_file.fetch(*fetch_args, **fetch_kwargs)
    else:
        iterator = alignment_file.fetch(until_eof=True)

    skipped = 0
    while True:
        try:
            yield next(iterator)
        except StopIteration:
            break
        except Exception as exc:
            if not is_corrupt_alignment_error(exc):
                raise
            skipped += 1
    if skipped:
        logger.warning(
            "Skipped %d alignment(s) with CIGAR/query-length mismatch",
            skipped,
        )
