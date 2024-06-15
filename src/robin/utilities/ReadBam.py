"""
Original code written by: Thomas Murray and taken from https://github.com/tom-murray98/minFQ_BAM/blob/master/BAM_RG_tag_extractor.py
"""

import pysam
import os
from typing import Optional, Tuple, Dict, Any, Generator, Set
import logging

# Configure logging
#logger = logging.getLogger(__name__)


class ReadBam:
    """
    A class to handle reading and extracting information from BAM files.

    Args:
        bam_file (Optional[str]): A BAM format alignment file.

    Attributes:
        sam_file (Optional[pysam.AlignmentFile]): The SAM/BAM file object.
        bam_file (Optional[str]): The BAM file path.
        mapped_reads (int): The number of mapped reads.
        unmapped_reads (int): The number of unmapped reads.
        state (str): The state of the BAM file, either 'pass' or 'fail'.
    """

    def __init__(self, bam_file: Optional[str] = None) -> None:
        self.sam_file: Optional[pysam.AlignmentFile] = None
        self.bam_file: Optional[str] = bam_file
        self.mapped_reads: int = 0
        self.unmapped_reads: int = 0
        self.yield_tracking: int = 0
        if self.bam_file and "pass" in self.bam_file:
            self.state = "pass"
        else:
            self.state = "fail"

    def summary(self) -> Dict[str, Any]:
        """
        Returns a summary of the BAM file.

        Returns:
            A dictionary containing the number of mapped and unmapped reads.
        """
        return {
            "mapped_reads": self.mapped_reads,
            "unmapped_reads": self.unmapped_reads,
            "yield_tracking": self.yield_tracking,
            "state": self.state,
        }

    def get_rg_tags(
        self,
    ) -> Optional[
        Tuple[
            Optional[str],
            Optional[str],
            Optional[str],
            Optional[str],
            Optional[str],
            Optional[str],
            Optional[str],
            Optional[str],
            Optional[str],
        ]
    ]:
        """
        Extracts RG (Read Group) tags from the BAM file header.

        Returns:
            A tuple containing RG tag information or None if no RG tags are found.
        """
        if self.sam_file:
            rg_tags = self.sam_file.header.get("RG", [])
            if not rg_tags:
                #logger.warning("This BAM file does not contain an @RG field")
                return None

            for rg_tag in rg_tags:
                id_tag = rg_tag.get("ID", None)
                dt_tag = rg_tag.get("DT", None)
                ds_tag = rg_tag.get("DS", None)
                ds_tags = ds_tag.split(" ") if ds_tag else ["", ""]
                basecall_model_tag = (
                    ds_tags[1].replace("basecall_model=", "")
                    if len(ds_tags) > 1
                    else None
                )
                runid_tag = (
                    ds_tags[0].replace("runid=", "") if len(ds_tags) > 0 else None
                )
                lb_tag = rg_tag.get("LB", None)
                pl_tag = rg_tag.get("PL", None)
                pm_tag = rg_tag.get("PM", None)
                pu_tag = rg_tag.get("PU", None)
                al_tag = rg_tag.get("al", None)

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
                )
        return None

    def read_bam(
        self,
    ) -> Generator[
        Tuple[
            Optional[str],
            Optional[str],
            Optional[str],
            Optional[str],
            Optional[str],
            Optional[str],
            Optional[str],
            Optional[str],
            Optional[str],
        ],
        None,
        None,
    ]:
        """
        Reads the BAM file and extracts RG tags.

        Yields:
            A tuple containing RG tag information for each read in the BAM file.
        """
        if self.bam_file:
            if not os.path.isfile(f"{self.bam_file}.bai"):
                #logger.warning(
                #    f"Index file for {self.bam_file} does not exist. Generating index file."
                #)
                pysam.index(self.bam_file)
            self.sam_file = pysam.AlignmentFile(self.bam_file, "rb", check_sq=False)

            rg_tags = self.get_rg_tags()
            if rg_tags:
                yield rg_tags

    def process_reads(self) -> Optional[Dict[str, Any]]:
        """
        Processes the reads in the BAM file and aggregates information.

        Returns:
            A dictionary containing aggregated read information.
        """
        for (
            id_tag,
            dt_tag,
            basecall_model_tag,
            runid_tag,
            lb_tag,
            pl_tag,
            pm_tag,
            pu_tag,
            al_tag,
        ) in self.read_bam():
            bam_read = {
                "ID": id_tag,
                "time_of_run": dt_tag,
                "sample_id": lb_tag,
                "basecall_model": basecall_model_tag,
                "runid": runid_tag,
                "platform": pl_tag,
                "flow_cell_id": pu_tag,
                "device_position": pm_tag,
                "al": al_tag,
                "state": self.state,
            }

            filtered_reads_count = 0
            readset: Set[str] = set()
            if self.sam_file:
                for read in self.sam_file.fetch():
                    # Check if read is not unmapped and not a secondary alignment
                    if not read.is_unmapped and not read.is_secondary:
                        filtered_reads_count += 1
                        if read.query_name not in readset:
                            readset.add(read.query_name)
                            if (
                                read.infer_query_length()
                                and read.infer_query_length() > 0
                            ):
                                self.yield_tracking += read.infer_query_length()

            self.mapped_reads = len(readset)
            self.unmapped_reads = self.sam_file.mapped if self.sam_file else 0

            #logger.debug(f"Mapped reads: {filtered_reads_count}")
            #logger.debug(f"Total reads: {filtered_reads_count + self.unmapped_reads}")

            return bam_read
        return None
