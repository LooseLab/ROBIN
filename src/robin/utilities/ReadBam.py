import pysam
import os
from typing import Optional, Tuple, Dict, Any, Generator, Set
from dataclasses import dataclass, field, asdict
import logging
from dateutil import parser
import re

# Create a logger for this module
logger = logging.getLogger(__name__)


def configure_logging(level=None):
    """
    Configure the logging for this module.

    Args:
        level (int, optional): The logging level to set. If None, uses the root logger's level.
    """
    # If no level specified, use the root logger's level
    if level is None:
        level = logging.getLogger().getEffectiveLevel()
    
    logger.setLevel(level)

    # Create a console handler if no handlers are set
    if not logger.handlers:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(level)
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)


@dataclass
class BamRead:
    ID: Optional[str] = None
    time_of_run: Optional[str] = None
    sample_id: Optional[str] = None
    basecall_model: Optional[str] = None
    runid: Optional[str] = None
    platform: Optional[str] = None
    flow_cell_id: Optional[str] = None
    device_position: Optional[str] = None
    al: Optional[str] = None
    state: str = "fail"
    last_start: Optional[int] = None
    elapsed_time: Optional[int] = None


@dataclass
class ReadBam:
    bam_file: Optional[str] = None
    sam_file: Optional[pysam.AlignmentFile] = None
    mapped_reads: int = 0
    unmapped_reads: int = 0
    yield_tracking: int = 0
    state: str = field(init=False)
    # New counters for read numbers
    mapped_reads_num: int = 0
    unmapped_reads_num: int = 0
    pass_mapped_reads_num: int = 0
    fail_mapped_reads_num: int = 0
    pass_unmapped_reads_num: int = 0
    fail_unmapped_reads_num: int = 0
    # New counters for bases
    mapped_bases: int = 0
    unmapped_bases: int = 0
    pass_mapped_bases: int = 0
    fail_mapped_bases: int = 0
    pass_unmapped_bases: int = 0
    fail_unmapped_bases: int = 0

    def __post_init__(self):
        self.state = "pass" if self.bam_file and "pass" in self.bam_file else "fail"

    def summary(self) -> Dict[str, Any]:
        """
        Returns a summary of the BAM file.

        Returns:
            dict: A dictionary containing read statistics including counts and mean lengths.
        """
        # Calculate mean lengths
        mean_mapped_length = self.mapped_bases / self.mapped_reads_num if self.mapped_reads_num > 0 else 0
        mean_unmapped_length = self.unmapped_bases / self.unmapped_reads_num if self.unmapped_reads_num > 0 else 0
        mean_pass_mapped_length = self.pass_mapped_bases / self.pass_mapped_reads_num if self.pass_mapped_reads_num > 0 else 0
        mean_fail_mapped_length = self.fail_mapped_bases / self.fail_mapped_reads_num if self.fail_mapped_reads_num > 0 else 0
        mean_pass_unmapped_length = self.pass_unmapped_bases / self.pass_unmapped_reads_num if self.pass_unmapped_reads_num > 0 else 0
        mean_fail_unmapped_length = self.fail_unmapped_bases / self.fail_unmapped_reads_num if self.fail_unmapped_reads_num > 0 else 0

        return {
            "mapped_reads": self.mapped_reads,
            "unmapped_reads": self.unmapped_reads,
            "yield_tracking": self.yield_tracking,
            "state": self.state,
            # Add read numbers
            "mapped_reads_num": self.mapped_reads_num,
            "unmapped_reads_num": self.unmapped_reads_num,
            "pass_mapped_reads_num": self.pass_mapped_reads_num,
            "fail_mapped_reads_num": self.fail_mapped_reads_num,
            "pass_unmapped_reads_num": self.pass_unmapped_reads_num,
            "fail_unmapped_reads_num": self.fail_unmapped_reads_num,
            # Add bases
            "mapped_bases": self.mapped_bases,
            "unmapped_bases": self.unmapped_bases,
            "pass_mapped_bases": self.pass_mapped_bases,
            "fail_mapped_bases": self.fail_mapped_bases,
            "pass_unmapped_bases": self.pass_unmapped_bases,
            "fail_unmapped_bases": self.fail_unmapped_bases,
            # Add mean lengths
            "mean_mapped_length": mean_mapped_length,
            "mean_unmapped_length": mean_unmapped_length,
            "mean_pass_mapped_length": mean_pass_mapped_length,
            "mean_fail_mapped_length": mean_fail_mapped_length,
            "mean_pass_unmapped_length": mean_pass_unmapped_length,
            "mean_fail_unmapped_length": mean_fail_unmapped_length
        }

    def get_rg_tags(self) -> Optional[Tuple[Optional[str], ...]]:
        """
        Extracts RG (Read Group) tags from the BAM file header.

        Returns:
            Optional[Tuple[Optional[str], ...]]: A tuple containing RG tag information
                                                 or None if no RG tags are found.
        """
        if not self.sam_file:
            logger.warning("SAM file is not initialized")
            return None

        rg_tags = self.sam_file.header.get("RG", [])
        if not rg_tags:
            logger.warning("This BAM file does not contain an @RG field")
            return None

        for rg_tag in rg_tags:
            id_tag = rg_tag.get("ID")
            dt_tag = rg_tag.get("DT")
            ds_tag = rg_tag.get("DS", "")
            ds_tags = ds_tag.split(" ")
            basecall_model_tag = (
                ds_tags[1].replace("basecall_model=", "") if len(ds_tags) > 1 else None
            )
            runid_tag = ds_tags[0].replace("runid=", "") if ds_tags else None
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
            )

        return None

    def read_bam(self) -> Generator[Tuple[Optional[str], ...], None, None]:
        """
        Reads the BAM file and extracts RG tags.

        Yields:
            Tuple[Optional[str], ...]: A tuple containing RG tag information for each read in the BAM file.
        """
        if not self.bam_file:
            logger.error("BAM file path is not set")
            return

        logger.debug(f"Read Bam {self.bam_file}")

        index_file = f"{self.bam_file}.bai"
        if not os.path.isfile(index_file):
            logger.warning(
                f"Index file for {self.bam_file} does not exist. Generating index file."
            )
            pysam.index(self.bam_file)

        try:
            self.sam_file = pysam.AlignmentFile(self.bam_file, "rb", check_sq=False)
        except (IOError, ValueError) as e:
            logger.error(f"Error opening BAM file: {e}")
            return

        rg_tags = self.get_rg_tags()
        if rg_tags:
            yield rg_tags

    def process_reads(self) -> Optional[Dict[str, Any]]:
        """
        Processes the reads in the BAM file and aggregates information.
        Also extracts and prints RG tags for each read.

        Returns:
            Optional[Dict[str, Any]]: A dictionary containing aggregated read information,
                                      or None if processing fails.
        """
        for rg_tags in self.read_bam():
            bam_read = BamRead(
                ID=rg_tags[0],
                time_of_run=rg_tags[1],
                sample_id=rg_tags[4],
                basecall_model=rg_tags[2],
                runid=rg_tags[3],
                platform=rg_tags[5],
                flow_cell_id=rg_tags[7],
                device_position=rg_tags[6],
                al=rg_tags[8],
                state=self.state,
                last_start=None,
                elapsed_time=None,
            )

            if not self.sam_file:
                logger.error("SAM file is not initialized")
                return None

            readset: Set[str] = set()
            mapped_readset: Set[str] = set()
            unmapped_readset: Set[str] = set()
            barcode_found = False

            for read in self.sam_file.fetch(until_eof=True):
                # Extract RG tag and check for barcode
                if read.has_tag('RG') and not barcode_found:
                    rg_tag = read.get_tag('RG')
                    # Check if RG tag ends with _barcodeNN where NN is 01-96
                    barcode_match = re.search(r'_barcode(\d{1,2})$', rg_tag)
                    if barcode_match and bam_read.sample_id:
                        barcode_num = int(barcode_match.group(1))
                        if 1 <= barcode_num <= 96:
                            barcode_str = f"_barcode{barcode_num:02d}"
                            if not bam_read.sample_id.endswith(barcode_str):
                                bam_read.sample_id = f"{bam_read.sample_id}{barcode_str}"
                                barcode_found = True
                                logger.info(f"Updated sample_id to: {bam_read.sample_id}")

                read_length = read.query_length if read.query_length else 0
                
                if not read.is_secondary:  # Only process primary alignments
                    if read.query_name not in readset:
                        readset.add(read.query_name)
                        
                        if not read.is_unmapped:
                            mapped_readset.add(read.query_name)
                            self.mapped_bases += read_length
                            if self.state == "pass":
                                self.pass_mapped_bases += read_length
                                self.pass_mapped_reads_num += 1
                            else:
                                self.fail_mapped_bases += read_length
                                self.fail_mapped_reads_num += 1
                        else:
                            unmapped_readset.add(read.query_name)
                            self.unmapped_bases += read_length
                            if self.state == "pass":
                                self.pass_unmapped_bases += read_length
                                self.pass_unmapped_reads_num += 1
                            else:
                                self.fail_unmapped_bases += read_length
                                self.fail_unmapped_reads_num += 1
                        
                        if read_length > 0:
                            self.yield_tracking += read_length

                if not bam_read.last_start:
                    bam_read.last_start = read.get_tag("st")
                if read.get_tag("st") > bam_read.last_start:
                    bam_read.last_start = read.get_tag("st")

            self.mapped_reads = len(mapped_readset)
            self.unmapped_reads = len(unmapped_readset)
            self.mapped_reads_num = len(mapped_readset)
            self.unmapped_reads_num = len(unmapped_readset)

            logger.info(f"Mapped reads: {self.mapped_reads}")
            logger.info(f"Total reads: {self.mapped_reads + self.unmapped_reads}")
            bam_read.elapsed_time = parser.parse(bam_read.last_start) - parser.parse(
                bam_read.time_of_run
            )

            return asdict(bam_read)

        return None

    def get_last_read(self) -> Optional[BamRead]:
        """
        Looks at the last read in the BAM file and returns a BamRead object.

        Returns:
            Optional[BamRead]: A BamRead object containing information about the last read,
                               or None if the file cannot be processed.
        """
        for rg_tags in self.read_bam():
            if not self.sam_file:
                logger.error("SAM file is not initialized")
                return None

            last_read = None
            for read in self.sam_file.fetch(until_eof=True):
                last_read = read

            if last_read is None:
                logger.warning("No reads found in the BAM file")
                return None

            bam_read = BamRead(
                ID=rg_tags[0],
                time_of_run=rg_tags[1],
                sample_id=rg_tags[4],
                basecall_model=rg_tags[2],
                runid=rg_tags[3],
                platform=rg_tags[5],
                flow_cell_id=rg_tags[7],
                device_position=rg_tags[6],
                al=rg_tags[8],
                state=self.state,
                last_start=last_read.get_tag("st") if last_read.has_tag("st") else None,
            )
            return asdict(bam_read)


# If this script is run directly, configure logging with default settings
if __name__ == "__main__":
    configure_logging()
