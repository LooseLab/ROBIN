from multiprocessing import Manager
from queue import Empty, Queue
from typing import Callable, Generator, List, Optional
from nicegui import background_tasks, run
import asyncio
import os
import logging
import subprocess
import time
import pickle
from datetime import datetime
from pathlib import Path
import pandas as pd
import polars as pl
import pyranges as pr
import pysam
from alive_progress import alive_bar
from robin.utilities.ReadBam import ReadBam
from robin.utilities.mnp_flex import APIClient as MnpFlexClient
from robin import resources
from robin.subpages.RandomForest_object import load_modkit_data


class Worker:

    def __init__(self) -> None:
        self._queue: Queue
        self.progress: float = 0.0
        self.is_running: bool = False
        self._create_queue()

    async def run(self, func: Callable[..., Generator[float, None, None]]) -> None:
        background_tasks.create(run.cpu_bound(self._run_generator, func, self._queue))
        background_tasks.create(self._consume_queue())

    @staticmethod
    def _run_generator(
        func: Callable[..., Generator[float, None, None]], queue: Queue
    ) -> None:
        for progress in func():
            queue.put({"progress": progress})
        queue.put({"progress": 1.0})

    def _create_queue(self) -> None:
        self._queue = Manager().Queue()

    async def _consume_queue(self) -> None:
        self.is_running = True
        self.progress = 0.0
        while self.progress < 1.0:
            try:
                msg = self._queue.get_nowait()
                self.progress = msg["progress"]
            except Empty:
                await asyncio.sleep(0.1)
        self.is_running = False


def merge_modkit_files(
    new_files: List[str],
    existing_file: str,
    output_file: str,
    filter_bed_file: str,
    sample_id: str,
    output_dir: str,
    mnpflex_config: Optional[dict],
) -> None:
    """
    Merge modkit files with improved caching and error handling.

    Args:
        new_files (List[str]): List of new modkit files to merge
        existing_file (str): Path to existing parquet file
        output_file (str): Path to output parquet file
        filter_bed_file (str): Path to BED file for filtering
        sample_id (str): Sample ID for organizing output files
        output_dir (str): Base output directory
        mnpflex_config (Optional[dict]): Configuration for MNP-FLEX integration
    """
    # Create sample-specific output directory
    sample_output_dir = os.path.join(output_dir, sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)

    # Define schema
    cols = [
        "chrom",
        "chromStart",
        "chromEnd",
        "mod_code",
        "score_bed",
        "strand",
        "thickStart",
        "thickEnd",
        "color",
        "valid_cov",
        "percent_modified",
        "n_mod",
        "n_canonical",
        "n_othermod",
        "n_delete",
        "n_fail",
        "n_diff",
        "n_nocall",
    ]
    categorical_cols = ["chrom", "mod_code", "strand", "color"]
    int_cols = ["thickStart", "thickEnd"]
    unsigned_int_cols = [
        "chromStart",
        "chromEnd",
        "valid_cov",
        "n_mod",
        "n_canonical",
        "n_othermod",
        "n_delete",
        "n_fail",
        "n_diff",
        "n_nocall",
    ]
    float_cols = ["score_bed"]

    try:
        # Enable StringCache for consistent categorical encoding
        pl.enable_string_cache()

        # Cache or build PyRanges filter with improved caching
        cache_path = os.path.join(
            sample_output_dir, f"{os.path.basename(filter_bed_file)}.pgr_cache"
        )
        if os.path.exists(cache_path):
            with open(cache_path, "rb") as f:
                filter_ranges = pickle.load(f)
        else:
            comp = "gzip" if filter_bed_file.endswith(".gz") else None
            bed_df = pd.read_csv(
                filter_bed_file,
                sep="\t",
                header=None,
                names=["Chromosome", "Start", "End", "cg_label"],
                compression=comp,
                dtype={"Chromosome": str},
            )
            filter_ranges = pr.PyRanges(bed_df[["Chromosome", "Start", "End"]])
            with open(cache_path, "wb") as f:
                pickle.dump(filter_ranges, f)

        # Process new files in chunks using Polars lazy evaluation
        new_frames = []
        for bed in new_files:
            try:
                # Read with pandas first to handle regex separator
                df = pd.read_csv(
                    bed,
                    sep="\s+",
                    header=None,
                    names=cols,
                    dtype={c: str for c in categorical_cols},
                )

                # Validate required columns
                missing_cols = set(cols) - set(df.columns)
                if missing_cols:
                    raise ValueError(f"Missing required columns: {missing_cols}")

                # Convert to Polars for efficient processing
                pl_df = pl.from_pandas(df)

                # Convert numeric columns with proper error handling
                for c in int_cols:
                    pl_df = pl_df.with_columns(pl.col(c).cast(pl.Int64, strict=False))
                for c in unsigned_int_cols:
                    pl_df = pl_df.with_columns(pl.col(c).cast(pl.UInt32, strict=False))
                for c in float_cols:
                    pl_df = pl_df.with_columns(pl.col(c).cast(pl.Float32, strict=False))

                # Filter using PyRanges
                # Rename columns using Polars syntax
                pr_df = pl_df.rename(
                    {"chrom": "Chromosome", "chromStart": "Start", "chromEnd": "End"}
                )
                # Convert to pandas for PyRanges
                pr_df_pandas = pr_df.to_pandas()
                gr = pr.PyRanges(pr_df_pandas[["Chromosome", "Start", "End"]])
                inter = gr.intersect(filter_ranges).df
                filt = pd.merge(
                    inter.rename(
                        columns={
                            "Chromosome": "chrom",
                            "Start": "chromStart",
                            "End": "chromEnd",
                        }
                    ),
                    pl_df.to_pandas(),
                    on=["chrom", "chromStart", "chromEnd"],
                    how="inner",
                )
                new_frames.append(filt)
            except Exception as e:
                logging.error(f"Error processing file {bed}: {str(e)}")
                continue

        if not new_frames:
            logging.warning("No valid data to merge after filtering")
            return

        # Combine new data
        new_df = pd.concat(new_frames, ignore_index=True)

        # If no existing file, just save the new data
        if not os.path.exists(existing_file):
            pl_df = pl.from_pandas(new_df)
            pl_df.write_parquet(output_file)
            return

        # Process existing data in chunks
        existing_df = pl.scan_parquet(existing_file)

        # Combine and aggregate using Polars lazy evaluation
        count_cols = [
            "n_mod",
            "n_canonical",
            "n_othermod",
            "n_delete",
            "n_fail",
            "n_diff",
            "n_nocall",
            "valid_cov",
        ]

        # Convert new data to Polars
        pl_new_df = pl.from_pandas(new_df)

        # Convert existing data to regular DataFrame for concatenation
        existing_df = existing_df.collect()

        # Ensure consistent data types and column order
        for c in categorical_cols:
            if c in existing_df.columns and c in pl_new_df.columns:
                existing_df = existing_df.with_columns(pl.col(c).cast(pl.Categorical))
                pl_new_df = pl_new_df.with_columns(pl.col(c).cast(pl.Categorical))

        # Ensure columns are in the same order
        pl_new_df = pl_new_df.select(existing_df.columns)

        # Validate column names match
        if set(existing_df.columns) != set(pl_new_df.columns):
            missing_cols = set(existing_df.columns) - set(pl_new_df.columns)
            extra_cols = set(pl_new_df.columns) - set(existing_df.columns)
            raise ValueError(
                f"Column mismatch: missing {missing_cols}, extra {extra_cols}"
            )

        # Combine existing and new data
        combined = pl.concat([existing_df, pl_new_df])

        # Define aggregation expressions
        exprs = [
            pl.first("mod_code"),
            pl.mean("score_bed").alias("score_bed"),
            pl.first("strand"),
            pl.first("thickStart"),
            pl.first("thickEnd"),
            pl.first("color"),
            *[pl.sum(c).alias(c) for c in count_cols],
        ]

        # Perform groupby and aggregation
        grouped = combined.group_by(["chrom", "chromStart", "chromEnd"]).agg(exprs)

        # Calculate percent_modified with error handling
        grouped = grouped.with_columns(
            [(pl.col("n_mod") / pl.col("valid_cov") * 100).alias("percent_modified")]
        )

        # Save using Polars' efficient parquet writer
        grouped.write_parquet(output_file)

        # If we are running with mnp_flex, send the file to the server
        if mnpflex_config["mnpuser"] and mnpflex_config["mnppass"]:
            logging.info("Prepare data for mnpflex")
            test_df = load_modkit_data(output_file)
            logging.info(f"Loaded modkit data with shape: {test_df.shape}")

            test_df.rename(
                columns={
                    "chromStart": "start_pos",
                    "chromEnd": "end_pos",
                    "mod_code": "mod",
                    "thickStart": "start_pos2",
                    "thickEnd": "end_pos2",
                    "color": "colour",
                },
                inplace=True,
            )
            test_df.rename(
                columns={
                    "n_canonical": "Ncanon",
                    "n_delete": "Ndel",
                    "n_diff": "Ndiff",
                    "n_fail": "Nfail",
                    "n_mod": "Nmod",
                    "n_nocall": "Nnocall",
                    "n_othermod": "Nother",
                    "valid_cov": "Nvalid",
                    "percent_modified": "score",
                },
                inplace=True,
            )

            # Log column names and data types
            logging.info("Column names after renaming:")
            logging.info(test_df.columns.tolist())
            logging.info("Data types:")
            logging.info(test_df.dtypes)

            savepath = os.path.join(sample_output_dir, f"{sample_id}.mnpflex.bed")
            test_df.to_csv(savepath, sep="\t", index=False, header=False)

            logging.info(f"Saved intermediate file to: {savepath}")
            logging.info(f"Intermediate file size: {os.path.getsize(savepath)} bytes")

            logging.info("Sending file to mnpflex")
            mnpFlex = MnpFlexClient(base_url="https://mnp-flex.org", verify_ssl=True)
            mnpFlex.authenticate(
                username=mnpflex_config["mnpuser"],
                password=mnpflex_config["mnppass"],
                client_id="ROBIN",
                client_secret="SECRET",
            )
            cpg_file = os.path.join(
                os.path.dirname(os.path.abspath(resources.__file__)),
                "mnp_flex_sample_clean.bed",
            )

            logging.info(f"Processing with reference file: {cpg_file}")
            mnpFlex.process_streaming(cpg_file, f"{savepath}", f"{savepath}.clean")

            # Print the length of the clean file
            with open(f"{savepath}.clean", "r") as f:
                clean_file_length = sum(1 for _ in f)
            logging.info(f"Length of clean file: {clean_file_length} lines")

            # Log first few lines of clean file
            with open(f"{savepath}.clean", "r") as f:
                first_lines = [next(f) for _ in range(5)]
                logging.info("First 5 lines of clean file:")
                for line in first_lines:
                    logging.info(line.strip())

            # Upload a sample file
            try:
                response = mnpFlex.upload_sample(
                    file_path=f"{savepath}.clean",
                    sample_name=f"{sample_id}.mnpFlex",
                    disclaimer_confirmed=True,
                )

                logging.info(f"MNP-FLEX upload successful: {response}")
                sample_id = response["id"]
                sample = mnpFlex.get_sample(sample_id)
                logging.info(f"Sample details: {sample}")
                result_status = sample["bed_file_sample"]["analysis_status"]
                
                # Add timeout for analysis completion
                max_wait_time = 300  # 5 minutes
                start_time = time.time()
                while result_status == "initialized" and (time.time() - start_time) < max_wait_time:
                    time.sleep(1)
                    try:
                        sample = mnpFlex.get_sample(sample_id)
                        result_status = sample["bed_file_sample"]["analysis_status"]
                    except Exception as e:
                        logging.error(f"Error checking analysis status: {str(e)}")
                        break

                if (time.time() - start_time) >= max_wait_time:
                    logging.error("Analysis timed out after 5 minutes")
                    mnpFlex.delete_sample(sample_id)
                elif result_status == "Analysis error":
                    logging.error(f"Analysis error for {sample_id}")
                    mnpFlex.delete_sample(sample_id)
                else:
                    try:
                        report_content = mnpFlex.get_sample_report(sample_id)

                        # If the response is binary content (e.g., a PDF), save it to a file
                        if isinstance(report_content, bytes):
                            report_path = os.path.join(
                                sample_output_dir,
                                f'{sample["sample_name"]}_{clean_file_length}_{datetime.now().strftime("%Y%m%d_%H%M%S")}.pdf',
                            )
                            with open(report_path, "wb") as file:
                                file.write(report_content)
                            logging.info(f"Report saved as {report_path}")

                            try:
                                # Extract and store data from the report
                                from robin.reporting.pdf_extractor import PDFExtractor

                                # Initialize PDF extractor with path in sample directory
                                extractor = PDFExtractor(
                                    os.path.join(sample_output_dir, "extracted_data")
                                )

                                # Extract data from the report
                                data = extractor.extract_from_pdf(report_path)
                                if data:
                                    # Save extracted data
                                    extractor.save_data(data)
                                    logging.info(
                                        f"Extracted and stored data from report: {report_path}"
                                    )
                            except Exception as e:
                                logging.error(f"Error extracting PDF data: {str(e)}")
                        else:
                            logging.info(f"Report content: {report_content}")

                    except Exception as e:
                        logging.error(f"Error getting/saving report: {str(e)}")
                    finally:
                        try:
                            mnpFlex.delete_sample(sample_id)
                        except Exception as e:
                            logging.error(f"Error deleting sample: {str(e)}")

            except Exception as e:
                # Log the error but don't let it crash the process
                logging.error("MNP-FLEX upload/analysis failed", exc_info=True)
                logging.info("Continuing with processing despite MNP-FLEX error")
                # Don't re-raise the exception - allow processing to continue

        logging.debug(
            f"✅ Merged with optimized Polars and cache saved to: {output_file}"
        )

    except Exception as e:
        logging.error(f"Error in merge_modkit_files: {str(e)}")
        raise
    finally:
        # Cleanup temporary files if needed
        if os.path.exists(cache_path) and not os.path.exists(filter_bed_file):
            try:
                os.remove(cache_path)
            except Exception as e:
                logging.error(f"Error removing cache file: {str(e)}")
        # Disable StringCache
        pl.disable_string_cache()


def run_modkit(sortfile: str, temp: str, threads: int) -> None:
    """
    Executes modkit on a bam file and extracts the methylation data.

    Args:
        sortfile (str): Path to the sorted BAM file.
        temp (str): Path to the temporary output file.
        threads (int): Number of threads to use.
    """
    cmd = [
        "modkit",
        "pileup",
        "-t",
        str(threads),  # Ensure threads is a string
        "--filter-threshold",
        "0.73",
        "--chunk-size",
        str(threads),
        "--interval-size",
        "25000000",
        "--combine-mods",
        sortfile,
        temp,
        # "--suppress-progress"
    ]
    
    # Run the command
    subprocess.run(cmd, capture_output=True, text=True)
    
def run_samtools_sort(
    file: str, tomerge: List[str], sortfile: str, threads: int
) -> None:
    """
    Sorts BAM files using Samtools.

    Args:
        file (str): Path to the output BAM file.
        tomerge (List[str]): List of BAM files to merge.
        sortfile (str): Path to the sorted BAM file.
        threads (int): Number of threads to use.
    """
    pysam.cat("-o", file, *tomerge)
    pysam.sort("-@", f"{threads}", "--write-index", "-o", sortfile, file)


def check_bam(bamfile):
    """
    Check a BAM file and return its attributes.

    :param bamfile: Path to the BAM file.
    :return: Tuple containing baminfo and bamdata.
    """
    logging.info(f"Checking BAM file: {bamfile}")

    # Check if file exists and is accessible
    if not os.path.exists(bamfile):
        logging.error(f"BAM file does not exist: {bamfile}")
        raise FileNotFoundError(f"BAM file not found: {bamfile}")

    # Check if file is being written to
    try:
        file_size = os.path.getsize(bamfile)
        if os.path.getsize(bamfile) != file_size:
            logging.warning(f"BAM file is still being written: {bamfile}")
            raise IOError("BAM file is still being written")
    except OSError as e:
        logging.error(f"Error checking BAM file size: {str(e)}")
        raise

    # Try to index the BAM file with retries
    max_retries = 3
    retry_delay = 1  # seconds

    for attempt in range(max_retries):
        try:
            # Check if index exists and is newer than BAM file
            bai_file = f"{bamfile}.bai"
            if os.path.exists(bai_file) and os.path.getmtime(
                bai_file
            ) > os.path.getmtime(bamfile):
                logging.info(f"Using existing index for {bamfile}")
            else:
                logging.info(
                    f"Indexing BAM file (attempt {attempt + 1}/{max_retries}): {bamfile}"
                )
                pysam.index(bamfile)
            break
        except pysam.utils.SamtoolsError as e:
            if "Resource temporarily unavailable" in str(e):
                if attempt < max_retries - 1:
                    logging.warning(
                        f"BAM indexing failed (attempt {attempt + 1}), retrying in {retry_delay}s: {str(e)}"
                    )
                    time.sleep(retry_delay)
                    retry_delay *= 2  # Exponential backoff
                    continue
                else:
                    logging.error(
                        f"Failed to index BAM file after {max_retries} attempts: {bamfile}"
                    )
                    raise
            else:
                logging.error(f"Error indexing BAM file: {str(e)}")
                raise

    # Read BAM file
    bam = ReadBam(bamfile)
    baminfo = bam.process_reads()
    bamdata = bam.summary()
    logging.info(f"BAM file processed successfully: {bamfile}")
    return baminfo, bamdata


def sort_bams(files_and_timestamps, watchfolder, file_endings, simtime):
    """
    Sort BAM files by timestamp.

    :param files_and_timestamps: List to store sorted files and timestamps
    :param watchfolder: Folder to watch for BAM files
    :param file_endings: Set of file endings to look for
    :param simtime: Whether to simulate time
    :return: Sorted list of files and timestamps
    """
    import bisect

    def insert_sorted(file_timestamp_tuple):
        file, timestamp, elapsed = file_timestamp_tuple
        datetime_obj = datetime.fromisoformat(timestamp)
        bisect.insort(files_and_timestamps, (datetime_obj, file, elapsed))

    for path, dirs, files in os.walk(watchfolder):
        with alive_bar(len(files)) as bar:
            for f in files:
                if "".join(Path(f).suffix) in file_endings:
                    logging.info(f"Reading BAM file: {os.path.join(path, f)}")
                    if simtime:
                        bam = ReadBam(os.path.join(path, f))
                        baminfo = bam.process_reads()
                        insert_sorted(
                            (
                                os.path.join(path, f),
                                baminfo["last_start"],
                                baminfo["elapsed_time"],
                            )
                        )
                    else:
                        filetime = datetime.fromtimestamp(
                            os.path.getctime(os.path.join(path, f))
                        )
                        elapsedtime = datetime.now() - filetime
                        insert_sorted(
                            (os.path.join(path, f), filetime.isoformat(), elapsedtime)
                        )
                bar()
    return files_and_timestamps
