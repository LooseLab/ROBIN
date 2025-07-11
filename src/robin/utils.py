from typing import List, Optional
import gc
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
from tqdm import tqdm
from robin.utilities.ReadBam import ReadBam
from robin.utilities.mnp_flex import APIClient as MnpFlexClient
from robin import resources
import json


def merge_modkit_files(
    new_files: List[str],
    existing_file: str,
    output_file: str,
    filter_bed_file: str,
    sample_id: str,
    output_dir: str,
    mnpflex_config: Optional[dict],
    num_bam_files_seen: int,
) -> None:
    """
    Merge modkit files with optimized column set and improved caching.

    This function uses only essential columns to reduce memory usage and processing time:
    - chrom, chromStart: Required for all classifiers
    - mod_code: Required for Sturgeon filtering
    - strand: Required for proper aggregation
    - valid_cov: Required for coverage calculation
    - percent_modified: Primary methylation data (required)
    - n_mod, n_canonical: Required for modification counts

    Args:
        new_files (List[str]): List of new modkit files to merge
        existing_file (str): Path to existing parquet file
        output_file (str): Path to output parquet file
        filter_bed_file (str): Path to BED file for filtering
        sample_id (str): Sample ID for organizing output files
        output_dir (str): Base output directory
        mnpflex_config (Optional[dict]): Configuration for MNP-FLEX integration
        num_bam_files_seen (int): Number of BAM files being processed
    """
    # Create sample-specific output directory
    sample_output_dir = os.path.join(output_dir, sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)

    # Define optimized schema with only essential columns
    essential_cols = [
        "chrom",
        "chromStart",
        "mod_code",
        "strand",
        "valid_cov",
        "percent_modified",
        "n_mod",
        "n_canonical",
    ]

    categorical_cols = ["chrom", "mod_code", "strand"]
    unsigned_int_cols = ["chromStart", "valid_cov", "n_mod", "n_canonical"]
    float_cols = ["percent_modified"]

    try:
        # Enable StringCache for consistent categorical encoding
        pl.enable_string_cache()

        # Track cumulative BAM file count
        cumulative_bam_file_count = 0

        # Check if existing file has metadata about BAM file count
        if os.path.exists(existing_file):
            try:
                # Read existing metadata
                metadata_file = existing_file.replace(".parquet", "_metadata.json")
                if os.path.exists(metadata_file):
                    with open(metadata_file, "r") as f:
                        metadata = json.load(f)
                        cumulative_bam_file_count = metadata.get("bam_file_count", 0)
                        logging.info(
                            f"Found existing metadata with {cumulative_bam_file_count} BAM files"
                        )
            except Exception as e:
                logging.warning(f"Could not read existing metadata: {str(e)}")
                cumulative_bam_file_count = 0

        # Add the number of new BAM files being processed
        cumulative_bam_file_count += num_bam_files_seen

        logging.info(
            f"Total cumulative BAM files contributing to parquet: {cumulative_bam_file_count} (added {num_bam_files_seen} new files)"
        )

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
                # Read all 18 columns from the input file
                full_cols = [
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

                df = pd.read_csv(
                    bed,
                    sep="\s+",
                    header=None,
                    names=full_cols,
                    dtype={c: str for c in ["chrom", "mod_code", "strand", "color"]},
                )

                # Validate required columns
                missing_cols = set(full_cols) - set(df.columns)
                if missing_cols:
                    raise ValueError(f"Missing required columns: {missing_cols}")

                # Extract only essential columns for processing
                df_essential = df[essential_cols].copy()

                # Convert to Polars for efficient processing
                pl_df = pl.from_pandas(df_essential)

                # Convert numeric columns with proper error handling
                for c in unsigned_int_cols:
                    pl_df = pl_df.with_columns(pl.col(c).cast(pl.UInt32, strict=False))
                for c in float_cols:
                    pl_df = pl_df.with_columns(pl.col(c).cast(pl.Float32, strict=False))

                # Filter using PyRanges
                # Rename columns using Polars syntax
                pr_df = pl_df.rename({"chrom": "Chromosome", "chromStart": "Start"})
                # Add chromEnd for filtering (calculate from chromStart)
                pr_df = pr_df.with_columns((pl.col("Start") + 1).alias("End"))

                # Convert to pandas for PyRanges
                pr_df_pandas = pr_df.to_pandas()
                gr = pr.PyRanges(pr_df_pandas[["Chromosome", "Start", "End"]])
                inter = gr.intersect(filter_ranges).df
                filt = pd.merge(
                    inter.rename(
                        columns={
                            "Chromosome": "chrom",
                            "Start": "chromStart",
                        }
                    ),
                    pl_df.to_pandas(),
                    on=["chrom", "chromStart"],
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
        if len(new_frames) > 1:
            new_df = pd.concat(new_frames, ignore_index=True)
        else:
            new_df = new_frames[0]

        # If no existing file, just save the new data
        if not os.path.exists(existing_file):
            pl_df = pl.from_pandas(new_df)
            pl_df.write_parquet(output_file)

            # Save metadata with cumulative BAM file count
            metadata = {
                "bam_file_count": cumulative_bam_file_count,
                "last_updated": datetime.now().isoformat(),
                "sample_id": sample_id,
                "files_added_in_this_update": num_bam_files_seen,
                "column_format": "optimized",  # Mark as optimized format
            }
            metadata_file = output_file.replace(".parquet", "_metadata.json")
            with open(metadata_file, "w") as f:
                json.dump(metadata, f, indent=2)

            logging.info(
                f"Created new optimized parquet file with {cumulative_bam_file_count} cumulative BAM files"
            )
            return

        # Process existing data in chunks
        existing_df = pl.scan_parquet(existing_file)

        # Convert new data to Polars
        pl_new_df = pl.from_pandas(new_df)

        # Convert existing data to regular DataFrame for concatenation
        existing_df = existing_df.collect()

        # Check if existing file is in old format (18 columns) or new format (8 columns)
        is_old_format = (
            len(existing_df.columns) > 10
        )  # More than 10 columns indicates old format

        if is_old_format:
            # Convert old format to new format by selecting only essential columns
            logging.info("Converting existing file from old format to optimized format")
            existing_df = existing_df.select(essential_cols)

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

        # Define aggregation expressions for essential columns only
        exprs = [
            pl.first("mod_code"),
            pl.first("strand"),
            pl.mean("percent_modified").alias("percent_modified"),
            *[pl.sum(c).alias(c) for c in ["valid_cov", "n_mod", "n_canonical"]],
        ]

        # Perform groupby and aggregation
        grouped = combined.group_by(["chrom", "chromStart"]).agg(exprs)

        # Save using Polars' efficient parquet writer
        grouped.write_parquet(output_file)

        # Save metadata with updated cumulative BAM file count
        metadata = {
            "bam_file_count": cumulative_bam_file_count,
            "last_updated": datetime.now().isoformat(),
            "sample_id": sample_id,
            "files_added_in_this_update": num_bam_files_seen,
            "column_format": "optimized",  # Mark as optimized format
        }
        metadata_file = output_file.replace(".parquet", "_metadata.json")
        with open(metadata_file, "w") as f:
            json.dump(metadata, f, indent=2)

        logging.info(
            f"Updated optimized parquet file with {cumulative_bam_file_count} cumulative BAM files (added {num_bam_files_seen} in this update)"
        )

        # If we are running with mnp_flex, send the file to the server
        if mnpflex_config["mnpuser"] and mnpflex_config["mnppass"]:
            logging.info("Prepare data for mnpflex")
            # Load the optimized data and reconstruct full format for MNP-FLEX compatibility
            test_df = reconstruct_full_bedmethyl_for_mnpflex(output_file)
            logging.info(
                f"Reconstructed full bedmethyl data with shape: {test_df.shape}"
            )

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
                while (
                    result_status == "initialized"
                    and (time.time() - start_time) < max_wait_time
                ):
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

            except Exception:
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
        gc.collect()


def reconstruct_full_bedmethyl_for_mnpflex(parquet_file_path: str) -> pd.DataFrame:
    """
    Reconstruct full 18-column bedmethyl format from optimized parquet file for MNP-FLEX compatibility.

    This function takes the optimized 8-column parquet file and reconstructs the full
    18-column format expected by MNP-FLEX integration.

    Args:
        parquet_file_path (str): Path to the optimized parquet file

    Returns:
        pd.DataFrame: Full 18-column bedmethyl data compatible with MNP-FLEX
    """
    try:
        # Read the optimized parquet file
        df = pd.read_parquet(parquet_file_path)

        # Check if this is already in full format
        if len(df.columns) >= 15:  # Full format has 18 columns
            logging.info("Parquet file already in full format, returning as-is")
            return df

        # Reconstruct missing columns
        logging.info("Reconstructing full bedmethyl format from optimized data")

        # Add calculated columns
        df["chromEnd"] = df["chromStart"] + 1
        df["score_bed"] = df["percent_modified"]
        df["thickStart"] = df["chromStart"]
        df["thickEnd"] = df["chromEnd"]
        df["color"] = "255,0,0"

        # Add columns with default values
        df["n_othermod"] = 0
        df["n_delete"] = 0
        df["n_fail"] = 0
        df["n_diff"] = 0
        df["n_nocall"] = 0

        # Ensure correct column order
        full_columns = [
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

        df = df[full_columns]

        logging.info(
            f"Reconstructed full bedmethyl data with {len(df)} rows and {len(df.columns)} columns"
        )
        return df

    except Exception as e:
        logging.error(f"Error reconstructing full bedmethyl data: {str(e)}")
        raise


# Pre-compile default count dictionary to avoid repeated creation
DEFAULT_COUNT_DICT = {
    "mod": 0,
    "canonical": 0,
    "other_mod": 0,
    "delete": 0,
    "fail": 0,
    "diff": 0,
    "nocall": 0,
    "total": 0,
    "max_prob_C": 0,
    "max_prob_m": 0,
    "max_prob_h": 0,
}

# Pre-define constants to avoid repeated string creation and calculations
PLUS_STRAND = "+"
MINUS_STRAND = "-"
COLOR_VALUE = "255,0,0"
CANONICAL_CODE = "C"


def run_matkit(sortfile: str, temp: str) -> None:
    """
    Executes modkit2-style processing on a bam file and extracts the methylation data.

    Args:
        sortfile (str): Path to the sorted BAM file.
        temp (str): Path to the temporary output file.
        threads (int): Number of threads to use (not used in current implementation).
    """
    logging.info(f"Processing BAM file with modkit2: {sortfile}")

    # Use the improved approach from modkit2.py with memory efficiency and better classification
    counts, debug_data = process_bam_counts_improved(
        sortfile,
        threshold=0.73,
        combine_mods=True,
        chrom_filter=None,
        ignore_supp=False,
        ref_fasta=None,
    )

    # Write output to BEDMethyl format
    write_bedmethyl_improved(temp, counts, debug_probs=False)
    # Clear memory
    del counts
    del debug_data
    logging.info(f"Modkit2 processing complete. Output written to: {temp}")
    gc.collect()


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
    subprocess.run(cmd)  # , capture_output=True, text=True)


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
    try:
        pysam.cat("-o", file, *tomerge)
    except Exception as e:
        print(f"Error merging BAM files: {e}")
        raise
    try:
        pysam.sort("-@", f"{threads}", "--write-index", "-o", sortfile, file)
    except Exception as e:
        print(f"Error merging BAM files: {e}")
        raise


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


def map_read_to_ref_positions(read):
    """
    Map read positions to reference positions for aligned bases.

    This function creates a mapping between read positions and reference positions
    for bases that are aligned (matches only). This is essential for correctly
    positioning modification calls in the reference coordinate system.

    According to modkit documentation, modification calls must be positioned
    relative to the reference genome, not the read sequence.

    Args:
        read: pysam.AlignedSegment object containing the read alignment

    Returns:
        dict: Mapping from read position to reference position for aligned bases

    References:
        - Modkit pileup documentation: https://nanoporetech.github.io/modkit/intro_pileup.html
        - BAM format specification for modified bases
    """
    # Optimized version: use direct iteration and avoid try/except overhead
    ref_map = {}
    try:
        # Try with sequence first (requires MD tag for base information)
        # This provides the most accurate mapping when MD tag is available
        for read_pos, ref_pos, base in read.get_aligned_pairs(
            with_seq=True, matches_only=True
        ):
            if read_pos is not None and ref_pos is not None:
                ref_map[read_pos] = ref_pos
    except ValueError:
        # Fall back to without sequence if MD tag is missing
        # This is less precise but still functional for modification calling
        for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
            if read_pos is not None and ref_pos is not None:
                ref_map[read_pos] = ref_pos
    return ref_map


def write_bedmethyl_improved(output_path, counts, debug_probs=False):
    """
    Write counts to BEDMethyl format with improved performance.

    This function generates the 18-column BEDMethyl format as specified in the
    modkit documentation with optimizations from modkit2.py. The format includes
    all the classification counts and follows the exact column order and naming conventions.

    BEDMethyl format columns (per modkit documentation):
    1. chrom - Reference sequence name
    2. start - 0-based start position
    3. end - 0-based exclusive end position
    4. modified base code - Single letter code for modified base
    5. score - Equal to Nvalid_cov
    6. strand - '+' for positive, '-' for negative, '.' when combined
    7. start position - Included for compatibility
    8. end position - Included for compatibility
    9. color - Always 255,0,0
    10. Nvalid_cov - Valid coverage
    11. fraction modified - Nmod / Nvalid_cov
    12. Nmod - Number of modified calls
    13. Ncanonical - Number of canonical calls
    14. Nother_mod - Number of other modification calls
    15. Ndelete - Number of deletions
    16. Nfail - Number of failed calls
    17. Ndiff - Number of different base calls
    18. Nnocall - Number of no-calls

    Args:
        output_path: Path to output BEDMethyl file
        counts: Dictionary containing classification counts with tuple keys
        debug_probs: Whether to include probability columns for debugging

    References:
        - BEDMethyl format specification: https://github.com/nanoporetech/modkit?tab=readme-ov-file#description-of-bedmethyl-output
        - Modkit pileup output documentation
    """
    with open(output_path, "w") as out:
        # Sort by chromosome, position, strand, and modification code for consistent output
        sorted_sites = sorted(counts.keys(), key=lambda x: (x[0], x[1], x[2], x[3]))

        # Pre-allocate string templates for better performance
        pos_str = str
        percent_str = "{:.2f}".format

        for site_key in sorted_sites:
            chrom, pos, strand, mod_code = site_key
            c = counts[site_key]
            total = c["total"]
            # Calculate fraction modified as Nmod / Nvalid_cov
            percent = (c["mod"] / total) * 100 if total else 0.0

            # Construct the 18-column BEDMethyl line using optimized string operations
            pos_plus_1 = pos + 1
            fields = [
                chrom,
                pos_str(pos),
                pos_str(pos_plus_1),  # Columns 1-3: chrom, start, end
                mod_code,
                pos_str(total),
                strand,  # Columns 4-6: mod_code, score, strand
                pos_str(pos),
                pos_str(pos_plus_1),
                COLOR_VALUE,  # Columns 7-9: start, end, color
                pos_str(total),
                percent_str(percent),  # Columns 10-11: Nvalid_cov, fraction_modified
                pos_str(c["mod"]),
                pos_str(c["canonical"]),  # Columns 12-13: Nmod, Ncanonical
                pos_str(c["other_mod"]),
                pos_str(c["delete"]),  # Columns 14-15: Nother_mod, Ndelete
                pos_str(c["fail"]),
                pos_str(c["diff"]),
                pos_str(c["nocall"]),  # Columns 16-18: Nfail, Ndiff, Nnocall
            ]

            # Add probability columns only if debug_probs is True
            if debug_probs:
                # Compute probabilities for canonical, 5mC, and 5hmC for each site
                # Output as 0-1 float (divide by 255)
                prob_canonical = c.get("max_prob_C", 0.0) / 255.0
                prob_5mc = c.get("max_prob_m", 0.0) / 255.0
                prob_5hmc = c.get("max_prob_h", 0.0) / 255.0
                fields.extend(
                    [f"{prob_canonical:.3f}", f"{prob_5mc:.3f}", f"{prob_5hmc:.3f}"]
                )

            out.write("\t".join(fields) + "\n")


def process_bam_counts_improved(
    bam_path,
    threshold=0.7,
    combine_mods=False,
    chrom_filter=None,
    ignore_supp=False,
    ref_fasta=None,
    debug_positions=None,
):
    """
    Process BAM file and return detailed counts for BEDMethyl format with improved memory efficiency.

    This function implements the detailed counting logic required for the full
    18-column BEDMethyl format with optimizations from modkit2.py. It tracks all
    the different classification categories (modified, canonical, other_mod, delete, fail, diff, nocall)
    for each genomic position using memory-efficient data structures.

    According to modkit documentation, the classification logic is:
    - Nmod: Calls passing filters that were classified as modified (prob >= threshold)
    - Ncanonical: Calls passing filters that were classified as canonical (prob == 0)
    - Nfail: Calls that were below the threshold but above zero (0 < prob < threshold)
    - Nother_mod: Calls classified as modified but with different modification type
    - Ndelete: Number of deletions at this position
    - Ndiff: Number of different base calls
    - Nnocall: Number of no-calls

    Only output sites with at least one MM/ML tag (i.e., modkit output sites).
    This matches modkit's behavior of only outputting sites with modification data.

    Args:
        bam_path: Path to BAM file containing modification data
        threshold: Probability threshold (0-1) for modification calling
        combine_mods: Whether to combine different modification types
        chrom_filter: Optional chromosome filter to limit processing
        ignore_supp: Whether to ignore supplementary alignments
        ref_fasta: Reference genome file for validation (optional)
        debug_positions: Set of (chrom, pos, strand) tuples to debug

    Returns:
        tuple: (counts_dict, debug_data) where counts_dict contains classification counts with tuple keys

    References:
        - BEDMethyl column descriptions: https://github.com/nanoporetech/modkit?tab=readme-ov-file#bedmethyl-column-descriptions
        - Modkit classification logic documentation
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    thresh = int(threshold * 255)  # Convert to 8-bit integer threshold
    iterable = (
        bam.fetch(contig=chrom_filter) if chrom_filter else bam.fetch(until_eof=True)
    )

    # Use a flat dictionary with tuple keys for memory efficiency
    # Structure: counts[(chrom, pos, strand, mod_code)] = count_dict
    counts = {}

    # Track all sites with at least one MM/ML tag (modkit output sites)
    # This ensures we only output sites that modkit would output
    mod_sites = set()

    # Load reference genome if provided for validation
    ref_fasta_obj = None
    if ref_fasta:
        ref_fasta_obj = pysam.FastaFile(ref_fasta)

    # Debug tracking for specific positions
    debug_data = {} if debug_positions else None

    for read in tqdm(iterable, desc="Processing reads", disable=True):
        # Skip unmapped, secondary, and optionally supplementary alignments
        # Note: modkit appears to filter out secondary alignments by default
        # This matches the behavior of the original modkit tool
        if (
            read.is_unmapped
            or read.is_secondary
            or (ignore_supp and read.is_supplementary)
        ):
            continue

        chrom = bam.get_reference_name(read.reference_id)
        ref_map = map_read_to_ref_positions(read)

        try:
            # Extract modification data from MM/ML tags
            mods = read.modified_bases or {}
        except AttributeError:
            mods = {}

        if combine_mods:
            # Combine all modification types into a single 'mod' category
            # This matches modkit's --combine-mods behavior
            # Use a more memory-efficient approach with direct processing

            # Collect all modification data for this read
            read_sites = {}  # (refpos, strand) -> {mod_code: [probs]}

            for key, values in mods.items():
                if not isinstance(key, tuple) or len(key) != 3:
                    continue
                _, strand_flag, mod_code = key
                strand = PLUS_STRAND if strand_flag == 0 else MINUS_STRAND
                for rpos, prob in values:
                    if rpos not in ref_map:
                        continue
                    refpos = ref_map[rpos]

                    # Store modification probabilities for this site
                    site_key = (refpos, strand)
                    if site_key not in read_sites:
                        read_sites[site_key] = {}
                    if mod_code not in read_sites[site_key]:
                        read_sites[site_key][mod_code] = []
                    read_sites[site_key][mod_code].append(prob)
                    # Track this as a modkit output site
                    mod_sites.add((chrom, refpos, strand, CANONICAL_CODE))

            # Process each site for this read
            for (refpos, strand), mod_probs in read_sites.items():
                # Update counts directly - each read contributes once per site
                site_key = (chrom, refpos, strand, CANONICAL_CODE)
                if site_key not in counts:
                    counts[site_key] = DEFAULT_COUNT_DICT.copy()

                c = counts[site_key]
                c["total"] += 1

                # Track max probabilities for C, m, h from the original mod_probs
                total_mod_prob = 0.0
                max_prob_m = c["max_prob_m"]
                max_prob_h = c["max_prob_h"]

                for mod_code, probs in mod_probs.items():
                    if mod_code != CANONICAL_CODE:  # Skip canonical base
                        max_mod_prob = max(probs)
                        max_mod_prob_01 = max_mod_prob / 255.0  # Convert to 0-1
                        total_mod_prob += max_mod_prob_01
                        if mod_code == "m":
                            if max_mod_prob > max_prob_m:
                                max_prob_m = max_mod_prob
                                c["max_prob_m"] = max_mod_prob
                        elif mod_code == "h":
                            if max_mod_prob > max_prob_h:
                                max_prob_h = max_mod_prob
                                c["max_prob_h"] = max_mod_prob

                # Canonical probability is the complement
                canonical_prob = max(0.0, 1.0 - total_mod_prob)
                canonical_prob_255 = int(canonical_prob * 255)
                if canonical_prob_255 > c["max_prob_C"]:
                    c["max_prob_C"] = canonical_prob_255

                # Use the highest probability among canonical, 5mC, and 5hmC for classification
                max_prob = max(canonical_prob_255, max_prob_m, max_prob_h)

                # Apply modkit classification logic based on the highest probability
                if max_prob >= thresh:
                    if canonical_prob_255 == max_prob:
                        c["canonical"] += 1  # Canonical call
                    else:
                        c["mod"] += 1  # Modified call
                elif max_prob == 0:
                    c["canonical"] += 1  # Canonical call
                else:
                    c["fail"] += 1  # Failed call (0 < prob < threshold)

                # Debug tracking for specific positions
                if debug_positions and (chrom, refpos, strand) in debug_positions:
                    debug_key = (chrom, refpos, strand)
                    if debug_key not in debug_data:
                        debug_data[debug_key] = {
                            "read_id": read.query_name,
                            "probs": {"C": [], "m": [], "h": []},
                            "classification": None,
                        }

                    # Store probabilities for this read
                    for mod_code, probs in mod_probs.items():
                        if mod_code in ["C", "m", "h"]:
                            debug_data[debug_key]["probs"][mod_code].extend(probs)

                    # Store classification
                    if max_prob >= thresh:
                        if canonical_prob_255 == max_prob:
                            debug_data[debug_key]["classification"] = "canonical"
                        else:
                            debug_data[debug_key]["classification"] = "modified"
                    elif max_prob == 0:
                        debug_data[debug_key]["classification"] = "canonical"
                    else:
                        debug_data[debug_key]["classification"] = "fail"
        else:
            # Process each modification type separately
            # This preserves the distinction between different modification types
            read_mod_calls = {}  # (refpos, strand, mod_code) -> [probs]

            for key, values in mods.items():
                if isinstance(key, tuple) and len(key) == 3:
                    _, strand_flag, mod_code = key
                    strand = PLUS_STRAND if strand_flag == 0 else MINUS_STRAND
                    for rpos, prob in values:
                        rp = ref_map.get(rpos)
                        if rp is not None:
                            # Validate reference base if available
                            # This ensures modifications are only called at appropriate reference positions
                            if ref_fasta_obj:
                                try:
                                    ref_base = ref_fasta_obj.fetch(
                                        chrom, rp, rp + 1
                                    ).upper()
                                    # Check if reference base matches expected canonical base
                                    expected_base = (
                                        "C"
                                        if mod_code in ["m", "h"]
                                        else "A" if mod_code == "a" else None
                                    )
                                    if expected_base and ref_base != expected_base:
                                        continue
                                except Exception:
                                    continue

                            site_key = (rp, strand, mod_code)
                            if site_key not in read_mod_calls:
                                read_mod_calls[site_key] = []
                            read_mod_calls[site_key].append(prob)
                            # Track this as a modkit output site
                            mod_sites.add((chrom, rp, strand, mod_code))

            # For each site, update counts using "most likely" classification
            # This implements the modkit logic for handling multiple modification types
            for (rp, strand, mod_code), probs in read_mod_calls.items():
                site_key = (chrom, rp, strand, mod_code)
                if site_key not in counts:
                    counts[site_key] = DEFAULT_COUNT_DICT.copy()

                c = counts[site_key]
                c["total"] += 1
                max_prob = max(probs)

                # Apply modkit classification logic
                if max_prob >= thresh:
                    c["mod"] += 1  # Modified call
                elif max_prob == 0:
                    # Check if there are other modification types at this site with prob >= threshold
                    # Only classify as other_mod if the other modification is above threshold
                    # Optimized: use early exit for better performance
                    has_other_mod_above_thresh = False
                    for (r, s, m), ps in read_mod_calls.items():
                        if (r, s) == (rp, strand) and m != mod_code:
                            for prob in ps:
                                if prob >= thresh:
                                    has_other_mod_above_thresh = True
                                    break
                            if has_other_mod_above_thresh:
                                break

                    if has_other_mod_above_thresh:
                        c["other_mod"] += 1  # Other modification call
                    else:
                        c["canonical"] += 1  # Canonical call
                else:
                    # Check if there are other modification types at this site with prob >= threshold
                    # If so, classify as other_mod instead of fail
                    # Optimized: use early exit for better performance
                    has_other_mod_above_thresh = False
                    for (r, s, m), ps in read_mod_calls.items():
                        if (r, s) == (rp, strand) and m != mod_code:
                            for prob in ps:
                                if prob >= thresh:
                                    has_other_mod_above_thresh = True
                                    break
                            if has_other_mod_above_thresh:
                                break

                    if has_other_mod_above_thresh:
                        c["other_mod"] += 1  # Other modification call
                    else:
                        c["fail"] += 1  # Failed call (0 < prob < threshold)

        # Process reads without modification data as canonical bases
        # Note: modkit only outputs positions with explicit modification data (MM tags)
        # Reads without modification data are not used to generate output positions
        # This matches modkit's behavior of only outputting sites with modification data
        if not mods:
            continue  # Skip reads without modification data

    bam.close()
    if ref_fasta_obj:
        ref_fasta_obj.close()

    # Filter counts to only include sites that match modkit's filtering logic
    # and have at least one passing call (Nmod > 0 or Ncanonical > 0)
    filtered_counts = {}
    for site_key, c in counts.items():
        chrom, pos, strand, mod_code = site_key
        # modkit outputs only positions with at least one passing call
        if site_key in mod_sites and (c["mod"] > 0 or c["canonical"] > 0):
            filtered_counts[site_key] = c

    # Force garbage collection to free memory
    gc.collect()

    return filtered_counts, debug_data


def get_bam_file_count(parquet_file_path: str) -> int:
    """
    Get the cumulative number of BAM files that have contributed to a parquet file.

    Args:
        parquet_file_path (str): Path to the parquet file

    Returns:
        int: Cumulative number of BAM files that have contributed to the parquet file
    """
    try:
        metadata_file = parquet_file_path.replace(".parquet", "_metadata.json")
        if os.path.exists(metadata_file):
            with open(metadata_file, "r") as f:
                metadata = json.load(f)
                return metadata.get("bam_file_count", 0)
        else:
            logging.warning(f"No metadata file found for {parquet_file_path}")
            return 0
    except Exception as e:
        logging.error(f"Error reading BAM file count from metadata: {str(e)}")
        return 0


def get_parquet_metadata(parquet_file_path: str) -> dict:
    """
    Get all metadata associated with a parquet file.

    Args:
        parquet_file_path (str): Path to the parquet file

    Returns:
        dict: Dictionary containing metadata (bam_file_count, last_updated, sample_id, files_added_in_this_update)
    """
    try:
        metadata_file = parquet_file_path.replace(".parquet", "_metadata.json")
        if os.path.exists(metadata_file):
            with open(metadata_file, "r") as f:
                return json.load(f)
        else:
            logging.warning(f"No metadata file found for {parquet_file_path}")
            return {}
    except Exception as e:
        logging.error(f"Error reading metadata: {str(e)}")
        return {}


def get_bam_file_history(parquet_file_path: str) -> dict:
    """
    Get detailed BAM file contribution history for a parquet file.

    Args:
        parquet_file_path (str): Path to the parquet file

    Returns:
        dict: Dictionary containing detailed history information
    """
    metadata = get_parquet_metadata(parquet_file_path)
    if not metadata:
        return {}

    return {
        "total_bam_files": metadata.get("bam_file_count", 0),
        "last_update": metadata.get("last_updated", "Unknown"),
        "sample_id": metadata.get("sample_id", "Unknown"),
        "files_in_last_update": metadata.get("files_added_in_this_update", 0),
        "metadata_file": parquet_file_path.replace(".parquet", "_metadata.json"),
    }
