import gzip
import shutil
import os
import logging

# Use the main logger configured in the main application
logger = logging.getLogger(__name__)


def decompress_gzip_file(file_path):
    # Ensure the file has a .gz extension
    if not file_path.endswith(".gz"):
        raise ValueError("The file is not a .gz file")

    decompressed_file_path = file_path[:-3]  # Remove .gz extension

    # Check if the decompressed file already exists
    if not os.path.exists(decompressed_file_path):
        print("First run. Preparing code.")
        logger.info(f"Decompressing {file_path}...")

        with gzip.open(file_path, "rb") as f_in:
            with open(decompressed_file_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        logger.info(f"Decompressed {file_path} to {decompressed_file_path}")
    else:
        logger.info(f"{decompressed_file_path} already exists, skipping decompression.")

    logger.info("Decompression check complete for file:", file_path)


# Example usage:
# decompress_gzip_file('/path/to/your/file1.gz')
# decompress_gzip_file('/path/to/your/file2.gz')
