"""
cli.py

Command-line interface for the ROBIN report generation tool.
"""

import sys
import click
import logging
from pathlib import Path
from typing import Optional

from robin.reporting.report import create_pdf


@click.command()
@click.argument("filename", type=str)
@click.argument("output", type=str)
@click.option(
    "--center",
    required=True,
    help="Center ID running the analysis (e.g., 'Oxford', 'Cambridge', 'London')",
)
@click.option(
    "--export-csv-dir",
    type=str,
    default=None,
    help="Directory to write CSV exports of report data",
)
@click.option(
    "--zip/--no-zip", default=False, help="Zip CSVs and manifest into a single archive"
)
@click.option("--debug", is_flag=True, help="Enable debug logging")
@click.option(
    "--cnv-normalized-difference",
    is_flag=True,
    help=(
        "Plot CNV summary and per-chromosome figures using log2(ploidy / expected copy number) "
        "(overrides the admin Plotting preference)"
    ),
)
def main(
    filename: str,
    output: str,
    center: str,
    export_csv_dir: Optional[str],
    zip: bool,
    debug: bool,
    cnv_normalized_difference: bool,
):
    """Create a PDF report from ROBIN analysis results.

    Args:
        FILENAME: The filename for the PDF report
        OUTPUT: The directory containing the analysis output files
    """
    # Configure logging
    log_level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        level=log_level, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    logger = logging.getLogger(__name__)

    try:
        # Validate input
        output_path = Path(output)
        if not output_path.exists():
            logger.error(f"Output directory does not exist: {output}")
            sys.exit(1)
        if not output_path.is_dir():
            logger.error(f"Output path is not a directory: {output}")
            sys.exit(1)

        # Create PDF
        logger.info(f"Creating PDF report {filename} from {output}")
        pdf_file = create_pdf(
            filename,
            output,
            center=center,
            export_csv_dir=export_csv_dir,
            export_xlsx=False,
            export_zip=zip,
            cnv_summary_normalized=True if cnv_normalized_difference else None,
        )
        logger.info(f"Successfully created PDF report: {pdf_file}")

    except Exception as e:
        logger.error(f"Error creating PDF report: {str(e)}", exc_info=debug)
        sys.exit(1)


if __name__ == "__main__":
    main()
