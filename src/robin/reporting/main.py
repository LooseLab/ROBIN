"""
main.py

This module is the entry point for generating the PDF report.
"""

import argparse

from robin.reporting.report import create_pdf

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a PDF report from a ROBIN analysis output directory."
    )
    parser.add_argument(
        "output_dir",
        help="Directory containing ROBIN analysis outputs for one sample/run.",
    )
    parser.add_argument(
        "--center",
        default="dev",
        help="Centre ID recorded in the report (default: dev).",
    )
    parser.add_argument(
        "--pdf",
        default="sample_report.pdf",
        help="Output PDF filename (default: sample_report.pdf).",
    )
    args = parser.parse_args()
    create_pdf(args.pdf, args.output_dir, args.center)
