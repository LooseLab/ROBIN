"""
main.py

This module is the entry point for generating the PDF report.
"""

import os
from robin.reporting.report import create_pdf

if __name__ == "__main__":
    output_dir = (
        "/Users/mattloose/GIT/niceGUI/cnsmeth/fusion_output/ds1305_Intraop003_a_NEB"
    )
    fonts_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fonts")
    images_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "images")

    create_pdf("sample_report.pdf", output_dir)
