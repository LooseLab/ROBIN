"""
ROBIN Report Generation Package

This package contains all the code needed to generate PDF reports from ROBIN analysis results.
"""

from .report import create_pdf, RobinReport

__all__ = ["create_pdf", "RobinReport"]
