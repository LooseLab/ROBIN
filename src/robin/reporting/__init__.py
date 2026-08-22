"""
ROBIN Report Generation Package

This package contains all the code needed to generate PDF reports from ROBIN analysis results.
"""

from .report import RobinReport, create_pdf

__all__ = ["create_pdf", "RobinReport"]
