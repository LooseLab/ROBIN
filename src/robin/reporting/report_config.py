"""
report_config.py

This module contains configuration classes and constants for the report generator.
"""

from enum import Flag, auto

class ReportSection(Flag):
    """Flags for configuring which sections to include in the report."""
    SUMMARY = auto()
    ANALYSIS = auto()
    SEQUENCING = auto()
    DISCLAIMER = auto()
    
    # Predefined configurations
    DEFAULT = SUMMARY | ANALYSIS | SEQUENCING | DISCLAIMER
    SUMMARY_ONLY = SUMMARY | DISCLAIMER
    ANALYSIS_REPORT = SUMMARY | ANALYSIS | DISCLAIMER 