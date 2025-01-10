"""
Report Sections Package

This package contains the individual sections that make up the ROBIN report.
"""

from .base import ReportSection
from .classification import ClassificationSection
from .cnv import CNVSection
from .fusion import FusionSection
from .coverage import CoverageSection
from .mgmt import MGMTSection
from .run_data import RunDataSection
from .disclaimer import DisclaimerSection
from .variants import VariantsSection

__all__ = [
    "ReportSection",
    "ClassificationSection",
    "CNVSection",
    "FusionSection",
    "CoverageSection",
    "MGMTSection",
    "RunDataSection",
    "DisclaimerSection",
    "VariantsSection"
]
