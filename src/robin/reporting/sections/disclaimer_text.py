"""
disclaimer_text.py

This module contains the centralized disclaimer text used across the ROBIN application.
"""

DISCLAIMER_TEXT = (
    "This report and the data contained within it are intended for research use only and "
    "should not be used for direct diagnostic purposes. The methylation-based classification "
    "and other analyses provided here may be considered by neuropathologists as supplementary "
    "information in the context of comprehensive diagnostic assessment, which should include "
    "clinical history, radiological findings, and complete histopathological and molecular evaluation. "
    "The final interpretation and diagnosis should always be made by qualified healthcare professionals "
    "based on all available information."
)

EXTENDED_DISCLAIMER_TEXT = (
    DISCLAIMER_TEXT + "\n\nThe authors of this software note that:\n"
    "Oxford Nanopore Technologies products are not intended for use for health assessment "
    "or to diagnose, treat, mitigate, cure, or prevent any disease or condition."
)
