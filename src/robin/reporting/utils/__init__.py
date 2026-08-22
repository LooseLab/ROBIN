"""
Report Utilities Package

This package contains utility functions used across the ROBIN report generation code.
"""

from .formatting import (
    convert_to_space_separated_string,
    format_number,
    format_timestamp,
    split_text,
)

__all__ = [
    "format_number",
    "format_timestamp",
    "convert_to_space_separated_string",
    "split_text",
]
