"""
formatting.py

This module contains utility functions for formatting data in the report.
"""

from datetime import datetime
import logging

logger = logging.getLogger(__name__)


def format_number(n):
    """Format number with commas for readability.

    Args:
        n: Number to format

    Returns:
        Formatted string with commas
    """
    try:
        return "{:,}".format(int(float(n)))
    except (ValueError, TypeError):
        return str(n)


def format_timestamp(timestamp_str):
    """Convert timestamp string to readable format.

    Args:
        timestamp_str: Timestamp string to format

    Returns:
        Formatted datetime string
    """
    clean_ts = str(timestamp_str).strip("[]'\"")
    try:
        dt = datetime.strptime(clean_ts.split("+")[0], "%Y-%m-%dT%H:%M:%S.%f")
        return dt.strftime("%B %d, %Y %H:%M:%S")
    except Exception as e:
        logger.debug(f"Error formatting timestamp {timestamp_str}: {e}")
        return str(timestamp_str)


def convert_to_space_separated_string(value):
    """Convert various data types to space-separated string.

    Args:
        value: Value to convert (can be list, tuple, or single value)

    Returns:
        Space-separated string
    """
    if isinstance(value, (list, tuple)):
        return " ".join(str(x) for x in value)
    return str(value)


def split_text(text, max_length=80):
    """Split text into lines of maximum length.

    Args:
        text: Text to split
        max_length: Maximum line length

    Returns:
        List of lines
    """
    words = text.split()
    lines = []
    current_line = []
    current_length = 0

    for word in words:
        word_length = len(word)
        if current_length + word_length + 1 <= max_length:
            current_line.append(word)
            current_length += word_length + 1
        else:
            lines.append(" ".join(current_line))
            current_line = [word]
            current_length = word_length + 1

    if current_line:
        lines.append(" ".join(current_line))

    return lines
