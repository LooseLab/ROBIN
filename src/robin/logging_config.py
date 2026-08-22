#!/usr/bin/env python3
"""
Logging configuration for robin

This module provides a configurable logging system that supports:
- Global log level configuration
- Per-job log level configuration
- Structured logging with job context
- Different log levels for different components
"""

import logging
import os
import sys
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any, Dict, Optional


@dataclass
class LogConfig:
    """Configuration for logging levels"""

    global_level: str = "INFO"
    job_levels: Dict[str, str] = field(default_factory=dict)
    component_levels: Dict[str, str] = field(default_factory=dict)

    def get_level(self, job_type: str = None, component: str = None) -> str:
        """Get the appropriate log level for a job or component"""
        # Check job-specific level first
        if job_type and job_type in self.job_levels:
            return self.job_levels[job_type]

        # Check component-specific level
        if component and component in self.component_levels:
            return self.component_levels[component]

        # Return global level
        return self.global_level


class JobLogger:
    """Logger that includes job context in log messages"""

    def __init__(
        self,
        logger: logging.Logger,
        job_id: str = None,
        job_type: str = None,
        filepath: str = None,
    ):
        self.logger = logger
        self.job_id = job_id
        self.job_type = job_type
        self.filepath = filepath

    def _format_message(self, level: str, message: Any) -> str:
        """Format log message with job context"""
        parts = []

        if self.job_type:
            parts.append(f"[{self.job_type}]")

        if self.job_id:
            parts.append(f"Job-{self.job_id}")

        if self.filepath:
            filename = os.path.basename(self.filepath)
            parts.append(f"({filename})")

        if parts:
            return f"{' '.join(parts)} {message}"
        return message

    def debug(self, message: Any, *args, **kwargs) -> None:
        """Log debug message"""
        self.logger.debug(self._format_message("DEBUG", message), *args, **kwargs)

    def info(self, message: Any, *args, **kwargs) -> None:
        """Log info message"""
        self.logger.info(self._format_message("INFO", message), *args, **kwargs)

    def warning(self, message: Any, *args, **kwargs) -> None:
        """Log warning message"""
        self.logger.warning(self._format_message("WARNING", message), *args, **kwargs)

    def error(self, message: Any, *args, **kwargs) -> None:
        """Log error message"""
        self.logger.error(self._format_message("ERROR", message), *args, **kwargs)

    def critical(self, message: Any, *args, **kwargs) -> None:
        """Log critical message"""
        self.logger.critical(self._format_message("CRITICAL", message), *args, **kwargs)

    def exception(self, message: Any, *args, **kwargs) -> None:
        """Log an exception message with job context."""
        self.logger.exception(self._format_message("ERROR", message), *args, **kwargs)


class LogManager:
    """Manages logging configuration and logger creation"""

    def __init__(self, config: LogConfig = None):
        self.config = config or LogConfig()
        self._loggers: Dict[str, logging.Logger] = {}
        self._setup_root_logger()

    def _setup_root_logger(self) -> None:
        """Setup the root logger with basic configuration"""
        root_logger = logging.getLogger()
        root_logger.setLevel(getattr(logging, self.config.global_level.upper()))

        # Remove existing handlers to avoid duplicates
        for handler in root_logger.handlers[:]:
            root_logger.removeHandler(handler)

        # Create console handler
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(getattr(logging, self.config.global_level.upper()))

        # Create formatter
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        console_handler.setFormatter(formatter)

        # Add handler to root logger
        root_logger.addHandler(console_handler)

        # Create file handler
        file_handler = logging.FileHandler("robin.log")
        file_handler.setLevel(getattr(logging, self.config.global_level.upper()))
        file_handler.setFormatter(formatter)

        # Add file handler to root logger
        root_logger.addHandler(file_handler)

    def get_logger(self, name: str, level: str = None) -> logging.Logger:
        """Get or create a logger with the specified name and level"""
        if name not in self._loggers:
            logger = logging.getLogger(name)
            if level:
                logger.setLevel(getattr(logging, level.upper()))
            self._loggers[name] = logger

        return self._loggers[name]

    def get_job_logger(
        self, job_id: str, job_type: str, filepath: str = None, level: str = None
    ) -> JobLogger:
        """Get a job-specific logger with context"""
        # Determine log level
        if level is None:
            level = self.config.get_level(job_type=job_type)

        # Get or create logger
        logger_name = f"robin.{job_type}"
        logger = self.get_logger(logger_name, level)

        return JobLogger(logger, job_id, job_type, filepath)

    def set_job_level(self, job_type: str, level: str) -> None:
        """Set log level for a specific job type"""
        self.config.job_levels[job_type] = level
        # Update existing logger if it exists
        logger_name = f"robin.{job_type}"
        if logger_name in self._loggers:
            self._loggers[logger_name].setLevel(getattr(logging, level.upper()))

    def set_component_level(self, component: str, level: str) -> None:
        """Set log level for a specific component"""
        self.config.component_levels[component] = level
        # Update existing logger if it exists
        logger_name = f"robin.{component}"
        if logger_name in self._loggers:
            self._loggers[logger_name].setLevel(getattr(logging, level.upper()))

    def set_global_level(self, level: str) -> None:
        """Set global log level"""
        self.config.global_level = level
        self._setup_root_logger()


# Global log manager instance
_log_manager: Optional[LogManager] = None


def get_log_manager() -> LogManager:
    """Get the global log manager instance"""
    global _log_manager
    if _log_manager is None:
        _log_manager = LogManager()
    return _log_manager


def configure_logging(
    global_level: str = "INFO",
    job_levels: Dict[str, str] = None,
    component_levels: Dict[str, str] = None,
) -> None:
    """Configure the global logging system"""
    global _log_manager

    config = LogConfig(
        global_level=global_level,
        job_levels=job_levels or {},
        component_levels=component_levels or {},
    )

    _log_manager = LogManager(config)


def get_job_logger(job_id: str, job_type: str, filepath: str = None) -> JobLogger:
    """Get a job-specific logger"""
    return get_log_manager().get_job_logger(job_id, job_type, filepath)


@contextmanager
def temporary_log_level(logger_name: str, level: str):
    """Temporarily change log level for a logger"""
    log_manager = get_log_manager()
    logger = log_manager.get_logger(logger_name)
    original_level = logger.level

    try:
        logger.setLevel(getattr(logging, level.upper()))
        yield logger
    finally:
        logger.setLevel(original_level)


# Convenience functions for common log levels
def set_debug_level(job_type: str = None, component: str = None) -> None:
    """Set debug level for a job type or component"""
    log_manager = get_log_manager()
    if job_type:
        log_manager.set_job_level(job_type, "DEBUG")
    if component:
        log_manager.set_component_level(component, "DEBUG")


def set_info_level(job_type: str = None, component: str = None) -> None:
    """Set info level for a job type or component"""
    log_manager = get_log_manager()
    if job_type:
        log_manager.set_job_level(job_type, "INFO")
    if component:
        log_manager.set_component_level(component, "INFO")


def set_warning_level(job_type: str = None, component: str = None) -> None:
    """Set warning level for a job type or component"""
    log_manager = get_log_manager()
    if job_type:
        log_manager.set_job_level(job_type, "WARNING")
    if component:
        log_manager.set_component_level(component, "WARNING")


def set_error_level(job_type: str = None, component: str = None) -> None:
    """Set error level for a job type or component"""
    log_manager = get_log_manager()
    if job_type:
        log_manager.set_job_level(job_type, "ERROR")
    if component:
        log_manager.set_component_level(component, "ERROR")
