#!/usr/bin/env python3
"""
Memory Management Module for Robin Ray Workflow

This module provides memory management utilities to help prevent memory leaks
in long-running Ray actors. It implements garbage collection triggers based
on both iteration count and RSS memory usage.

Features:
- Automatic garbage collection based on iteration count
- RSS memory monitoring with configurable thresholds
- Optional malloc_trim() for returning freed memory to OS
- Configurable parameters for different workload types
- Integration with Ray actors and analysis handlers

Classes
-------
MemoryManager
    Main memory management class that monitors and triggers cleanup.

Usage
-----
.. code-block:: python

    from robin.memory_manager import MemoryManager

    # Initialize with default settings
    memory_manager = MemoryManager()

    # In your processing loop
    for batch in batches:
        result = process_batch(batch)
        memory_manager.check_and_cleanup()
        return result

    # Or with custom settings
    memory_manager = MemoryManager(
        gc_every=25,           # GC every 25 iterations
        rss_trigger_mb=1024,  # Trigger GC at 1GB RSS
        enable_malloc_trim=True
    )

Authors
-------
Matt Loose
"""

import gc
import logging
import os
from typing import Any, Dict, Optional

import psutil


class MemoryManager:
    """
    Memory management utility for long-running Ray actors.

    Monitors memory usage and triggers garbage collection based on:
    1. Iteration count (every N iterations)
    2. RSS memory threshold (when memory exceeds limit)

    Also optionally attempts to return freed memory to the OS using malloc_trim().
    """

    def __init__(
        self,
        gc_every: int = 50,
        rss_trigger_mb: int = 2048,
        enable_malloc_trim: bool = True,
        restart_every: int = 10000,
        restart_rss_trigger_mb: int = 4096,
        logger: Optional[logging.Logger] = None,
    ):
        """
        Initialize the memory manager.

        Args:
            gc_every: Trigger garbage collection every N iterations
            rss_trigger_mb: Trigger garbage collection when RSS exceeds this many MB
            enable_malloc_trim: Whether to attempt malloc_trim() after GC
            restart_every: Trigger actor restart every N iterations
            restart_rss_trigger_mb: Trigger actor restart when RSS exceeds this many MB
            logger: Optional logger for memory management events
        """
        self._iteration_count = 0
        self._gc_every = max(1, gc_every)
        self._rss_trigger_bytes = rss_trigger_mb * 1024 * 1024
        self._enable_malloc_trim = enable_malloc_trim
        self._logger = logger or logging.getLogger(__name__)

        # Actor restart mechanism disabled - keeping only memory cleanup
        self._restart_every = restart_every  # Keep for stats but not used
        self._restart_rss_trigger_bytes = (
            restart_rss_trigger_mb * 1024 * 1024
        )  # Keep for stats but not used
        self._restart_requested = False
        self._restart_reason = None

        # Get process handle for memory monitoring
        try:
            self._process = psutil.Process(os.getpid())
        except Exception as e:
            self._logger.warning(f"Failed to initialize psutil process handle: {e}")
            self._process = None

        # Initialize malloc_trim if enabled
        self._malloc_trim_func = None
        if self._enable_malloc_trim:
            try:
                import ctypes
                import ctypes.util

                libc = ctypes.CDLL(ctypes.util.find_library("c"))
                self._malloc_trim_func = libc.malloc_trim
                self._logger.debug("malloc_trim() available for memory optimization")
            except Exception as e:
                self._logger.debug(f"malloc_trim() not available: {e}")
                # Don't change the setting, just disable the function
                self._malloc_trim_func = None

        self._logger.info(
            f"MemoryManager initialized: gc_every={self._gc_every}, "
            f"rss_trigger={rss_trigger_mb}MB, malloc_trim={self._enable_malloc_trim}, "
            f"actor_restarts=DISABLED (memory cleanup only)"
        )

    def check_and_cleanup(
        self, force: bool = False, is_idle: bool = True
    ) -> Dict[str, Any]:
        """
        Check memory usage and trigger cleanup if needed.

        Args:
            force: Force garbage collection regardless of thresholds
            is_idle: Whether the actor is currently idle (no active jobs)

        Returns:
            Dictionary with cleanup statistics
        """
        self._iteration_count += 1

        # Get current memory info
        memory_info = self._get_memory_info()
        rss_bytes = memory_info.get("rss_bytes", 0)

        # Determine if cleanup is needed
        # Only perform cleanup if actor is idle (unless forced)
        needs_gc = force or (
            is_idle
            and (
                (self._iteration_count % self._gc_every == 0)
                or (rss_bytes > self._rss_trigger_bytes)
            )
        )

        cleanup_stats = {
            "iteration": self._iteration_count,
            "rss_mb": memory_info.get("rss_mb", 0),
            "gc_triggered": False,
            "malloc_trim_attempted": False,
            "trigger_reason": None,
            "is_idle": is_idle,
            "skipped_due_to_busy": False,
            "restart_requested": False,
            "restart_reason": None,
        }

        if needs_gc:
            cleanup_stats["gc_triggered"] = True

            # Determine trigger reason
            if force:
                cleanup_stats["trigger_reason"] = "forced"
            elif self._iteration_count % self._gc_every == 0:
                cleanup_stats["trigger_reason"] = "iteration_count"
            elif rss_bytes > self._rss_trigger_bytes:
                cleanup_stats["trigger_reason"] = "rss_threshold"

            # Perform garbage collection
            try:
                collected = gc.collect()
                cleanup_stats["objects_collected"] = collected
                self._logger.debug(f"Garbage collection: {collected} objects collected")
            except Exception as e:
                self._logger.warning(f"Garbage collection failed: {e}")
                cleanup_stats["gc_error"] = str(e)

            # Attempt malloc_trim if enabled
            if self._enable_malloc_trim and self._malloc_trim_func:
                try:
                    # malloc_trim(0) returns freed memory to OS
                    result = self._malloc_trim_func(0)
                    cleanup_stats["malloc_trim_attempted"] = True
                    cleanup_stats["malloc_trim_result"] = result
                    self._logger.debug(f"malloc_trim() returned: {result}")
                except Exception as e:
                    self._logger.debug(f"malloc_trim() failed: {e}")
                    cleanup_stats["malloc_trim_error"] = str(e)

            # Log cleanup event
            self._logger.info(
                f"Memory cleanup triggered ({cleanup_stats['trigger_reason']}): "
                f"RSS={cleanup_stats['rss_mb']:.1f}MB, "
                f"iterations={self._iteration_count}, "
                f"idle={is_idle}"
            )
        elif not is_idle and not force:
            # Log that cleanup was skipped due to busy state
            cleanup_stats["skipped_due_to_busy"] = True
            self._logger.debug(
                f"Memory cleanup skipped: actor busy (RSS={cleanup_stats['rss_mb']:.1f}MB, "
                f"iterations={self._iteration_count})"
            )

        # Actor restart mechanism disabled - keeping only memory cleanup
        # This allows actors to run indefinitely while still managing memory
        cleanup_stats["restart_requested"] = False
        cleanup_stats["restart_reason"] = None

        return cleanup_stats

    def _get_memory_info(self) -> Dict[str, Any]:
        """Get current memory usage information."""
        if self._process is None:
            return {"rss_bytes": 0, "rss_mb": 0, "error": "process_not_available"}

        try:
            memory_info = self._process.memory_info()
            rss_bytes = memory_info.rss
            rss_mb = rss_bytes / (1024 * 1024)

            return {
                "rss_bytes": rss_bytes,
                "rss_mb": rss_mb,
                "vms_bytes": memory_info.vms,
                "vms_mb": memory_info.vms / (1024 * 1024),
            }
        except Exception as e:
            self._logger.warning(f"Failed to get memory info: {e}")
            return {"rss_bytes": 0, "rss_mb": 0, "error": str(e)}

    def get_stats(self) -> Dict[str, Any]:
        """Get current memory manager statistics."""
        memory_info = self._get_memory_info()

        return {
            "iteration_count": self._iteration_count,
            "gc_every": self._gc_every,
            "rss_trigger_mb": self._rss_trigger_bytes / (1024 * 1024),
            "enable_malloc_trim": self._enable_malloc_trim,
            "restart_every": self._restart_every,
            "restart_rss_trigger_mb": self._restart_rss_trigger_bytes / (1024 * 1024),
            "restart_requested": self._restart_requested,
            "restart_reason": self._restart_reason,
            "current_memory": memory_info,
            "next_gc_at": self._gc_every - (self._iteration_count % self._gc_every),
            "next_restart_at": self._restart_every
            - (self._iteration_count % self._restart_every),
        }

    def reset_iteration_count(self):
        """Reset the iteration counter."""
        self._iteration_count = 0
        self._logger.debug("Iteration count reset")

    def is_restart_requested(self) -> bool:
        """Check if actor restart has been requested."""
        return self._restart_requested

    def get_restart_reason(self) -> Optional[str]:
        """Get the reason for restart request."""
        return self._restart_reason

    def clear_restart_request(self):
        """Clear the restart request (call after restart)."""
        self._restart_requested = False
        self._restart_reason = None
        self._logger.info("Actor restart request cleared")

    def update_settings(
        self,
        gc_every: Optional[int] = None,
        rss_trigger_mb: Optional[int] = None,
        enable_malloc_trim: Optional[bool] = None,
        restart_every: Optional[int] = None,
        restart_rss_trigger_mb: Optional[int] = None,
    ):
        """Update memory manager settings."""
        if gc_every is not None:
            self._gc_every = max(1, gc_every)
            self._logger.info(f"Updated gc_every to {self._gc_every}")

        if rss_trigger_mb is not None:
            self._rss_trigger_bytes = rss_trigger_mb * 1024 * 1024
            self._logger.info(f"Updated rss_trigger to {rss_trigger_mb}MB")

        if enable_malloc_trim is not None:
            self._enable_malloc_trim = enable_malloc_trim
            self._logger.info(
                f"Updated enable_malloc_trim to {self._enable_malloc_trim}"
            )

        if restart_every is not None:
            self._restart_every = max(1, restart_every)
            self._logger.info(f"Updated restart_every to {self._restart_every}")

        if restart_rss_trigger_mb is not None:
            self._restart_rss_trigger_bytes = restart_rss_trigger_mb * 1024 * 1024
            self._logger.info(
                f"Updated restart_rss_trigger to {restart_rss_trigger_mb}MB"
            )

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit with forced cleanup."""
        self.check_and_cleanup(force=True, is_idle=True)


# Convenience function for quick integration
def create_memory_manager(
    gc_every: int = 50,
    rss_trigger_mb: int = 2048,
    enable_malloc_trim: bool = True,
    restart_every: int = 10000,
    restart_rss_trigger_mb: int = 4096,
    logger: Optional[logging.Logger] = None,
) -> MemoryManager:
    """
    Create a MemoryManager instance with the specified settings.

    Args:
        gc_every: Trigger garbage collection every N iterations
        rss_trigger_mb: Trigger garbage collection when RSS exceeds this many MB
        enable_malloc_trim: Whether to attempt malloc_trim() after GC
        restart_every: Trigger actor restart every N iterations
        restart_rss_trigger_mb: Trigger actor restart when RSS exceeds this many MB
        logger: Optional logger for memory management events

    Returns:
        Configured MemoryManager instance
    """
    return MemoryManager(
        gc_every=gc_every,
        rss_trigger_mb=rss_trigger_mb,
        enable_malloc_trim=enable_malloc_trim,
        restart_every=restart_every,
        restart_rss_trigger_mb=restart_rss_trigger_mb,
        logger=logger,
    )


# Context manager for automatic cleanup
class MemoryManagedContext:
    """Context manager that automatically triggers memory cleanup on exit."""

    def __init__(self, memory_manager: MemoryManager, force_cleanup: bool = True):
        self.memory_manager = memory_manager
        self.force_cleanup = force_cleanup

    def __enter__(self):
        return self.memory_manager

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.force_cleanup:
            self.memory_manager.check_and_cleanup(force=True)
