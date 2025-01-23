from typing import Set, Union
from pathlib import Path
import logging


class WatchlistManager:
    def __init__(self):
        self._folders = set()
        self._callbacks = []

    def register_callback(self, callback):
        """Register a callback to be called when the watchlist changes"""
        if callback not in self._callbacks:
            self._callbacks.append(callback)

    def add_folder(self, folder):
        """Add a folder to the watchlist and notify callbacks"""
        if folder and folder not in self._folders:
            self._folders.add(folder)
            self._notify_callbacks()

    def remove_folder(self, folder):
        """Remove a folder from the watchlist and notify callbacks"""
        if folder in self._folders:
            self._folders.remove(folder)
            self._notify_callbacks()

    def _notify_callbacks(self):
        """Notify all registered callbacks of the current folder list"""
        for callback in self._callbacks:
            try:
                callback(self._folders)
            except Exception as e:
                logging.error(f"Error in watchlist callback: {e}")

    @property
    def folders(self):
        """Get the current set of watched folders"""
        return self._folders.copy()
