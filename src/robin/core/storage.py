from typing import Optional, Dict, Any
from nicegui import app


class StorageManager:
    """Manages storage operations for the application.

    This class handles all storage-related operations including initialization,
    configuration, and state management.
    """

    def __init__(self, mainuuid: str):
        """Initialize the storage manager.

        Args:
            mainuuid: Unique identifier for the current session
        """
        self.mainuuid = mainuuid
        self._initialize_storage()

    def _initialize_storage(self):
        """Initialize the storage structure."""
        if self.mainuuid not in app.storage.general:
            app.storage.general[self.mainuuid] = {}

    def configure_storage(self, sample_id: str) -> None:
        """Configure storage for a specific sample.

        Args:
            sample_id: Identifier for the sample
        """
        if self.mainuuid in app.storage.general:
            if "samples" not in app.storage.general[self.mainuuid]:
                app.storage.general[self.mainuuid]["samples"] = {}
            if sample_id not in app.storage.general[self.mainuuid]["samples"]:
                app.storage.general[self.mainuuid]["samples"][sample_id] = {
                    "files": [],
                    "metrics": {},
                    "status": "initializing",
                }

    def get_sample_data(self, sample_id: str) -> Optional[Dict[str, Any]]:
        """Get data for a specific sample.

        Args:
            sample_id: Identifier for the sample

        Returns:
            Dictionary containing sample data if found, None otherwise
        """
        if self.mainuuid in app.storage.general:
            if "samples" in app.storage.general[self.mainuuid]:
                return app.storage.general[self.mainuuid]["samples"].get(sample_id)
        return None

    def update_sample_status(self, sample_id: str, status: str) -> None:
        """Update the status of a specific sample.

        Args:
            sample_id: Identifier for the sample
            status: New status for the sample
        """
        if self.mainuuid in app.storage.general:
            if "samples" in app.storage.general[self.mainuuid]:
                if sample_id in app.storage.general[self.mainuuid]["samples"]:
                    app.storage.general[self.mainuuid]["samples"][sample_id][
                        "status"
                    ] = status

    def add_file_to_sample(self, sample_id: str, file_path: str) -> None:
        """Add a file to a sample's file list.

        Args:
            sample_id: Identifier for the sample
            file_path: Path to the file to add
        """
        if self.mainuuid in app.storage.general:
            if "samples" in app.storage.general[self.mainuuid]:
                if sample_id in app.storage.general[self.mainuuid]["samples"]:
                    if (
                        file_path
                        not in app.storage.general[self.mainuuid]["samples"][sample_id][
                            "files"
                        ]
                    ):
                        app.storage.general[self.mainuuid]["samples"][sample_id][
                            "files"
                        ].append(file_path)

    def clear_storage(self) -> None:
        """Clear all storage data for the current session."""
        if self.mainuuid in app.storage.general:
            del app.storage.general[self.mainuuid]
