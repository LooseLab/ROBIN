from typing import Dict, Any, Optional
from dataclasses import dataclass
from nicegui import app

@dataclass
class MetricsState:
    """State container for various metrics in the application."""
    bam_count: Dict[str, Any]
    file_metrics: Dict[str, Any]
    read_metrics: Dict[str, Any]

class MetricsManager:
    """Manages metrics and counters for the application.
    
    This class handles all metrics-related operations including initialization,
    updates, and state management.
    """
    
    def __init__(self, mainuuid: str):
        """Initialize the metrics manager.
        
        Args:
            mainuuid: Unique identifier for the current session
        """
        self.mainuuid = mainuuid
        self._initialize_storage()
    
    def _initialize_storage(self):
        """Initialize the storage structure for metrics."""
        if self.mainuuid not in app.storage.general:
            app.storage.general[self.mainuuid] = {}
        if "metrics" not in app.storage.general[self.mainuuid]:
            app.storage.general[self.mainuuid]["metrics"] = MetricsState(
                bam_count={"counter": 0, "total_files": 0, "file": {}},
                file_metrics={},
                read_metrics={}
            )
    
    def update_counter(self, counter_name: str, value: Any) -> None:
        """Update a specific counter in the metrics state.
        
        Args:
            counter_name: Name of the counter to update
            value: New value for the counter
        """
        if self.mainuuid in app.storage.general:
            if "metrics" not in app.storage.general[self.mainuuid]:
                self._initialize_storage()
            app.storage.general[self.mainuuid]["metrics"].bam_count[counter_name] = value
    
    def get_counter(self, counter_name: str) -> Any:
        """Get the value of a specific counter.
        
        Args:
            counter_name: Name of the counter to retrieve
            
        Returns:
            The current value of the counter
        """
        if self.mainuuid in app.storage.general:
            if "metrics" in app.storage.general[self.mainuuid]:
                return app.storage.general[self.mainuuid]["metrics"].bam_count.get(counter_name)
        return None
    
    def update_file_metrics(self, sample_id: str, metrics: Dict[str, Any]) -> None:
        """Update metrics for a specific sample.
        
        Args:
            sample_id: Identifier for the sample
            metrics: Dictionary of metrics to update
        """
        if self.mainuuid in app.storage.general:
            if "metrics" not in app.storage.general[self.mainuuid]:
                self._initialize_storage()
            app.storage.general[self.mainuuid]["metrics"].file_metrics[sample_id] = metrics
    
    def update_read_metrics(self, sample_id: str, metrics: Dict[str, Any]) -> None:
        """Update read metrics for a specific sample.
        
        Args:
            sample_id: Identifier for the sample
            metrics: Dictionary of read metrics to update
        """
        if self.mainuuid in app.storage.general:
            if "metrics" not in app.storage.general[self.mainuuid]:
                self._initialize_storage()
            app.storage.general[self.mainuuid]["metrics"].read_metrics[sample_id] = metrics 