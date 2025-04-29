from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List


@dataclass
class BrainMethConfig:
    """Configuration for the BrainMeth application.

    This dataclass holds all configuration parameters for the BrainMeth application,
    providing a centralized way to manage application settings.
    """

    threads: int = 4
    force_sampleid: Optional[str] = None
    kit: Optional[str] = None
    centreID: Optional[str] = None
    simtime: bool = False
    watchfolder: Optional[Path] = None
    output: Optional[Path] = None
    sequencing_summary: Optional[Path] = None
    target_panel: Optional[str] = None
    showerrors: bool = False
    browse: bool = False
    exclude: List[str] = None
    reference: Optional[Path] = None
    bed_file: Optional[Path] = None
    mainuuid: Optional[str] = None
    readfish_toml: Optional[Path] = None
    basecall_config: Optional[Path] = None
    experiment_duration: Optional[int] = None
    unique_id: Optional[str] = None
    telemetry_instance: Optional[object] = None

    def __post_init__(self):
        """Initialize default values for mutable attributes."""
        if self.exclude is None:
            self.exclude = []
