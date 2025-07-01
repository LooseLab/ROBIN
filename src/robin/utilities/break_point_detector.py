import numpy as np
from collections import defaultdict
from typing import Dict, List, Tuple
import mmap
import tempfile
import struct
import array
from pathlib import Path
import json
import shutil
import pickle


class MemoryMappedStorage:
    """
    Handles efficient storage of genomic position data using memory-mapped files.
    This class provides a way to store large arrays of genomic positions on disk
    while maintaining efficient access patterns.
    """

    def __init__(self, temp_dir: str = None):
        self.temp_dir = temp_dir or tempfile.gettempdir()
        self.storage_dir = Path(self.temp_dir) / "robin_cnv_storage"
        self.storage_dir.mkdir(exist_ok=True)
        self.mapped_files = {}
        self.position_arrays = {}

    def _get_file_path(self, chromosome: str) -> Path:
        return self.storage_dir / f"{chromosome}_positions.bin"

    def store_positions(
        self, chromosome: str, positions: List[Tuple[int, int]]
    ) -> None:
        """Store position data for a chromosome using memory mapping."""
        if not positions:
            return

        file_path = self._get_file_path(chromosome)
        # Convert positions to a flat array of integers
        flat_data = array.array("Q")  # Using unsigned long long (8 bytes) for positions
        for pos, val in positions:
            flat_data.extend([pos, val])

        # Write to file
        with open(file_path, "wb") as f:
            # Write number of positions as header
            f.write(struct.pack("Q", len(positions)))
            flat_data.tofile(f)

    def load_positions(self, chromosome: str) -> List[Tuple[int, int]]:
        """Load position data for a chromosome using memory mapping."""
        file_path = self._get_file_path(chromosome)
        if not file_path.exists():
            return []

        with open(file_path, "rb") as f:
            # Read number of positions from header
            num_positions = struct.unpack("Q", f.read(8))[0]
            # Memory map the file
            mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            # Read the data
            positions = []
            for _ in range(num_positions):
                pos = struct.unpack("Q", mm.read(8))[0]
                val = struct.unpack("Q", mm.read(8))[0]
                positions.append((pos, val))
            mm.close()
        return positions

    def clear(self, chromosome: str = None) -> None:
        """Clear stored data for a specific chromosome or all chromosomes."""
        if chromosome:
            file_path = self._get_file_path(chromosome)
            if file_path.exists():
                file_path.unlink()
        else:
            for file in self.storage_dir.glob("*_positions.bin"):
                file.unlink()

    def save_state(self, output_path: str) -> None:
        """Save the current state of all stored data to a directory."""
        output_dir = Path(output_path)
        output_dir.mkdir(exist_ok=True)

        # Save metadata about stored chromosomes
        metadata = {
            "temp_dir": str(self.temp_dir),
            "chromosomes": [
                f.stem.split("_")[0] for f in self.storage_dir.glob("*_positions.bin")
            ],
        }

        with open(output_dir / "storage_metadata.pkl", "wb") as f:
            pickle.dump(metadata, f)

        # Copy all position files to the output directory
        for file in self.storage_dir.glob("*_positions.bin"):
            shutil.copy2(file, output_dir / file.name)

    def load_state(self, input_path: str) -> None:
        """Load state from a saved directory."""
        input_dir = Path(input_path)
        if not input_dir.exists():
            raise FileNotFoundError(f"State directory not found: {input_path}")

        # Load metadata
        try:
            with open(input_dir / "storage_metadata.pkl", "rb") as f:
                metadata = pickle.load(f)
        except (FileNotFoundError, pickle.UnpicklingError):
            # Fallback to JSON if pickle file doesn't exist or is corrupted
            try:
                with open(input_dir / "storage_metadata.json", "r") as f:
                    metadata = json.load(f)
            except (FileNotFoundError, json.JSONDecodeError):
                metadata = {"temp_dir": str(self.temp_dir), "chromosomes": []}
                print(f"Error loading state from {input_path}")
                print(f"Metadata: {metadata}")

        # Clear existing data
        self.clear()

        # Copy position files to storage directory
        for file in input_dir.glob("*_positions.bin"):
            shutil.copy2(file, self.storage_dir / file.name)


class CNVChangeDetectorTracker:
    """
    Tracks and analyzes copy number variation (CNV) breakpoints in genomic data.
    Uses memory-mapped storage for efficient memory usage.

    This class handles the detection and tracking of CNV breakpoints across different chromosomes.
    It uses a binning approach to manage large genomic data efficiently, allowing for the analysis
    of CNV patterns by calculating a threshold based on scaled genomic proportions.

    Attributes:
        base_proportion (float): The base proportion of the genome considered for CNV analysis.
        bin_width (int): The bin width used for segmenting the genome.
        max_bin_width (int): The maximum bin width used in the analysis.
        coordinates (Dict[str, Dict[str, any]]): A dictionary storing breakpoint data and analysis results.
    """

    def __init__(
        self,
        base_proportion: float = 0.01,
        bin_width: int = 1_000_000,
        temp_dir: str = None,
    ):
        """
        Initializes the CNVChangeDetectorTracker with a specified base proportion and bin width.

        Args:
            base_proportion (float): The proportion of the genome to be used as a base for CNV analysis.
            bin_width (int): The width of the bins to be used for genomic segmentation.
        """
        self.base_proportion: float = base_proportion
        self.bin_width: int = bin_width
        self.max_bin_width: int = (
            bin_width  # Used to track the largest bin width during updates
        )

        self.storage = MemoryMappedStorage(temp_dir)

        # Store only essential metadata in memory
        self.coordinates: Dict[str, Dict[str, any]] = (
            {}
        )  # Stores CNV data per chromosome
        self.ruptures_breakpoints = {}

    def _calculate_scaled_proportion(self, length: int) -> int:
        """
        Calculates a scaled proportion of the genome based on the current bin width.

        The proportion scales with the bin width, allowing for adaptive analysis depending on the
        granularity of the data.

        Args:
            length (int): The length of the chromosome or genomic region.

        Returns:
            int: The scaled proportion of the genome.
        """
        # The proportion is scaled by a factor dependent on the bin width
        proportion = self.base_proportion * ((self.bin_width / 1_000_000) ** (1 / 4))
        return int(length * proportion)

    def add_breakpoints(
        self,
        chromosome: str,
        breakpoints: List[Dict[str, int]],
        local_breakpoints: List[Dict[str, int]],
        length: int,
        bin_width: int,
    ) -> None:
        """
        Adds detected CNV breakpoints to the tracker and updates internal structures accordingly.

        This method handles new breakpoints for a given chromosome, ensuring that internal structures
        are updated and the analysis is recalculated as necessary.

        Args:
            chromosome (str): The name of the chromosome (e.g., 'chr1').
            breakpoints (List[Dict[str, int]]): A list of breakpoint dictionaries with 'start' and 'end' positions.
            length (int): The total length of the chromosome or genomic region.
            bin_width (int): The width of the bins used for this set of breakpoints.
        """
        if bin_width != self.bin_width:
            self._update_bin_width(bin_width)

        increments = defaultdict(int)
        for entry in breakpoints:
            increments[entry["start"]] += 1
            increments[entry["end"]] -= 1

        breakpointlist = []
        for entry in local_breakpoints:
            breakpointlist.append((entry["start"], entry["end"]))

        self.ruptures_breakpoints[chromosome] = breakpointlist

        # Process positions and store in memory-mapped file
        positions_array = np.array(list(increments.keys()))
        values_array = np.array(list(increments.values()))
        sorted_indices = np.argsort(positions_array)
        sorted_positions = positions_array[sorted_indices]
        sorted_values = values_array[sorted_indices]

        # Store positions in memory-mapped file
        positions = list(zip(sorted_positions, np.cumsum(sorted_values)))
        self.storage.store_positions(chromosome, positions)

        # Initialize chromosome metadata if needed
        if chromosome not in self.coordinates:
            self._initialize_chromosome(chromosome, length)

        # Recalculate threshold and update BED data
        self._calculate_threshold(chromosome)
        self.coordinates[chromosome]["bed_data"] = self.get_bed_targets(chromosome)
        self.coordinates[chromosome]["bed_data_breakpoints"] = (
            self.get_bed_targets_breakpoints(chromosome)
        )

    def _initialize_chromosome(self, chromosome: str, length: int) -> None:
        """
        Initializes internal data structures for tracking a new chromosome.

        Args:
            chromosome (str): The name of the chromosome to initialize.
            length (int): The total length of the chromosome.
        """
        self.coordinates[chromosome] = {
            "length": length,
            "target": self._calculate_scaled_proportion(length),
            "threshold": 0,
            "current_value": 0,
            "flush": True,
            "previous_value": 0,
            "start_value": 0,
            "end_value": 0,
            "cumulative_sum": 0,
            "start_positions": [],
            "end_positions": [],
            "bed_data": "",
        }

    def _update_chromosome_target(self, chromosome: str) -> None:
        """
        Updates the target proportion for a chromosome based on the current bin width.

        Args:
            chromosome (str): The name of the chromosome.
        """
        if chromosome in self.coordinates:
            length = self.coordinates[chromosome]["length"]
            self.coordinates[chromosome]["target"] = self._calculate_scaled_proportion(
                length
            )
            self._calculate_threshold(chromosome)

    def _update_bin_width(self, new_bin_width: int) -> None:
        """
        Updates the bin width and recalculates the target proportion for all chromosomes.

        Args:
            new_bin_width (int): The new bin width to be applied.
        """
        if new_bin_width < self.max_bin_width:
            self.bin_width = new_bin_width
            for chromosome in self.coordinates:
                self._update_chromosome_target(chromosome)

    def _calculate_threshold(self, chromosome: str) -> None:
        """
        Calculates the threshold for detecting significant CNVs for a specific chromosome.

        The threshold is incrementally adjusted until the cumulative sum of positions meets the target proportion.

        Args:
            chromosome (str): The name of the chromosome.
        """
        data = self.coordinates[chromosome]
        length = data["length"]
        target = data["target"]

        while length > target:
            self._init_positions(chromosome)
            # Load positions from memory-mapped storage
            positions = self.storage.load_positions(chromosome)
            for pos, depth in positions:
                self.update(chromosome, 1 if depth > data["threshold"] else 0, pos)
            length = self.get_cumulative_sum(chromosome)
            data["threshold"] += 1

    def _init_positions(self, chromosome: str) -> None:
        """
        Resets internal tracking variables for a chromosome before recalculating thresholds.
        """
        data = self.coordinates[chromosome]
        data.update(
            {
                "current_value": 0,
                "flush": True,
                "previous_value": 0,
                "start_value": 0,
                "end_value": 0,
                "cumulative_sum": 0,
                "start_positions": [],
                "end_positions": [],
            }
        )

    def get_threshold(self, chromosome: str) -> int:
        """
        Retrieves the current threshold used for detecting significant CNVs on a chromosome.

        Args:
            chromosome (str): The name of the chromosome.

        Returns:
            int: The current threshold value for CNV detection.
        """
        return self.coordinates[chromosome]["threshold"]

    def get_cumulative_sum(self, chromosome: str) -> int:
        """
        Calculates the cumulative sum of CNVs for a chromosome, which reflects the total length of regions with CNVs.

        Args:
            chromosome (str): The name of the chromosome.

        Returns:
            int: The cumulative sum of detected CNVs.
        """
        data = self.coordinates[chromosome]
        return (
            data["cumulative_sum"]
            if data["flush"]
            else data["cumulative_sum"] + (data["end_value"] - data["start_value"] + 1)
        )

    def update(self, chromosome: str, new_value: int, position: int) -> None:
        """
        Updates the internal tracking variables for a chromosome based on the latest position and CNV state.

        Args:
            chromosome (str): The name of the chromosome.
            new_value (int): The new CNV state at the given position (0 or 1).
            position (int): The genomic position where the update occurs.
        """
        data = self.coordinates[chromosome]

        if new_value != data["current_value"]:
            if new_value == 1 and data["current_value"] == 0:
                # Starting a new CNV region
                data["start_value"] = position
                data["end_value"] = position
                data["start_positions"].append(position)
                data["flush"] = False
            elif new_value == 0 and data["current_value"] == 1:
                # Ending the current CNV region
                data["end_value"] = position
                data["cumulative_sum"] += data["end_value"] - data["start_value"] + 1
                data["end_positions"].append(position)
                data["flush"] = False
            data["previous_value"] = data["current_value"]
            data["current_value"] = new_value
        elif new_value == 1:
            # Extending the current CNV region
            data["end_value"] = position

    def get_breakpoints(self, chromosome: str) -> Tuple[List[int], List[int]]:
        """
        Retrieves the start and end positions of all detected CNV breakpoints for a chromosome.

        Args:
            chromosome (str): The name of the chromosome.

        Returns:
            Tuple[List[int], List[int]]: Two lists containing start and end positions of CNV breakpoints.
        """
        data = self.coordinates[chromosome]
        return data["start_positions"], data["end_positions"]

    def get_bed_targets(self, chromosome: str) -> str:
        """
        Generates BED format data for detected CNV regions on a chromosome.

        Args:
            chromosome (str): The name of the chromosome.

        Returns:
            str: A string containing BED format lines for each CNV region.
        """
        if chromosome in self.coordinates:
            data = self.coordinates[chromosome]
            bedlist = [
                (chromosome, start, end)
                for start, end in zip(data["start_positions"], data["end_positions"])
            ]
            bed_lines = [
                f"{chrom}\t{start - self.bin_width}\t{end}\.\t.\t+"
                f"\n{chrom}\t{start}\t{end + self.bin_width}\.\t.\t-"
                for chrom, start, end in bedlist
            ]
            return "\n".join(bed_lines)
        return ""

    def get_bed_targets_breakpoints(self, chromosome: str) -> str:
        """
        Generates BED format data for detected CNV regions on a chromosome.

        Args:
            chromosome (str): The name of the chromosome.

        Returns:
            str: A string containing BED format lines for each CNV region.
        """
        if chromosome in self.ruptures_breakpoints:
            data = self.ruptures_breakpoints[chromosome]
            bedlist = [(chromosome, start, end) for start, end in data]
            bed_lines = [
                f"{chrom}\t{start - int(1.5*self.bin_width)}\t{end}\tCNV_detected\t.\t+"
                f"\n{chrom}\t{start}\t{end + int(1.5*self.bin_width)}\tCNV_detected\t.\t-"
                for chrom, start, end in bedlist
            ]
            return "\n".join(bed_lines)
        return ""

    def __del__(self):
        """Clean up memory-mapped files when the object is destroyed."""
        self.storage.clear()

    def save_state(self, output_path: str) -> None:
        """Save the current state of the tracker."""
        # Save the memory-mapped storage state
        self.storage.save_state(output_path)

        # Convert coordinates to a serializable format
        serializable_coordinates = {}
        for chrom, data in self.coordinates.items():
            serializable_coordinates[chrom] = {
                k: float(v) if isinstance(v, np.float32) else v
                for k, v in data.items()
                if k
                not in [
                    "positions"
                ]  # Skip positions as they're stored in memory-mapped files
            }

        # Save the metadata using pickle for better handling of complex types
        metadata = {
            "base_proportion": self.base_proportion,
            "bin_width": self.bin_width,
            "max_bin_width": self.max_bin_width,
            "coordinates": serializable_coordinates,
            "ruptures_breakpoints": dict(
                self.ruptures_breakpoints
            ),  # Convert to regular dict
        }

        with open(Path(output_path) / "tracker_metadata.pkl", "wb") as f:
            pickle.dump(metadata, f)

    def load_state(self, input_path: str) -> None:
        """Load state from a saved directory."""
        # Load the memory-mapped storage state
        self.storage.load_state(input_path)

        # Load the metadata
        try:
            with open(Path(input_path) / "tracker_metadata.pkl", "rb") as f:
                metadata = pickle.load(f)
        except (FileNotFoundError, pickle.UnpicklingError):
            # Fallback to JSON if pickle file doesn't exist or is corrupted
            try:
                with open(Path(input_path) / "tracker_metadata.json", "r") as f:
                    metadata = json.load(f)
            except (FileNotFoundError, json.JSONDecodeError):
                # If both files are missing or corrupted, initialize with defaults
                metadata = {
                    "base_proportion": self.base_proportion,
                    "bin_width": self.bin_width,
                    "max_bin_width": self.max_bin_width,
                    "coordinates": {},
                    "ruptures_breakpoints": {},
                }

        self.base_proportion = metadata["base_proportion"]
        self.bin_width = metadata["bin_width"]
        self.max_bin_width = metadata["max_bin_width"]
        self.coordinates = metadata["coordinates"]
        self.ruptures_breakpoints = metadata["ruptures_breakpoints"]
