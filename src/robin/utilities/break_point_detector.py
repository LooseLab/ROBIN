import numpy as np
from collections import defaultdict
from typing import Dict, List, Tuple


class CNVChangeDetectorTracker:
    """
    Tracks and analyzes copy number variation (CNV) breakpoints in genomic data.

    This class handles the detection and tracking of CNV breakpoints across different chromosomes.
    It uses a binning approach to manage large genomic data efficiently, allowing for the analysis
    of CNV patterns by calculating a threshold based on scaled genomic proportions.

    Attributes:
        base_proportion (float): The base proportion of the genome considered for CNV analysis.
        bin_width (int): The bin width used for segmenting the genome.
        max_bin_width (int): The maximum bin width used in the analysis.
        coordinates (Dict[str, Dict[str, any]]): A dictionary storing breakpoint data and analysis results.
    """

    def __init__(self, base_proportion: float = 0.01, bin_width: int = 1_000_000):
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
        self.coordinates: Dict[str, Dict[str, any]] = (
            {}
        )  # Stores CNV data per chromosome

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

        # Convert the increments dictionary to NumPy arrays for efficient processing
        positions_array = np.array(list(increments.keys()))
        values_array = np.array(list(increments.values()))

        # Sort positions to ensure proper ordering of breakpoints
        sorted_indices = np.argsort(positions_array)
        sorted_positions = positions_array[sorted_indices]
        sorted_values = values_array[sorted_indices]

        # Initialize chromosome data if not already present
        self._initialize_chromosome(chromosome, length)

        # Process positions and calculate cumulative values
        positions = self.coordinates[chromosome]["positions"]
        cumulative_values = np.cumsum(sorted_values)
        final_positions = np.vstack((sorted_positions, cumulative_values)).T

        # Store the processed positions and values
        positions.extend(map(tuple, final_positions))

        # Recalculate the threshold and update BED data
        self._calculate_threshold(chromosome)
        self.coordinates[chromosome]["bed_data"] = self.get_bed_targets(chromosome)

    def _initialize_chromosome(self, chromosome: str, length: int) -> None:
        """
        Initializes internal data structures for tracking a new chromosome.

        Args:
            chromosome (str): The name of the chromosome to initialize.
            length (int): The total length of the chromosome.
        """
        self.coordinates[chromosome] = {
            "positions": [],
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
            for pos, depth in data["positions"]:
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
                f"{chrom}\t{start}\t{end}\t.\t.\t+"
                f"\n{chrom}\t{start}\t{end}\t.\t.\t-"
                for chrom, start, end in bedlist
            ]
            return "\n".join(bed_lines)
        return ""
