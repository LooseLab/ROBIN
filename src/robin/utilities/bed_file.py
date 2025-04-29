"""
File included for development purposes.

Needs to be refactored.

"""

# Python imports.
from __future__ import annotations
from nicegui import ui, app

from robin import theme
from collections import defaultdict
import numpy as np
import os

import csv
from natsort import natsorted
from io import StringIO
import copy
import hashlib
from datetime import datetime
import pandas as pd

from typing import Optional
import tomli
import tomli_w
from pathlib import Path


class MasterBedTree:
    """
    A class for managing and processing BED (Browser Extensible Data) files.

    This class is used to manage multiple BedTree instances.
    It provides functionality to:
    - Create and add a new BedTree instances
    - Remove BedTree instances
    - Get the results from a specific BedTree instance for a specific sample
    """

    def __init__(
        self,
        default_preserve_original_tree: bool = False,
        default_reference_file: str = None,
        default_output_location: str = None,
        default_readfish_toml: Optional[Path] = None,
    ):
        """
        Initialize a new MasterBedTree instance.

        Args:
            default_preserve_original_tree (bool, optional): Default value for preserve_original_tree. Defaults to False.
            default_reference_file (str, optional): Default path to reference genome FASTA index file. Defaults to None.
            default_output_location (str, optional): Default directory for output files. Defaults to None.
            default_readfish_toml (Optional[Path], optional): Default path to readfish TOML configuration. Defaults to None.
        """
        self.bed_trees = defaultdict(lambda: None, key=None)
        self.default_preserve_original_tree = default_preserve_original_tree
        self.default_reference_file = default_reference_file
        self.default_output_location = default_output_location
        self.default_readfish_toml = default_readfish_toml

    def __getitem__(self, sample_id: str) -> Optional[BedTree]:
        """
        Get a BedTree instance for a specific sample.

        Args:
            sample_id (str): The ID of the sample to get the BedTree for

        Returns:
            Optional[BedTree]: The BedTree instance if found, None otherwise
        """
        return self.bed_trees[sample_id]

    def add_bed_tree(
        self,
        sample_id: str,
        preserve_original_tree: bool = None,
        reference_file: str = None,
        output_location: str = None,
        readfish_toml: Optional[Path] = None,
    ):
        """
        Add a new BedTree instance for a specific sample.

        Args:
            sample_id (str): The ID of the sample to create the BedTree for
            preserve_original_tree (bool, optional): Whether to preserve the original tree structure.
                If None, uses the default from MasterBedTree initialization. Defaults to None.
            reference_file (str, optional): Path to reference genome FASTA index file.
                If None, uses the default from MasterBedTree initialization. Defaults to None.
            output_location (str, optional): Directory for output files.
                If None, uses the default from MasterBedTree initialization. Defaults to None.
            readfish_toml (Optional[Path], optional): Path to readfish TOML configuration.
                If None, uses the default from MasterBedTree initialization. Defaults to None.
        """
        # Use provided values or fall back to defaults
        preserve_original_tree = (
            preserve_original_tree
            if preserve_original_tree is not None
            else self.default_preserve_original_tree
        )
        reference_file = (
            reference_file
            if reference_file is not None
            else self.default_reference_file
        )
        output_location = (
            output_location
            if output_location is not None
            else self.default_output_location
        )
        readfish_toml = (
            readfish_toml if readfish_toml is not None else self.default_readfish_toml
        )

        bed_tree = BedTree(
            preserve_original_tree=preserve_original_tree,
            reference_file=reference_file,
            output_location=output_location,
            readfish_toml=readfish_toml,
        )
        bed_tree.sample_id = sample_id
        self.bed_trees[sample_id] = bed_tree

    def remove_bed_tree(self, sample_id: str):
        del self.bed_trees[sample_id]


class BedTree:
    """
    A class for managing and processing BED (Browser Extensible Data) files.

    This class provides functionality for:
    - Loading and processing BED files
    - Managing target regions and their metadata
    - Tracking proportions and statistics
    - Writing updated BED files
    - Preserving original tree structures

    Attributes:
        total_count (int): Total number of target regions
        total_range_sum (int): Total number of bases covered by targets
        tree_data (dict): Hierarchical data structure for targets
        tree_dict (dict): Dictionary representation of the tree structure
        preserve_original_tree (bool): Whether to preserve the original tree structure
        reference_tree (dict): Copy of the original tree structure
        file_counter (int): Counter for file naming
        previous_targets_hash (str): Hash of previous targets for change detection
        proportions_df (pd.DataFrame): DataFrame storing proportions and timestamps
        output_location (str): Location for output files
        readfish_toml (Optional[Path]): Path to readfish TOML configuration
        toml_dict (dict): Parsed TOML configuration
        chromosome_lengths (dict): Dictionary of chromosome lengths
        total_length (int): Total length of all chromosomes

    Example:
        >>> bed_tree = BedTree(preserve_original_tree=True, reference_file="reference.fa.fai")
        >>> bed_tree.load_from_file("targets.bed")
    """

    def __init__(
        self,
        preserve_original_tree=False,
        reference_file=None,
        output_location=None,
        readfish_toml: Optional[Path] = None,
    ):
        """
        Initialize a new BedTree instance.

        Args:
            preserve_original_tree (bool, optional): Whether to preserve the original tree structure. Defaults to False.
            reference_file (str, optional): Path to reference genome FASTA index file. Defaults to None.
            output_location (str, optional): Directory for output files. Defaults to None.
            readfish_toml (Optional[Path], optional): Path to readfish TOML configuration. Defaults to None.

        Returns:
            None

        Example:
            >>> bed_tree = BedTree(
            ...     preserve_original_tree=True,
            ...     reference_file="reference.fa.fai",
            ...     output_location="output/"
            ... )
        """
        self.total_count = 0
        self.total_range_sum = 0
        self.tree_data = {
            "id": "Chromosomes",
            "description": "Summary of currently used targets.",
            "children": [],
            "count": 0,
            "range_sum": 0,
            "proportion": 0,
        }
        self.tree_dict = {}
        self.preserve_original_tree = preserve_original_tree
        self.reference_tree = None
        self.file_counter = 0  # Counter for the file naming
        self.previous_targets_hash = None  # To track changes in targets
        self.proportions_df = (
            pd.DataFrame()
        )  # DataFrame to store the proportions and timestamps
        self.output_location = output_location
        self.readfish_toml: Optional[Path] = readfish_toml
        if self.readfish_toml:
            with open(self.readfish_toml, "rb") as f:
                self.toml_dict = tomli.load(f)
        else:
            self.toml_dict = None
        if reference_file:
            self.chromosome_lengths = self._get_chromosome_lengths(reference_file)
            self.total_length = sum(self.chromosome_lengths.values())
        else:
            self.chromosome_lengths = None
            self.total_length = None

    def _hash_current_targets(self):
        """
        Generate a hash for the current targets to detect changes.

        Returns:
            str: MD5 hash of the current target structure

        Example:
            >>> hash_value = bed_tree._hash_current_targets()
            >>> print(f"Current hash: {hash_value}")
        """
        m = hashlib.md5()
        for chromosome, chromosome_data in sorted(self.tree_dict.items()):
            m.update(chromosome.encode())
            for strand_group in sorted(
                chromosome_data["children"], key=lambda x: x["id"]
            ):
                m.update(strand_group["id"].encode())
                for target in sorted(strand_group["children"], key=lambda x: x["id"]):
                    m.update(target["id"].encode())
        return m.hexdigest()

    def _check_and_create_folder(self, path, folder_name=None):
        """
        Check if a path exists and optionally create a subfolder.

        Args:
            path (str): Base path to check
            folder_name (str, optional): Name of subfolder to create. Defaults to None.

        Returns:
            str: Full path to the folder

        Raises:
            FileNotFoundError: If the base path does not exist

        Example:
            >>> folder_path = bed_tree._check_and_create_folder("/output", "results")
            >>> print(f"Created folder at: {folder_path}")
        """
        # Check if the path exists
        if not os.path.exists(path):
            raise FileNotFoundError(f"The specified path does not exist: {path}")

        # If folder_name is provided
        if folder_name:
            full_path = os.path.join(path, folder_name)
            # Create the folder if it doesn't exist
            if not os.path.exists(full_path):
                os.makedirs(full_path)
            return full_path
        else:
            return path

    def _write_bed_file(self):
        """
        Write the current targets to a new BED file if they have changed.

        This method:
        1. Checks if targets have changed using hash comparison
        2. Creates a new BED file with updated targets
        3. Updates the readfish TOML configuration if present

        Returns:
            None

        Example:
            >>> bed_tree._write_bed_file()
            >>> # A new BED file will be created if targets have changed
        """
        current_hash = self._hash_current_targets()
        if current_hash != self.previous_targets_hash:
            self.file_counter += 1
            if self.output_location:
                filename = os.path.join(
                    self._check_and_create_folder(
                        self.output_location, folder_name="bed_files"
                    ),
                    f"new_file_{self.file_counter}.bed",
                )
            else:
                filename = f"new_file_{self.file_counter}.bed"
            with open(filename, "w") as f:
                for chromosome, chromosome_data in sorted(self.tree_dict.items()):
                    for strand_group in sorted(
                        chromosome_data["children"], key=lambda x: x["id"]
                    ):
                        for target in sorted(
                            strand_group["children"], key=lambda x: x["id"]
                        ):
                            try:
                                if target["id"] and "-" in target["id"]:
                                    start, end = map(int, target["id"].split("-"))
                                    # Check if this target is from the original preserved tree
                                    is_preserved = (
                                        self.preserve_original_tree
                                        and self.reference_tree
                                        and any(
                                            any(
                                                any(
                                                    t["id"] == target["id"]
                                                    for t in sg["children"]
                                                )
                                                for sg in chrom["children"]
                                            )
                                            for chrom in self.reference_tree.values()
                                            if chrom["id"] == chromosome
                                        )
                                    )
                                    # Use "fixed_target" for preserved targets, keep original name for others
                                    name = (
                                        "fixed_target"
                                        if is_preserved
                                        else target.get("name", "unknown")
                                    )
                                    strand = strand_group["id"].split()[
                                        -1
                                    ]  # Get the strand from the group ID
                                    f.write(
                                        f"{chromosome}\t{start}\t{end}\t{name}\t.\t{strand}\n"
                                    )
                            except (ValueError, IndexError):
                                # Skip invalid entries
                                continue
                self.previous_targets_hash = current_hash
                self._update_proportions_df()
                if self.toml_dict:
                    for region in self.toml_dict["regions"]:
                        region["targets"] = filename
                    with open(f"{self.readfish_toml}_live", "wb") as f:
                        tomli_w.dump(self.toml_dict, f)

    def _get_chromosome_lengths(self, fai_file):
        """
        Read chromosome lengths from a FASTA index file.

        Args:
            fai_file (str): Path to the FASTA index file

        Returns:
            dict: Dictionary mapping chromosome names to their lengths

        Example:
            >>> lengths = bed_tree._get_chromosome_lengths("reference.fa.fai")
            >>> print(f"Chromosome lengths: {lengths}")
        """
        chromosome_lengths = {}
        with open(fai_file, "r") as f:
            for line in f:
                parts = line.split("\t")
                chromosome = parts[0]
                length = int(parts[1])
                chromosome_lengths[chromosome] = length
        return chromosome_lengths

    def _merge_ranges(self, existing_ranges, new_range):
        """
        Merge ranges while preserving those with different names.

        Args:
            existing_ranges (list): List of tuples (start, end, name)
            new_range (tuple): Tuple (start, end, name)

        Returns:
            list: List of merged ranges sorted by start position

        Example:
            >>> existing = [(100, 200, "gene1"), (300, 400, "gene2")]
            >>> new = (150, 250, "gene1")
            >>> merged = bed_tree._merge_ranges(existing, new)
            >>> print(f"Merged ranges: {merged}")
        """
        start, end, name = new_range
        merged_ranges = []
        i = 0

        while i < len(existing_ranges):
            existing_start, existing_end, existing_name = existing_ranges[i]

            # Only merge if names match or both are default (".")
            if name == existing_name or (name == "." and existing_name == "."):
                if not (end < existing_start or start > existing_end):
                    start = min(existing_start, start)
                    end = max(existing_end, end)
                else:
                    merged_ranges.append((existing_start, existing_end, existing_name))
            else:
                # Keep ranges separate if names differ
                merged_ranges.append((existing_start, existing_end, existing_name))
            i += 1

        merged_ranges.append((start, end, name))
        # Sort by start position, then end position if starts are equal
        merged_ranges.sort(key=lambda x: (x[0], x[1]))

        return merged_ranges

    def _process_bed_data(self, reader):
        """
        Process BED data from a CSV reader.

        Args:
            reader (csv.reader): CSV reader containing BED data

        Returns:
            None

        Example:
            >>> with open("targets.bed") as f:
            ...     reader = csv.reader(f, delimiter="\t")
            ...     bed_tree._process_bed_data(reader)
        """
        for row in reader:
            if len(row) > 0:
                chromosome = row[0]
                start = int(row[1])
                end = int(row[2])
                name = row[3]
                strand = row[5] if len(row) > 5 else "."

                if chromosome not in self.tree_dict:
                    self.tree_dict[chromosome] = {
                        "id": chromosome,
                        "children": [],
                        "count": 0,
                        "range_sum": 0,
                        "description": "chromosome",
                    }

                strand_id = f"{chromosome} {strand}"
                strand_group = next(
                    (
                        child
                        for child in self.tree_dict[chromosome]["children"]
                        if child["id"] == strand_id
                    ),
                    None,
                )
                if not strand_group:
                    strand_group = {
                        "id": strand_id,
                        "children": [],
                        "count": 0,
                        "range_sum": 0,
                    }
                    self.tree_dict[chromosome]["children"].append(strand_group)

                # Include name in the range information
                existing_ranges = []
                for r in strand_group.get("children", []):
                    try:
                        if r["id"] and "-" in r["id"]:
                            start_pos, end_pos = map(int, r["id"].split("-"))
                            existing_ranges.append((start_pos, end_pos, r["name"]))
                    except (ValueError, IndexError):
                        # Skip invalid entries
                        continue

                merged_ranges = self._merge_ranges(existing_ranges, (start, end, name))

                strand_group["children"] = [
                    {
                        "id": f"{s}-{e}",
                        "description": "target region",
                        "chromosome_length": self.chromosome_lengths.get(chromosome, 0),
                        "range_sum": e - s + 1,
                        "proportion": (e - s + 1)
                        / self.chromosome_lengths.get(chromosome, 0)
                        * 100,
                        "name": n,
                    }
                    for s, e, n in merged_ranges
                ]

                self.tree_dict[chromosome]["count"] = sum(
                    c["count"] for c in self.tree_dict[chromosome]["children"]
                )
                self.tree_dict[chromosome]["range_sum"] = sum(
                    c["range_sum"] for c in self.tree_dict[chromosome]["children"]
                )

        self.total_count = sum(v["count"] for v in self.tree_dict.values())
        self.total_range_sum = sum(v["range_sum"] for v in self.tree_dict.values())

    def _add_if_not_exists(self, entries, new_item, description, chromosome):
        """
        Add a new entry if it doesn't already exist in the entries list.

        Args:
            entries (list): List of existing entries
            new_item (str): ID of the new item
            description (str): Description of the new item
            chromosome (str): Chromosome name

        Returns:
            None

        Example:
            >>> bed_tree._add_if_not_exists(
            ...     entries=existing_entries,
            ...     new_item="chr1",
            ...     description="Chromosome 1",
            ...     chromosome="chr1"
            ... )
        """
        if not any(entry["id"] == new_item for entry in entries):
            entries.append(
                {
                    "id": new_item,
                    "children": [],
                    "description": description,
                    "count": 0,
                    "range_sum": 0,
                    "proportion": 0,
                    "chromosome_length": self.chromosome_lengths.get(chromosome, 0),
                }
            )

    def _propagate_values(self, node):
        """
        Recursively propagate count and range sum values through the tree.

        Args:
            node (dict): Current node in the tree

        Returns:
            tuple: (total_count, total_range_sum) for the node

        Example:
            >>> count, range_sum = bed_tree._propagate_values(root_node)
            >>> print(f"Total count: {count}, Total range: {range_sum}")
        """
        if (
            "children" not in node or not node["children"]
        ):  # Here we are looking at an individual target on a chromosome.
            node["count"] = node.get("count", 1)
            node["range_sum"] = node.get("range_sum", 0)
            return node["count"], node["range_sum"]

        total_count = 0
        total_range_sum = 0

        for child in node[
            "children"
        ]:  # Here we are not in a tiny node - we are in a bigger node.
            child_count, child_range_sum = self._propagate_values(child)
            total_count += child_count
            total_range_sum += child_range_sum
        node["count"] = total_count
        node["range_sum"] = total_range_sum

        if (
            node["description"] == "Stranded Targets"
        ):  # We are looking at a stranded chromosome.
            node["proportion"] = (
                node["range_sum"] / node["chromosome_length"] * 100
                if self.total_length
                else 0
            )
        elif node["description"] == "chromosome":
            node["proportion"] = (
                node["range_sum"] / (2 * node["chromosome_length"]) * 100
                if self.total_length
                else 0
            )
        elif node["description"] == "Summary of currently used targets.":
            node["proportion"] = (
                node["range_sum"] / (2 * self.total_length) * 100
                if self.total_length
                else 0
            )

        return node["count"], node["range_sum"]

    def _build_tree(self):
        """
        Build the hierarchical tree structure from the tree dictionary.

        This method:
        1. Adds chromosomes to the tree
        2. Sorts strand groups and their children
        3. Updates strands and targets
        4. Propagates values through the tree

        Returns:
            None

        Example:
            >>> bed_tree._build_tree()
            >>> # The tree structure will be built and values propagated
        """
        for chromosome_data in self.tree_dict.values():
            self._add_if_not_exists(
                self.tree_data["children"],
                chromosome_data["id"],
                "chromosome",
                chromosome_data["id"],
            )

            for strand_group in natsorted(chromosome_data["children"]):
                try:
                    strand_group["children"] = natsorted(
                        strand_group["children"],
                        key=lambda x: (
                            int(x["id"].split("-")[0])
                            if x["id"] and "-" in x["id"]
                            else 0
                        ),
                    )
                except (ValueError, IndexError):
                    # If conversion fails, sort by the original ID
                    strand_group["children"] = natsorted(
                        strand_group["children"], key=lambda x: x["id"]
                    )
                self._update_strands(chromosome_data["id"], strand_group["id"])
                self._update_targets(
                    chromosome_data["id"], strand_group["id"], strand_group
                )

        self.tree_data["children"] = natsorted(
            self.tree_data["children"], key=lambda x: x["id"]
        )

        self._propagate_values(self.tree_data)

    def _update_strands(self, chromosome, strand):
        """
        Update the strands for a given chromosome.

        Args:
            chromosome (str): Chromosome name
            strand (str): Strand identifier

        Returns:
            None

        Example:
            >>> bed_tree._update_strands("chr1", "chr1 +")
        """
        chromosome_entry = next(
            item for item in self.tree_data["children"] if item["id"] == chromosome
        )
        self._add_if_not_exists(
            chromosome_entry["children"], strand, "Stranded Targets", chromosome
        )

    def _update_targets(self, chromosome, strand, strand_group):
        """
        Update the targets for a given chromosome and strand.

        Args:
            chromosome (str): Chromosome name
            strand (str): Strand identifier
            strand_group (dict): Strand group data

        Returns:
            None

        Example:
            >>> bed_tree._update_targets("chr1", "chr1 +", strand_data)
        """
        chromosome_entry = next(
            item for item in self.tree_data["children"] if item["id"] == chromosome
        )
        strand_entry = next(
            orientation
            for orientation in chromosome_entry["children"]
            if orientation["id"] == strand
        )
        strand_entry["children"] = strand_group["children"]
        strand_entry["count"] = strand_group["count"]
        strand_entry["range_sum"] = strand_group["range_sum"]
        strand_entry["chromosome_length"] = self.chromosome_lengths.get(
            chromosome, 0
        )  # Add chromosome length

    def load_from_file(self, file_path, merge=False):
        """
        Load BED data from a file.

        Args:
            file_path (str): Path to the BED file
            merge (bool, optional): Whether to merge with existing data. Defaults to False.

        Returns:
            None

        Example:
            >>> bed_tree.load_from_file("targets.bed")
            >>> # The BED file will be loaded and processed
        """
        if not merge:
            if "children" in self.tree_data.keys():
                self.tree_data["children"] = []
            self.tree_dict = {}
        with open(file_path, "r") as bed_file:
            reader = csv.reader(bed_file, delimiter="\t")
            self._process_bed_data(reader)
        self._build_tree()
        self._propagate_values(self.tree_data)
        if self.preserve_original_tree:
            self.reference_tree = copy.deepcopy(self.tree_dict)
        self._propagate_values(self.tree_data)

    def load_from_string(
        self,
        bed_string,
        merge=False,
        write_files=False,
        output_location=None,
        source_type=None,
    ):
        """
        Load BED data from a string while preserving entries from other sources.

        Args:
            bed_string (str): The BED format string to load
            merge (bool, optional): Whether to merge with existing data. Defaults to False.
            write_files (bool, optional): Whether to write output files. Defaults to False.
            output_location (str, optional): Where to write output files. Defaults to None.
            source_type (str, optional): The type of source updating the BED ("CNV", "FUSION", etc.). Defaults to None.

        Returns:
            None

        Example:
            >>> bed_data = "chr1\t100\t200\tgene1\t.\t+\nchr1\t300\t400\tgene2\t.\t-"
            >>> bed_tree.load_from_string(bed_data, source_type="FUSION")
        """
        if output_location:
            self.output_location = output_location

        if not merge:
            if "children" in self.tree_data.keys():
                if source_type:
                    # First pass: Remove all entries of the type we're updating
                    for chromosome in self.tree_data["children"][:]:
                        for strand_group in chromosome.get("children", [])[:]:
                            if source_type == "FUSION":
                                # Remove all entries that are NOT CNV_detected
                                to_remove = []
                                for i, target in enumerate(
                                    strand_group.get("children", [])
                                ):
                                    if target.get("name") != "CNV_detected":
                                        to_remove.append(i)
                                # Remove entries in reverse order to maintain correct indices
                                for i in reversed(to_remove):
                                    del strand_group["children"][i]

                            elif source_type == "CNV":
                                # Remove all CNV_detected entries
                                to_remove = []
                                for i, target in enumerate(
                                    strand_group.get("children", [])
                                ):
                                    if target.get("name") == "CNV_detected":
                                        to_remove.append(i)
                                # Remove entries in reverse order to maintain correct indices
                                for i in reversed(to_remove):
                                    del strand_group["children"][i]

                            # Remove empty strand groups
                            if not strand_group.get("children", []):
                                chromosome["children"].remove(strand_group)

                    # Remove empty chromosomes
                    self.tree_data["children"] = [
                        chrom
                        for chrom in self.tree_data["children"]
                        if chrom.get("children")
                    ]
                else:
                    # If no source_type specified, clear everything as before
                    self.tree_data["children"] = []

            # Handle the tree_dict preservation
            if self.preserve_original_tree and self.reference_tree:
                if not self.tree_dict:
                    # If tree_dict is empty, simply copy the reference tree
                    self.tree_dict = copy.deepcopy(self.reference_tree)
                else:
                    # Merge reference tree entries with existing tree_dict
                    reference_copy = copy.deepcopy(self.reference_tree)
                    for chrom_key, chrom_value in reference_copy.items():
                        if chrom_key not in self.tree_dict:
                            # If chromosome doesn't exist, add it completely
                            self.tree_dict[chrom_key] = chrom_value
                        else:
                            # Chromosome exists, need to merge strand groups
                            existing_chrom = self.tree_dict[chrom_key]
                            if "children" in chrom_value:
                                if "children" not in existing_chrom:
                                    existing_chrom["children"] = []

                                # For each strand group in reference
                                for ref_strand_group in chrom_value["children"]:
                                    strand_group_id = ref_strand_group["id"]
                                    # Find matching strand group in existing chromosome
                                    existing_strand_group = next(
                                        (
                                            group
                                            for group in existing_chrom["children"]
                                            if group["id"] == strand_group_id
                                        ),
                                        None,
                                    )
                                    if existing_strand_group is None:
                                        # If strand group doesn't exist, add it
                                        existing_chrom["children"].append(
                                            ref_strand_group
                                        )
                                    else:
                                        # Merge children of strand groups
                                        if "children" in ref_strand_group:
                                            if "children" not in existing_strand_group:
                                                existing_strand_group["children"] = []
                                            # Add any missing targets
                                            existing_ids = {
                                                child["id"]
                                                for child in existing_strand_group[
                                                    "children"
                                                ]
                                            }
                                            for ref_target in ref_strand_group[
                                                "children"
                                            ]:
                                                if ref_target["id"] not in existing_ids:
                                                    existing_strand_group[
                                                        "children"
                                                    ].append(ref_target)
            elif not self.preserve_original_tree:
                self.tree_dict = {}

        bed_file = StringIO(bed_string)
        reader = csv.reader(bed_file, delimiter="\t")
        self._process_bed_data(reader)
        self._build_tree()
        self._propagate_values(self.tree_data)
        if write_files:
            self._write_bed_file()

    def _update_proportions_df(self):
        """
        Update the proportions DataFrame with the current data and timestamp.

        This method:
        1. Creates a new DataFrame with current proportions
        2. Appends it to the existing DataFrame
        3. Writes the updated DataFrame to a CSV file if output_location is set

        Returns:
            None

        Example:
            >>> bed_tree._update_proportions_df()
            >>> # The proportions DataFrame will be updated and written to file
        """
        timestamp = datetime.now()

        # Top-level (genome-wide) proportion
        top_level_proportion = {
            "timestamp": timestamp,
            "chromosome": "Genome-wide",
            "proportion": self.tree_data["proportion"],
        }

        # List to hold the rows
        rows = [top_level_proportion]

        # Proportions for each chromosome
        for chromosome in self.tree_data["children"]:
            chromosome_proportion = {
                "timestamp": timestamp,
                "chromosome": chromosome["id"],
                "proportion": chromosome["proportion"],
            }
            rows.append(chromosome_proportion)

        # Create a DataFrame from the rows and append it to the existing DataFrame
        new_df = pd.DataFrame(rows)
        self.proportions_df = pd.concat(
            [self.proportions_df, new_df], ignore_index=True
        )
        if self.output_location:
            self.proportions_df.to_csv(
                os.path.join(self.output_location, "bedranges.csv")
            )

    def update_target_table(self) -> None:
        """
        Update the target panel information table with current BedTree data.

        This method:
        1. Checks if gene_bed data is available
        2. Updates the target table with current status and source information
        3. Updates the UI if a target_table exists

        Returns:
            None

        Example:
            >>> bed_tree.update_target_table()
            >>> # The target table will be updated with current information
        """
        print("Updating target table")  # Basic debug print
        logger.debug("Starting target table update")
        if not hasattr(self, "gene_bed") or self.gene_bed.empty:
            logger.warning("No gene_bed data available")
            return

        table_rows = []
        for _, row in self.gene_bed.iterrows():
            target_status = "Original"
            target_source = "Panel"

            if self.NewBed and self.NewBed.tree_dict:
                chrom_data = self.NewBed.tree_dict.get(row.chrom, {})
                if chrom_data:
                    # Check for any overlap
                    for strand_group in chrom_data.get("children", []):
                        for target in strand_group.get("children", []):
                            target_start, target_end = map(int, target["id"].split("-"))
                            # Check for any overlap between intervals
                            if (
                                target_start <= row.end_pos
                                and target_end >= row.start_pos
                            ):
                                target_status = "Active"
                                source_name = target.get("name", "")
                                if "CNV" in source_name:
                                    target_source = "CNV_detected"
                                elif source_name:
                                    target_source = source_name
                                break
                        if target_status == "Active":
                            break

            table_rows.append(
                {
                    "chrom": str(row.chrom),
                    "gene": str(row.gene),
                    "start_pos": int(row.start_pos),
                    "end_pos": int(row.end_pos),
                    "size": int(row.end_pos - row.start_pos),
                    "status": target_status,
                    "source": target_source,
                }
            )

        if self.target_table:
            self.target_table.rows = table_rows
            ui.update(self.target_table)

    def get_unique_entries_count(self) -> dict:
        """
        Count the number of unique entries in the tree_dict.

        Returns:
            dict: A dictionary containing:
                - total_targets: Total number of unique target regions
                - by_chromosome: Dictionary of counts per chromosome
                - by_source: Dictionary of counts by source type (Panel, CNV_detected, etc.)

        Example:
            >>> counts = bed_tree.get_unique_entries_count()
            >>> print(f"Total targets: {counts['total_targets']}")
            >>> print(f"By chromosome: {counts['by_chromosome']}")
            >>> print(f"By source: {counts['by_source']}")
        """
        if not self.tree_dict:
            return {"total_targets": 0, "by_chromosome": {}, "by_source": {}}

        total_targets = 0
        by_chromosome = {}
        by_source = {}

        for chrom, chrom_data in self.tree_dict.items():
            chrom_count = 0
            for strand_group in chrom_data.get("children", []):
                for target in strand_group.get("children", []):
                    chrom_count += 1
                    # Count by source
                    source = target.get("name", "unknown")
                    by_source[source] = by_source.get(source, 0) + 1

            by_chromosome[chrom] = chrom_count
            total_targets += chrom_count

        return {
            "total_targets": total_targets,
            "by_chromosome": by_chromosome,
            "by_source": by_source,
        }


@ui.page("/", response_timeout=30)
def index_page() -> None:
    # from check_connection import ConnectionDialog
    # initial_ip = "127.0.0.1"
    # my_connection = ConnectionDialog(initial_ip)
    # my_connection = None
    existing_bed = "/Users/mattloose/GIT/robin_dev/robin/src/robin/resources/panel_11092024_5kb_pad.bed"
    existing_bed2 = "/Users/mattloose/GIT/robin_dev/robin/src/robin/resources/AML_panel_name_uniq.bed"

    master_bed_tree = MasterBedTree()

    master_bed_tree.add_bed_tree(
        sample_id="ds1305_CNVDetection_0048_a",
        preserve_original_tree=True,
        reference_file="/Users/mattloose/references/hg38_simple.fa.fai",
    )

    master_bed_tree.add_bed_tree(
        sample_id="ds1305_CNVDetection_40_b",
        preserve_original_tree=True,
        reference_file="/Users/mattloose/references/hg38_simple.fa.fai",
    )

    # original_bed_tree = BedTree(
    #    preserve_original_tree=True,
    #    reference_file="/Users/mattloose/references/hg38_simple.fa.fai",
    # )
    original_bed_tree = master_bed_tree.bed_trees["ds1305_CNVDetection_0048_a"]

    original_bed_tree.load_from_file(existing_bed)

    original_bed_tree2 = master_bed_tree.bed_trees["ds1305_CNVDetection_40_b"]
    original_bed_tree2.load_from_file(existing_bed2)

    original_bed_tree3 = master_bed_tree.bed_trees["camel"]

    if original_bed_tree3:
        print("original_bed_tree3 is not None")
    else:
        print("original_bed_tree3 is None")

    output = "/Users/mattloose/GIT/robin_dev/robin/single_sample_results/ds1305_CNVDetection_0048_a/"
    output2 = "/Users/mattloose/GIT/robin_dev/robin/single_sample_results/ds1305_CNVDetection_40_b/"
    CNVResults = np.load(
        os.path.join(output, "ruptures.npy"), allow_pickle="TRUE"
    ).item()

    CNVResults2 = np.load(
        os.path.join(output2, "ruptures.npy"), allow_pickle="TRUE"
    ).item()

    with theme.frame(
        "<strong><font color='#000000'>R</font></strong>apid nanop<strong><font color='#000000'>O</font></strong>re <strong><font color='#000000'>B</font></strong>rain intraoperat<strong><font color='#000000'>I</font></strong>ve classificatio<strong><font color='#000000'>N</font></strong>",
        smalltitle="<strong><font color='#000000'>R.O.B.I.N</font></strong>",
    ):
        with ui.card().classes("w-full"):

            ui.label("New Target Information")
            ui.label(f"Sample:ID = {original_bed_tree.sample_id}")
            orig_tree = ui.tree(
                [original_bed_tree.tree_data],
                label_key="id",
            )
            orig_tree.add_slot(
                "default-body",
                """
                <div v-if="props.node.description">
                        <span class="text-weight-bold">{{ props.node.description }}</span>
                    </div>
                <div v-if="props.node.range_sum">
                        <span class="text-weight-bold">{{props.node.range_sum}} bases</span>
                </div>
                <div v-if="props.node.count">
                        <span class="text-weight-bold">{{ props.node.count }} targets</span>
                </div>
                <div v-if="props.node.proportion">
                        <span class="text-weight-bold">{{ props.node.proportion }} proportion</span>
                </div>
                <div v-if="props.node.chromosome_length">
                        <span class="text-weight-bold">{{ props.node.chromosome_length }} chromosome length</span>
                </div>
        
            """,
            )

        with ui.card().classes("w-full"):

            ui.label("New Target Information")
            ui.label(f"Sample:ID = {original_bed_tree2.sample_id}")
            orig_tree2 = ui.tree(
                [original_bed_tree2.tree_data],
                label_key="id",
            )
            orig_tree2.add_slot(
                "default-body",
                """
                <div v-if="props.node.description">
                        <span class="text-weight-bold">{{ props.node.description }}</span>
                    </div>
                <div v-if="props.node.range_sum">
                        <span class="text-weight-bold">{{props.node.range_sum}} bases</span>
                </div>
                <div v-if="props.node.count">
                        <span class="text-weight-bold">{{ props.node.count }} targets</span>
                </div>
                <div v-if="props.node.proportion">
                        <span class="text-weight-bold">{{ props.node.proportion }} proportion</span>
                </div>
                <div v-if="props.node.chromosome_length">
                        <span class="text-weight-bold">{{ props.node.chromosome_length }} chromosome length</span>
                </div>
        
            """,
            )


def run_class(port: int, reload: bool):
    """
    Helper function to run the app.
    :param port: The port to serve the app on.
    :param reload: Should we reload the app on changes.
    :return:
    """
    ui.add_css(
        """
        .shadows-into light-regular {
            font-family: "Shadows Into Light", cursive;
            font-weight: 400;
            font-style: normal;
        }
    """
    )
    app.add_static_files("/fonts", str(Path(__file__).parent / "../fonts"))
    # app.on_startup(mainpage.index_page)
    # index_page()
    ui.run(
        port=port,
        reload=reload,
        title="MethClass NiceGUI",
        storage_secret="slartibartfast",
        native=True,
        window_size=(1200, 800),
        fullscreen=False,
    )  # , native=True, fullscreen=False, window_size=(1200, 1000))


def main():  # , threads, simtime, watchfolder, output, sequencing_summary):
    """
    Entrypoint for when GUI is launched directly.
    :return: None
    """
    run_class(port=12398, reload=False)


# Entrypoint for when GUI is launched by the CLI.
# e.g.: python my_app/my_cli.py
if __name__ in {"__main__", "__mp_main__"}:
    """
    Entrypoint for when GUI is launched by the CLI
    :return: None
    """
    if __name__ == "__mp_main__":
        print("GUI launched by auto-reload")

    run_class(port=12398, reload=True)
