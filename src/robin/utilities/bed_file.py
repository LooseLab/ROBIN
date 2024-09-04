"""
File included for development purposes.

Needs to be refactored.

"""

# Python imports.
from __future__ import annotations
from nicegui import ui, app
from pathlib import Path

import theme
import numpy as np
import os
from natsort import natsorted
import csv

from io import StringIO


class BedTree:
    def __init__(self):
        self.tree_dict = {}
        self.total_count = 0
        self.total_range_sum = 0
        self.tree_data = []
        self.original_tree_data = None

    def _reset(self):
        """Resets the tree data, count, and range sum."""
        self.tree_dict = {}
        self.total_count = 0
        self.total_range_sum = 0
        self.tree_data = []

    def _merge_ranges(self, existing_ranges, new_range, chromosome, strand):
        """Merges a new range with existing ranges if they overlap, with logging."""
        start, end = new_range
        merged_ranges = []
        i = 0

        while i < len(existing_ranges):
            existing_start, existing_end = existing_ranges[i]
            if not (end < existing_start or start > existing_end):  # Check if overlapping
                # Merge the ranges
                print(
                    f"Merging range {start}-{end} with existing range {existing_start}-{existing_end} on {chromosome} {strand}")
                start = min(existing_start, start)
                end = max(existing_end, end)
            else:
                # If there's no overlap, keep the existing range
                merged_ranges.append((existing_start, existing_end))
            i += 1

        # Add the new merged range
        merged_ranges.append((start, end))

        # Sort the merged ranges by start position
        merged_ranges.sort()

        return merged_ranges

    def _process_bed_data(self, reader):
        """Processes the BED data from a CSV reader."""
        for row in reader:
            chromosome = row[0]
            start = int(row[1])
            end = int(row[2])
            strand = row[5] if len(row) > 5 else '.'  # If strand is not available, use '.' as default

            # Calculate the range (end - start)
            target_range = end - start

            # Create the ID for the range including the strand
            range_id = f'{start}-{end}'

            # Check if the chromosome is already in the dictionary
            if chromosome not in self.tree_dict:
                self.tree_dict[chromosome] = {'id': chromosome, 'children': [], 'count': 0, 'range_sum': 0}

            # Create the id for the strand with the chromosome name
            strand_id = f'{chromosome} {strand}'

            # Check if the strand is already a child of the chromosome
            strand_group = next((child for child in self.tree_dict[chromosome]['children'] if child['id'] == strand_id),
                                None)
            if not strand_group:
                strand_group = {'id': strand_id, 'children': [], 'count': 0, 'range_sum': 0}
                self.tree_dict[chromosome]['children'].append(strand_group)

            # Merge the new range into the existing ranges within the specific chromosome and strand
            existing_ranges = [(int(r['id'].split('-')[0]), int(r['id'].split('-')[1])) for r in
                               strand_group['children']]
            merged_ranges = self._merge_ranges(existing_ranges, (start, end), chromosome, strand)

            # Update the strand group with merged ranges
            strand_group['children'] = [{'id': f'{s}-{e}', 'range_sum': e - s} for s, e in merged_ranges]

            # Recalculate the counts and range sums for this specific chromosome and strand
            strand_group['count'] = len(merged_ranges)
            strand_group['range_sum'] = sum(e - s for s, e in merged_ranges)
            self.tree_dict[chromosome]['count'] = sum(c['count'] for c in self.tree_dict[chromosome]['children'])
            self.tree_dict[chromosome]['range_sum'] = sum(
                c['range_sum'] for c in self.tree_dict[chromosome]['children'])

            # Recalculate the total count and range sums across all chromosomes and strands
            self.total_count = sum(v['count'] for v in self.tree_dict.values())
            self.total_range_sum = sum(v['range_sum'] for v in self.tree_dict.values())

    def load_from_file(self, file_path, merge=False):
        """Loads and processes a BED file from the given file path."""
        if not merge:
            self._reset()
        else:
            self._preserve_original_tree()

        with open(file_path, 'r') as bed_file:
            reader = csv.reader(bed_file, delimiter='\t')
            self._process_bed_data(reader)
        self._build_tree()

    def load_from_string(self, bed_string, merge=False):
        """Loads and processes a BED file from a string."""
        if not merge:
            self._reset()
        else:
            self._preserve_original_tree()

        bed_file = StringIO(bed_string)
        reader = csv.reader(bed_file, delimiter='\t')
        self._process_bed_data(reader)
        self._build_tree()

    def merge_with_original(self, bed_string):
        """Merges new BED data with the original tree only."""
        if self.original_tree_data is None:
            raise ValueError("No original tree data available to merge with. Load initial data first.")

        # Temporarily store the current tree
        current_tree_dict = self.tree_dict

        # Restore the original tree structure
        self.tree_dict = self._restore_original_tree_dict()

        # Process the new BED data and merge it with the original tree
        bed_file = StringIO(bed_string)
        reader = csv.reader(bed_file, delimiter='\t')
        self._process_bed_data(reader)

        # Build the tree with the merged data
        self._build_tree()

        # Optionally, you could restore the current tree if you want to keep it
        # self.tree_dict = current_tree_dict

    def _restore_original_tree_dict(self):
        """Helper method to restore the original tree dictionary."""
        restored_tree_dict = {}
        for chromosome_data in self.original_tree_data[0]['children']:
            chromosome = chromosome_data['id']
            restored_tree_dict[chromosome] = {
                'id': chromosome,
                'children': [],
                'count': chromosome_data['count'],
                'range_sum': chromosome_data['range_sum']
            }
            for strand_group in chromosome_data['children']:
                restored_tree_dict[chromosome]['children'].append({
                    'id': strand_group['id'],
                    'children': list(strand_group['children']),
                    'count': strand_group['count'],
                    'range_sum': strand_group['range_sum']
                })
        return restored_tree_dict

    def _preserve_original_tree(self):
        """Preserve the original tree before any updates."""
        if self.original_tree_data is None:
            self.original_tree_data = self.get_tree().copy()

    def _build_tree(self):
        """Builds the tree structure from the processed data."""
        for chromosome_data in self.tree_dict.values():
            for strand_group in chromosome_data['children']:
                # Sort children within each strand group by the start position
                strand_group['children'] = sorted(strand_group['children'], key=lambda x: int(x['id'].split('-')[0]))

        original_bed_tree = [
            {
                'id': k,
                'children': v['children'],
                'count': v['count'],
                'range_sum': v['range_sum']
            }
            for k, v in natsorted(self.tree_dict.items())
        ]

        # Add the top-level node
        self.tree_data = [
            {
                "id": "chromosome",
                "description": "Targets On Each Chromosome",
                "count": self.total_count,
                "range_sum": self.total_range_sum,
                "children": original_bed_tree,
            }
        ]

    def get_tree(self):
        """Returns the current tree structure with sorted coordinates."""
        self._build_tree()  # Ensure the tree is built with sorted data
        return self.tree_data

    def save_to_bed_file(self, file_path):
        """Saves the current tree structure to a sorted BED file."""
        bed_entries = []

        # Collect all BED entries
        for chromosome, data in self.tree_dict.items():
            for strand_group in data['children']:
                strand = strand_group['id'].split()[1]
                for range_item in strand_group['children']:
                    start, end = map(int, range_item['id'].split('-'))
                    bed_entries.append((chromosome, start, end, '.', '.', strand))

        # Sort the entries by chromosome, start position, and strand
        sorted_bed_entries = natsorted(bed_entries, key=lambda x: (x[0], x[1], x[5]))

        # Write the sorted entries to the BED file
        with open(file_path, 'w', newline='') as bed_file:
            writer = csv.writer(bed_file, delimiter='\t')
            for entry in sorted_bed_entries:
                writer.writerow(entry)


def bed_to_tree_dict_with_ranges_and_chromosome_strand(bed_file_path):
    # Initialize an empty dictionary to hold the tree structure
    tree_dict = {}
    total_count = 0  # To keep track of the total count of all ranges
    total_range_sum = 0  # To keep track of the total range sum of all ranges

    # Open and parse the BED file
    with open(bed_file_path, 'r') as bed_file:
        reader = csv.reader(bed_file, delimiter='\t')
        for row in reader:
            chromosome = row[0]
            start = int(row[1])
            end = int(row[2])
            strand = row[5] if len(row) > 5 else '.'  # If strand is not available, use '.' as default

            # Calculate the range (end - start)
            target_range = end - start

            # Create the ID for the range including the strand
            range_id = f'{start}-{end}'

            # Check if the chromosome is already in the dictionary
            if chromosome not in tree_dict:
                tree_dict[chromosome] = {'id': chromosome, 'children': [], 'count': 0, 'range_sum': 0}

            # Create the id for the strand with the chromosome name
            strand_id = f'{chromosome} {strand}'

            # Check if the strand is already a child of the chromosome
            strand_group = next((child for child in tree_dict[chromosome]['children'] if child['id'] == strand_id),
                                None)
            if not strand_group:
                strand_group = {'id': strand_id, 'children': [], 'count': 0, 'range_sum': 0}
                tree_dict[chromosome]['children'].append(strand_group)

            # Add the range to the strand group
            strand_group['children'].append({'id': range_id, 'range_sum': target_range})
            # Increment the count and range sum for both the strand group and the chromosome
            strand_group['count'] += 1
            strand_group['range_sum'] += target_range
            tree_dict[chromosome]['count'] += 1
            tree_dict[chromosome]['range_sum'] += target_range
            total_count += 1  # Increment the total count
            total_range_sum += target_range  # Increment the total range sum

    # Convert the dictionary to a list of dictionaries as required by NiceGUI
    # Sort the chromosomes in natural order
    original_bed_tree = [
        {'id': k, 'children': v['children'], 'count': v['count'], 'range_sum': v['range_sum']}
        for k, v in natsorted(tree_dict.items())
    ]

    # Add the top-level node
    tree_data = [
        {
            "id": "chromosome",
            "description": "Targets On Each Chromosome",
            "count": total_count,
            "range_sum": total_range_sum,
            "children": original_bed_tree,
        }
    ]

    return tree_data

def add_child(data_structure, strand_id, new_child):
    """
    Adds a new child to the children list of the specified strand (+ or -) in the given data structure.

    Args:
        data_structure (list): The list containing dictionaries for each strand.
        strand_id (str): The strand identifier, either '+' or '-'.
        new_child (dict): The new child item to be added to the children list.

    Returns:
        bool: True if the item was added successfully, False if the strand_id was not found.
    """
    for item in data_structure:
        if item["id"] == strand_id:
            item["children"].append(new_child)
            return True
    return False


@ui.page("/", response_timeout=30)
def index_page() -> None:
    # from check_connection import ConnectionDialog
    # initial_ip = "127.0.0.1"
    # my_connection = ConnectionDialog(initial_ip)
    # my_connection = None
    existing_bed = "/Users/mattloose/GIT/niceGUI/cnsmeth/panel_adaptive_nogenenames_20122021_hg38_VGL_CHD7_5kb_pad.bed"

    original_bed_tree = bed_to_tree_dict_with_ranges_and_chromosome_strand(existing_bed)

    new_bed_tree = BedTree()
    new_bed_tree.load_from_file(existing_bed)

    new_bed_tree.load_from_string("chr1\t1\t1000000\t.\t.\t+", merge=True)

    new_bed_tree.save_to_bed_file("testingbedoutput.bed")

    output = "/Users/mattloose/GIT/niceGUI/cnsmeth/robin_output_test/_ds1305_Intraop0051_c_ULK_2_highFRA"
    CNVResults = np.load(
        os.path.join(output, "ruptures.npy"), allow_pickle="TRUE"
    ).item()
    with theme.frame(
        "<strong><font color='#000000'>R</font></strong>apid nanop<strong><font color='#000000'>O</font></strong>re <strong><font color='#000000'>B</font></strong>rain intraoperat<strong><font color='#000000'>I</font></strong>ve classificatio<strong><font color='#000000'>N</font></strong>",
        smalltitle="<strong><font color='#000000'>R.O.B.I.N</font></strong>",
    ):
        with ui.card().classes("w-full"):

            ui.label("New Target Information")
            orig_tree = ui.tree(
                new_bed_tree.get_tree()
               ,label_key="id",
                            )
            orig_tree.add_slot(
                "default-body",
                """
                <div v-if="props.node.description">
                        <span class="text-weight-bold">{{ props.node.description }}.</span>
                    </div>
                <div v-if="props.node.count">
                        <span class="text-weight-bold">{{ props.node.count }} targets covering {{props.node.range_sum}} bases.</span>
                </div>
                    
            """,
            )
            orig_tree.add_slot('default-header', r'''
                    <div class="text-weight-bold text-primary">{{ props.node.id }}</div>
                                    
            ''')

            ui.label("Control Target Information")
            cont_tree = ui.tree(
                original_bed_tree
                , label_key="id",
            )
            cont_tree.add_slot(
                "default-body",
                """
                <div v-if="props.node.description">
                        <span class="text-weight-bold">{{ props.node.description }}.</span>
                    </div>
                <div v-if="props.node.count">
                        <span class="text-weight-bold">{{ props.node.count }} targets covering {{props.node.range_sum}} bases.</span>
                </div>

            """,
            )
            cont_tree.add_slot('default-header', r'''
                                <div class="text-weight-bold text-primary">{{ props.node.id }}</div>

                        ''')

        # my_connection.connect_to_minknow()
        with ui.card().classes("w-full"):

            ui.label("Target Information")

            trees = []
            for key in natsorted(CNVResults):
                tree_dict = {}
                # print(key)
                tree_dict["id"] = key
                tree_dict["description"] = f"Targets On Chromosome {key}"
                tree_dict["children"] = [
                    {"id": "+", "description": "Forward Strand", "children": []},
                    {"id": "-", "description": "Reverse Strand", "children": []},
                ]
                # print(CNVResults[key])
                for line in CNVResults[key]["bed_data"].split("\n"):
                    # print(line)
                    # print(tree_dict['children'])
                    # Split each line into its components
                    chrom, start, end, _, _, strand = line.split("\t")
                    if strand == "+":
                        add_child(
                            tree_dict["children"],
                            strand,
                            {"id": f"{chrom}:{start}-{end}"},
                        )
                    else:
                        add_child(
                            tree_dict["children"],
                            strand,
                            {"id": f"{chrom}:{start}-{end}"},
                        )

                # print(dir(tree))
                # print(tree.id)
                trees.append(tree_dict)
            with ui.card().classes("w-full border-[1px]"):
                tree = ui.tree(
                    [
                        {
                            "id": "chromosome",
                            "description": "Targets On Each Chromosome",
                            "children": trees,
                        }
                    ],
                    label_key="id",
                    on_select=lambda e: ui.notify(e.value),
                )

                tree.add_slot(
                    "default-header",
                    """
                    <span :props="props"><strong>{{ props.node.id }}</strong></span>
                """,
                )
                tree.add_slot(
                    "default-body",
                    """
                    <span :props="props">Description: "{{ props.node.description }}"</span>
                """,
                )
                # ui.label(f"{CNVResults}")

                with ui.row():
                    ui.button("+ all", on_click=tree.expand)
                    ui.button("- all", on_click=tree.collapse)
                    ui.button("open + strand", on_click=lambda: tree.expand(["+"]))
                    ui.button("close + strand", on_click=lambda: tree.collapse(["+"]))

        new_bed_tree.load_from_string("chr2\t1\t1000000\t.\t.\t+", merge=True)
        new_bed_tree.load_from_string("chr1\t1\t100000000\t.\t.\t+", merge=True)
        # my_object = MinknowHistograms(my_connection.positions[0])


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
    app.add_static_files("/fonts", str(Path(__file__).parent / "fonts"))
    # app.on_startup(mainpage.index_page)
    # index_page()
    ui.run(
        port=port,
        reload=reload,
        title="MethClass NiceGUI",
        storage_secret="slartibartfast",
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
