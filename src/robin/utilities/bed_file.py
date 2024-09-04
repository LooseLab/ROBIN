"""
File included for development purposes.

Needs to be refactored.

"""

# Python imports.
from __future__ import annotations
from nicegui import ui, app
from pathlib import Path

from robin import theme
import numpy as np
import os

from nicegui import ui
import csv
from natsort import natsorted
from io import StringIO
import copy
import hashlib
import logging
from datetime import datetime
import pandas as pd


class BedTree:
    def __init__(self, preserve_original_tree=False, reference_file=None, output_location=None):
        self.total_count = 0
        self.total_range_sum = 0
        self.tree_data = {'id': "Chromosomes", 'description': 'Summary of currently used targets.', 'children': [], 'count': 0, 'range_sum': 0, 'proportion': 0}
        self.tree_dict = {}
        self.preserve_original_tree = preserve_original_tree
        self.reference_tree = None
        self.file_counter = 0  # Counter for the file naming
        self.previous_targets_hash = None  # To track changes in targets
        self.proportions_df = pd.DataFrame()  # DataFrame to store the proportions and timestamps
        self.output_location = output_location
        if reference_file:
            self.chromosome_lengths = self._get_chromosome_lengths(reference_file)
            self.total_length = sum(self.chromosome_lengths.values())
        else:
            self.chromosome_lengths = None
            self.total_length = None

    def _hash_current_targets(self):
        """Generates a hash for the current targets to detect changes."""
        m = hashlib.md5()
        for chromosome, chromosome_data in sorted(self.tree_dict.items()):
            m.update(chromosome.encode())
            for strand_group in sorted(chromosome_data['children'], key=lambda x: x['id']):
                m.update(strand_group['id'].encode())
                for target in sorted(strand_group['children'], key=lambda x: x['id']):
                    m.update(target['id'].encode())
        return m.hexdigest()

    def _check_and_create_folder(self, path, folder_name=None):
        # Check if the path exists
        if not os.path.exists(path):
            raise FileNotFoundError(f"The specified path does not exist: {path}")

        # If folder_name is provided
        if folder_name:
            full_path = os.path.join(path, folder_name)
            # Create the folder if it doesn't exist
            if not os.path.exists(full_path):
                os.makedirs(full_path)
                #logger.info(f"Folder created: {full_path}")
            return full_path
        else:
            return path

    def _write_bed_file(self):
        """Writes the current targets to a new BED file if they have changed."""
        current_hash = self._hash_current_targets()
        if current_hash != self.previous_targets_hash:
            self.file_counter += 1
            if self.output_location:
                filename = os.path.join(self._check_and_create_folder(self.output_location, folder_name="bed_files"), f"new_file_{self.file_counter}.bed")
            else:
                filename = f"new_file_{self.file_counter}.bed"
            with open(filename, 'w') as f:
                for chromosome, chromosome_data in sorted(self.tree_dict.items()):
                    for strand_group in sorted(chromosome_data['children'], key=lambda x: x['id']):
                        for target in sorted(strand_group['children'], key=lambda x: x['id']):
                            start, end = map(int, target['id'].split('-'))
                            f.write(f"{chromosome}\t{start}\t{end}\t.\t.\t{strand_group['id'].split()[-1]}\n")
            self.previous_targets_hash = current_hash
            self._update_proportions_df()
            #print(f"Wrote new BED file: {filename}")

    def _get_chromosome_lengths(self, fai_file):
        chromosome_lengths = {}
        with open(fai_file, 'r') as f:
            for line in f:
                parts = line.split('\t')
                chromosome = parts[0]
                length = int(parts[1])
                chromosome_lengths[chromosome] = length
        return chromosome_lengths

    def _merge_ranges(self, existing_ranges, new_range):
        start, end = new_range
        merged_ranges = []
        i = 0

        while i < len(existing_ranges):
            existing_start, existing_end = existing_ranges[i]
            if not (end < existing_start or start > existing_end):
                start = min(existing_start, start)
                end = max(existing_end, end)
            else:
                merged_ranges.append((existing_start, existing_end))
            i += 1

        merged_ranges.append((start, end))
        merged_ranges.sort()

        return merged_ranges

    def _process_bed_data(self, reader):
        for row in reader:
            if len(row)>0:
                chromosome = row[0]
                start = int(row[1])
                end = int(row[2])
                strand = row[5] if len(row) > 5 else '.'

                if chromosome not in self.tree_dict:
                    self.tree_dict[chromosome] = {'id': chromosome, 'children': [], 'count': 0, 'range_sum': 0, 'description': "chromosome"}

                strand_id = f'{chromosome} {strand}'
                strand_group = next((child for child in self.tree_dict[chromosome]['children'] if child['id'] == strand_id), None)
                if not strand_group:
                    strand_group = {'id': strand_id, 'children': [], 'count': 0, 'range_sum': 0}
                    self.tree_dict[chromosome]['children'].append(strand_group)

                existing_ranges = [(int(r['id'].split('-')[0]), int(r['id'].split('-')[1])) for r in strand_group['children']]
                merged_ranges = self._merge_ranges(existing_ranges, (start, end))

                strand_group['children'] = [{'id': f'{s}-{e}', 'description': 'target region', "chromosome_length": self.chromosome_lengths.get(chromosome, 0), 'range_sum': e - s + 1, "proportion": (e - s + 1 )/self.chromosome_lengths.get(chromosome, 0) * 100  } for s, e in merged_ranges]
                #strand_group['count'] = len(merged_ranges)
                #strand_group['range_sum'] = sum(e - s + 1 for s, e in merged_ranges)
                #strand_group['proportion'] = sum(e - s + 1 for s, e in merged_ranges)/self.chromosome_lengths.get(chromosome, 0) * 100

                self.tree_dict[chromosome]['count'] = sum(c['count'] for c in self.tree_dict[chromosome]['children'])
                self.tree_dict[chromosome]['range_sum'] = sum(c['range_sum'] for c in self.tree_dict[chromosome]['children'])

        self.total_count = sum(v['count'] for v in self.tree_dict.values())
        self.total_range_sum = sum(v['range_sum'] for v in self.tree_dict.values())

    def _add_if_not_exists(self, entries, new_item, description, chromosome):
        if not any(entry['id'] == new_item for entry in entries):
            entries.append({'id': new_item, 'children': [], 'description': description, 'count': 0, 'range_sum': 0, 'proportion': 0, "chromosome_length":self.chromosome_lengths.get(chromosome, 0)})

    def _propagate_values(self, node):
        if 'children' not in node or not node['children']: # Here we are looking at an individual target on a chromosome.
            node['count'] = node.get('count', 1)
            node['range_sum'] = node.get('range_sum', 0)
            #node['proportion'] = node['range_sum'] / self.total_length if self.total_length else 0
            return node['count'], node['range_sum']

        total_count = 0
        total_range_sum = 0

        for child in node['children']: # Here we are not in a tiny node - we are in a bigger node.
            child_count, child_range_sum = self._propagate_values(child)
            total_count += child_count
            total_range_sum += child_range_sum

        node['count'] = total_count
        node['range_sum'] = total_range_sum

        if node['description'] == "Stranded Targets": # We are looking at a stranded chromosome.
            node['proportion'] = node['range_sum']/node['chromosome_length']* 100 if self.total_length else 0
        elif node['description'] == "chromosome":
            node['proportion'] = node['range_sum'] / (2*node['chromosome_length']) * 100 if self.total_length else 0
        elif node['description'] == "Summary of currently used targets.":
            node['proportion'] = node['range_sum'] / (2 * self.total_length) * 100 if self.total_length else 0

        return node['count'], node['range_sum']

    def _build_tree(self):
        for chromosome_data in self.tree_dict.values():
            self._add_if_not_exists(self.tree_data['children'], chromosome_data['id'], "chromosome",chromosome_data['id'])

            for strand_group in natsorted(chromosome_data['children']):
                strand_group['children'] = natsorted(strand_group['children'], key=lambda x: int(x['id'].split('-')[0]))
                self._update_strands(chromosome_data['id'], strand_group['id'])
                self._update_targets(chromosome_data['id'], strand_group['id'], strand_group)

        self.tree_data["children"] = natsorted(self.tree_data["children"], key=lambda x: x['id'])

        self._propagate_values(self.tree_data)

    def _update_strands(self, chromosome, strand):
        chromosome_entry = next(item for item in self.tree_data['children'] if item["id"] == chromosome)
        self._add_if_not_exists(chromosome_entry['children'], strand, "Stranded Targets", chromosome)

    def _update_targets(self, chromosome, strand, strand_group):
        chromosome_entry = next(item for item in self.tree_data['children'] if item["id"] == chromosome)
        strand_entry = next(orientation for orientation in chromosome_entry['children'] if orientation['id'] == strand)
        strand_entry['children'] = strand_group['children']
        strand_entry['count'] = strand_group['count']
        strand_entry['range_sum'] = strand_group['range_sum']
        strand_entry['chromosome_length'] = self.chromosome_lengths.get(chromosome, 0)  # Add chromosome length

    def load_from_file(self, file_path, merge=False):
        if not merge:
            if 'children' in self.tree_data.keys():
                self.tree_data['children'] = []
            self.tree_dict = {}
        with open(file_path, 'r') as bed_file:
            reader = csv.reader(bed_file, delimiter='\t')
            self._process_bed_data(reader)
        self._build_tree()
        self._propagate_values(self.tree_data)
        if self.preserve_original_tree:
            self.reference_tree = copy.deepcopy(self.tree_dict)
        self._propagate_values(self.tree_data)
        #self._write_bed_file()

    def load_from_string(self, bed_string, merge=False, write_files=False, output_location=None):
        if output_location:
            self.output_location = output_location
        if not merge:
            if 'children' in self.tree_data.keys():
                self.tree_data['children'] = []
            self.tree_dict = {}
            if self.preserve_original_tree:
                if self.reference_tree:
                    self.tree_dict = copy.deepcopy(self.reference_tree)

        bed_file = StringIO(bed_string)
        reader = csv.reader(bed_file, delimiter='\t')
        self._process_bed_data(reader)
        self._build_tree()
        self._propagate_values(self.tree_data)
        if write_files:
            self._write_bed_file()


    def _update_proportions_df(self):
        """Updates the proportions DataFrame with the current data and timestamp."""
        timestamp = datetime.now()

        # Top-level (genome-wide) proportion
        top_level_proportion = {
            'timestamp': timestamp,
            'chromosome': 'Genome-wide',
            'proportion': self.tree_data['proportion']
        }

        # List to hold the rows
        rows = [top_level_proportion]

        # Proportions for each chromosome
        for chromosome in self.tree_data['children']:
            chromosome_proportion = {
                'timestamp': timestamp,
                'chromosome': chromosome['id'],
                'proportion': chromosome['proportion']
            }
            rows.append(chromosome_proportion)

        # Create a DataFrame from the rows and append it to the existing DataFrame
        new_df = pd.DataFrame(rows)
        self.proportions_df = pd.concat([self.proportions_df, new_df], ignore_index=True)
        if self.output_location:
            self.proportions_df.to_csv(
                os.path.join(self.output_location,"bedranges.csv")
                )



        #print(self.proportions_df)

    #def _log_proportions(self):
    #    """Logs the top-level proportion and the proportion for each chromosome along with the time."""
    #    # Log the top-level proportion (for the whole genome)
    #    print(f"Total genome proportion covered: {self.tree_data['proportion']:.6f}")#

        # Log proportions for each chromosome
    #    for chromosome in self.tree_data['children']:
    #        print(f"Chromosome {chromosome['id']} proportion covered: {chromosome['proportion']:.6f}")

@ui.page("/", response_timeout=30)
def index_page() -> None:
    # from check_connection import ConnectionDialog
    # initial_ip = "127.0.0.1"
    # my_connection = ConnectionDialog(initial_ip)
    # my_connection = None
    existing_bed = "/Users/mattloose/GIT/niceGUI/cnsmeth/panel_adaptive_nogenenames_20122021_hg38_VGL_CHD7_5kb_pad.bed"

    original_bed_tree = BedTree(preserve_original_tree=True, reference_file="/Users/mattloose/references/hg38_simple.fa.fai")

    original_bed_tree.load_from_file(existing_bed)

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
                [original_bed_tree.tree_data]
               ,label_key="id",
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
        native=True, window_size=(1200, 800), fullscreen=False
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
