import numpy as np
import pandas as pd
import os
from collections import defaultdict


class CNVChangeDetectorTracker:
    """
    This class is designed to store breakpoints detected from copy
    number variation analysis in real time. It will identify candidate
    start and end points for targets that are likely to include structural
    variants of importance.
    The class will automatically limit the range of targets detected to a given
    proportion of the target genome.
    """
    def __init__(self, proportion = 0.01):
        """
        The proportion value will be used to scale the targets provided back to the
        user.
        """
        self.proportion = proportion
        self.coordinates = {}

    def add_breakpoints(self, chromosome, breakpoints, length):
        increments = defaultdict(int)
        for entry in breakpoints:
            increments[entry['start']] += 1
            increments[entry['end']] -= 1

        current_value = 0
        if chromosome not in self.coordinates:
            self.coordinates[chromosome] = {}
            self.coordinates[chromosome]["positions"] = []
            self.coordinates[chromosome]["length"] = length
            self.coordinates[chromosome]['target'] = int(length * self.proportion)

        sorted_positions = sorted(increments.keys())
        for position in sorted_positions:
            current_value += increments[position]
            self.coordinates[chromosome]["positions"].append((position, current_value))
        self._calculate_threshold(chromosome)

    def _init_positions(self, chromosome):
        self.coordinates[chromosome]["current_value"] = 0
        self.coordinates[chromosome]["flush"] = True
        self.coordinates[chromosome]["current_value"] = 0
        self.coordinates[chromosome]["previous_value"] = 0
        self.coordinates[chromosome]["start_value"] = 0
        self.coordinates[chromosome]["end_value"] = 0
        self.coordinates[chromosome]["cumulative_sum"] = 0
        self.coordinates[chromosome]["start_positions"] = []
        self.coordinates[chromosome]["end_positions"] = []

    def _calculate_threshold(self, chromosome):
        self.coordinates[chromosome]["threshold"] = 0
        length = self.coordinates[chromosome]["length"]

        while length > self.coordinates[chromosome]["target"]:
            self._init_positions(chromosome)
            for (pos, depth) in self.coordinates[chromosome]["positions"]:
                if depth > self.coordinates[chromosome]["threshold"]:
                    self.update(chromosome, 1, pos)
                else:
                    self.update(chromosome,0, pos)
            length = self.get_cumulative_sum(chromosome)
            self.coordinates[chromosome]["threshold"] += 1

    def get_threshold(self,chromosome):
        return self.coordinates[chromosome]["threshold"]

    def get_cumulative_sum(self, chromosome):
        if self.coordinates[chromosome]["flush"]:
            return self.coordinates[chromosome]["cumulative_sum"]
        else:
            return self.coordinates[chromosome]["cumulative_sum"] + (self.coordinates[chromosome]["end_value"] - self.coordinates[chromosome]["start_value"] + 1)

    def update(self, chromosome, new_value, position):
        if new_value != self.coordinates[chromosome]["current_value"]:
            if new_value == 1 and self.coordinates[chromosome]["current_value"] == 0:
                #print("Change detected: 0 -> 1")
                self.coordinates[chromosome]["start_value"] = position
                self.coordinates[chromosome]["end_value"] = position
                self.coordinates[chromosome]["start_positions"].append(position)
                self.coordinates[chromosome]["flush"] = False
            elif new_value == 0 and self.coordinates[chromosome]["current_value"] == 1:
                #print("Change detected: 1 -> 0")
                self.coordinates[chromosome]["end_value"] = position
                self.coordinates[chromosome]["cumulative_sum"] += (self.coordinates[chromosome]["end_value"] - self.coordinates[chromosome]["start_value"] + 1)
                self.coordinates[chromosome]["end_positions"].append(position)
                self.coordinates[chromosome]["flush"] = False
            self.coordinates[chromosome]["previous_value"] = self.coordinates[chromosome]["current_value"]
            self.coordinates[chromosome]["current_value"] = new_value
        else:
            #print("No change detected")
            if new_value == 1:
                self.coordinates[chromosome]["end_value"] = position

    def get_breakpoints(self, chromosome):
        return (self.coordinates[chromosome]["start_positions"], self.coordinates[chromosome]["end_positions"])

    def get_bed_targets(self, chromosome):
        if chromosome in self.coordinates.keys():
            bedlist = [(chromosome, first, second) for first, second in zip(self.coordinates[chromosome]["start_positions"],self.coordinates[chromosome]["end_positions"])]
            # Convert the list to BED format
            bed_lines = []
            for idx, (chrom, start, end) in enumerate(bedlist):
                name = "." #f"region_{idx + 1}"  # Create a name for each region
                score = "." #0  # You can adjust the score as needed
                bed_lines.append(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t+")
                bed_lines.append(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t-")

            return "\n".join(bed_lines)
        #return (list(zip(*(self.coordinates[chromosome]["start_positions"],self.coordinates[chromosome]["end_positions"]))))

class ChangeDetector:
    def __init__(self, breakpoints, proportion = 0.01, length = 1_000_000_000):
        self.length = length
        self.target = int(self.length * proportion)
        self._init_positions()
        self.flush = True
        self.threshold = 0
        self.increments = defaultdict(int)
        for entry in breakpoints:
            self.increments[entry['start']] += 1
            self.increments[entry['end'] + 1] -= 1
        self.coordinates = []
        self.current_value = 0
        self.sorted_positions = sorted(self.increments.keys())
        for position in self.sorted_positions:
            self.current_value += self.increments[position]
            self.coordinates.append((position, self.current_value))

    def update(self, new_value, position):
        if new_value != self.current_value:
            if new_value == 1 and self.current_value == 0:
                #print("Change detected: 0 -> 1")
                self.start_value = position
                self.end_value = position
                self.start_positions.append(position)
                self.flush = False
            elif new_value == 0 and self.current_value == 1:
                #print("Change detected: 1 -> 0")
                self.end_value = position
                self.cumulative_sum += (self.end_value - self.start_value + 1)
                self.end_positions.append(position)
                self.flush = False
            self.previous_value = self.current_value
            self.current_value = new_value
        else:
            #print("No change detected")
            if new_value == 1:
                self.end_value = position

    def _init_positions(self):
        self.flush = True
        self.current_value = 0
        self.previous_value = 0
        self.start_value = 0
        self.end_value = 0
        self.cumulative_sum = 0
        self.start_positions = []
        self.end_positions = []

    def calculate_threshold(self):
        length = self.length

        while length > self.target:
            self._init_positions()
            for (pos, depth) in self.coordinates:
                if depth > self.threshold:
                    self.update(1, pos)
                else:
                    self.update(0, pos)
            length = self.get_cumulative_sum()
            self.threshold += 1
        return self.threshold

    def get_cumulative_sum(self):
        if self.flush:
            return self.cumulative_sum
        else:
            return self.cumulative_sum + (self.end_value - self.start_value + 1)

    def get_breakpoints(self):
        return (self.start_positions, self.end_positions)

    def get_bed_targets(self):
        return (list(zip(*(self.start_positions,self.end_positions))))



if __name__ == "__main__":
    output = "/Users/mattloose/GIT/niceGUI/cnsmeth/robin_output_test/_ds1305_Intraop0051_c_ULK_2_highFRA"
    DATA_ARRAY = np.load(os.path.join(output, "ruptures.npy"), allow_pickle="TRUE").item()
    #print(type(DATA_ARRAY))
    #print(DATA_ARRAY)
    #print(DATA_ARRAY["chr9"])
    #print(DATA_ARRAY["chr9"].keys())
    #for item in list(range(1, 23)) + ['X', 'Y']:
    #    contig = f"chr{item}"
    #    breakpoints = DATA_ARRAY[contig]
    #    detector = ChangeDetector(breakpoints)
    #    print(f"Chromosome {contig}")
    #    print(f"Threshold now is {detector.calculate_threshold()}")
    #    print(f"Breakpoints are {detector.get_breakpoints()}")
    #    print(f"Bed targets are {detector.get_bed_targets()}")

    """
    detector2 = CNVChangeDetectorTracker()
    for item in list(range(1, 23)) + ['X', 'Y']:
        contig = f"chr{item}"
        breakpoints = DATA_ARRAY[contig]
        detector2.add_breakpoints(contig,breakpoints,1_000_000_000)


    for item in list(range(1, 23)) + ['X', 'Y']:
        contig = f"chr{item}"
        print (f"contig - {contig}, {detector2.get_threshold(contig)}, {detector2.get_cumulative_sum(contig)}")
        print(detector2.get_bed_targets(contig))
    """
    detector2 = CNVChangeDetectorTracker()
    detector2.coordinates = DATA_ARRAY
    for item in list(range(1, 23)) + ['X', 'Y']:
        contig = f"chr{item}"
        result = detector2.get_bed_targets(contig)
        if result:
            print(result)
