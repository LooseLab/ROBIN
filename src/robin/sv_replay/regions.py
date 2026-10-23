"""Indexed genomic context; masks annotate retained evidence, never delete it."""
from bisect import bisect_right
from pathlib import Path

RESOURCE_DIR = Path(__file__).resolve().parent.parent / 'resources'


class IntervalIndex:
    def __init__(self, path):
        intervals = {}
        with open(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith(('#', 'track', 'browser')):
                    continue
                chrom, start, end, *_ = line.split()
                start, end = int(start), int(end)
                if start < 0 or end <= start:
                    raise ValueError(f'Invalid BED interval: {line.strip()}')
                intervals.setdefault(chrom, []).append((start, end))
        self.intervals = intervals
        self.index = {}
        for chrom, rows in intervals.items():
            rows.sort()
            merged = []
            for start, end in rows:
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
                else:
                    merged.append((start, end))
            rows = merged
            starts, ends, maximum = [], [], -1
            for start, end in rows:
                maximum = max(maximum, end)
                starts.append(start)
                ends.append(maximum)
            self.index[chrom] = starts, ends

    def overlap_context(self, chrom, pos):
        """Containing repeat union and distance to its nearest edge (BED coordinates)."""
        starts, ends = self.index.get(chrom, ([], []))
        i = bisect_right(starts, pos) - 1
        if i < 0 or pos >= ends[i]:
            return None
        end = ends[i]
        while i > 0 and starts[i] <= ends[i-1]:
            i -= 1
        start = starts[i]
        return dict(start=start, end=end, edge_distance=min(pos-start, end-1-pos))

    def contains(self, chrom, pos):
        starts, ends = self.index.get(chrom, ([], []))
        i = bisect_right(starts, pos) - 1
        return i >= 0 and pos < ends[i]


def call_status(support, flags, min_support):
    if 'mitochondrial_reference_boundary' in flags:
        return 'reference_boundary'
    if support < min_support:
        return 'candidate'
    return 'masked_candidate' if flags else 'provisional'


class RegionContext:
    def __init__(self, lengths, tandem_repeats, centromeres, extra_beds=(),
                 end_margin=1_000_000, origin_margin=100):
        self.lengths = lengths
        self.indexes = [('tandem_repeat', IntervalIndex(tandem_repeats)),
                        ('centromere_satellite', IntervalIndex(centromeres))]
        self.indexes.extend(('custom_mask', IntervalIndex(p)) for p in extra_beds)
        self.end_margin = end_margin
        self.origin_margin = origin_margin

    def flags(self, sig):
        a, b = sig['a'], sig['b']
        mitochondrial = {'chrM', 'MT', 'M'}
        if (sig['kind'] == 'BND' and a[0] == b[0] and a[0] in mitochondrial
                and a[2] == '-' and b[2] == '+'
                and a[0] in self.lengths
                and 0 <= a[1] <= self.origin_margin
                and 0 <= self.lengths[a[0]] - b[1] <= self.origin_margin):
            return {'mitochondrial_reference_boundary'}
        flags = set()
        for chrom, pos, orientation in (a, b):
            # Annotate the retained flank base, not the entire deleted interval.
            base = max(0, pos - 1 if orientation == '+' else pos)
            for name, index in self.indexes:
                if index.contains(chrom, base):
                    flags.add(name)
            if chrom not in mitochondrial and self.end_margin:
                length = self.lengths.get(chrom)
                if base < self.end_margin or (length is not None and
                                               base >= length - self.end_margin):
                    flags.add('chromosome_end_margin')
        return flags
