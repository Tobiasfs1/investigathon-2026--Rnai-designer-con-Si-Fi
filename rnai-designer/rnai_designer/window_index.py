"""Read-only indexes for exact interval aggregates; floats retain source sum order."""
from bisect import bisect_left, bisect_right
from itertools import accumulate
from bisect import insort
from collections import Counter


class WindowIndex:
    def __init__(self, sirnas, records, duplicate_extras, metrics, gc_percent):
        self.sirnas = sirnas
        self.starts = [s['start_1based'] for s in sirnas]
        self.prefixes = {}
        flags = {}
        for m in metrics:
            flags[m + '_true'] = [s[m] is True for s in sirnas]
            flags[m + '_false'] = [s[m] is False for s in sirnas]
            flags[m + '_known'] = [s[m] is not None and not s[m+'_conflict'] for s in sirnas]
        flags['known_cells'] = [sum(s[m] is not None and not s[m+'_conflict'] for m in metrics) for s in sirnas]
        flags['complete'] = [all(s[m] is not None and not s[m+'_conflict'] for m in metrics) for s in sirnas]
        flags['conflicts'] = [any(s[m+'_conflict'] for m in metrics) for s in sirnas]
        flags['matching'] = [s['reference_sequence_match'] for s in sirnas]
        gc = [gc_percent(s['sequence']) for s in sirnas]
        flags['gc_valid'] = [v is not None for v in gc]
        flags['gc_target'] = [v is not None and 30 <= v <= 52 for v in gc]
        flags['position_change'] = [i == 0 or self.starts[i] != self.starts[i-1] for i in range(len(sirnas))]
        flags['covered_step'] = [0 if i == 0 else min(len(sirnas[i]['sequence']), self.starts[i]-self.starts[i-1]) for i in range(len(sirnas))]
        for key, values in flags.items():
            self.prefixes[key] = list(accumulate(values, initial=0))
        self.raw_starts = [r['sirna_position'] for r in records]
        positions = sorted(set(self.raw_starts))
        self.duplicate_positions = positions
        self.duplicates = list(accumulate((duplicate_extras.get(p, 0) for p in positions), initial=0))
        self.access = [s['accessibility_value'] if not s['accessibility_value_conflict'] else None for s in sirnas]
        self.valid_before = list(accumulate((v is not None for v in self.access), initial=0))
        self.valid_access = [v for v in self.access if v is not None]
        self.sequences = [s['sequence'] for s in sirnas]
        self.previous = None
        self.ordered = []
        self.sequence_counts = Counter()

    def range_data(self, first, last):
        previous = self.previous
        if (previous is None or first < previous[0] or last < previous[1] or
                first > previous[1] or (first-previous[0]+last-previous[1]) > last-first):
            access = self.valid_access[self.valid_before[first]:self.valid_before[last]]
            self.ordered = sorted(access)
            self.sequence_counts = Counter(self.sequences[first:last])
        else:
            for i in range(previous[0], first):
                seq = self.sequences[i]
                self.sequence_counts[seq] -= 1
                if self.sequence_counts[seq] == 0:
                    del self.sequence_counts[seq]
                value = self.access[i]
                if value is not None:
                    self.ordered.pop(bisect_left(self.ordered, value))
            for i in range(previous[1], last):
                self.sequence_counts[self.sequences[i]] += 1
                value = self.access[i]
                if value is not None:
                    insort(self.ordered, value)
            access = self.valid_access[self.valid_before[first]:self.valid_before[last]]
        self.previous = first, last
        return access, self.ordered, len(self.sequence_counts)

    def set_reference(self, sequence):
        self.gc_bases = list(accumulate((b in 'GC' for b in sequence.upper()), initial=0))
        self.canonical_bases = list(accumulate((b in 'ACGT' for b in sequence.upper()), initial=0))

    def gc(self, start, end):
        canonical = self.canonical_bases[end] - self.canonical_bases[start-1]
        return 100.0 * (self.gc_bases[end]-self.gc_bases[start-1])/canonical if canonical else None

    def count(self, key, first, last):
        values = self.prefixes[key]
        return values[last] - values[first]

    def duplicate_count(self, start, max_start):
        if max_start < start:
            return 0
        return (self.duplicates[bisect_right(self.duplicate_positions, max_start)] -
                self.duplicates[bisect_left(self.duplicate_positions, start)])

    def coverage(self, first, last, length, size):
        if first == last:
            return 0, 0.0, 0, 0.0
        represented = 1 + self.count('position_change', first+1, last)
        covered = size + self.count('covered_step', first+1, last)
        possible = max(0, length-size+1)
        return represented, represented/possible if possible else 0.0, covered, covered/length


class OfftargetIndex:
    """Overlap counts via sorted starts/ends, without materializing alignment copies."""
    def __init__(self, evaluations, query_length):
        from .window_explorer import _offtarget_counts
        self.statuses = _offtarget_counts([], evaluations)
        self.intervals = {}
        self.total_starts, self.total_ends = [], []
        for e in evaluations:
            for a in e.get('criterion', {}).get('alignments', []):
                organism = a.get('organism', e['organism'])
                species = ({'Homo sapiens': 'human', 'Zea mays': 'maize', 'Apis mellifera': 'bee'}.get(organism)
                           if isinstance(organism, str) else None)
                start, end = a.get('query_start_1based'), a.get('query_end_1based_inclusive')
                if end is None and type(start) is int and type(a.get('length')) is int:
                    end = start + a['length'] - 1
                if type(start) is not int or type(end) is not int or start < 1 or end < start or end > query_length:
                    continue
                self.total_starts.append(start)
                self.total_ends.append(end)
                if species:
                    suffixes = []
                    if a.get('contains_exact_21nt') is True:
                        suffixes.append('_exact_21nt_overlapping_rows')
                    if a.get('match_type') == 'mismatch_alignment':
                        suffixes.append('_mismatch_overlapping_rows')
                    if a.get('match_type') == 'unverified':
                        suffixes.append('_unverified_overlapping_rows')
                    for suffix in suffixes:
                        starts, ends = self.intervals.setdefault(species+suffix, ([], []))
                        starts.append(start)
                        ends.append(end)
        self.total_starts.sort()
        self.total_ends.sort()
        for starts, ends in self.intervals.values():
            starts.sort()
            ends.sort()

    def counts(self, start, end):
        result = self.statuses.copy()
        for key, (starts, ends) in self.intervals.items():
            result[key] = bisect_right(starts, end) - bisect_left(ends, start)
        return bisect_right(self.total_starts, end)-bisect_left(self.total_ends, start), result
