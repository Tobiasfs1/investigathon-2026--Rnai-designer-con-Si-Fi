"""Column storage with the existing iterable-of-dictionaries Python interface."""
from array import array
from collections.abc import Sequence
import math
from .priority_scoring import priority_key


class CompactWindows(Sequence):
    def __init__(self):
        self.columns = {}
        self.kinds = {}
        self.categories = {}
        self.nullable = set()
        self.n = 0
        self.best = None
        self.best_index = None
        self.best_key = None
        self.best_score = None
        self.best_score_tie_count = 0
        self.short_windows = []
        self._schema = []
        self._row_size = None

    def append(self, row):
        if row.get('window_class') == 'short_green_region':
            self.short_windows.append(row)
        score = row.get('adjusted_priority_score', row['efficacy_priority_score'])
        if score is not None:
            key = priority_key(row)
            if self.best_key is None or key < self.best_key:
                self.best_key, self.best = key, row
                self.best_index = self.n
            if self.best_score is None or score > self.best_score:
                self.best_score, self.best_score_tie_count = score, 1
            elif score == self.best_score:
                self.best_score_tie_count += 1
        if not self.columns:
            self._row_size = len(row)
            for parent, value in row.items():
                if parent in ('window_id', 'offtarget_by_species'):
                    continue
                children = value.keys() if isinstance(value, dict) else (None,)
                for child in children:
                    key = parent+'.'+child if child is not None else parent
                    column = self.columns[key] = array('d')
                    self.kinds[key] = None
                    self._schema.append((key, parent, child, column))
        if len(row) != self._row_size:
            raise ValueError('Inconsistent window schema')
        kinds, categories, nullable = self.kinds, self.categories, self.nullable
        for key, parent, child, column in self._schema:
            value = row[parent] if child is None else row[parent][child]
            if value is None:
                encoded = math.nan
                nullable.add(key)
            elif isinstance(value, str):
                entry = categories.get(key)
                if entry is None:
                    entry = categories[key] = ({}, [])
                    kinds[key] = str
                codes, values = entry
                encoded = codes.get(value)
                if encoded is None:
                    encoded = codes[value] = len(values)
                    values.append(value)
            else:
                kind = kinds[key]
                if kind is None:
                    kinds[key] = type(value)
                elif kind is int and type(value) is float:
                    # A later integer must never truncate an earlier fractional value.
                    kinds[key] = float
                encoded = value
            column.append(encoded)
        self.n += 1

    def ranked_indices(self):
        """Try the known best first; sort scalar columns only if it fails."""
        if self.best_index is None:
            return
        yield self.best_index
        columns = self.columns
        score = columns['adjusted_priority_score']
        component_names = ('sifi_is_efficient_fraction', 'sirna_gc_30_52_fraction', 'unique_sirna_sequence_fraction')
        components = [columns['efficacy_priority_components.'+name] for name in component_names]
        regions = self.categories['region_id'][1]
        def key(i):
            return (-score[i], *(-column[i] for column in components), -columns['quality_score'][i],
                    columns['start_1based'][i], columns['length_nt'][i], regions[int(columns['region_id'][i])])
        yield from sorted((i for i in range(self.n) if i != self.best_index and not math.isnan(score[i])), key=key)

    def eligible_count(self):
        return sum(not math.isnan(value) for value in self.columns.get('adjusted_priority_score', ()))

    def __len__(self):
        return self.n

    def __getitem__(self, index):
        if isinstance(index, slice):
            return [self[i] for i in range(*index.indices(self.n))]
        if index < 0:
            index += self.n
        if not 0 <= index < self.n:
            raise IndexError(index)
        row = {}
        for key, column in self.columns.items():
            value = column[index]
            kind = self.kinds[key]
            value = (None if math.isnan(value) else self.categories[key][1][int(value)]
                     if kind is str else kind(value))
            if '.' in key:
                parent, child = key.split('.', 1)
                row.setdefault(parent, {})[child] = value
            else:
                row[key] = value
        row['window_id'] = f"{row['region_id']}:{row['start_1based']}-{row['end_1based_inclusive']}"
        row['offtarget_by_species'] = {key: value for key, value in row.items()
                                     if key.startswith(('human_', 'maize_', 'bee_'))}
        return row

    def __iter__(self):
        for i in range(self.n):
            yield self[i]

    def scalar_rows(self):
        """Flatten nested scalar components; never serialize objects per window."""
        fields = list(self.columns)
        iterators = []
        for key in fields:
            column, kind = self.columns[key], self.kinds[key]
            if key in self.nullable:
                categories = self.categories[key][1] if kind is str else None
                def decode(value, kind=kind, categories=categories):
                    return (None if math.isnan(value) else categories[int(value)]
                            if kind is str else kind(value))
                values = map(decode, column)
            elif kind is str:
                values = map(self.categories[key][1].__getitem__, map(int, column))
            elif kind is float:
                values = iter(column)
            else:
                values = map(kind, column)
            iterators.append(values)
        region, start, end = (fields.index(key) for key in ('region_id', 'start_1based', 'end_1based_inclusive'))
        for values in zip(*iterators):
            yield (f'{values[region]}:{values[start]}-{values[end]}', *values)
