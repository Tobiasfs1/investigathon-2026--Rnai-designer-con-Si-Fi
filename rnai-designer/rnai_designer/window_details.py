"""Recover source detail on demand, without embedding it in each table row."""
from bisect import bisect_left, bisect_right
from .window_explorer import explore_windows, _sirnas, _overlapping_offtargets


def resolve_window_details(report, window_id, sirna_size=None):
    region_id, coordinates = window_id.rsplit(':', 1)
    start, end = map(int, coordinates.split('-'))
    if start < 1 or end < start:
        raise ValueError('Coordenadas de ventana inválidas')
    analysis = explore_windows(report, end-start+1, sirna_size=sirna_size, validation_only=True)
    region = next((r for r in analysis['green_regions'] if r['region_id'] == region_id), None)
    if region is None or start < region['start_1based'] or end > region['end_1based_inclusive']:
        raise ValueError('La ventana no pertenece a la región verde indicada')
    sirnas, _ = _sirnas(report['sifi']['records'])
    starts = [s['start_1based'] for s in sirnas]
    last_start = end-analysis['sirna_size']+1
    contained = (sirnas[bisect_left(starts, start):bisect_right(starts, last_start)]
                 if last_start >= start else [])
    sequence = report['candidate']['reference']['sequence']
    return {'window_id': window_id, 'region_id': region_id,
            'start_1based': start, 'end_1based_inclusive': end, 'length_nt': end-start+1,
            'sequence': sequence[start-1:end], 'sirnas': contained,
            'offtarget_matches': _overlapping_offtargets(report.get('offtarget_evaluations', []), start, end, len(sequence))}
