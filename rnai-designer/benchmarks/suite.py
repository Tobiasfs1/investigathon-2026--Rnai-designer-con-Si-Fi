"""Final comparison: one worker at a time, fresh process per size and version."""
import csv
import hashlib
import json
import os
import platform
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT/'benchmarks'/'results'


def latest_export(mode, n):
    # Directory mtime changes when obsolete large CSVs are removed. The export
    # identifier encodes creation order; incomplete/cleaned exports are ineligible.
    exports = sorted((RESULTS/mode/str(n)).glob('export_*'),
                     key=lambda p: int(p.name.removeprefix('export_')), reverse=True)
    for directory in exports:
        try:
            if (directory/'all_windows.csv').is_file() and (directory/'selected_window.json').is_file():
                return directory
        except OSError:
            continue
    raise FileNotFoundError(f'No complete export for {mode}/{n}')


def compare_tables(n):
    before, after = (latest_export(mode, n) for mode in ('before', 'after'))
    csv.field_size_limit(2**30)
    count = 0
    with (before/'all_windows.csv').open(encoding='utf-8', newline='') as a, (after/'all_windows.csv').open(encoding='utf-8', newline='') as b:
        for old, new in zip(csv.DictReader(a), csv.DictReader(b), strict=True):
            old.pop('offtarget_matches')
            old.pop('offtarget_by_species')
            for key in ('quality_components', 'efficacy_priority_components'):
                for k, value in json.loads(old.pop(key)).items():
                    old[key+'.'+k] = '' if value is None else str(value)
            assert old == new, (n, count, [(k, v, new.get(k)) for k, v in old.items() if new.get(k) != v])
            count += 1
            if count % 50000 == 0:
                print(f'Checked {count} rows for {n} nt', flush=True)
    old_winner, new_winner = (json.loads((d/'selected_window.json').read_text(encoding='utf-8')) for d in (before, after))
    # Paths/provenance refer to different benchmark input files, scientific values must match.
    old_winner.pop('input_report')
    new_winner.pop('input_report')
    assert old_winner == new_winner
    return count


if __name__ == '__main__':
    RESULTS.mkdir(parents=True, exist_ok=True)
    environment = {'python': sys.version, 'platform': platform.platform(), 'processor': platform.processor(),
                   'cpu_count': os.cpu_count(), 'seed': 817, 'stride': 1, 'lengths': [180, 250],
                   'note': 'Sequential workers, one run per size/version; no biological databases. Time excludes interpreter/module imports.'}
    environment['reference_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'benchmarks'/'reference').glob('*.py')}
    (RESULTS/'environment.json').write_text(json.dumps(environment, indent=2))
    comparison = []
    for n in (500, 1000, 2500, 5000):
        for mode in ('before', 'after'):
            subprocess.run([sys.executable, str(ROOT/'benchmarks'/'benchmark.py'), '--mode', mode, '--lengths', str(n), '--worker'], cwd=ROOT, check=True)
        old, new = (json.loads((RESULTS/mode/str(n)/'measurement.json').read_text()) for mode in ('before', 'after'))
        verified = compare_tables(n)
        row = {'nt': n, 'verified_identical_rows': verified, 'before': old, 'after': new,
               'speedup': old['total_s']/new['total_s'],
               'peak_rss_reduction_percent': 100*(1-new['peak_rss_bytes']/old['peak_rss_bytes'])}
        comparison.append(row)
        (RESULTS/'comparison.json').write_text(json.dumps(comparison, indent=2))
        print('VERIFIED ' + json.dumps({'nt': n, 'rows': verified, 'speedup': row['speedup']}), flush=True)
