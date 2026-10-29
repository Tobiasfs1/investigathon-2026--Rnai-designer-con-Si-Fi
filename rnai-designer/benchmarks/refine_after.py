"""Measure final column refinement; compare exact bytes with baseline-verified CSVs."""
import json
import subprocess
import sys
import time
from pathlib import Path
from suite import latest_export

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT/'benchmarks'/'results'


if __name__ == '__main__':
    while not (RESULTS/'instrumented_complete.json').is_file():
        time.sleep(1)
    records = json.loads((RESULTS/'comparison.json').read_text())
    for record in records:
        n = record['nt']
        previous = latest_export('after', n)
        previous_winner = json.loads((previous/'selected_window.json').read_text(encoding='utf-8'))
        subprocess.run([sys.executable, 'benchmarks/benchmark.py', '--mode', 'after', '--lengths', str(n)], cwd=ROOT, check=True)
        current = latest_export('after', n)
        with (previous/'all_windows.csv').open('rb') as a, (current/'all_windows.csv').open('rb') as b:
            while True:
                left, right = a.read(1024*1024), b.read(1024*1024)
                assert left == right, (n, 'CSV bytes changed')
                if not left:
                    break
        winner = json.loads((current/'selected_window.json').read_text(encoding='utf-8'))
        previous_winner.pop('input_report')
        winner.pop('input_report')
        assert previous_winner == winner
        record.setdefault('previous_optimized_measurements', []).append(record['after'])
        record['after'] = json.loads((RESULTS/'after'/str(n)/'measurement.json').read_text())
        assert record['after']['windows'] == record['verified_identical_rows']
        record['final_verification'] = {'method': 'Byte-for-byte CSV equality with previously baseline-verified export, plus exact winner equality',
                                        'previous_export': str(previous), 'final_export': str(current)}
        record['speedup'] = record['before']['total_s']/record['after']['total_s']
        record['peak_rss_reduction_percent'] = 100*(1-record['after']['peak_rss_bytes']/record['before']['peak_rss_bytes'])
        (RESULTS/'comparison.json').write_text(json.dumps(records, indent=2))
        print(f'Final column implementation verified: {n} nt', flush=True)
    subprocess.run([sys.executable, 'benchmarks/benchmark.py', '--mode', 'after', '--lengths', '500', '--profile'], cwd=ROOT, check=True)
    (RESULTS/'final_after_complete.json').write_text(json.dumps({'completed': True}))
