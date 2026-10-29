"""Complete instrumented probes after the primary sequential suite has finished."""
import json
import subprocess
import sys
import time
import os
import hashlib
from uuid import uuid4
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


if __name__ == '__main__':
    comparison = ROOT/'benchmarks'/'results'/'comparison.json'
    while True:
        try:
            if len(json.loads(comparison.read_text())) == 4:
                break
        except (OSError, ValueError):
            pass
        time.sleep(1)
    # Re-run every optimized case with the final source and the same fresh-worker
    # protocol, after the entire reference run has released its memory.
    from suite import compare_tables
    records = json.loads(comparison.read_text())
    for record in records:
        n = record['nt']
        subprocess.run([sys.executable, 'benchmarks/benchmark.py', '--mode', 'after', '--lengths', str(n)], cwd=ROOT, check=True)
        record['after'] = json.loads((ROOT/'benchmarks'/'results'/'after'/str(n)/'measurement.json').read_text())
        record['verified_identical_rows'] = compare_tables(n)
        record['speedup'] = record['before']['total_s']/record['after']['total_s']
        record['peak_rss_reduction_percent'] = 100*(1-record['after']['peak_rss_bytes']/record['before']['peak_rss_bytes'])
        comparison.write_text(json.dumps(records, indent=2))
        print('Final optimized version verified: ' + str(n), flush=True)
    environment_path = ROOT/'benchmarks'/'results'/'environment.json'
    environment = json.loads(environment_path.read_text())
    environment['optimized_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                       for p in (ROOT/'rnai_designer').glob('*.py')}
    environment_path.write_text(json.dumps(environment, indent=2))
    commands = [
        ['benchmarks/benchmark.py', '--mode', 'after', '--lengths', '500', '--profile'],
        ['benchmarks/benchmark.py', '--mode', 'before', '--lengths', '500', '--trace'],
        ['benchmarks/benchmark.py', '--mode', 'after', '--lengths', '500', '--trace'],
        ['benchmarks/parallel_probe.py'],
    ]
    for command in commands:
        print('Running instrumented probe: ' + ' '.join(command), flush=True)
        subprocess.run([sys.executable, *command], cwd=ROOT, check=True)
    env = os.environ.copy()
    env['SIFI_BOWTIE'] = str(ROOT.parent/'siFi21--1.2.3'/'ToCopy'/'Bowtie'/'bowtie.exe')
    env['SIFI_BOWTIE_INSPECT'] = str(ROOT.parent/'siFi21--1.2.3'/'ToCopy'/'Bowtie'/'bowtie-inspect.exe')
    with (ROOT/'benchmarks'/'results'/'tests.txt').open('w', encoding='utf-8') as handle:
        subprocess.run([sys.executable, '-m', 'pytest', '-q', '-rs', '--basetemp',
                        str(ROOT/'benchmarks'/('test-temp-final-'+uuid4().hex))], cwd=ROOT, env=env,
                       stdout=handle, stderr=subprocess.STDOUT, check=True)
    (ROOT/'benchmarks'/'results'/'instrumented_complete.json').write_text(json.dumps({'completed': True}))
