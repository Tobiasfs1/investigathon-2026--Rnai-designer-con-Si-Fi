"""Synthetic only. Run from project root: python benchmarks/benchmark.py --mode before."""
import argparse
import cProfile
import importlib
import json
import random
import sys
import time
import tracemalloc
import os
import ctypes
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'benchmarks'))


def peak_rss(detail=False):
    if os.name != 'nt':
        import resource
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024)
        return {'rss': rss, 'commit': None} if detail else rss
    from ctypes import wintypes
    class Counters(ctypes.Structure):
        _fields_ = [('cb', wintypes.DWORD), ('PageFaultCount', wintypes.DWORD),
                    *[(name, ctypes.c_size_t) for name in ('PeakWorkingSetSize', 'WorkingSetSize',
                      'QuotaPeakPagedPoolUsage', 'QuotaPagedPoolUsage', 'QuotaPeakNonPagedPoolUsage',
                      'QuotaNonPagedPoolUsage', 'PagefileUsage', 'PeakPagefileUsage')]]
    counters = Counters()
    counters.cb = ctypes.sizeof(counters)
    kernel = ctypes.windll.kernel32
    kernel.GetCurrentProcess.restype = wintypes.HANDLE
    get_info = ctypes.windll.psapi.GetProcessMemoryInfo
    get_info.argtypes = [wintypes.HANDLE, ctypes.c_void_p, wintypes.DWORD]
    get_info.restype = wintypes.BOOL
    if not get_info(kernel.GetCurrentProcess(), ctypes.byref(counters), counters.cb):
        raise ctypes.WinError()
    return ({'rss': counters.PeakWorkingSetSize, 'commit': counters.PeakPagefileUsage}
            if detail else counters.PeakWorkingSetSize)


def fixture(n):
    rng = random.Random(817)
    sequence = ''.join(rng.choices('ACGT', k=n))
    candidate = {'gene_id': 'SYNTHETIC_ONLY', 'priority_score': None,
                 'metadata': {'simulated': True},
                 'reference': {'sequence': sequence, 'source': 'synthetic', 'version': '817'}}
    records = [dict(query_name='SYNTHETIC_ONLY', sirna_name=f'sirna{i}',
                    sirna_position=i, sirna_sequence=sequence[i-1:i+20],
                    is_off_target=False, hit_name='synthetic-target',
                    is_efficient=i % 3 != 0, strand_selection=True, end_stability=False,
                    target_site_accessibility=True, accessibility_value=(i % 31)/32)
               for i in range(1, n-19)]
    evaluations = [{'organism': 'Homo sapiens', 'status': 'imported_unvalidated',
                    'criterion': {'status': 'exact_match_detected', 'alignments': [
                        {'query_start_1based': i, 'query_end_1based_inclusive': i+20,
                         'length': 21, 'contains_exact_21nt': True, 'match_type': 'exact_21nt',
                         'reference_id': f'SYNTHETIC_{i}'} for i in range(1, n-19, 10)]}}]
    return {'candidate': candidate, 'sifi': {'records': records, 'status': 'evaluated',
            'source': {'sha256': 'synthetic-only'}, 'warnings': [],
            'run': {'query_id': 'SYNTHETIC_ONLY', 'parameters': {'sirna_size': 21}}},
            'offtarget_evaluations': evaluations}


def run(mode, n, trace=False, profile=False):
    package = 'reference' if mode == 'before' else 'rnai_designer'
    selection_module = importlib.import_module(package + '.automatic_window_selection')
    reporting = importlib.import_module(package + '.automatic_window_reporting')
    out = ROOT / 'benchmarks' / 'results' / mode / str(n)
    out.mkdir(parents=True, exist_ok=True)
    stages = {}
    t = time.perf_counter()
    report = fixture(n)
    stages['synthetic_preparation_s'] = time.perf_counter()-t
    report_path = out / 'input.json'
    report_path.write_text(json.dumps(report), encoding='utf-8')
    t = time.perf_counter()
    data_models = importlib.import_module(package + '.data_models')
    adapter = importlib.import_module(package + '.sifi_adapter')
    efficiency = importlib.import_module(package + '.efficiency_analysis')
    candidate = data_models.Candidate.from_fields('SYNTHETIC_ONLY', report['candidate']['reference']['sequence'], 'synthetic', '817')
    candidate.original['metadata']['simulated'] = True
    stages['sequence_validation_s'] = time.perf_counter()-t
    t = time.perf_counter()
    data_models.read_json(report_path)
    stages['json_read_s'] = time.perf_counter()-t
    raw_path = out/'raw_sifi.json'
    raw_path.write_text(json.dumps(report['sifi']['records']))
    t = time.perf_counter()
    imported = adapter.import_sifi(raw_path, candidate, report['sifi']['run'])
    stages['sifi_import_s'] = time.perf_counter()-t
    t = time.perf_counter()
    efficiency.summarize(imported['records'])
    stages['efficiency_summary_s'] = time.perf_counter()-t
    plot_time = [0.0]
    original_plot = reporting._plot_selected
    def timed_plot(*args, **kwargs):
        start = time.perf_counter()
        try:
            return original_plot(*args, **kwargs)
        finally:
            plot_time[0] += time.perf_counter()-start
    reporting._plot_selected = timed_plot
    if trace:
        tracemalloc.start()
    profiler = cProfile.Profile()
    if profile:
        profiler.enable()
    started = time.perf_counter()
    selection = selection_module.select_windows(report)
    computed = time.perf_counter()
    # Fresh destination required by the public writer.
    destination = out / ('export_' + str(time.time_ns()))
    reporting.write_selection_outputs(report, selection, report_path, destination)
    ended = time.perf_counter()
    reporting._plot_selected = original_plot
    if profile:
        profiler.disable()
        profiler.dump_stats(str(out / 'profile.pstats'))
        import pstats
        with (out/'profile.txt').open('w', encoding='utf-8') as handle:
            pstats.Stats(profiler, stream=handle).strip_dirs().sort_stats('cumulative').print_stats(60)
    peak = tracemalloc.get_traced_memory()[1] if trace else None
    if trace:
        snapshot = tracemalloc.take_snapshot()
        (out/'allocations.txt').write_text('\n'.join(str(stat) for stat in snapshot.statistics('lineno')[:30]), encoding='utf-8')
        tracemalloc.stop()
    memory = peak_rss(detail=True)
    row = {'mode': mode, 'nt': n, 'windows': selection['evaluated_window_count'],
           'compute_s': computed-started, 'export_s': ended-computed, 'total_s': ended-started,
           'tracemalloc_peak_bytes': peak, 'traced': trace, 'profiled': profile,
           'peak_rss_bytes': memory['rss'], 'peak_commit_bytes': memory['commit'],
           'stages': stages, 'plot_s': plot_time[0],
           'table_and_json_export_s': ended-computed-plot_time[0],
           'output_bytes': sum(p.stat().st_size for p in destination.iterdir()),
           'csv_bytes': (destination/'all_windows.csv').stat().st_size}
    name = 'measurement_traced.json' if trace else 'measurement_profiled.json' if profile else 'measurement.json'
    (out / name).write_text(json.dumps(row, indent=2))
    print(json.dumps(row), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['before', 'after'], required=True)
    parser.add_argument('--lengths', nargs='+', type=int, default=[500, 1000, 2500, 5000])
    parser.add_argument('--trace', action='store_true')
    parser.add_argument('--profile', action='store_true')
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()
    if not args.worker:
        import subprocess
        for n in args.lengths:
            subprocess.run([sys.executable, __file__, '--mode', args.mode, '--lengths', str(n),
                            '--worker', *(['--trace'] if args.trace else []),
                            *(['--profile'] if args.profile else [])], check=True)
        sys.exit(0)
    for n in args.lengths:
        run(args.mode, n, args.trace, args.profile)
