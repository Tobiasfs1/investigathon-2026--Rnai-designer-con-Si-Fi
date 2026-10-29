"""Measure Windows spawn plus transfer of complete compact results (not just workers)."""
from concurrent.futures import ProcessPoolExecutor
import json
import time
from benchmark import fixture, ROOT, peak_rss
from rnai_designer.automatic_window_selection import select_windows


def work(n):
    result = select_windows(fixture(n))
    return result, peak_rss()


if __name__ == '__main__':
    start = time.perf_counter()
    sequential = [work(500), work(500)]
    sequential_s = time.perf_counter()-start
    start = time.perf_counter()
    with ProcessPoolExecutor(max_workers=2) as pool:
        parallel = list(pool.map(work, [500, 500]))
    parallel_s = time.perf_counter()-start
    for (a, _), (b, _) in zip(sequential, parallel):
        assert a['selected_window'] == b['selected_window']
        assert list(a['windows'].scalar_rows()) == list(b['windows'].scalar_rows())
    result = {'nt_per_reference': 500, 'references': 2, 'sequential_s': sequential_s,
              'spawn_and_transfer_s': parallel_s, 'speedup': sequential_s/parallel_s,
              'child_peak_rss_bytes': [rss for _, rss in parallel], 'parent_peak_rss_bytes': peak_rss(),
              'scope': 'Python window calculation only, not concurrent si-Fi/Conda or Bowtie'}
    (ROOT/'benchmarks'/'results'/'parallel_probe.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(result))
