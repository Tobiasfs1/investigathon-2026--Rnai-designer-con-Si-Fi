"""Profile additional pipeline stages using explicit synthetic data, no installed databases."""
import importlib
import json
import sys
import tempfile
import time
from pathlib import Path
from benchmark import fixture, ROOT


def run(package, n):
    def mod(name):
        return importlib.import_module(package+'.'+name)
    measurements = {'nt': n, 'package': package,
                    'headless_conda_s': None, 'external_bowtie_s': None,
                    'external_note': 'External execution is measured separately by native_tools.py; Conda headless is unavailable.'}
    with tempfile.TemporaryDirectory(prefix='rnai-synthetic-', dir=ROOT/'benchmarks') as temporary:
        root = Path(temporary)
        report = fixture(n)
        t = time.perf_counter()
        candidate = mod('data_models').Candidate.from_fields('SYNTHETIC_ONLY', report['candidate']['reference']['sequence'], 'synthetic', '817')
        candidate.original['metadata']['simulated'] = True
        measurements['sequence_preparation_validation_s'] = time.perf_counter()-t
        raw = root/'sifi.json'
        raw.write_text(json.dumps(report['sifi']['records']))
        t = time.perf_counter()
        imported = mod('sifi_adapter').import_sifi(raw, candidate, report['sifi']['run'])
        measurements['sifi_json_import_validation_s'] = time.perf_counter()-t
        t = time.perf_counter()
        summary = mod('efficiency_analysis').summarize(imported['records'])
        measurements['descriptive_summary_s'] = time.perf_counter()-t
        t = time.perf_counter()
        mod('window_explorer').green_positions(imported['records'], 21)
        measurements['green_mask_s'] = time.perf_counter()-t
        alignment = root/'synthetic.tsv'
        with alignment.open('w') as handle:
            for i in range(1, n-19, 10):
                handle.write(f'sirna{i}\t+\tSYNTHETIC_REFERENCE\t{i}\t{candidate.sequence[i-1:i+20]}\tIIIIIIIIIIIIIIIIIIIII\t0\t\n')
        t = time.perf_counter()
        parsed = mod('bowtie_adapter').parse_bowtie(alignment, candidate, {'query_id': candidate.gene_id})
        measurements['bowtie_result_parsing_s'] = time.perf_counter()-t
        t = time.perf_counter()
        mod('offtarget_assessment').assess({'organism': 'Homo sapiens', 'status': 'imported_unvalidated',
                                          'results_format': 'bowtie_legacy_tab', 'result': parsed}, candidate)
        measurements['offtarget_assessment_s'] = time.perf_counter()-t
        t = time.perf_counter()
        with (root/'serialized.json').open('w') as handle:
            json.dump(report, handle)
        measurements['json_serialization_s'] = time.perf_counter()-t
        evaluations = mod('offtarget_reports').load_evaluations(None, candidate)
        t = time.perf_counter()
        mod('reporting').write_reports(root/'integrated', candidate, imported, summary, evaluations)
        measurements['integrated_reports_csv_json_figures_s'] = time.perf_counter()-t
        # Consolidation only: the same compact per-reference table stream, two independent targets.
        t = time.perf_counter()
        import csv
        count = 0
        with (root/'references.csv').open('w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(['main_target', 'start', 'end'])
            for target in ('SYNTHETIC_A', 'SYNTHETIC_B'):
                for length in range(180, 251):
                    for start in range(1, n-length+2):
                        writer.writerow([target, start, start+length-1])
                        count += 1
        measurements['reference_consolidation_scalar_io_s'] = time.perf_counter()-t
        measurements['consolidated_rows'] = count
    print(json.dumps(measurements), flush=True)
    return measurements


if __name__ == '__main__':
    records = [run(package, n) for package in ('reference', 'rnai_designer') for n in (500, 1000, 2500, 5000)]
    (ROOT/'benchmarks'/'results'/'stages.json').write_text(json.dumps(records, indent=2))
