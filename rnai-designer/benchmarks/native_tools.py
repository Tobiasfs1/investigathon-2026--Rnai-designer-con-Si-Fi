"""Run packaged, unchanged scientific executables against temporary synthetic data only."""
from contextlib import nullcontext
import importlib
import json
import random
import subprocess
import tempfile
import time
from pathlib import Path
from benchmark import ROOT, fixture


TOOLS = ROOT.parent/'siFi21--1.2.3'/'ToCopy'
BOWTIE = TOOLS/'Bowtie'/'bowtie.exe'
BUILDER = TOOLS/'Bowtie'/'bowtie-build.exe'
RNAPLFOLD = TOOLS/'RNAplfold'/'RNAplfold.exe'


def run(n):
    from rnai_designer.data_models import Candidate
    from rnai_designer.config import Settings
    from rnai_designer.reuse import offtarget_reuse_session
    sequence = fixture(n)['candidate']['reference']['sequence']
    candidate = Candidate.from_fields('SYNTHETIC_ONLY', sequence, 'synthetic', '817')
    candidate.original['metadata']['simulated'] = True
    result = {'nt': n, 'headless_conda_s': None, 'note': 'Conda unavailable; native tools measured directly, without changing their algorithms.'}
    with tempfile.TemporaryDirectory(prefix='native-synthetic-', dir=ROOT/'benchmarks') as temporary:
        root = Path(temporary)
        db = root/'synthetic-db'
        db.mkdir()
        bindings = {}
        t = time.perf_counter()
        for i, (species, slug) in enumerate((('Homo sapiens', 'human'), ('Zea mays', 'maize'), ('Apis mellifera', 'bee'))):
            rng = random.Random(817+i)
            reference = ''.join(rng.choices('ACGT', k=n))
            if i != 1:
                reference += sequence[:71] if i == 0 else sequence[-71:]
            fasta = root/(slug+'.fa')
            fasta.write_text('>SYNTHETIC_'+slug+'\n'+reference+'\n', encoding='ascii')
            subprocess.run([str(BUILDER), str(fasta), str(db/slug)], capture_output=True, check=True,
                           creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
            bindings[species] = {'prefix': slug, 'source': 'explicitly synthetic temporary FASTA', 'version': '817'}
        result['synthetic_index_build_s'] = time.perf_counter()-t
        settings = Settings(db, BOWTIE, BOWTIE.with_name('bowtie-inspect.exe'), RNAPLFOLD, bindings)
        criteria = {}
        for package, context in (('reference', nullcontext), ('rnai_designer', offtarget_reuse_session)):
            module = importlib.import_module(package+'.automation')
            actual_execution_times = []
            original_execute = module.execute
            def timed_execute(*args, _execute=original_execute, **kwargs):
                started = time.perf_counter()
                try:
                    return _execute(*args, **kwargs)
                finally:
                    actual_execution_times.append(time.perf_counter()-started)
            module.execute = timed_execute
            measurements = []
            with context():
                for reference in (1, 2):
                    actual_execution_times.clear()
                    t = time.perf_counter()
                    evaluations = module.run_offtargets(candidate, settings, root/package/str(reference))
                    measurements.append(time.perf_counter()-t)
                    result[package+f'_actual_bowtie_execution_{reference}_s'] = sum(actual_execution_times)
                    result[package+f'_actual_bowtie_execution_{reference}_calls'] = len(actual_execution_times)
                    assert all(e['status'] == 'imported_unvalidated' for e in evaluations)
                    science = [e['criterion'] for e in evaluations]
                    if package not in criteria:
                        criteria[package] = science
                    assert science == criteria[package]
                    if package == 'rnai_designer' and reference == 2:
                        assert all('reused_from' in e for e in evaluations)
            result[package+'_three_species_first_s'] = measurements[0]
            module.execute = original_execute
            result[package+'_same_gene_second_reference_s'] = measurements[1]
            result[package+'_output_bytes'] = sum(p.stat().st_size for p in (root/package).rglob('*') if p.is_file())
        assert criteria['reference'] == criteria['rnai_designer']
        result['scientific_criteria_identical'] = True
        rna_out = root/'rna'
        rna_out.mkdir()
        t = time.perf_counter()
        subprocess.run([str(RNAPLFOLD), '-W', '80', '-L', '40', '-u', '21', '-T', '22.00'],
                       input='>SYNTHETIC_ONLY\n'+sequence+'\n\n', text=True, capture_output=True,
                       cwd=rna_out, check=True, creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
        result['rnaplfold_direct_s'] = time.perf_counter()-t
        assert (rna_out/'SYNTHETIC_ONLY_lunp').is_file()
        result['rnaplfold_output_bytes'] = sum(p.stat().st_size for p in rna_out.iterdir() if p.is_file())
    print(json.dumps(result), flush=True)
    return result


if __name__ == '__main__':
    records = [run(n) for n in (500, 1000, 2500, 5000)]
    (ROOT/'benchmarks'/'results'/'native_tools.json').write_text(json.dumps(records, indent=2))
