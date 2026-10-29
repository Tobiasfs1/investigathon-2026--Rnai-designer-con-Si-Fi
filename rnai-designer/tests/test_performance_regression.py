"""Differential checks against the frozen pre-optimization implementation."""
import copy
import csv
import json
import random
from pathlib import Path

import pytest
from benchmarks.benchmark import fixture
from benchmarks.reference.automatic_window_selection import select_windows as original
from rnai_designer.automatic_window_selection import select_windows
from rnai_designer.automatic_window_reporting import write_selection_outputs
from rnai_designer.atomic_io import atomic_text
from rnai_designer.window_explorer import _overlapping_offtargets


@pytest.mark.parametrize('variant', ['normal', 'missing', 'conflicts', 'short', 'gaps', 'empty', 'ambiguous', 'mixed_numeric'])
def test_every_metric_and_winner_identical(variant):
    report = fixture(270)
    rows = report['sifi']['records']
    rng = random.Random(81)
    for row in rows:
        row['accessibility_value'] = rng.random()  # Non-binary fractions: exact sum order matters.
    if variant == 'missing':
        for row in rows[::3]:
            row['is_efficient'] = row['accessibility_value'] = None
    elif variant == 'conflicts':
        rows.extend([dict(rows[0]), dict(rows[3], accessibility_value=.999, hit_name='other')])
    elif variant == 'ambiguous':
        seq = 'NURY' * 67 + 'NN'
        report['candidate']['reference']['sequence'] = seq
        for row in rows:
            i = row['sirna_position']-1
            row['sirna_sequence'] = seq[i:i+21]
    elif variant == 'mixed_numeric':
        for i, row in enumerate(rows):
            row['accessibility_value'] = 1 if i % 13 == 0 else .25 if i % 2 else 0
    elif variant in ('short', 'gaps', 'empty'):
        positions = list(range(1, 100)) if variant == 'short' else list(range(1, 200)) + list(range(220, 271))
        report['sifi']['green_mask'] = {'positions_1based': positions, 'sirna_size': 21,
                                      'source_function': 'general_helpers.get_target_data',
                                      'coordinate_system': '1-based positions'}
        if variant == 'empty':
            report['sifi']['records'] = []
    old, new = original(report), select_windows(report)
    assert old['evaluated_window_count'] == new['evaluated_window_count']
    for a, b in zip(old['windows'], new['windows'], strict=True):
        details = a.pop('offtarget_matches')
        assert a == {k: v for k, v in b.items() if k not in ("length_penalty", "adjusted_priority_score")}  # Original science remains exact.
        assert details == _overlapping_offtargets(report['offtarget_evaluations'], b['start_1based'], b['end_1based_inclusive'], 270)
    assert old['selected_window'] == ({k: v for k, v in new['selected_window'].items() if k not in ('length_penalty', 'adjusted_priority_score')} if new['selected_window'] else None)
    assert old['warnings'] == new['warnings']


def test_export_is_scalar_and_failure_never_publishes_directory(tmp_path, monkeypatch):
    import rnai_designer.automatic_window_reporting as module
    report = fixture(181)
    source = tmp_path/'input.json'
    source.write_text(json.dumps(report))
    selection = select_windows(report)
    original_plot = module._plot_selected
    def interrupt(*args):
        raise KeyboardInterrupt('synthetic interruption')
    monkeypatch.setattr(module, '_plot_selected', interrupt)
    with pytest.raises(KeyboardInterrupt):
        write_selection_outputs(report, selection, source, tmp_path/'export')
    assert not (tmp_path/'export').exists()
    assert not list(tmp_path.glob('*.partial'))
    monkeypatch.setattr(module, '_plot_selected', original_plot)
    write_selection_outputs(report, selection, source, tmp_path/'export')
    with (tmp_path/'export'/'all_windows.csv').open(newline='') as handle:
        records = list(csv.DictReader(handle))
    assert len(records) == selection['evaluated_window_count']
    assert all(not value.startswith(('{', '[')) for r in records for value in r.values())
    assert 'offtarget_matches' not in records[0]


def test_atomic_write_preserves_previous_valid_file(tmp_path):
    path = tmp_path/'result.json'
    path.write_text('previous')
    with pytest.raises(OSError):
        with atomic_text(path) as handle:
            handle.write('unfinished')
            raise OSError('disk full')
    assert path.read_text() == 'previous'


def test_cache_keys_invalidate_inputs_parameters_versions_and_simulation(tmp_path):
    from rnai_designer.cache import content_key, file_digest, ArtifactCache
    source = tmp_path/'source.json'
    source.write_text('[1]')
    params = {'simulated': True, 'version': '1', 'stride': 1}
    key = content_key(params, [source])
    for change in ({'simulated': False}, {'version': '2'}, {'stride': 2}):
        assert key != content_key({**params, **change}, [source])
    cache = ArtifactCache(tmp_path/'cache', key)
    cache.store({'source.json': source})
    target = tmp_path/'restored.json'
    assert cache.restore({'source.json': target})
    assert target.read_bytes() == source.read_bytes()
    (cache.directory/'source.json').write_text('[9]')
    assert not cache.restore({'source.json': target})
    source.write_text('[2, 3]')
    assert key != content_key(params, [source])
    assert file_digest(source) != file_digest(target)


def test_stale_or_incomplete_cache_never_counts_as_result(tmp_path):
    from rnai_designer.cache import ArtifactCache
    cache = ArtifactCache(tmp_path, 'incomplete')
    cache.directory.mkdir()
    (cache.directory/'source.json').write_text('[]')
    assert not cache.restore({'source.json': tmp_path/'restored'})
    assert not (tmp_path/'restored').exists()


def test_order_sensitive_green_mask_matches_original_for_shuffled_hits():
    from rnai_designer.window_explorer import green_positions
    from benchmarks.reference.window_explorer import green_positions as old_green
    for seed in range(30):
        report = fixture(120)
        rng = random.Random(seed)
        records = report['sifi']['records']
        for row in records:
            row['is_off_target'] = rng.choice([True, False, None])
            row['hit_name'] = rng.choice(['A', 'B', 'C'])
        records.extend(copy.deepcopy(records[::9]))
        rng.shuffle(records)
        assert green_positions(records, 21) == old_green(records, 21)


def test_tiny_region_with_duplicates_before_start_and_detail_lookup():
    from rnai_designer.window_details import resolve_window_details
    report = fixture(270)
    report['sifi']['records'].extend(copy.deepcopy(report['sifi']['records'][70:100]))
    report['sifi']['green_mask'] = {'positions_1based': list(range(100, 106)), 'sirna_size': 21,
                                  'source_function': 'general_helpers.get_target_data',
                                  'coordinate_system': '1-based positions'}
    old, new = original(report), select_windows(report)
    a, b = old['windows'][0], new['windows'][0]
    detail = a.pop('offtarget_matches')
    assert a == {k: v for k, v in b.items() if k not in ("length_penalty", "adjusted_priority_score")}
    assert b['exact_duplicate_record_count'] == 0
    resolved = resolve_window_details(report, b['window_id'])
    assert resolved['sirnas'] == []
    assert resolved['offtarget_matches'] == detail
    with pytest.raises(ValueError):
        resolve_window_details(report, 'green-001:90-105')


def test_completed_search_reuse_is_scoped_verified_and_never_crosses_simulation(tmp_path, monkeypatch):
    import rnai_designer.automation as automation
    from rnai_designer.reuse import offtarget_reuse_session
    from rnai_designer.config import Settings, index_files
    from rnai_designer.data_models import Candidate
    from types import SimpleNamespace
    db = tmp_path/'db'
    db.mkdir()
    tool = tmp_path/'bowtie.exe'
    tool.write_bytes(b'synthetic-test-placeholder-not-an-executable')
    bindings = {species: {'prefix': slug, 'source': 'synthetic', 'version': '1'}
                for species, slug in automation.SPECIES_SLUGS.items()}
    for binding in bindings.values():
        for path in index_files(db, binding['prefix']):
            path.write_bytes(b'synthetic-index-placeholder')
    settings = Settings(db, tool, None, None, bindings)
    candidate = Candidate.from_fields('SYNTHETIC_ONLY', 'ACGT'*45, 'synthetic', '1')
    candidate.original['metadata']['simulated'] = True
    calls = []
    monkeypatch.setattr(automation.subprocess, 'run', lambda *a, **kw:
                        SimpleNamespace(returncode=0, stdout='bowtie version 0.11.3', stderr=''))
    def execute(command, directory, timeout):
        calls.append(command)
        (directory/'alignments.tsv').write_text('')
        (directory/'stdout.log').write_text('')
        (directory/'stderr.log').write_text('# reads processed: 160\n')
        return {'status': 'completed', 'returncode': 0, 'command': command}
    monkeypatch.setattr(automation, 'execute', execute)
    with offtarget_reuse_session():
        a = automation.run_offtargets(candidate, settings, tmp_path/'first')
        b = automation.run_offtargets(candidate, settings, tmp_path/'second')
        assert len(calls) == 3
        assert all('reused_from' in e for e in b)
        assert [e['criterion'] for e in a] == [e['criterion'] for e in b]
        assert a[0]['result']['records'] is b[0]['result']['records']
        assert (tmp_path/'second'/'reuse_manifest.json').is_file()
        # A real-mode invocation can never consume an explicitly synthetic cached search.
        real_mode = copy.deepcopy(candidate)
        real_mode.original['metadata']['simulated'] = False
        automation.run_offtargets(real_mode, settings, tmp_path/'real-mode-contract-test')
        assert len(calls) == 6
        # Changed database bytes invalidate all previous results.
        index_files(db, 'human')[0].write_bytes(b'changed synthetic index')
        automation.run_offtargets(candidate, settings, tmp_path/'changed-db')
        assert len(calls) == 9
        # Corrupted source evidence invalidates a previously completed entry.
        (tmp_path/'changed-db'/'human'/'alignments.tsv').write_text('corrupt')
        automation.run_offtargets(candidate, settings, tmp_path/'corrupt-evidence')
        assert len(calls) == 12
    automation.run_offtargets(candidate, settings, tmp_path/'outside-session')
    assert len(calls) == 15


def test_headless_cache_reuses_across_output_directories_and_invalidates_parameters(tmp_path, monkeypatch):
    import rnai_designer.__main__ as cli
    from rnai_designer.config import index_files
    from rnai_designer.data_models import Candidate
    from types import SimpleNamespace
    report = fixture(21)
    candidate = Candidate.from_fields('SYNTHETIC_ONLY', report['candidate']['reference']['sequence'], 'synthetic', '1')
    candidate.original['metadata']['simulated'] = True
    db, engine, prefix = (tmp_path/name for name in ('db', 'engine', 'env'))
    for path in (db, engine, prefix):
        path.mkdir()
    (engine/'original.py').write_text('# synthetic placeholder')
    conda, tool = tmp_path/'conda.exe', tmp_path/'bowtie.exe'
    for path in (conda, tool, prefix/'python.exe'):
        path.write_bytes(b'synthetic-test-placeholder')
    for path in index_files(db, 'synthetic'):
        path.write_bytes(b'synthetic-index-placeholder')
    runtime = {'prefix': str(prefix), 'conda': str(conda), 'python': '2.7', 'bits': 64}
    settings = SimpleNamespace(database_dir=db, databases={'Dalbulus maidis': {'prefix': 'synthetic'}},
                               bowtie=tool, rnaplfold=tool)
    args = SimpleNamespace(cache_dir=tmp_path/'cache', sifi_root=engine, legacy_sirna_size=21,
                           sifi_mismatches=0, accessibility_window=8, end_stability_threshold=1.,
                           target_accessibility_threshold=.1, sifi_timeout=10)
    def context(args, candidate, out, settings):
        inputs = out/'inputs'
        inputs.mkdir(parents=True)
        query = out/'designer_handoff'/'query.fasta'
        query.parent.mkdir()
        query.write_text('>SYNTHETIC_ONLY\n'+candidate.sequence+'\n')
        source, mask = inputs/'designer_original.json', inputs/'green_mask.json'
        return runtime, query, source, mask, ['synthetic-conda', '--out-json', str(source), '--out-mask', str(mask)]
    monkeypatch.setattr(cli, '_legacy_sifi_context', context)
    calls = []
    def run(command, **kwargs):
        calls.append(command)
        Path(command[command.index('--out-json')+1]).write_text(json.dumps(report['sifi']['records']))
        Path(command[command.index('--out-json')+1]+'.bowtie.tmp').write_text('SYNTHETIC CACHE TRANSPORT TEST')
        Path(command[command.index('--out-mask')+1]).write_text(json.dumps({
            'positions_1based': list(range(1, 22)), 'sirna_size': 21,
            'source_function': 'general_helpers.get_target_data', 'coordinate_system': '1-based positions'}))
        return SimpleNamespace(returncode=0, stdout='', stderr='')
    monkeypatch.setattr(cli.subprocess, 'run', run)
    a = cli._run_legacy_sifi(args, candidate, tmp_path/'a', settings, 'synthetic-target')
    b = cli._run_legacy_sifi(args, candidate, tmp_path/'b', settings, 'synthetic-target')
    assert len(calls) == 1
    assert a[1]['records'] == b[1]['records']
    assert Path(str(a[0])+'.bowtie.tmp').read_bytes() == Path(str(b[0])+'.bowtie.tmp').read_bytes()
    assert (tmp_path/'b'/'headless_performance.json').is_file()
    args.accessibility_window = 9
    cli._run_legacy_sifi(args, candidate, tmp_path/'c', settings, 'synthetic-target')
    assert len(calls) == 2


def test_batched_plot_preserves_png_pixels(tmp_path):
    from PIL import Image, ImageChops
    from benchmarks.reference.automatic_window_reporting import _plot_selected as old_plot
    from rnai_designer.automatic_window_reporting import _plot_selected as new_plot
    report = fixture(500)
    selection = select_windows(report)
    old_plot(tmp_path/'before.svg', report, selection)
    new_plot(tmp_path/'after.svg', report, selection)
    with Image.open(tmp_path/'before.png') as a, Image.open(tmp_path/'after.png') as b:
        assert a.size == b.size
        assert ImageChops.difference(a.convert('RGB'), b.convert('RGB')).getbbox() is None


def test_partial_reference_csv_is_rolled_back_before_next_reference(tmp_path, monkeypatch):
    import rnai_designer.__main__ as cli
    from rnai_designer.data_models import Candidate
    from types import SimpleNamespace
    candidate = Candidate.from_fields('SYNTHETIC_ONLY', 'ACGT'*50, 'synthetic', '1')
    candidate.original['metadata']['simulated'] = True
    args = SimpleNamespace(timeout=10, min_window=180, max_window=250, sirna_size=21, _sifi_runtime={})
    def fake_sifi(args, candidate, out, settings, target):
        source = out/'inputs'/'source.json'
        source.parent.mkdir()
        source.write_text('[]')
        return source, {'records': []}, {}
    def fake_pipeline(*args):
        return {'summary': {'green_regions': [], 'evaluated_window_count': 2, 'selected_window': None,
                            'selected_window_offtarget_filter': {'passes': False}, 'warnings': []},
                'selection': {'selection_method': {}, 'windows': [
                    {'window_id': 'green-001:1-180', 'start_1based': 1},
                    {'window_id': 'green-001:2-181', 'start_1based': 2}]}}
    monkeypatch.setattr(cli, '_run_legacy_sifi', fake_sifi)
    monkeypatch.setattr(cli, 'run_integrated_workflow', fake_pipeline)
    real_writer = csv.writer
    failed = [False]
    class Writer:
        def __init__(self, handle, *a, **kw):
            self.writer = real_writer(handle, *a, **kw)
        def writerow(self, row):
            return self.writer.writerow(row)
        def writerows(self, rows):
            if not failed[0]:
                self.writer.writerow(next(iter(rows)))
                failed[0] = True
                raise OSError('synthetic transient write failure after first row')
            return self.writer.writerows(rows)
    monkeypatch.setattr(csv, 'writer', Writer)
    out = tmp_path/'multi'
    out.mkdir()
    report = cli._run_multi_reference_workflow(args, candidate, object(), out, ['A', 'B'])
    assert report['failed_reference_count'] == 1
    assert report['completed_reference_count'] == 1
    assert report['combined_window_count'] == 2
    with (out/'all_references_windows.csv').open(newline='') as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2
    assert {r['main_target'] for r in rows} == {'B'}


def test_integrated_report_overwrite_preserves_evidence_and_rolls_back_publish_failure(tmp_path, monkeypatch):
    from rnai_designer.reporting import write_reports
    from rnai_designer.efficiency_analysis import summarize
    from rnai_designer.offtarget_reports import load_evaluations
    from rnai_designer.data_models import Candidate
    import rnai_designer.atomic_io as atomic
    report = fixture(21)
    candidate = Candidate.from_fields('SYNTHETIC_ONLY', report['candidate']['reference']['sequence'], 'synthetic', '1')
    candidate.original['metadata']['simulated'] = True
    imported = report['sifi']
    summary = summarize(imported['records'])
    evaluations = load_evaluations(None, candidate)
    out = tmp_path/'report'
    evidence = out/'offtarget_runs'/'keep.txt'
    evidence.parent.mkdir(parents=True)
    evidence.write_text('original scientific evidence')
    write_reports(out, candidate, imported, summary, evaluations)
    old_report = (out/'report.json').read_bytes()
    old_records = (out/'records.csv').read_bytes()
    replace = atomic.os.replace
    failed = [False]
    def fail_one_publish(source, destination):
        if Path(source).name == 'records.csv' and Path(destination) == out/'records.csv' and not failed[0]:
            failed[0] = True
            raise OSError('synthetic publication failure')
        return replace(source, destination)
    monkeypatch.setattr(atomic.os, 'replace', fail_one_publish)
    with pytest.raises(OSError):
        write_reports(out, candidate, imported, summary, evaluations)
    assert (out/'report.json').read_bytes() == old_report
    assert (out/'records.csv').read_bytes() == old_records
    assert evidence.read_text() == 'original scientific evidence'
    write_reports(out, candidate, imported, summary, evaluations)
    assert evidence.read_text() == 'original scientific evidence'


def test_legacy_list_container_exports_identical_compact_csv(tmp_path):
    report = fixture(181)
    source = tmp_path/'input.json'
    source.write_text(json.dumps(report))
    selection = select_windows(report)
    write_selection_outputs(report, selection, source, tmp_path/'columnar')
    legacy_container = {**selection, 'windows': list(selection['windows'])}
    write_selection_outputs(report, legacy_container, source, tmp_path/'list')
    assert (tmp_path/'columnar'/'all_windows.csv').read_bytes() == (tmp_path/'list'/'all_windows.csv').read_bytes()


def test_search_reuse_memory_is_bounded_and_released_between_workflows(tmp_path):
    from rnai_designer.reuse import offtarget_reuse_session, search_cache, remember_search
    with offtarget_reuse_session():
        cache = search_cache.get()
        for i in range(12):
            source = tmp_path/str(i)
            source.mkdir()
            (source/'synthetic.txt').write_text('SYNTHETIC CACHE LIFETIME TEST')
            remember_search(cache, str(i), source, [{'status': 'imported_unvalidated', 'result': {'records': []}}])
        assert len(cache) == 4
        assert '0' not in cache
        assert '11' in cache
    assert search_cache.get() is None
