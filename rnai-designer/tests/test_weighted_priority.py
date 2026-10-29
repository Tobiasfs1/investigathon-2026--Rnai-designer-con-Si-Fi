"""Small scoring checks only: no si-Fi, Bowtie, RNAplfold or real batch."""
import csv
import json

import pytest

from rnai_designer.automatic_window_selection import _efficacy_priority, select_windows
from rnai_designer.automatic_window_reporting import write_selection_outputs
from rnai_designer.batch import result_row, save_reports
from rnai_designer.priority_scoring import length_scores, priority_key
from test_automatic_window_selection import make_report

CASES = [(1, 0, 75), (0, 1, 25), (1, 1, 100), (0, 0, 0)]


def single_window(efficiency, gc):
    sequence = "G"*7 + "A"*14 if gc else "A"*21
    report = make_report([(1, 21)], sequence_length=21)
    report['candidate']['reference']['sequence'] = sequence
    report['sifi']['records'][0].update(sirna_sequence=sequence, is_efficient=bool(efficiency))
    return report, select_windows(report)


@pytest.mark.parametrize('efficiency,gc,expected', CASES)
def test_requested_scores_in_helper_and_indexed_selection(efficiency, gc, expected):
    report, selection = single_window(efficiency, gc)
    window = selection['windows'][0]
    score, components = _efficacy_priority(
        [{'sequence': report['candidate']['reference']['sequence']}], efficiency, 1)
    assert score == expected
    assert components['sifi_is_efficient_fraction'] == efficiency
    assert components['sirna_gc_30_52_fraction'] == gc
    assert window['efficacy_priority_score'] == expected
    assert window['efficacy_priority_components'] == components
    assert window['length_penalty'] == 26.5
    assert window['adjusted_priority_score'] == max(0, expected-26.5)
    assert selection['selected_window']['efficacy_priority_score'] == expected
    assert set(components) == {'sifi_is_efficient_fraction', 'sirna_gc_30_52_fraction',
                               'unique_sirna_sequence_fraction'}
    assert '0.75' in selection['selection_method']['formula']
    assert '0.25' in selection['selection_method']['formula']


def test_diversity_does_not_change_score_but_still_breaks_ties():
    sequence = 'G'*7 + 'A'*14
    repeated = [{'sequence': sequence} for _ in range(4)]
    distinct = [{'sequence': sequence[i:]+sequence[:i]} for i in range(4)]
    a, ac = _efficacy_priority(repeated, 2, 4)
    b, bc = _efficacy_priority(distinct, 2, 4)
    assert a == b == 62.5
    assert ac['unique_sirna_sequence_fraction'] == .25
    assert bc['unique_sirna_sequence_fraction'] == 1
    assert ac['sifi_is_efficient_fraction'] == bc['sifi_is_efficient_fraction'] == .5
    assert ac['sirna_gc_30_52_fraction'] == bc['sirna_gc_30_52_fraction'] == 1
    row = dict(adjusted_priority_score=a, quality_score=90, start_1based=1,
               length_nt=180, region_id='green-001')
    assert priority_key(dict(row, efficacy_priority_components=bc)) < priority_key(
        dict(row, efficacy_priority_components=ac))


@pytest.mark.parametrize('length,penalty', [(27, 25.5), (100, 40/3), (179, 1/6),
                                           (180, 0), (200, 0), (250, 0)])
def test_length_penalty_is_unchanged(length, penalty):
    scores = length_scores(75, length)
    assert scores['length_penalty'] == pytest.approx(penalty)
    assert scores['adjusted_priority_score'] == pytest.approx(max(0, 75-penalty))
    assert length_scores(None, length)['adjusted_priority_score'] is None


def test_missing_original_metrics_remain_unscored():
    assert _efficacy_priority([{'sequence': 'A'*21}], 0, 0)[0] is None
    assert _efficacy_priority([{'sequence': 'N'*21}], 1, 1)[0] is None


def test_individual_export_and_batch_preserve_the_same_scores(tmp_path, monkeypatch):
    import rnai_designer.automatic_window_reporting as reporting
    monkeypatch.setattr(reporting, '_plot_selected', lambda *args: None)
    entries = []
    for index, (efficiency, gc, expected) in enumerate(CASES, 1):
        report, selection = single_window(efficiency, gc)
        source = tmp_path / f'input_{index}.json'
        source.write_text(json.dumps(report), encoding='utf-8')
        exported = write_selection_outputs(report, selection, source, tmp_path / f'gene_{index}')
        window = exported['selected_window']
        saved = json.loads((tmp_path / f'gene_{index}' / 'selected_window.json').read_text(encoding='utf-8'))
        with (tmp_path / f'gene_{index}' / 'all_windows.csv').open(encoding='utf-8', newline='') as handle:
            csv_row = next(csv.DictReader(handle))
        assert saved['efficacy_priority_score'] == float(csv_row['efficacy_priority_score']) == expected
        assert saved['adjusted_priority_score'] == float(csv_row['adjusted_priority_score']) == max(0, expected-26.5)
        assert saved['efficacy_priority_components']['unique_sirna_sequence_fraction'] == 1
        assert float(csv_row['efficacy_priority_components.unique_sirna_sequence_fraction']) == 1
        entry = dict(input_index=index, gene_id=f'gene{index}', status='complete', result={
            'status': 'complete', 'selected_window': window,
            'selected_window_offtarget_filter': {'passes': True, 'status': 'criterion_met'}})
        assert result_row(entry)['efficacy_priority_score'] == expected
        assert result_row(entry)['adjusted_priority_score'] == window['adjusted_priority_score']
        entries.append(entry)
    ranked = save_reports(tmp_path, {'status': 'complete', 'genes': entries})
    assert [row['input_index'] for row in ranked] == [3, 1, 2, 4]
    with (tmp_path/'ranking_genes.csv').open(encoding='utf-8', newline='') as handle:
        assert [float(row['adjusted_priority_score']) for row in csv.DictReader(handle)] == [73.5, 48.5, 0, 0]
