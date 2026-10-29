import copy
import io
import pytest
from rnai_designer.priority_scoring import length_scores, priority_key
from rnai_designer.console_output import print_workflow_summary
from rnai_designer.automatic_window_selection import select_windows
from test_automatic_window_selection import make_report

@pytest.mark.parametrize('length', [27,100,179,180,200,250])
def test_lengths_preserve_original_science(length):
    report = make_report([(1,length)])
    before = copy.deepcopy(report)
    result = select_windows(report)
    baseline = select_windows(report, length_penalty_weight=0)
    for row, old in zip(result['windows'], baseline['windows']):
        assert row['length_penalty'] == pytest.approx(30*max(0,1-length/180))
        assert row['adjusted_priority_score'] == pytest.approx(max(0,row['efficacy_priority_score']-row['length_penalty']))
        assert {k:v for k,v in row.items() if k not in ('length_penalty','adjusted_priority_score')} == {k:v for k,v in old.items() if k not in ('length_penalty','adjusted_priority_score')}
    assert report == before

def test_missing_and_bicc_example():
    assert length_scores(None,27) == {'length_penalty':25.5,'adjusted_priority_score':None}
    assert length_scores(66.6667,27)['adjusted_priority_score'] == pytest.approx(41.1667)
    assert select_windows(make_report([(1,27)],missing=True))['selected_window'] is None

@pytest.mark.parametrize('weight',[-1,float('nan'),float('inf')])
def test_invalid_weight(weight):
    with pytest.raises(ValueError): length_scores(50,100,weight)

@pytest.mark.parametrize('status',['complete','partial','error','no_reference_completed'])
def test_console_old_reports_and_failed_references(status):
    summary={'candidate_id':'test','status':status,'references':[
        {'status':'error','main_target':'failed','error':'incomplete'},
        {'status':'complete','main_target':'ok','selected_window':None}]}
    original=copy.deepcopy(summary)
    stream=io.StringIO()
    print_workflow_summary(summary,'.',stream=stream,verbose=True)
    assert 'Sin ganador' in stream.getvalue()
    assert 'incomplete' in stream.getvalue()
    assert summary == original

def test_console_tie_and_old_score():
    w={'efficacy_priority_score':70,'start_1based':1,'end_1based_inclusive':27,'length_nt':27}
    summary={'references':[{'main_target':r,'status':'complete','selected_window':dict(w)} for r in ['B','A']]}
    stream=io.StringIO()
    print_workflow_summary(summary,'.',sequence='A'*27,stream=stream)
    assert any(line.startswith('Referencia ganadora') and line.endswith(': A') for line in stream.getvalue().splitlines())
    assert 'A'*27 not in stream.getvalue()  # Missing authorization must retain sequence.
    assert priority_key(w) == priority_key(dict(w))


def test_existing_bicc_regression():
    import csv
    from pathlib import Path
    rows=[]
    for source in (Path(__file__).resolve().parents[1]/'output'/'BicC_Dm').glob('references/*/window_selection/all_windows.csv'):
        for raw in csv.DictReader(source.open(encoding='utf-8-sig',newline='')):
            score=float(raw['efficacy_priority_score']) if raw.get('efficacy_priority_score') else None
            row={**raw,'length_nt':int(raw['length_nt']), 'efficacy_priority_score':score}
            row.update(length_scores(score,row['length_nt']))
            rows.append(row)
    if not rows: pytest.skip('Existing BicC_Dm report unavailable')
    short=[r for r in rows if r['length_nt']==27 and r['efficacy_priority_score'] is not None]
    assert any(r['length_penalty']==25.5 and r['adjusted_priority_score']==pytest.approx(41.1666666667) for r in short)
    winner=max((r for r in rows if r['adjusted_priority_score'] is not None),key=lambda r:r['adjusted_priority_score'])
    assert winner['length_nt']==166
    assert int(winner['start_1based'])==1170
    assert winner['adjusted_priority_score']==pytest.approx(59.99543378995434)
