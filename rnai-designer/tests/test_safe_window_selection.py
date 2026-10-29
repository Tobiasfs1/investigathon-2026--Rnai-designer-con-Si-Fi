import copy
from pathlib import Path
import pytest

from rnai_designer.automatic_window_selection import select_windows
from rnai_designer.data_models import Candidate
from rnai_designer.safe_window_selection import choose_safe_window, evaluate_filter
from rnai_designer.priority_scoring import priority_key
from test_automatic_window_selection import make_report


def evaluations(status='criterion_met'):
    return [dict(organism=o,status='imported_unvalidated',criterion=dict(status=status,
                 exact_21nt_rows=1 if status=='exact_match_detected' else 0,reasons=[],alignments=[]))
            for o in ['Homo sapiens','Zea mays','Apis mellifera']]


def run(tmp_path,results,report=None):
    report=report or make_report([(1,182)])
    candidate=Candidate.from_fields('test',report['candidate']['reference']['sequence'])
    selection=select_windows(report,180,180)
    original=copy.deepcopy(list(selection['windows']))
    calls=[]
    def off(scoped,settings,out,timeout):
        calls.append((scoped.sequence,out))
        result=results[len(calls)-1]
        if isinstance(result,Exception):raise result
        return result
    result=choose_safe_window(candidate,report,selection,None,tmp_path,300,off,lambda *a:None)
    assert list(selection['windows'])==original
    return selection,result,calls,original


@pytest.mark.parametrize('failed_first', ['exact_match_detected','not_evaluated','indeterminate','error'])
def test_fallback_to_second(tmp_path,failed_first):
    selection,(ev,filt,audit),calls,original=run(tmp_path,[evaluations(failed_first),evaluations()])
    expected=sorted(original,key=priority_key)[1]
    winner=selection['selected_window']
    assert winner['window_id']==expected['window_id']
    assert winner['efficacy_priority_score']==expected['efficacy_priority_score']
    assert winner['adjusted_priority_score']==expected['adjusted_priority_score']
    assert winner['length_penalty']==expected['length_penalty']
    assert winner['sequence_release_status']=='released_after_filter'
    assert winner['offtarget_safety_filter']['passes'] is True
    assert audit['selected_priority_rank']==2 and audit['attempted_window_count']==2
    assert audit['attempts'][0]['filter']['status']==failed_first
    assert calls[0][1]!=calls[1][1]
    assert winner['sequence']==calls[1][0]


def test_third_and_early_stop(tmp_path):
    selection,(_,filt,audit),calls,_=run(tmp_path,[evaluations('exact_match_detected'),evaluations('indeterminate'),evaluations()])
    assert selection['selected_window']['selected_priority_rank']==3
    assert len(calls)==3 and filt['passes']


def test_first_pass_stops_immediately(tmp_path):
    selection,(_,filt,audit),calls,rows=run(tmp_path,[evaluations()])
    assert len(calls)==1 and audit['selected_priority_rank']==1
    assert selection['selected_window']['window_id']==min(rows,key=priority_key)['window_id']


def test_no_winner_when_all_rejected(tmp_path):
    selection,(_,filt,audit),calls,_=run(tmp_path,[evaluations('exact_match_detected')]*3)
    assert selection['selected_window'] is None
    assert audit['attempted_window_count']==3 and audit['status']=='no_window_passed'
    assert not filt['passes'] and len(calls)==3


def test_exception_does_not_authorize_and_next_is_tested(tmp_path):
    selection,(_,filt,audit),calls,_=run(tmp_path,[OSError('tool failed'),evaluations()])
    assert filt['passes'] and audit['selected_priority_rank']==2
    assert audit['attempts'][0]['filter']['status']=='error'
    assert audit['attempts'][0]['filter']['reasons']==['tool failed']


def test_missing_metrics_not_eligible(tmp_path):
    selection,(_,filt,audit),calls,_=run(tmp_path,[],make_report([(1,182)],missing=True))
    assert selection['selected_window'] is None and not calls
    assert audit['status']=='no_eligible_windows' and not filt['passes']


def test_too_short_never_passes(tmp_path):
    report=make_report([(1,20)])
    report['sifi']['green_mask']={'positions_1based':list(range(1,21)), 'sirna_size':21, 'source_function':'general_helpers.get_target_data', 'coordinate_system':'1-based positions'}
    selection,(_,filt,audit),calls,_=run(tmp_path,[],report)
    assert selection['selected_window'] is None and not calls
    assert not filt['passes']


@pytest.mark.parametrize('variant',['missing_species','duplicate_species','positive_hit','error_run','not_evaluated_run','unvalidated_criterion'])
def test_filter_requires_all_species_and_valid_criteria(variant):
    ev=evaluations()
    if variant=='missing_species':ev.pop()
    if variant=='duplicate_species':ev.append(copy.deepcopy(ev[0]))
    if variant=='positive_hit':ev[0]['criterion']['exact_21nt_rows']=1
    if variant=='error_run':ev[0]['status']='error'
    if variant=='not_evaluated_run':ev[0]['status']='not_evaluated'
    if variant=='unvalidated_criterion':ev[0]['criterion']['status']='indeterminate'
    assert not evaluate_filter(ev)['passes']


def test_ties_use_existing_order_and_short_window_can_win(tmp_path):
    report=make_report([(1,27),(50,230)])
    selection,(_,filt,audit),calls,_=run(tmp_path,[evaluations('exact_match_detected'),evaluations('exact_match_detected'),evaluations()],report)
    # The two 180-nt windows outrank the 27-nt region under the existing penalty.
    assert selection['selected_window']['window_class']=='short_green_region'
    assert selection['selected_window']['length_nt']==27
    assert audit['selected_priority_rank']==3


def test_empty_green_regions(tmp_path):
    report=make_report([])
    report['sifi']['green_mask']={'positions_1based':[], 'sirna_size':21, 'source_function':'general_helpers.get_target_data', 'coordinate_system':'1-based positions'}
    selection,(_,filt,audit),calls,_=run(tmp_path,[],report)
    assert selection['selected_window'] is None and not calls


@pytest.mark.parametrize('passes_second',[True,False])
def test_pipeline_exports_only_approved_winner(tmp_path,monkeypatch,passes_second):
    import csv
    import json
    import rnai_designer.pipeline_workflow as pipeline
    report=make_report([(1,182)])
    candidate=Candidate(report['candidate'],report['candidate']['gene_id'],report['candidate']['reference']['sequence'])
    calls=[]
    def off(candidate,settings,out,timeout):
        if out.name=='offtarget_runs': ev=evaluations('not_evaluated')
        else:
            calls.append(out)
            ev=evaluations() if passes_second and len(calls)==2 else evaluations('exact_match_detected')
        for item in ev:item['criterion']['label']=item['criterion']['status']
        return ev
    def write_reports(out,*args):
        report['offtarget_criterion']={'status':'not_evaluated'}
        out.mkdir(parents=True)
        (out/'report.json').write_text(json.dumps(report),encoding='utf-8')
        return report
    monkeypatch.setattr(pipeline,'run_offtargets',off)
    monkeypatch.setattr(pipeline,'write_reports',write_reports)
    result=pipeline.run_integrated_workflow(candidate,dict(records=report['sifi']['records']),None,tmp_path,max_length=180)
    summary=result['summary']
    persisted=json.loads((tmp_path/'pipeline_result.json').read_text(encoding='utf-8'))
    selection=json.loads((tmp_path/'window_selection'/'selection_report.json').read_text(encoding='utf-8'))
    assert persisted==summary
    assert summary['offtarget_selection']['attempted_window_count']==(2 if passes_second else 3)
    assert len(list(csv.DictReader((tmp_path/'window_selection'/'all_windows.csv').open())))==3
    if passes_second:
        assert summary['selected_window']['selected_priority_rank']==2
        assert summary['selected_window_offtarget_filter']['passes'] is True
        assert selection['selection_reasons']['best_score_tie_count']==3
        assert (tmp_path/'window_selection'/'selected_window.json').exists()
    else:
        assert summary['selected_window'] is None
        assert not (tmp_path/'window_selection'/'selected_window.json').exists()
        assert selection['selected_window_id'] is None
