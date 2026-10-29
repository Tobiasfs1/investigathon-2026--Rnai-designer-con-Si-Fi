import copy
import io
import json
from pathlib import Path
import pytest
from rnai_designer.console_output import print_workflow_summary


def fixture():
    w = dict(efficacy_priority_score=70, adjusted_priority_score=65, length_penalty=5,
             quality_score=92, start_1based=1, end_1based_inclusive=166, length_nt=166,
             window_class="short_green_region", sequence="ACGT"*41+"AA",
             is_efficient_true_count=19, is_efficient_evaluable_count=146,
             sirna_unique_sequences_count=100, sirna_unique_count=146, gc_percent=42,
             sequence_release_status="released_after_filter")
    filt = dict(passes=True,status="criterion_met",by_species=[
        dict(organism=o,status=s,reasons=[]) for o,s in [
        ("Homo sapiens","criterion_met"),("Zea mays","not_evaluated"),("Apis mellifera","imported_unvalidated")]])
    return dict(candidate_id="test",status="complete",reference_count=2,combined_window_count=100,
                references=[dict(status="complete",main_target="winner",selected_window=w,
                                 selected_window_offtarget_filter=filt,warnings=["duplicate"]*20),
                            dict(status="error",main_target="failed",error="technical error")])


def render(s, **kwargs):
    stream=io.StringIO()
    print_workflow_summary(s,"out",stream=stream,**kwargs)
    return stream.getvalue()


def test_compact_format_and_no_mutation():
    s=fixture();before=copy.deepcopy(s);text=render(s)
    seq=s['references'][0]['selected_window']['sequence']
    lines=text.splitlines();sequence_lines=[line for line in lines if line and set(line)<=set('ACGT')]
    assert ''.join(sequence_lines)==seq
    assert [len(line) for line in sequence_lines]==[60,60,46]
    assert 20<=len(lines)-len(sequence_lines)<=25
    assert '19 / 146 (13.01%) evaluables' in text
    assert '100 / 146 (68.49%)' in text
    assert 'no evaluado' in text and 'importado sin validar' in text
    assert text.count('Advertencia:')==2
    assert '180 nt' in text and 'eficacia experimental' in text
    assert 'technical error' not in text and '[]' not in text
    assert '.json' not in text and '.csv' not in text
    assert 'failed' not in text and 'Carpeta de salida' in text
    assert s==before


@pytest.mark.parametrize('verbose',[False,True])
@pytest.mark.parametrize('change', ['missing','failed','withheld','unvalidated'])
def test_sequence_authorization(verbose,change):
    s=fixture();r=s['references'][0];w=r['selected_window']
    if change=='missing':r.pop('selected_window_offtarget_filter')
    if change=='failed':r['selected_window_offtarget_filter']['passes']=False
    if change=='withheld':w['sequence_release_status']='withheld_filter_failed_or_indeterminate'
    if change=='unvalidated':r['selected_window_offtarget_filter']['status']='imported_unvalidated'
    text=render(s,verbose=verbose,sequence='A'*500)
    assert w['sequence'][:60] not in text
    assert 'A'*60 not in text


@pytest.mark.parametrize('status',['partial','error','no_reference_completed'])
def test_no_eligible_and_incomplete(status):
    s=fixture();s['status']=status
    for r in s['references']:r['selected_window']=None
    text=render(s)
    assert 'no hubo resultado elegible' in text
    assert 'incompleta' in text
    assert text.count('Advertencia:')<=2


def test_verbose_retains_details():
    text=render(fixture(),verbose=True)
    assert 'technical error' in text and 'selection_report.json' in text


def test_bicc_existing_output(tmp_path):
    root=Path(__file__).resolve().parents[1]/'output'
    source=root/'BicC_Dm4'/'pipeline_result.json'
    if not source.exists():pytest.skip('BicC_Dm4 unavailable')
    raw=source.read_bytes();s=json.loads(raw)
    text=render(s)
    assert 'GITV01000968_1' in text and '59.995' not in text
    assert '60.00 /100' in text and '1170 - 1335' in text
    assert source.read_bytes()==raw
