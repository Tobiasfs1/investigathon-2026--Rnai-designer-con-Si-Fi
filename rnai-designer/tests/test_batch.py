import copy
import csv
import json
from pathlib import Path
from types import SimpleNamespace
import pytest

import rnai_designer.__main__ as cli
import rnai_designer.simple_run as simple
from rnai_designer.batch import batch_command, parse_fasta, result_row


@pytest.fixture
def setup(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(simple, 'MAX_PATH_CHARS', 2000)
    return tmp_path


def write_fasta(root, names):
    p=root/'genes.fasta'
    p.write_text(''.join(f'>{name}\nacgtacgt\n' for name in names),encoding='utf-8')
    return p


def summary(gene, score=75):
    window=dict(start_1based=1,end_1based_inclusive=180,length_nt=180,
                window_class='conventional_180_250nt',efficacy_priority_score=score,
                adjusted_priority_score=score,length_penalty=0,quality_score=92,gc_percent=45,
                is_efficient_true_count=80,is_efficient_evaluable_count=100,
                sirna_unique_count=100,sirna_unique_sequences_count=90,
                sequence_release_status='withheld_filter_failed_or_indeterminate')
    return dict(candidate_id=gene,status='complete',main_target='ref',selected_window=window,
                selected_window_offtarget_filter=dict(status='criterion_met',passes=True,
                by_species=[dict(organism='Homo sapiens',status='imported_unvalidated')]))


def fake_workflow(calls, failures=(), interrupt=None):
    def run(args,**kwargs):
        gene=args[args.index('--gene-id')+1];calls.append(gene)
        out=Path(args[args.index('--out')+1])
        assert not list(out.iterdir())
        assert '--run-sifi-headless' in args
        assert kwargs['input_metadata']['original_gene_name']==gene
        print('extensive report and repeated warning')
        if gene==interrupt: raise KeyboardInterrupt()
        if gene in failures: raise ValueError('synthetic error')
        payload=summary(gene,60+len(calls))
        (out/'pipeline_result.json').write_text(json.dumps(payload),encoding='utf-8')
        return 0
    return run


def manifest(root):
    batch=next((root/'output').glob('batch_*'))
    return batch,json.loads((batch/'resumen_batch.json').read_text(encoding='utf-8'))


def test_twenty_genes_and_export(setup,capsys):
    source=write_fasta(setup,[f'Gene{i}' for i in range(20)])
    calls=[]
    assert batch_command([str(source)],fake_workflow(calls))==0
    root,state=manifest(setup)
    assert len(calls)==20 and state['counts']==dict(total=20,completed=20,failed=0,evaluable=20)
    results=list(csv.DictReader((root/'resultados_genes.csv').open(encoding='utf-8')))
    ranked=list(csv.DictReader((root/'ranking_genes.csv').open(encoding='utf-8')))
    assert len(results)==len(ranked)==20
    assert [float(r['adjusted_priority_score']) for r in ranked]==sorted([float(r['adjusted_priority_score']) for r in results],reverse=True)
    assert results[0]['efficient_percent']=='80.0' and results[0]['diversity_percent']=='90.0'
    assert results[0]['human_offtarget']=='imported_unvalidated'
    assert results[0]['maize_offtarget']=='not_evaluated'
    assert all((root/e['output_directory']/'pipeline_result.json').exists() for e in state['genes'])
    text=capsys.readouterr().out
    assert 'Analizando gen 4/20: Gene3' in text and 'extensive report' not in text
    assert 'heuristico' in text and 'conventional_180_250nt' in text


def test_partial_duplicate_names_invalid_and_no_score(setup):
    source=setup/'genes.fasta';source.write_text('>same\nACGT\n>same\nACGT\n>bad\nACGT\n>invalid\nXYZ!\n>no_score\nACGT\n')
    calls=[];base=fake_workflow(calls,failures={'bad'})
    def run(args,**kwargs):
        code=base(args,**kwargs)
        if args[1]=='no_score':
            out=Path(args[args.index('--out')+1]);s=summary('no_score');s['selected_window']=None
            (out/'pipeline_result.json').write_text(json.dumps(s))
        return code
    assert batch_command([str(source)],run)==2
    root,s=manifest(setup)
    assert calls.count('same')==2
    assert len(set(e.get('output_directory') for e in s['genes'] if e.get('output_directory')))==4
    assert s['counts']==dict(total=5,completed=3,failed=2,evaluable=2)
    rows=list(csv.DictReader((root/'ranking_genes.csv').open()))
    assert all(r['adjusted_priority_score']=='' and r['rank']=='' for r in rows[2:])
    assert any(e.get('error')=='synthetic error' for e in s['genes'])
    assert s['genes'][0]['gene_id']==s['genes'][1]['gene_id']=='same'


def test_resume_skips_completed_and_retries_failed_without_overwrite(setup):
    source=write_fasta(setup,['good','bad']);calls=[]
    assert batch_command([str(source)],fake_workflow(calls,failures={'bad'}))==2
    root,state=manifest(setup);good=root/state['genes'][0]['output_directory']/'pipeline_result.json';raw=good.read_bytes()
    old_bad=state['genes'][1]['output_directory']
    calls=[]
    assert batch_command(['--resume',str(root)],fake_workflow(calls))==0
    state=json.loads((root/'resumen_batch.json').read_text())
    assert calls==['bad'] and good.read_bytes()==raw
    assert state['genes'][1]['output_directory']!=old_bad
    assert len(state['genes'][1]['attempts'])==2
    assert (root/old_bad).exists()
    assert batch_command(['--resume',str(root)],lambda *a,**k:pytest.fail('completed genes repeated'))==0


def test_interrupt_and_resume(setup):
    source=write_fasta(setup,['first','second','third']);calls=[]
    assert batch_command([str(source)],fake_workflow(calls,interrupt='second'))==130
    root,state=manifest(setup)
    assert state['status']=='interrupted' and state['genes'][1]['status']=='interrupted'
    calls=[]
    assert batch_command(['--resume',str(root)],fake_workflow(calls))==0
    assert calls==['second','third']


def test_recover_report_written_before_interruption(setup):
    source=write_fasta(setup,['gene']);calls=[];base=fake_workflow(calls)
    def interrupted(*args,**kwargs):
        base(*args,**kwargs);raise KeyboardInterrupt()
    assert batch_command([str(source)],interrupted)==130
    root,_=manifest(setup)
    assert batch_command(['--resume',str(root)],lambda *a,**k:pytest.fail('repeated completed analysis'))==0


def test_source_mismatch_and_new_batch_collision(setup):
    source=write_fasta(setup,['a']);calls=[]
    assert batch_command([str(source)],fake_workflow(calls))==0
    root,_=manifest(setup)
    source.write_text('>different\nACGT\n')
    assert batch_command([str(source),'--resume',str(root)],fake_workflow(calls))==2
    assert batch_command([str(source)],fake_workflow(calls))==0
    assert len(list((setup/'output').glob('batch_*')))==2


@pytest.mark.parametrize('text',['','ACGT\n>name\nACGT'])
def test_invalid_structure_before_creation(setup,text):
    source=setup/'genes.fasta';source.write_text(text)
    assert batch_command([str(source)],lambda *a,**k:pytest.fail('not valid'))==2
    assert not (setup/'output').exists()


def test_multi_fasta_lowercase_and_whitespace(setup):
    source=setup/'genes.fasta';source.write_text('>full gene name\nac gt\n tgca\n>other\nACGT\n')
    rows=parse_fasta(source)
    assert rows[0]['gene_id']=='full gene name' and rows[0]['sequence']=='ACGTTGCA'


def test_individual_and_batch_same_workflow_result(setup,monkeypatch):
    db=setup/'db';db.mkdir()
    monkeypatch.setattr(cli.Settings,'load',lambda *a:SimpleNamespace(database_dir=db))
    monkeypatch.setattr(cli,'_legacy_sifi_context',lambda args,*a:setattr(args,'_sifi_runtime',{}))
    monkeypatch.setattr(cli,'_list_legacy_targets',lambda *a:{'available_targets':[{'hit_name':'ref'}]})
    calls=[]
    def analysis(args,candidate,settings,out,targets):
        calls.append((candidate.gene_id,candidate.sequence,targets,args.min_window,args.max_window,args.length_penalty_weight))
        report=summary(candidate.gene_id)
        report.update(completed_reference_count=1,failed_reference_count=0)
        (out/'pipeline_result.json').write_text(json.dumps(report))
        return report
    monkeypatch.setattr(cli,'_run_multi_reference_workflow',analysis)
    source=write_fasta(setup,['gene'])
    assert cli.main(['batch',str(source)])==0
    assert cli.main(['gene','acgtacgt'])==0
    root,state=manifest(setup)
    individual=next(p for p in (setup/'output').glob('gene_*/pipeline_result.json'))
    assert state['genes'][0]['result']==json.loads(individual.read_text())
    assert calls[0]==calls[1]


def test_no_scores_invented():
    s=summary('missing');s['selected_window'].pop('adjusted_priority_score')
    row=result_row(dict(input_index=1,gene_id='missing',status='complete',result=s))
    assert row['adjusted_priority_score'] is None and row['ranking_status']=='no_evaluable_result'


def test_zero_scores_ties_short_classes_and_verification(setup):
    source=write_fasta(setup,['first','second','missing'])
    def run(args,**kwargs):
        gene=args[1];out=Path(args[args.index('--out')+1]);s=summary(gene,0)
        if gene=='second':
            s['selected_window'].update(length_nt=27,window_class='short_green_region')
        if gene=='missing':s['selected_window']=None
        (out/'pipeline_result.json').write_text(json.dumps(s))
        return 0
    assert batch_command([str(source)],run)==0
    root,s=manifest(setup)
    rows=list(csv.DictReader((root/'ranking_genes.csv').open()))
    assert [r['gene_id'] for r in rows]==['first','second','missing']
    assert rows[0]['rank']=='1' and rows[1]['rank']=='2' and rows[2]['rank']==''
    assert rows[1]['window_class']=='short_green_region'
    assert rows[1]['adjusted_priority_score']=='0'
    assert rows[1]['offtarget_filter_status']=='criterion_met'


def test_no_report_and_return_code_failure(setup):
    source=write_fasta(setup,['missing'])
    assert batch_command([str(source)],lambda *a,**k:2)==2
    _,state=manifest(setup)
    assert state['genes'][0]['status']=='failed' and state['counts']['evaluable']==0


def test_path_check_precedes_creation(setup,monkeypatch):
    monkeypatch.setattr(simple,'MAX_PATH_CHARS',10)
    source=write_fasta(setup,['gene'])
    assert batch_command([str(source)],lambda *a,**k:pytest.fail('long path'))==2
    assert not (setup/'output').exists()


def test_changed_completed_report_is_not_reused(setup):
    source=write_fasta(setup,['gene']);calls=[]
    assert batch_command([str(source)],fake_workflow(calls))==0
    root,state=manifest(setup)
    old=root/state['genes'][0]['output_directory']/'pipeline_result.json'
    old.write_text(json.dumps(summary('modified',1)))
    calls=[]
    assert batch_command(['--resume',str(root)],fake_workflow(calls))==0
    assert calls==['gene']
    assert json.loads(old.read_text())['candidate_id']=='modified'


def test_concurrent_resume_lock(setup):
    from rnai_designer.batch import batch_lock
    root=setup/'lock';root.mkdir()
    with batch_lock(root):
        with pytest.raises(ValueError,match='ya esta en ejecucion'):
            with batch_lock(root):pytest.fail('lock not exclusive')


@pytest.mark.parametrize('error',[RuntimeError('engine failed'),SystemExit(2)])
def test_unexpected_engine_failure_is_isolated(setup,error):
    source=write_fasta(setup,['bad','good']);calls=[];base=fake_workflow(calls)
    def run(args,**kwargs):
        if args[1]=='bad':raise error
        return base(args,**kwargs)
    assert batch_command([str(source)],run)==2
    _,state=manifest(setup)
    assert calls==['good'] and state['counts']['completed']==1
    assert state['genes'][0]['error_type']==type(error).__name__


def test_old_failed_winner_cannot_enter_batch_ranking():
    report=summary('old')
    report['selected_window_offtarget_filter'].update(status='indeterminate',passes=False)
    row=result_row(dict(input_index=1,gene_id='old',status='complete',result=report))
    assert row['ranking_status']=='no_evaluable_result'
    assert row['adjusted_priority_score'] is None
