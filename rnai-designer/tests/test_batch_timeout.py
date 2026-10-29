"""Simulated Python workers only; never run si-Fi, Bowtie, RNAplfold or real genes."""
import csv
import io
import json
import os
from pathlib import Path
import sys
from time import monotonic

import pytest

import rnai_designer.batch as batch
import rnai_designer.batch_process as supervision
from rnai_designer.batch_process import run_worker, completed_worker
from rnai_designer.atomic_io import write_json

SIMULATOR = r'''
import json,os,subprocess,sys,time
from pathlib import Path
bootstrap=json.loads(sys.stdin.readline())
r=json.loads(Path(bootstrap['request_file']).read_text(encoding='utf-8'))
r.update(started=bootstrap['started'],deadline=bootstrap.get('deadline'))
out=Path(r['argv'][r['argv'].index('--out')+1])
gene=r['metadata']['original_gene_name']
Path(r['progress']).write_text(json.dumps({'stage':'evaluando_ventanas','reference_id':'GITV_TEST'}))
(out/'temporary_path.txt').write_text(os.environ['TMP'])
w={'start_1based':1,'end_1based_inclusive':21,'length_nt':21,'sequence':'A'*21,
   'window_class':'short_green_region','efficacy_priority_score':75,'length_penalty':26.5,
   'adjusted_priority_score':48.5,'gc_percent':0,'quality_score':100,
   'is_efficient_true_count':1,'is_efficient_evaluable_count':1,
   'sirna_unique_count':1,'sirna_unique_sequences_count':1}
s={'candidate_id':gene,'status':'complete','selected_window':w,
   'selected_window_offtarget_filter':{'status':'criterion_met','passes':True}}
(out/'pipeline_result.json').write_text(json.dumps(s))
if gene in ('slow','cpu') and os.environ.get('SIM_FAST')!='1':
    grand="import time;time.sleep(60)"
    child="import subprocess,sys,time;from pathlib import Path;p=subprocess.Popen([sys.executable,'-B','-c',%r],creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0));Path(sys.argv[1]).write_text(str(p.pid));time.sleep(60)" % grand
    process=subprocess.Popen([sys.executable,'-B','-c',child,str(out/'grandchild.pid')],creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
    (out/'child.pid').write_text(str(process.pid))
    if gene=='cpu':
        while True:pass
    time.sleep(60)
Path(r['receipt']).write_text(json.dumps({'return_code':0,'within_deadline':r['deadline'] is None or time.monotonic()<=r['deadline']}))
'''


@pytest.fixture
def simulator(tmp_path,monkeypatch):
    import rnai_designer.simple_run as simple
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(simple,'MAX_PATH_CHARS',2000)
    (tmp_path/'fake_batch_worker.py').write_text(SIMULATOR,encoding='utf-8')
    monkeypatch.setattr(supervision,'_WORKER_MODULE','fake_batch_worker')
    monkeypatch.setenv('PYTHONPATH',str(tmp_path)+os.pathsep+os.environ.get('PYTHONPATH',''))
    return tmp_path


def fasta(root,names):
    p=root/'genes.fasta'
    p.write_text(''.join('>'+name+'\n'+'A'*21+'\n' for name in names))
    return p


def state(root):
    folder=next((root/'output').glob('batch_*'))
    return folder,json.loads((folder/'resumen_batch.json').read_text(encoding='utf-8'))


def must_not_run_in_parent(*a,**k):
    pytest.fail('Workflow must execute in isolated worker')


def running(pid):
    if os.name=='nt':
        import ctypes as c
        from ctypes import wintypes as w
        api=c.WinDLL('kernel32',use_last_error=True)
        api.OpenProcess.argtypes=[w.DWORD,w.BOOL,w.DWORD];api.OpenProcess.restype=w.HANDLE
        api.GetExitCodeProcess.argtypes=[w.HANDLE,c.POINTER(w.DWORD)];api.GetExitCodeProcess.restype=w.BOOL
        api.CloseHandle.argtypes=[w.HANDLE]
        handle=api.OpenProcess(0x1000,False,pid)
        if not handle:return False
        try:
            code=w.DWORD()
            return bool(api.GetExitCodeProcess(handle,c.byref(code))) and code.value==259
        finally:api.CloseHandle(handle)
    stat=Path(f'/proc/{pid}/stat')
    if stat.exists() and stat.read_text().split()[2]=='Z':return False
    try:os.kill(pid,0);return True
    except ProcessLookupError:return False


def test_worker_finishes_normally(tmp_path):
    code="import sys;sys.stdin.readline();print('simulated success')"
    result=run_worker([sys.executable,'-B','-c',code],{},tmp_path/'normal.log',5)
    assert result.returncode==0 and not result.timed_out and result.elapsed_s<5
    assert result.cleanup_confirmed and not running(result.pid)


@pytest.mark.parametrize('gene',['slow','cpu'])
def test_timeout_cleans_worker_children_and_private_temp(simulator,gene):
    source=fasta(simulator,[gene])
    assert batch.batch_command([str(source),'--gene-timeout','0.8'],must_not_run_in_parent)==2
    root,s=state(simulator);entry=s['genes'][0];attempt=entry['attempts'][0]
    assert entry['status']=='incomplete_timeout'
    assert entry['result'] is None and entry['result_sha256'] is None
    assert entry['elapsed_s']>=.8 and entry['timeout_seconds']==.8
    assert entry['timeout_stage']=='evaluando_ventanas' and entry['timeout_reference_id']=='GITV_TEST'
    assert attempt['cleanup_confirmed'] is True
    assert not running(attempt['worker_pid'])
    out=root/entry['output_directory']
    for file in ['child.pid','grandchild.pid']:
        assert not running(int((out/file).read_text()))
    assert not Path((out/'temporary_path.txt').read_text()).exists()
    # A report published just before cancellation must not authorize completion.
    assert json.loads((out/'pipeline_result.json').read_text())['status']=='complete'
    assert not completed_worker(root,attempt)
    for name in ['resultados_genes.csv','ranking_genes.csv']:
        row=next(csv.DictReader((root/name).open(encoding='utf-8')))
        assert row['status']=='incomplete_timeout' and row['adjusted_priority_score']==''
    stale={**entry,'result':json.loads((out/'pipeline_result.json').read_text())}
    assert batch.result_row(stale)['ranking_status']=='no_evaluable_result'


def test_next_gene_continues_and_resume_retries_only_timeout(simulator,monkeypatch):
    source=fasta(simulator,['slow','normal'])
    assert batch.batch_command([str(source),'--gene-timeout','0.7'],must_not_run_in_parent)==2
    root,s=state(simulator)
    assert [e['status'] for e in s['genes']]==['incomplete_timeout','complete']
    normal=root/s['genes'][1]['output_directory']/'pipeline_result.json';original=normal.read_bytes()
    old_timeout=s['genes'][0]['output_directory']
    monkeypatch.setenv('SIM_FAST','1')
    assert batch.batch_command(['--resume',str(root),'--gene-timeout','5'],must_not_run_in_parent)==0
    s=json.loads((root/'resumen_batch.json').read_text())
    assert [e['status'] for e in s['genes']]==['complete','complete']
    assert s['genes'][0]['output_directory']!=old_timeout and (root/old_timeout).is_dir()
    assert len(s['genes'][0]['attempts'])==2 and len(s['genes'][1]['attempts'])==1
    assert normal.read_bytes()==original and s['gene_timeout_seconds']==5
    assert batch.batch_command(['--resume',str(root)],must_not_run_in_parent)==0
    assert normal.read_bytes()==original


@pytest.mark.parametrize('options,expected',[([],360),(['--gene-timeout','300'],300),(['--gene-timeout','0'],0)])
def test_default_custom_and_disabled_timeout(simulator,options,expected):
    assert batch.batch_command([str(fasta(simulator,['normal'])),*options],must_not_run_in_parent)==0
    root,s=state(simulator)
    assert s['gene_timeout_seconds']==expected and s['genes'][0]['timeout_seconds']==expected
    assert s['genes'][0]['status']=='complete'
    assert s['genes'][0]['elapsed_s']>0


def test_disabled_supervision_allows_work_to_finish(tmp_path):
    code="import sys,time;sys.stdin.readline();time.sleep(.15)"
    r=run_worker([sys.executable,'-B','-c',code],{},tmp_path/'disabled.log',0)
    assert not r.timed_out and r.elapsed_s>=.15 and r.returncode==0


@pytest.mark.parametrize('value',['-1','nan','inf'])
def test_invalid_timeout_creates_no_batch(simulator,value):
    with pytest.raises(SystemExit):
        batch.batch_command([str(fasta(simulator,['normal'])),'--gene-timeout',value],must_not_run_in_parent)
    assert not (simulator/'output').exists()


def test_preparation_time_counts_toward_deadline(tmp_path):
    code="import sys;sys.stdin.readline();raise Exception('must not start')"
    result=run_worker([sys.executable,'-B','-c',code],{},tmp_path/'expired.log',.1,started=monotonic()-1)
    assert result.timed_out and result.returncode==124


def test_unfinished_receipt_never_recovers_complete_report(tmp_path):
    attempt={'worker_state':'x.worker.json','worker_receipt':'x.done.json'}
    write_json(tmp_path/attempt['worker_state'],{'status':'running'})
    assert not completed_worker(tmp_path,attempt)
    write_json(tmp_path/attempt['worker_receipt'],{'return_code':0,'within_deadline':False})
    assert not completed_worker(tmp_path,attempt)
    write_json(tmp_path/attempt['worker_receipt'],{'return_code':0,'within_deadline':True})
    assert completed_worker(tmp_path,attempt)
    write_json(tmp_path/attempt['worker_state'],{'status':'incomplete_timeout'})
    assert not completed_worker(tmp_path,attempt)


def test_worker_context_tracks_existing_progress_without_changing_workflow(tmp_path,monkeypatch):
    from rnai_designer.batch_worker import main
    request={'argv':[],'metadata':{},'progress':str(tmp_path/'progress.json'),
             'receipt':str(tmp_path/'receipt.json'),'started':monotonic(),'deadline':None}
    monkeypatch.setattr(sys,'stdin',io.StringIO(json.dumps(request)+'\n'))
    def workflow(argv,**kwargs):
        print('[si-Fi] GITV0001: ejecutando Conda',file=sys.stderr)
        print('[rnai-designer] Filtro off-target: ventana 2/5',file=sys.stderr)
        return 0
    assert main(workflow)==0
    context=json.loads((tmp_path/'progress.json').read_text())
    assert context=={'reference_id':'GITV0001','stage':'Filtro off-target: ventana 2/5'}
    assert json.loads((tmp_path/'receipt.json').read_text())['return_code']==0


def test_normal_worker_does_not_leave_descendants(tmp_path):
    marker=tmp_path/'child.pid'
    code=("import json,sys,subprocess;from pathlib import Path;r=json.loads(sys.stdin.readline());"
          "p=subprocess.Popen([sys.executable,'-B','-c','import time;time.sleep(60)'],"
          "creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0));Path(r['marker']).write_text(str(p.pid))")
    result=run_worker([sys.executable,'-B','-c',code],{'marker':str(marker)},tmp_path/'normal_child.log',5)
    assert result.returncode==0 and not result.timed_out
    assert not running(int(marker.read_text()))


def test_late_receipt_is_timeout_even_if_process_exit_is_observed_as_zero(tmp_path,monkeypatch):
    def simulated(command,request,log_path,timeout,**kwargs):
        payload=json.loads(Path(request['request_file']).read_text())
        write_json(payload['receipt'],{'return_code':0,'within_deadline':False})
        return supervision.WorkerResult(0,1)
    monkeypatch.setattr(supervision,'run_worker',simulated)
    result=supervision.run_gene_worker(None,[],{},tmp_path/'late.log',.5)
    assert result.timed_out and result.returncode==124
    assert json.loads((tmp_path/'late.worker.json').read_text())['status']=='incomplete_timeout'


def test_worker_sidecars_are_not_overwritten(tmp_path):
    state=tmp_path/'same.worker.json';state.write_text('previous evidence')
    with pytest.raises(FileExistsError):
        supervision.run_gene_worker(None,[],{},tmp_path/'same.log',1)
    assert state.read_text()=='previous evidence'


def test_timeout_marker_is_durable_before_tree_cleanup(tmp_path,monkeypatch):
    from rnai_designer.process_tree import ProcessTree
    state=tmp_path/'supervision.json';write_json(state,{'status':'running'})
    observed=[];original=ProcessTree.stop
    def stop(self):
        observed.append(json.loads(state.read_text())['status'])
        return original(self)
    monkeypatch.setattr(ProcessTree,'stop',stop)
    code="import sys,time;sys.stdin.readline();time.sleep(60)"
    result=run_worker([sys.executable,'-B','-c',code],{'supervision_state':str(state)},tmp_path/'durable.log',.15)
    assert result.timed_out and observed==['incomplete_timeout']
