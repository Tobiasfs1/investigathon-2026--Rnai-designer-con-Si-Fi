"""Small operation-count/differential checks; no real engines or gene batches."""
import copy
import json
from pathlib import Path
import sys
from types import ModuleType, SimpleNamespace
import textwrap

import pytest
from rnai_designer.legacy_bridge import run_sifi_headless as bridge
from rnai_designer.legacy_session import shared_sifi_session, shared_directory
from rnai_designer.automatic_window_selection import select_windows
from rnai_designer.priority_scoring import priority_key
from test_automatic_window_selection import make_report

ROOT=Path(__file__).resolve().parents[1]
LEGACY=ROOT.parent/'siFi21--1.2.3'


def original_function(filename, name, namespace):
    text=(LEGACY/filename).read_text(encoding='utf-8')
    start=text.index('    def '+name+'(') if name=='data_to_json' else text.index('def '+name+'(')
    following=text.find('\n    def ',start+8) if name=='data_to_json' else text.find('\ndef ',start+4)
    code=textwrap.dedent(text[start:following if following!=-1 else len(text)])
    # Python 2's print has no effect on return values; suppress only that output.
    code=code.replace("print 'ok'",'pass')
    exec(code,namespace)
    return namespace[name]


def test_bridge_shared_results_equal_original_per_target(tmp_path,monkeypatch):
    counts=dict(bowtie=0,rnaplfold=0,efficiency=0)
    # Execute the original data_to_json and original get_target_data on tiny data.
    convert=original_function('sifi_pipeline.py','data_to_json',{'np':SimpleNamespace(float=float)})
    helper=original_function('general_helpers.py','get_target_data',{
        'prepare_json_data':lambda p:json.loads(Path(p).read_text()),
        'range':lambda *a:list(range(*a))})
    class Array:
        def __getitem__(self,key):return self
        def astype(self,kind):return self
        def tolist(self):return [0.25]*10
    class Pipeline:
        data_to_json=convert
        def __init__(self,*args):
            self.mode=0;self.accessibility_window=args[6];self.directory=Path(args[14])
        def create_single_fasta_file(self,*args):return 'unused'
        def create_sirnas(self,sequence,size):
            tab=self.directory/'tab.txt';tab.write_text('sirna1\t'+sequence[:size]+'\nsirna2\t'+sequence[1:size+1]+'\n')
            return 'unused',str(tab)
        def bowtie_to_lst(self,lines):return [line.strip().split('\t') for line in lines]
        def run_rnaplfold(self,*args):counts['rnaplfold']+=1;return Array()
        def calculate_efficiency(self,*args):counts['efficiency']+=1;return (True,True,False,True,False)
    rows=[['sirna1','+','A','0','A'*21,'x','0','0'],
          ['sirna1','+','B','0','A'*21,'x','0','0'],
          ['sirna2','+','B','1','A'*21,'x','0','0'],
          ['sirna2','-','A','1','T'*21,'x','0','0']]
    class Process:
        returncode=0
        def __init__(self,command,**kw):
            counts['bowtie']+=1
            Path(command[-1]).write_text(''.join('\t'.join(r)+'\n' for r in rows))
        def communicate(self):return b'',b''
    monkeypatch.setitem(sys.modules,'Bio',SimpleNamespace(SeqIO=SimpleNamespace(parse=lambda *a:[SimpleNamespace(id='query',seq='A'*22)])))
    monkeypatch.setitem(sys.modules,'sifi_pipeline',SimpleNamespace(SifiPipeline=Pipeline))
    monkeypatch.setitem(sys.modules,'general_helpers',SimpleNamespace(get_target_data=helper))
    monkeypatch.setattr(bridge.subprocess,'Popen',Process)
    tool=tmp_path/'bowtie.exe';tool.write_text('test')
    query=tmp_path/'q.fasta';query.write_text('>query\n'+'A'*22)
    shared=tmp_path/'shared';shared.mkdir()
    def run(name,target=None,use_shared=False,list_only=False):
        out=tmp_path/name;out.mkdir()
        args=['--sifi-root',str(LEGACY),'--query-fasta',str(query),'--database-prefix','db',
              '--database-dir',str(tmp_path),'--bowtie',str(tool),'--rnaplfold-dir',str(tmp_path),
              '--out-json',str(out/'r.json'),'--out-mask',str(out/'m.json')]
        if target:args+=['--main-target',target]
        if use_shared:args+=['--shared-dir',str(shared)]
        if list_only:args+=['--list-targets']
        assert bridge.main(args)==0
        if not list_only:
            return json.loads((out/'r.json').read_text()),json.loads((out/'m.json').read_text()),Path(str(out/'r.json')+'.bowtie.tmp').read_bytes()
    a=run('original_A','A');b=run('original_B','B')
    counts.update(bowtie=0,rnaplfold=0,efficiency=0)
    run('list',use_shared=True,list_only=True)
    assert run('shared_A','A',True)==a
    assert run('shared_B','B',True)==b
    assert counts==dict(bowtie=1,rnaplfold=1,efficiency=1)
    # Changing a dependency cannot silently reuse the shared data.
    tool.write_text('changed tool')
    with pytest.raises(ValueError,match='inputs changed'):run('changed','A',True)


def test_session_is_scoped_and_removed():
    assert shared_directory.get() is None
    with shared_sifi_session():
        path=Path(shared_directory.get());assert path.is_dir()
        with shared_sifi_session():assert Path(shared_directory.get())!=path
        assert Path(shared_directory.get())==path
    assert not path.exists() and shared_directory.get() is None


def test_scalar_lazy_ranking_matches_existing_tie_break(monkeypatch):
    windows=select_windows(make_report([(1,185),(201,227)]),180,180)['windows']
    expected=sorted(range(len(windows)),key=lambda i:priority_key(windows[i]))
    accesses=[];original=type(windows).__getitem__
    def get(self,index):accesses.append(index);return original(self,index)
    monkeypatch.setattr(type(windows),'__getitem__',get)
    order=windows.ranked_indices()
    assert next(order)==expected[0]
    assert accesses==[]  # No row decoding or sorting to try the first candidate.
    assert list(order)==expected[1:]
    assert accesses==[]


def test_no_green_skips_indices_and_preserves_warnings(monkeypatch):
    import rnai_designer.automatic_window_selection as module
    report=make_report([(1,30)])
    for row in report['sifi']['records']:row['is_off_target']=True
    def forbidden(*a,**kw):pytest.fail('unnecessary index with no green regions')
    monkeypatch.setattr(module,'WindowIndex',forbidden)
    monkeypatch.setattr(module,'OfftargetIndex',forbidden)
    result=module.select_windows(report)
    assert result['evaluated_window_count']==0 and result['selected_window'] is None
    assert result['warnings']


def test_pin_gene_search_but_bound_window_cache(tmp_path):
    from rnai_designer.reuse import offtarget_reuse_session, search_cache, remember_search
    with offtarget_reuse_session():
        cache=search_cache.get()
        gene=tmp_path/'offtarget_runs';gene.mkdir();(gene/'r').write_text('gene')
        ev=[{'status':'imported_unvalidated'}]
        remember_search(cache,'gene',gene,ev)
        for i in range(10):
            path=tmp_path/str(i);path.mkdir();(path/'r').write_text('window')
            remember_search(cache,str(i),path,ev)
        assert 'gene' in cache and len(cache)==5
    assert search_cache.get() is None


@pytest.mark.parametrize('variant',['normal','missing','duplicates','no_green'])
def test_small_scientific_differential(variant):
    from benchmarks.reference.automatic_window_selection import select_windows as original
    report=make_report([(1,185)])
    if variant=='missing':
        for row in report['sifi']['records'][::2]:row['accessibility_value']=None
    if variant=='duplicates':report['sifi']['records'].append(copy.deepcopy(report['sifi']['records'][0]))
    if variant=='no_green':
        for row in report['sifi']['records']:row['is_off_target']=True
    old=original(copy.deepcopy(report),180,180)
    new=select_windows(copy.deepcopy(report),180,180)
    assert old['evaluated_window_count']==new['evaluated_window_count']
    assert old['warnings']==new['warnings']
    for a,b in zip(old['windows'],new['windows'],strict=True):
        a.pop('offtarget_matches')
        c=a['efficacy_priority_components']
        if a['efficacy_priority_score'] is not None:
            a['efficacy_priority_score']=100*(0.75*c['sifi_is_efficient_fraction']+0.25*c['sirna_gc_30_52_fraction'])
        assert a=={k:v for k,v in b.items() if k not in ('length_penalty','adjusted_priority_score')}


def test_efficiency_cache_uses_all_original_inputs():
    calls=[]
    def calculate(*args):calls.append(args);return args
    pipeline=SimpleNamespace(calculate_efficiency=calculate)
    bridge.memoize_efficiency(pipeline)
    for args in [('AAA',None,0.2),('AAA',None,0.2),('AAA','TTT',0.2),('AAA','TTT',0.3),('CCC','TTT',0.3)]:
        assert pipeline.calculate_efficiency(*args)==args
    assert len(calls)==4
