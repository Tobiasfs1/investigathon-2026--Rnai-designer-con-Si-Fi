import copy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
import rnai_designer.__main__ as cli
import rnai_designer.simple_run as simple
from rnai_designer.data_models import Candidate


@pytest.fixture
def cwd(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # pytest's descriptive temp directories are longer than the user's project.
    monkeypatch.setattr(simple, "MAX_PATH_CHARS", 1000)
    return tmp_path


@pytest.mark.parametrize("name", ["BicC_Dm", "Rpt5", "Mi gen muy largo / con caracteres", "CON", "NUL", "COM1", "../../escape", "\u57fa\u56e0"])
def test_names_and_normalization(cwd, name):
    calls=[]
    def workflow(args, **kwargs):calls.append((args,kwargs));return 0
    assert simple.simple_command([name,"ac gt\nAT\tgc"],workflow)==0
    fasta=next((cwd/'inputs').glob('*.fasta'));out=next((cwd/'output').iterdir())
    assert fasta.stem==out.name
    assert len(out.name.rsplit('_',1)[0])<=12
    assert fasta.read_text().splitlines()[1]=='ACGTATGC'
    assert calls[0][0]==['--gene-id',name,'--sequence-file',str(fasta),'--out',str(out),'--run-sifi-headless']
    assert calls[0][1]['input_metadata']['original_gene_name']==name


def test_existing_files_and_truncated_collisions(cwd):
    (cwd/'inputs').mkdir();(cwd/'output').mkdir()
    (cwd/'output'/'BicC_Dm_1').mkdir()
    old=cwd/'inputs'/'BicC_Dm_2.fasta';old.write_text('do not overwrite')
    calls=[]
    for _ in range(2):
        assert simple.simple_command(['BicC_Dm','ACGT'],lambda args,**kw: calls.append(args) or 0)==0
    assert (cwd/'output'/'BicC_Dm_3').is_dir()
    assert (cwd/'output'/'BicC_Dm_4').is_dir()
    assert old.read_text()=='do not overwrite'
    for name in ['abcdefghijkl111','abcdefghijkl222']:
        simple.simple_command([name,'ACGT'],lambda *a,**k:0)
    assert (cwd/'inputs'/'abcdefghijkl_1.fasta').exists()
    assert (cwd/'inputs'/'abcdefghijkl_2.fasta').exists()


@pytest.mark.parametrize('name,seq',[('', 'ACGT'),('  ','ACGT'),('gene',''),('gene','ACGT!'),('gene','1234')])
def test_invalid_input_creates_nothing(cwd,name,seq):
    assert simple.simple_command([name,seq],lambda *a,**k:pytest.fail('must not run'))==2
    assert not (cwd/'inputs').exists() and not (cwd/'output').exists()


@pytest.mark.parametrize('args',[[],['gene'],['gene','ACGT','extra'],['gene','ACGT','--out','x']])
def test_exactly_two_arguments(cwd,args):
    with pytest.raises(SystemExit):simple.simple_command(args,lambda *a,**k:None)
    assert not (cwd/'inputs').exists()


def test_long_sequence(cwd):
    sequence='acgt'*25000
    assert simple.simple_command(['long',sequence],lambda *a,**k:0)==0
    lines=(cwd/'inputs'/'long_1.fasta').read_text().splitlines()
    assert ''.join(lines[1:])==sequence.upper()
    assert max(map(len,lines[1:]))<=60


def test_path_length_checked_before_creation(cwd,monkeypatch,capsys):
    monkeypatch.setattr(simple,'MAX_PATH_CHARS',30)
    assert simple.simple_command(['gene','ACGT'],lambda *a,**k:pytest.fail('must not run'))==2
    assert 'Ruta demasiado larga' in capsys.readouterr().err
    assert not (cwd/'inputs').exists()


def test_failure_preserves_inputs_and_returns_code(cwd):
    assert simple.simple_command(['gene','ACGT'],lambda *a,**k:2)==2
    assert (cwd/'inputs'/'gene_1.fasta').is_file()
    assert (cwd/'output'/'gene_1').is_dir()


def test_simple_and_advanced_use_identical_analysis(cwd,monkeypatch,capsys):
    database_dir=cwd/'db';database_dir.mkdir()
    monkeypatch.setattr(cli.Settings,'load',lambda *a:SimpleNamespace(database_dir=database_dir))
    monkeypatch.setattr(cli,'_legacy_sifi_context',lambda args,*a:setattr(args,'_sifi_runtime',{}))
    monkeypatch.setattr(cli,'_list_legacy_targets',lambda *a:{'available_targets':[{'hit_name':'ref_A'},{'hit_name':'ref_B'}]})
    calls=[]
    def analyze(args,candidate,settings,out,targets):
        calls.append((copy.deepcopy(candidate),targets,(args.min_window,args.max_window,args.length_penalty_weight,args.legacy_sirna_size,args.sifi_mismatches,args.timeout,args.run_sifi_headless,args.verbose_report)))
        return {'candidate_id':candidate.gene_id,'status':'complete','reference_count':2,
                'completed_reference_count':2,'failed_reference_count':0,'combined_window_count':100}
    monkeypatch.setattr(cli,'_run_multi_reference_workflow',analyze)
    name='Original gene name over twelve characters';seq='ac gt'*100
    assert cli.main([name,seq])==0
    assert cli.main(['workflow','--gene-id',name,'--sequence',seq,'--out',str(cwd/'advanced'),'--run-sifi-headless'])==0
    a,b=calls
    assert (a[0].gene_id,a[0].sequence,a[1:])==(b[0].gene_id,b[0].sequence,b[1:])
    assert a[0].original['reference']==b[0].original['reference']
    stored=json.loads(next((cwd/'output').glob('*/candidate.json')).read_text(encoding='utf-8'))
    assert stored['metadata']['original_gene_name']==name
    assert stored['gene_id']==name
    assert 'Ventanas evaluadas' in capsys.readouterr().out


def test_advanced_dispatch_remains(monkeypatch):
    monkeypatch.setattr(cli,'workflow_command',lambda args: ('workflow',args))
    monkeypatch.setattr(cli,'utility',lambda command,args:(command,args))
    assert cli.main(['workflow','--resume','x'])==('workflow',['--resume','x'])
    assert cli.main(['doctor'])==('doctor',[])
