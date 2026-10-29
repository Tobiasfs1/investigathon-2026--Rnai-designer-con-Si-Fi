import json
from pathlib import Path
from types import SimpleNamespace
import pytest
from rnai_designer import legacy_runtime as runtime
from rnai_designer.config import Settings, doctor
from rnai_designer.__main__ import main


def source(path):
    path.mkdir(parents=True)
    for name in runtime.SOURCE_FILES:
        (path / name).touch()
    return path


@pytest.fixture(autouse=True)
def clean_env(monkeypatch):
    for name in ('CONDA_EXE','CONDA_PREFIX','SIFI_ROOT','SIFI_CONDA_ENV','SIFI_HOME','SIFI_DB_DIR','SIFI_BOWTIE','SIFI_BOWTIE_INSPECT','SIFI_RNAPLFOLD'):
        monkeypatch.delenv(name, raising=False)


@pytest.mark.parametrize('drive,user', [('C','Ana'),('D','Other User'),('Z','Usuario con espacios')])
def test_project_sibling_and_external(tmp_path, drive, user):
    project = tmp_path / drive / user / 'project' / 'rnai-designer'
    project.mkdir(parents=True)
    sibling = source(project.parent / 'original renombrado')
    assert runtime.find_sifi_root(base=project) == sibling
    external = source(tmp_path / 'external')
    assert runtime.find_sifi_root(external, project) == external
    with pytest.raises(ValueError, match='no encontrado'):
        runtime.find_sifi_root(tmp_path / 'missing', project)


def test_conda_path_has_priority(tmp_path, monkeypatch):
    first, second = tmp_path/'conda.exe', tmp_path/'other.exe'
    first.touch(); second.touch()
    monkeypatch.setattr(runtime.shutil, 'which', lambda name: str(first))
    monkeypatch.setenv('CONDA_EXE', str(second))
    assert runtime.find_conda() == first


@pytest.mark.parametrize('location', ['home','local','programs','prefix','explicit'])
def test_conda_fallbacks(tmp_path, monkeypatch, location):
    monkeypatch.setattr(runtime.shutil, 'which', lambda name: None)
    monkeypatch.setattr(Path, 'home', classmethod(lambda cls: tmp_path/'user'))
    monkeypatch.setenv('LOCALAPPDATA', str(tmp_path/'local'))
    roots = {'home':tmp_path/'user'/'anaconda3', 'local':tmp_path/'local'/'miniforge3',
             'programs':tmp_path/'local'/'Programs'/'miniconda3', 'prefix':tmp_path/'custom', 'explicit':tmp_path/'another'}
    exe = roots[location]/'Scripts'/'conda.exe'
    exe.parent.mkdir(parents=True); exe.touch()
    if location == 'prefix': monkeypatch.setenv('CONDA_PREFIX', str(roots[location]/'envs'/'sifi_32bit'))
    if location == 'explicit': monkeypatch.setenv('CONDA_EXE', str(exe))
    assert runtime.find_conda() == exe


@pytest.mark.parametrize('selection', ['sifi_32bit','custom', 'prefix'])
def test_select_environment(tmp_path, monkeypatch, selection):
    root = source(tmp_path/'original')
    prefix = tmp_path/'envs'/('custom' if selection != 'sifi_32bit' else selection)
    monkeypatch.setattr(runtime, 'find_conda', lambda: tmp_path/'conda.exe')
    monkeypatch.setattr(runtime, '_run', lambda *a, **k: SimpleNamespace(returncode=0, stdout=json.dumps({'envs':[str(prefix)]})))
    seen=[]
    def probe(conda, selected, root, timeout):
        seen.append(selected)
        return {'available':True, 'python':'2.7.18'}
    monkeypatch.setattr(runtime, '_probe', probe)
    result=runtime.detect_sifi_runtime(root, str(prefix) if selection=='prefix' else selection)
    assert result['available'] and seen==[str(prefix)]


def test_missing_environment_and_conda(tmp_path, monkeypatch):
    root=source(tmp_path/'source')
    monkeypatch.setattr(runtime,'find_conda',lambda:None)
    assert 'Conda no encontrado' in runtime.detect_sifi_runtime(root)['reason']
    monkeypatch.setattr(runtime,'find_conda',lambda:tmp_path/'conda.exe')
    monkeypatch.setattr(runtime,'_run',lambda *a, **k:SimpleNamespace(returncode=0,stdout='{"envs":[]}'))
    assert 'Entorno no encontrado' in runtime.detect_sifi_runtime(root)['reason']


@pytest.mark.parametrize('kind,label', [('python_incompatible','Python incompatible'),('dependency_missing','Dependencia faltante'),('import_error','Error al importar')])
def test_probe_diagnostics(tmp_path, monkeypatch, kind, label):
    data={'available':False,'errors':[{'kind':kind,'module':'numpy','detail':'failed'}]}
    monkeypatch.setattr(runtime,'_run',lambda *a, **k:SimpleNamespace(returncode=0,stdout=json.dumps(data)))
    assert label in runtime._probe(tmp_path/'conda',str(tmp_path/'env'),tmp_path)['reason']


def test_real_probe_rejects_python3(tmp_path, monkeypatch):
    import subprocess, sys
    monkeypatch.setattr(runtime,'_run',lambda command, timeout:subprocess.run([sys.executable,'-c',command[-1]],capture_output=True,text=True))
    result=runtime._probe(tmp_path/'conda',str(tmp_path),tmp_path)
    assert not result['available'] and 'Python incompatible' in result['reason']


def test_configuration_and_overrides(tmp_path, monkeypatch):
    config=tmp_path/'config.json'
    config.write_text(json.dumps({'sifi_root':'../original','legacy_conda_env':'custom','sifi_home':'portable','bowtie':'tools/bowtie.exe'}))
    settings=Settings.load(config)
    assert settings.sifi_root==(tmp_path.parent/'original').resolve()
    assert settings.legacy_conda_env=='custom'
    assert settings.database_dir==tmp_path/'portable'/'Databases'
    assert settings.bowtie==tmp_path/'tools'/'bowtie.exe'
    monkeypatch.setenv('SIFI_DB_DIR',str(tmp_path/'override'))
    monkeypatch.setenv('SIFI_CONDA_ENV','override')
    assert Settings.load(config).database_dir==tmp_path/'override'
    assert Settings.load(config).legacy_conda_env=='override'


def test_localappdata_tools_and_missing_doctor(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv('LOCALAPPDATA',str(tmp_path/'another user'))
    settings=Settings.load()
    assert settings.database_dir==tmp_path/'another user'/'siFi2015'/'Databases'
    monkeypatch.setattr(runtime,'detect_sifi_runtime',lambda *a: {'available':False,'reason':'Conda no encontrado'})
    result=doctor(settings)
    assert all(not entry['available'] for entry in result['databases']['configured_indexes'])
    assert main(['doctor'])==0
    assert 'configured_indexes' in capsys.readouterr().out


def test_probe_uses_configured_source_path(tmp_path, monkeypatch):
    import subprocess, sys
    root=source(tmp_path/'original with spaces')
    (root/'free_energy.py').touch()
    monkeypatch.setattr(runtime,'_run',lambda command, timeout:subprocess.run([sys.executable,'-c',command[-1]],capture_output=True,text=True))
    result=runtime._probe(tmp_path/'conda',str(tmp_path/'env'),root)
    assert set(result['source_imports'])=={'sifi_pipeline','general_helpers','show_plot','free_energy'}
    assert Path(result['module_paths']['sifi_pipeline']).parent==root


def test_config_relative_environment_prefix(tmp_path):
    config=tmp_path/'config.json'
    config.write_text(json.dumps({'legacy_conda_env':'../envs/custom','database_dir':'db'}))
    assert Settings.load(config).legacy_conda_env==str((tmp_path.parent/'envs'/'custom').resolve())


def test_source_protected_outside_project(tmp_path):
    from rnai_designer.__main__ import protect_output
    root=source(tmp_path/'external')
    with pytest.raises(ValueError):
        protect_output(root/'outputs',sifi_root=root)
    assert not (root/'outputs').exists()


def test_environment_probe_failure_is_preserved(tmp_path, monkeypatch):
    root=source(tmp_path/'source')
    monkeypatch.setattr(runtime,'find_conda',lambda:tmp_path/'conda.exe')
    monkeypatch.setattr(runtime,'_run',lambda *a, **k:SimpleNamespace(returncode=0,stdout=json.dumps({'envs':[str(tmp_path/'sifi_32bit')]})))
    monkeypatch.setattr(runtime,'_probe',lambda *a:{'available':False,'reason':'Dependencia faltante: numpy'})
    result=runtime.detect_sifi_runtime(root)
    assert 'numpy' in result['reason'] and len(result['probed_environments'])==1
