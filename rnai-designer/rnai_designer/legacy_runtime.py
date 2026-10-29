"""Portable, read-only discovery of the Python 2 si-Fi runtime."""
from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess

SOURCE_FILES = ('sifi_pipeline.py', 'general_helpers.py', 'show_plot.py')


def find_sifi_root(explicit=None, base=None):
    configured_base = Path(base or Path.cwd()).resolve()
    base = Path(base or Path(__file__).resolve().parents[1]).resolve()
    configured = explicit or os.environ.get('SIFI_ROOT')
    if configured:
        path = Path(os.path.expandvars(str(configured))).expanduser()
        path = (path if path.is_absolute() else configured_base / path).resolve()
        if all((path / name).is_file() for name in SOURCE_FILES):
            return path
        raise ValueError('Código original de si-Fi no encontrado en %s; se requieren %s' % (path, ', '.join(SOURCE_FILES)))
    candidates = [base, base.parent]
    for parent in (base, base.parent):
        try:
            candidates.extend(sorted(p for p in parent.iterdir() if p.is_dir()))
        except OSError:
            pass
    for path in candidates:
        if all((path / name).is_file() for name in SOURCE_FILES):
            return path.resolve()
    raise ValueError('Código original de si-Fi no encontrado junto al proyecto; configurar SIFI_ROOT o --sifi-root (%s).' % ', '.join(SOURCE_FILES))


def find_conda() -> Path | None:
    candidates = []
    for name in ('conda.exe', 'conda'):
        located = shutil.which(name)
        if located:
            candidates.append(Path(located))
    if os.environ.get('CONDA_EXE'):
        candidates.append(Path(os.path.expandvars(os.environ['CONDA_EXE'])))
    roots = []
    if os.environ.get('CONDA_PREFIX'):
        prefix = Path(os.environ['CONDA_PREFIX'])
        roots.extend((prefix, prefix.parent.parent))
    parents = [Path.home()]
    parents.extend(Path(os.environ[name]) for name in ('LOCALAPPDATA', 'PROGRAMDATA', 'ProgramFiles') if os.environ.get(name))
    for parent in parents:
        for name in ('miniconda3', 'miniforge3', 'anaconda3', 'mambaforge', 'Miniconda3', 'Anaconda3'):
            roots.extend((parent / name, parent / 'Programs' / name))
    for root in roots:
        candidates.extend((root / 'Scripts' / 'conda.exe', root / 'condabin' / 'conda.bat', root / 'bin' / 'conda'))
    return next((p.resolve() for p in candidates if p.is_file()), None)


def _run(command, timeout=30):
    return subprocess.run(command, capture_output=True, text=True, errors='replace', timeout=timeout,
                          creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))


def _probe(conda: Path, prefix: str, sifi_root: Path, timeout=60) -> dict:
    # Python 2 compatible; import only the modules used by the headless adapters.
    code = '''import sys, struct, json, importlib
result = {'python':sys.version.split()[0], 'bits':struct.calcsize('P')*8, 'packages':{}, 'package_versions':{}, 'source_imports':[], 'module_paths':{}, 'errors':[]}
if sys.version_info[:2] != (2,7):
    result['errors'].append({'kind':'python_incompatible','detail':'Se requiere Python 2.7'})
sys.path.insert(0, %s)
for name in ['PyQt4','Bio.Alphabet','numpy','matplotlib','seaborn','scipy','sifi_pipeline','general_helpers','show_plot','free_energy']:
    try:
        module = importlib.import_module(name)
        result['packages'][name] = True
        result['package_versions'][name] = getattr(module,'__version__',None)
        if name in ['sifi_pipeline','general_helpers','show_plot','free_energy']:
            result['source_imports'].append(name)
            result['module_paths'][name] = getattr(module, '__file__', None)
    except Exception as error:
        result['packages'][name] = False
        result['errors'].append({'kind':'dependency_missing' if name in ['PyQt4','Bio.Alphabet','numpy','matplotlib','seaborn','scipy'] else 'import_error','module':name,'detail':str(error)})
result['available'] = not result['errors']
print(json.dumps(result))
''' % ('u' + ascii(str(sifi_root)))
    try:
        process = _run([str(conda), 'run', '-p', prefix, 'python', '-c', 'exec(' + ascii(code) + ')'], timeout)
    except (OSError, subprocess.TimeoutExpired) as error:
        return {'available': False, 'reason': 'No se pudo ejecutar Python del entorno: ' + str(error)}
    for line in reversed(process.stdout.splitlines()):
        try:
            result = json.loads(line)
            if not isinstance(result, dict) or 'available' not in result:
                continue
            if process.returncode:
                result['available'] = False
            labels = {'python_incompatible':'Python incompatible', 'dependency_missing':'Dependencia faltante', 'import_error':'Error al importar módulos'}
            result['reason'] = '; '.join(labels[e['kind']] + ': ' + e.get('module','') + ' ' + e['detail'] for e in result['errors'])
            return result
        except ValueError:
            continue
    return {'available': False, 'reason': 'Error al ejecutar/importar módulos: ' + (process.stderr or process.stdout)[:4000]}


def detect_sifi_runtime(sifi_root=None, env_name: str | None = None, timeout=60) -> dict:
    try:
        sifi_root = find_sifi_root(sifi_root)
    except ValueError as error:
        return {'available': False, 'reason': str(error)}
    conda = find_conda()
    if conda is None:
        return {'available': False, 'reason': 'Conda no encontrado; instalar Conda o configurar PATH/CONDA_EXE.'}
    try:
        process = _run([str(conda), 'env', 'list', '--json'], timeout)
        if process.returncode:
            raise ValueError((process.stderr or process.stdout)[-2000:])
        envs = json.loads(process.stdout)['envs']
        if not isinstance(envs, list):
            raise ValueError('envs debe ser una lista')
    except (OSError, ValueError, KeyError, subprocess.TimeoutExpired) as error:
        return {'available': False, 'conda': str(conda), 'reason': 'conda env list --json falló: ' + str(error)}
    selected = env_name or os.environ.get('SIFI_CONDA_ENV')
    active = os.environ.get('CONDA_PREFIX')
    if active and active not in envs:
        envs.append(active)
    if selected and (Path(selected).is_absolute() or '/' in selected or '\\' in selected):
        envs = [str(Path(os.path.expandvars(selected)).expanduser().resolve())]
    else:
        envs = [p for p in envs if (Path(p).name.casefold() == selected.casefold() if selected else 'sifi' in Path(p).name.casefold())]
    envs = sorted(dict.fromkeys(envs), key=lambda p: Path(p).name.casefold() != 'sifi_32bit')
    probes = []
    for prefix in envs:
        probe = _probe(conda, prefix, sifi_root, timeout)
        probe.update(name=Path(prefix).name, prefix=str(Path(prefix).resolve()), sifi_root=str(sifi_root), interpreter=str(Path(prefix) / ('python.exe' if os.name == 'nt' else 'bin/python')))
        probes.append(probe)
        if probe.get('available'):
            probe.update(conda=str(conda), launch_command=[str(conda), 'run', '-p', probe['prefix'], 'python', str(sifi_root / 'main.py')])
            return probe
    return {'available': False, 'conda': str(conda), 'sifi_root': str(sifi_root), 'probed_environments': probes,
            'reason': '; '.join(p.get('reason', 'Error al validar entorno') for p in probes) if probes else 'Entorno no encontrado: %s. Configurar SIFI_CONDA_ENV o --legacy-conda-env (nombre o prefijo).' % (selected or 'sifi_32bit')}
