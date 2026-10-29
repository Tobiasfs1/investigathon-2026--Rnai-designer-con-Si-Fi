"""Read-only discovery and invocation helpers for the legacy si-Fi Conda runtime."""
from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess


def find_conda() -> Path | None:
    candidates = []
    configured = os.environ.get("CONDA_EXE")
    if configured:
        candidates.append(Path(configured))
    located = shutil.which("conda")
    if located:
        candidates.append(Path(located))
    local = os.environ.get("LOCALAPPDATA")
    program_data = os.environ.get("PROGRAMDATA", r"C:\ProgramData")
    for root in (Path(local) / "miniforge3" if local else None,
                 Path(local) / "miniconda3" if local else None,
                 Path(program_data) / "miniforge3", Path(program_data) / "miniconda3"):
        if root:
            candidates.extend((root / "Scripts" / "conda.exe", root / "condabin" / "conda.bat"))
    return next((path.resolve() for path in candidates if path.is_file()), None)


def _run(command, timeout=30):
    return subprocess.run(command, capture_output=True, text=True, errors="replace", timeout=timeout,
                          creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))


def _probe(conda: Path, prefix: str, sifi_root: Path, timeout=60) -> dict:
    code = ("import sys,struct,json,importlib; "
            "names=['PyQt4','Bio.Alphabet','numpy','matplotlib','seaborn','scipy']; "
            "result={'python':sys.version.split()[0],'bits':struct.calcsize('P')*8,'packages':{}}; "
            "[(result['packages'].__setitem__(n, bool(importlib.import_module(n)))) for n in names]; "
            "version_names=['Bio','numpy','matplotlib','seaborn','scipy']; "
            "result['package_versions']={n:getattr(importlib.import_module(n),'__version__',None) for n in version_names}; "
            "from PyQt4 import QtCore; result['package_versions']['PyQt4']=QtCore.PYQT_VERSION_STR; "
            "result['package_versions']['Qt']=QtCore.QT_VERSION_STR; "
            "sys.path.insert(0,%r); "
            "mods=['main','sifi_pipeline','general_helpers','show_plot','free_energy']; "
            "[(importlib.import_module(n)) for n in mods]; "
            "result.update(source_imports=mods); print(json.dumps(result))" % str(sifi_root))
    process = _run([str(conda), "run", "-p", prefix, "python", "-c", code], timeout)
    if process.returncode:
        return {"available": False, "probe_error": (process.stderr or process.stdout)[-3000:]}
    # conda may print activation notices; consume the last JSON object line.
    for line in reversed(process.stdout.splitlines()):
        try:
            result = json.loads(line)
            result["available"] = True
            return result
        except ValueError:
            continue
    return {"available": False, "probe_error": "El entorno no devolvió el resultado JSON esperado."}


def detect_sifi_runtime(sifi_root: Path, env_name: str | None = None, timeout=60) -> dict:
    """Find an existing Python 2 si-Fi environment without installing or changing anything."""
    conda = find_conda()
    if conda is None:
        return {"available": False, "reason": "No se encontró conda.exe."}
    process = _run([str(conda), "env", "list", "--json"], timeout=timeout)
    if process.returncode:
        return {"available": False, "conda": str(conda),
                "reason": "conda env list falló: " + (process.stderr or process.stdout)[-2000:]}
    try:
        envs = json.loads(process.stdout).get("envs", [])
    except (ValueError, AttributeError):
        return {"available": False, "conda": str(conda), "reason": "No se pudo leer conda env list --json."}
    candidates = []
    for prefix in envs:
        name = Path(prefix).name
        if env_name and name.casefold() != env_name.casefold():
            continue
        if not env_name and "sifi" not in name.casefold():
            continue
        candidates.append((prefix, name))
    probes = []
    for prefix, name in candidates:
        probe = _probe(conda, prefix, sifi_root.resolve(), timeout)
        probe.update(name=name, prefix=str(Path(prefix).resolve()))
        probes.append(probe)
        if probe.get("available"):
            probe.update(conda=str(conda.resolve()),
                         launch_command=[str(conda.resolve()), "run", "-p", str(Path(prefix).resolve()),
                                         "python", str((sifi_root / "main.py").resolve())])
            return probe
    return {"available": False, "conda": str(conda.resolve()), "probed_environments": probes,
            "reason": ("No hay un entorno Conda con las dependencias y módulos originales de si-Fi."
                       if candidates else "No se encontró un entorno cuyo nombre contenga 'sifi'.")}
