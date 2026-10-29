"""Publish files/directories only after successful completion, including interrupts."""
from contextlib import contextmanager
from functools import wraps
import json
import os
from pathlib import Path
import shutil
import tempfile
import inspect
from uuid import uuid4


def new_staging_directory(parent, prefix):
    # mkdtemp creates owner-only ACLs on recent Windows/Python. A directory that
    # will be published must inherit its intended parent's access permissions.
    path = Path(parent) / (prefix + uuid4().hex + '.partial')
    path.mkdir()
    return path


@contextmanager
def temporary_directory(parent, prefix):
    path = new_staging_directory(parent, prefix)
    try:
        yield path
    finally:
        if path.exists():
            shutil.rmtree(path)


@contextmanager
def atomic_text(path, **kwargs):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.' + path.name + '.', suffix='.partial', dir=path.parent)
    try:
        with os.fdopen(fd, 'w', encoding='utf-8', **kwargs) as handle:
            yield handle
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def write_json(path, value):
    with atomic_text(path) as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, allow_nan=False)


def _publish_existing(staging, out):
    """Replace generated files transactionally; preserve unrelated input/evidence."""
    backup = new_staging_directory(out.parent, '.' + out.name + '.previous.')
    published = []
    clean_backup = False
    marker = 'report.json'
    try:
        files = sorted(staging.iterdir(), key=lambda p: (p.name == marker, p.name))
        if (out/marker).exists():
            if not (out/marker).is_file():
                raise IsADirectoryError(out/marker)
            os.replace(out/marker, backup/marker)
        for source in files:
            destination = out/source.name
            if destination.exists():
                if not destination.is_file():
                    raise IsADirectoryError(destination)
                os.replace(destination, backup/source.name)
            os.replace(source, destination)
            published.append(destination)
        clean_backup = True
    except BaseException:
        for destination in reversed(published):
            destination.unlink()
        # Restore the previous completion marker last, after all prior data.
        for original in sorted(backup.iterdir(), key=lambda p: (p.name == marker, p.name)):
            os.replace(original, out/original.name)
        clean_backup = True
        raise
    finally:
        # A failed rollback leaves recoverable originals in the backup directory.
        if clean_backup:
            shutil.rmtree(backup)


def staged_directory(function=None, *, allow_existing=False):
    """Publish the bound 'out' directory only after all artifacts are complete."""
    if function is None:
        return lambda target: staged_directory(target, allow_existing=allow_existing)
    @wraps(function)
    def wrapped(*args, **kwargs):
        bound = inspect.signature(function).bind(*args, **kwargs)
        out = Path(bound.arguments['out'])
        if out.exists() and any(out.iterdir()) and not allow_existing:
            raise ValueError('La carpeta de salida ya contiene archivos')
        out.parent.mkdir(parents=True, exist_ok=True)
        staging = new_staging_directory(out.parent, '.' + out.name + '.')
        try:
            bound.arguments['out'] = staging
            result = function(*bound.args, **bound.kwargs)
            if allow_existing and out.exists():
                _publish_existing(staging, out)
                return result
            if out.exists():
                out.rmdir()
            os.replace(staging, out)
            return result
        finally:
            if staging.exists():
                shutil.rmtree(staging)
    return wrapped
