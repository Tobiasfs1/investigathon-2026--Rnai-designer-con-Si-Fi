"""Content-addressed optional headless cache; digest reuse is process-local only."""
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import os

_digests = {}


def signature(path):
    s = Path(path).stat()
    return (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)


def file_digest(path):
    path = Path(path).resolve()
    before = signature(path)
    cached = _digests.get(path)
    if cached and cached[0] == before:
        return cached[1]
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            digest.update(block)
    if signature(path) != before:
        raise ValueError(f'El archivo cambió durante el cálculo del hash: {path}')
    result = digest.hexdigest()
    _digests[path] = (before, result)
    return result


def content_key(parameters, files):
    payload = {'cache_schema': 2, 'parameters': parameters,
               'files': [(str(Path(p).resolve()), file_digest(p)) for p in sorted(files, key=str)]}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


class ArtifactCache:
    def __init__(self, root, key):
        self.directory = Path(root) / key

    def restore(self, destinations):
        try:
            manifest = json.loads((self.directory/'complete.json').read_text(encoding='utf-8'))
            if set(manifest) != set(destinations):
                return False
            # Verify every byte before copying any artifact. Never trust mere file existence.
            for name, digest in manifest.items():
                if file_digest(self.directory/name) != digest:
                    return False
            parent = Path(next(iter(destinations.values()))).parent
            parent.mkdir(parents=True, exist_ok=True)
            from .atomic_io import temporary_directory
            with temporary_directory(parent, '.restore.') as temporary:
                for name in destinations:
                    copied = Path(temporary)/name
                    shutil.copyfile(self.directory/name, copied)
                    if file_digest(copied) != manifest[name]:
                        return False
                for name, destination in destinations.items():
                    Path(destination).parent.mkdir(parents=True, exist_ok=True)
                    os.replace(Path(temporary)/name, destination)
            return True
        except (OSError, ValueError, TypeError):
            return False

    def store(self, sources):
        if self.directory.exists():
            return
        self.directory.parent.mkdir(parents=True, exist_ok=True)
        from .atomic_io import new_staging_directory
        staging = new_staging_directory(self.directory.parent, '.cache.')
        try:
            manifest = {}
            for name, source in sources.items():
                shutil.copyfile(source, staging/name)
                manifest[name] = file_digest(staging/name)
            from .atomic_io import write_json
            write_json(staging/'complete.json', manifest)
            try:
                os.rename(staging, self.directory)
            except FileExistsError:
                pass  # Another completed writer already published the same key.
        finally:
            if staging.exists():
                shutil.rmtree(staging)
