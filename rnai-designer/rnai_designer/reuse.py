"""Bounded-lifetime reuse of completed external searches across main targets."""
from contextlib import contextmanager
from contextvars import ContextVar
from collections import OrderedDict
from pathlib import Path
import json
import hashlib
import os
import shutil
import tempfile
from .cache import file_digest
from .atomic_io import write_json, new_staging_directory

search_cache = ContextVar('rnai_completed_searches', default=None)


@contextmanager
def offtarget_reuse_session():
    token = search_cache.set(OrderedDict())
    try:
        yield
    finally:
        search_cache.reset(token)


def search_key(candidate, settings, fingerprint):
    tool = settings.bowtie
    if tool is None:
        return None
    payload = {'schema': 1, 'candidate': candidate.original,
               'tool': (str(tool.resolve()), file_digest(tool)),
               'tool_dependencies': [(str(p), file_digest(p)) for p in sorted(tool.parent.glob('*.dll'))],
               'databases': [fingerprint(settings, species) for species in
                             ('Homo sapiens', 'Zea mays', 'Apis mellifera')],
               'parameters': {'sirna_size': 21, 'mismatches': 0, 'all_alignments': True, 'both_strands': True}}
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def restore_search(entry, root, key):
    source, evaluations, digests = entry
    try:
        if any(file_digest(source/name) != digest for name, digest in digests.items()):
            return None
        staging = new_staging_directory(root.parent, '.search.')
        try:
            for name, digest in digests.items():
                destination = staging/name
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source/name, destination)
                if file_digest(destination) != digest:
                    return None
            write_json(staging/'reuse_manifest.json', {'source_directory': str(source), 'content_key': key,
                                                      'note': 'Original completed execution and evidence copied byte-for-byte; no new external run.'})
            root.rmdir()  # Only an empty, freshly validated destination is eligible.
            os.replace(staging, root)
            return [{**evaluation, 'reused_from': str(source)} for evaluation in evaluations]
        finally:
            if staging.exists():
                shutil.rmtree(staging)
    except (OSError, ValueError):
        return None


def remember_search(cache, key, root, evaluations):
    if key is not None and all(e['status'] == 'imported_unvalidated' for e in evaluations):
        digests = {p.relative_to(root): file_digest(p) for p in root.rglob('*') if p.is_file()}
        cache[key] = root, evaluations, digests
        cache.move_to_end(key)
        # Keep the gene-wide run: otherwise several rejected windows evict it
        # and the identical whole-gene Bowtie searches repeat for the next GITV.
        window_keys = [k for k, item in cache.items() if item[0].name != "offtarget_runs"]
        for old_key in window_keys[:-4]:
            del cache[old_key]
