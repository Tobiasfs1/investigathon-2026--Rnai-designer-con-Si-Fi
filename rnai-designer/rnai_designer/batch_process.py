"""Wall-clock supervision for one batch gene in a separate process."""
from dataclasses import dataclass, asdict
from tempfile import TemporaryDirectory
import json
import os
from pathlib import Path
import subprocess
import sys
from time import monotonic

from .atomic_io import write_json
from .data_models import read_json
from .process_tree import ProcessTree

DEFAULT_GENE_TIMEOUT = 360.0
_WORKER_MODULE = "rnai_designer.batch_worker"


@dataclass
class WorkerResult:
    returncode: int
    elapsed_s: float
    timed_out: bool = False
    cleanup_confirmed: bool = True
    pid: int | None = None
    stage: str | None = None
    reference_id: str | None = None


class CleanupFailure(OSError):
    def __init__(self, error, timed_out, elapsed_s):
        super().__init__("No se pudo confirmar la limpieza de procesos: " + str(error))
        self.timed_out, self.elapsed_s = timed_out, elapsed_s


def run_worker(command, request, log_path, timeout, *, started=None):
    """stdin handshake prevents descendants escaping containment at startup."""
    started = monotonic() if started is None else started
    deadline = started+timeout if timeout else None
    request = {**request, 'started': started, 'deadline': deadline}
    process = None
    tree = ProcessTree()
    timed_out = False
    private_temp = None
    try:
        private_temp = TemporaryDirectory(prefix="rg_", ignore_cleanup_errors=True)
        with Path(log_path).open('x',encoding='utf-8') as log:
            env = os.environ.copy()
            env.update(PYTHONIOENCODING='utf-8', PYTHONUNBUFFERED='1',
                       TMP=private_temp.name, TEMP=private_temp.name, TMPDIR=private_temp.name)
            package_root = str(Path(__file__).resolve().parents[1])
            env['PYTHONPATH'] = package_root + os.pathsep + env.get('PYTHONPATH','')
            flags = getattr(subprocess,'CREATE_NO_WINDOW',0)
            process = subprocess.Popen(command,stdin=subprocess.PIPE,stdout=log,stderr=log,
                          env=env,creationflags=flags,start_new_session=os.name != 'nt')
            try:
                tree.attach(process)
            except BaseException:
                # No request was sent, so this uncontained worker has no children.
                process.kill();process.wait(timeout=5)
                raise
            try:
                remaining = None if deadline is None else max(0,deadline-monotonic())
                if remaining == 0:
                    timed_out = True
                else:
                    process.stdin.write((json.dumps(request)+'\n').encode('utf-8'))
                    process.stdin.close()
                    remaining = None if deadline is None else max(0,deadline-monotonic())
                    if remaining == 0:
                        timed_out = True
                    else:
                        process.wait(timeout=remaining)
            except subprocess.TimeoutExpired:
                timed_out = True
            finally:
                try:
                    if timed_out and request.get('supervision_state'):
                        # Publish cancellation BEFORE killing: a parent crash during
                        # cleanup must never turn a timeout into recovered success.
                        write_json(request['supervision_state'],{'status':'incomplete_timeout',
                                   'timeout_seconds':timeout,'elapsed_s':monotonic()-started})
                finally:
                    # Also remove descendants after a normal/failed worker exit.
                    try:
                        tree.stop()
                    except (OSError,subprocess.TimeoutExpired) as error:
                        raise CleanupFailure(error,timed_out,monotonic()-started) from error
    finally:
        try:
            if process and process.stdin and not process.stdin.closed:
                process.stdin.close()
        finally:
            try:
                tree.close()  # Also cancels descendants if the parent dies.
            finally:
                if private_temp:
                    private_temp.cleanup()
    progress = {}
    if request.get('progress'):
        try:
            progress = read_json(request['progress'])
        except (OSError, ValueError):
            pass
    return WorkerResult(124 if timed_out else process.returncode,monotonic()-started,
                        timed_out,pid=process.pid,stage=progress.get('stage'),
                        reference_id=progress.get('reference_id'))


def run_gene_worker(workflow, argv, metadata, log_path, timeout, *, started=None):
    # workflow remains an injectable callback at the batch boundary for small
    # unit tests; the production worker imports the same CLI in a new interpreter.
    log_path = Path(log_path)
    state = log_path.with_suffix('.worker.json')
    receipt = log_path.with_suffix('.done.json')
    progress = log_path.with_suffix('.progress.json')
    request_file = log_path.with_suffix('.request.json')
    for path in (log_path,state,receipt,progress,request_file):
        if path.exists() or path.is_symlink():
            raise FileExistsError('Ya existe un archivo del intento: ' + str(path))
    write_json(state, {'status':'running','timeout_seconds':timeout})
    write_json(request_file,{'argv':argv,'metadata':metadata,'receipt':str(receipt),'progress':str(progress)})
    # Only a short path crosses stdin; even unusually long FASTA headers cannot
    # block the supervising parent on a full pipe.
    request = {'request_file':str(request_file),'progress':str(progress),'supervision_state':str(state)}
    try:
        result = run_worker([sys.executable,'-B','-m',_WORKER_MODULE],
                            request,log_path,timeout,started=started)
    except CleanupFailure as error:
        write_json(state,{'status':'incomplete_timeout' if error.timed_out else 'cleanup_error',
                         'cleanup_confirmed':False,'timeout_seconds':timeout,'error':str(error)})
        raise
    except KeyboardInterrupt:
        write_json(state,{'status':'interrupted','timeout_seconds':timeout})
        raise
    try:
        confirmation = read_json(receipt)
    except (OSError,ValueError):
        confirmation = {}
    if timeout and confirmation.get('within_deadline') is False:
        result.timed_out = True
        result.returncode = 124
    elif result.returncode == 0 and confirmation.get('return_code') != 0:
        result.returncode = 2  # No successful completion receipt, no complete gene.
    write_json(state,{'status':'incomplete_timeout' if result.timed_out else 'finished',
                      'timeout_seconds':timeout,**asdict(result)})
    return result


def completed_worker(root, attempt):
    """Legacy attempts retain recovery; managed attempts need a successful receipt."""
    if not attempt.get('worker_state'):
        return True
    try:
        state = read_json(Path(root)/attempt['worker_state'])
        if state.get('status') in ('incomplete_timeout','cleanup_error') or state.get('cleanup_confirmed') is False:
            return False
        receipt = read_json(Path(root)/attempt['worker_receipt'])
        return (receipt.get('return_code') == 0 and receipt.get('within_deadline') is True
                and (state.get('status') != 'finished' or
                     state.get('returncode') == 0 and not state.get('timed_out')))
    except (OSError,ValueError,KeyError):
        return False
