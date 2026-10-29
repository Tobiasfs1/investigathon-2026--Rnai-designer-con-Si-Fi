"""Resumable batch orchestration; scientific work stays in workflow."""
import argparse
from contextlib import contextmanager
import csv
import hashlib
import math
from pathlib import Path
import sys
from time import monotonic

from .atomic_io import atomic_text, write_json
from .data_models import Candidate, read_json
from .simple_run import reserve_run, check_paths
from .priority_scoring import priority_key
from .batch_process import DEFAULT_GENE_TIMEOUT, run_gene_worker, completed_worker, CleanupFailure

NOTE = "Ranking heuristico; no es una prediccion validada de eficacia biologica."
FIELDS = ["input_index", "gene_id", "status", "error", "ranking_status", "rank",
          "reference_id", "start_1based", "end_1based_inclusive", "length_nt", "window_class",
          "efficacy_priority_score", "length_penalty", "adjusted_priority_score", "quality_score",
          "data_sufficiency_status", "gc_percent", "sirna_record_count", "sirna_unique_count",
          "is_efficient_true_count", "is_efficient_evaluable_count", "efficient_percent",
          "sirna_unique_sequences_count", "diversity_percent", "human_offtarget", "maize_offtarget",
          "bee_offtarget", "offtarget_filter_status", "sequence_release_status", "output_directory"]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_fasta(path):
    records = []
    name, parts = None, []
    def append():
        entry = {"input_index": len(records)+1, "gene_id": name, "sequence": "".join(parts),
                 "status": "pending", "attempts": []}
        try:
            candidate = Candidate.from_fields(name, entry["sequence"])
            entry["sequence"] = candidate.sequence
        except ValueError as error:
            entry.update(status="invalid_input", error=str(error))
        records.append(entry)
    for line_no, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                append()
            name, parts = line[1:].strip(), []
        elif name is None:
            raise ValueError(f"FASTA: secuencia sin encabezado en linea {line_no}")
        else:
            parts.append("".join(line.split()))
    if name is not None:
        append()
    if not records:
        raise ValueError("FASTA vacio: no contiene genes")
    return records


@contextmanager
def batch_lock(root):
    # The OS releases the lock on process death, so an interrupted run is resumable.
    with (root / ".lock").open("a+b") as handle:
        handle.seek(0, 2)
        if not handle.tell():
            handle.write(b"0")
            handle.flush()
        handle.seek(0)
        try:
            if sys.platform == "win32":
                import msvcrt
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            raise ValueError("Este batch ya esta en ejecucion") from error
        try:
            yield
        finally:
            handle.seek(0)
            if sys.platform == "win32":
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def approved(reference, window):
    if not window:
        return False
    filt = reference.get("selected_window_offtarget_filter") or window.get("offtarget_safety_filter") or {}
    return filt.get("passes") is True and filt.get("status") == "criterion_met"


def selected(summary):
    refs = summary.get("references")
    if refs is None:
        window = summary.get("selected_window")
        return (summary, window) if approved(summary, window) else ({}, None)
    candidates = [(r, r["selected_window"]) for r in refs if r.get("status") == "complete"
                  and approved(r, r.get("selected_window")) and r["selected_window"].get("adjusted_priority_score") is not None]
    return min(candidates, key=lambda pair: (priority_key(pair[1]), str(pair[0].get("main_target", "")))) if candidates else ({}, None)


def result_row(entry):
    row = {key: None for key in FIELDS}
    row.update({key: entry.get(key) for key in ("input_index", "gene_id", "status", "error", "output_directory")})
    ref, window = ({}, None) if entry.get("status") == "incomplete_timeout" else selected(entry.get("result") or {})
    score = (window or {}).get("adjusted_priority_score")
    evaluable = type(score) in (int, float) and math.isfinite(score)
    row["ranking_status"] = "evaluable" if evaluable else "no_evaluable_result"
    if not window:
        for slug in ("human", "maize", "bee"):
            row[slug+"_offtarget"] = "not_evaluated"
        return row
    for key in FIELDS:
        if key in window:
            row[key] = window[key]
    # Missing scores remain null: neither fallback to original nor zero imputation.
    row["adjusted_priority_score"] = score if evaluable else None
    row["reference_id"] = ref.get("main_target") or ref.get("reference_id") or window.get("main_target")
    def percent(count, total):
        return 100*count/total if count is not None and total else None
    row["efficient_percent"] = percent(window.get("is_efficient_true_count"), window.get("is_efficient_evaluable_count"))
    row["diversity_percent"] = percent(window.get("sirna_unique_sequences_count"), window.get("sirna_unique_count"))
    filt = ref.get("selected_window_offtarget_filter") or window.get("offtarget_safety_filter") or {}
    row["offtarget_filter_status"] = filt.get("status", "not_evaluated")
    species = {r.get("organism"): r for r in filt.get("by_species", [])}
    evaluations = {r.get("organism"): r for r in ref.get("selected_window_offtarget_evaluations", [])}
    for organism, slug in [("Homo sapiens", "human"), ("Zea mays", "maize"), ("Apis mellifera", "bee")]:
        evaluation = evaluations.get(organism, {})
        state = evaluation.get("status")
        row[slug+"_offtarget"] = (state if state in ("imported_unvalidated", "not_evaluated", "error")
                                 else species.get(organism, {}).get("status", "not_evaluated"))
    return row


def save_reports(root, manifest):
    rows = [result_row(e) for e in manifest["genes"]]
    ranked = sorted(rows, key=lambda r: (r["ranking_status"] != "evaluable",
                    -(r["adjusted_priority_score"] or 0) if r["ranking_status"] == "evaluable" else 0, r["input_index"]))
    rank = 0
    for row in ranked:
        if row["ranking_status"] == "evaluable":
            rank += 1
            row["rank"] = rank
    manifest["counts"] = {"total": len(rows), "completed": sum(r["status"] == "complete" for r in rows),
                          "failed": sum(r["status"] in ("failed", "invalid_input", "incomplete_timeout") for r in rows),
                          "evaluable": rank}
    for name, data in [("resultados_genes.csv", rows), ("ranking_genes.csv", ranked)]:
        with atomic_text(root / name, newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(data)
    write_json(root / "resumen_batch.json", manifest)
    return ranked


def new_batch(records):
    index = 1
    while True:
        root = Path.cwd() / "output" / f"batch_{index}"
        check_paths(root / "inputs" / "g1_1.fasta", root / "output" / f"g{len(records)}_1")
        root.parent.mkdir(parents=True, exist_ok=True)
        try:
            root.mkdir()
            return root
        except FileExistsError:
            index += 1


def run_genes(root, manifest, workflow, timeout=DEFAULT_GENE_TIMEOUT):
    for entry in manifest["genes"]:
        if entry["status"] == "invalid_input":
            continue
        # Recover an output published before a crash in the batch bookkeeping.
        if entry.get("output_directory"):
            report = root / entry["output_directory"] / "pipeline_result.json"
            if report.is_file():
                try:
                    result = read_json(report)
                    verified = not entry.get("result_sha256") or digest(report) == entry["result_sha256"]
                    if (verified and result.get("status") == "complete"
                            and entry["status"] in ("complete", "running", "interrupted")
                            and completed_worker(root, (entry.get("attempts") or [{}])[-1])):
                        entry.update(status="complete", result=result, result_sha256=digest(report), error=None)
                        if entry.get("attempts"):
                            entry["attempts"][-1]["status"] = "complete"
                        continue
                except (ValueError, OSError):
                    pass
        started = monotonic()
        print(f"Analizando gen {entry['input_index']}/{len(manifest['genes'])}: {entry['gene_id']}", flush=True)
        try:
            fasta, out = reserve_run(root, f"g{entry['input_index']}", entry["sequence"])
            entry.update(status="running", error=None, error_type=None, result=None, result_sha256=None,
                         elapsed_s=None, timeout_seconds=timeout, timeout_stage=None, timeout_reference_id=None,
                         output_directory=str(out.relative_to(root)))
            attempt = {"output_directory": entry["output_directory"], "log_file": out.name + ".log",
                       "worker_state": out.name + ".worker.json", "worker_receipt": out.name + ".done.json",
                       "timeout_seconds": timeout, "status": "running"}
            entry["attempts"].append(attempt)
            save_reports(root, manifest)
            result = run_gene_worker(workflow,
                ["--gene-id", entry["gene_id"], "--sequence-file", str(fasta),
                 "--out", str(out), "--run-sifi-headless"],
                {"input_mode":"rnai-designer-batch", "original_gene_name":entry["gene_id"],
                 "batch_input_index":entry["input_index"], "input_fasta":str(fasta)},
                root/(out.name+".log"), timeout, started=started)
            code = result.returncode
            entry.update(elapsed_s=result.elapsed_s, return_code=code)
            attempt.update(elapsed_s=result.elapsed_s, return_code=code,
                           worker_pid=result.pid, cleanup_confirmed=result.cleanup_confirmed)
            if result.timed_out:
                message = (f"Timeout de {timeout:g} s: {entry['gene_id']}; transcurrido {result.elapsed_s:.2f} s; "
                           f"etapa={result.stage or 'no disponible'}; GITV={result.reference_id or 'no disponible'}")
                entry.update(status="incomplete_timeout", error=message, error_type="GeneTimeout",
                             result=None, result_sha256=None, timeout_stage=result.stage,
                             timeout_reference_id=result.reference_id)
                attempt.update(status="incomplete_timeout", error=message, stage=result.stage,
                               reference_id=result.reference_id)
                with (root/(out.name+".log")).open('a',encoding='utf-8') as log:
                    log.write("\n[batch] "+message+"; arbol de procesos cancelado.\n")
            else:
                report = out / "pipeline_result.json"
                if report.is_file():
                    entry.update(result=read_json(report), result_sha256=digest(report))
                complete = code == 0 and (entry.get("result") or {}).get("status") == "complete"
                entry.update(status="complete" if complete else "failed",
                             error=None if complete else f"Workflow incompleto o fallido (codigo {code}); ver log del intento")
                attempt.update(status=entry["status"], error=entry["error"])
        except CleanupFailure as error:
            # Do not start another gene if tree termination could not be verified.
            entry.update(status="incomplete_timeout" if error.timed_out else "failed",
                         result=None,result_sha256=None,error=str(error),elapsed_s=error.elapsed_s)
            if entry["attempts"]:
                entry["attempts"][-1].update(status=entry["status"],cleanup_confirmed=False,error=str(error))
            manifest["status"]="partial"
            save_reports(root,manifest)
            raise
        except KeyboardInterrupt:
            entry.update(status="interrupted", error="Ejecucion interrumpida", elapsed_s=monotonic()-started)
            if entry["attempts"]:
                entry["attempts"][-1].update(status="interrupted", error=entry["error"])
            manifest["status"] = "interrupted"
            save_reports(root, manifest)
            raise
        except (Exception, SystemExit) as error:
            # A batch job is an isolation boundary: one engine/parser failure
            # must not prevent unrelated genes from running. Ctrl+C stays separate.
            entry.update(status="failed", error=str(error), error_type=type(error).__name__)
            if entry["attempts"]:
                entry["attempts"][-1].update(status="failed", error=str(error))
        save_reports(root, manifest)
    manifest["status"] = "complete" if all(e["status"] == "complete" for e in manifest["genes"]) else "partial"
    return save_reports(root, manifest)


def print_batch(manifest, ranked, root):
    counts = manifest["counts"]
    print(f"Genes: {counts['total']} | Completados: {counts['completed']} | Fallidos: {counts['failed']} | Evaluables: {counts['evaluable']}")
    timeouts = sum(e["status"] == "incomplete_timeout" for e in manifest["genes"])
    if timeouts:
        print(f"Incompletos por timeout: {timeouts} (reintentables con --resume)")
    print(NOTE)
    print("Gen                          | Ajustado /100 | nt   | Clase / estado")
    for row in ranked:
        score = f"{row['adjusted_priority_score']:.2f}" if row["ranking_status"] == "evaluable" else "sin resultado"
        print(f"{row['gene_id']:<28} | {score:>13} | {str(row['length_nt'] or '-'):>4} | {row['window_class'] or row['ranking_status']} / {row['status']}")
    print("Carpeta de salida: " + str(root))


def batch_command(argv, workflow):
    parser = argparse.ArgumentParser(prog="rnai-designer batch")
    parser.add_argument("fasta", type=Path, nargs="?")
    parser.add_argument("--resume", type=Path, help="Carpeta batch existente")
    parser.add_argument("--gene-timeout", type=float, default=None, metavar="SEGUNDOS",
                        help="Limite total por gen: 360 s por defecto; 0 desactiva el limite. En resume se conserva el valor guardado si no se indica.")
    args = parser.parse_args(argv)
    if args.gene_timeout is not None and (not math.isfinite(args.gene_timeout) or args.gene_timeout < 0):
        parser.error("--gene-timeout debe ser finito y no negativo; usar 0 para desactivar")
    if args.fasta is None and args.resume is None:
        parser.error("Se requiere genes.fasta o --resume CARPETA")
    try:
        if args.resume:
            root = args.resume.resolve()
            # Acquire lock before reading state, not after (concurrent resume).
            if not (root / "resumen_batch.json").is_file():
                raise ValueError("No existe resumen_batch.json en la carpeta indicada")
            manifest = None
        else:
            source = args.fasta.resolve()
            records = parse_fasta(source)
            root = new_batch(records)
            manifest = {"schema_version": "1.0", "source_fasta": str(source), "source_sha256": digest(source),
                        "status": "running", "interpretation": NOTE, "genes": records}
        with batch_lock(root):
            if manifest is None:
                manifest = read_json(root / "resumen_batch.json")
                if manifest.get("schema_version") != "1.0" or not isinstance(manifest.get("genes"), list):
                    raise ValueError("Manifiesto batch incompatible")
                if args.fasta and digest(args.fasta) != manifest["source_sha256"]:
                    raise ValueError("El FASTA no coincide con el batch original")
            timeout = args.gene_timeout if args.gene_timeout is not None else manifest.get("gene_timeout_seconds", DEFAULT_GENE_TIMEOUT)
            if type(timeout) not in (int,float) or not math.isfinite(timeout) or timeout < 0:
                raise ValueError("gene_timeout_seconds invalido en el resumen batch")
            manifest["gene_timeout_seconds"] = timeout
            manifest["status"] = "running"
            save_reports(root, manifest)
            try:
                ranked = run_genes(root, manifest, workflow, timeout)
            except KeyboardInterrupt:
                print('Batch interrumpido. Reanudar con: python -m rnai_designer batch --resume "' + str(root) + '"', file=sys.stderr)
                return 130
            print_batch(manifest, ranked, root)
            return 0 if manifest["counts"]["failed"] == 0 else 2
    except (OSError, ValueError, TypeError, KeyError) as error:
        print("Error batch: " + str(error), file=sys.stderr)
        return 2
