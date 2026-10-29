"""Búsquedas externas con el Bowtie de si-Fi. No invoca su GUI ni predice eficiencia."""
from pathlib import Path
import json
import re
import subprocess
from .cache import file_digest
from .reuse import search_cache, search_key, restore_search, remember_search
from .atomic_io import write_json

from .bowtie_adapter import parse_bowtie, sequence_hash, write_queries
from .config import Settings, index_files
from .data_models import Candidate
from .external import execute
from .offtarget_assessment import assess
from .sifi_adapter import provenance

SPECIES_SLUGS = {"Homo sapiens": "human", "Zea mays": "maize", "Apis mellifera": "bee"}


def hash_file(path: Path) -> str:
    return file_digest(path)


def database_fingerprint(settings: Settings, organism: str) -> dict:
    binding = settings.databases[organism]
    files = index_files(settings.database_dir, binding["prefix"])
    missing = [p.name for p in files if not p.is_file() or p.stat().st_size == 0]
    if missing:
        raise ValueError(f"Índice ausente/incompleto de {organism}: {', '.join(missing)}")
    return {"prefix": binding["prefix"], "declared_source": binding.get("source"),
            "declared_version": binding.get("version"),
            "identity_note": "Organismo asociado por configuración; versión biológica no inferida de las fechas",
            "files": [{"name": p.name, "size": p.stat().st_size, "mtime_ns": p.stat().st_mtime_ns,
                       "sha256": hash_file(p)} for p in files]}


def run_offtargets(candidate: Candidate, settings: Settings, out: Path, timeout: float = 300) -> list[dict]:
    root = out.resolve()
    if root.is_relative_to(settings.database_dir.resolve()):
        raise ValueError("La salida debe estar fuera del directorio de bases originales")
    if root.exists() and any(root.iterdir()):
        raise ValueError("La carpeta de búsqueda ya contiene archivos; usar una salida nueva para conservar evidencia")
    root.mkdir(parents=True, exist_ok=True)
    cache = search_cache.get()
    key = None
    if cache is not None:
        try:
            key = search_key(candidate, settings, database_fingerprint)
            if key in cache:
                restored = restore_search(cache[key], root, key)
                if restored is not None:
                    cache.move_to_end(key)
                    return restored
        except (OSError, ValueError):
            pass  # Missing/invalid databases keep the original per-species error behavior.
    evaluations, manifest_entries = [], []
    for organism, slug in SPECIES_SLUGS.items():
        directory = root / slug
        directory.mkdir(parents=True, exist_ok=True)
        entry = {"organism": organism, "status": "error", "results_format": "bowtie_legacy_tab",
                 "safety_conclusion": "undetermined"}
        run = {"organism": organism, "tool": "Bowtie (motor externo utilizado por si-Fi)",
               "tool_version": "unknown", "database": settings.databases[organism]["prefix"],
               "database_version": settings.databases[organism].get("version") or "unknown",
               "source": "rnai-designer: proceso externo, sin ejecutar RNAi Designer",
               "query_kind": "synthetic" if candidate.original["metadata"].get("simulated") else "provided_reference",
               "query_id": candidate.gene_id, "query_sha256": sequence_hash(candidate),
               "parameters": {"sirna_size": 21, "mismatches": 0, "alignment_policy": "-v",
                              "gapped": False, "all_alignments": True, "both_strands": True}}
        try:
            if settings.bowtie is None:
                raise ValueError("No se encontró bowtie; configurar ruta del ejecutable")
            run["database_fingerprint"] = database_fingerprint(settings, organism)
            version = subprocess.run([str(settings.bowtie), "--version"], capture_output=True,
                                     text=True, errors="replace", timeout=10,
                                     creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            if version.returncode != 0:
                raise ValueError("Bowtie --version falló")
            version_text = version.stdout + version.stderr
            detected = re.search(r"\bversion\s+([0-9.]+)", version_text)
            if not detected or "bowtie" not in version_text.lower():
                raise ValueError("El ejecutable no identifica una versión compatible de Bowtie")
            run["tool_version"] = detected[1]
            count = write_queries(directory / "queries.fasta", candidate)
            target = directory / "alignments.tsv"
            command = [str(settings.bowtie), "-a", "-v", "0", "-y",
                       str(settings.database_dir / run["database"]), "-f",
                       str(directory / "queries.fasta"), str(target)]
            # Cada corrida tiene una carpeta nueva; un archivo previo nunca cuenta como resultado nuevo.
            if target.exists():
                raise ValueError("Ya existe alignments.tsv: usar un directorio de salida nuevo")
            run["completion"] = execute(command, directory, timeout)
            stderr = (directory / "stderr.log").read_text(encoding="utf-8", errors="replace")
            processed = re.search(r"# reads processed:\s*([0-9]+)", stderr)
            run["coverage"] = {"expected_windows": len(candidate.sequence)-20, "submitted_windows": count,
                               "processed_windows": int(processed[1]) if processed else None}
            run["logs"] = {name: provenance(directory / name) for name in ("stdout.log", "stderr.log")}
            if target.exists():
                run["result_sha256"] = hash_file(target)
                entry["result"] = parse_bowtie(target, candidate, run)
                entry["hit_rows"] = len(entry["result"]["records"])
            if run["completion"]["status"] != "completed":
                raise ValueError(f"Bowtie no terminó correctamente: {run['completion']}")
            if not target.is_file():
                raise ValueError("Bowtie no generó archivo de resultados")
            for file in run["database_fingerprint"]["files"]:
                stat = (settings.database_dir / file["name"]).stat()
                if stat.st_size != file["size"] or stat.st_mtime_ns != file["mtime_ns"]:
                    raise ValueError("La base cambió durante la ejecución")
            entry["status"] = "imported_unvalidated"
        except (OSError, ValueError, subprocess.TimeoutExpired) as error:
            entry["error"] = str(error)
        entry["manifest"] = run
        entry["criterion"] = assess(entry, candidate)
        write_json(directory / 'run.json', run)
        manifest_entries.append({**run, "status": entry["status"], "error": entry.get("error"),
                                 "results_format": "bowtie_legacy_tab", "results_file": f"{slug}/alignments.tsv"})
        evaluations.append(entry)
    write_json(root / 'offtargets.json', {'schema_version': '1.1', 'evaluations': manifest_entries})
    if cache is not None:
        remember_search(cache, key, root, evaluations)
    return evaluations
