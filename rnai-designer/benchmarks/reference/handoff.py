"""Entrada para la GUI legada y recuperación inequívoca de JSON existentes."""
from pathlib import Path
import hashlib
import json
import tempfile

from .data_models import Candidate, read_json
from .sifi_adapter import import_sifi, provenance


def prepare(candidate: Candidate, out: Path) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    alias = "rnai_" + hashlib.sha256((candidate.gene_id + candidate.sequence).encode("utf-8")).hexdigest()[:16]
    (out / "query.fasta").write_text(f">{alias}\n{candidate.sequence}\n", encoding="utf-8")
    request = {"schema_version": "1.0", "query_id": alias, "gene_id": candidate.gene_id,
               "tool": "si-Fi", "tool_version": None, "database": None, "database_version": None,
               "source": "GUI de si-Fi; completar después de ejecutar", "parameters": {},
               "status": "not_evaluated",
               "instructions": ["Abrir query.fasta en RNAi Designer y seleccionar explícitamente la base Dalbulus",
                                "Seleccionar los main targets; registrar base, parámetros y fecha en run.json",
                                "Ejecutar collect o importar el JSON con --run run.json; no usar el CSV resumen"]}
    # No reemplazar un manifiesto que el usuario haya completado.
    if not (out / "run.json").exists():
        (out / "run.json").write_text(json.dumps(request, ensure_ascii=False, indent=2), encoding="utf-8")
    return {"directory": str(out), "query_id": alias, "designer_status": "manual_required"}


def collect(candidate: Candidate, run: dict | None, search_dirs: list[Path] | None = None):
    directories = search_dirs or [Path(tempfile.gettempdir())]
    matches, skipped = {}, 0
    for directory in directories:
        if not directory.is_dir():
            raise ValueError(f"Directorio de búsqueda no disponible: {directory}")
        for path in sorted(directory.glob("*.json")):
            # El productor escribe archivos en la raíz temporal. No recorrer otras aplicaciones ni leer enormes JSON.
            try:
                if path.stat().st_size > 100 * 1024 * 1024:
                    skipped += 1
                    continue
                rows = read_json(path)
                if not isinstance(rows, list) or not rows:
                    continue
                imported = import_sifi(path, candidate, run)
                if any("secuencia no coincide con referencia" in warning for warning in imported["warnings"]):
                    continue
                matches.setdefault(provenance(path)["sha256"], (path, imported))
            except (OSError, ValueError, TypeError, KeyError):
                continue
    if not matches:
        raise ValueError(f"No hay JSON no vacío compatible con query y referencia; archivos >100MiB omitidos: {skipped}")
    if len(matches) > 1:
        names = [str(path) for path, _ in matches.values()]
        raise ValueError(f"Varios resultados distintos compatibles; indicar --sifi explícitamente: {names}")
    path, imported = next(iter(matches.values()))
    imported["warnings"].append("Recuperado automáticamente por ID/secuencia; completar y verificar parámetros de ejecución.")
    return path, imported
