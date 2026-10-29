"""Importa exclusivamente la lista JSON verificada en SifiPipeline.data_to_json."""
from pathlib import Path
import hashlib
from .data_models import Candidate, ImportedResult, read_json, require_text

BOOL_FIELDS = ("is_efficient", "strand_selection", "end_stability", "target_site_accessibility", "is_off_target")
FIELDS = ("query_name", "sirna_name", "sirna_position", "sirna_sequence", *BOOL_FIELDS,
          "accessibility_value", "hit_name", "reference_strand_pos", "strand", "mismatches")


def provenance(path: Path | str) -> dict[str, str]:
    path = Path(path)
    return {"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def import_sifi(path: Path | str, candidate: Candidate, run: dict | None = None) -> ImportedResult:
    rows = read_json(path)
    if not isinstance(rows, list):
        raise ValueError("si-Fi: se esperaba una lista JSON; CSV resumen no contiene posiciones")
    warnings = ["Predicción computacional: no demuestra eficacia biológica ni seguridad ambiental.",
                "La copia inspeccionada exporta solo hebra +; ausencia de hits no demuestra ausencia de off-targets."]
    if not run:
        warnings.append("Falta manifiesto: parámetros, base y versión de la ejecución desconocidos.")
    else:
        for key in ("query_id", "tool", "tool_version", "source", "database", "database_version"):
            if run.get(key) is None:
                warnings.append(f"Manifiesto: falta {key}")
            else:
                require_text(run, key)
        if run.get("parameters") is None:
            warnings.append("Manifiesto: parameters desconocidos")
        elif not isinstance(run["parameters"], dict):
            raise ValueError("Manifiesto: parameters debe ser objeto")
        elif not run["parameters"]:
            warnings.append("Manifiesto: parameters vacío; configuración desconocida")
    if candidate.original["metadata"].get("simulated") or ((run or {}).get("parameters") or {}).get("simulated"):
        warnings.append("DATOS FICTICIOS: no son resultados de una ejecución real ni evidencia experimental.")
    if not rows:
        warnings.append("Resultado vacío: no se interpreta como evaluación negativa ni como seguridad.")
    query_id = (run or {}).get("query_id") or candidate.gene_id
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            raise ValueError(f"Registro {index}: debe ser objeto")
        for key in FIELDS:
            if key not in row or row[key] is None:
                warnings.append(f"Registro {index}: {key} ausente o null")
        require_text(row, "sirna_name")
        sequence = require_text(row, "sirna_sequence")
        if set(sequence.upper()) - set("ACGTURYSWKMBDHVN"):
            raise ValueError(f"Registro {index}: sirna_sequence no es IUPAC")
        start = row.get("sirna_position")
        if type(start) is not int or start < 1 or start + len(sequence) - 1 > len(candidate.sequence):
            raise ValueError(f"Registro {index}: intervalo fuera de la referencia")
        if sequence.upper() != candidate.sequence[start-1:start-1+len(sequence)].upper():
            warnings.append(f"Registro {index}: secuencia no coincide con referencia")
        if row.get("query_name") != query_id:
            raise ValueError(f"Registro {index}: query_name no coincide con query_id/gene_id")
        for key in BOOL_FIELDS:
            if row.get(key) is not None and type(row[key]) is not bool:
                raise ValueError(f"Registro {index}: {key} debe ser booleano o null")
        value = row.get("accessibility_value")
        if value is not None and (type(value) not in (int, float) or not 0 <= value <= 1):
            raise ValueError(f"Registro {index}: accessibility_value debe estar entre 0 y 1")
        if row.get("is_efficient") is False and all(row.get(k) is None for k in BOOL_FIELDS[1:4]):
            warnings.append(f"Registro {index}: eficiencia posiblemente desactivada; False no equivale a ineficacia")
        pos = row.get("reference_strand_pos")
        if pos is not None and (type(pos) is not int or pos < 0):
            raise ValueError(f"Registro {index}: reference_strand_pos inválida")
        if row.get("strand") not in (None, "+", "-"):
            raise ValueError(f"Registro {index}: strand inválida")
        if row.get("strand") == "-":
            warnings.append(f"Registro {index}: hebra - inesperada para la copia inspeccionada; verificar origen")
    sizes = {len(row["sirna_sequence"]) for row in rows}
    if len(sizes) == 1:
        size = next(iter(sizes))
        expected = len(candidate.sequence) - size + 1
        observed = len({row["sirna_position"] for row in rows})
        if observed < expected:
            warnings.append(f"Solo {observed} de {expected} inicios posibles están representados; "
                            "posiciones ausentes no son mediciones de eficiencia cero.")
    elif len(sizes) > 1:
        warnings.append("Resultado con longitudes de siRNA distintas; verificar si combina ejecuciones.")
    return {"source": provenance(path), "run": run, "records": rows, "warnings": warnings}
