"""Organiza evaluaciones externas; no ejecuta búsquedas ni infiere seguridad."""
from pathlib import Path
from .data_models import Candidate
from .data_models import read_json, require_text
from .sifi_adapter import import_sifi, provenance
from .bowtie_adapter import parse_bowtie
from .offtarget_assessment import assess

SPECIES = ("Homo sapiens", "Zea mays", "Apis mellifera")


def load_evaluations(path: Path | None, candidate: Candidate, tolerant: bool = False) -> list[dict]:
    if path is None:
        entries = [{"organism": name, "status": "not_evaluated"} for name in SPECIES]
        for entry in entries:
            entry["criterion"] = assess(entry, candidate)
        return entries
    manifest = read_json(path)
    if not isinstance(manifest, dict) or manifest.get("schema_version") not in ("1.0", "1.1"):
        raise ValueError("Off-target: schema_version debe ser 1.0 o 1.1")
    entries = manifest.get("evaluations")
    if not isinstance(entries, list):
        raise ValueError("evaluations debe ser lista")
    result, seen = [], set()
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("Evaluación debe ser objeto")
        name = require_text(entry, "organism")
        if name not in SPECIES or name in seen:
            raise ValueError("Organismo desconocido o repetido")
        seen.add(name)
        for key in ("tool", "tool_version", "database", "database_version", "source", "query_id"):
            require_text(entry, key)
        if not isinstance(entry.get("parameters"), dict):
            raise ValueError("parameters debe ser objeto (puede estar vacío si se desconoce)")
        relative = Path(require_text(entry, "results_file"))
        if relative.is_absolute():
            raise ValueError("results_file debe ser relativo al manifiesto")
        result_path = path.parent / relative
        item = {"organism": name, "status": "imported_unvalidated", "manifest": entry,
                "manifest_source": provenance(path), "safety_conclusion": "undetermined",
                "results_format": entry.get("results_format", "sifi_json")}
        try:
            if item["results_format"] == "sifi_json":
                imported = import_sifi(result_path, candidate, entry)
            elif item["results_format"] == "bowtie_legacy_tab":
                imported = parse_bowtie(result_path, candidate, entry)
            else:
                raise ValueError(f"Formato no soportado: {item['results_format']}")
            item.update(result=imported, hit_rows=len(imported["records"]))
        except (OSError, ValueError) as error:
            if not tolerant and entry.get("status") != "error":
                raise
            item.update(status="error", error=str(error))
        if entry.get("status") == "error":
            item.update(status="error", error=entry.get("error") or "Proceso externo incompleto")
        item["criterion"] = assess(item, candidate)
        result.append(item)
    result.extend({"organism": name, "status": "not_evaluated"} for name in SPECIES if name not in seen)
    for item in result:
        if "criterion" not in item:
            item["criterion"] = assess(item, candidate)
    return result
