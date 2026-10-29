"""Interpretación descriptiva de evidencia de alineamientos de herramientas externas."""
import re

from .bowtie_adapter import reverse_complement, sequence_hash
from .data_models import Candidate

LABELS = {"criterion_met": "Criterio computacional cumplido",
          "exact_match_detected": "Coincidencia exacta detectada", "indeterminate": "Indeterminado",
          "not_evaluated": "No evaluado", "error": "Error"}


def mismatch_count(value, length: int) -> int | None:
    if type(value) is int and value == 0:
        return 0
    if value == "":
        return 0
    if not isinstance(value, str):
        return None
    descriptors = value.split(",")
    positions = []
    for descriptor in descriptors:
        match = re.fullmatch(r"([0-9]+):([ACGTN])>([ACGTN])", descriptor, re.IGNORECASE)
        if not match or int(match[1]) >= length or match[2].upper() == match[3].upper():
            return None
        positions.append(int(match[1]))
    return len(positions) if len(set(positions)) == len(positions) else None


def describe_alignment(row: dict, candidate: Candidate) -> dict:
    sequence = row.get("sirna_sequence", "").upper()
    length = len(sequence)
    count = mismatch_count(row.get("mismatches"), length)
    start, reference_start = row.get("sirna_position"), row.get("reference_strand_pos")
    result = {"reference_id": row.get("hit_name"), "query_start_1based": start,
              "reference_start_0based": reference_start, "strand": row.get("strand"),
              "length": length, "mismatch_count": count, "match_type": "unverified",
              "contains_exact_21nt": None}
    if type(start) is not int or start < 1 or start+length-1 > len(candidate.sequence):
        return result
    expected = candidate.sequence[start-1:start-1+length].upper()
    if row.get("strand") == "-":
        expected = reverse_complement(expected)
    reference_id = row.get("hit_name")
    if (not isinstance(reference_id, str) or not reference_id or type(reference_start) is not int
            or reference_start < 0 or row.get("strand") not in ("+", "-")
            or not sequence or set(sequence) - set("ACGT") or sequence != expected
            or row.get("cigar") not in (None, f"{length}M", f"{length}=")):
        return result
    result["query_end_1based_inclusive"] = start + length - 1
    result["reference_end_0based_exclusive"] = reference_start + length
    if count == 0:
        result.update(match_type="exact_21nt" if length == 21 else "exact_other_length",
                      contains_exact_21nt=length >= 21)
    elif count is not None:
        result.update(match_type="mismatch_alignment", contains_exact_21nt=False if length <= 21 else None)
    return result


def assess(entry: dict, candidate: Candidate) -> dict:
    imported = entry.get("result", {})
    run = imported.get("run") or entry.get("manifest", {})
    alignments = [describe_alignment(row, candidate) for row in imported.get("records", [])]
    exact = sum(a["contains_exact_21nt"] is True for a in alignments)
    mismatched = sum(a["match_type"] == "mismatch_alignment" for a in alignments)
    reasons = []
    status = "indeterminate"
    if entry.get("status") == "not_evaluated":
        status = "not_evaluated"
    elif exact:
        status = "exact_match_detected"
    elif entry.get("status") == "error":
        status = "error"
        reasons.append(entry.get("error", "La evaluación no terminó correctamente"))
    else:
        params = run.get("parameters") if isinstance(run.get("parameters"), dict) else {}
        completion = run.get("completion") if isinstance(run.get("completion"), dict) else {}
        coverage = run.get("coverage") if isinstance(run.get("coverage"), dict) else {}
        expected = max(0, len(candidate.sequence) - 20)
        checks = {
            "Se requiere Bowtie original completo; el JSON de si-Fi filtra la hebra negativa":
                entry.get("results_format") == "bowtie_legacy_tab",
            "La ejecución no acredita finalización correcta":
                completion.get("status") == "completed" and type(completion.get("returncode")) is int
                and completion["returncode"] == 0,
            "Falta búsqueda exhaustiva de 21nt, cero mismatches, sin gaps y ambas hebras":
                type(params.get("sirna_size")) is int and params["sirna_size"] == 21
                and type(params.get("mismatches")) is int and params["mismatches"] == 0
                and params.get("alignment_policy") == "-v" and params.get("gapped") is False
                and params.get("all_alignments") is True and params.get("both_strands") is True,
            "No se acredita procesamiento de todas las ventanas de la query":
                expected > 0 and all(type(coverage.get(key)) is int and coverage[key] == expected
                                     for key in ("submitted_windows", "processed_windows", "expected_windows")),
            "La query contiene bases ambiguas o no coincide su hash":
                not (set(candidate.sequence.upper()) - set("ACGT"))
                and run.get("query_sha256") == sequence_hash(candidate),
            "El archivo de resultados no coincide con el hash de la ejecución":
                bool(imported.get("source", {}).get("sha256"))
                and run.get("result_sha256") == imported.get("source", {}).get("sha256"),
            "Hay alineamientos de interpretación insuficiente o longitud distinta de 21nt":
                all(a["contains_exact_21nt"] is False and a["length"] == 21 for a in alignments),
        }
        reasons.extend(reason for reason, valid in checks.items() if not valid)
        if not reasons:
            status = "criterion_met"
    return {"status": status, "label": LABELS[status], "exact_21nt_rows": exact,
            "mismatch_rows": mismatched, "unverified_rows": sum(a["match_type"] == "unverified" for a in alignments),
            "alignments": alignments, "reasons": reasons,
            "scope": "Coincidencias en las secuencias representadas por la base y bajo los parámetros registrados",
            "limitation": "Este criterio computacional no garantiza seguridad biológica ni ambiental"}


def overall(evaluations: list[dict]) -> dict:
    statuses = [entry["criterion"]["status"] for entry in evaluations]
    if "exact_match_detected" in statuses:
        status = "exact_match_detected"
    elif statuses and all(s == "criterion_met" for s in statuses):
        status = "criterion_met"
    elif statuses and all(s == "not_evaluated" for s in statuses):
        status = "not_evaluated"
    elif "error" in statuses:
        status = "error"
    else:
        status = "indeterminate"
    return {"status": status, "label": LABELS[status], "safety_conclusion": "undetermined"}
