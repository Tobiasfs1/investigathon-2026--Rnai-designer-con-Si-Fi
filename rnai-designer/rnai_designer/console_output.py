"""Dependency-free workflow console report, also accepts legacy summaries."""
from pathlib import Path
import sys
from .priority_scoring import priority_key


def _print_detailed(summary, out, sequence=None, stream=None):
    stream = stream or sys.stdout
    def table(rows):
        for label, value in rows:
            if isinstance(value, float):
                value = f"{value:.4f}"
            print(f"{label:<28} | {value if value is not None else 'no disponible'}", file=stream)
    refs = summary.get("references")
    if refs is None:
        refs = [summary]
    complete = [r for r in refs if r.get("status", "complete") == "complete"]
    table([("Gen", summary.get("candidate_id")), ("Estado final", summary.get("status")),
           ("Referencias completas", summary.get("completed_reference_count", len(complete))),
           ("Referencias fallidas", summary.get("failed_reference_count", len(refs)-len(complete))),
           ("Ventanas evaluadas", summary.get("combined_window_count", summary.get("evaluated_window_count")))])
    candidates = [(r, r["selected_window"]) for r in complete if r.get("selected_window") and
                  r["selected_window"].get("adjusted_priority_score", r["selected_window"].get("efficacy_priority_score")) is not None]
    for r in refs:
        table([("Referencia", r.get("main_target", r.get("reference_id"))), ("Estado referencia", r.get("status")), ("Error", r.get("error"))])
    if not candidates:
        print("Sin ganador: no hay ventanas elegibles con puntaje y datos suficientes en referencias completas.", file=stream)
    else:
        r, w = min(candidates, key=lambda pair: (priority_key(pair[1]), str(pair[0].get("main_target", ""))))
        table([("Referencia ganadora", r.get("main_target", w.get("reference_source"))),
               ("Coordenadas (1-based)", f"{w.get('start_1based')}-{w.get('end_1based_inclusive')}"),
               ("Longitud nt", w.get("length_nt")), ("Clase", w.get("window_class")),
               ("Puntaje original", w.get("efficacy_priority_score")), ("Penalizacion", w.get("length_penalty")),
               ("Puntaje ajustado", w.get("adjusted_priority_score")), ("Calidad de evidencia", w.get("quality_score")),
               ("Suficiencia de datos", w.get("data_sufficiency_status")), ("GC %", w.get("gc_percent")),
               ("siRNAs totales", w.get("sirna_record_count")), ("siRNAs unicos", w.get("sirna_unique_count")),
               ("Fraccion eficientes si-Fi", w.get("is_efficient_true_fraction_among_evaluable"))])
        filt = r.get("selected_window_offtarget_filter", w.get("offtarget_safety_filter", {}))
        table([("Off-target ventana", filt.get("status", "not_evaluated")), ("Cumple filtro", filt.get("passes"))])
        for species in filt.get("by_species", []):
            table([(species.get("organism", "Especie"), species.get("status")), ("Coincidencias exactas 21nt", species.get("exact_21nt_rows")), ("Razones", species.get("reasons"))])
        seq = w.get("sequence") if _sequence_authorized(r, w) else None
        if _sequence_authorized(r, w) and not seq and sequence and w.get("start_1based") and w.get("end_1based_inclusive"):
            seq = sequence[w["start_1based"]-1:w["end_1based_inclusive"]]
        print("Secuencia ganadora (ver estado del filtro):", file=stream)
        print(seq or "no disponible o retenida por el filtro", file=stream)
    print("Estados off-target: criterion_met=cumple criterio; failed=no cumple; not_evaluated=sin evaluar; imported_unvalidated=importado sin validar.", file=stream)
    print("Limitaciones: prioridad heuristica, no eficacia validada; filtro computacional acotado, no certificacion de seguridad. Regiones cortas y convencionales conservan su clase y no son equivalentes.", file=stream)
    for r in refs:
        for warning in r.get("warnings", []):
            print("Advertencia: " + str(warning), file=stream)
        safety = r.get("offtarget_21nt_filter", {})
        if safety:
            table([("Off-target gen completo", safety.get("overall", {}).get("status"))])
            for item in safety.get("by_species", []):
                table([(item.get("organism", "Especie"), item.get("criterion_status")), ("Razones", item.get("reasons"))])
        base = Path(out) / r.get("output_directory", "")
        for name in ("pipeline_result.json", "window_selection/selection_report.json", "window_selection/all_windows.csv"):
            print("Reporte: " + str(base / name), file=stream)
    if summary.get("references") is not None:
        for name in ("reference_comparison.json", "all_references_windows.csv"):
            print("Reporte: " + str(Path(out) / name), file=stream)


def _filter(reference, window):
    return reference.get("selected_window_offtarget_filter") or window.get("offtarget_safety_filter") or {}


def _sequence_authorized(reference, window):
    filt = _filter(reference, window)
    return (filt.get("passes") is True and filt.get("status") == "criterion_met"
            and not window.get("sequence_withheld")
            and window.get("sequence_release_status") in (None, "released_after_filter"))


def print_workflow_summary(summary, out, sequence=None, stream=None, verbose=False):
    """Render only; never infer missing scores or alter persisted report data."""
    if verbose:
        return _print_detailed(summary, out, sequence, stream)
    import textwrap
    stream = stream or sys.stdout
    def emit(label, value):
        print(f"{label:<24} : {value}", file=stream)
    def number(value):
        return "no disponible" if value is None else f"{value:.2f}"
    def ratio(count, total):
        percent = number(100 * count / total) if count is not None and total else "no disponible"
        return f"{count if count is not None else '?'} / {total if total is not None else '?'} ({percent}%)"
    refs = summary.get("references")
    refs = [summary] if refs is None else refs
    candidates = [(r, r["selected_window"]) for r in refs
                  if r.get("status", "complete") == "complete" and r.get("selected_window")
                  and r["selected_window"].get("adjusted_priority_score", r["selected_window"].get("efficacy_priority_score")) is not None]
    emit("Gen", summary.get("candidate_id", "no disponible"))
    emit("Estado", summary.get("status", "no disponible"))
    emit("Referencias procesadas", summary.get("reference_count", len(refs)))
    emit("Ventanas evaluadas", summary.get("combined_window_count", summary.get("evaluated_window_count", "no disponible")))
    warnings = ["Puntuacion heuristica; no representa eficacia experimental."]
    if candidates:
        r, w = min(candidates, key=lambda pair: (priority_key(pair[1]), str(pair[0].get("main_target", ""))))
        print("Ventana destacada", file=stream)
        emit("Referencia ganadora", r.get("main_target") or r.get("reference_id") or w.get("main_target") or w.get("reference_source") or "no disponible")
        emit("Inicio - fin (1-based)", f"{w.get('start_1based', '?')} - {w.get('end_1based_inclusive', '?')}")
        emit("Longitud", f"{w.get('length_nt', '?')} nt")
        emit("Clase", w.get("window_class", "no disponible"))
        emit("Puntaje ajustado", number(w.get("adjusted_priority_score")) + " /100")
        emit("Puntaje original", number(w.get("efficacy_priority_score")) + " /100")
        emit("Penalizacion longitud", number(w.get("length_penalty")))
        emit("Calidad de evidencia", number(w.get("quality_score")) + " /100")
        emit("GC", number(w.get("gc_percent")) + "%")
        emit("siRNAs eficientes si-Fi", ratio(w.get("is_efficient_true_count"), w.get("is_efficient_evaluable_count")) + " evaluables")
        emit("Diversidad siRNAs", ratio(w.get("sirna_unique_sequences_count"), w.get("sirna_unique_count")))
        seq = None
        if _sequence_authorized(r, w):
            seq = w.get("sequence")
            if not seq and sequence and w.get("start_1based") and w.get("end_1based_inclusive"):
                seq = sequence[w["start_1based"]-1:w["end_1based_inclusive"]]
        print("Secuencia:" if seq else "Secuencia: no disponible o retenida por el filtro.", file=stream)
        if seq:
            print("\n".join(textwrap.wrap(seq, 60)), file=stream)
        filt = _filter(r, w)
        species = {x.get("organism"): x for x in filt.get("by_species", [])}
        evaluations = {x.get("organism"): x for x in r.get("selected_window_offtarget_evaluations", [])}
        labels = {"criterion_met": "criterio cumplido", "failed": "criterio no cumplido",
                  "not_evaluated": "no evaluado", "imported_unvalidated": "importado sin validar",
                  "error": "evaluacion fallida"}
        for organism, label in [("Homo sapiens", "Humano"), ("Zea mays", "Maiz"), ("Apis mellifera", "Abeja")]:
            item = species.get(organism, {})
            evaluation = evaluations.get(organism, {})
            status = evaluation.get("status") if evaluation.get("status") in ("imported_unvalidated", "not_evaluated", "error") else item.get("status", "not_evaluated")
            emit("Off-target " + label, labels.get(status, "no validado (" + str(status) + ")"))
        if isinstance(w.get("length_nt"), (int, float)) and w["length_nt"] < 180:
            warnings.insert(0, "Ventana menor de 180 nt; region corta.")
        extra = r.get("warnings", [])
    else:
        if any((r.get("offtarget_selection") or {}).get("status") == "no_window_passed" for r in refs):
            print("Sin ganador: ninguna ventana acredito el filtro off-target en las tres especies.", file=stream)
        else:
            print("Sin ganador: no hubo resultado elegible en referencias completas.", file=stream)
        extra = summary.get("warnings", [])
    if summary.get("status") in ("partial", "error", "no_reference_completed"):
        warnings[0] += " Ejecucion incompleta; revisar reportes."
    for warning in extra:
        clean = " ".join(str(warning).split())
        if clean and clean not in warnings:
            warnings.append(clean)
    for warning in list(dict.fromkeys(warnings))[:2]:
        print("Advertencia: " + warning, file=stream)
    emit("Carpeta de salida", Path(out))
