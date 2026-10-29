"""Run the available post-si-Fi stages as one reproducible workflow."""
from pathlib import Path

from .automation import run_offtargets
from .efficiency_analysis import summarize
from .reporting import write_reports
from .sifi_adapter import provenance
from .automatic_window_selection import select_windows
from .automatic_window_reporting import write_selection_outputs
from .safe_window_selection import choose_safe_window
from .data_models import Candidate
from .atomic_io import write_json
from time import perf_counter
import sys
import shutil
from .cache import file_digest


def run_integrated_workflow(candidate, imported: dict, settings, out: Path,
                            sifi_path: Path | None = None, run_path: Path | None = None,
                            timeout: float = 300, min_length: int = 180,
                            max_length: int = 250, sirna_size: int | None = None, length_penalty_weight: float = 30.0) -> dict:
    """Run Bowtie, integrated reporting, exhaustive windows, and the final summary."""
    out.mkdir(parents=True, exist_ok=True)
    started = perf_counter()
    stages = {}
    def progress(name):
        print('[rnai-designer] ' + name, file=sys.stderr, flush=True)
    if imported.get("run") is None and run_path is not None:
        from .data_models import read_json
        imported["run"] = read_json(run_path)
    if sifi_path is not None:
        imported["source"] = provenance(sifi_path)
        inputs = out / "inputs"
        inputs.mkdir(parents=True, exist_ok=True)
        target = inputs / "designer_original.json"
        if sifi_path.resolve() != target.resolve():
            if target.exists() and file_digest(target) != file_digest(sifi_path):
                raise ValueError("Ya existe un JSON si-Fi distinto en inputs; usar una carpeta nueva")
            if not target.exists():
                shutil.copyfile(sifi_path, target)
        imported["recovered_copy"] = "inputs/designer_original.json"
    if run_path is not None:
        imported["run_source"] = provenance(run_path)
        inputs = out / "inputs"
        inputs.mkdir(parents=True, exist_ok=True)
        run_copy = inputs / "run.json"
        if run_path.resolve() != run_copy.resolve():
            source_bytes = run_path.read_bytes()
            if run_copy.exists() and run_copy.read_bytes() != source_bytes:
                raise ValueError("Ya existe un manifiesto distinto en inputs; usar una carpeta nueva")
            if not run_copy.exists():
                run_copy.write_bytes(source_bytes)
    progress('Off-targets del gen completo')
    t = perf_counter()
    evaluations = run_offtargets(candidate, settings, out / "offtarget_runs", timeout)
    stages['full_gene_offtargets_s'] = perf_counter()-t
    progress('Reporte integrado y gráficos')
    t = perf_counter()
    integrated = write_reports(out / "integrated", candidate, imported,
                               summarize(imported["records"]), evaluations)
    report_path = out / "integrated" / "report.json"
    stages['integrated_reporting_s'] = perf_counter()-t
    progress('Enumeración exhaustiva y métricas de ventanas')
    t = perf_counter()
    selection = select_windows(integrated, min_length, max_length, sirna_size, length_penalty_weight)
    stages['window_calculation_s'] = perf_counter()-t
    t = perf_counter()
    selected_window_evaluations, selected_window_filter, offtarget_selection = choose_safe_window(
        candidate, integrated, selection, settings, out, timeout, run_offtargets, progress)
    stages['selected_window_offtargets_s'] = perf_counter()-t
    progress(f'Exportando {selection["evaluated_window_count"]} ventanas y gráficos')
    t = perf_counter()
    selection = write_selection_outputs(integrated, selection, report_path, out / "window_selection")
    stages['window_export_s'] = perf_counter()-t

    overall = integrated["offtarget_criterion"]
    by_species = [{"organism": entry["organism"], "run_status": entry["status"],
                   "criterion_status": entry["criterion"]["status"],
                   "exact_21nt_rows": entry["criterion"]["exact_21nt_rows"],
                   "label": entry["criterion"]["label"],
                   "database": (entry.get("result", {}).get("run") or entry.get("manifest", {})).get("database"),
                   "database_version": (entry.get("result", {}).get("run") or entry.get("manifest", {})).get("database_version"),
                   "database_fingerprint": (entry.get("result", {}).get("run") or entry.get("manifest", {})).get("database_fingerprint"),
                   "reasons": entry["criterion"]["reasons"]}
                  for entry in evaluations]
    winner = selection["selected_window"]
    selected_summary = None
    if winner is not None:
        selected_summary = dict(winner)
        selected_summary["score_name"] = "adjusted_priority_score"
        selected_summary["score_interpretation"] = (
            "Prioridad heuristica exploratoria basada en si-Fi, GC y diversidad; no es eficacia medida ni validada.")
    summary = {
        "schema_version": "1.0", "candidate_id": candidate.gene_id,
        "main_target": ((imported.get("run") or {}).get("main_target") or
                        ((imported.get("run") or {}).get("parameters") or {}).get("main_targets")),
        "status": "complete", "integrated_report": "integrated/report.json",
        "sifi_result": "inputs/designer_original.json" if sifi_path else imported.get("source"),
        "sifi_run": imported.get("run"),
        "green_regions": selection["green_regions"],
        "green_mask_source": selection["green_mask_source"],
        "green_mask_provenance": imported.get("green_mask_provenance"),
        "evaluated_window_count": selection["evaluated_window_count"],
        "selected_window": selected_summary,
        "selection_status": selection["selection_status"],
        "offtarget_selection": offtarget_selection,
        "selected_window_offtarget_filter": selected_window_filter,
        "selected_window_offtarget_evaluations": [
            {"organism": item["organism"], "status": item["status"],
             "criterion": item["criterion"], "run": item.get("result", {}).get("run") or item.get("manifest")}
            for item in selected_window_evaluations],
        "offtarget_21nt_filter": {"overall": overall, "by_species": by_species,
                                  "passes_all_species": overall["status"] == "criterion_met",
                                  "scope": "Todas las ventanas de 21 nt de la secuencia de referencia, no solo la ventana destacada.",
                                  "interpretation": "Solo criterio computacional de ausencia de coincidencias exactas 21nt bajo el alcance documentado; no declara seguridad."},
        "artifacts": {"selection_report": "window_selection/selection_report.json",
                      "all_windows": "window_selection/all_windows.csv",
                      "selected_window": "window_selection/selected_window.json" if winner else None,
                      "selected_visualization": "window_selection/selected_window.svg"},
        "warnings": list(dict.fromkeys(integrated["warnings"] + selection["warnings"])),
    }
    stages['total_post_sifi_s'] = perf_counter()-started
    write_json(out / 'performance.json', stages)
    write_json(out / 'pipeline_result.json', summary)
    return {"report": integrated, "selection": selection, "summary": summary}
