"""Run the available post-si-Fi stages as one reproducible workflow."""
import json
from pathlib import Path

from .automation import run_offtargets
from .efficiency_analysis import summarize
from .reporting import write_reports
from .sifi_adapter import provenance
from .automatic_window_selection import select_windows
from .automatic_window_reporting import write_selection_outputs
from .offtarget_assessment import overall as overall_offtarget
from .data_models import Candidate


def run_integrated_workflow(candidate, imported: dict, settings, out: Path,
                            sifi_path: Path | None = None, run_path: Path | None = None,
                            timeout: float = 300, min_length: int = 180,
                            max_length: int = 250, sirna_size: int | None = None) -> dict:
    """Run Bowtie, integrated reporting, exhaustive windows, and the final summary."""
    out.mkdir(parents=True, exist_ok=True)
    if imported.get("run") is None and run_path is not None:
        from .data_models import read_json
        imported["run"] = read_json(run_path)
    if sifi_path is not None:
        imported["source"] = provenance(sifi_path)
        inputs = out / "inputs"
        inputs.mkdir(parents=True, exist_ok=True)
        target = inputs / "designer_original.json"
        if sifi_path.resolve() != target.resolve():
            source_bytes = sifi_path.read_bytes()
            if target.exists() and target.read_bytes() != source_bytes:
                raise ValueError("Ya existe un JSON si-Fi distinto en inputs; usar una carpeta nueva")
            if not target.exists():
                target.write_bytes(source_bytes)
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
    evaluations = run_offtargets(candidate, settings, out / "offtarget_runs", timeout)
    integrated = write_reports(out / "integrated", candidate, imported,
                               summarize(imported["records"]), evaluations)
    report_path = out / "integrated" / "report.json"
    selection = select_windows(integrated, min_length, max_length, sirna_size)
    selected_window_evaluations = []
    selected_window_filter = {"status": "not_evaluated", "passes": False,
                              "interpretation": "La secuencia no se libera sin tres corridas completas y cero coincidencias exactas de 21 nt."}
    winner = selection["selected_window"]
    if winner is not None:
        scoped = Candidate.from_fields(
            f'{candidate.gene_id}_window_{winner["start_1based"]}_{winner["end_1based_inclusive"]}',
            winner["sequence"], candidate.original["reference"].get("source"),
            candidate.original["reference"].get("version"))
        if winner["length_nt"] < 21:
            selected_overall = {"status": "not_evaluated"}
            selected_window_evaluations = []
        else:
            selected_window_evaluations = run_offtargets(
                scoped, settings, out / "offtarget_selected_window", timeout)
            selected_overall = overall_offtarget(selected_window_evaluations)
        passes = (selected_overall["status"] == "criterion_met" and
                  all(item["criterion"]["status"] == "criterion_met" and
                      item["criterion"]["exact_21nt_rows"] == 0
                      for item in selected_window_evaluations) and
                  {item["organism"] for item in selected_window_evaluations} ==
                  {"Homo sapiens", "Zea mays", "Apis mellifera"})
        selected_window_filter = {
            "status": selected_overall["status"], "passes": passes,
            "criterion": "Cero coincidencias exactas de 21 nt en las tres especies; Bowtie completo, cero mismatches, sin gaps, todas las alineaciones y ambas hebras.",
            "by_species": [{"organism": item["organism"], "status": item["criterion"]["status"],
                            "exact_21nt_rows": item["criterion"]["exact_21nt_rows"],
                            "reasons": item["criterion"]["reasons"]}
                           for item in selected_window_evaluations],
            "scope": "Todos los siRNAs de 21 nt contenidos en la ventana destacada; contra las secuencias representadas por cada indice Bowtie.",
            "interpretation": "Criterio computacional acotado; no es una certificacion general de seguridad.",
        }
        winner["offtarget_safety_filter"] = selected_window_filter
        winner["sequence_release_status"] = "released_after_filter" if passes else "withheld_filter_failed_or_indeterminate"
        if not passes:
            winner.pop("sequence", None)
            winner["sequence_withheld"] = True
            selection["warnings"].append("La secuencia destacada fue retenida: el filtro off-target no acredito cero coincidencias exactas en las tres especies.")
    selection = write_selection_outputs(integrated, selection, report_path, out / "window_selection")

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
        selected_summary["score_name"] = "efficacy_priority_score"
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
    (out / "pipeline_result.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    return {"report": integrated, "selection": selection, "summary": summary}
