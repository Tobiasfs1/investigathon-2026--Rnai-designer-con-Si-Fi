"""Choose the highest-ranked window that meets the existing off-target criterion."""
from pathlib import Path
import subprocess

from .data_models import Candidate
from .offtarget_assessment import overall
from .priority_scoring import priority_key
from .window_explorer import _overlapping_offtargets

SPECIES = {"Homo sapiens", "Zea mays", "Apis mellifera"}
POLICY = "ranked_until_offtarget_pass_v1"


def evaluate_filter(evaluations):
    status = overall(evaluations)["status"]
    passes = (status == "criterion_met" and len(evaluations) == 3
              and {item["organism"] for item in evaluations} == SPECIES
              and all(item.get("status") not in ("error", "not_evaluated")
                      and item["criterion"]["status"] == "criterion_met"
                      and item["criterion"]["exact_21nt_rows"] == 0 for item in evaluations))
    return {"status": status if passes or status != "criterion_met" else "indeterminate", "passes": passes,
            "criterion": "Cero coincidencias exactas de 21 nt en las tres especies; Bowtie completo, cero mismatches, sin gaps, todas las alineaciones y ambas hebras.",
            "by_species": [{"organism": item["organism"], "status": item["criterion"]["status"],
                            "run_status": item.get("status"),
                            "exact_21nt_rows": item["criterion"]["exact_21nt_rows"],
                            "reasons": item["criterion"].get("reasons", [])} for item in evaluations],
            "scope": "Todos los siRNAs de 21 nt contenidos en esta ventana, contra los indices documentados.",
            "interpretation": "Criterio computacional acotado; no es una certificacion general de seguridad."}


def choose_safe_window(candidate, integrated, selection, settings, out, timeout, run_offtargets, progress):
    windows = selection["windows"]
    if hasattr(windows, 'ranked_indices'):
        eligible_count = windows.eligible_count()
        indices = windows.ranked_indices()
    else:
        indices = [i for i, row in enumerate(windows) if row.get("adjusted_priority_score", row.get("efficacy_priority_score")) is not None]
        indices.sort(key=lambda i: priority_key(windows[i]))
        eligible_count = len(indices)
    original = selection["selected_window"]
    selection["selected_window"] = None
    attempts = []
    last_evaluations = []
    last_filter = {"status": "not_evaluated", "passes": False, "by_species": []}
    for rank, index in enumerate(indices, 1):
        row = windows[index]
        start, end = row["start_1based"], row["end_1based_inclusive"]
        sequence = candidate.sequence[start-1:end]
        directory = Path(out) / ("offtarget_selected_window" if rank == 1 else f"offtarget_w{rank}")
        evaluations = []
        filt = {"status": "not_evaluated", "passes": False, "by_species": [],
                "reasons": ["Ventana menor de 21 nt: no contiene siRNAs completos para este filtro."]}
        if row["length_nt"] >= 21:
            progress(f'Filtro off-target: ventana {rank}/{eligible_count}')
            try:
                scoped = Candidate.from_fields(f"{candidate.gene_id}_window_{start}_{end}", sequence,
                           candidate.original["reference"].get("source"), candidate.original["reference"].get("version"))
                scoped.original["metadata"]["simulated"] = bool(candidate.original["metadata"].get("simulated"))
                evaluations = run_offtargets(scoped, settings, directory, timeout)
                filt = evaluate_filter(evaluations)
            except (OSError, ValueError, TypeError, KeyError, subprocess.SubprocessError) as error:
                filt = {"status": "error", "passes": False, "by_species": [], "reasons": [str(error)]}
        attempt = {"rank": rank, "window_id": row["window_id"], "start_1based": start,
                   "end_1based_inclusive": end, "length_nt": row["length_nt"],
                   "window_class": row.get("window_class"), "efficacy_priority_score": row["efficacy_priority_score"],
                   "length_penalty": row.get("length_penalty"), "adjusted_priority_score": row.get("adjusted_priority_score"),
                   "evaluation_directory": str(directory.relative_to(out)) if row["length_nt"] >= 21 else None,
                   "filter": filt, "outcome": "accepted" if filt["passes"] else "rejected_or_unverified"}
        attempts.append(attempt)
        last_evaluations, last_filter = evaluations, filt
        if not filt["passes"]:
            continue
        winner = dict(original) if original and row["window_id"] == original["window_id"] else {
            **row, "sequence": sequence,
            "reference_source": candidate.original["reference"].get("source"),
            "reference_version": candidate.original["reference"].get("version"),
            "green_region": next(r for r in selection["green_regions"] if r["region_id"] == row["region_id"]),
            "offtarget_matches": _overlapping_offtargets(integrated.get("offtarget_evaluations", []), start, end, len(candidate.sequence)),
            "selection_interpretation": "Mayor prioridad heuristica entre las ventanas que acreditaron el filtro computacional off-target; no eficacia validada."}
        winner.update(offtarget_safety_filter=filt, sequence_release_status="released_after_filter",
                      selected_priority_rank=rank)
        selection["selected_window"] = winner
        break
    winner = selection["selected_window"]
    audit = {"policy": POLICY, "eligible_window_count": eligible_count, "attempted_window_count": len(attempts),
             "selected_priority_rank": winner.get("selected_priority_rank") if winner else None,
             "status": "passed" if winner else "no_window_passed" if eligible_count else "no_eligible_windows",
             "attempts": attempts,
             "interpretation": "Se prueban en orden de prioridad; fallos y evaluaciones indeterminadas no autorizan una ganadora. No se modifican los puntajes."}
    selection["offtarget_selection"] = audit
    selection["selection_method"]["offtarget_selection_policy"] = POLICY
    selection["selection_method"]["final_selection_rule"] = "Primera ventana del ranking que cumple el filtro off-target en las tres especies."
    if winner:
        selection["selection_status"] = "selected_after_offtarget_filter"
    elif eligible_count:
        selection["selection_status"] = "no_window_passed_offtarget_filter"
        selection["warnings"].append("Sin ganadora: ninguna ventana evaluada acredito el filtro off-target en las tres especies.")
    return last_evaluations, last_filter, audit
