"""Files and full-gene visualization for automatic window selection."""
import csv
import json
from pathlib import Path

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from .sifi_adapter import provenance


def _jsonable_row(row: dict) -> dict:
    return {key: json.dumps(value, ensure_ascii=False, sort_keys=True) if isinstance(value, (dict, list)) else value
            for key, value in row.items()}


def write_selection_outputs(report: dict, selection: dict, report_path: Path, out: Path) -> dict:
    if out.exists() and any(out.iterdir()):
        raise ValueError("La carpeta de selección ya contiene archivos; elegir una salida nueva")
    out.mkdir(parents=True, exist_ok=True)
    input_info = provenance(report_path)
    selection = {**selection, "input_report": input_info}
    if selection["selected_window"] is not None:
        selection["selected_window"] = {
            **selection["selected_window"],
            "candidate_id": selection["candidate_id"],
            "priority_score": report["candidate"].get("priority_score"),
            "sifi_source": report["sifi"].get("source"),
            "sifi_run": report["sifi"].get("run"),
            "input_report": input_info,
            "coordinate_system": selection["coordinate_system"],
        }
    (out / "all_windows.csv").write_text("", encoding="utf-8")
    if selection["windows"]:
        fields = list(dict.fromkeys(key for row in selection["windows"] for key in row))
        with (out / "all_windows.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(_jsonable_row(row) for row in selection["windows"])
    else:
        (out / "all_windows.csv").write_text("window_id,region_id,start_1based,end_1based_inclusive,length_nt,efficacy_priority_score,quality_score\n",
                                             encoding="utf-8")
    if selection["selected_window"] is not None:
        (out / "selected_window.json").write_text(
            json.dumps(selection["selected_window"], ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    report_payload = {key: value for key, value in selection.items() if key not in ("windows", "selected_window")}
    report_payload["selected_window_id"] = selection["selected_window"]["window_id"] if selection["selected_window"] else None
    report_payload["selected_quality_score"] = selection["selected_window"]["quality_score"] if selection["selected_window"] else None
    report_payload["selected_efficacy_priority_score"] = (selection["selected_window"]["efficacy_priority_score"]
                                                           if selection["selected_window"] else None)
    if selection["selected_window"]:
        best_score = max(window["efficacy_priority_score"] for window in selection["windows"]
                         if window["efficacy_priority_score"] is not None)
        tied_count = sum(window["efficacy_priority_score"] == best_score for window in selection["windows"])
        report_payload["selection_reasons"] = {
            "efficacy_priority_components": selection["selected_window"]["efficacy_priority_components"],
            "quality_score": selection["selected_window"]["quality_score"],
            "best_score_tie_count": tied_count,
            "tie_break_applied": tied_count > 1,
            "reason": "Mayor prioridad heuristica; calidad de evidencia se conserva como auditoria independiente.",
            "biological_interpretation": "Prioridad exploratoria basada en datos originales, GC y diversidad; eficacia real no validada.",
        }
    else:
        report_payload["selection_reasons"] = []
    (out / "selection_report.json").write_text(
        json.dumps(report_payload, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    _plot_selected(out / "selected_window.svg", report, selection)
    return selection


def _plot_selected(path: Path, report: dict, selection: dict) -> None:
    reference = report["candidate"]["reference"]["sequence"]
    selected = selection["selected_window"]
    figure = Figure(figsize=(14, 4.5), layout="constrained")
    FigureCanvasAgg(figure)
    gene, tracks = figure.subplots(2, 1, gridspec_kw={"height_ratios": [2, 1]}, sharex=True)
    for region in selection["green_regions"]:
        gene.axvspan(region["start_1based"] - .5, region["end_1based_inclusive"] + .5,
                     color="#53a653", alpha=.35)
    gene.plot([1, len(reference)], [0, 0], color="#555555", linewidth=3)
    if selected:
        gene.axvspan(selected["start_1based"] - .5, selected["end_1based_inclusive"] + .5,
                     color="#f2b632", alpha=.75, label="Ventana destacada")
        gene.annotate(f'{selected["start_1based"]}–{selected["end_1based_inclusive"]} nt',
                      ((selected["start_1based"] + selected["end_1based_inclusive"]) / 2, .08),
                      ha="center", va="bottom")
    gene.set_yticks([])
    gene.set_ylabel("Gen / regiones verdes")
    gene.set_title(f'{selection["candidate_id"]}: ventana priorizada por metrica exploratoria')
    evaluations = report.get("offtarget_evaluations", [])
    species = [("Homo sapiens", "human", "#7a5195"), ("Zea mays", "maize", "#ef5675"),
               ("Apis mellifera", "bee", "#ffa600")]
    offsets = {name: i for i, (name, _, _) in enumerate(species)}
    for evaluation in evaluations:
        name = evaluation.get("organism")
        if name not in offsets:
            continue
        color = next(c for n, _, c in species if n == name)
        xs, xe = [], []
        for alignment in evaluation.get("criterion", {}).get("alignments", []):
            start = alignment.get("query_start_1based")
            end = alignment.get("query_end_1based_inclusive")
            if end is None and type(start) is int and type(alignment.get("length")) is int:
                end = start + alignment["length"] - 1
            if type(start) is int and type(end) is int and 1 <= start <= end <= len(reference):
                xs.append(start)
                xe.append(end)
                tracks.plot([start, end], [offsets[name], offsets[name]], color=color, linewidth=3)
    tracks.set_yticks(list(offsets.values()), [name for name, _, _ in species])
    tracks.set_ylim(-.6, 2.6)
    tracks.set_xlabel("Coordenada en referencia (1-based, inclusiva)")
    tracks.set_ylabel("Off-target por especie")
    tracks.set_xlim(.5, len(reference) + .5)
    if len(reference) > 10000:
        tracks.ticklabel_format(style="plain", axis="x")
    figure.savefig(path, format="svg", metadata={"Date": None})
    figure.savefig(path.with_suffix(".png"), dpi=160)
    figure.clear()
