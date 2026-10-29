"""Exporta tablas y comparativas descriptivas del explorador de ventanas."""
import csv
import json
import math
from pathlib import Path

from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from .sifi_adapter import provenance
from .window_explorer import _overlapping_offtargets
from .atomic_io import staged_directory

WINDOW_FIELDS = (
    "window_id", "region_id", "start_1based", "end_1based_inclusive", "length_nt",
    "sirna_start_count", "is_efficient_true_count", "is_efficient_false_count",
    "is_efficient_unknown_count", "is_efficient_classified_count",
    "is_efficient_true_fraction_among_classified", "strand_selection_true_count",
    "strand_selection_false_count", "strand_selection_unknown_count", "end_stability_true_count",
    "end_stability_false_count", "end_stability_unknown_count",
    "target_site_accessibility_true_count", "target_site_accessibility_false_count",
    "target_site_accessibility_unknown_count", "accessibility_value_count",
    "accessibility_value_missing_count", "accessibility_value_mean", "accessibility_value_median",
    "accessibility_value_min", "accessibility_value_max", "offtarget_alignment_count",
    "human_evaluation_status", "human_criterion_status", "human_exact_21nt_overlapping_rows",
    "human_mismatch_overlapping_rows", "human_unverified_overlapping_rows",
    "maize_evaluation_status", "maize_criterion_status", "maize_exact_21nt_overlapping_rows",
    "maize_mismatch_overlapping_rows", "maize_unverified_overlapping_rows",
    "bee_evaluation_status", "bee_criterion_status", "bee_exact_21nt_overlapping_rows",
    "bee_mismatch_overlapping_rows", "bee_unverified_overlapping_rows",
)


@staged_directory
def write_window_outputs(report: dict, analysis: dict, report_path: Path, out: Path) -> dict:
    if out.exists() and any(out.iterdir()):
        raise ValueError("La carpeta de exploración ya contiene archivos; elegir una salida nueva")
    out.mkdir(parents=True, exist_ok=True)
    analysis = {**analysis, "input_report": provenance(report_path)}
    if analysis["selected_window"] is not None:
        selected = {**analysis["selected_window"], "input_report": analysis["input_report"],
                    "candidate_id": analysis["candidate_id"], "priority_score": analysis["priority_score"],
                    "sifi_source": analysis["sifi_source"], "sifi_run": analysis["sifi_run"]}
        analysis["selected_window"] = selected
    with (out / 'window_analysis.json').open('w', encoding='utf-8') as handle:
        json.dump(analysis, handle, ensure_ascii=False, indent=2, allow_nan=False)
    with (out / "windows.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=WINDOW_FIELDS)
        writer.writeheader()
        writer.writerows({field: window.get(field) for field in WINDOW_FIELDS}
                         for window in analysis["windows"])
    fields = ["window_id", "candidate_id", "priority_score", "organism", "evaluation_status",
              "criterion_status", "database", "database_fingerprint_json", "result_sha256",
              "window_start_1based", "window_end_1based_inclusive",
              "window_overlap_start_1based", "window_overlap_end_1based_inclusive", "window_overlap_nt",
              "reference_id", "query_start_1based", "query_end_1based_inclusive",
              "reference_start_0based", "reference_end_0based_exclusive", "strand", "length",
              "mismatch_count", "match_type", "contains_exact_21nt"]
    with (out / "window_offtargets.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        candidate_data = report["candidate"]
        for window in analysis["windows"]:
            for match in _overlapping_offtargets(report.get("offtarget_evaluations", []),
                                                  window["start_1based"],
                                                  window["end_1based_inclusive"],
                                                  len(candidate_data["reference"]["sequence"])):
                writer.writerow({"window_id": window["window_id"],
                                 "candidate_id": analysis["candidate_id"],
                                 "priority_score": analysis["priority_score"],
                                 "window_start_1based": window["start_1based"],
                                 "window_end_1based_inclusive": window["end_1based_inclusive"],
                                 "database_fingerprint_json": json.dumps(match.get("database_fingerprint"), ensure_ascii=False),
                                 **{key: match.get(key) for key in fields if key in match}})
    if analysis["selected_window"] is not None:
        (out / "selected_window.json").write_text(
            json.dumps(analysis["selected_window"], ensure_ascii=False, indent=2, allow_nan=False),
            encoding="utf-8")
    _plot_windows(out / "windows.svg", analysis)
    return analysis


def _plot_windows(path: Path, analysis: dict) -> None:
    windows = analysis["windows"]
    figure = Figure(figsize=(12, 7), layout="constrained")
    FigureCanvasAgg(figure)
    access, efficiency = figure.subplots(2, 1, sharex=True)
    if windows:
        x = [w["start_1based"] for w in windows]
        known = [w for w in windows if w["accessibility_value_mean"] is not None]
        access.errorbar([w["start_1based"] for w in known],
                        [w["accessibility_value_mean"] for w in known],
                        yerr=[[w["accessibility_value_mean"] - w["accessibility_value_min"] for w in known],
                              [w["accessibility_value_max"] - w["accessibility_value_mean"] for w in known]],
                        fmt="o", markersize=3, capsize=2, color="#355f9c",
                        label="Media de accesibilidad; barras: mínimo–máximo")
        access.set_ylabel("accessibility_value")
        access.set_title("Resumen descriptivo de valores importados por ventana")
        access.legend(loc="best", fontsize=8)
        width = max(0.25, min(4, analysis["step"] * 0.75))
        true = [w["is_efficient_true_count"] for w in windows]
        false = [w["is_efficient_false_count"] for w in windows]
        unknown = [w["is_efficient_unknown_count"] for w in windows]
        efficiency.bar(x, true, width=width, color="#16853b", label="si-Fi true")
        efficiency.bar(x, false, width=width, bottom=true, color="#c83232", label="si-Fi false")
        bottoms = [a + b for a, b in zip(true, false)]
        efficiency.bar(x, unknown, width=width, bottom=bottoms, color="#888888", label="null/conflicto")
        efficiency.set_ylabel("Inicios de siRNA")
        efficiency.set_xlabel("Inicio de ventana (1-based)")
        ratio_ax = efficiency.twinx()
        ratio = [w["is_efficient_true_fraction_among_classified"] for w in windows]
        ratio_x = [start for start, value in zip(x, ratio) if value is not None]
        ratio_y = [value for value in ratio if value is not None]
        ratio_ax.scatter(ratio_x, ratio_y, color="#222222", marker=".",
                         label="Fracción etiquetada true entre etiquetas conocidas")
        ratio_ax.set_ylim(-0.05, 1.05)
        ratio_ax.set_ylabel("Fracción descriptiva")
        handles, labels = efficiency.get_legend_handles_labels()
        handles2, labels2 = ratio_ax.get_legend_handles_labels()
        efficiency.legend(handles + handles2, labels + labels2, loc="upper left", fontsize=8)
        selected = analysis["selected_window"]
        if selected:
            for axis in (access, efficiency):
                axis.axvspan(selected["start_1based"] - 0.5,
                             selected["end_1based_inclusive"] + 0.5,
                             color="#f2c94c", alpha=0.2)
        access.set_xlim(min(x) - 0.5, max(x) + analysis["window_length"] - 0.5)
        if len(x) > 25:
            stride = math.ceil(len(x) / 20)
            efficiency.set_xticks(x[::stride])
            efficiency.tick_params(axis="x", labelrotation=45)
    else:
        access.text(0.5, 0.5, "No hay ventanas comparables en las regiones verdes.",
                    ha="center", va="center", transform=access.transAxes)
        access.set_axis_off()
        efficiency.set_axis_off()
    figure.suptitle(f'{analysis["candidate_id"]} — ventanas dentro de bandas verdes de si-Fi\n'
                    "Estadísticas descriptivas; no son predicciones de eficacia ni seguridad")
    figure.savefig(path, format="svg", metadata={"Date": None})
    figure.savefig(path.with_suffix(".png"), dpi=160)
    figure.clear()
