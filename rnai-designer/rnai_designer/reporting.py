"""Exportaciones completas y trazables."""
import csv
import json
from . import __version__
from .visualization import write_plot, write_offtarget_plot
from .offtarget_assessment import overall
from .sifi_adapter import FIELDS
from .data_models import Candidate, ImportedResult
from .atomic_io import staged_directory


@staged_directory(allow_existing=True)
def write_reports(out, candidate: Candidate, imported: ImportedResult, summary: dict,
                  evaluations: list[dict]) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    warnings = imported["warnings"] + summary["warnings"]
    for entry in evaluations:
        warnings.extend(entry.get("result", {}).get("warnings", []))
        if entry["status"] == "not_evaluated":
            warnings.append(f'{entry["organism"]}: no evaluado')
        if entry["status"] == "error":
            warnings.append(f'{entry["organism"]}: error: {entry.get("error")}')
        warnings.extend(f'{entry["organism"]}: {reason}' for reason in entry["criterion"]["reasons"])
    report = {"schema_version": "1.1", "producer": {"name": "rnai-designer", "version": __version__},
              "csv_schema_version": "2.0",
              "candidate": candidate.original, "sifi": imported,
              "efficiency_summary": summary, "offtarget_evaluations": evaluations,
              "offtarget_criterion": overall(evaluations),
              "warnings": warnings, "safety_conclusion": "undetermined"}
    with (out / 'report.json').open('w', encoding='utf-8') as handle:
        json.dump(report, handle, ensure_ascii=False, indent=2, allow_nan=False)
    with (out / "records.csv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["gene_id", "priority_score", "context", "organism", *FIELDS,
                  "record_json", "record_index", "dataset_ref"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        datasets = [("primary", "Dalbulus maidis", imported, 'report.json#/sifi')]
        datasets += [("offtarget", e["organism"], e["result"], f'report.json#/offtarget_evaluations/{i}/result')
                     for i, e in enumerate(evaluations) if "result" in e]
        for context, organism, dataset, dataset_ref in datasets:
            for index, row in enumerate(dataset["records"]):
                writer.writerow({**{key: row.get(key) for key in FIELDS},
                    "gene_id": candidate.gene_id, "priority_score": candidate.original["priority_score"],
                    "context": context, "organism": organism, "record_json": json.dumps(row, ensure_ascii=False),
                    "record_index": index, "dataset_ref": dataset_ref})
    with (out / "evaluations.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["gene_id", "organism", "status", "hit_rows", "safety_conclusion", "evaluation_ref",
                         "criterion_status", "criterion_label", "exact_21nt_rows", "mismatch_rows"])
        for index, entry in enumerate(evaluations):
            criterion = entry["criterion"]
            writer.writerow([candidate.gene_id, entry["organism"], entry["status"], entry.get("hit_rows"),
                             "undetermined", f'report.json#/offtarget_evaluations/{index}', criterion["status"],
                             criterion["label"], criterion["exact_21nt_rows"], criterion["mismatch_rows"]])
    with (out / "offtarget_matches.csv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["gene_id", "organism", "database", "criterion_status", "reference_id", "query_start_1based",
                  "query_end_1based_inclusive", "reference_start_0based", "reference_end_0based_exclusive",
                  "strand", "length", "mismatch_count", "match_type", "contains_exact_21nt"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for entry in evaluations:
            for alignment in entry["criterion"]["alignments"]:
                writer.writerow({"gene_id": candidate.gene_id, "organism": entry["organism"],
                                 "database": entry.get("manifest", {}).get("database"),
                                 "criterion_status": entry["criterion"]["status"], **alignment})
    with (out / "positions.csv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["sirna_name", "start_1based", "end_1based_inclusive", "accessibility_value",
                  "is_efficient", "hit_rows", "conflict"]
        writer = csv.DictWriter(handle, fieldnames=["gene_id", *fields])
        writer.writeheader()
        writer.writerows({"gene_id": candidate.gene_id, **point} for point in summary["points"])
    simulated = any("DATOS FICTICIOS" in warning for warning in warnings)
    write_plot(out / "efficiency.svg", summary["points"], candidate.gene_id,
               len(candidate.sequence), simulated=simulated, evaluated=imported.get("status") != "not_evaluated")
    write_offtarget_plot(out / "offtargets.svg", evaluations, candidate.gene_id,
                         len(candidate.sequence), simulated=simulated)
    return report
