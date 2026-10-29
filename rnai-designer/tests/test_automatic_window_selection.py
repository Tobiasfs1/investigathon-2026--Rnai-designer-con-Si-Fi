"""Exhaustive automatic selection is a transparent data-quality review aid."""
import copy
import csv
import json
import pytest

from rnai_designer.automatic_window_selection import select_windows
from rnai_designer.automatic_window_reporting import write_selection_outputs
from rnai_designer.data_models import read_json
from rnai_designer.__main__ import main


DEMO = read_json(__import__("pathlib").Path(__file__).resolve().parents[1] / "output" / "demo" / "report.json")


def make_report(regions, sequence_length=600, missing=False):
    report = copy.deepcopy(DEMO)
    sequence = ("ACGT" * ((sequence_length + 3) // 4))[:sequence_length]
    report["candidate"]["reference"]["sequence"] = sequence
    report["candidate"]["reference"]["source"] = "synthetic-test"
    report["candidate"]["reference"]["version"] = "v1"
    report["sifi"]["run"]["query_id"] = "DM_FICTICIO_001"
    report["sifi"]["run"]["parameters"]["sirna_size"] = 21
    rows = []
    for region_start, region_end in regions:
        for start in range(region_start, region_end - 20 + 1):
            rows.append({
                "query_name": "DM_FICTICIO_001", "sirna_name": f"sirna-{start}",
                "sirna_position": start, "sirna_sequence": sequence[start - 1:start + 20],
                "is_off_target": False, "hit_name": f"hit-{region_start}",
                "is_efficient": None if missing else True,
                "strand_selection": None if missing else True,
                "end_stability": None if missing else False,
                "target_site_accessibility": None if missing else True,
                "accessibility_value": None if missing else .25,
            })
    report["sifi"]["records"] = rows
    report["sifi"]["status"] = "evaluated"
    report["sifi"]["source"] = {"sha256": "a" * 64, "path": "synthetic.json"}
    report["offtarget_evaluations"] = []
    return report


def test_exhaustive_range_includes_180_and_250_and_all_integer_lengths():
    result = select_windows(make_report([(1, 250)]))
    assert result["evaluated_window_count"] == sum(251 - length for length in range(180, 251))
    assert {w["length_nt"] for w in result["windows"]} == set(range(180, 251))
    assert any(w["start_1based"] == 1 and w["length_nt"] == 180 for w in result["windows"])
    assert any(w["start_1based"] == 1 and w["length_nt"] == 250 for w in result["windows"])
    assert all(w["end_1based_inclusive"] - w["start_1based"] + 1 == w["length_nt"] for w in result["windows"])


def test_one_nt_stride_and_no_windows_cross_green_region_gaps():
    shifted = select_windows(make_report([(1, 181)]), 180, 180)
    assert [(w["start_1based"], w["end_1based_inclusive"]) for w in shifted["windows"]] == [(1, 180), (2, 181)]
    separated = select_windows(make_report([(1, 200), (221, 500)]), 180, 250)
    assert separated["evaluated_window_count"] > 0
    for window in separated["windows"]:
        region = next(r for r in separated["green_regions"] if r["region_id"] == window["region_id"])
        assert region["start_1based"] <= window["start_1based"]
        assert window["end_1based_inclusive"] <= region["end_1based_inclusive"]
        assert not (window["start_1based"] <= 220 and window["end_1based_inclusive"] >= 201)


def test_empty_range_and_length_validation():
    empty_mask = make_report([(1, 21)])
    empty_mask["sifi"]["records"][0]["is_off_target"] = True
    result = select_windows(empty_mask)
    assert result["selection_status"] == "no_windows_available"
    assert result["selected_window"] is None and result["evaluated_window_count"] == 0
    try:
        select_windows(make_report([(1, 250)]), 250, 180)
    except ValueError as error:
        assert "min_length" in str(error)
    else:
        raise AssertionError("Expected invalid range to fail")


def test_short_green_region_is_one_exact_window_and_does_not_extend():
    result = select_windows(make_report([(1, 179)]), 180, 250)
    assert len(result["windows"]) == 1
    window = result["windows"][0]
    assert (window["start_1based"], window["end_1based_inclusive"], window["length_nt"]) == (1, 179, 179)
    assert window["window_class"] == "short_green_region"
    assert window["conventional_range_comparable"] is False


def test_short_regions_keep_gaps_and_report_insufficient_sirna_evidence():
    report = make_report([(1, 179), (201, 379)], sequence_length=400)
    result = select_windows(report, 180, 250)
    short = [row for row in result["windows"] if row["window_class"] == "short_green_region"]
    assert [(row["start_1based"], row["end_1based_inclusive"]) for row in short] == [(1, 179), (201, 379)]
    assert all(row["end_1based_inclusive"] - row["start_1based"] + 1 == 179 for row in short)
    assert all(row["conventional_range_comparable"] is False for row in short)


def test_short_region_from_explicit_mask_is_kept_when_siRNA_records_are_missing():
    report = make_report([(1, 21)])
    report["sifi"]["records"] = []
    report["sifi"]["green_mask"] = {
        "source_function": "general_helpers.get_target_data",
        "coordinate_system": "1-based positions", "sirna_size": 21,
        "positions_1based": list(range(1, 26)),
    }
    result = select_windows(report, 180, 250)
    assert result["evaluated_window_count"] == 1
    window = result["windows"][0]
    assert (window["start_1based"], window["end_1based_inclusive"], window["length_nt"]) == (1, 25, 25)
    assert window["data_sufficiency_status"] == "insufficient_no_complete_sirna_records"
    assert result["selected_window"] is None
    assert any("evidencia insuficiente" in warning for warning in result["warnings"])


def test_tie_break_is_earliest_start_and_missing_values_are_not_imputed():
    result = select_windows(make_report([(1, 181)], missing=True), 180, 180)
    assert result["selected_window"] is None
    assert result["selection_status"] == "no_evaluable_efficacy_metrics"
    assert result["windows"][0]["sirna_complete_count"] == 0
    assert result["windows"][0]["missing_data_percent"] == 100
    assert result["windows"][0]["is_efficient_true_count"] == 0
    assert result["windows"][0]["is_efficient_evaluable_count"] == 0
    assert result["windows"][0]["accessibility_value_mean"] is None


def test_metrics_quality_and_duplicate_audit_are_separate():
    report = make_report([(1, 180)])
    first = report["sifi"]["records"][0]
    duplicate = dict(first)
    report["sifi"]["records"].append(duplicate)
    result = select_windows(report, 180, 180)
    selected = result["selected_window"]
    assert selected["length_nt"] == 180
    assert selected["gc_percent"] == 50
    assert selected["is_efficient_true_count"] == 160
    assert selected["is_efficient_true_fraction_among_evaluable"] == 1
    assert selected["sirna_record_count"] == selected["sirna_unique_count"] + 1
    assert selected["exact_duplicate_record_count"] == 1
    assert selected["quality_score"] == 100 * sum(selected["quality_components"].values()) / 5
    assert selected["efficacy_priority_score"] == pytest.approx(87.5)
    assert selected["efficacy_priority_components"]["sirna_gc_30_52_fraction"] == pytest.approx(.5)
    assert selected["efficacy_priority_components"]["unique_sirna_sequence_fraction"] == pytest.approx(.025)
    assert "is_efficient" in result["selection_method"]["formula"]


def test_efficacy_priority_uses_si_fi_efficiency_gc_band_and_reports_diversity():
    report = make_report([(1, 180), (181, 360)], sequence_length=360)
    report["candidate"]["reference"]["sequence"] = "ACGT" * 45 + "A" * 180
    for row in report["sifi"]["records"]:
        start = row["sirna_position"]
        row["sirna_sequence"] = report["candidate"]["reference"]["sequence"][start - 1:start + 20]
    result = select_windows(report, 180, 180)
    by_start = {row["start_1based"]: row for row in result["windows"]}
    assert by_start[1]["efficacy_priority_components"]["sirna_gc_30_52_fraction"] > \
           by_start[181]["efficacy_priority_components"]["sirna_gc_30_52_fraction"]
    assert by_start[1]["efficacy_priority_components"]["unique_sirna_sequence_fraction"] > \
           by_start[181]["efficacy_priority_components"]["unique_sirna_sequence_fraction"]
    assert result["selected_window"]["efficacy_priority_score"] == max(
        row["efficacy_priority_score"] for row in result["windows"])


def test_conflicting_metric_reduces_evidence_and_offtargets_stay_by_species():
    report = make_report([(1, 180)])
    conflict = dict(report["sifi"]["records"][0], hit_name="second-hit", accessibility_value=.9)
    report["sifi"]["records"].append(conflict)
    report["offtarget_evaluations"] = [{
        "organism": "Homo sapiens", "status": "imported_unvalidated",
        "criterion": {"status": "exact_match_detected", "alignments": [{
            "query_start_1based": 20, "query_end_1based_inclusive": 40,
            "match_type": "exact_21nt", "contains_exact_21nt": True,
        }]},
    }]
    result = select_windows(report, 180, 180)["selected_window"]
    assert result["sirnas_with_metric_conflicts_count"] == 1
    assert result["accessibility_value_count"] == 159
    assert result["human_exact_21nt_overlapping_rows"] == 1
    assert result["maize_exact_21nt_overlapping_rows"] == 0
    assert result["bee_exact_21nt_overlapping_rows"] == 0
    assert result["offtarget_by_species"]["human_evaluation_status"] != "not_evaluated"


def test_cli_writes_auditable_outputs_and_full_gene_figure(tmp_path):
    report = make_report([(1, 250)])
    report_path = tmp_path / "report.json"
    report_path.write_text(json.dumps(report), encoding="utf-8")
    out = tmp_path / "automatic"
    assert main(["select-windows", "--report", str(report_path), "--out", str(out)]) == 0
    selected = read_json(out / "selected_window.json")
    selection_report = read_json(out / "selection_report.json")
    assert selected["window_id"] == selection_report["selected_window_id"]
    assert selection_report["evaluated_window_count"] == sum(251 - n for n in range(180, 251))
    assert (out / "selected_window.svg").is_file()
    assert (out / "selected_window.png").read_bytes().startswith(b"\x89PNG")
    with (out / "all_windows.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == selection_report["evaluated_window_count"]
    assert rows[0]["start_1based"] == "1"
