"""Reconstrucción de bandas verdes, resumen de ventanas y asociación off-target."""
import copy
import csv
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from rnai_designer.data_models import read_json
from rnai_designer.window_explorer import _overlapping_offtargets, explore_windows, green_positions
from rnai_designer.window_reporting import write_window_outputs
from rnai_designer.__main__ import main

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
DEMO_REPORT = Path(__file__).resolve().parents[1] / "output" / "demo" / "report.json"


def row(name, start, off_target, efficient=True, accessibility=.2, hit="main"):
    return {"query_name": "DM_FICTICIO_001", "sirna_name": name, "sirna_position": start,
            "sirna_sequence": "ACGT" * 5 + "A", "is_efficient": efficient,
            "strand_selection": True, "end_stability": False,
            "target_site_accessibility": True, "accessibility_value": accessibility,
            "is_off_target": off_target, "hit_name": hit}


@pytest.fixture
def report():
    return read_json(DEMO_REPORT)


def test_green_definition_is_main_target_hit_coverage_and_order_sensitive():
    rows = [row("sirna1", 1, False, hit="gene"),
            row("sirna10", 10, True, hit="other"),
            row("sirna20", 20, False, hit="gene")]
    positions, warnings = green_positions(rows, 21)
    assert positions == set(range(1, 10)) | set(range(31, 41))
    assert warnings == []
    # si-Fi guarda snapshots por hit: un off-target posterior no reescribe
    # posiciones verdes que la fila main-target previa ya había registrado.
    positions_late_off = green_positions([rows[0], rows[2], rows[1]], 21)[0]
    assert positions_late_off == set(range(1, 41))


def test_green_is_not_efficiency_and_null_follows_original_python_truthiness():
    rows = [row("sirna1", 1, False, efficient=False), row("sirna2", 2, None)]
    positions, warnings = green_positions(rows, 21)
    # get_target_data uses `if data['is_off_target']`; null is falsey and enters main-target.
    assert positions == set(range(1, 23))
    assert len(warnings) == 1 and "null" in warnings[0]


def test_demo_windows_statistics_and_manual_detail(report):
    analysis = explore_windows(report, window_length=21, step=2, selected_start=1)
    assert analysis["coordinate_system"] == "1-based inclusive"
    assert analysis["green_regions"] == [{"region_id": "green-001", "start_1based": 1,
                                           "end_1based_inclusive": 25, "length_nt": 25}]
    assert len(analysis["windows"]) == 3
    selected = analysis["selected_window"]
    assert selected["sequence"] == report["candidate"]["reference"]["sequence"][:21]
    assert selected["sirna_start_count"] == 1
    assert selected["is_efficient_true_count"] == 1
    assert selected["is_efficient_false_count"] == 0
    assert selected["is_efficient_unknown_count"] == 0
    assert selected["accessibility_value_mean"] == pytest.approx(.2)
    assert selected["accessibility_value_min"] == pytest.approx(.2)
    assert selected["accessibility_value_max"] == pytest.approx(.2)
    assert selected["is_efficient_true_fraction_among_classified"] == 1


def test_explicit_mask_uses_original_helper_positions_and_preserves_separated_runs(report):
    report = copy.deepcopy(report)
    report["sifi"]["green_mask"] = {
        "schema_version": "1.0", "source_function": "general_helpers.get_target_data",
        "coordinate_system": "1-based positions", "sirna_size": 21,
        "positions_1based": [1, 2, 3, 10, 11],
    }
    analysis = explore_windows(report, 1)
    assert analysis["green_mask_source"] == "general_helpers.get_target_data"
    assert [(region["start_1based"], region["end_1based_inclusive"])
            for region in analysis["green_regions"]] == [(1, 3), (10, 11)]
    assert len(analysis["windows"]) == 5
    assert any("difiere de la reconstrucción" in warning for warning in analysis["warnings"])


def test_offtarget_alignments_link_only_when_their_query_interval_overlaps(report):
    report = copy.deepcopy(report)
    evaluation = report["offtarget_evaluations"][0]
    evaluation["criterion"] = {"status": "exact_match_detected", "alignments": [
        {"query_start_1based": 10, "query_end_1based_inclusive": 30,
         "reference_id": "NM_DEMO", "reference_start_0based": 8,
         "reference_end_0based_exclusive": 29, "strand": "+", "length": 21,
         "mismatch_count": 0, "match_type": "exact_21nt", "contains_exact_21nt": True},
        {"query_start_1based": 23, "query_end_1based_inclusive": None, "length": 21,
         "match_type": "unverified", "contains_exact_21nt": None},
        {"query_start_1based": None, "query_end_1based_inclusive": None,
         "match_type": "unverified", "contains_exact_21nt": None}]}
    analysis = explore_windows(report, 21, selected_start=1)
    assert analysis["windows"][0]["offtarget_alignment_count"] == 1
    overlapping = next(w for w in analysis["windows"] if w["start_1based"] == 5)
    assert overlapping["offtarget_alignment_count"] == 2
    match = analysis["selected_window"]["offtarget_matches"]
    assert len(match) == 1  # La ventana 1-21 se superpone al hit 10-30.
    assert any("carecen de coordenadas query válidas" in warning for warning in analysis["warnings"])
    linked = [w for w in analysis["windows"] if w["start_1based"] == 5][0]
    assert linked["human_exact_21nt_overlapping_rows"] == 1
    assert linked["human_unverified_overlapping_rows"] == 1
    unverified = _overlapping_offtargets(report["offtarget_evaluations"], 5, 25, 80)[1]
    assert unverified["query_interval_derived_from_start_and_length"] is True
    assert unverified["query_end_1based_inclusive_derived"] == 43


def test_conflicting_values_are_excluded_from_statistics(report):
    report = copy.deepcopy(report)
    original = report["sifi"]["records"][0]
    conflict = dict(original, hit_name="second_hit", accessibility_value=.8)
    report["sifi"]["records"].append(conflict)
    analysis = explore_windows(report, 21, selected_start=1)
    sirna = next(item for item in analysis["selected_window"]["sirnas"]
                 if item["sirna_name"] == original["sirna_name"])
    assert sirna["accessibility_value"] is None
    assert sirna["accessibility_value_conflict"] is True
    assert any("discrepa" in warning for warning in analysis["warnings"])


def test_missing_zone_data_errors_and_empty_zone_warns(report):
    no_data = copy.deepcopy(report)
    no_data["sifi"] = {"status": "not_evaluated", "records": []}
    with pytest.raises(ValueError, match="no tiene resultados"):
        explore_windows(no_data, 21)
    no_main = copy.deepcopy(report)
    for entry in no_main["sifi"]["records"]:
        entry["is_off_target"] = True
    analysis = explore_windows(no_main, 21)
    assert analysis["green_regions"] == [] and analysis["windows"] == []
    assert any("No se reconstruyeron" in warning for warning in analysis["warnings"])


@pytest.mark.parametrize("length,step", [(0, 1), (2, 0), (True, 1)])
def test_invalid_window_parameters(report, length, step):
    with pytest.raises(ValueError):
        explore_windows(report, length, step)


def test_out_of_reference_position_is_rejected(report):
    report["sifi"]["records"][0]["sirna_position"] = 1000
    with pytest.raises(ValueError, match="fuera de la referencia"):
        explore_windows(report, 21)


def test_window_cli_exports_comparison_table_plot_and_detail(report, tmp_path):
    report_path = tmp_path / "report.json"
    report_path.write_text(json.dumps(report), encoding="utf-8")
    output = tmp_path / "windows"
    assert main(["windows", "--report", str(report_path), "--window-length", "21",
                 "--step", "2", "--select-start", "1", "--out", str(output)]) == 0
    result = read_json(output / "window_analysis.json")
    assert result["selected_window"]["window_id"] == "green-001:1-21"
    assert (output / "selected_window.json").is_file()
    ET.parse(output / "windows.svg")
    assert (output / "windows.png").read_bytes().startswith(b"\x89PNG")
    with (output / "windows.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 3 and rows[0]["start_1based"] == "1"
    with (output / "window_offtargets.csv").open(encoding="utf-8", newline="") as handle:
        assert list(csv.DictReader(handle)) == []


def test_window_output_refuses_to_overwrite_existing_files(report, tmp_path):
    analysis = explore_windows(report, 21)
    output = tmp_path / "windows"
    output.mkdir()
    sentinel = output / "keep.txt"
    sentinel.write_text("keep", encoding="utf-8")
    with pytest.raises(ValueError, match="ya contiene archivos"):
        write_window_outputs(report, analysis, DEMO_REPORT, output)
    assert sentinel.read_text(encoding="utf-8") == "keep"
