import csv
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

import rnai_designer.__main__ as cli
from rnai_designer.data_models import Candidate, read_json


@pytest.mark.parametrize("targets,failed", [
    (["ref_A"], set()),
    (["ref_A", "ref_B"], set()),
    (["ref_A", "ref_B"], {"ref_A"}),
])
def test_multi_reference_report_keeps_each_target_separate(tmp_path, monkeypatch, targets, failed):
    candidate = Candidate.from_fields("ATPase", "ACGT" * 100, "synthetic", "v1")
    args = SimpleNamespace(timeout=10, min_window=180, max_window=250, sirna_size=21,
                           _sifi_runtime={"name": "test-runtime"})
    calls = []

    def fake_sifi(target_args, candidate, reference_out, settings, target):
        if target in failed:
            raise ValueError("synthetic si-Fi failure")
        calls.append(target)
        source = reference_out / "inputs" / "designer_original.json"
        source.parent.mkdir(parents=True, exist_ok=True)
        source.write_text("[]", encoding="utf-8")
        return source, {"records": []}, {"main_target": target}

    def fake_pipeline(candidate, imported, settings, reference_out, source, run, timeout,
                      min_length, max_length, sirna_size):
        target = read_json(reference_out / "designer_handoff" / "run.json")["query_id"]
        start = 11 if "ref_A" in reference_out.name else 51
        row = {"region_id": "green_1", "start_1based": start, "end_1based_inclusive": start + 179,
               "length_nt": 180, "window_class": "conventional_180_250nt", "quality_score": 90,
               "efficacy_priority_score": 80, "sirna_record_count": 160,
               "gc_percent": 50, "sirna_unique_sequences_count": 100}
        summary = {"green_regions": [{"region_id": "green_1", "start_1based": start,
                                      "end_1based_inclusive": start + 179, "length_nt": 180}],
                   "evaluated_window_count": 1, "selected_window": {"start_1based": start,
                   "end_1based_inclusive": start + 179, "length_nt": 180,
                   "efficacy_priority_score": 80, "quality_score": 90},
                   "selected_window_offtarget_filter": {"passes": False}, "warnings": []}
        return {"summary": summary, "selection": {
            "windows": [row], "selection_method": {"score_name": "efficacy_priority_score",
            "formula": "test", "evidence_quality_score_formula": "test"}}}

    monkeypatch.setattr(cli, "_run_legacy_sifi", fake_sifi)
    monkeypatch.setattr(cli, "run_integrated_workflow", fake_pipeline)
    out = tmp_path / "workflow"
    out.mkdir()
    report = cli._run_multi_reference_workflow(args, candidate, object(), out, targets)

    assert report["reference_count"] == len(targets)
    assert report["completed_reference_count"] == len(set(targets) - failed)
    assert report["failed_reference_count"] == len(failed)
    assert report["status"] == ("partial" if failed and len(failed) < len(targets) else "complete")
    assert len(report["references"]) == len(targets)
    assert all(item["main_target"] == target for item, target in zip(report["references"], targets))
    assert (out / "reference_comparison.json").is_file()
    with (out / "all_references_windows.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert {row["main_target"] for row in rows} == set(targets) - failed
    if len(rows) > 1:
        assert rows[0]["start_1based"] != rows[1]["start_1based"]
    for entry in report["references"]:
        ref_dir = out / entry["output_directory"]
        if entry["status"] == "error":
            assert (ref_dir / "reference_error.json").is_file()
        else:
            assert (ref_dir / "designer_handoff" / "run.json").is_file()


def test_workflow_auto_detects_and_processes_every_target_from_one_input(tmp_path, monkeypatch):
    database_dir = tmp_path / "databases"
    database_dir.mkdir()
    monkeypatch.setattr(cli.Settings, "load", lambda config, db_dir:
                        SimpleNamespace(database_dir=database_dir))
    monkeypatch.setattr(cli, "_legacy_sifi_context", lambda args, candidate, out, settings:
                        setattr(args, "_sifi_runtime", {"name": "test-runtime"}))
    monkeypatch.setattr(cli, "_list_legacy_targets", lambda *args: {
        "available_targets": [{"hit_name": "GITV01001458_1"}, {"hit_name": "GITV01002133_1"}]})
    captured = {}
    def run_all(args, candidate, settings, out, targets):
        captured["targets"] = targets
        return {"reference_count": 2, "completed_reference_count": 2,
                "failed_reference_count": 0, "combined_window_count": 3}
    monkeypatch.setattr(cli, "_run_multi_reference_workflow", run_all)
    out = tmp_path / "ATPase"
    code = cli.main(["workflow", "--gene-id", "ATPase", "--sequence", "ACGT" * 50,
                     "--out", str(out), "--run-sifi-headless"])
    assert code == 0
    assert captured["targets"] == ["GITV01001458_1", "GITV01002133_1"]
