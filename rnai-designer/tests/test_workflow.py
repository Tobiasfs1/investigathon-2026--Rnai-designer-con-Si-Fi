"""The workflow command prepares input and integrates the available analysis stages."""
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace
import pytest

from rnai_designer.data_models import Candidate, read_json
from rnai_designer.offtarget_reports import load_evaluations
from rnai_designer.__main__ import main
import rnai_designer.__main__ as cli
import rnai_designer.pipeline_workflow as pipeline


def test_candidate_builder_accepts_fasta_and_unknown_provenance():
    candidate = Candidate.from_fields("DM_test", ">DM_test\nacgt\nacgt", None, None)
    assert candidate.sequence == "ACGTACGT"
    assert candidate.original["reference"]["source"] is None
    assert candidate.original["tribolium_homolog_id"] is None
    assert candidate.original["priority_score"] is None


def test_candidate_builder_rejects_multiple_fasta_records():
    try:
        Candidate.from_fields("DM_test", ">gene1\nACGT\n>gene2\nTGCA")
    except ValueError as error:
        assert "una sola" in str(error)
    else:
        raise AssertionError("Expected multiple FASTA records to be rejected")


def test_python2_compatible_mask_bridge_exports_disconnected_runs(tmp_path):
    helper_root = tmp_path / "sifi"
    helper_root.mkdir()
    (helper_root / "general_helpers.py").write_text(
        "def get_target_data(path, size):\n"
        "    assert size == 21\n"
        "    return {}, {'target': set([1, 2, 5, 6])}, {}, []\n",
        encoding="utf-8")
    sifi_json = tmp_path / "sifi.json"
    sifi_json.write_text("[]", encoding="utf-8")
    mask = tmp_path / "mask.json"
    bridge = Path(__file__).resolve().parents[1] / "rnai_designer" / "legacy_bridge" / "export_green_mask.py"
    process = subprocess.run([sys.executable, str(bridge), "--sifi-json", str(sifi_json),
                              "--sirna-size", "21", "--sifi-root", str(helper_root), "--out", str(mask)],
                             capture_output=True, text=True)
    assert process.returncode == 0, process.stderr
    result = read_json(mask)
    assert result["source_function"] == "general_helpers.get_target_data"
    assert result["positions_1based"] == [1, 2, 5, 6]
    assert [(region["start_1based"], region["end_1based_inclusive"])
            for region in result["green_regions"]] == [(1, 2), (5, 6)]


def test_workflow_prepare_only_creates_query_and_resume_instructions(tmp_path, capsys):
    out = tmp_path / "prepared"
    assert main(["workflow", "--gene-id", "DM_test", "--sequence", "ACGT" * 30,
                 "--out", str(out)]) == 0
    candidate = read_json(out / "candidate.json")
    handoff = read_json(out / "designer_handoff" / "run.json")
    fasta = (out / "designer_handoff" / "query.fasta").read_text(encoding="utf-8")
    assert candidate["gene_id"] == "DM_test"
    assert handoff["query_id"] in fasta
    assert "--resume" in capsys.readouterr().out
    assert not (out / "pipeline_result.json").exists()


@pytest.mark.parametrize("window_filter_passes", [False, True])
def test_workflow_integrates_sifi_offtargets_windows_and_filter_status(tmp_path, monkeypatch, window_filter_passes):
    sequence = ("ACGT" * 45)[:180]
    sifi_path = tmp_path / "original_sifi.json"
    records = []
    for start in range(1, 161):
        records.append({
            "query_name": "DM_test", "sirna_name": f"sirna{start}", "sirna_position": start,
            "sirna_sequence": sequence[start - 1:start + 20], "is_efficient": True,
            "strand_selection": True, "end_stability": True, "target_site_accessibility": True,
            "accessibility_value": .2, "is_off_target": False, "hit_name": "DM_test",
            "reference_strand_pos": start - 1, "strand": "+", "mismatches": 0,
        })
    sifi_path.write_text(json.dumps(records), encoding="utf-8")
    run_path = tmp_path / "run.json"
    run_path.write_text(json.dumps({
        "query_id": "DM_test", "tool": "si-Fi", "tool_version": "v1.2.3-0008",
        "source": "local si-Fi", "database": "Dalbulus_maidis", "database_version": "unknown",
        "parameters": {"sirna_size": 21},
    }), encoding="utf-8")
    mask_path = tmp_path / "green_mask.json"
    mask_path.write_text(json.dumps({
        "schema_version": "1.0", "source_function": "general_helpers.get_target_data",
        "coordinate_system": "1-based positions", "sirna_size": 21,
        "positions_1based": list(range(1, 181)),
    }), encoding="utf-8")

    monkeypatch.setattr(cli.Settings, "load", lambda config, db_dir: SimpleNamespace(database_dir=tmp_path / "db"))
    def mocked_offtargets(candidate, settings, out, timeout):
        evaluations = load_evaluations(None, candidate)
        if window_filter_passes and out.name == "offtarget_selected_window":
            for item in evaluations:
                item["status"] = "imported_unvalidated"
                item["criterion"].update(status="criterion_met", label="Criterio computacional cumplido",
                                         exact_21nt_rows=0, reasons=[], alignments=[])
        return evaluations
    monkeypatch.setattr(pipeline, "run_offtargets", mocked_offtargets)
    out = tmp_path / "workflow"
    assert main(["workflow", "--gene-id", "DM_test", "--sequence", sequence,
                 "--reference-source", "test-fasta", "--reference-version", "v1",
                 "--sifi", str(sifi_path), "--run", str(run_path), "--green-mask", str(mask_path),
                 "--min-window", "180", "--max-window", "180", "--out", str(out)]) == 0

    result = read_json(out / "pipeline_result.json")
    assert result["evaluated_window_count"] == 1
    assert result["selected_window_offtarget_filter"]["passes"] is window_filter_passes
    if window_filter_passes:
        assert result["selected_window"]["start_1based"] == 1
        assert result["selected_window"]["end_1based_inclusive"] == 180
        assert result["selected_window"]["quality_score"] == 100
        assert result["selected_window"]["sequence"] == sequence
        assert result["selected_window"]["sequence_release_status"] == "released_after_filter"
    else:
        assert result["selected_window"] is None
        assert result["selection_status"] == "no_window_passed_offtarget_filter"
        assert result["offtarget_selection"]["attempted_window_count"] == 1
        assert result["selected_window_offtarget_filter"]["status"] == "not_evaluated"
    assert {entry["organism"] for entry in result["selected_window_offtarget_evaluations"]} == {
        "Homo sapiens", "Zea mays", "Apis mellifera"}
    assert result["green_mask_source"] == "general_helpers.get_target_data"
    assert result["offtarget_21nt_filter"]["overall"]["status"] == "not_evaluated"
    assert result["offtarget_21nt_filter"]["passes_all_species"] is False
    assert len(result["offtarget_21nt_filter"]["by_species"]) == 3
    assert (out / "integrated" / "report.json").is_file()
    assert (out / "window_selection" / "all_windows.csv").is_file()
    assert (out / "inputs" / "designer_original.json").is_file()
