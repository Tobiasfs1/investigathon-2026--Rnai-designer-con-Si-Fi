"""Opt-in integration with the user's existing legacy si-Fi Conda runtime."""
import json
from pathlib import Path
import random
import subprocess

import pytest

from rnai_designer.config import Settings
from rnai_designer.legacy_runtime import detect_sifi_runtime
from rnai_designer.__main__ import main
from rnai_designer.offtarget_reports import load_evaluations
import rnai_designer.pipeline_workflow as pipeline


ROOT = Path(__file__).resolve().parents[2]
SIFI_ROOT = ROOT / "siFi21--1.2.3"
RUNTIME = detect_sifi_runtime(SIFI_ROOT, timeout=60)
pytestmark = [pytest.mark.integration,
              pytest.mark.skipif(not RUNTIME.get("available"), reason="si-Fi Conda runtime unavailable")]


def test_workflow_runs_original_sifi_headless_and_exports_mask(tmp_path, monkeypatch, capsys):
    settings = Settings.load()
    if not settings.bowtie or not settings.rnaplfold:
        pytest.skip("Bowtie/RNAplfold unavailable; see doctor")
    builder = settings.bowtie.with_name("bowtie-build.exe")
    if not builder.is_file():
        pytest.skip("bowtie-build.exe unavailable")
    rng = random.Random(137)
    sequence = "".join(rng.choice("ACGT") for _ in range(180))
    reference = tmp_path / "synthetic.fa"
    reference.write_text(">synthetic_target\n" + sequence + "\n>synthetic_target_2\n" + sequence + "\n",
                         encoding="ascii")
    database_dir = tmp_path / "db"
    database_dir.mkdir()
    built = subprocess.run([str(builder), str(reference), str(database_dir / "Dalbulus_maidis")],
                           capture_output=True, text=True, timeout=120)
    assert built.returncode == 0, built.stderr

    monkeypatch.setattr(pipeline, "run_offtargets",
                        lambda candidate, config, out, timeout: load_evaluations(None, candidate))
    out = tmp_path / "integrated-workflow"
    code = main(["workflow", "--gene-id", "DM_synthetic", "--sequence", sequence,
                 "--reference-source", "pytest synthetic fixture", "--out", str(out),
                 "--db-dir", str(database_dir), "--run-sifi-headless", "--sifi-timeout", "600"])
    assert code == 0, capsys.readouterr()

    ref_out = out / "references" / "01_synthetic_target"
    second_ref_out = out / "references" / "02_synthetic_target_2"
    result = json.loads((ref_out / "inputs" / "designer_original.json").read_text(encoding="utf-8"))
    mask = json.loads((ref_out / "inputs" / "green_mask.json").read_text(encoding="utf-8"))
    manifest = json.loads((ref_out / "designer_handoff" / "run.json").read_text(encoding="utf-8"))
    second_result = json.loads((second_ref_out / "inputs" / "designer_original.json").read_text(encoding="utf-8"))
    second_mask = json.loads((second_ref_out / "inputs" / "green_mask.json").read_text(encoding="utf-8"))
    assert len(result) == len(second_result) == 320
    assert {"is_efficient", "strand_selection", "end_stability", "target_site_accessibility"} <= set(result[0])
    assert mask["source_function"] == "general_helpers.get_target_data"
    assert mask["main_targets_selected"] == ["synthetic_target"]
    assert all(1 <= position <= 180 for position in mask["positions_1based"])
    assert second_mask["main_targets_selected"] == ["synthetic_target_2"]
    assert all(1 <= position <= 180 for position in second_mask["positions_1based"])
    assert mask["main_targets_selected"] != second_mask["main_targets_selected"]
    assert mask["positions_1based"] != second_mask["positions_1based"]
    assert mask["main_offtarget_overlap_start_count"] == 160
    assert second_mask["main_offtarget_overlap_start_count"] == 160
    assert mask["warnings"] and second_mask["warnings"]
    assert manifest["headless"] is True
    assert manifest["runtime"]["python"].startswith("2.7")
    assert (ref_out / "pipeline_result.json").is_file()
    assert (out / "reference_comparison.json").is_file()
    assert (out / "all_references_windows.csv").is_file()
    consolidated = json.loads((out / "reference_comparison.json").read_text(encoding="utf-8"))
    assert consolidated["completed_reference_count"] == 2
