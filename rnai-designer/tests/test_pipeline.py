"""Pruebas de integridad de contratos, evidencia faltante y exportación."""
from pathlib import Path
import csv
import json
import xml.etree.ElementTree as ET

import pytest

from rnai_designer.data_models import Candidate, read_json
from rnai_designer.sifi_adapter import import_sifi
from rnai_designer.efficiency_analysis import summarize
from rnai_designer.offtarget_reports import load_evaluations
from rnai_designer.__main__ import main

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


@pytest.fixture
def candidate():
    return Candidate.load(EXAMPLES / "candidate.json")


@pytest.fixture
def rows():
    return read_json(EXAMPLES / "sifi_result.json")


def save(tmp_path, value, name="input.json"):
    path = tmp_path / name
    path.write_text(json.dumps(value), encoding="utf-8")
    return path


def test_complete_cli(tmp_path):
    out = tmp_path / "report"
    assert main(["--candidate", str(EXAMPLES / "candidate.json"),
                 "--sifi", str(EXAMPLES / "sifi_result.json"),
                 "--run", str(EXAMPLES / "run.json"),
                 "--offtargets", str(EXAMPLES / "offtargets.json"), "--out", str(out)]) == 0
    report = read_json(out / "report.json")
    assert report["candidate"] == read_json(EXAMPLES / "candidate.json")
    assert report["sifi"]["records"] == read_json(EXAMPLES / "sifi_result.json")
    assert len(report["offtarget_evaluations"]) == 3
    assert report["safety_conclusion"] == "undetermined"
    assert any("FICTICIOS" in w for w in report["warnings"])
    assert len(report["sifi"]["source"]["sha256"]) == 64
    assert report["sifi"]["run_source"]["file"] == "run.json"
    ET.parse(out / "efficiency.svg")
    assert "FICTICIOS" in (out / "efficiency.svg").read_text(encoding="utf-8")
    assert (out / "efficiency.png").read_bytes().startswith(b"\x89PNG")
    with (out / "records.csv").open(encoding="utf-8", newline="") as handle:
        records = list(csv.DictReader(handle))
    assert len(records) == 7
    assert json.loads(records[0]["record_json"])["mismatches"] == 0
    assert records[0]["sirna_position"] == "1"
    with (out / "positions.csv").open(encoding="utf-8", newline="") as handle:
        assert len(list(csv.DictReader(handle))) == 4


def test_no_double_count_and_conflicts(rows):
    result = summarize(rows + [dict(rows[0], hit_name="second")])
    assert result["hit_rows"] == 5 and result["unique_sirnas"] == 4
    result = summarize([rows[0], dict(rows[0], accessibility_value=.9)])
    assert result["points"][0]["accessibility_value"] is None
    assert result["points"][0]["is_efficient"] is None
    assert result["warnings"]


@pytest.mark.parametrize("change", [
    {"sirna_position": 0}, {"sirna_position": True}, {"sirna_position": 79},
    {"is_efficient": "false"}, {"accessibility_value": 2}, {"query_name": "OTHER"},
    {"sirna_sequence": "XYZ"}, {"reference_strand_pos": -1}, {"strand": "invalid"},
])
def test_reject_invalid(candidate, rows, tmp_path, change):
    with pytest.raises(ValueError):
        import_sifi(save(tmp_path, [dict(rows[0], **change)]), candidate)


def test_empty_is_unknown(candidate, tmp_path):
    assert import_sifi(save(tmp_path, []), candidate)["warnings"]
    assert all(e["status"] == "not_evaluated" for e in load_evaluations(None, candidate))


def test_raw_and_missing_preserved(candidate, rows, tmp_path):
    row = rows[0]
    row.pop("accessibility_value")
    row.update(mismatches="3:A>G", extra={"original": True})
    result = import_sifi(save(tmp_path, [row]), candidate)
    assert result["records"] == [row]
    assert any("accessibility_value" in w for w in result["warnings"])
    assert any("inicios posibles" in w for w in result["warnings"])


@pytest.mark.parametrize("text", ['{"wrong": 1}', '[', '[NaN]', '[1e999]',
                                  '[{"key":1,"key":2}]'])
def test_malformed_json(candidate, tmp_path, text):
    path = tmp_path / "bad.json"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError):
        import_sifi(path, candidate)


def test_missing_file_is_error(candidate, tmp_path):
    manifest = read_json(EXAMPLES / "offtargets.json")
    with pytest.raises(OSError):
        load_evaluations(save(tmp_path, manifest), candidate)


def test_duplicate_species_rejected(candidate, tmp_path):
    manifest = read_json(EXAMPLES / "offtargets.json")
    entry = manifest["evaluations"][0]
    entry["results_file"] = "hits.json"
    manifest["evaluations"] = [entry, entry]
    save(tmp_path, read_json(EXAMPLES / "human.json"), "hits.json")
    with pytest.raises(ValueError, match="repetido"):
        load_evaluations(save(tmp_path, manifest), candidate)


def test_empty_offtarget_preserves_evidence(candidate, tmp_path):
    entry = read_json(EXAMPLES / "offtargets.json")["evaluations"][0]
    entry["results_file"] = "empty.json"
    save(tmp_path, [], "empty.json")
    path = save(tmp_path, {"schema_version": "1.0", "evaluations": [entry]})
    results = load_evaluations(path, candidate)
    assert results[0]["status"] == "imported_unvalidated"
    assert results[0]["safety_conclusion"] == "undetermined"
    assert results[0]["manifest"] == entry
    assert results[0]["hit_rows"] == 0 and results[0]["result"]["warnings"]
    assert results[1]["status"] == "not_evaluated"


def test_disabled_efficiency_is_flagged(candidate, rows, tmp_path):
    row = dict(rows[0], is_efficient=False, strand_selection=None,
               end_stability=None, target_site_accessibility=None)
    result = import_sifi(save(tmp_path, [row]), candidate)
    assert result["records"][0]["is_efficient"] is False
    assert any("desactivada" in w for w in result["warnings"])


def test_reference_mismatch_warns_without_rewriting(candidate, rows, tmp_path):
    row = dict(rows[0], sirna_sequence="A" * 21)
    result = import_sifi(save(tmp_path, [row]), candidate)
    assert result["records"][0] == row
    assert any("no coincide con referencia" in w for w in result["warnings"])


def test_query_alias_is_explicit(candidate, rows, tmp_path):
    row = dict(rows[0], query_name="FASTA_ALIAS")
    result = import_sifi(save(tmp_path, [row]), candidate, {"query_id": "FASTA_ALIAS"})
    assert result["records"][0]["query_name"] == "FASTA_ALIAS"


@pytest.mark.parametrize("change", [{"priority_score": True}, {"metadata": []},
                                    {"evidence": {}}, {"schema_version": "2.0"}])
def test_candidate_contract(tmp_path, change):
    data = read_json(EXAMPLES / "candidate.json")
    with pytest.raises(ValueError):
        Candidate.load(save(tmp_path, dict(data, **change)))


def test_empty_cli_reports_missing_species(tmp_path):
    empty = save(tmp_path, [])
    out = tmp_path / "output"
    assert main(["--candidate", str(EXAMPLES / "candidate.json"), "--sifi", str(empty),
                 "--out", str(out)]) == 0
    report = read_json(out / "report.json")
    assert report["efficiency_summary"]["unique_sirnas"] == 0
    assert all(e["status"] == "not_evaluated" for e in report["offtarget_evaluations"])
    with (out / "evaluations.csv").open(encoding="utf-8", newline="") as handle:
        assert len(list(csv.DictReader(handle))) == 3


def test_cli_error(tmp_path):
    assert main(["--candidate", str(tmp_path / "missing"), "--sifi", "missing",
                 "--out", str(tmp_path / "out")]) == 2
    assert not (tmp_path / "out").exists()
