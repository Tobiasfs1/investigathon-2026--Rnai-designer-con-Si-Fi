"""Pruebas de contratos de ejecución y, cuando está instalado, Bowtie real."""
from pathlib import Path
from dataclasses import replace
import copy
import json
import subprocess
import sys

import pytest

from rnai_designer.automation import run_offtargets
from rnai_designer.bowtie_adapter import parse_bowtie, sequence_hash, write_queries, reverse_complement
from rnai_designer.config import Settings, inventory, INDEX_PARTS
from rnai_designer.data_models import Candidate, read_json
from rnai_designer.external import execute
from rnai_designer.handoff import prepare, collect
from rnai_designer.offtarget_assessment import assess, mismatch_count, describe_alignment
from rnai_designer.offtarget_reports import load_evaluations
from rnai_designer.sifi_adapter import provenance

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


@pytest.fixture
def candidate():
    return Candidate.load(EXAMPLES / "candidate.json")


def complete_empty_evaluation(candidate):
    checksum = "example-sha256-fixture"
    run = {"query_sha256": sequence_hash(candidate), "result_sha256": checksum,
           "completion": {"status": "completed", "returncode": 0},
           "parameters": {"sirna_size": 21, "mismatches": 0, "alignment_policy": "-v",
                          "gapped": False, "all_alignments": True, "both_strands": True},
           "coverage": {"expected_windows": 60, "submitted_windows": 60, "processed_windows": 60}}
    return {"status": "imported_unvalidated", "results_format": "bowtie_legacy_tab",
            "result": {"records": [], "run": run, "source": {"sha256": checksum}}}


def test_verified_negative_is_conditional(candidate):
    result = assess(complete_empty_evaluation(candidate), candidate)
    assert result["status"] == "criterion_met"
    assert "no garantiza" in result["limitation"]


@pytest.mark.parametrize("mutation", [
    lambda e: e.update(results_format="sifi_json"),
    lambda e: e["result"]["run"].update(completion={"status": "timeout", "returncode": None}),
    lambda e: e["result"]["run"]["parameters"].update(both_strands=False),
    lambda e: e["result"]["run"]["parameters"].update(all_alignments=False),
    lambda e: e["result"]["run"]["parameters"].update(mismatches=1),
    lambda e: e["result"]["run"]["coverage"].update(processed_windows=59),
    lambda e: e["result"]["run"].update(query_sha256="OTHER"),
    lambda e: e["result"]["source"].update(sha256="OTHER"),
    lambda e: e["result"]["run"].update(completion=[]),
    lambda e: e["result"]["run"].pop("coverage"),
])
def test_incomplete_never_negative(candidate, mutation):
    entry = complete_empty_evaluation(candidate)
    mutation(entry)
    assert assess(entry, candidate)["status"] == "indeterminate"


@pytest.mark.parametrize("descriptor,expected", [(0, 0), ("", 0), (None, None), ("0", None),
                                                  ("3:A>G", 1), ("3:A>G,7:T>C", 2),
                                                  ("23:A>G", None), ("3:A>A", None), (False, None)])
def test_mismatch_descriptors(descriptor, expected):
    assert mismatch_count(descriptor, 21) == expected


def test_exact_does_not_require_negative_proof(candidate):
    row = read_json(EXAMPLES / "human.json")[0]
    entry = {"status": "imported_unvalidated", "result": {"records": [row]}}
    assert assess(entry, candidate)["status"] == "exact_match_detected"
    row["mismatches"] = "3:A>G"
    result = assess(entry, candidate)
    assert result["status"] == "indeterminate" and result["mismatch_rows"] == 1


def test_invalid_sequence_cannot_prove_exact(candidate):
    row = read_json(EXAMPLES / "human.json")[0]
    row["sirna_sequence"] = "A" * 21
    assert describe_alignment(row, candidate)["match_type"] == "unverified"
    row["sirna_sequence"] = candidate.sequence[:21]
    row["cigar"] = "10M1I10M"
    assert describe_alignment(row, candidate)["match_type"] == "unverified"


def test_failed_run_and_missing_analysis(candidate):
    assert assess({"status": "error", "error": "failed"}, candidate)["status"] == "error"
    assert assess({"status": "not_evaluated"}, candidate)["status"] == "not_evaluated"


def test_raw_reverse_strand_and_positions(candidate, tmp_path):
    path = tmp_path / "alignments.tsv"
    reverse = reverse_complement(candidate.sequence[:21])
    path.write_text(f"sirna1\t-\tref\t10\t{reverse}\t{'I'*21}\t0\t\n", encoding="utf-8")
    result = parse_bowtie(path, candidate, {})
    alignment = describe_alignment(result["records"][0], candidate)
    assert alignment["contains_exact_21nt"] is True
    assert alignment["reference_end_0based_exclusive"] == 31
    assert result["records"][0]["mismatches"] == ""
    assert result["source"] == provenance(path)


def test_truncated_bowtie_rejected(candidate, tmp_path):
    path = tmp_path / "alignments.tsv"
    path.write_text("sirna1\t+\tref", encoding="utf-8")
    with pytest.raises(ValueError):
        parse_bowtie(path, candidate, {})


def test_queries_and_ambiguous_input(candidate, tmp_path):
    assert write_queries(tmp_path / "reads.fasta", candidate) == 60
    assert (tmp_path / "reads.fasta").read_text().count(">") == 60
    with pytest.raises(ValueError):
        write_queries(tmp_path / "bad.fasta", replace(candidate, sequence="N" * 80))


def test_configuration_environment_and_overrides(tmp_path, monkeypatch):
    monkeypatch.setenv("LOCALAPPDATA", str(tmp_path))
    monkeypatch.delenv("SIFI_HOME", raising=False)
    monkeypatch.delenv("SIFI_DB_DIR", raising=False)
    settings = Settings.load()
    assert settings.database_dir == tmp_path / "siFi2015" / "Databases"
    assert not inventory(settings)["exists"]
    settings = Settings.load(database_dir=tmp_path / "alternate")
    assert settings.database_dir == tmp_path / "alternate"


def test_inventory_preserves_dotted_prefix_and_incomplete_index(tmp_path):
    settings = Settings(tmp_path, None, None, None, {})
    for part in INDEX_PARTS[:-1]:
        (tmp_path / ("reference.v1" + part + ".ebwt")).write_bytes(b"fixture")
    result = inventory(settings)
    assert result["indexes"][0]["prefix"] == "reference.v1"
    assert result["indexes"][0]["complete"] is False
    assert result["indexes"][0]["missing"] == ["reference.v1.rev.2.ebwt"]


def test_inventory_reports_nested_indexes_without_binding_them(tmp_path):
    settings = Settings(tmp_path, None, None, None,
                       {"Homo sapiens": {"prefix": "human"}})
    nested = tmp_path / "Off_Targets"
    nested.mkdir()
    for part in INDEX_PARTS:
        (nested / ("human" + part + ".ebwt")).write_bytes(b"fixture")
    result = inventory(settings)
    nested_index = next(index for index in result["indexes"] if index["directory"] == "Off_Targets")
    assert nested_index["complete"] is True
    assert nested_index["organism_labels"] == []
    assert "Off_Targets/human.1.ebwt" in [entry["name"] for entry in result["files"]]


def test_prepare_does_not_overwrite_user_manifest(candidate, tmp_path):
    first = prepare(candidate, tmp_path)
    (tmp_path / "run.json").write_text('{"completed_by_user":true}', encoding="utf-8")
    assert prepare(candidate, tmp_path) == first
    assert read_json(tmp_path / "run.json") == {"completed_by_user": True}
    assert first["query_id"] in (tmp_path / "query.fasta").read_text()


def test_collect_checks_identity_and_ambiguity(candidate, tmp_path):
    rows = read_json(EXAMPLES / "sifi_result.json")
    (tmp_path / "one.json").write_text(json.dumps(rows), encoding="utf-8")
    (tmp_path / "unrelated.json").write_text('{"other":true}', encoding="utf-8")
    path, imported = collect(candidate, None, [tmp_path])
    assert path.name == "one.json" and imported["records"] == rows
    different = copy.deepcopy(rows)
    different[0]["accessibility_value"] = .9
    (tmp_path / "two.json").write_text(json.dumps(different), encoding="utf-8")
    with pytest.raises(ValueError, match="Varios"):
        collect(candidate, None, [tmp_path])


def test_collect_cannot_associate_empty_json(candidate, tmp_path):
    (tmp_path / "empty.json").write_text("[]", encoding="utf-8")
    with pytest.raises(ValueError, match="No hay JSON"):
        collect(candidate, None, [tmp_path])


def test_external_timeout(tmp_path):
    result = execute([sys.executable, "-c", "import time; time.sleep(10)"], tmp_path, .05)
    assert result["status"] == "timeout" and result["returncode"] is None
    assert (tmp_path / "stderr.log").exists()


def test_tolerant_import_keeps_species_error(candidate, tmp_path):
    manifest = read_json(EXAMPLES / "offtargets.json")
    path = tmp_path / "offtargets.json"
    path.write_text(json.dumps(manifest), encoding="utf-8")
    evaluations = load_evaluations(path, candidate, tolerant=True)
    assert all(e["status"] == "error" for e in evaluations)
    assert all(e["criterion"]["status"] == "error" for e in evaluations)


@pytest.mark.integration
def test_real_bowtie_with_synthetic_indexes(candidate, tmp_path):
    try:
        installed = Settings.load()
    except ValueError:
        pytest.skip("Sin instalación local de Bowtie")
    if installed.bowtie is None or not installed.bowtie.is_file():
        pytest.skip("Bowtie no disponible")
    build = installed.bowtie.with_name("bowtie-build" + installed.bowtie.suffix)
    if not build.is_file():
        pytest.skip("bowtie-build no disponible")
    databases = tmp_path / "databases"
    databases.mkdir()
    # Índices minúsculos de prueba; jamás escribir en las bases instaladas.
    sequence = "ACGTTGCAACGATCGTACGTA"
    positive = replace(candidate, sequence=sequence)
    for prefix, reference in (("positive", sequence + "NNNNN" + reverse_complement(sequence)),
                              ("negative", "A" * 80)):
        fasta = tmp_path / (prefix + ".fasta")
        fasta.write_text(f">synthetic_reference\n{reference}\n", encoding="ascii")
        process = subprocess.run([str(build), str(fasta), str(databases / prefix)], capture_output=True, timeout=30)
        assert process.returncode == 0, process.stderr.decode(errors="replace")
        bindings = {name: {"prefix": prefix, "source": "SYNTHETIC TEST", "version": "fixture-1"}
                    for name in installed.databases}
        settings = replace(installed, database_dir=databases, databases=bindings)
        before = {p.name: (p.stat().st_size, p.stat().st_mtime_ns) for p in databases.iterdir()}
        evaluations = run_offtargets(positive, settings, tmp_path / prefix, timeout=30)
        expected = "exact_match_detected" if prefix == "positive" else "criterion_met"
        assert all(e["criterion"]["status"] == expected for e in evaluations), evaluations
        if prefix == "positive":
            assert {r["strand"] for r in evaluations[0]["result"]["records"]} == {"+", "-"}
        after = {p.name: (p.stat().st_size, p.stat().st_mtime_ns) for p in databases.iterdir()}
        assert before == after
        reimported = load_evaluations(tmp_path / prefix / "offtargets.json", positive)
        assert [e["criterion"]["status"] for e in reimported] == [expected] * 3
        with pytest.raises(ValueError, match="salida nueva"):
            run_offtargets(positive, settings, tmp_path / prefix, timeout=30)
