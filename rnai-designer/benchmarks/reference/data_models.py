"""Contratos de entrada y lectura JSON estricta."""
from dataclasses import dataclass
from pathlib import Path
import json
import math
from typing import Any, TypedDict


class ImportContext(TypedDict, total=False):
    candidate_source: dict[str, str]
    run_source: dict[str, str]


class ImportedResult(ImportContext):
    source: dict[str, str]
    run: dict[str, Any] | None
    records: list[dict[str, Any]]
    warnings: list[str]


def read_json(path: Path | str) -> Any:
    def invalid(value):
        raise ValueError(f"Constante JSON no finita: {value}")
    def finite_float(value):
        number = float(value)
        if not math.isfinite(number):
            raise ValueError(f"Número fuera de rango JSON: {value}")
        return number

    def unique_keys(pairs):
        data = {}
        for key, value in pairs:
            if key in data:
                raise ValueError(f"Clave JSON repetida: {key}")
            data[key] = value
        return data

    return json.loads(Path(path).read_text(encoding="utf-8-sig"), parse_constant=invalid,
                      parse_float=finite_float, object_pairs_hook=unique_keys)


def require_text(data: dict, key: str) -> str:
    value = data.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Campo requerido de texto: {key}")
    return value


@dataclass
class Candidate:
    original: dict[str, Any]
    gene_id: str
    sequence: str

    @classmethod
    def from_fields(cls, gene_id: str, sequence: str, reference_source: str | None = None,
                    reference_version: str | None = None, homolog_id: str | None = None,
                    priority_score: float | None = None):
        gene_id = require_text({"gene_id": gene_id}, "gene_id")
        if not isinstance(sequence, str):
            raise ValueError("sequence debe ser texto")
        lines = sequence.splitlines()
        headers = [line for line in lines if line.lstrip().startswith(">")]
        if len(headers) > 1:
            raise ValueError("El workflow acepta una sola secuencia FASTA por corrida")
        if headers:  # also accept a pasted single-record FASTA
            sequence = "".join(line.strip() for line in lines if not line.lstrip().startswith(">"))
        else:
            sequence = "".join(sequence.split())
        sequence = sequence.upper()
        if not sequence or set(sequence) - set("ACGTURYSWKMBDHVN"):
            raise ValueError("Secuencia invÃ¡lida: se requiere ADN IUPAC sin espacios ni caracteres ambiguos no permitidos")
        if reference_source is not None:
            require_text({"reference_source": reference_source}, "reference_source")
        if reference_version is not None:
            require_text({"reference_version": reference_version}, "reference_version")
        if homolog_id is not None:
            require_text({"homolog_id": homolog_id}, "homolog_id")
        if priority_score is not None and (type(priority_score) not in (int, float)):
            raise ValueError("priority_score debe ser número o null")
        original = {
            "schema_version": "1.0", "gene_id": gene_id,
            "tribolium_homolog_id": homolog_id,
            "reference": {"sequence": sequence, "source": reference_source, "version": reference_version},
            "priority_score": priority_score,
            "metadata": {"simulated": False, "input_mode": "rnai-designer-workflow"},
            "evidence": [],
        }
        return cls(original, gene_id, sequence)

    @classmethod
    def load(cls, path):
        data = read_json(path)
        if not isinstance(data, dict) or data.get("schema_version") != "1.0":
            raise ValueError("Candidato: schema_version debe ser 1.0")
        gene_id = require_text(data, "gene_id")
        if data.get("tribolium_homolog_id") is not None:
            require_text(data, "tribolium_homolog_id")
        reference = data.get("reference")
        if not isinstance(reference, dict):
            raise ValueError("reference debe ser objeto")
        sequence = require_text(reference, "sequence")
        if set(sequence.upper()) - set("ACGTURYSWKMBDHVN"):
            raise ValueError("Secuencia inválida: se requiere IUPAC sin espacios")
        for key in ("source", "version"):
            if reference.get(key) is not None:
                require_text(reference, key)
        if "priority_score" not in data:
            raise ValueError("Falta priority_score del módulo anterior")
        if type(data["priority_score"]) not in (int, float) and data["priority_score"] is not None:
            raise ValueError("priority_score debe ser número o null")
        for key in ("metadata", "evidence"):
            if key not in data:
                raise ValueError(f"Falta {key}")
        if not isinstance(data["metadata"], dict) or not isinstance(data["evidence"], list):
            raise ValueError("metadata debe ser objeto y evidence debe ser lista")
        return cls(data, gene_id, sequence)
