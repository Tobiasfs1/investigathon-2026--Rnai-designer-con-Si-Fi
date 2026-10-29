"""Adaptador al Bowtie instalado de si-Fi; no calcula eficiencia ni optimiza secuencias."""
from pathlib import Path
import hashlib
import re

from .data_models import Candidate, ImportedResult
from .sifi_adapter import provenance


def sequence_hash(candidate: Candidate) -> str:
    return hashlib.sha256(candidate.sequence.upper().encode("ascii")).hexdigest()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def write_queries(path: Path, candidate: Candidate) -> int:
    """Enumera todas las ventanas 21nt como entrada del alineador existente."""
    sequence = candidate.sequence.upper()
    if set(sequence) - set("ACGT") or len(sequence) < 21:
        raise ValueError("Búsqueda exacta automática: requiere referencia A/C/G/T de al menos 21 nt")
    count = len(sequence) - 20
    with path.open("w", encoding="ascii") as handle:
        for start in range(count):
            handle.write(f">sirna{start+1}\n{sequence[start:start+21]}\n")
    return count


def parse_bowtie(path: Path, candidate: Candidate, run: dict) -> ImportedResult:
    rows, warnings = [], []
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            columns = line.rstrip("\r\n").split("\t")
            if len(columns) not in (7, 8):
                raise ValueError(f"Bowtie línea {number}: se requieren 7 u 8 columnas originales")
            name, strand, reference, offset, sequence, quality, other = columns[:7]
            match = re.fullmatch(r"sirna([1-9][0-9]*)", name)
            if not match or strand not in ("+", "-") or not reference or not offset.isdigit():
                raise ValueError(f"Bowtie línea {number}: ID, hebra, referencia o posición inválidos")
            start = int(match.group(1))
            sequence = sequence.upper()
            if set(sequence) - set("ACGT") or not sequence or start + len(sequence) - 1 > len(candidate.sequence):
                raise ValueError(f"Bowtie línea {number}: secuencia/intervalo inválidos")
            expected = candidate.sequence[start-1:start-1+len(sequence)].upper()
            if strand == "-":
                expected = reverse_complement(expected)
            if sequence != expected:
                raise ValueError(f"Bowtie línea {number}: read no corresponde a la referencia candidata")
            descriptor = columns[7] if len(columns) == 8 else ""
            rows.append({"query_name": run.get("query_id", candidate.gene_id), "sirna_name": name,
                         "sirna_position": start, "sirna_sequence": sequence,
                         "hit_name": reference, "reference_strand_pos": int(offset), "strand": strand,
                         "mismatches": descriptor, "quality": quality, "other_alignments": other,
                         "raw_bowtie_line": line.rstrip("\r\n")})
    warnings.append("Alineamientos Bowtie: no incluyen métricas de eficiencia de RNAi Designer.")
    if run.get("database_version") in (None, "unknown"):
        warnings.append("Versión biológica de la base desconocida; el fingerprint identifica los índices locales.")
    if not run.get("database_fingerprint", {}).get("declared_source"):
        warnings.append("Procedencia original de la base no documentada; verificar accesiones y taxonomía.")
    warnings.append("El criterio se limita a las secuencias indexadas y a bases resueltas A/C/G/T; no cubre referencias ausentes/ambiguas.")
    return {"source": provenance(path), "run": run, "records": rows, "warnings": warnings}
