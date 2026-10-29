"""Two-argument entry point; all analysis is delegated to workflow."""
import argparse
from pathlib import Path
import re
import sys
import unicodedata

from .data_models import Candidate

# Legacy Windows tools still need MAX_PATH even when Python supports long paths.
MAX_PATH_CHARS = 259
WORKFLOW_PATH_RESERVE = 140  # Nested reference, staging and off-target artifacts.


def short_identifier(gene_id):
    text = unicodedata.normalize("NFKD", gene_id).encode("ascii", "ignore").decode("ascii")
    text = re.sub(r"[^A-Za-z0-9_-]+", "_", text).strip("_-")[:12] or "gen"
    if text.upper() in {"CON", "PRN", "AUX", "NUL", *[f"COM{i}" for i in range(1, 10)], *[f"LPT{i}" for i in range(1, 10)]}:
        text = "g_" + text[:10]
    return text


def check_paths(fasta, out):
    paths = [(fasta.resolve(), 0), (out.resolve(), WORKFLOW_PATH_RESERVE)]
    for path, reserve in paths:
        if len(str(path).encode("utf-16-le")) // 2 + reserve > MAX_PATH_CHARS:
            raise ValueError("Ruta demasiado larga para herramientas Windows; ejecutar desde una carpeta mas corta: " + str(path))


def reserve_run(root, identifier, sequence):
    """Exclusive creation protects both namespaces, including concurrent invocations."""
    inputs, outputs = root / "inputs", root / "output"
    index = 1
    while True:
        name = f"{identifier}_{index}"
        fasta, out = inputs / (name + ".fasta"), outputs / name
        check_paths(fasta, out)
        if fasta.exists() or fasta.is_symlink() or out.exists() or out.is_symlink():
            index += 1
            continue
        inputs.mkdir(parents=True, exist_ok=True)
        outputs.mkdir(parents=True, exist_ok=True)
        try:
            out.mkdir()
        except FileExistsError:
            index += 1
            continue
        try:
            with fasta.open("x", encoding="ascii", newline="\n") as handle:
                handle.write(f">{name}\n")
                for offset in range(0, len(sequence), 60):
                    handle.write(sequence[offset:offset+60] + "\n")
        except FileExistsError:
            out.rmdir()  # Only our newly reserved, empty directory.
            index += 1
            continue
        except BaseException:
            out.rmdir()
            raise
        return fasta, out


def simple_command(argv, workflow):
    parser = argparse.ArgumentParser(prog="rnai-designer", description="Analizar un gen sin preparar archivos ni rutas")
    parser.add_argument("gene_id", help="Nombre original del gen")
    parser.add_argument("sequence", help="Secuencia IUPAC completa, entre comillas")
    args = parser.parse_args(argv)
    try:
        # Use the same IUPAC contract and normalization as the advanced workflow.
        candidate = Candidate.from_fields(args.gene_id, "".join(args.sequence.split()))
        fasta, out = reserve_run(Path.cwd(), short_identifier(candidate.gene_id), candidate.sequence)
    except (ValueError, OSError) as error:
        print("Error: " + str(error), file=sys.stderr)
        return 2
    return workflow(["--gene-id", candidate.gene_id, "--sequence-file", str(fasta),
                     "--out", str(out), "--run-sifi-headless"], input_metadata={
                         "input_mode": "rnai-designer-simple",
                         "original_gene_name": candidate.gene_id,
                         "input_fasta": str(fasta), "run_identifier": out.name})
