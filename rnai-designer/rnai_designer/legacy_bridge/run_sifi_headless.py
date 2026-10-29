# -*- coding: utf-8 -*-
"""Python 2.7 adapter around si-Fi's existing design methods; no GUI is created."""
from __future__ import print_function

import argparse
import collections
import json
import os
import subprocess
import sys
import shutil
import glob


def contiguous_runs(positions):
    ordered = sorted(positions)
    if not ordered:
        return []
    runs = []
    start = previous = ordered[0]
    for position in ordered[1:]:
        if position != previous + 1:
            runs.append({"start_1based": start, "end_1based_inclusive": previous,
                         "length_nt": previous - start + 1})
            start = position
        previous = position
    runs.append({"start_1based": start, "end_1based_inclusive": previous,
                 "length_nt": previous - start + 1})
    return runs


def parse_tab(path):
    rows = []
    with open(path, "r") as handle:
        for line in handle:
            rows.append(line.rstrip("\r\n").split("\t"))
    return rows


def relabel_records(records, selected):
    # The original design-mode data_to_json uses main_targets only for this flag.
    selected = set(selected)
    for row in records:
        row["is_off_target"] = row["hit_name"] not in selected
    return records


def memoize_efficiency(pipeline):
    original = pipeline.calculate_efficiency
    cache = {}
    def calculate(sequence, previous, accessibility):
        key = (sequence, previous, accessibility)
        if key not in cache:
            cache[key] = original(sequence, previous, accessibility)
        return cache[key]
    pipeline.calculate_efficiency = calculate


def bind_shared(directory, args, query_name, sequence):
    # The directory is private and lives only for this workflow. Check identities
    # and file signatures before any reuse; changed inputs fail closed.
    files = [args.bowtie, args.query_fasta, __file__]
    files += glob.glob(os.path.join(args.database_dir, args.database_prefix + ".*"))
    files += glob.glob(os.path.join(args.sifi_root, "*.py"))
    files += glob.glob(os.path.join(args.rnaplfold_dir, "*.exe"))
    files += glob.glob(os.path.join(args.rnaplfold_dir, "*.dll"))
    files += glob.glob(os.path.join(os.path.dirname(args.bowtie), "*.dll"))
    # Query FASTA paths differ between per-reference handoffs: sequence and ID
    # bind their content instead of their filesystem identity.
    files = [p for p in files if p != args.query_fasta]
    signatures = []
    for path in sorted(set(files)):
        stat = os.stat(path)
        signatures.append([os.path.abspath(path), stat.st_size, stat.st_mtime, stat.st_ctime])
    binding = {"query_name": query_name, "sequence": sequence, "files": signatures,
               "database": os.path.abspath(args.database_dir), "prefix": args.database_prefix,
               "sifi_root": os.path.abspath(args.sifi_root), "rnaplfold_dir": os.path.abspath(args.rnaplfold_dir),
               "sirna_size": args.sirna_size, "mismatches": args.mismatches,
               "accessibility_window": args.accessibility_window,
               "end_stability_threshold": args.end_stability_threshold,
               "target_accessibility_threshold": args.target_accessibility_threshold}
    manifest = os.path.join(directory, "binding.json")
    if os.path.isfile(manifest):
        with open(manifest) as handle:
            if json.load(handle) != binding:
                raise ValueError("Shared si-Fi inputs changed during workflow")
    else:
        with open(manifest + ".tmp", "w") as handle:
            json.dump(binding, handle)
        os.rename(manifest + ".tmp", manifest)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Ejecuta el cálculo si-Fi original sin popup ni gráfico")
    parser.add_argument("--sifi-root", required=True)
    parser.add_argument("--query-fasta", required=True)
    parser.add_argument("--database-prefix", required=True)
    parser.add_argument("--database-dir", required=True)
    parser.add_argument("--bowtie", required=True)
    parser.add_argument("--rnaplfold-dir", required=True)
    parser.add_argument("--out-json", required=True)
    parser.add_argument("--out-mask", required=True)
    parser.add_argument("--sirna-size", type=int, default=21)
    parser.add_argument("--mismatches", type=int, default=0)
    parser.add_argument("--accessibility-window", type=int, default=8)
    parser.add_argument("--end-stability-threshold", type=float, default=1.0)
    parser.add_argument("--target-accessibility-threshold", type=float, default=0.1)
    parser.add_argument("--main-target", action="append", default=[])
    parser.add_argument("--list-targets", action="store_true")
    parser.add_argument("--shared-dir", help="Internal per-workflow temporary shared results")
    args = parser.parse_args(argv)

    if args.sirna_size < 1 or args.mismatches < 0:
        parser.error("Tamaño siRNA y mismatches fuera de rango")
    sys.path.insert(0, os.path.abspath(args.sifi_root))
    from Bio import SeqIO
    import general_helpers
    import sifi_pipeline

    fasta_records = list(SeqIO.parse(args.query_fasta, "fasta"))
    if len(fasta_records) != 1:
        raise ValueError("Se requiere exactamente un registro FASTA")
    query = fasta_records[0]
    query_name = query.id
    query_sequence = str(query.seq)
    if len(query_sequence) < args.sirna_size:
        raise ValueError("La query es más corta que sirna_size")

    database_dir = os.path.abspath(args.database_dir)
    if not database_dir.endswith(os.sep):
        database_dir += os.sep
    bowtie_dir = os.path.dirname(os.path.abspath(args.bowtie))
    pipeline = sifi_pipeline.SifiPipeline(
        args.database_prefix, database_dir, args.query_fasta, args.sirna_size, args.mismatches,
        2, args.accessibility_window, args.rnaplfold_dir, bowtie_dir, 0, 2, 2,
        args.end_stability_threshold, args.target_accessibility_threshold,
        os.path.dirname(os.path.abspath(args.out_json)), 2, False)
    pipeline.len_seq = len(query_sequence)
    single_query = pipeline.create_single_fasta_file(query_name, query_sequence)
    sirna_fasta, tab_file = pipeline.create_sirnas(query_sequence, args.sirna_size)
    tab_rows = parse_tab(tab_file)

    alignment_file = os.path.abspath(args.out_json) + ".bowtie.tmp"
    database_index = database_dir + args.database_prefix
    command = [os.path.abspath(args.bowtie), "-a", "-v", str(args.mismatches), "-y",
               database_index, "-f", sirna_fasta, alignment_file]
    shared = args.shared_dir
    if shared:
        bind_shared(shared, args, query_name, query_sequence)
    cached_alignment = os.path.join(shared, "alignments.tsv") if shared else None
    counts_file = os.path.join(shared, "counts.json") if shared else None
    if shared and os.path.isfile(cached_alignment) and os.path.isfile(counts_file):
        shutil.copyfile(cached_alignment, alignment_file)
        with open(counts_file) as handle:
            counts = collections.Counter(json.load(handle))
        bowtie_rows = None
        bowtie_returncode = 0
    else:
        process = subprocess.Popen(command, cwd=bowtie_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            raise RuntimeError("Bowtie failed (%d): %s" % (process.returncode, stderr[-4000:]))
        with open(alignment_file, "r") as handle:
            bowtie_rows = pipeline.bowtie_to_lst(handle.readlines())
        counts = collections.Counter(row[2] for row in bowtie_rows)
        bowtie_returncode = process.returncode
        if shared:
            shutil.copyfile(alignment_file, cached_alignment)
            with open(counts_file + ".tmp", "w") as handle:
                json.dump(dict(counts), handle)
            os.rename(counts_file + ".tmp", counts_file)
    alignment_count = sum(counts.values())
    targets = sorted(counts.keys())

    if args.list_targets:
        print(json.dumps({"status": "targets_available", "query_id": query_name,
                          "alignment_rows": alignment_count,
                          "available_targets": [{"hit_name": name, "hit_rows": counts[name]}
                                                for name in targets]}, ensure_ascii=True))
        return 0

    selected = list(dict.fromkeys(args.main_target))
    unknown = sorted(set(selected) - set(targets))
    if unknown:
        raise ValueError("Main targets no encontrados en Bowtie: %s" % ", ".join(unknown))
    template = os.path.join(shared, "records.json") if shared else None
    if template and os.path.isfile(template):
        with open(template) as handle:
            records = relabel_records(json.load(handle), selected) if alignment_count else json.load(handle)
    else:
        pipeline.sirna_l = tab_rows
        lunp_data = pipeline.run_rnaplfold(query_name, single_query)
        if bowtie_rows is None:
            with open(alignment_file) as handle:
                bowtie_rows = pipeline.bowtie_to_lst(handle.readlines())
        input_rows = bowtie_rows if bowtie_rows else tab_rows
        # Repeated Bowtie hits use identical efficiency arguments. The original
        # method is pure for a fixed pipeline configuration; preserve every result.
        memoize_efficiency(pipeline)
        records = pipeline.data_to_json(query_name, input_rows, not bool(bowtie_rows), lunp_data, selected)
        if template:
            with open(template + ".tmp", "w") as handle:
                json.dump(records, handle)
            os.rename(template + ".tmp", template)
    with open(args.out_json, "w") as handle:
        json.dump(records, handle, indent=2)

    unused_off_targets, main_target_dict, unused_efficient, unused_histogram = \
        general_helpers.get_target_data(args.out_json, args.sirna_size)
    positions = set()
    for target_positions in main_target_dict.values():
        positions.update(target_positions)
    # DrawPlot.plot_design explicitly clears green positions when no main targets were selected.
    mask_cleared = not selected
    if mask_cleared:
        positions = set()
    main_starts = set(row.get("sirna_position") for row in records if not row.get("is_off_target"))
    off_target_starts = set(row.get("sirna_position") for row in records if row.get("is_off_target"))
    cross_target_starts = sorted(main_starts & off_target_starts)
    mask_warnings = []
    if cross_target_starts:
        mask_warnings.append("Main-target siRNA positions also have off-target alignments; the original get_target_data mask subtracts off-target positions and may shorten/suppress green runs.")
    mask = {"schema_version": "1.0", "source_function": "general_helpers.get_target_data",
            "coordinate_system": "1-based positions", "sirna_size": args.sirna_size,
            "positions_1based": sorted(positions), "green_regions": contiguous_runs(positions),
            "main_targets_selected": selected,
            "main_offtarget_overlap_start_count": len(cross_target_starts),
            "main_offtarget_overlap_starts_1based": cross_target_starts,
            "warnings": mask_warnings,
            "plot_mask_empty_due_to_no_main_targets": mask_cleared}
    with open(args.out_mask, "w") as handle:
        json.dump(mask, handle, indent=2, sort_keys=True)
    print(json.dumps({"status": "completed", "query_id": query_name,
                      "records": len(records), "alignment_rows": alignment_count,
                      "main_targets": selected, "green_regions": mask["green_regions"],
                      "bowtie_returncode": bowtie_returncode}, ensure_ascii=True))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print("ERROR: %s" % error, file=sys.stderr)
        sys.exit(2)
