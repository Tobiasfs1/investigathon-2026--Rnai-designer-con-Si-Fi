"""Reconstruye las bandas verdes de si-Fi y resume observaciones en ventanas."""
from collections import defaultdict
from bisect import bisect_left, bisect_right
from statistics import mean, median

from .data_models import Candidate

METRICS = ("is_efficient", "strand_selection", "end_stability",
           "target_site_accessibility", "accessibility_value")
FLAG_METRICS = METRICS[:-1]


def green_positions(records: list[dict], sirna_size: int) -> tuple[set[int], list[str]]:
    """Reproduce show_plot.plot_design + general_helpers.get_target_data.

    Mantiene el orden de filas porque el código original resta únicamente los
    off-targets ya vistos al guardar snapshots de main_target_dict. El helper
    original evalÃºa is_off_target por truthiness, asÃ­ que null entra en main.
    """
    # Integer bitsets preserve ordered snapshots without repeatedly copying large sets.
    off_positions = 0
    main_positions = 0
    main_by_hit = {}
    main_seen: set[tuple[str, int]] = set()
    warnings = []
    for index, row in enumerate(records):
        flag = row.get("is_off_target")
        if flag is None:
            warnings.append(f"Fila {index}: is_off_target es null; si-Fi lo trata como false por truthiness de Python.")
        elif type(flag) is not bool:
            raise ValueError(f"Fila {index}: is_off_target debe ser booleano o null")
        if flag is True:
            start = row["sirna_position"]
            off_positions |= ((1 << sirna_size)-1) << (start-1)
        else:
            start = row["sirna_position"]
            key = (row["sirna_name"], start)
            if key in main_seen:
                continue
            main_seen.add(key)
            main_positions |= ((1 << sirna_size)-1) << (start-1)
            if "hit_name" not in row:
                warnings.append(f"Fila {index}: falta hit_name que requiere el helper original.")
                continue
            hit = row["hit_name"]
            try:
                hash(hit)
            except TypeError:
                warnings.append(f"Fila {index}: hit_name no es una clave válida para el helper original.")
                continue
            main_by_hit[hit] = main_positions & ~off_positions
    bits = 0
    for snapshot in main_by_hit.values():
        bits |= snapshot
    green = set()
    while bits:
        lowest = bits & -bits
        green.add(lowest.bit_length())
        bits ^= lowest
    return green, warnings


def _regions(positions: set[int]) -> list[dict]:
    if not positions:
        return []
    ordered = sorted(positions)
    regions = []
    start = previous = ordered[0]
    for position in ordered[1:]:
        if position != previous + 1:
            regions.append({"region_id": f"green-{len(regions)+1:03d}",
                            "start_1based": start, "end_1based_inclusive": previous,
                            "length_nt": previous - start + 1})
            start = position
        previous = position
    regions.append({"region_id": f"green-{len(regions)+1:03d}",
                    "start_1based": start, "end_1based_inclusive": previous,
                    "length_nt": previous - start + 1})
    return regions


def _explicit_green_positions(mask: dict, sequence_length: int) -> tuple[set[int], list[str]]:
    """Validate a lossless mask exported from original get_target_data()."""
    positions = mask.get("positions_1based")
    if not isinstance(positions, list):
        raise ValueError("green_mask debe contener positions_1based como lista")
    if any(type(position) is not int or position < 1 or position > sequence_length for position in positions):
        raise ValueError("green_mask contiene coordenadas fuera de la referencia o no enteras")
    if len(positions) != len(set(positions)):
        raise ValueError("green_mask contiene posiciones duplicadas")
    warnings = []
    if mask.get("source_function") != "general_helpers.get_target_data":
        warnings.append("green_mask explícita: source_function no acredita general_helpers.get_target_data.")
    return set(positions), warnings


def _sirnas(records: list[dict]) -> tuple[list[dict], list[str]]:
    grouped = defaultdict(list)
    for row in records:
        key = (row["sirna_name"], row["sirna_position"], row["sirna_sequence"])
        grouped[key].append(row)
    results, warnings = [], []
    for (name, start, sequence), rows in grouped.items():
        item = {"sirna_name": name, "start_1based": start,
                "end_1based_inclusive": start + len(sequence) - 1,
                "sequence": sequence, "hit_names": sorted({r["hit_name"] for r in rows
                                                              if isinstance(r.get("hit_name"), str)}),
                "source_row_count": len(rows), "source_rows": rows}
        for metric in METRICS:
            values = {r.get(metric) for r in rows if r.get(metric) is not None}
            item[metric] = next(iter(values)) if len(values) == 1 else None
            item[metric + "_conflict"] = len(values) > 1
            item[metric + "_missing_source_rows"] = sum(r.get(metric) is None for r in rows)
            if len(values) > 1:
                warnings.append(f"{name}: {metric} discrepa entre filas de hits; se excluye del resumen.")
            if item[metric + "_missing_source_rows"]:
                warnings.append(f"{name}: {metric} falta en {item[metric + '_missing_source_rows']} de {len(rows)} filas de hit.")
        results.append(item)
    return sorted(results, key=lambda item: (item["start_1based"], item["sirna_name"])), warnings


def _metric_summary(sirnas: list[dict]) -> dict:
    summary = {"sirna_start_count": len(sirnas)}
    for metric in FLAG_METRICS:
        summary[metric + "_true_count"] = sum(s[metric] is True for s in sirnas)
        summary[metric + "_false_count"] = sum(s[metric] is False for s in sirnas)
        summary[metric + "_unknown_count"] = sum(s[metric] is None for s in sirnas)
        summary[metric + "_missing_source_rows"] = sum(s[metric + "_missing_source_rows"] for s in sirnas)
    classified = summary["is_efficient_true_count"] + summary["is_efficient_false_count"]
    summary["is_efficient_classified_count"] = classified
    summary["is_efficient_true_fraction_among_classified"] = (
        summary["is_efficient_true_count"] / classified if classified else None)
    values = [s["accessibility_value"] for s in sirnas if s["accessibility_value"] is not None]
    summary["accessibility_value_count"] = len(values)
    summary["accessibility_value_missing_count"] = len(sirnas) - len(values)
    summary["accessibility_value_mean"] = mean(values) if values else None
    summary["accessibility_value_median"] = median(values) if values else None
    summary["accessibility_value_min"] = min(values) if values else None
    summary["accessibility_value_max"] = max(values) if values else None
    return summary


def _overlapping_offtargets(evaluations: list[dict], start: int, end: int,
                           query_length: int | None = None) -> list[dict]:
    matches = []
    for evaluation in evaluations:
        result = evaluation.get("result", {})
        source = result.get("source", {})
        run = result.get("run") or evaluation.get("manifest", {})
        criterion = evaluation.get("criterion", {})
        for alignment in criterion.get("alignments", []):
            a = alignment.get("query_start_1based")
            b = alignment.get("query_end_1based_inclusive")
            derived_end = False
            if b is None and type(a) is int and type(alignment.get("length")) is int:
                b = a + alignment["length"] - 1
                derived_end = True
            if (type(a) is not int or type(b) is not int or a < 1 or b < a
                    or (query_length is not None and b > query_length) or b < start or a > end):
                continue
            matches.append({"organism": evaluation["organism"],
                            "evaluation_status": evaluation["status"],
                            "criterion_status": criterion.get("status", "indeterminate"),
                            "database": run.get("database"),
                            "database_fingerprint": run.get("database_fingerprint"),
                            "result_sha256": run.get("result_sha256") or source.get("sha256"),
                            "window_overlap_start_1based": max(a, start),
                            "window_overlap_end_1based_inclusive": min(b, end),
                            "window_overlap_nt": min(b, end) - max(a, start) + 1,
                            **alignment,
                            "query_interval_derived_from_start_and_length": derived_end,
                            **({"query_end_1based_inclusive_derived": b} if derived_end else {})})
    return matches


def _offtarget_counts(matches: list[dict], evaluations: list[dict]) -> dict:
    from .offtarget_reports import SPECIES

    result = {}
    for organism in SPECIES:
        key = {"Homo sapiens": "human", "Zea mays": "maize", "Apis mellifera": "bee"}[organism]
        rows = [m for m in matches if m["organism"] == organism]
        entry = next((e for e in evaluations if e["organism"] == organism), None)
        result[key + "_evaluation_status"] = entry["status"] if entry else "not_evaluated"
        criterion = entry.get("criterion", {}) if entry else {}
        result[key + "_criterion_status"] = criterion.get("status", "not_evaluated")
        result[key + "_exact_21nt_overlapping_rows"] = sum(m.get("contains_exact_21nt") is True for m in rows)
        result[key + "_mismatch_overlapping_rows"] = sum(m.get("match_type") == "mismatch_alignment" for m in rows)
        result[key + "_unverified_overlapping_rows"] = sum(m.get("match_type") == "unverified" for m in rows)
    return result


def _unlinked_offtarget_counts(evaluations: list[dict], query_length: int) -> list[tuple[str, int]]:
    counts = []
    for evaluation in evaluations:
        missing = 0
        for alignment in evaluation.get("criterion", {}).get("alignments", []):
            start = alignment.get("query_start_1based")
            end = alignment.get("query_end_1based_inclusive")
            if end is None and type(start) is int and type(alignment.get("length")) is int:
                end = start + alignment["length"] - 1
            if (type(start) is not int or type(end) is not int or start < 1
                    or end < start or end > query_length):
                missing += 1
        if missing:
            counts.append((evaluation["organism"], missing))
    return counts


def explore_windows(report: dict, window_length: int, step: int = 1,
                    sirna_size: int | None = None, selected_start: int | None = None,
                    validation_only: bool = False, _include_prepared: bool = False) -> dict:
    """Resume observaciones si-Fi en ventanas enteras contenidas en zonas verdes.

    Las métricas se asignan al inicio 1-based de su siRNA. Una observación
    participa cuando todo su intervalo cabe en la ventana; no se interpola.
    """
    if type(window_length) is not int or window_length < 1:
        raise ValueError("window_length debe ser un entero positivo")
    if type(step) is not int or step < 1:
        raise ValueError("step debe ser un entero positivo")
    candidate_data = report.get("candidate")
    if not isinstance(candidate_data, dict):
        raise ValueError("El reporte no contiene candidate")
    reference = candidate_data.get("reference", {})
    candidate = Candidate(candidate_data, candidate_data.get("gene_id", "unknown"),
                          reference.get("sequence", ""))
    if not candidate.sequence:
        raise ValueError("El reporte no contiene secuencia de referencia")
    sifi = report.get("sifi", {})
    records = sifi.get("records", [])
    explicit_mask = sifi.get("green_mask")
    if (sifi.get("status") == "not_evaluated" or not records) and explicit_mask is None:
        raise ValueError("RNAi Designer no tiene resultados posicionales para reconstruir zonas verdes")
    lengths = {len(row.get("sirna_sequence", "")) for row in records
               if isinstance(row.get("sirna_sequence"), str) and row["sirna_sequence"]}
    run = sifi.get("run") or {}
    parameters = run.get("parameters") if isinstance(run.get("parameters"), dict) else {}
    configured_size = parameters.get("sirna_size")
    if type(configured_size) is not int and isinstance(explicit_mask, dict):
        configured_size = explicit_mask.get("sirna_size")
    warnings = []
    if sirna_size is None:
        if type(configured_size) is int:
            sirna_size = configured_size
        elif len(lengths) == 1:
            sirna_size = next(iter(lengths))
            warnings.append("El JSON de si-Fi no guarda sirna_size; se infirió de la longitud uniforme de sirna_sequence.")
        else:
            raise ValueError("Indicar --sirna-size; el manifiesto no lo registra y no hay longitud uniforme")
    if type(sirna_size) is not int or sirna_size < 1:
        raise ValueError("sirna_size debe ser un entero positivo")
    if lengths and lengths != {sirna_size}:
        raise ValueError(f"sirna_size={sirna_size} no coincide con las secuencias importadas: {sorted(lengths)}")
    query_id = run.get("query_id") or candidate.gene_id
    for i, row in enumerate(records):
        start = row.get("sirna_position")
        sequence = row.get("sirna_sequence")
        if row.get("query_name") != query_id:
            raise ValueError(f"Fila {i}: query_name no coincide con la ejecución")
        if type(start) is not int or start < 1 or start + sirna_size - 1 > len(candidate.sequence):
            raise ValueError(f"Fila {i}: coordenadas fuera de la referencia")
        if not isinstance(sequence, str) or len(sequence) != sirna_size:
            raise ValueError(f"Fila {i}: sirna_sequence no coincide con sirna_size")
    if explicit_mask is not None:
        if not isinstance(explicit_mask, dict):
            raise ValueError("sifi.green_mask debe ser objeto")
        if type(explicit_mask.get("sirna_size")) is int and explicit_mask["sirna_size"] != sirna_size:
            raise ValueError("sirna_size del green_mask no coincide con el manifiesto/importación")
        green, green_warnings = _explicit_green_positions(explicit_mask, len(candidate.sequence))
        green_mask_source = explicit_mask.get("source_function", "explicit sidecar")
        if explicit_mask.get("coordinate_system") != "1-based positions":
            green_warnings.append("green_mask no declara coordinate_system='1-based positions'.")
        reconstructed, _ = green_positions(records, sirna_size)
        if records and green != reconstructed:
            green_warnings.append(
                f"green_mask explícita difiere de la reconstrucción de los registros en "
                f"{len(green.symmetric_difference(reconstructed))} posiciones; verificar que ambos archivos sean de la misma corrida.")
        elif not records:
            green_warnings.append("No hay registros siRNA para contrastar las coordenadas de la mascara explicita; metricas de ventana incompletas.")
    else:
        green, green_warnings = green_positions(records, sirna_size)
        green_mask_source = "compatible reconstruction of general_helpers.get_target_data from source JSON rows"
    warnings.extend(green_warnings)
    warnings.extend(f"RNAi Designer: {warning}" for warning in sifi.get("warnings", []))
    warnings.append("La máscara verde reproduce la acumulación ordenada de main_target_dict de la versión original de si-Fi.")
    regions = _regions(green)
    sirnas, metric_warnings = _sirnas(records)
    warnings.extend(metric_warnings)
    starts = [s["start_1based"] for s in sirnas]
    windows = []
    evaluations = report.get("offtarget_evaluations", [])
    for organism, count in _unlinked_offtarget_counts(evaluations, len(candidate.sequence)):
        warnings.append(f"{organism}: {count} alineamientos carecen de coordenadas query válidas; no se asignaron a ventanas.")
    for region in (() if validation_only else regions):
        last_start = region["end_1based_inclusive"] - window_length + 1
        for start in range(region["start_1based"], last_start + 1, step):
            end = start + window_length - 1
            first = bisect_left(starts, start)
            last = bisect_right(starts, end - sirna_size + 1)
            observed = sirnas[first:last]
            matches = _overlapping_offtargets(evaluations, start, end, len(candidate.sequence))
            item = {"window_id": f'{region["region_id"]}:{start}-{end}',
                    "region_id": region["region_id"], "start_1based": start,
                    "end_1based_inclusive": end, "length_nt": window_length,
                    **_metric_summary(observed), **_offtarget_counts(matches, evaluations),
                    "offtarget_alignment_count": len(matches)}
            windows.append(item)
    selected = None
    if selected_start is not None:
        candidates = [w for w in windows if w["start_1based"] == selected_start]
        if not candidates:
            raise ValueError("--select-start no coincide con el inicio de una ventana generada; consultar windows.csv")
        selected = candidates[0]
        first = bisect_left(starts, selected_start)
        last = bisect_right(starts, selected["end_1based_inclusive"] - sirna_size + 1)
        selected = {**selected, "sirnas": sirnas[first:last],
                    "offtarget_matches": _overlapping_offtargets(evaluations, selected_start,
                                                                  selected["end_1based_inclusive"],
                                                                  len(candidate.sequence)),
                    "sequence": candidate.sequence[selected_start-1:selected["end_1based_inclusive"]],
                    "reference_source": reference.get("source"),
                    "reference_version": reference.get("version")}
    if not regions:
        warnings.append("No se reconstruyeron regiones verdes: confirmar que el JSON es de RNAi design y conserva is_off_target.")
    if (not validation_only and not windows) or (validation_only and not any(r['length_nt'] >= window_length for r in regions)):
        warnings.append("No hay ventanas de esa longitud completamente contenidas en las regiones verdes.")
    if not report.get("offtarget_evaluations"):
        warnings.append("No hay evaluaciones off-target vinculadas al reporte.")
    return {**({"_sirnas": sirnas, "_metric_warnings": metric_warnings} if _include_prepared else {}),
            "schema_version": "1.0", "candidate_id": candidate.gene_id,
            "homolog_id": candidate_data.get("tribolium_homolog_id"),
            "priority_score": candidate_data.get("priority_score"),
            "reference_source": reference.get("source"), "reference_version": reference.get("version"),
            "sifi_source": sifi.get("source"), "sifi_run": run,
            "definition": "Reconstrucción de show_plot.DrawPlot.plot_design y general_helpers.get_target_data; ver documentación.",
            "coordinate_system": "1-based inclusive",
            "green_mask_source": green_mask_source,
            "window_metric_assignment": "Observaciones si-Fi asignadas por sirna_position, inicio 1-based; una observación por siRNA. Sin interpolación.",
            "sirna_size": sirna_size, "window_length": window_length, "step": step,
            "green_regions": regions, "windows": windows, "selected_window": selected,
            "warnings": list(dict.fromkeys(warnings)),
            "offtarget_global_criterion": report.get("offtarget_criterion", {"status": "indeterminate"})}
