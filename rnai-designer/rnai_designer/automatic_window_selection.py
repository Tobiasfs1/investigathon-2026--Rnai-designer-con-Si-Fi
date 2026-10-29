"""Exhaustive, auditable selection by evidence quality (not biological efficacy)."""
from bisect import bisect_left, bisect_right
import json
from functools import lru_cache
from .compact_windows import CompactWindows
from .priority_scoring import length_scores
from .window_index import WindowIndex, OfftargetIndex

from .window_explorer import (
    METRICS, _offtarget_counts, _overlapping_offtargets, _sirnas, explore_windows,
)


@lru_cache(maxsize=65536)
def _gc_percent(sequence: str) -> float | None:
    canonical = [base for base in sequence.upper() if base in "ACGT"]
    if not canonical:
        return None
    return 100.0 * sum(base in "GC" for base in canonical) / len(canonical)


def _efficacy_priority(contained: list[dict], efficient_count: int,
                       evaluable_count: int) -> tuple[float | None, dict]:
    """Exploratory rank: 75% original si-Fi efficiency, 25% GC; diversity is descriptive."""
    gc_values = [_gc_percent(item["sequence"]) for item in contained]
    gc_values = [value for value in gc_values if value is not None]
    gc_fraction = (sum(30.0 <= value <= 52.0 for value in gc_values) / len(gc_values)
                   if gc_values else None)
    diversity = (len({item["sequence"] for item in contained}) / len(contained)
                 if contained else None)
    sifi_fraction = efficient_count / evaluable_count if evaluable_count else None
    components = {"sifi_is_efficient_fraction": sifi_fraction,
                  "sirna_gc_30_52_fraction": gc_fraction,
                  "unique_sirna_sequence_fraction": diversity}
    if None in components.values():
        return None, components
    return 100.0 * (0.75 * sifi_fraction + 0.25 * gc_fraction), components


def _coverage(sirnas: list[dict], start: int, end: int, sirna_size: int) -> tuple[int, float, int, float]:
    possible_starts = max(0, end - start - sirna_size + 2)
    positions = {s["start_1based"] for s in sirnas}
    represented = len(positions)
    start_coverage = represented / possible_starts if possible_starts else 0.0
    intervals = sorted((max(start, s["start_1based"]), min(end, s["end_1based_inclusive"])) for s in sirnas)
    covered_nt = 0
    left = right = None
    for a, b in intervals:
        if left is None:
            left, right = a, b
        elif a > right + 1:
            covered_nt += right - left + 1
            left, right = a, b
        else:
            right = max(right, b)
    if left is not None:
        covered_nt += right - left + 1
    return represented, start_coverage, covered_nt, covered_nt / (end - start + 1)


def _window_metrics(region: dict, start: int, length: int, sequence: str,
                    raw_records: list[dict], raw_starts: list[int], sirnas: list[dict],
                    sirna_starts: list[int], sirna_size: int, evaluations: list[dict],
                    traceability: float, duplicate_extras: dict[int, int], index: WindowIndex | None = None,
                    off_index: OfftargetIndex | None = None) -> dict:
    end = start + length - 1
    max_start = end - sirna_size + 1
    first = bisect_left(sirna_starts, start)
    last = bisect_right(sirna_starts, max_start) if max_start >= start else first
    count = lambda key: index.count(key, first, last)
    # All source lengths were validated once; interval containment follows from start.
    raw_count = bisect_right(raw_starts, max_start)-bisect_left(raw_starts, start) if max_start >= start else 0
    unique = last-first
    cells = unique * len(METRICS)
    known_cells = count('known_cells')
    complete = count('complete')
    missing_pct = 100.0 * (cells - known_cells) / cells if cells else 100.0
    classified = count("is_efficient_true") + count("is_efficient_false")
    efficient = count("is_efficient_true")
    access, ordered_access, unique_sequences = index.range_data(first, last)
    exact_duplicates = index.duplicate_count(start, max_start)
    conflict_count = count('conflicts')
    represented, start_coverage, covered_nt, base_coverage = index.coverage(first, last, length, sirna_size)
    # Coordinates are syntactically checked by explore_windows; this also checks each imported
    # sequence against the corresponding reference interval, which the input format does not require.
    matching = count("matching")
    coordinate_consistency = matching / unique if unique else 0.0
    metric_completeness = known_cells / cells if cells else 0.0
    duplicate_integrity = (raw_count - exact_duplicates) / raw_count if raw_count else 0.0
    # Equal-weight mean of five evidence audit dimensions. No si-Fi efficacy/accessibility
    # value is part of the score. Each component is exposed for independent review.
    components = {
        "record_metric_completeness": metric_completeness,
        "coordinate_sequence_consistency": coordinate_consistency,
        "position_start_coverage": start_coverage,
        "exact_duplicate_integrity": duplicate_integrity,
        "traceability": traceability,
    }
    quality_score = 100.0 * sum(components.values()) / len(components)
    gc_valid = count('gc_valid')
    efficacy_priority_components = {
        'sifi_is_efficient_fraction': efficient / classified if classified else None,
        'sirna_gc_30_52_fraction': count('gc_target') / gc_valid if gc_valid else None,
        'unique_sirna_sequence_fraction': unique_sequences / unique if unique else None,
    }
    efficacy_priority_score = (None if None in efficacy_priority_components.values() else
                              100.0 * (0.75 * efficacy_priority_components['sifi_is_efficient_fraction'] +
                                       0.25 * efficacy_priority_components['sirna_gc_30_52_fraction']))
    alignment_count, species_counts = off_index.counts(start, end)
    row = {
        "window_id": f'{region["region_id"]}:{start}-{end}', "region_id": region["region_id"],
        "start_1based": start, "end_1based_inclusive": end, "length_nt": length,
        "gc_percent": index.gc(start, end),
        "sirna_record_count": raw_count, "sirna_unique_count": unique,
        "sirna_unique_sequences_count": unique_sequences,
        "is_efficient_true_count": efficient, "is_efficient_false_count": count("is_efficient_false"),
        "is_efficient_evaluable_count": classified,
        "is_efficient_true_fraction_among_evaluable": efficient / classified if classified else None,
        "sirna_complete_count": complete, "missing_data_percent": missing_pct,
        "metric_cells_expected": cells, "metric_cells_available": known_cells,
        "represented_start_positions": represented, "possible_start_positions": max(0, length - sirna_size + 1),
        "position_start_coverage_fraction": start_coverage,
        "covered_nucleotides": covered_nt, "nucleotide_coverage_fraction": base_coverage,
        "accessibility_value_count": len(access),
        "accessibility_value_mean": sum(access) / len(access) if access else None,
        "accessibility_value_median": ordered_access[len(access)//2] if len(access) % 2 else
            ((ordered_access[len(access)//2 - 1] + ordered_access[len(access)//2]) / 2 if access else None),
        "accessibility_value_min": ordered_access[0] if access else None,
        "accessibility_value_max": ordered_access[-1] if access else None,
        "exact_duplicate_record_count": exact_duplicates,
        "repeated_hit_record_count": max(0, raw_count - unique),
        "sirnas_with_metric_conflicts_count": conflict_count,
        "coordinate_sequence_match_count": matching,
        "coordinate_sequence_mismatch_count": unique - matching,
        "quality_score": quality_score, "quality_components": components,
        "efficacy_priority_score": efficacy_priority_score,
        "efficacy_priority_components": efficacy_priority_components,
        "data_sufficiency_status": ("insufficient_no_complete_sirna_records" if unique == 0 else
                                     "insufficient_no_evaluable_sifi_classifications" if classified == 0 else
                                     "partial_original_metrics" if complete < unique else
                                     "priority_components_available" if efficacy_priority_score is not None else
                                     "insufficient_gc_or_diversity_data"),
        "offtarget_alignment_count": alignment_count,
        "offtarget_by_species": species_counts,
    }
    row.update(row["offtarget_by_species"])
    for metric in METRICS:
        if metric != 'accessibility_value':
            row[metric + '_true_count'] = count(metric + '_true')
            row[metric + '_false_count'] = count(metric + '_false')
            row[metric + '_unknown_count'] = unique - count(metric + '_known')
    return row


def select_windows(report: dict, min_length: int = 180, max_length: int = 250,
                   sirna_size: int | None = None, length_penalty_weight: float = 30.0) -> dict:
    """Enumerate every contained window with a one-nucleotide stride."""
    if type(min_length) is not int or type(max_length) is not int or min_length < 1 or max_length < min_length:
        raise ValueError("Se requiere 1 <= min_length <= max_length")
    length_scores(None, 180, length_penalty_weight)
    # Reuse the existing strict input validation and authoritative green-mask reconstruction.
    baseline = explore_windows(report, min_length, 1, sirna_size, validation_only=True, _include_prepared=True)
    candidate_data = report["candidate"]
    reference = candidate_data["reference"]
    sequence = reference["sequence"]
    records = report["sifi"]["records"]
    sirnas, metric_warnings = baseline.pop("_sirnas"), baseline.pop("_metric_warnings")
    for item in (sirnas if baseline["green_regions"] else []):
        item["reference_sequence_match"] = (
            sequence[item["start_1based"] - 1:item["end_1based_inclusive"]].upper() == item["sequence"].upper())
    sirna_starts = [s["start_1based"] for s in sirnas]
    raw_records = sorted(records, key=lambda row: row["sirna_position"]) if baseline["green_regions"] else []
    raw_starts = [row["sirna_position"] for row in raw_records]
    signature_counts = {}
    for row in (records if baseline["green_regions"] else []):
        signature = json.dumps(row, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
        signature_counts[signature] = signature_counts.get(signature, 0) + 1
    duplicate_extras = {}
    for signature, count in signature_counts.items():
        if count > 1:
            position = json.loads(signature)["sirna_position"]
            duplicate_extras[position] = duplicate_extras.get(position, 0) + count - 1
    size = baseline["sirna_size"]
    evaluations = report.get("offtarget_evaluations", [])
    run = baseline["sifi_run"] or {}
    source = report["sifi"].get("source", {})
    trace_fields = [reference.get("source"), reference.get("version"),
                    source.get("sha256") if isinstance(source, dict) else None,
                    run.get("query_id"), (run.get("parameters") or {}).get("sirna_size")]
    traceability = sum(value is not None and value != "" for value in trace_fields) / len(trace_fields)
    index = off_index = None
    if baseline["green_regions"]:
        index = WindowIndex(sirnas, raw_records, duplicate_extras, METRICS, _gc_percent)
        index.set_reference(sequence)
        off_index = OfftargetIndex(evaluations, len(sequence))
    windows = CompactWindows()
    counts_by_region = {}
    for region in baseline["green_regions"]:
        count = 0
        region_length = region["end_1based_inclusive"] - region["start_1based"] + 1
        if region_length < min_length:
            window = _window_metrics(region, region["start_1based"], region_length, sequence,
                                     raw_records, raw_starts, sirnas, sirna_starts, size,
                                     evaluations, traceability, duplicate_extras, index, off_index)
            window["window_class"] = "short_green_region"
            window["conventional_range_comparable"] = False
            window.update(length_scores(window["efficacy_priority_score"], window["length_nt"], length_penalty_weight))
            windows.append(window)
            count += 1
        for length in range(min_length, max_length + 1):
            last_start = region["end_1based_inclusive"] - length + 1
            for start in range(region["start_1based"], last_start + 1):
                window = _window_metrics(region, start, length, sequence, raw_records, raw_starts,
                                         sirnas, sirna_starts, size, evaluations, traceability,
                                         duplicate_extras, index, off_index)
                is_conventional = 180 <= length <= 250
                window["window_class"] = "conventional_180_250nt" if is_conventional else "configured_range"
                window["conventional_range_comparable"] = is_conventional
                window.update(length_scores(window["efficacy_priority_score"], window["length_nt"], length_penalty_weight))
                windows.append(window)
                count += 1
        counts_by_region[region["region_id"]] = count
    # Keep the documentation-quality audit separate from this exploratory priority rank.
    selected = windows.best
    if selected:
        start, end = selected["start_1based"], selected["end_1based_inclusive"]
        chosen_region = next(r for r in baseline["green_regions"] if r["region_id"] == selected["region_id"])
        selected = {**selected, "sequence": sequence[start - 1:end],
                    "offtarget_matches": _overlapping_offtargets(evaluations, start, end, len(sequence)),
                    "reference_source": reference.get("source"),
                    "reference_version": reference.get("version"),
                    "selection_interpretation": "Prioridad heuristica: 75% eficiencia original si-Fi y 25% GC; diversidad descriptiva y de desempate, sin peso en el puntaje. No es eficacia validada.",
                    "green_region": chosen_region}
    formula = ("efficacy_priority_score = 100 * (0.75 * si-Fi is_efficient True / evaluable + "
               "0.25 * siRNAs con GC 30-52% / siRNAs con GC valido). Diversidad: 0% de peso; se conserva como metrica y desempate. "
               "Se conserva la validacion de disponibilidad de las metricas; no se imputan. GC 30-52%: Reynolds et al. 2004. "
               "Es una heuristica exploratoria, no eficacia validada. quality_score audita la documentacion "
               "por separado y no participa en esta puntuacion.")
    warnings = list(dict.fromkeys(baseline["warnings"] + metric_warnings))
    if not windows:
        warnings.append(f"Ninguna región verde contiene ventanas de {min_length}–{max_length} nt.")
    if selected:
        warnings.append("La ventana priorizada usa una heuristica exploratoria no validada y requiere evaluacion experimental.")
    elif windows:
        warnings.append("No hay ventanas con clasificaciones si-Fi evaluables; no se calculo prioridad ni se selecciono una.")
    for row in windows.short_windows:
        if row["window_class"] == "short_green_region" and row["data_sufficiency_status"] not in (
                "priority_components_available", "partial_original_metrics"):
            warnings.append("La region corta {0} ({1}-{2}) tiene evidencia insuficiente: {3}.".format(
                row["region_id"], row["start_1based"], row["end_1based_inclusive"],
                row["data_sufficiency_status"]))
        elif row["window_class"] == "short_green_region" and row["data_sufficiency_status"] == "partial_original_metrics":
            warnings.append("La region corta {0} ({1}-{2}) tiene metricas originales incompletas: {3}% faltante.".format(
                row["region_id"], row["start_1based"], row["end_1based_inclusive"], row["missing_data_percent"]))
    return {
        "schema_version": "1.0", "candidate_id": baseline["candidate_id"],
        "coordinate_system": "1-based inclusive; source si-Fi sirna_position is 1-based. Python sequence slices are converted internally to 0-based half-open.",
        "window_length_range_nt": [min_length, max_length], "stride_nt": 1,
        "green_regions": baseline["green_regions"], "green_region_window_counts": counts_by_region,
        "green_mask_source": baseline["green_mask_source"],
        "evaluated_window_count": len(windows), "windows": windows, "selected_window": selected,
        "selection_status": ("selected" if selected else
                             "no_windows_available" if not windows else "no_evaluable_efficacy_metrics"),
        "selection_method": {"score_name": "adjusted_priority_score", "formula": formula,
                              "length_penalty_weight": length_penalty_weight,
                              "adjustment_formula": "length_penalty = weight * max(0, 1 - length_nt / 180); adjusted_priority_score = max(0, efficacy_priority_score - length_penalty)",
                              "components": ["sifi_is_efficient_fraction", "sirna_gc_30_52_fraction",
                                             "unique_sirna_sequence_fraction"],
                              "evidence_quality_score_formula": "100 * mean(record_metric_completeness, coordinate_sequence_consistency, position_start_coverage, exact_duplicate_integrity, traceability)",
                              "evidence_quality_score_components": ["record_metric_completeness",
                                                                     "coordinate_sequence_consistency",
                                                                     "position_start_coverage",
                                                                     "exact_duplicate_integrity", "traceability"],
                              "tie_break": "Mayor prioridad; empate: mayor fraccion is_efficient, fraccion GC objetivo, diversidad, calidad de evidencia, inicio menor, longitud menor y region_id.",
                              "interpretation": "Heuristica exploratoria; no equivale a eficacia medida ni validada."},
        "input_provenance": {"sifi_source": source, "sifi_run": run,
                             "reference_source": reference.get("source"),
                             "reference_version": reference.get("version")},
        "warnings": list(dict.fromkeys(warnings)),
        "offtarget_evaluations": [{"organism": e.get("organism"), "status": e.get("status"),
                                   "criterion_status": e.get("criterion", {}).get("status"),
                                   "database": (e.get("result", {}).get("run") or e.get("manifest", {})).get("database"),
                                   "database_fingerprint": (e.get("result", {}).get("run") or e.get("manifest", {})).get("database_fingerprint"),
                                   "result_sha256": (e.get("result", {}).get("run") or {}).get("result_sha256") or
                                                    e.get("result", {}).get("source", {}).get("sha256")}
                                  for e in evaluations],
    }
