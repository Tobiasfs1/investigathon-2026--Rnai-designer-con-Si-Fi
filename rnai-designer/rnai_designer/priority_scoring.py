"""Presentation priority; original scientific metrics are never changed."""
import math


def length_scores(score, length_nt, weight=30.0):
    if not math.isfinite(weight) or weight < 0:
        raise ValueError("length penalty weight must be finite and nonnegative")
    penalty = weight * max(0, 1 - length_nt / 180)
    return {"length_penalty": penalty,
            "adjusted_priority_score": None if score is None else max(0, score - penalty)}


def priority_key(row):
    score = row.get("adjusted_priority_score", row.get("efficacy_priority_score"))
    c = row.get("efficacy_priority_components", {})
    return (-(score if score is not None else -1),
            *(-float(c.get(k) or 0) for k in ("sifi_is_efficient_fraction", "sirna_gc_30_52_fraction", "unique_sirna_sequence_fraction")),
            -float(row.get("quality_score") or 0), row.get("start_1based", 0),
            row.get("length_nt", 0), row.get("region_id", ""))
