"""Resumen descriptivo sin recalcular ni optimizar actividad biológica."""
from collections import defaultdict


def summarize(records: list[dict]) -> dict:
    grouped = defaultdict(list)
    for row in records:
        grouped[(row.get("query_name"), row["sirna_name"], row["sirna_position"], row["sirna_sequence"])].append(row)
    points, warnings = [], []
    for key, rows in sorted(grouped.items(), key=lambda item: item[0][2]):
        values = {r.get("accessibility_value") for r in rows}
        labels = {r.get("is_efficient") for r in rows}
        conflict = len(values) > 1 or len(labels) > 1
        if conflict:
            warnings.append(f"Métricas discordantes entre hits de {key[1]}; no se promedian")
        points.append({"sirna_name": key[1], "start_1based": key[2],
                       "end_1based_inclusive": key[2]+len(key[3])-1,
                       "accessibility_value": next(iter(values)) if not conflict else None,
                       "is_efficient": next(iter(labels)) if not conflict else None,
                       "hit_rows": len(rows), "conflict": conflict})
    return {"hit_rows": len(records), "unique_sirnas": len(points),
            "sifi_true": sum(p["is_efficient"] is True for p in points),
            "sifi_false": sum(p["is_efficient"] is False for p in points),
            "unknown_or_conflicting": sum(p["is_efficient"] is None for p in points),
            "points": points, "warnings": warnings}
