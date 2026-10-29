"""Python 2.7 bridge: export exact si-Fi mask positions without rendering a plot.

Run with the Python/Biopython environment used by si-Fi. This file lives outside
the original si-Fi source tree and only calls its read-only helper.
"""
from __future__ import print_function

import argparse
import json
import os
import sys


def runs(positions):
    if not positions:
        return []
    ordered = sorted(positions)
    result = []
    start = previous = ordered[0]
    for position in ordered[1:]:
        if position != previous + 1:
            result.append({"start_1based": start, "end_1based_inclusive": previous,
                           "length_nt": previous - start + 1})
            start = position
        previous = position
    result.append({"start_1based": start, "end_1based_inclusive": previous,
                   "length_nt": previous - start + 1})
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description="Exporta la máscara main_target_dict de si-Fi sin leer la imagen")
    parser.add_argument("--sifi-json", required=True)
    parser.add_argument("--sirna-size", type=int, required=True)
    parser.add_argument("--sifi-root", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    if args.sirna_size < 1:
        parser.error("--sirna-size debe ser positivo")
    sys.path.insert(0, os.path.abspath(args.sifi_root))
    import general_helpers
    unused_off_targets, main_target_dict, unused_efficient, unused_histogram = \
        general_helpers.get_target_data(args.sifi_json, args.sirna_size)
    positions = set()
    for target_positions in main_target_dict.values():
        positions.update(target_positions)
    payload = {
        "schema_version": "1.0",
        "source_function": "general_helpers.get_target_data",
        "coordinate_system": "1-based positions",
        "sirna_size": args.sirna_size,
        "positions_1based": sorted(positions),
        "green_regions": runs(positions),
        "source_json": os.path.abspath(args.sifi_json),
        "main_target_dict_keys": sorted(str(key) for key in main_target_dict.keys()),
    }
    with open(args.out, "w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
    print("Exportadas %d posiciones en %d runs verdes: %s" %
          (len(positions), len(payload["green_regions"]), args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
