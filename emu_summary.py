#!/usr/bin/env python3
"""
emu_summary.py — Summarise Emu abundance output with optional comparison
to a reference mock community.

Usage:
    python3 emu_summary.py --result <path/to/*_rel-abundance.tsv>
    python3 emu_summary.py --result <...> --reference references/zymobiomics_d6305.json

The Emu TSV has these columns (header row present):
    tax_id  abundance  ...  species  ...  estimated counts
Special rows: tax_id == "unmapped" and tax_id == "mapped_unclassified".

Reference JSON format:
    {
        "name": "ZymoBIOMICS D6305",
        "species": {
            "Bacillus subtilis": 17.4,
            ...
        },
        "aliases": {
            "Lactobacillus fermentum": "Limosilactobacillus fermentum"
        }
    }
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--result", required=True, type=Path,
                   help="Path to Emu *_rel-abundance.tsv")
    p.add_argument("--reference", type=Path, default=None,
                   help="Path to reference mock-community JSON (optional)")
    p.add_argument("--top-n", type=int, default=20,
                   help="When no reference is given, show top N species (default: 20)")
    return p.parse_args()


def load_emu_tsv(path: Path) -> tuple[dict[str, dict], float, float]:
    """
    Returns:
        species_rows: dict mapping species name -> {"abundance_pct": float, "counts": float}
        unmapped_counts: estimated counts for tax_id == "unmapped"
        unclassified_counts: estimated counts for tax_id == "mapped_unclassified"
    """
    species_rows: dict[str, dict] = {}
    unmapped = 0.0
    unclassified = 0.0

    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        # Emu's columns vary slightly between versions; locate by header name.
        if "abundance" not in reader.fieldnames or "species" not in reader.fieldnames:
            sys.exit(f"ERROR: unexpected Emu header in {path}: {reader.fieldnames}")
        count_col = "estimated counts"
        if count_col not in reader.fieldnames:
            # Fall back to last column if name differs
            count_col = reader.fieldnames[-1]

        for row in reader:
            tax_id = (row.get("tax_id") or "").strip()
            try:
                counts = float(row.get(count_col) or 0)
            except ValueError:
                counts = 0.0

            if tax_id == "unmapped":
                unmapped = counts
                continue
            if tax_id == "mapped_unclassified":
                unclassified = counts
                continue

            species = (row.get("species") or "").strip()
            if not species:
                continue
            try:
                abund_pct = float(row["abundance"]) * 100.0
            except (ValueError, KeyError):
                continue
            species_rows[species] = {"abundance_pct": abund_pct, "counts": counts}

    return species_rows, unmapped, unclassified


def pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    if n == 0:
        return 0.0
    sx = sum(xs); sy = sum(ys)
    sxx = sum(x * x for x in xs); syy = sum(y * y for y in ys)
    sxy = sum(x * y for x, y in zip(xs, ys))
    denom = math.sqrt((n * sxx - sx * sx) * (n * syy - sy * sy))
    return (n * sxy - sx * sy) / denom if denom > 0 else 0.0


def report_with_reference(
    species_rows: dict[str, dict],
    unmapped: float,
    unclassified: float,
    reference: dict,
) -> None:
    ref_name = reference.get("name", "reference")
    expected: dict[str, float] = reference["species"]
    aliases: dict[str, str] = reference.get("aliases", {})

    # Resolve aliases — if the observed table contains an alias, treat it as the canonical name.
    resolved: dict[str, dict] = {}
    for sp, vals in species_rows.items():
        canonical = aliases.get(sp, sp)
        if canonical in resolved:
            # Combine if both forms appear (rare but safe)
            resolved[canonical]["abundance_pct"] += vals["abundance_pct"]
            resolved[canonical]["counts"] += vals["counts"]
        else:
            resolved[canonical] = dict(vals)

    # Expected species table
    bar = "=" * 60
    print(bar)
    print(f" EXPECTED SPECIES ({ref_name})")
    print(bar)
    print(f" {'Species':<32} {'Theor.%':>10} {'Obs.%':>10} {'Delta_pp':>10}")
    print(f" {'-' * 32:<32} {'-' * 7:>10} {'-' * 6:>10} {'-' * 8:>10}")

    theor_vals: list[float] = []
    obs_vals: list[float] = []
    for sp, theor in expected.items():
        obs = resolved.get(sp, {}).get("abundance_pct", 0.0)
        delta = obs - theor
        print(f" {sp:<32} {theor:>10.2f} {obs:>10.2f} {delta:>+10.2f}")
        theor_vals.append(theor)
        obs_vals.append(obs)

    r = pearson(theor_vals, obs_vals)
    print()
    print(f" Pearson r (observed vs theoretical):  {r:.4f}")
    print()

    # Unexpected species
    print(bar)
    print(" UNEXPECTED SPECIES (contaminants / classifier noise)")
    print(bar)
    print(f" {'Species':<45} {'Obs.%':>10} {'Est. counts':>12}")
    print(f" {'-' * 45:<45} {'-' * 6:>10} {'-' * 11:>12}")

    unexpected = [(sp, v) for sp, v in resolved.items() if sp not in expected]
    unexpected.sort(key=lambda kv: kv[1]["abundance_pct"], reverse=True)

    if not unexpected:
        print(" (none above min-abundance threshold)")
    else:
        for sp, v in unexpected:
            print(f" {sp:<45} {v['abundance_pct']:>10.4f} {v['counts']:>12.1f}")
    print()

    # Read fate
    print(bar)
    print(" READ FATE")
    print(bar)
    print(f" {'Category':<25} {'Est. counts':>12}")
    print(f" {'-' * 25:<25} {'-' * 11:>12}")
    print(f" {'unmapped':<25} {unmapped:>12.1f}")
    print(f" {'mapped_unclassified':<25} {unclassified:>12.1f}")
    print()


def report_without_reference(
    species_rows: dict[str, dict],
    unmapped: float,
    unclassified: float,
    top_n: int,
) -> None:
    bar = "=" * 60
    items = sorted(species_rows.items(), key=lambda kv: kv[1]["abundance_pct"], reverse=True)

    print(bar)
    print(f" TOP {min(top_n, len(items))} SPECIES BY ABUNDANCE")
    print(bar)
    print(f" {'Species':<45} {'Obs.%':>10} {'Est. counts':>12}")
    print(f" {'-' * 45:<45} {'-' * 6:>10} {'-' * 11:>12}")
    for sp, v in items[:top_n]:
        print(f" {sp:<45} {v['abundance_pct']:>10.4f} {v['counts']:>12.1f}")
    print()

    print(bar)
    print(" READ FATE")
    print(bar)
    print(f" {'Category':<25} {'Est. counts':>12}")
    print(f" {'-' * 25:<25} {'-' * 11:>12}")
    print(f" {'unmapped':<25} {unmapped:>12.1f}")
    print(f" {'mapped_unclassified':<25} {unclassified:>12.1f}")
    print()


def main() -> int:
    args = parse_args()

    if not args.result.is_file():
        print(f"ERROR: result file not found: {args.result}", file=sys.stderr)
        return 1

    species_rows, unmapped, unclassified = load_emu_tsv(args.result)

    if args.reference is None:
        report_without_reference(species_rows, unmapped, unclassified, args.top_n)
        return 0

    if not args.reference.is_file():
        print(f"ERROR: reference file not found: {args.reference}", file=sys.stderr)
        return 1

    with args.reference.open() as fh:
        reference = json.load(fh)

    report_with_reference(species_rows, unmapped, unclassified, reference)
    return 0


if __name__ == "__main__":
    sys.exit(main())
