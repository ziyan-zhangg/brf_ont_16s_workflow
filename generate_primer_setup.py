#!/usr/bin/env python3
"""
generate_primer_setup.py

Generate a primer setup file for the ONT 16S workflow.

For each sample in the input sample sheet, look up its barcode in the
Twist 384 barcode reference and write out the corresponding indices and
primers in the format required downstream.

Usage:
    python3 generate_primer_setup.py <samplesheet.csv> [-o OUTPUT_DIR] [-b BARCODE_FILE]

Default paths assume execution on Gadi under /g/data/vz35/ONT_16s_workflow.
"""

import argparse
import csv
import re
import sys
from datetime import datetime
from pathlib import Path


# --- Defaults for Gadi ---------------------------------------------------------
DEFAULT_BARCODE_FILE = Path(
    "/g/data/vz35/ONT_16s_workflow/tools/Twist_16S_384_barcode.txt"
)
DEFAULT_OUTPUT_DIR = Path("/g/data/vz35/ONT_16s_workflow/sample_sheet")

OUTPUT_HEADER = ["SampleID", "FwIndex", "FwPrimer", "RvIndex", "RvPrimer"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate primer_setup_<date>.txt from a sample sheet."
    )
    parser.add_argument(
        "samplesheet",
        type=Path,
        help="Input sample sheet CSV (columns: Client, Sample_ID, Barcode).",
    )
    parser.add_argument(
        "-b",
        "--barcode-file",
        type=Path,
        default=DEFAULT_BARCODE_FILE,
        help=f"Twist 384 barcode reference TSV (default: {DEFAULT_BARCODE_FILE}).",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--date",
        default=None,
        help="Date stamp for output filename (default: today, YYYYMMDD).",
    )
    return parser.parse_args()


def load_barcode_reference(barcode_file: Path) -> dict:
    """
    Load the Twist 384 barcode reference into a dict keyed by Barcode ID.

    Expected columns:
        Barcode  Tem_index_F  Tem_index_R  FwIndex  FwPrimer  RvIndex  RvPrimer
    """
    if not barcode_file.exists():
        sys.exit(f"ERROR: Barcode file not found: {barcode_file}")

    barcode_map = {}
    with barcode_file.open("r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"Barcode", "FwIndex", "FwPrimer", "RvIndex", "RvPrimer"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            sys.exit(
                f"ERROR: Barcode file {barcode_file} is missing column(s): "
                f"{', '.join(sorted(missing))}"
            )
        for row in reader:
            bc = row["Barcode"].strip()
            if not bc:
                continue
            barcode_map[bc] = {
                "FwIndex": row["FwIndex"].strip(),
                "FwPrimer": row["FwPrimer"].strip(),
                "RvIndex": row["RvIndex"].strip(),
                "RvPrimer": row["RvPrimer"].strip(),
            }
    return barcode_map


def load_samplesheet(samplesheet: Path) -> list:
    """
    Load the sample sheet CSV. Expected columns: Client, Sample_ID, Barcode.
    Returns a list of (sample_id, barcode) tuples in the original order.
    """
    if not samplesheet.exists():
        sys.exit(f"ERROR: Sample sheet not found: {samplesheet}")

    samples = []
    with samplesheet.open("r", newline="") as f:
        reader = csv.DictReader(f)
        required = {"Sample_ID", "Barcode"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            sys.exit(
                f"ERROR: Sample sheet {samplesheet} is missing column(s): "
                f"{', '.join(sorted(missing))}"
            )
        for i, row in enumerate(reader, start=2):  # start=2 -> account for header
            sample_id = (row.get("Sample_ID") or "").strip()
            barcode = (row.get("Barcode") or "").strip()
            if not sample_id and not barcode:
                continue  # skip blank lines
            if not sample_id or not barcode:
                sys.exit(
                    f"ERROR: Sample sheet line {i} is missing Sample_ID or Barcode."
                )
            sanitized = re.sub(r"[^A-Za-z0-9_-]", "_", sample_id)
            if sanitized != sample_id:
                print(
                    f"WARNING: Sample_ID {sample_id!r} (line {i}) rewritten to {sanitized!r}.",
                    file=sys.stderr,
                )
            samples.append((sanitized, barcode))
    if not samples:
        sys.exit(f"ERROR: No samples found in {samplesheet}.")
    return samples


def main():
    args = parse_args()

    barcode_map = load_barcode_reference(args.barcode_file)
    samples = load_samplesheet(args.samplesheet)

    missing_barcodes = [(sid, bc) for sid, bc in samples if bc not in barcode_map]
    if missing_barcodes:
        print(
            f"ERROR: The following sample(s) had barcodes not found in {args.barcode_file}:",
            file=sys.stderr,
        )
        for sid, bc in missing_barcodes:
            print(f"  {sid}\t{bc}", file=sys.stderr)
        sys.exit(1)

    date_stamp = args.date or datetime.now().strftime("%Y%m%d")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_path = args.output_dir / f"16S_primer_setup_{date_stamp}.txt"

    with output_path.open("w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_HEADER)
        for sample_id, barcode in samples:
            entry = barcode_map[barcode]
            writer.writerow([
                sample_id,
                entry["FwIndex"],
                entry["FwPrimer"],
                entry["RvIndex"],
                entry["RvPrimer"],
            ])

    print(f"Wrote {len(samples)} samples to {output_path}")


if __name__ == "__main__":
    main()
