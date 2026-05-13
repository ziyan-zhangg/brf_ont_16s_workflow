#!/usr/bin/env python3
"""
filter_chopper_demux_minibar.py

Python rewrite of the core function of Filter_Chopper_Demux_Minibar.qsub.

Pipeline:
    Step 1 - Chopper:    filter each fastq.gz for Q >= 15, length 1kb-2kb
    Step 2 - Minibar:    demultiplex each filtered file but keep the primer
             usage: python3 minibar.py -e 1 -E 5 -l 200 -M 2 -F
    Step 2.5 - Cutadapt: two-pass per-sample orientation + 5'/3' primer trim
             Pass 1: -g FWD --revcomp --rename={header} --discard-untrimmed
             Pass 2: -a REV_RC  (tolerant; keeps truncated reads)
             Replaces sample_<SampleID>.fastq with <SampleID>.fastq.gz.
    Step 3 - Organise:   group by client, write summaries
    Final  - Logs:       move all *.log/*.txt logs to run_log_<date>/

This script assumes the primer setup file has already been generated
(via generate_primer_setup.py) before it is invoked. A separate qsub
wrapper is responsible for module loading, PBS directives, and calling
this script.

Usage:
    python3 filter_chopper_demux_minibar.py \\
        --raw-input-dir   /path/to/fastq_pass \\
        --output-dir      /path/to/run_output \\
        --primer-file     /path/to/16S_primer_setup_YYYYMMDD.txt \\
        --samplesheet     /path/to/16s_samplesheet.csv \\
        [--chopper /path/to/chopper] \\
        [--minibar  /path/to/minibar.py] \\
        [--min-quality 15] [--min-length 1000] [--max-length 2000]
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

# Local import: two-pass cutadapt step (file lives in the same tools/ dir).
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cutadapt_2pass import run_cutadapt_step, DEFAULT_CUTADAPT  # noqa: E402


# --- Defaults (Gadi) ---------------------------------------------------------
DEFAULT_CHOPPER = Path("/g/data/vz35/zpfeng/tools/chopper/chopper-linux-musl")
DEFAULT_MINIBAR = Path("/g/data/vz35/ONT_16s_workflow/tools/minibar/minibar.py")


# =============================================================================
# Helpers
# =============================================================================
def log(msg: str = "") -> None:
    """Flushed print so the PBS job log stays in order."""
    print(msg, flush=True)


def banner(title: str) -> None:
    log("========================================")
    log(f" {title}: {datetime.now():%Y-%m-%d %H:%M:%S}")
    log("========================================")


def count_fastq_reads(path: Path) -> int:
    """Count reads (4 lines per record) in a plain or gzipped fastq."""
    opener = gzip.open if path.suffix == ".gz" else open
    lines = 0
    with opener(path, "rb") as fh:
        for _ in fh:
            lines += 1
    return lines // 4


def sanitise(name: str) -> str:
    """Match the awk gsub(/[^A-Za-z0-9._-]/, "_") used in the bash version."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def sanitise_sample_id(sid: str) -> str:
    """Match the rule applied by generate_primer_setup.py — must agree with it
    so we look for the same filename minibar wrote based on the primer setup."""
    return re.sub(r"[^A-Za-z0-9_-]", "_", sid)


# =============================================================================
# Step 1 - Chopper filtering
# =============================================================================
def run_chopper(
    raw_input_dir: Path,
    filtered_dir: Path,
    chopper_prog: Path,
    min_quality: int,
    min_length: int,
    max_length: int,
) -> list[Path]:
    """
    Pipe each fastq.gz through chopper and gzip the result back.

    Equivalent to:
        zcat IN | chopper -q Q --minlength MIN --maxlength MAX | gzip > OUT
    """
    banner("STEP 1: Chopper filtering")
    log(f" Input dir:    {raw_input_dir}")
    log(f" Filtered dir: {filtered_dir}")
    log("")

    filtered_dir.mkdir(parents=True, exist_ok=True)
    inputs = sorted(raw_input_dir.glob("*.fastq.gz"))
    if not inputs:
        sys.exit(f"ERROR: no *.fastq.gz files found in {raw_input_dir}")

    outputs: list[Path] = []
    for f in inputs:
        sample = f.name.removesuffix(".fastq.gz")
        out_file = filtered_dir / f"{sample}_filtered.fastq.gz"

        log(f"  Filtering: {sample}")
        log(f"    Input:  {f}")
        log(f"    Output: {out_file}")
        log(f"    Start:  {datetime.now():%Y-%m-%d %H:%M:%S}")

        # zcat | chopper | gzip   -- assemble as a process pipeline
        with gzip.open(f, "rb") as zin, gzip.open(out_file, "wb") as zout:
            chopper = subprocess.Popen(
                [
                    str(chopper_prog),
                    "-q", str(min_quality),
                    "--minlength", str(min_length),
                    "--maxlength", str(max_length),
                ],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
            )

            # Feed decompressed bytes into chopper, stream chopper stdout to gzip out.
            # Use small chunks so we don't buffer entire fastq files in memory.
            assert chopper.stdin is not None and chopper.stdout is not None
            try:
                # writer thread substitute: alternate reads from input and writes to chopper
                # Simpler: use shutil.copyfileobj on each side via a thread.
                import threading

                def feed():
                    try:
                        shutil.copyfileobj(zin, chopper.stdin, length=1 << 20)
                    finally:
                        chopper.stdin.close()

                t = threading.Thread(target=feed, daemon=True)
                t.start()
                shutil.copyfileobj(chopper.stdout, zout, length=1 << 20)
                t.join()
            finally:
                rc = chopper.wait()
            if rc != 0:
                sys.exit(f"ERROR: chopper failed on {f} (exit {rc})")

        log(f"    Done:   {datetime.now():%Y-%m-%d %H:%M:%S}")
        log("")
        outputs.append(out_file)

    log(f" Chopper filtering complete: {datetime.now():%Y-%m-%d %H:%M:%S}")
    log("========================================")
    return outputs


def count_filtered_total(filtered_files: Iterable[Path]) -> int:
    log("")
    log(f" Counting total filtered input reads: {datetime.now():%Y-%m-%d %H:%M:%S}")
    total = 0
    for f in filtered_files:
        n = count_fastq_reads(f)
        log(f"  {f.name}: {n} reads")
        total += n
    log(f" Total filtered input reads: {total}")
    log("========================================")
    return total


# =============================================================================
# Step 2 - Minibar demultiplexing
# =============================================================================
def run_minibar(
    filtered_files: list[Path],
    output_dir: Path,
    primer_file: Path,
    minibar_prog: Path,
) -> list[Path]:
    """
    Run minibar on each filtered file in its own subdirectory.

    Equivalent to:
        cd OUTDIR/<base> && python3 minibar.py -e 1 -E 5 -l 200 -M 2 -T -F PRIMER IN
    """
    log("")
    banner("STEP 2: Minibar demultiplexing")
    log(f" Files found:")
    for f in filtered_files:
        log(f"   {f}")
    log("")

    per_file_dirs: list[Path] = []
    for input_file in filtered_files:
        base = input_file.name.removesuffix(".fastq.gz")
        outd = output_dir / base
        outd.mkdir(parents=True, exist_ok=True)

        log(f"--- Processing: {base} ---")
        log(f"  Input:  {input_file}")
        log(f"  Outdir: {outd}")
        log(f"  Start:  {datetime.now():%Y-%m-%d %H:%M:%S}")

        cmd = [
            sys.executable, str(minibar_prog),
            "-e", "1",
            "-E", "5",
            "-l", "200",
            "-M", "2",
            "-F",
            str(primer_file),
            str(input_file),
        ]
        # minibar writes per-sample fastqs into its current working directory
        rc = subprocess.call(cmd, cwd=outd)
        if rc != 0:
            sys.exit(f"ERROR: minibar failed on {input_file} (exit {rc})")

        log(f"  Done:   {datetime.now():%Y-%m-%d %H:%M:%S}")
        log("")
        per_file_dirs.append(outd)

    log(f" All files demultiplexed: {datetime.now():%Y-%m-%d %H:%M:%S}")
    log("========================================")
    return per_file_dirs


def merge_per_file_outputs(
    per_file_dirs: list[Path],
    integrated_dir: Path,
) -> None:
    """Concatenate matching sample_*.fastq across all per-file subdirs."""
    log("")
    banner("Merging outputs")
    integrated_dir.mkdir(parents=True, exist_ok=True)

    # collect every distinct sample_*.fastq filename
    sample_names: set[str] = set()
    for d in per_file_dirs:
        for f in d.glob("sample_*.fastq"):
            sample_names.add(f.name)

    for name in sorted(sample_names):
        log(f"  Merging: {name}")
        dest = integrated_dir / name
        with dest.open("ab") as out:
            for d in per_file_dirs:
                src = d / name
                if src.is_file():
                    with src.open("rb") as inp:
                        shutil.copyfileobj(inp, out, length=1 << 20)

    log(f" Merging complete: {datetime.now():%Y-%m-%d %H:%M:%S}")
    log("========================================")


def cleanup(per_file_dirs: list[Path], filtered_dir: Path) -> None:
    log("")
    banner("Cleaning up per-file subdirectories")
    for d in per_file_dirs:
        log(f"  Removing: {d}/")
        shutil.rmtree(d, ignore_errors=True)
    log(f" Cleanup complete: {datetime.now():%Y-%m-%d %H:%M:%S}")
    log("========================================")

    log("")
    banner("Deleting Chopper filtered files")
    shutil.rmtree(filtered_dir, ignore_errors=True)
    log(f" Deleted: {filtered_dir}")
    log(f" Cleanup complete: {datetime.now():%Y-%m-%d %H:%M:%S}")
    log("========================================")


# =============================================================================
# Step 3 - Organise by client + summaries
# =============================================================================
@dataclass
class ClientSample:
    client: str
    sample_id: str


def load_client_map(samplesheet: Path) -> list[ClientSample]:
    """Read Client + Sample_ID columns; sanitise client names like the awk does."""
    if not samplesheet.exists():
        sys.exit(f"ERROR: sample sheet not found: {samplesheet}")
    out: list[ClientSample] = []
    with samplesheet.open("r", newline="") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames or "Client" not in reader.fieldnames \
                or "Sample_ID" not in reader.fieldnames:
            sys.exit(
                f"ERROR: {samplesheet} must have Client and Sample_ID columns"
            )
        for row in reader:
            client = (row.get("Client") or "").strip()
            sid = (row.get("Sample_ID") or "").strip()
            if not client or not sid:
                continue
            out.append(ClientSample(sanitise(client), sanitise_sample_id(sid)))
    return out


def organise_by_client(
    integrated_dir: Path,
    samples: list[ClientSample],
) -> None:
    log("")
    banner("STEP 3: Organising reads by client")

    for cs in samples:
        client_dir = integrated_dir / cs.client
        client_dir.mkdir(parents=True, exist_ok=True)
        # Cutadapt step renamed sample_<id>.fastq -> <id>.fastq.gz.
        # Fall back to the pre-cutadapt name if cutadapt skipped the sample
        # (e.g. no primers in the setup file).
        src = integrated_dir / f"{cs.sample_id}.fastq.gz"
        legacy = integrated_dir / f"sample_{cs.sample_id}.fastq"
        if src.is_file():
            shutil.move(str(src), str(client_dir / src.name))
            log(f"  Moved: {src.name} -> {cs.client}/")
        elif legacy.is_file():
            shutil.move(str(legacy), str(client_dir / legacy.name))
            log(f"  Moved (uncut): {legacy.name} -> {cs.client}/")
        else:
            log(f"  WARNING: neither {src.name} nor {legacy.name} found")

    log("")
    log(" Generating per-client summaries...")
    for client_subdir in sorted(p for p in integrated_dir.iterdir() if p.is_dir()):
        fastqs = sorted(client_subdir.glob("*.fastq.gz")) \
                 + sorted(client_subdir.glob("sample_*.fastq"))
        client_total = sum(count_fastq_reads(f) for f in fastqs)

        lines = [
            "========================================",
            f" Client: {client_subdir.name}",
            f" Date:   {datetime.now():%Y-%m-%d %H:%M:%S}",
            "========================================",
            f"{'Sample':<40} {'Reads':>10}",
            "----------------------------------------",
        ]
        for f in fastqs:
            lines.append(f"{f.stem:<40} {count_fastq_reads(f):>10d}")
        lines += [
            "----------------------------------------",
            f"{'CLIENT TOTAL':<40} {client_total:>10d}",
            "========================================",
        ]
        body = "\n".join(lines)
        log(body)
        summary = client_subdir / "summary.txt"
        summary.write_text(body + "\n")
        log(f"  Saved: {summary}")
        log("")

    log(f" Organisation complete: {datetime.now():%Y-%m-%d %H:%M:%S}")
    log("========================================")


def write_run_summary(
    integrated_dir: Path,
    summary_file: Path,
    raw_input_dir: Path,
    filtered_dir: Path,
    output_dir: Path,
    total_input: int,
) -> None:
    log("")
    banner("Generating read count summary")

    # collect all per-sample fastqs (top-level + one client level deep).
    # Cutadapt step renames sample_<id>.fastq -> <id>.fastq.gz; the minibar
    # catch-all bins (sample_unk, sample_Multiple_Matches) keep the old name.
    fastqs: list[Path] = []
    fastqs += sorted(integrated_dir.glob("sample_*.fastq"))
    fastqs += sorted(integrated_dir.glob("*/sample_*.fastq"))
    fastqs += sorted(integrated_dir.glob("*.fastq.gz"))
    fastqs += sorted(integrated_dir.glob("*/*.fastq.gz"))

    lines = [
        "========================================",
        " Minibar Demultiplexing Read Count Summary",
        f" Raw input dir:     {raw_input_dir}",
        f" Filtered dir:      {filtered_dir}",
        f" Out dir:           {output_dir}",
        f" Date:              {datetime.now():%Y-%m-%d %H:%M:%S}",
        "========================================",
        f"{'Sample':<40} {'Reads':>10} {'% of Input':>10}",
        "----------------------------------------",
    ]

    total_demux = 0
    multi_match_reads = 0
    unk_reads = 0
    for f in fastqs:
        # .fastq.gz files have suffixes ".fastq" + ".gz"; .stem strips only ".gz".
        sample = f.name.removesuffix(".fastq.gz") if f.name.endswith(".fastq.gz") else f.stem
        reads = count_fastq_reads(f)
        pct = (reads / total_input * 100) if total_input else 0.0
        lines.append(f"{sample:<40} {reads:>10d} {pct:>9.2f}%")
        total_demux += reads
        if sample == "sample_Multiple_Matches":
            multi_match_reads = reads
        elif sample == "sample_unk":
            unk_reads = reads

    success_reads = total_input - multi_match_reads - unk_reads
    success_pct = (success_reads / total_input * 100) if total_input else 0.0
    total_pct = (total_demux / total_input * 100) if total_input else 0.0

    lines += [
        "----------------------------------------",
        f"{'TOTAL DEMULTIPLEXED':<40} {total_demux:>10d} {total_pct:>9.2f}%",
        f"{'SUCCESSFULLY DEMULTIPLEXED':<40} {success_reads:>10d} {success_pct:>9.2f}%",
        f"{'TOTAL FILTERED INPUT':<40} {total_input:>10d} {100.00:>9.2f}%",
        "========================================",
    ]

    body = "\n".join(lines)
    log(body)
    summary_file.write_text(body + "\n")
    log("")
    log(f" Summary saved to: {summary_file}")
    log(f" Pipeline complete: {datetime.now():%Y-%m-%d %H:%M:%S}")
    log("========================================")


# =============================================================================
# CLI
# =============================================================================
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--raw-input-dir", type=Path, required=True,
                   help="Directory containing raw *.fastq.gz (e.g. fastq_pass).")
    p.add_argument("--output-dir", type=Path, required=True,
                   help="Run output directory (will be created).")
    p.add_argument("--primer-file", type=Path, required=True,
                   help="Primer setup file produced by generate_primer_setup.py.")
    p.add_argument("--samplesheet", type=Path, required=True,
                   help="Sample sheet CSV (Client, Sample_ID, Barcode).")
    p.add_argument("--chopper", type=Path, default=DEFAULT_CHOPPER,
                   help=f"Path to chopper binary (default: {DEFAULT_CHOPPER}).")
    p.add_argument("--minibar", type=Path, default=DEFAULT_MINIBAR,
                   help=f"Path to minibar.py (default: {DEFAULT_MINIBAR}).")
    p.add_argument("--min-quality", type=int, default=15)
    p.add_argument("--min-length", type=int, default=1000)
    p.add_argument("--max-length", type=int, default=2000)
    p.add_argument("--cutadapt", type=Path, default=DEFAULT_CUTADAPT,
                   help=f"Path to cutadapt binary (default: {DEFAULT_CUTADAPT}).")
    p.add_argument("--cutadapt-threads", type=int, default=4,
                   help="Threads for cutadapt (default: 4).")
    p.add_argument("--cutadapt-error-rate", type=float, default=0.2,
                   help="Cutadapt error rate -e (default: 0.2).")
    return p.parse_args()


def collect_logs(output_dir: Path) -> Path:
    """Move all *.log and *.txt summaries into a dated run_log_<date> directory."""
    log("")
    banner("Collecting logs")
    date_stamp = datetime.now().strftime("%Y%m%d")
    log_dir = output_dir / f"run_log_{date_stamp}"
    log_dir.mkdir(parents=True, exist_ok=True)

    moved = 0
    # Top-level log/summary files in output_dir
    for pattern in ("*.log", "*.txt"):
        for src in sorted(output_dir.glob(pattern)):
            if src.is_file() and src.parent == output_dir:
                dest = log_dir / src.name
                shutil.move(str(src), str(dest))
                moved += 1
                log(f"  Moved: {src.name}")

    # Cutadapt per-sample logs live in output_dir/cutadapt_logs/
    cutadapt_logs = output_dir / "cutadapt_logs"
    if cutadapt_logs.is_dir():
        dest = log_dir / cutadapt_logs.name
        # If a previous run already populated run_log_<date>/cutadapt_logs/, merge
        if dest.exists():
            for f in cutadapt_logs.iterdir():
                shutil.move(str(f), str(dest / f.name))
            cutadapt_logs.rmdir()
        else:
            shutil.move(str(cutadapt_logs), str(dest))
        log(f"  Moved: cutadapt_logs/ -> {log_dir.name}/cutadapt_logs/")
        moved += 1

    log(f" Collected {moved} log items into {log_dir}")
    log("========================================")
    return log_dir


def main() -> None:
    args = parse_args()

    if not args.raw_input_dir.is_dir():
        sys.exit(f"ERROR: raw input dir not found: {args.raw_input_dir}")
    if not args.primer_file.is_file():
        sys.exit(f"ERROR: primer file not found: {args.primer_file}")
    if not args.chopper.is_file():
        sys.exit(f"ERROR: chopper not found: {args.chopper}")
    if not args.minibar.is_file():
        sys.exit(f"ERROR: minibar.py not found: {args.minibar}")
    if not args.cutadapt.is_file():
        sys.exit(f"ERROR: cutadapt not found: {args.cutadapt}")

    output_dir = args.output_dir
    filtered_dir = output_dir / "chopper_filtered"
    integrated_dir = output_dir / "integrated_demultiplexing"
    summary_file = output_dir / "read_counts_summary.txt"

    output_dir.mkdir(parents=True, exist_ok=True)
    integrated_dir.mkdir(parents=True, exist_ok=True)

    log(f" Primer file: {args.primer_file}")

    # ---- Step 1 ----
    filtered = run_chopper(
        args.raw_input_dir, filtered_dir, args.chopper,
        args.min_quality, args.min_length, args.max_length,
    )
    total_input = count_filtered_total(filtered)

    # ---- Step 2 ----
    per_file_dirs = run_minibar(filtered, output_dir, args.primer_file, args.minibar)
    merge_per_file_outputs(per_file_dirs, integrated_dir)
    cleanup(per_file_dirs, filtered_dir)

    # ---- Step 2.5: two-pass cutadapt ----
    run_cutadapt_step(
        integrated_dir=integrated_dir,
        primer_file=args.primer_file,
        output_dir=output_dir,
        cutadapt=args.cutadapt,
        threads=args.cutadapt_threads,
        error_rate=args.cutadapt_error_rate,
    )

    # ---- Step 3 ----
    samples = load_client_map(args.samplesheet)
    organise_by_client(integrated_dir, samples)
    write_run_summary(
        integrated_dir, summary_file,
        args.raw_input_dir, filtered_dir, output_dir,
        total_input,
    )

    # ---- Final: gather logs into run_log_<date>/ ----
    collect_logs(output_dir)


if __name__ == "__main__":
    main()
