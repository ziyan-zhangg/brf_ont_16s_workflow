#!/usr/bin/env python3
"""
filter_chopper_demux_minibar.py

Python rewrite of the core function of Filter_Chopper_Demux_Minibar.qsub.

Two filter orderings are supported via --filter-stage; both are runnable on
the same input so they can be compared directly.

The step numbering below is shared by both modes -- step N always means the
same thing -- so the job logs of a pre run and a post run line up. Each mode
skips exactly one step: pre skips step 4, post skips step 1.

--filter-stage pre (default, unchanged legacy behaviour):
    Step 1 - Chopper:  filter each raw fastq.gz for Q/length before demux
    Step 2 - Minibar:  demultiplex each filtered file, merge, keep primers
    Step 3 - Demux QC: per-sample stats (see note below)
    Step 4 - skipped.  Chopper already ran at step 1.
    Step 5 - Cutadapt: two-pass orientation + 5'/3' primer trim
    Step 6 - Organise: group by client, write summaries
    Step 7 - Logs:     move all *.log/*.txt logs to run_log_<date>/

--filter-stage post:
    Step 1 - skipped.  Raw fastq_pass/*.fastq.gz go straight to Minibar;
             fastq_pass is already MinKNOW min-qscore filtered, so there
             is no pre-filter at all.
    Step 2 - Minibar:  demultiplex the raw files, merge, keep primers
    Step 3 - Demux QC: per-sample length/quality profile of the raw
             demultiplexed reads, before any filtering
    Step 4 - Chopper:  filter each merged sample_<SampleID>.fastq
    Step 5 - Cutadapt: as above
    Step 6 - Organise: group by client, write summaries
    Step 7 - Logs:     move all *.log/*.txt logs to run_log_<date>/

Minibar flags are identical in both modes: -e 1 -E 5 -l 200 -M 2 -F.

Demux QC runs in both modes but means different things: in post mode it
profiles genuinely raw demultiplexed reads, in pre mode the reads have
already been through Chopper. The stage is recorded in every QC output.
Only post mode can report a raw_demuxed count, so pre mode leaves that
column of read_counts_summary.txt blank rather than filling it with a
post-filter number.

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
        [--filter-stage pre|post] \\
        [--chopper /path/to/chopper] \\
        [--minibar  /path/to/minibar.py] \\
        [--min-quality 15] [--minlength 1000] [--maxlength 2000] \\
        [--minibar-work-dir $PBS_JOBFS] \\
        [--diagnostic]

    --diagnostic is shorthand for
        --filter-stage post --minlength 100 --maxlength 5000
    so a first pass observes the real length distribution instead of
    assuming the 1-2 kb window.
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
import threading
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import BinaryIO, Iterable

# Local imports: sibling modules in the same tools/ dir.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cutadapt_2pass import (  # noqa: E402
    CutadaptStats,
    DEFAULT_CUTADAPT,
    run_cutadapt_step,
)
from demux_qc import SampleQC, run_demux_qc  # noqa: E402


# --- Defaults (Gadi) ---------------------------------------------------------
DEFAULT_CHOPPER = Path("/g/data/vz35/ONT_16s_workflow/tools/chopper/chopper-linux-musl")
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
# Chopper filtering
#
# Used at two different points depending on --filter-stage:
#   pre  -> Step 1, whole raw fastq.gz files, before demultiplexing
#   post -> Step 4, per-sample merged fastq, after demultiplexing and QC
# Both go through _chopper_stream so the filtering itself is identical.
# =============================================================================
def _chopper_stream(
    chopper_prog: Path,
    src: BinaryIO,
    dst: BinaryIO,
    min_quality: int,
    min_length: int,
    max_length: int,
    what: str,
) -> None:
    """
    Stream src -> chopper -> dst.

    Equivalent to:  cat SRC | chopper -q Q --minlength MIN --maxlength MAX > DST
    The caller owns src/dst, so it decides whether either side is gzipped.
    """
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
    assert chopper.stdin is not None and chopper.stdout is not None
    try:
        # Feed chopper from a thread while the main thread drains its stdout,
        # so neither pipe can fill up and deadlock. Small chunks keep whole
        # fastq files out of memory.
        def feed() -> None:
            try:
                shutil.copyfileobj(src, chopper.stdin, length=1 << 20)
            finally:
                chopper.stdin.close()

        t = threading.Thread(target=feed, daemon=True)
        t.start()
        shutil.copyfileobj(chopper.stdout, dst, length=1 << 20)
        t.join()
    finally:
        rc = chopper.wait()
    if rc != 0:
        sys.exit(f"ERROR: chopper failed on {what} (exit {rc})")


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
    Only used by --filter-stage pre.

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
            _chopper_stream(
                chopper_prog, zin, zout,
                min_quality, min_length, max_length,
                what=str(f),
            )

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


def count_raw_total(raw_files: Iterable[Path]) -> int:
    """Total reads handed to Minibar in post mode (no pre-filter applied)."""
    log("")
    log(f" Counting total raw input reads: {datetime.now():%Y-%m-%d %H:%M:%S}")
    total = 0
    for f in raw_files:
        n = count_fastq_reads(f)
        log(f"  {f.name}: {n} reads")
        total += n
    log(f" Total raw input reads: {total}")
    log("========================================")
    return total


def run_chopper_per_sample(
    integrated_dir: Path,
    chopper_prog: Path,
    min_quality: int,
    min_length: int,
    max_length: int,
    qc_by_sid: dict[str, SampleQC] | None = None,
) -> dict[str, tuple[int, int]]:
    """
    Post-demux Chopper (--filter-stage post, Step 4).

    Filters each merged sample_<SampleID>.fastq in place. Returns
    {sample_id: (reads_before, reads_after)} so the demux loss and the
    filter loss stay separable in the run summary.

    Minibar's catch-all bins (sample_unk, sample_Multiple_Matches) are left
    unfiltered on purpose: their whole diagnostic value is the unfiltered
    length profile of what failed to demultiplex.

    Pre-filter counts are taken from the demux QC pass when it has already
    read the file, so no extra pass over the data is needed.
    """
    log("")
    banner("STEP 4: Chopper filtering (per sample, post-demux)")
    log(f" Integrated dir: {integrated_dir}")
    log(f" Filter:         -q {min_quality} --minlength {min_length} "
        f"--maxlength {max_length}")
    log("")

    counts: dict[str, tuple[int, int]] = {}
    sample_fastqs = sorted(integrated_dir.glob("sample_*.fastq"))
    if not sample_fastqs:
        log("  WARNING: no sample_*.fastq found; nothing to filter.")
        log("========================================")
        return counts

    for fq in sample_fastqs:
        sid = fq.stem.removeprefix("sample_")
        if sid in ("unk", "Multiple_Matches"):
            log(f"  Skipping catch-all bin: {fq.name} (left unfiltered)")
            continue

        qc = (qc_by_sid or {}).get(sid)
        before = qc.reads if qc is not None else count_fastq_reads(fq)

        # Deliberately not *.fastq: a crash mid-filter must not leave behind
        # something the sample_*.fastq globs downstream would pick up.
        tmp_out = fq.with_name(fq.name + ".chopped.tmp")
        log(f"  Filtering: {sid}")
        log(f"    Start:  {datetime.now():%Y-%m-%d %H:%M:%S}")
        with fq.open("rb") as fin, tmp_out.open("wb") as fout:
            _chopper_stream(
                chopper_prog, fin, fout,
                min_quality, min_length, max_length,
                what=str(fq),
            )
        tmp_out.replace(fq)

        after = count_fastq_reads(fq)
        counts[sid] = (before, after)
        pct = (100.0 * after / before) if before else 0.0
        log(f"    Done:   {before:,} -> {after:,} reads ({pct:.2f}% retained)")

    tot_before = sum(b for b, _ in counts.values())
    tot_after = sum(a for _, a in counts.values())
    tot_pct = (100.0 * tot_after / tot_before) if tot_before else 0.0
    log("")
    log(f" Per-sample Chopper complete: {tot_before:,} -> {tot_after:,} reads "
        f"({tot_pct:.2f}% retained) at {datetime.now():%Y-%m-%d %H:%M:%S}")
    log("========================================")
    return counts


# =============================================================================
# Step 2 - Minibar demultiplexing + merge
# =============================================================================
def run_minibar(
    filtered_files: list[Path],
    work_dir: Path,
    primer_file: Path,
    minibar_prog: Path,
) -> list[Path]:
    """
    Run minibar on each input file in its own subdirectory of work_dir.

    Equivalent to:
        cd WORKDIR/<base> && python3 minibar.py -e 1 -E 5 -l 200 -M 2 -T -F PRIMER IN

    work_dir defaults to the run output dir (pre mode, matching the legacy
    layout). In post mode Minibar sees unfiltered reads, so these per-file
    trees are several times larger and the caller points work_dir at
    $PBS_JOBFS instead; only the merged result lands in the output dir.
    """
    log("")
    banner("STEP 2: Minibar demultiplexing")
    log(f" Work dir: {work_dir}")
    log(f" Files found:")
    for f in filtered_files:
        log(f"   {f}")
    log("")

    work_dir.mkdir(parents=True, exist_ok=True)
    per_file_dirs: list[Path] = []
    for input_file in filtered_files:
        base = input_file.name.removesuffix(".fastq.gz")
        outd = work_dir / base
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


def cleanup(per_file_dirs: list[Path], filtered_dir: Path | None) -> None:
    """Remove the per-file Minibar trees, and the Chopper scratch dir if the
    run made one (post mode has no pre-filter, so filtered_dir is None)."""
    log("")
    banner("Cleaning up per-file subdirectories")
    for d in per_file_dirs:
        log(f"  Removing: {d}/")
        shutil.rmtree(d, ignore_errors=True)
    log(f" Cleanup complete: {datetime.now():%Y-%m-%d %H:%M:%S}")
    log("========================================")

    if filtered_dir is None:
        return

    log("")
    banner("Deleting Chopper filtered files")
    shutil.rmtree(filtered_dir, ignore_errors=True)
    log(f" Deleted: {filtered_dir}")
    log(f" Cleanup complete: {datetime.now():%Y-%m-%d %H:%M:%S}")
    log("========================================")


# =============================================================================
# Step 6 - Organise by client + summaries
# =============================================================================
@dataclass
class ClientSample:
    client: str
    sample_id: str
    comment: str = ""  # optional, from samplesheet 'Comment' column


# Threshold below which a low-read-count comment is shown in the per-client
# summary. Set to None to disable, or change to suit future runs.
LOW_READ_THRESHOLD = 15000


def load_client_map(samplesheet: Path) -> list[ClientSample]:
    """Read Client + Sample_ID columns; sanitise client names like the awk does.

    Also reads an optional 'Comment' column if present in the samplesheet.
    Extra columns are silently ignored.
    """
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
        has_comment = "Comment" in reader.fieldnames
        for row in reader:
            client = (row.get("Client") or "").strip()
            sid = (row.get("Sample_ID") or "").strip()
            if not client or not sid:
                continue
            comment = (row.get("Comment") or "").strip() if has_comment else ""
            out.append(ClientSample(
                client=sanitise(client),
                sample_id=sanitise_sample_id(sid),
                comment=comment,
            ))
    return out


def organise_by_client(
    integrated_dir: Path,
    samples: list[ClientSample],
) -> None:
    log("")
    banner("STEP 6: Organising reads by client")

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

    # Build sample_id -> comment lookup once. Use the sanitised ID so it
    # matches what's in the per-client filenames.
    comment_by_sid: dict[str, str] = {cs.sample_id: cs.comment for cs in samples}

    # Column widths
    W_SAMPLE = 40
    W_READS = 10
    W_COMMENT = 50

    for client_subdir in sorted(p for p in integrated_dir.iterdir() if p.is_dir()):
        fastqs = sorted(client_subdir.glob("*.fastq.gz")) \
                 + sorted(client_subdir.glob("sample_*.fastq"))
        client_total = sum(count_fastq_reads(f) for f in fastqs)

        sep_short = "-" * (W_SAMPLE + 1 + W_READS)
        sep_long = "-" * (W_SAMPLE + 1 + W_READS + 1 + W_COMMENT)
        edge_long = "=" * (W_SAMPLE + 1 + W_READS + 1 + W_COMMENT)

        lines = [
            edge_long,
            f" Client: {client_subdir.name}",
            f" Date:   {datetime.now():%Y-%m-%d %H:%M:%S}",
            f" Low-read comment threshold: < {LOW_READ_THRESHOLD:,} reads",
            edge_long,
            f"{'Sample':<{W_SAMPLE}} {'Reads':>{W_READS}} {'Comment':<{W_COMMENT}}",
            sep_long,
        ]
        for f in fastqs:
            reads = count_fastq_reads(f)
            # Show the .fastq.gz extension explicitly (Path.stem strips only .gz).
            display_name = f.name
            # Look up the matching ClientSample via the file stem-without-extensions.
            sid_key = f.name.removesuffix(".fastq.gz") if f.name.endswith(".fastq.gz") \
                else f.stem.removeprefix("sample_")
            # Only annotate when below threshold; otherwise leave blank.
            comment = ""
            if reads < LOW_READ_THRESHOLD:
                comment = comment_by_sid.get(sid_key, "")
            lines.append(
                f"{display_name:<{W_SAMPLE}} {reads:>{W_READS}d} {comment:<{W_COMMENT}}"
            )
        lines += [
            sep_long,
            f"{'CLIENT TOTAL':<{W_SAMPLE}} {client_total:>{W_READS}d}",
            edge_long,
        ]
        body = "\n".join(lines)
        log(body)
        summary = client_subdir / "summary.txt"
        summary.write_text(body + "\n")
        log(f"  Saved: {summary}")
        log("")

    log(f" Organisation complete: {datetime.now():%Y-%m-%d %H:%M:%S}")
    log("========================================")


# Column widths for the extended read count table.
_W = {
    "sample": 40, "reads": 10, "pct_input": 10, "raw_demuxed": 12,
    "post_chopper": 13, "pct_retained": 13, "pass1": 20, "pct_pass1": 15,
}
_TABLE_WIDTH = sum(_W.values()) + len(_W) - 1


def write_run_summary(
    integrated_dir: Path,
    summary_file: Path,
    raw_input_dir: Path,
    filtered_dir: Path | None,
    output_dir: Path,
    total_input: int,
    filter_stage: str,
    qc_by_sid: dict[str, SampleQC] | None = None,
    chopper_by_sid: dict[str, tuple[int, int]] | None = None,
    cutadapt_by_sid: dict[str, CutadaptStats] | None = None,
) -> None:
    """
    Write read_counts_summary.txt.

    Beyond the legacy Reads / % of Input columns this tracks the read count
    through each stage:

        raw_demuxed -> post_chopper -> pct_retained
                    -> post_cutadapt_pass1 -> pct_pass1_kept

    pct_pass1_kept is the non-specificity signal. A barcode with a healthy
    raw_demuxed but a low pct_pass1_kept means Minibar found the right
    barcode on reads that are not the expected product -- correct barcode,
    wrong amplicon -- which a post_chopper count alone would hide.

    In --filter-stage pre the reads were filtered before demultiplexing, so
    no raw_demuxed count exists; that column and pct_retained are left blank
    rather than being filled with a post-filter number.
    """
    log("")
    banner("Generating read count summary")

    qc_by_sid = qc_by_sid or {}
    chopper_by_sid = chopper_by_sid or {}
    cutadapt_by_sid = cutadapt_by_sid or {}

    # collect all per-sample fastqs (top-level + one client level deep).
    # Cutadapt step renames sample_<id>.fastq -> <id>.fastq.gz; the minibar
    # catch-all bins (sample_unk, sample_Multiple_Matches) keep the old name.
    fastqs: list[Path] = []
    fastqs += sorted(integrated_dir.glob("sample_*.fastq"))
    fastqs += sorted(integrated_dir.glob("*/sample_*.fastq"))
    fastqs += sorted(integrated_dir.glob("*.fastq.gz"))
    fastqs += sorted(integrated_dir.glob("*/*.fastq.gz"))

    input_label = ("TOTAL RAW INPUT" if filter_stage == "post"
                   else "TOTAL FILTERED INPUT")
    sep = "-" * _TABLE_WIDTH

    lines = [
        "========================================",
        " Minibar Demultiplexing Read Count Summary",
        f" Raw input dir:     {raw_input_dir}",
        f" Filtered dir:      {filtered_dir if filtered_dir else '<none: post-demux filtering>'}",
        f" Out dir:           {output_dir}",
        f" Filter stage:      {filter_stage}",
        f" Date:              {datetime.now():%Y-%m-%d %H:%M:%S}",
        "========================================",
        (
            f"{'Sample':<{_W['sample']}} {'Reads':>{_W['reads']}} "
            f"{'% of Input':>{_W['pct_input']}} {'raw_demuxed':>{_W['raw_demuxed']}} "
            f"{'post_chopper':>{_W['post_chopper']}} {'pct_retained':>{_W['pct_retained']}} "
            f"{'post_cutadapt_pass1':>{_W['pass1']}} {'pct_pass1_kept':>{_W['pct_pass1']}}"
        ),
        sep,
    ]

    def fmt(value: int | None, width: int) -> str:
        """Blank rather than 0 when a stage genuinely has no number."""
        return f"{'':>{width}}" if value is None else f"{value:>{width},}"

    def fmt_pct(value: float | None, width: int) -> str:
        return f"{'':>{width}}" if value is None else f"{value:>{width - 1}.2f}%"

    total_demux = 0
    multi_match_reads = 0
    unk_reads = 0
    # tot_raw is the column total (every bin). tot_raw_chopped is the subset
    # that actually went through Chopper, and is the only honest denominator
    # for the TOTAL pct_retained -- the catch-all bins are left unfiltered.
    tot_raw = tot_raw_chopped = tot_chop = tot_pass1 = tot_pass1_in = 0
    for f in fastqs:
        # .fastq.gz files have suffixes ".fastq" + ".gz"; .stem strips only ".gz".
        sample = f.name.removesuffix(".fastq.gz") if f.name.endswith(".fastq.gz") else f.stem
        reads = count_fastq_reads(f)
        pct = (reads / total_input * 100) if total_input else 0.0

        # Sample ID as the QC/chopper/cutadapt maps key it (no sample_ prefix).
        sid = sample.removeprefix("sample_")
        chop = chopper_by_sid.get(sid)
        cut = cutadapt_by_sid.get(sid)

        if filter_stage == "post":
            raw_demuxed = chop[0] if chop else (
                qc_by_sid[sid].reads if sid in qc_by_sid else None
            )
            post_chopper = chop[1] if chop else None
        else:
            # No raw demuxed count exists in pre mode -- reads were already
            # filtered when Minibar saw them.
            raw_demuxed = None
            post_chopper = qc_by_sid[sid].reads if sid in qc_by_sid else None

        pct_retained = (
            100.0 * post_chopper / raw_demuxed
            if raw_demuxed and post_chopper is not None else None
        )
        pass1 = cut.pass1_kept if cut and not cut.error else None
        pct_pass1 = (
            100.0 * cut.pass1_kept / cut.input_reads
            if cut and not cut.error and cut.input_reads else None
        )

        lines.append(
            f"{sample:<{_W['sample']}} {reads:>{_W['reads']},} "
            f"{pct:>{_W['pct_input'] - 1}.2f}% {fmt(raw_demuxed, _W['raw_demuxed'])} "
            f"{fmt(post_chopper, _W['post_chopper'])} "
            f"{fmt_pct(pct_retained, _W['pct_retained'])} "
            f"{fmt(pass1, _W['pass1'])} {fmt_pct(pct_pass1, _W['pct_pass1'])}"
        )

        total_demux += reads
        tot_raw += raw_demuxed or 0
        tot_chop += post_chopper or 0
        if raw_demuxed and post_chopper is not None:
            tot_raw_chopped += raw_demuxed
        if cut and not cut.error:
            tot_pass1 += cut.pass1_kept
            tot_pass1_in += cut.input_reads
        if sample == "sample_Multiple_Matches":
            multi_match_reads = reads
        elif sample == "sample_unk":
            unk_reads = reads

    success_reads = total_input - multi_match_reads - unk_reads
    success_pct = (success_reads / total_input * 100) if total_input else 0.0
    total_pct = (total_demux / total_input * 100) if total_input else 0.0

    tot_raw_disp = tot_raw if filter_stage == "post" and tot_raw else None
    tot_chop_disp = tot_chop or None
    tot_retained = (
        100.0 * tot_chop / tot_raw_chopped if tot_raw_chopped and tot_chop else None
    )
    tot_pass1_disp = tot_pass1 or None
    tot_pct_pass1 = (100.0 * tot_pass1 / tot_pass1_in) if tot_pass1_in else None

    lines += [
        sep,
        (
            f"{'TOTAL DEMULTIPLEXED':<{_W['sample']}} {total_demux:>{_W['reads']},} "
            f"{total_pct:>{_W['pct_input'] - 1}.2f}% {fmt(tot_raw_disp, _W['raw_demuxed'])} "
            f"{fmt(tot_chop_disp, _W['post_chopper'])} "
            f"{fmt_pct(tot_retained, _W['pct_retained'])} "
            f"{fmt(tot_pass1_disp, _W['pass1'])} "
            f"{fmt_pct(tot_pct_pass1, _W['pct_pass1'])}"
        ),
        f"{'SUCCESSFULLY DEMULTIPLEXED':<{_W['sample']}} {success_reads:>{_W['reads']},} "
        f"{success_pct:>{_W['pct_input'] - 1}.2f}%",
        f"{input_label:<{_W['sample']}} {total_input:>{_W['reads']},} "
        f"{100.00:>{_W['pct_input'] - 1}.2f}%",
        "=" * _TABLE_WIDTH,
        "",
        " Columns:",
        "   Reads                final read count in the delivered fastq",
        f"   % of Input           share of {input_label.lower()}",
        "   raw_demuxed          reads Minibar assigned, before filtering"
        + ("" if filter_stage == "post"
           else " (n/a: filtered pre-demux)"),
        "   post_chopper         reads surviving the Chopper length/Q filter",
        "   pct_retained         post_chopper / raw_demuxed",
        "   post_cutadapt_pass1  reads with a findable forward primer",
        "   pct_pass1_kept       post_cutadapt_pass1 / cutadapt input;",
        "                        low here with a healthy raw_demuxed means",
        "                        correct barcode, wrong product",
    ]
    if filter_stage == "post":
        lines += [
            "",
            " Note: sample_unk and sample_Multiple_Matches are left unfiltered --",
            "       their raw length profile is what makes them diagnostic. They",
            "       are therefore excluded from the TOTAL pct_retained denominator.",
        ]
    lines.append("=" * _TABLE_WIDTH)

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
    # default None (not "pre") so an explicit --filter-stage pre can be told
    # apart from an unset one when --diagnostic is also given.
    p.add_argument("--filter-stage", choices=("pre", "post"), default=None,
                   help="Run Chopper before demultiplexing (pre, default, "
                        "legacy behaviour) or after it, per sample (post). "
                        "In post mode raw fastq_pass files go straight to "
                        "Minibar with no pre-filter at all.")
    p.add_argument("--minibar-work-dir", type=Path, default=None,
                   help="Where Minibar's per-file working directories are "
                        "created (default: the output dir). Point this at "
                        "$PBS_JOBFS in post mode -- unfiltered input makes "
                        "these trees several times larger, and only the "
                        "merged result needs to reach the output dir.")
    p.add_argument("--chopper", type=Path, default=DEFAULT_CHOPPER,
                   help=f"Path to chopper binary (default: {DEFAULT_CHOPPER}).")
    p.add_argument("--minibar", type=Path, default=DEFAULT_MINIBAR,
                   help=f"Path to minibar.py (default: {DEFAULT_MINIBAR}).")

    # Chopper filter parameters. --minlength/--maxlength match chopper's own
    # flag names; --min-length/--max-length are kept as aliases so existing
    # invocations do not break.
    p.add_argument("--min-quality", "-q", type=int, default=None,
                   help="Chopper -q minimum mean read quality (default: 15).")
    p.add_argument("--minlength", "--min-length", dest="minlength",
                   type=int, default=None,
                   help="Chopper --minlength (default: 1000).")
    p.add_argument("--maxlength", "--max-length", dest="maxlength",
                   type=int, default=None,
                   help="Chopper --maxlength (default: 2000).")

    p.add_argument("--diagnostic", action="store_true",
                   help="Shorthand for --filter-stage post --minlength 100 "
                        "--maxlength 5000, so a first pass observes the real "
                        "length distribution instead of assuming the window. "
                        "Explicitly given flags still win.")
    p.add_argument("--no-demux-qc", action="store_true",
                   help="Skip the per-sample demux QC pass (Step 3).")

    p.add_argument("--cutadapt", type=Path, default=DEFAULT_CUTADAPT,
                   help=f"Path to cutadapt binary (default: {DEFAULT_CUTADAPT}).")
    p.add_argument("--cutadapt-threads", type=int, default=4,
                   help="Threads for cutadapt (default: 4).")
    p.add_argument("--cutadapt-error-rate", type=float, default=0.2,
                   help="Cutadapt error rate -e (default: 0.2).")

    args = p.parse_args()
    _resolve_filter_settings(args, p)
    return args


# Standard filter window, and the wide window --diagnostic swaps in.
_DEFAULTS = {"min_quality": 15, "minlength": 1000, "maxlength": 2000}
_DIAGNOSTIC = {"filter_stage": "post", "minlength": 100, "maxlength": 5000}


def _resolve_filter_settings(args: argparse.Namespace,
                             parser: argparse.ArgumentParser) -> None:
    """
    Apply precedence: explicit flag > --diagnostic > default.

    The filter args default to None so an explicitly passed value is
    distinguishable from an unset one, which is what lets --diagnostic
    supply values without silently overriding the user.
    """
    if args.diagnostic:
        if args.filter_stage == "pre":
            parser.error(
                "--diagnostic implies --filter-stage post, but "
                "--filter-stage pre was given explicitly"
            )
        args.filter_stage = _DIAGNOSTIC["filter_stage"]
        for key in ("minlength", "maxlength"):
            if getattr(args, key) is None:
                setattr(args, key, _DIAGNOSTIC[key])

    if args.filter_stage is None:
        args.filter_stage = "pre"

    for key, value in _DEFAULTS.items():
        if getattr(args, key) is None:
            setattr(args, key, value)

    if args.minlength > args.maxlength:
        parser.error(
            f"--minlength ({args.minlength}) exceeds --maxlength ({args.maxlength})"
        )


def collect_logs(output_dir: Path) -> Path:
    """Move all *.log and *.txt summaries into a dated run_log_<date> directory."""
    log("")
    banner("STEP 7: Collecting logs")
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
    integrated_dir = output_dir / "integrated_demultiplexing"
    qc_dir = output_dir / "demux_qc"
    summary_file = output_dir / "read_counts_summary.txt"
    minibar_work_dir = args.minibar_work_dir or output_dir

    output_dir.mkdir(parents=True, exist_ok=True)
    integrated_dir.mkdir(parents=True, exist_ok=True)

    log(f" Primer file:  {args.primer_file}")
    log(f" Filter stage: {args.filter_stage}")
    log(f" Chopper:      -q {args.min_quality} --minlength {args.minlength} "
        f"--maxlength {args.maxlength}"
        + ("   [--diagnostic]" if args.diagnostic else ""))

    # ---- Step 1: Chopper, pre-demux only ----
    if args.filter_stage == "pre":
        filtered_dir = output_dir / "chopper_filtered"
        minibar_inputs = run_chopper(
            args.raw_input_dir, filtered_dir, args.chopper,
            args.min_quality, args.minlength, args.maxlength,
        )
        total_input = count_filtered_total(minibar_inputs)
    else:
        # No pre-filter: fastq_pass is already MinKNOW min-qscore filtered,
        # and the whole point of post mode is to let Minibar see everything.
        filtered_dir = None
        minibar_inputs = sorted(args.raw_input_dir.glob("*.fastq.gz"))
        if not minibar_inputs:
            sys.exit(f"ERROR: no *.fastq.gz files found in {args.raw_input_dir}")
        log("")
        banner("STEP 1: Chopper filtering SKIPPED (--filter-stage post)")
        log(f" Feeding {len(minibar_inputs)} raw file(s) straight to Minibar.")
        log("========================================")
        total_input = count_raw_total(minibar_inputs)

    # ---- Step 2: Minibar ----
    per_file_dirs = run_minibar(
        minibar_inputs, minibar_work_dir, args.primer_file, args.minibar,
    )
    merge_per_file_outputs(per_file_dirs, integrated_dir)
    cleanup(per_file_dirs, filtered_dir)

    # ---- Step 3: per-sample QC of the demultiplexed reads ----
    # In post mode this is the only place the pre-filter distribution is
    # observable, so it has to run before Step 4 touches the files.
    qc_by_sid: dict[str, SampleQC] = {}
    if not args.no_demux_qc:
        qc_by_sid = run_demux_qc(
            integrated_dir=integrated_dir,
            qc_dir=qc_dir,
            minlength=args.minlength,
            maxlength=args.maxlength,
            filter_stage=args.filter_stage,
            log=log,
        )

    # ---- Step 4: Chopper per sample, post-demux only ----
    chopper_by_sid: dict[str, tuple[int, int]] = {}
    if args.filter_stage == "post":
        chopper_by_sid = run_chopper_per_sample(
            integrated_dir=integrated_dir,
            chopper_prog=args.chopper,
            min_quality=args.min_quality,
            min_length=args.minlength,
            max_length=args.maxlength,
            qc_by_sid=qc_by_sid,
        )
    else:
        # Announced rather than silently absent, so the step numbering lines
        # up when a pre log and a post log are read side by side.
        log("")
        banner("STEP 4: Per-sample Chopper SKIPPED (--filter-stage pre)")
        log(" Chopper already ran at step 1, before demultiplexing.")
        log("========================================")

    # ---- Step 5: two-pass cutadapt ----
    _, cutadapt_stats = run_cutadapt_step(
        integrated_dir=integrated_dir,
        primer_file=args.primer_file,
        output_dir=output_dir,
        cutadapt=args.cutadapt,
        threads=args.cutadapt_threads,
        error_rate=args.cutadapt_error_rate,
    )
    cutadapt_by_sid = {s.sample_id: s for s in cutadapt_stats}

    # ---- Step 6 ----
    samples = load_client_map(args.samplesheet)
    organise_by_client(integrated_dir, samples)
    write_run_summary(
        integrated_dir, summary_file,
        args.raw_input_dir, filtered_dir, output_dir,
        total_input,
        filter_stage=args.filter_stage,
        qc_by_sid=qc_by_sid,
        chopper_by_sid=chopper_by_sid,
        cutadapt_by_sid=cutadapt_by_sid,
    )

    # ---- Step 7: gather logs into run_log_<date>/ ----
    collect_logs(output_dir)


if __name__ == "__main__":
    main()
