#!/usr/bin/env python3
"""
cutadapt_2pass.py

Two-pass cutadapt for ONT 16S amplicon reads, per sample.

This module is invoked from filter_chopper_demux_minibar.py after minibar
demultiplexing and before "organise by client". For each demultiplexed
sample_<SampleID>.fastq in the integrated_demultiplexing directory:

    Pass 1 (orient + 5' trim, STRICT)
        -g FWD  --revcomp  --rename={header}  --discard-untrimmed  -e 0.2
        Reads where FWD primer cannot be found are dropped.

    Pass 2 (3' trim, TOLERANT)
        -a REV_RC  -e 0.2
        Trims the reverse-complemented reverse construct from the 3' end.
        Reads where REV_RC is not found are kept (truncated reads retained).

Per-sample primers come from the primer setup file produced by
generate_primer_setup.py (columns: SampleID, FwIndex, FwPrimer,
RvIndex, RvPrimer). FwPrimer / RvPrimer already encode the
adapter+index+primer construct, so we use them as-is.

Output naming:
    integrated_demultiplexing/sample_<SampleID>.fastq      (input)
    integrated_demultiplexing/<SampleID>.fastq.gz          (output, replaces input)

Stats: a per-sample summary (kept / discarded / reverse-complemented /
REV_RC trimmed) is written to <output_dir>/cutadapt_summary.txt.

Logs:
    <output_dir>/cutadapt_logs/<SampleID>.pass1_orient.log
    <output_dir>/cutadapt_logs/<SampleID>.pass2_trim3.log
    (JSON files are read for stats then deleted.)
"""

from __future__ import annotations

import csv
import gzip
import json
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


# --- Defaults (Gadi) ---------------------------------------------------------
DEFAULT_CUTADAPT = Path("/g/data/vz35/ONT_16s_workflow/tools/cutadapt-env/bin/cutadapt")

# IUPAC complement table for reverse-complementing primer sequences
_IUPAC_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)


# =============================================================================
# Helpers (kept self-contained so this module is import-friendly)
# =============================================================================
def _log(msg: str = "") -> None:
    print(msg, flush=True)


def _banner(title: str) -> None:
    _log("========================================")
    _log(f" {title}: {datetime.now():%Y-%m-%d %H:%M:%S}")
    _log("========================================")


def reverse_complement(seq: str) -> str:
    """Reverse complement a DNA sequence supporting IUPAC ambiguity codes."""
    return seq.translate(_IUPAC_COMPLEMENT)[::-1]


# =============================================================================
# Primer setup file
# =============================================================================
@dataclass
class SamplePrimers:
    sample_id: str
    fw_primer: str   # full forward construct: adapter + FwIdx + 27F (etc.)
    rv_primer: str   # full reverse construct: adapter + RvIdx + 1492R (etc.)

    @property
    def rv_primer_rc(self) -> str:
        return reverse_complement(self.rv_primer)


def load_primer_setup(primer_file: Path) -> dict[str, SamplePrimers]:
    """
    Load the primer setup TSV produced by generate_primer_setup.py.
    Returns dict keyed by SampleID.
    Required columns: SampleID, FwPrimer, RvPrimer (FwIndex, RvIndex ignored).
    """
    if not primer_file.is_file():
        sys.exit(f"ERROR: primer setup file not found: {primer_file}")

    out: dict[str, SamplePrimers] = {}
    with primer_file.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"SampleID", "FwPrimer", "RvPrimer"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            sys.exit(
                f"ERROR: {primer_file} is missing column(s): "
                f"{', '.join(sorted(missing))}"
            )
        for row in reader:
            sid = (row.get("SampleID") or "").strip()
            fw = (row.get("FwPrimer") or "").strip()
            rv = (row.get("RvPrimer") or "").strip()
            if not sid:
                continue
            if not fw or not rv:
                sys.exit(
                    f"ERROR: primer file row for {sid} has empty FwPrimer / RvPrimer"
                )
            out[sid] = SamplePrimers(sample_id=sid, fw_primer=fw, rv_primer=rv)
    if not out:
        sys.exit(f"ERROR: no samples loaded from {primer_file}")
    return out


# =============================================================================
# Per-sample two-pass cutadapt
# =============================================================================
@dataclass
class CutadaptStats:
    sample_id: str
    input_reads: int = 0
    pass1_kept: int = 0
    pass1_discarded: int = 0
    pass1_revcomped: int = 0
    pass2_input: int = 0
    pass2_trimmed3: int = 0  # reads where REV_RC was found and trimmed
    pass2_output: int = 0
    bp_input: int = 0
    bp_output: int = 0
    error: str = ""  # populated on failure; sample is skipped from output then

    # display percentages
    def pct_kept(self) -> float:
        return 100.0 * self.pass2_output / self.input_reads if self.input_reads else 0.0

    def pct_discarded(self) -> float:
        return 100.0 * self.pass1_discarded / self.input_reads if self.input_reads else 0.0

    def pct_revcomped(self) -> float:
        return 100.0 * self.pass1_revcomped / self.pass1_kept if self.pass1_kept else 0.0

    def pct_trimmed3(self) -> float:
        return 100.0 * self.pass2_trimmed3 / self.pass2_input if self.pass2_input else 0.0


def _parse_pass1_json(p: Path) -> tuple[int, int, int, int, int]:
    """Return (input, output, discarded, revcomped, bp_input)."""
    with p.open() as fh:
        d = json.load(fh)
    rc = d["read_counts"]
    return (
        rc["input"],
        rc["output"],
        rc["filtered"].get("discard_untrimmed", 0),
        rc.get("reverse_complemented", 0),
        d["basepair_counts"]["input"],
    )


def _parse_pass2_json(p: Path) -> tuple[int, int, int, int]:
    """Return (input, output, trimmed3, bp_output)."""
    with p.open() as fh:
        d = json.load(fh)
    rc = d["read_counts"]
    trimmed = 0
    for a in d.get("adapters_read1", []) or []:
        trimmed += a.get("total_matches", 0)
    return rc["input"], rc["output"], trimmed, d["basepair_counts"]["output"]


def _run_pass1(
    cutadapt: Path,
    fw_primer: str,
    input_fq: Path,
    out_fq: Path,
    log_file: Path,
    json_file: Path,
    threads: int,
    error_rate: float,
) -> None:
    cmd = [
        str(cutadapt),
        "-g", fw_primer,
        "--revcomp",
        "--rename={header}",
        "--discard-untrimmed",
        "-e", str(error_rate),
        "-j", str(threads),
        "--json", str(json_file),
        "-o", str(out_fq),
        str(input_fq),
    ]
    with log_file.open("w") as lf:
        rc = subprocess.call(cmd, stdout=lf, stderr=subprocess.STDOUT)
    if rc != 0:
        raise RuntimeError(f"cutadapt pass 1 failed (exit {rc}); see {log_file}")


def _run_pass2(
    cutadapt: Path,
    rv_primer_rc: str,
    input_fq: Path,
    out_fq: Path,
    log_file: Path,
    json_file: Path,
    threads: int,
    error_rate: float,
) -> None:
    cmd = [
        str(cutadapt),
        "-a", rv_primer_rc,
        "-e", str(error_rate),
        "-j", str(threads),
        "--json", str(json_file),
        "-o", str(out_fq),
        str(input_fq),
    ]
    with log_file.open("w") as lf:
        rc = subprocess.call(cmd, stdout=lf, stderr=subprocess.STDOUT)
    if rc != 0:
        raise RuntimeError(f"cutadapt pass 2 failed (exit {rc}); see {log_file}")


def process_sample(
    primers: SamplePrimers,
    sample_fq: Path,
    final_fq: Path,
    log_dir: Path,
    tmp_dir: Path,
    cutadapt: Path,
    threads: int,
    error_rate: float,
) -> CutadaptStats:
    """Run two-pass cutadapt on a single sample. Returns stats."""
    stats = CutadaptStats(sample_id=primers.sample_id)

    # intermediate files
    oriented_fq = tmp_dir / f"{primers.sample_id}.oriented.fastq.gz"
    p1_json = tmp_dir / f"{primers.sample_id}.pass1_orient.json"
    p2_json = tmp_dir / f"{primers.sample_id}.pass2_trim3.json"
    p1_log = log_dir / f"{primers.sample_id}.pass1_orient.log"
    p2_log = log_dir / f"{primers.sample_id}.pass2_trim3.log"

    try:
        # ---- Pass 1: orient + 5' trim --------------------------------------
        _run_pass1(
            cutadapt=cutadapt,
            fw_primer=primers.fw_primer,
            input_fq=sample_fq,
            out_fq=oriented_fq,
            log_file=p1_log,
            json_file=p1_json,
            threads=threads,
            error_rate=error_rate,
        )
        in1, out1, disc1, revc1, bp_in = _parse_pass1_json(p1_json)
        stats.input_reads = in1
        stats.pass1_kept = out1
        stats.pass1_discarded = disc1
        stats.pass1_revcomped = revc1
        stats.bp_input = bp_in

        # ---- Pass 2: 3' trim, tolerant -------------------------------------
        _run_pass2(
            cutadapt=cutadapt,
            rv_primer_rc=primers.rv_primer_rc,
            input_fq=oriented_fq,
            out_fq=final_fq,
            log_file=p2_log,
            json_file=p2_json,
            threads=threads,
            error_rate=error_rate,
        )
        in2, out2, trim3, bp_out = _parse_pass2_json(p2_json)
        stats.pass2_input = in2
        stats.pass2_output = out2
        stats.pass2_trimmed3 = trim3
        stats.bp_output = bp_out

    except Exception as e:
        stats.error = str(e)

    # ---- Always clean up intermediates ------------------------------------
    for p in (oriented_fq, p1_json, p2_json):
        try:
            p.unlink(missing_ok=True)
        except Exception:
            pass

    return stats


# =============================================================================
# Driver: run cutadapt over all sample_*.fastq in integrated_demultiplexing
# =============================================================================
def run_cutadapt_step(
    integrated_dir: Path,
    primer_file: Path,
    output_dir: Path,
    cutadapt: Path = DEFAULT_CUTADAPT,
    threads: int = 4,
    error_rate: float = 0.2,
) -> Path:
    """
    For every sample_<SampleID>.fastq in integrated_dir whose SampleID is in
    the primer setup file:
        - run two-pass cutadapt
        - write <SampleID>.fastq.gz alongside (replacing the sample_*.fastq)
        - remove the source sample_*.fastq on success
    sample_unk.fastq and sample_Multiple_Matches.fastq are left untouched.

    Returns the path to the cutadapt summary file written to output_dir.
    """
    _log("")
    _banner("STEP 2.5: Cutadapt two-pass (orient + trim)")

    if not cutadapt.is_file():
        sys.exit(f"ERROR: cutadapt not found at {cutadapt}")
    primers_by_id = load_primer_setup(primer_file)
    _log(f" Loaded primers for {len(primers_by_id)} samples from {primer_file}")
    _log(f" Cutadapt binary:  {cutadapt}")
    _log(f" Integrated dir:   {integrated_dir}")
    _log("")

    log_dir = output_dir / "cutadapt_logs"
    tmp_dir = output_dir / "cutadapt_tmp"
    log_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    all_stats: list[CutadaptStats] = []
    missing_primers: list[str] = []
    skipped_files: list[str] = []

    # Iterate over sample_*.fastq currently in integrated_demultiplexing.
    # NOTE: at this point in the pipeline reads have not yet been moved into
    # per-client subdirectories, so all per-sample files sit at the top level.
    sample_fastqs = sorted(integrated_dir.glob("sample_*.fastq"))
    if not sample_fastqs:
        _log("  WARNING: no sample_*.fastq files found; nothing to do.")
        return _write_summary(output_dir, all_stats, missing_primers, skipped_files)

    for fq in sample_fastqs:
        stem = fq.stem  # e.g. "sample_PPCS-4"
        if not stem.startswith("sample_"):
            continue
        sample_id = stem[len("sample_"):]

        # Skip minibar's catch-all bins
        if sample_id in ("unk", "Multiple_Matches"):
            skipped_files.append(fq.name)
            continue

        primers = primers_by_id.get(sample_id)
        if primers is None:
            _log(f"  WARNING: no primers in setup for {sample_id}; leaving {fq.name} untouched.")
            missing_primers.append(sample_id)
            continue

        final_fq = integrated_dir / f"{sample_id}.fastq.gz"

        _log(f"--- {sample_id} ---")
        _log(f"  Input:  {fq}")
        _log(f"  Output: {final_fq}")
        _log(f"  Start:  {datetime.now():%Y-%m-%d %H:%M:%S}")

        stats = process_sample(
            primers=primers,
            sample_fq=fq,
            final_fq=final_fq,
            log_dir=log_dir,
            tmp_dir=tmp_dir,
            cutadapt=cutadapt,
            threads=threads,
            error_rate=error_rate,
        )

        if stats.error:
            _log(f"  ERROR: {stats.error}")
        else:
            _log(
                f"  Done:   in={stats.input_reads:,} "
                f"kept={stats.pass2_output:,} ({stats.pct_kept():.2f}%) "
                f"discarded={stats.pass1_discarded:,} ({stats.pct_discarded():.2f}%) "
                f"revcomped={stats.pass1_revcomped:,} ({stats.pct_revcomped():.2f}%) "
                f"3'trim={stats.pass2_trimmed3:,} ({stats.pct_trimmed3():.2f}%)"
            )
            # Replace the demultiplexed plaintext fastq with the trimmed gz
            try:
                fq.unlink()
            except Exception as e:
                _log(f"  WARNING: could not remove source {fq}: {e}")

        all_stats.append(stats)
        _log("")

    # Remove now-empty tmp dir
    try:
        tmp_dir.rmdir()
    except OSError:
        pass  # not empty, leave it

    summary_path = _write_summary(output_dir, all_stats, missing_primers, skipped_files)
    _log(f" Cutadapt step complete: {datetime.now():%Y-%m-%d %H:%M:%S}")
    _log("========================================")
    return summary_path


def _write_summary(
    output_dir: Path,
    stats_list: list[CutadaptStats],
    missing_primers: list[str],
    skipped_files: list[str],
) -> Path:
    summary_path = output_dir / "cutadapt_summary.txt"

    lines = [
        "========================================",
        " Cutadapt Two-Pass Summary",
        f" Date:  {datetime.now():%Y-%m-%d %H:%M:%S}",
        "========================================",
        f"{'Sample':<30} {'Input':>10} {'Kept':>10} {'%Kept':>8} "
        f"{'Disc':>8} {'%Disc':>8} {'RevC':>8} {'%RevC':>8} "
        f"{'Trim3':>8} {'%Trim3':>8}",
        "-" * 120,
    ]
    tot_in = tot_kept = tot_disc = tot_revc = tot_trim3 = 0
    for s in stats_list:
        if s.error:
            lines.append(f"{s.sample_id:<30}  ERROR: {s.error}")
            continue
        lines.append(
            f"{s.sample_id:<30} "
            f"{s.input_reads:>10,} {s.pass2_output:>10,} {s.pct_kept():>7.2f}% "
            f"{s.pass1_discarded:>8,} {s.pct_discarded():>7.2f}% "
            f"{s.pass1_revcomped:>8,} {s.pct_revcomped():>7.2f}% "
            f"{s.pass2_trimmed3:>8,} {s.pct_trimmed3():>7.2f}%"
        )
        tot_in += s.input_reads
        tot_kept += s.pass2_output
        tot_disc += s.pass1_discarded
        tot_revc += s.pass1_revcomped
        tot_trim3 += s.pass2_trimmed3

    pct = lambda n, d: (100.0 * n / d) if d else 0.0
    lines += [
        "-" * 120,
        f"{'TOTAL':<30} "
        f"{tot_in:>10,} {tot_kept:>10,} {pct(tot_kept, tot_in):>7.2f}% "
        f"{tot_disc:>8,} {pct(tot_disc, tot_in):>7.2f}% "
        f"{tot_revc:>8,} {pct(tot_revc, tot_kept):>7.2f}% "
        f"{tot_trim3:>8,} {pct(tot_trim3, tot_kept):>7.2f}%",
        "========================================",
    ]
    if missing_primers:
        lines.append(" Samples with no primers in setup (left untouched):")
        for s in missing_primers:
            lines.append(f"   - {s}")
        lines.append("")
    if skipped_files:
        lines.append(" Files skipped (minibar catch-all bins):")
        for s in skipped_files:
            lines.append(f"   - {s}")
        lines.append("")

    body = "\n".join(lines)
    _log(body)
    summary_path.write_text(body + "\n")
    _log(f" Cutadapt summary saved to: {summary_path}")
    return summary_path
