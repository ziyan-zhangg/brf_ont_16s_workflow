#!/usr/bin/env python3
"""
demux_qc.py

Per-sample QC of freshly demultiplexed reads, run *before* any length/quality
filtering. This is the step that makes the two-stage read loss separable: it
records what Minibar actually assigned to each barcode, so a later Chopper
pass can be scored against a real distribution instead of an assumed one.

Invoked from filter_chopper_demux_minibar.py between the merge step and the
per-sample Chopper step, but also runnable standalone on an existing run
directory (see --help).

For every sample_<SampleID>.fastq in the integrated demultiplexing directory
-- including Minibar's catch-all bins sample_unk and sample_Multiple_Matches
-- it reports:

    * raw demuxed read count
    * length histogram in 50 bp bins, plus N50 and modal length
    * mean read Q (per-read mean error probability -> Phred, averaged)
    * % of reads falling inside [minlength, maxlength]

Outputs, all under <output_dir>/demux_qc/:

    <SampleID>.tsv            one row per 50 bp bin, scalar stats in a
                              '#'-prefixed header block (pandas: comment='#')
    demux_qc_summary.tsv      one row per sample, all scalar stats

Note on --filter-stage pre: the reads reaching this step have already been
through Chopper, so the distributions describe *filtered* reads. The header
block and the summary table both record which stage produced them so the two
orderings can be compared directly on the same input.
"""

from __future__ import annotations

import argparse
import gzip
import math
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import BinaryIO, Callable

# Width of the length histogram bins, in bp.
BIN_WIDTH = 50

# Minibar's catch-all output bins. These are not real samples but are still
# reported: an inflated sample_unk is the first sign of a demux problem, and
# in post mode it is the only place the unassigned length profile survives.
CATCH_ALL = {
    "unk": "unassigned",
    "Multiple_Matches": "multiple_matches",
}

# Error probability for every possible Phred+33 byte, indexed by byte value.
# Built once so the per-read loop is a table lookup rather than a pow().
_ERR_PROB = [10.0 ** (-(b - 33) / 10.0) for b in range(256)]


# =============================================================================
# Helpers
# =============================================================================
def _log(msg: str = "") -> None:
    print(msg, flush=True)


def _open_fastq(path: Path) -> BinaryIO:
    return gzip.open(path, "rb") if path.suffix == ".gz" else open(path, "rb")


def sample_id_from_path(path: Path) -> str:
    """sample_<SampleID>.fastq[.gz] -> <SampleID>."""
    name = path.name
    name = name.removesuffix(".gz").removesuffix(".fastq")
    return name.removeprefix("sample_")


# =============================================================================
# Per-sample statistics
# =============================================================================
@dataclass
class SampleQC:
    sample_id: str
    category: str = "sample"          # sample | unassigned | multiple_matches
    reads: int = 0
    total_bp: int = 0
    min_length: int = 0
    max_length: int = 0
    mean_length: float = 0.0
    median_length: int = 0
    n50: int = 0
    modal_length: int = 0             # most frequent exact read length
    modal_bin_start: int = 0          # most populous 50 bp bin
    modal_bin_count: int = 0
    mean_read_q: float = 0.0
    window_min: int = 0
    window_max: int = 0
    reads_in_window: int = 0
    # exact length -> count; the 50 bp histogram is derived from this so that
    # N50/median/mode stay exact without holding one int per read in memory.
    length_counts: dict[int, int] = field(default_factory=dict, repr=False)

    @property
    def pct_in_window(self) -> float:
        return 100.0 * self.reads_in_window / self.reads if self.reads else 0.0

    def histogram(self) -> list[tuple[int, int]]:
        """(bin_start, count) for every 50 bp bin from the shortest to the
        longest read, including empty bins so the profile stays plottable."""
        if not self.length_counts:
            return []
        binned: dict[int, int] = {}
        for length, n in self.length_counts.items():
            binned[(length // BIN_WIDTH) * BIN_WIDTH] = \
                binned.get((length // BIN_WIDTH) * BIN_WIDTH, 0) + n
        lo = min(binned)
        hi = max(binned)
        return [(b, binned.get(b, 0)) for b in range(lo, hi + BIN_WIDTH, BIN_WIDTH)]


def scan_fastq(
    path: Path,
    sample_id: str,
    minlength: int,
    maxlength: int,
) -> SampleQC:
    """Single streaming pass over a fastq, accumulating every reported stat."""
    qc = SampleQC(
        sample_id=sample_id,
        category=CATCH_ALL.get(sample_id, "sample"),
        window_min=minlength,
        window_max=maxlength,
    )

    lengths: dict[int, int] = {}
    q_sum = 0.0
    err_prob = _ERR_PROB.__getitem__

    with _open_fastq(path) as fh:
        readline = fh.readline
        while True:
            header = readline()
            if not header:
                break
            seq = readline().rstrip()
            readline()                      # '+' separator
            qual = readline().rstrip()
            if not qual:
                # truncated final record; ignore it rather than miscounting
                break

            length = len(seq)
            qc.reads += 1
            qc.total_bp += length
            lengths[length] = lengths.get(length, 0) + 1
            if minlength <= length <= maxlength:
                qc.reads_in_window += 1

            if qual:
                mean_err = sum(map(err_prob, qual)) / len(qual)
                # Phred caps at ~93; a perfect-quality read would divide by 0.
                q_sum += -10.0 * math.log10(mean_err) if mean_err > 0 else 93.0

    if qc.reads == 0:
        return qc

    qc.length_counts = lengths
    qc.min_length = min(lengths)
    qc.max_length = max(lengths)
    qc.mean_length = qc.total_bp / qc.reads
    qc.mean_read_q = q_sum / qc.reads

    # Median length: walk the exact-length counter in ascending order.
    half_reads = qc.reads / 2.0
    seen = 0
    for length in sorted(lengths):
        seen += lengths[length]
        if seen >= half_reads:
            qc.median_length = length
            break

    # N50: shortest length L such that reads >= L hold half the total bases.
    half_bp = qc.total_bp / 2.0
    seen_bp = 0
    for length in sorted(lengths, reverse=True):
        seen_bp += length * lengths[length]
        if seen_bp >= half_bp:
            qc.n50 = length
            break

    qc.modal_length = max(lengths, key=lambda L: (lengths[L], -L))
    hist = qc.histogram()
    qc.modal_bin_start, qc.modal_bin_count = max(hist, key=lambda bc: (bc[1], -bc[0]))
    return qc


# =============================================================================
# Output
# =============================================================================
_SUMMARY_COLUMNS = [
    "sample", "category", "reads", "total_bp", "min_length", "max_length",
    "mean_length", "median_length", "n50", "modal_length", "modal_bin",
    "modal_bin_count", "mean_read_q", "window_min", "window_max",
    "reads_in_window", "pct_in_window",
]


def _summary_row(qc: SampleQC) -> list[str]:
    return [
        qc.sample_id,
        qc.category,
        str(qc.reads),
        str(qc.total_bp),
        str(qc.min_length),
        str(qc.max_length),
        f"{qc.mean_length:.1f}",
        str(qc.median_length),
        str(qc.n50),
        str(qc.modal_length),
        f"{qc.modal_bin_start}-{qc.modal_bin_start + BIN_WIDTH - 1}",
        str(qc.modal_bin_count),
        f"{qc.mean_read_q:.2f}",
        str(qc.window_min),
        str(qc.window_max),
        str(qc.reads_in_window),
        f"{qc.pct_in_window:.2f}",
    ]


def write_sample_tsv(qc: SampleQC, dest: Path, filter_stage: str, source: Path) -> None:
    """Histogram table with the scalar stats in a '#' header block."""
    stage_note = (
        "reads as demultiplexed, before any filtering"
        if filter_stage == "post"
        else "reads already Chopper-filtered before demultiplexing"
    )
    header = [
        f"# sample\t{qc.sample_id}",
        f"# category\t{qc.category}",
        f"# source\t{source}",
        f"# filter_stage\t{filter_stage}",
        f"# stage_note\t{stage_note}",
        f"# generated\t{datetime.now():%Y-%m-%d %H:%M:%S}",
        f"# bin_width_bp\t{BIN_WIDTH}",
        f"# reads\t{qc.reads}",
        f"# total_bp\t{qc.total_bp}",
        f"# min_length\t{qc.min_length}",
        f"# max_length\t{qc.max_length}",
        f"# mean_length\t{qc.mean_length:.1f}",
        f"# median_length\t{qc.median_length}",
        f"# n50\t{qc.n50}",
        f"# modal_length\t{qc.modal_length}",
        f"# modal_bin\t{qc.modal_bin_start}-{qc.modal_bin_start + BIN_WIDTH - 1}",
        f"# mean_read_q\t{qc.mean_read_q:.2f}",
        f"# length_window\t{qc.window_min}-{qc.window_max}",
        f"# reads_in_window\t{qc.reads_in_window}",
        f"# pct_in_window\t{qc.pct_in_window:.2f}",
        "bin_start\tbin_end\tcount\tpct_of_reads\tcum_pct",
    ]

    lines = list(header)
    cum = 0
    for bin_start, count in qc.histogram():
        cum += count
        pct = 100.0 * count / qc.reads if qc.reads else 0.0
        cum_pct = 100.0 * cum / qc.reads if qc.reads else 0.0
        lines.append(
            f"{bin_start}\t{bin_start + BIN_WIDTH - 1}\t{count}\t{pct:.3f}\t{cum_pct:.3f}"
        )

    dest.write_text("\n".join(lines) + "\n")


def write_summary_tsv(qcs: list[SampleQC], dest: Path) -> None:
    lines = ["\t".join(_SUMMARY_COLUMNS)]
    # real samples first, catch-all bins last
    for qc in sorted(qcs, key=lambda q: (q.category != "sample", q.sample_id)):
        lines.append("\t".join(_summary_row(qc)))
    dest.write_text("\n".join(lines) + "\n")


# =============================================================================
# Driver
# =============================================================================
def run_demux_qc(
    integrated_dir: Path,
    qc_dir: Path,
    minlength: int,
    maxlength: int,
    filter_stage: str = "post",
    log: Callable[[str], None] = _log,
) -> dict[str, SampleQC]:
    """
    QC every sample_*.fastq in integrated_dir. Returns {sample_id: SampleQC}.

    Read-only with respect to the fastq files -- it never rewrites reads, so
    it is safe to run in either filter stage.
    """
    log("")
    log("========================================")
    log(f" STEP 3: Per-sample demux QC: {datetime.now():%Y-%m-%d %H:%M:%S}")
    log("========================================")
    log(f" Input dir:     {integrated_dir}")
    log(f" QC dir:        {qc_dir}")
    log(f" Filter stage:  {filter_stage}")
    log(f" Length window: {minlength}-{maxlength} bp")
    log("")

    sample_fastqs = sorted(integrated_dir.glob("sample_*.fastq"))
    if not sample_fastqs:
        log("  WARNING: no sample_*.fastq found; skipping QC.")
        return {}

    qc_dir.mkdir(parents=True, exist_ok=True)

    results: dict[str, SampleQC] = {}
    for fq in sample_fastqs:
        sid = sample_id_from_path(fq)
        qc = scan_fastq(fq, sid, minlength, maxlength)
        results[sid] = qc

        write_sample_tsv(qc, qc_dir / f"{sid}.tsv", filter_stage, fq)
        log(
            f"  {sid:<30} reads={qc.reads:>9,}  N50={qc.n50:>6}  "
            f"mode={qc.modal_bin_start}-{qc.modal_bin_start + BIN_WIDTH - 1}  "
            f"meanQ={qc.mean_read_q:>5.2f}  "
            f"in {minlength}-{maxlength}bp={qc.pct_in_window:>6.2f}%"
        )

    summary = qc_dir / "demux_qc_summary.tsv"
    write_summary_tsv(list(results.values()), summary)
    log("")
    log(f" Per-sample TSVs: {qc_dir}/<SampleID>.tsv")
    log(f" Run summary:     {summary}")
    log(f" Demux QC complete: {datetime.now():%Y-%m-%d %H:%M:%S}")
    log("========================================")
    return results


# =============================================================================
# CLI (standalone use on an existing run directory)
# =============================================================================
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--integrated-dir", type=Path, required=True,
                   help="Directory holding the merged sample_*.fastq files.")
    p.add_argument("--qc-dir", type=Path, default=None,
                   help="Output directory (default: <integrated-dir>/../demux_qc).")
    p.add_argument("--minlength", type=int, default=1000,
                   help="Lower bound of the reported length window (default: 1000).")
    p.add_argument("--maxlength", type=int, default=2000,
                   help="Upper bound of the reported length window (default: 2000).")
    p.add_argument("--filter-stage", choices=("pre", "post"), default="post",
                   help="Recorded in the output so stages stay comparable.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if not args.integrated_dir.is_dir():
        sys.exit(f"ERROR: integrated dir not found: {args.integrated_dir}")
    qc_dir = args.qc_dir or (args.integrated_dir.parent / "demux_qc")
    run_demux_qc(
        integrated_dir=args.integrated_dir,
        qc_dir=qc_dir,
        minlength=args.minlength,
        maxlength=args.maxlength,
        filter_stage=args.filter_stage,
    )


if __name__ == "__main__":
    main()
