"""Stand-in for tools/chopper/chopper-linux-musl, for tests only.

Reads fastq on stdin, writes fastq on stdout, dropping records outside the
--minlength/--maxlength window or with mean read Q below -q. Mean Q is the
mean per-base error probability converted back to Phred, as the real chopper
computes it, so qc_report.py's Chopper cross-check can be exercised.
"""
import argparse
import math
import sys


def mean_q(qual: bytes) -> float:
    if not qual:
        return 0.0
    err = sum(10.0 ** (-(b - 33) / 10.0) for b in qual) / len(qual)
    return -10.0 * math.log10(err) if err > 0 else 93.0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-q", type=int, default=0)
    ap.add_argument("--minlength", type=int, default=0)
    ap.add_argument("--maxlength", type=int, default=10**9)
    args = ap.parse_args()

    inp, out = sys.stdin.buffer, sys.stdout.buffer
    kept = total = 0
    while True:
        hdr = inp.readline()
        if not hdr:
            break
        seq, plus, qual = inp.readline(), inp.readline(), inp.readline()
        total += 1
        if (args.minlength <= len(seq.strip()) <= args.maxlength
                and mean_q(qual.strip()) >= args.q):
            out.write(hdr + seq + plus + qual)
            kept += 1
    out.flush()
    print(f"Kept {kept} reads out of {total} reads", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
