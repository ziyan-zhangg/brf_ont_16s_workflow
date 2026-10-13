"""Stand-in for tools/chopper/chopper-linux-musl, for tests only.

Reads fastq on stdin, writes fastq on stdout, dropping records outside the
--minlength/--maxlength window. Enough to exercise run_chopper()'s plumbing;
the real chopper also filters on mean quality.
"""
import argparse
import sys


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
        if args.minlength <= len(seq.strip()) <= args.maxlength:
            out.write(hdr + seq + plus + qual)
            kept += 1
    out.flush()
    print(f"Kept {kept} reads out of {total} reads", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
