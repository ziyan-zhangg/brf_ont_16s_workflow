"""Stand-in for tools/minibar/minibar.py, for tests only.

Mimics the parts the pipeline depends on: reads one fastq (gzipped or plain),
bins each read by a sample id derived from its header, and writes
sample_<id>.fastq into the current working directory. Crashes the same way
the real minibar does on an empty input -- indexing the first line -- so the
zero-read skip is exercised against real behaviour rather than a guess.
"""
import gzip
import sys
from pathlib import Path


def main() -> int:
    path = Path(sys.argv[-1])
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        firstln = fh.readline()
        # The real minibar does exactly this, hence the IndexError on empty.
        if firstln[0] != "@":
            sys.exit("not a fastq")
        fh.seek(0)
        bins: dict[str, list[str]] = {}
        while True:
            hdr = fh.readline()
            if not hdr:
                break
            rec = hdr + fh.readline() + fh.readline() + fh.readline()
            # header looks like "@read<N> bc=<id>"
            sid = hdr.split("bc=")[1].strip() if "bc=" in hdr else "unk"
            bins.setdefault(sid, []).append(rec)
    for sid, recs in bins.items():
        Path(f"sample_{sid}.fastq").write_text("".join(recs))
    print(f"fake_minibar: {path.name} -> {len(bins)} bin(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
