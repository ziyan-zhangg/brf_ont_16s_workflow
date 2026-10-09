"""Fixture test for Step 1: Chopper writes plain FASTQ, and that is equivalent.

The point of Change 2 is that Minibar's result must not depend on whether its
input arrived gzipped. This builds both forms from the same Chopper output and
compares the merged per-sample digests.

Run:  python tests/test_chopper_plain.py
"""
import gzip
import hashlib
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
HERE = Path(__file__).resolve().parent
FAKE_MINIBAR = HERE / "fake_minibar.py"
FAKE_CHOPPER = HERE / "fake_chopper.py"

spec = importlib.util.spec_from_file_location(
    "core", ROOT / "filter_chopper_demux_minibar.py")
core = importlib.util.module_from_spec(spec)
sys.modules["core"] = core
spec.loader.exec_module(core)

# The real chopper is a binary; our stand-in is a .py, so it needs an
# interpreter in front of it. Test-harness only.
_real_popen = subprocess.Popen


def _popen(cmd, **kw):
    if cmd and str(cmd[0]).endswith(".py"):
        cmd = [sys.executable, *map(str, cmd)]
    return _real_popen(cmd, **kw)


subprocess.Popen = _popen


def md5(p: Path) -> str:
    return hashlib.md5(p.read_bytes()).hexdigest()


def digests(d: Path) -> dict:
    return {p.name: md5(p) for p in sorted(d.glob("sample_*.fastq"))}


def main() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="chopper_plain_"))
    raw = tmp / "fastq_pass"
    raw.mkdir()
    (tmp / "primer.txt").write_text("unused\n")
    os.environ["PBS_NCPUS"] = "2"
    failures = []

    # Two chunks of in-window reads, one chunk whose reads are all out of
    # window (so Chopper keeps 0 and emits an empty file), interleaved
    # samples so merge order is observable.
    def chunk(name, reads):
        with gzip.open(raw / name, "wt") as fh:
            for i, sid, length in reads:
                fh.write(f"@read{i} bc={sid}\n{'A' * length}\n+\n{'I' * length}\n")

    chunk("c1.fastq.gz", [(1, "S1", 1500), (2, "S2", 1500)])
    chunk("c2.fastq.gz", [(3, "S1", 1500), (4, "S3", 1500)])
    chunk("c3_allshort.fastq.gz", [(5, "S1", 10), (6, "S2", 10)])

    # --- Step 1 as it now runs: plain output ---
    filt = tmp / "chopper_filtered"
    outs = core.run_chopper(raw, filt, FAKE_CHOPPER, 15, 1000, 2000)

    if not all(o.suffix == ".fastq" for o in outs):
        failures.append(f"expected plain .fastq outputs, got {[o.name for o in outs]}")
    else:
        print(f"PASS  plain output: {[o.name for o in outs]}")

    # Chopper emits a file per chunk even when it keeps nothing.
    empty = [o.name for o in outs if not core.fastq_has_reads(o)]
    if empty != ["c3_allshort_filtered.fastq"]:
        failures.append(f"expected the all-short chunk to be empty, got {empty}")
    else:
        print(f"PASS  zero-read chunk detected on a plain file: {empty}")

    # count_filtered_total must read plain files correctly
    total = core.count_filtered_total(outs)
    if total != 4:
        failures.append(f"count_filtered_total: want 4 reads, got {total}")
    else:
        print(f"PASS  count_filtered_total on plain files: {total}")

    non_empty = [o for o in outs if core.fastq_has_reads(o)]

    # --- demux the plain files ---
    d_plain = core.run_minibar(non_empty, tmp / "w_plain", tmp / "primer.txt",
                               FAKE_MINIBAR, tmp / "l_plain")
    i_plain = tmp / "i_plain"
    core.merge_per_file_outputs(d_plain, i_plain)

    # --- same Chopper output, gzipped, as the old code wrote it ---
    gz_dir = tmp / "as_gz"
    gz_dir.mkdir()
    gz_inputs = []
    for o in non_empty:
        g = gz_dir / (o.name + ".gz")
        with o.open("rb") as fin, gzip.open(g, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        gz_inputs.append(g)
    d_gz = core.run_minibar(gz_inputs, tmp / "w_gz", tmp / "primer.txt",
                            FAKE_MINIBAR, tmp / "l_gz")
    i_gz = tmp / "i_gz"
    core.merge_per_file_outputs(d_gz, i_gz)

    if digests(i_plain) != digests(i_gz):
        failures.append(f"plain vs gzipped input changed the merged output:\n"
                        f"  plain {digests(i_plain)}\n  gz    {digests(i_gz)}")
    else:
        print(f"PASS  plain == gzipped merged output: {digests(i_plain)}")

    # --- chunk naming must agree across both forms ---
    if [d.name for d in d_plain] != [d.name for d in d_gz]:
        failures.append(f"work dir names differ: {[d.name for d in d_plain]} "
                        f"vs {[d.name for d in d_gz]}")
    else:
        print(f"PASS  chunk names match: {[d.name for d in d_plain]}")

    # --- cleanup() removes the plain filtered dir ---
    core.cleanup(d_plain, filt)
    if filt.exists():
        failures.append("cleanup() left the filtered dir behind")
    else:
        print("PASS  cleanup() removed the filtered dir")

    shutil.rmtree(tmp, ignore_errors=True)
    print()
    if failures:
        for f in failures:
            print("FAIL  " + f)
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
