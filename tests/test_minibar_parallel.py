"""Fixture test for Step 2: parallel Minibar, zero-read chunks, failures.

Checks that the parallel chunk loop produces byte-identical merged output to
a serial one, that empty chunks are skipped rather than crashing Minibar,
and that a malformed non-empty chunk is still fatal.

Run:  python tests/test_minibar_parallel.py
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
FAKE_MINIBAR = Path(__file__).resolve().parent / "fake_minibar.py"


def load_core():
    spec = importlib.util.spec_from_file_location(
        "core", ROOT / "filter_chopper_demux_minibar.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["core"] = mod
    spec.loader.exec_module(mod)
    return mod


core = load_core()


def write_chunk(path: Path, reads, gz=True):
    """reads = [(read_index, sample_id)]"""
    body = "".join(
        f"@read{i} bc={sid}\nACGTACGTAC\n+\nIIIIIIIIII\n" for i, sid in reads)
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(body)
    else:
        path.write_text(body)


def md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


def merged_digests(integrated: Path) -> dict:
    return {p.name: md5(p) for p in sorted(integrated.glob("sample_*.fastq"))}


def demux(tmp: Path, inputs, workers, tag):
    """Skip empty chunks, run Step 2, merge -- the main() path for Step 2."""
    os.environ["PBS_NCPUS"] = str(workers)
    work = tmp / f"work_{tag}"
    logs = tmp / f"logs_{tag}"
    integrated = tmp / f"integrated_{tag}"

    non_empty = [f for f in inputs if core.fastq_has_reads(f)]
    skipped = [f.name for f in inputs if f not in non_empty]
    dirs = core.run_minibar(non_empty, work, tmp / "primer.txt",
                            FAKE_MINIBAR, logs)
    core.merge_per_file_outputs(dirs, integrated)
    return integrated, skipped, logs


def main() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="minibar_test_"))
    (tmp / "primer.txt").write_text("unused by fake_minibar\n")
    failures = []

    # Three normal chunks, interleaved samples so merge order is observable,
    # plus one chunk with zero reads.
    write_chunk(tmp / "c1.fastq.gz", [(1, "S1"), (2, "S2")])
    write_chunk(tmp / "c2.fastq.gz", [(3, "S1"), (4, "S3")])
    write_chunk(tmp / "c3.fastq.gz", [(5, "S2"), (6, "S1")])
    write_chunk(tmp / "c0_empty.fastq.gz", [])
    inputs = sorted(tmp.glob("c*.fastq.gz"))

    # --- 1. serial vs parallel must be byte-identical ---
    ser, ser_skip, _ = demux(tmp, inputs, 1, "serial")
    par, par_skip, par_logs = demux(tmp, inputs, 2, "parallel")
    if merged_digests(ser) != merged_digests(par):
        failures.append(f"merged output differs:\n  serial   "
                        f"{merged_digests(ser)}\n  parallel {merged_digests(par)}")
    else:
        print(f"PASS  serial == parallel: {merged_digests(ser)}")

    # --- 2. the empty chunk was skipped, not fed to minibar ---
    if par_skip != ["c0_empty.fastq.gz"]:
        failures.append(f"empty chunk not skipped, got {par_skip}")
    else:
        print(f"PASS  skipped {par_skip}")

    # --- 3. per-chunk minibar logs were copied off the work dir ---
    got = sorted(p.name for p in par_logs.glob("*.log"))
    want = ["c1.log", "c2.log", "c3.log"]
    if got != want:
        failures.append(f"minibar logs: want {want}, got {got}")
    else:
        print(f"PASS  minibar logs saved: {got}")

    # --- 4. plain .fastq inputs work too (Change 2 writes these) ---
    write_chunk(tmp / "p1.fastq", [(1, "S1"), (2, "S2")], gz=False)
    write_chunk(tmp / "p2.fastq", [(3, "S1"), (4, "S3")], gz=False)
    write_chunk(tmp / "p0_empty.fastq", [], gz=False)
    plain = sorted(tmp.glob("p*.fastq"))
    pl, pl_skip, _ = demux(tmp, plain, 2, "plain")
    if pl_skip != ["p0_empty.fastq"]:
        failures.append(f"empty plain chunk not skipped, got {pl_skip}")
    elif not (pl / "sample_S1.fastq").is_file():
        failures.append("plain fastq produced no merged sample_S1.fastq")
    else:
        print(f"PASS  plain .fastq handled, skipped {pl_skip}")

    # --- 5. a malformed NON-empty chunk must still be fatal ---
    bad = tmp / "bad.fastq.gz"
    with gzip.open(bad, "wt") as fh:
        fh.write("this is not a fastq at all\n")
    rc = subprocess.call(
        [sys.executable, "-c",
         "import sys; sys.argv=['x']; "
         f"sys.path.insert(0, r'{ROOT}'); "
         "import importlib.util, pathlib; "
         f"spec=importlib.util.spec_from_file_location('c', r'{ROOT / 'filter_chopper_demux_minibar.py'}'); "
         "m=importlib.util.module_from_spec(spec); sys.modules['c']=m; "
         "spec.loader.exec_module(m); "
         f"m.run_minibar([pathlib.Path(r'{bad}')], pathlib.Path(r'{tmp / 'work_bad'}'), "
         f"pathlib.Path(r'{tmp / 'primer.txt'}'), pathlib.Path(r'{FAKE_MINIBAR}'), "
         f"pathlib.Path(r'{tmp / 'logs_bad'}'))"],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if rc == 0:
        failures.append("malformed non-empty chunk did NOT fail the run")
    else:
        print(f"PASS  malformed non-empty chunk is fatal (exit {rc})")

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
