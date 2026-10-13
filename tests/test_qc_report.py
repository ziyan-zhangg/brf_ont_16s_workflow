"""Fixture test for qc_report.py: slice profiling, counts parser, flagging.

Builds a small fake run (raw fastq_pass + read_counts_summary.txt written by
the production writer), runs the report end to end with fake Minibar and
Chopper, then checks the categories, the flags, the parser's failure modes,
run-level relative mode, and that --report-only rebuilds the same table.

Needs PyYAML (the config is YAML). matplotlib is optional; without it the
plot checks are skipped.

Run:  python tests/test_qc_report.py
"""
import contextlib
import gzip
import importlib.util
import io
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
HERE = Path(__file__).resolve().parent
FAKE_MINIBAR = HERE / "fake_minibar.py"
FAKE_CHOPPER = HERE / "fake_chopper.py"

sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location("qc_report", ROOT / "qc_report.py")
qc = importlib.util.module_from_spec(spec)
sys.modules["qc_report"] = qc
spec.loader.exec_module(qc)
core = qc.core
from cutadapt_2pass import CutadaptStats  # noqa: E402

# The real chopper is a binary; our stand-in is a .py, so it needs an
# interpreter in front of it. Test-harness only.
_real_popen = subprocess.Popen


def _popen(cmd, **kw):
    if cmd and str(cmd[0]).endswith(".py"):
        cmd = [sys.executable, *map(str, cmd)]
    return _real_popen(cmd, **kw)


subprocess.Popen = _popen

HIGH_Q = "I"    # Q40
LOW_Q = "+"     # Q10, below the Q15 Chopper floor

# sample -> per raw file: list of (n_reads, length, qual_char)
GOOD = [(30, 1500, HIGH_Q)]
DESIGN = {
    **{f"S{i:02d}": GOOD for i in range(1, 11)},
    "BP_short": [(15, 1500, HIGH_Q), (15, 500, HIGH_Q)],      # 50% short
    "BP_lowq": [(12, 1500, HIGH_Q), (18, 1500, LOW_Q)],       # 60% low Q
    "BP_long": [(27, 1500, HIGH_Q), (3, 2600, HIGH_Q)],       # 10% long
    "NTC1": [(2, 800, HIGH_Q)],
    "Tiny": [(1, 1500, HIGH_Q)],                              # n<200 in slice
    "unk": [(5, 300, LOW_Q)],
}
N_FILES = 20
EVERY_NTH = 2


def write_raw(raw: Path) -> dict:
    """Write N_FILES raw chunks; return reads per sample in the slice."""
    raw.mkdir(parents=True)
    picked = set(range(0, N_FILES, EVERY_NTH))
    in_slice: dict = {}
    for f in range(N_FILES):
        lines = []
        k = 0
        for sid, groups in DESIGN.items():
            for n, length, q in groups:
                for _ in range(n):
                    k += 1
                    tag = "" if sid == "unk" else f" bc={sid}"
                    lines.append(f"@r{f}_{k}{tag}\n{'A' * length}\n+\n{q * length}\n")
                    if f in picked:
                        in_slice[sid] = in_slice.get(sid, 0) + 1
        with gzip.open(raw / f"PAQ_pass_abc_{f}.fastq.gz", "wt") as fh:
            fh.write("".join(lines))
    return in_slice


def write_counts(out: Path, delivered: dict, pass1_pct: dict,
                 filter_stage="pre", long_name=None) -> Path:
    """read_counts_summary.txt via the production writer, as Step 6 does."""
    integ = out / "integrated_demultiplexing"
    for sid, n in delivered.items():
        d = integ / "ClientA"
        d.mkdir(parents=True, exist_ok=True)
        with gzip.open(d / f"{sid}.fastq.gz", "wt") as fh:
            fh.write("@r\nA\n+\nI\n" * n)
    (integ / "sample_unk.fastq").write_text("@r\nA\n+\nI\n" * 7)
    if long_name:
        with gzip.open(integ / "ClientA" / f"{long_name}.fastq.gz", "wt") as fh:
            fh.write("@r\nA\n+\nI\n" * 12345)
    cut = {sid: CutadaptStats(sid, input_reads=1000,
                              pass1_kept=int(10 * pass1_pct[sid]))
           for sid in pass1_pct}
    qc_by = {sid: type("Q", (), {"reads": n + 10})() for sid, n in delivered.items()}
    log_dir = out / "run_log_20261013"
    log_dir.mkdir(parents=True, exist_ok=True)
    path = log_dir / "read_counts_summary.txt"
    with contextlib.redirect_stdout(io.StringIO()):
        core.write_run_summary(integ, path, out / "raw", None, out,
                               total_input=10**6, filter_stage=filter_stage,
                               qc_by_sid=qc_by, cutadapt_by_sid=cut)
    return path


def expect_exit(fn, *args) -> str:
    try:
        fn(*args)
    except SystemExit as e:
        return str(e)
    return ""


def table(path: Path) -> dict:
    rows = [ln.split("\t") for ln in path.read_text().splitlines()
            if not ln.startswith("#")]
    head = rows[0]
    return {r[1]: dict(zip(head, r)) for r in rows[1:]}


def main() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="qc_report_test_"))
    failures = []

    def check(ok, msg_pass, msg_fail):
        if ok:
            print(f"PASS  {msg_pass}")
        else:
            failures.append(msg_fail)

    # ---- fake run ----------------------------------------------------------
    raw = tmp / "fastq_pass"
    in_slice = write_raw(raw)
    out = tmp / "run"
    samples = [s for s in DESIGN if s != "unk"]
    delivered = {s: 20000 for s in samples}
    delivered.update({"NTC1": 40, "Tiny": 120})
    pass1 = {s: 97.0 for s in samples}
    counts_path = write_counts(out, delivered, pass1)

    sheet = tmp / "sheet.csv"
    sheet.write_text("Client,Sample_ID,Location,Barcode,Sample_Type\n" + "".join(
        f"ClientA,{s},A{i + 1},TBC{i + 1:03d},{'NTC' if s == 'NTC1' else ''}\n"
        for i, s in enumerate(samples)))
    (tmp / "primer.txt").write_text("unused by fake_minibar\n")

    cfg_path = tmp / "qc_config.yaml"
    cfg_path.write_text((ROOT / "qc_config.yaml").read_text().replace(
        "every_nth_file: 10", f"every_nth_file: {EVERY_NTH}"))

    argv = ["--raw-input-dir", str(raw), "--output-dir", str(out),
            "--primer-file", str(tmp / "primer.txt"),
            "--samplesheet", str(sheet), "--config", str(cfg_path),
            "--minibar", str(FAKE_MINIBAR), "--chopper", str(FAKE_CHOPPER),
            "--work-dir", str(tmp / "work"), "--workers", "2"]
    with contextlib.redirect_stdout(io.StringIO()):
        rc = qc.main(argv)
    qc_dir = out / "qc_report"
    check(rc == 0, "end-to-end run completed", f"qc.main returned {rc}")

    # ---- 1. categories sum to slice reads, fractions to 100% ---------------
    profiles, meta = qc.read_bins(qc_dir / "qc_bins.tsv.gz")
    bad = {s: (profiles[s].n, in_slice[s]) for s in in_slice
           if s not in profiles or profiles[s].n != in_slice[s]}
    check(not bad, f"category counts == slice reads for {len(in_slice)} bins",
          f"category counts != slice reads: {bad}")
    off = [s for s, p in profiles.items()
           if abs(sum(p.frac(c) for c in qc.CATEGORIES) - 1) > 1e-9]
    check(not off, "short+long+lowQ+pass = 100% for every bin",
          f"fractions do not sum to 1: {off}")
    len_ok = all(sum(p.len_bins.values()) == p.n for p in profiles.values())
    check(len_ok, "length histogram totals == slice reads",
          "length histogram totals differ from slice reads")
    check(meta["slice_files"] == N_FILES // EVERY_NTH
          and meta["total_files"] == N_FILES,
          f"slice = {meta['slice_files']} of {meta['total_files']} files",
          f"unexpected slice meta: {meta}")

    # ---- 2. Chopper cross-check agrees ------------------------------------
    run_qc = (qc_dir / "run_qc.txt").read_text()
    check("-> OK" in run_qc and meta["chopper_in"] == sum(in_slice.values()),
          "Chopper cross-check OK and saw every slice read",
          f"cross-check not OK:\n{run_qc}")

    # ---- 3. flags ----------------------------------------------------------
    t = table(qc_dir / "per_sample_qc.tsv")
    want = {"S01": "PASS", "BP_short": "FAIL", "BP_lowq": "FAIL",
            "BP_long": "CHECK", "NTC1": "CTRL", "Tiny": "FAIL"}
    got = {s: t[s]["flag"] for s in want}
    check(got == want, f"flags {got}", f"flags: want {want}, got {got}")
    check(t["BP_short"]["reasons"].startswith("FAIL: SHORT(50%)"),
          f"BP_short reasons: {t['BP_short']['reasons']}",
          f"BP_short reasons: {t['BP_short']['reasons']}")
    check("LOWQ(60%)" in t["BP_lowq"]["reasons"]
          and "RETENTION(40%)" in t["BP_lowq"]["reasons"],
          f"BP_lowq reasons: {t['BP_lowq']['reasons']}",
          f"BP_lowq reasons: {t['BP_lowq']['reasons']}")
    check(t["Tiny"]["retention"] == "NA (n<200)"
          and t["Tiny"]["reasons"] == "FAIL: READS(120)",
          "n<200: estimates NA, exact metrics still flagged",
          f"Tiny row: {t['Tiny']}")
    check(t["S01"]["primer_found"] == "0.9700"
          and t["S01"]["delivered_reads"] == "20000",
          "exact metrics parsed from read_counts_summary.txt",
          f"S01 row: {t['S01']}")
    check((qc_dir / "per_client" / "ClientA_qc.tsv").is_file(),
          "per-client table written", "per_client/ClientA_qc.tsv missing")

    try:
        import matplotlib  # noqa: F401
        pngs = {p.name for p in qc_dir.glob("*.png")}
        flagged = {p.stem for p in (qc_dir / "flagged").glob("*.png")}
        check(pngs == {"length_heatmap.png", "q_heatmap.png"}
              and flagged == {"BP_short", "BP_lowq", "BP_long", "Tiny"},
              f"plots: {sorted(pngs)} + flagged {sorted(flagged)}",
              f"plots: {pngs}, flagged {flagged}")
    except ImportError:
        print("SKIP  plots (matplotlib not installed)")

    # ---- 4. --report-only rebuilds the same table --------------------------
    before = [ln for ln in (qc_dir / "per_sample_qc.tsv").read_text().splitlines()
              if not ln.startswith("# generated")]
    with contextlib.redirect_stdout(io.StringIO()):
        qc.main(["--report-only", "--output-dir", str(out),
                 "--samplesheet", str(sheet), "--config", str(cfg_path)])
    after = [ln for ln in (qc_dir / "per_sample_qc.tsv").read_text().splitlines()
             if not ln.startswith("# generated")]
    check(before == after, "--report-only rebuilds an identical table",
          "--report-only changed per_sample_qc.tsv")

    # ---- 5. counts parser --------------------------------------------------
    rows = qc.parse_read_counts(counts_path)
    check(rows["S01"].reads == 20000 and rows["S01"].pct_pass1_kept == 97.0
          and rows["S01"].post_chopper == 20010 and "unk" in rows,
          "parser: pre-mode table with blank raw_demuxed/pct_retained",
          f"parser pre-mode: {rows.get('S01')}")
    check("TOTAL DEMULTIPLEXED" not in rows and len(rows) == len(samples) + 1,
          "parser: TOTAL rows ignored", f"parser rows: {sorted(rows)}")

    long_name = "X" * 55
    out2 = tmp / "run_long"
    p2 = write_counts(out2, {"S01": 20000}, {"S01": 88.5, long_name: 91.2},
                      long_name=long_name)
    r2 = qc.parse_read_counts(p2)
    check(r2[long_name].reads == 12345 and r2[long_name].pct_pass1_kept == 91.2
          and r2["S01"].pct_pass1_kept == 88.5,
          "parser: over-long sample name shifts columns correctly",
          f"parser long name: {r2.get(long_name)}")

    out3 = tmp / "run_post"
    p3 = write_counts(out3, {"S01": 20000}, {"S01": 50.0}, filter_stage="post")
    r3 = qc.parse_read_counts(p3)
    check(r3["S01"].pct_pass1_kept == 50.0, "parser: post-mode table",
          f"parser post-mode: {r3.get('S01')}")

    broken = tmp / "broken.txt"
    broken.write_text(counts_path.read_text().replace("pct_pass1_kept", "pct_kept"))
    msg = expect_exit(qc.parse_read_counts, broken)
    check("pct_pass1_kept" in msg, "parser: missing column fails loudly by name",
          f"parser missing column message: {msg!r}")

    old = tmp / "old.txt"
    old.write_text("Sample                                        Reads % of Input\n"
                   "----------------------------------------\n"
                   "sample_S01                                    45249      4.00%\n")
    msg = expect_exit(qc.parse_read_counts, old)
    check("missing column" in msg, "parser: legacy 3-column table rejected",
          f"parser legacy message: {msg!r}")

    # ---- 6. run-level relative mode ----------------------------------------
    cfg = qc.load_config(cfg_path)
    sheet_rows = [qc.SheetSample("C", f"R{i:02d}", order=i) for i in range(12)]
    counts = {s.sample_id: qc.CountsRow(s.sample_id, 30000, None, None, 97.0)
              for s in sheet_rows}
    profs = {}
    for i, s in enumerate(sheet_rows):
        lowq = 80 if i == 0 else 38 + i % 5      # run median ~40% lowQ
        p = qc.SliceProfile(s.sample_id)
        p.cats.update(short=0, long=0, lowq=lowq * 10, **{"pass": (100 - lowq) * 10})
        profs[s.sample_id] = p
    warns = []
    res, runs = qc.evaluate(sheet_rows, counts, profs, cfg, warns)
    lowq_run = next(m for m in runs if m.metric == "lowq_frac")
    flags = {r.sample.sample_id: r.flag for r in res}
    check(lowq_run.run_level == 1 and lowq_run.relative
          and any("RUN: LOW_Q (median 40%" in ln for ln in qc.run_flag_lines(runs, cfg)),
          "run flag RUN: LOW_Q raised and metric switched to relative mode",
          f"run-level lowq: {lowq_run}")
    check(flags["R00"] == "FAIL"
          and all(f == "PASS" for s, f in flags.items() if s != "R00"),
          "relative mode flags only the outlier, not the whole run",
          f"relative flags: {flags}")
    check("z=" in res[0].reasons, f"outlier reason: {res[0].reasons}",
          f"outlier reason lacks z: {res[0].reasons}")

    # delivered_reads stays absolute even when the whole run is low
    for c in counts.values():
        c.reads = 9000
    res, runs = qc.evaluate(sheet_rows, counts, profs, cfg, [])
    rd = next(m for m in runs if m.metric == "delivered_reads")
    check(rd.run_level == 1 and not rd.relative
          and all(r.levels["delivered_reads"] == 1 for r in res),
          "delivered_reads: run flag raised but samples stay absolute",
          f"delivered_reads run: {rd}")

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
