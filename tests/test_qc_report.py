"""Fixture test for qc_report.py (v2) and the demux_qc.py Q histogram.

Builds a small fake run: per-sample demultiplexed FASTQs profiled by the real
demux_qc.run_demux_qc, and read_counts_summary.txt written by the production
writer, in post and in pre mode. Then runs the report end to end and checks
the Q bins, the metrics, the flags, the parser (incl. pre-mode NA handling),
run-level relative mode and controls.

Needs PyYAML (the config is YAML). matplotlib is optional; without it the
plot checks are skipped.

Run:  python tests/test_qc_report.py
"""
import contextlib
import gzip
import importlib.util
import io
import math
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

sys.path.insert(0, str(ROOT))
spec = importlib.util.spec_from_file_location("qc_report", ROOT / "qc_report.py")
qc = importlib.util.module_from_spec(spec)
sys.modules["qc_report"] = qc
spec.loader.exec_module(qc)
core = qc.core
from cutadapt_2pass import CutadaptStats  # noqa: E402
from demux_qc import run_demux_qc  # noqa: E402

Q40, Q20, Q17, Q10 = "I", "5", "2", "+"


def truncated(n):
    """n partial amplicons spread over 1,050-1,225 bp, so none of them
    becomes the modal length."""
    return [(1, 1050 + 25 * (i % 8), Q40) for i in range(n)]


# sample -> list of (n_reads, length, qual_char)
DESIGN = {
    **{f"S{i:02d}": [(300, 1500, Q40)] for i in range(1, 10)},
    "S10": [(300, 1500, Q20)],                         # exactly Q20: counts
    "BP_core": [(140, 1500, Q40), *truncated(160)],    # core 47%
    "BP_long": [(270, 1500, Q40), (30, 2600, Q40)],    # 10% > 2 kb
    "BP_q": [(150, 1500, Q40), (150, 1500, Q17)],      # 50% Q20+
    "ITS": [(300, 1350, Q40)],                         # modal 150 bp off
    "BP_primer": [(300, 1500, Q40)],                   # 70% primer found
    "NTC1": [(20, 1500, Q40)],
    "Tiny": [(100, 1500, Q40)],                        # n < 200
    "unk": [(50, 1500, Q10)],
}
SAMPLES = [s for s in DESIGN if s != "unk"]
PASS1 = {s: 97.0 for s in SAMPLES} | {"BP_primer": 70.0}
DELIVERED = {s: 20000 for s in SAMPLES} | {"NTC1": 40, "Tiny": 100}


def mean_q(qual: str) -> float:
    err = math.fsum(10 ** (-(ord(c) - 33) / 10) for c in qual) / len(qual)
    return -10 * math.log10(err)


def reads_of(sid):
    for n, length, q in DESIGN[sid]:
        for _ in range(n):
            yield length, q


def write_demux_fastqs(d: Path) -> None:
    d.mkdir(parents=True)
    for sid in DESIGN:
        with (d / f"sample_{sid}.fastq").open("w") as fh:
            for k, (length, q) in enumerate(reads_of(sid)):
                fh.write(f"@r{k}\n{'A' * length}\n+\n{q * length}\n")


def chopper_counts(sid):
    """(raw, kept) as production Chopper (-q 15, 1000-2000 bp) would."""
    raw = kept = 0
    for length, q in reads_of(sid):
        raw += 1
        kept += 1000 <= length <= 2000 and mean_q(q * 10) >= 15
    return raw, kept


def build_run(tmp: Path, name: str, stage: str, long_name=None) -> Path:
    """Run dir with demux_qc/ and run_log_*/read_counts_summary.txt."""
    out = tmp / name
    raw_demux = tmp / f"{name}_demux"
    write_demux_fastqs(raw_demux)
    with contextlib.redirect_stdout(io.StringIO()):
        qc_by = run_demux_qc(raw_demux, out / "demux_qc", 1000, 2000,
                             filter_stage=stage)

    integ = out / "integrated_demultiplexing"
    (integ / "ClientA").mkdir(parents=True)
    delivered = dict(DELIVERED)
    if long_name:
        delivered[long_name] = 12345
    for sid, n in delivered.items():
        with gzip.open(integ / "ClientA" / f"{sid}.fastq.gz", "wt") as fh:
            fh.write("@r\nA\n+\nI\n" * n)
    (integ / "sample_unk.fastq").write_text("@r\nA\n+\nI\n" * 7)

    pass1 = dict(PASS1)
    if long_name:
        pass1[long_name] = 91.2
    cut = {sid: CutadaptStats(sid, input_reads=1000,
                              pass1_kept=int(10 * pct))
           for sid, pct in pass1.items()}
    chop = ({sid: chopper_counts(sid) for sid in SAMPLES}
            if stage == "post" else None)
    log_dir = out / "run_log_20261013"
    log_dir.mkdir(parents=True)
    with contextlib.redirect_stdout(io.StringIO()):
        core.write_run_summary(integ, log_dir / "read_counts_summary.txt",
                               out / "raw", None, out, total_input=10**6,
                               filter_stage=stage, qc_by_sid=qc_by,
                               chopper_by_sid=chop, cutadapt_by_sid=cut)
    return out


def write_sheet(path: Path) -> None:
    path.write_text("Client,Sample_ID,Location,Barcode,Sample_Type\n" + "".join(
        f"ClientA,{s},A{i + 1},TBC{i + 1:03d},{'NTC' if s == 'NTC1' else ''}\n"
        for i, s in enumerate(SAMPLES)))


def expect_exit(fn, *args) -> str:
    try:
        fn(*args)
    except SystemExit as e:
        return str(e)
    return ""


def table(path: Path, key: int = 1) -> dict:
    rows = [ln.split("\t") for ln in path.read_text().splitlines()
            if not ln.startswith("#")]
    head = rows[0]
    return {r[key]: dict(zip(head, r)) for r in rows[1:]}


def run_report(out: Path, sheet: Path, *extra) -> int:
    with contextlib.redirect_stdout(io.StringIO()):
        return qc.main(["--output-dir", str(out), "--samplesheet", str(sheet),
                        "--config", str(ROOT / "qc_config.yaml"), *extra])


def main() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="qc_report_test_"))
    failures = []

    def check(ok, msg_pass, msg_fail):
        if ok:
            print(f"PASS  {msg_pass}")
        else:
            failures.append(msg_fail)

    sheet = tmp / "sheet.csv"
    write_sheet(sheet)
    cfg = qc.load_config(ROOT / "qc_config.yaml")

    # ---- 1. demux_qc Q bins ------------------------------------------------
    out = build_run(tmp, "run_post", "post")
    dq = out / "demux_qc"
    bad = []
    for sid in DESIGN:
        p = qc.read_demux_profile(dq, sid)
        n = sum(1 for _ in reads_of(sid))
        n_win = sum(1 for L, _ in reads_of(sid) if 1000 <= L <= 2000)
        hi = sum(1 for L, q in reads_of(sid)
                 if 1000 <= L <= 2000 and mean_q(q * L) >= 20 - 1e-9)
        direct = hi / n_win if n_win else None
        if (sum(p.q_all.values()) != n or sum(p.q_win.values()) != n_win
                or qc.pct_q20(p) != direct):
            bad.append((sid, sum(p.q_all.values()), n, qc.pct_q20(p), direct))
    check(not bad, f"Q bins sum to reads and pct_q20 matches a direct "
                   f"calculation for {len(DESIGN)} bins",
          f"Q bins / pct_q20 mismatch: {bad}")
    summ = table(dq / "demux_qc_summary.tsv", key=0)
    check(summ["S10"]["pct_q20"] == "100.00" and summ["BP_q"]["pct_q20"] == "50.00",
          "demux_qc_summary.tsv pct_q20 column (uniform Q20 counts as >= 20)",
          f"summary pct_q20: S10={summ['S10'].get('pct_q20')} "
          f"BP_q={summ['BP_q'].get('pct_q20')}")

    # ---- 2. post-mode end to end -------------------------------------------
    rc = run_report(out, sheet)
    qd = out / "qc_report"
    check(rc == 0, "post-mode end-to-end run completed", f"qc.main returned {rc}")
    t = table(qd / "per_sample_qc.tsv")
    want = {"S01": "PASS", "S10": "PASS", "BP_core": "FAIL", "BP_long": "CHECK",
            "BP_q": "CHECK", "ITS": "FAIL", "BP_primer": "FAIL",
            "NTC1": "CTRL", "Tiny": "FAIL"}
    got = {s: t[s]["flag"] for s in want}
    check(got == want, f"flags {got}", f"flags: want {want}, got {got}")
    expect_reasons = {
        "BP_core": "FAIL: LOW_CORE(47%)",
        "BP_long": "CHECK: LONG(10%)",
        "BP_q": "CHECK: LOW_Q20(50%)",
        "ITS": "FAIL: MODAL(150bp)",
        "BP_primer": "FAIL: PRIMER(70%)",
        "Tiny": "FAIL: READS(100)",
    }
    wrong = {s: t[s]["reasons"] for s, r in expect_reasons.items()
             if t[s]["reasons"] != r}
    check(not wrong, "reason codes as designed", f"reasons: {wrong}")
    s01 = t["S01"]
    check(s01["raw_demuxed"] == "300" and s01["delivered_reads"] == "20000"
          and s01["primer_found"] == "0.9700" and s01["long_frac"] == "0.0000"
          and s01["core_frac"] == "1.0000" and s01["modal_length"] == "1500"
          and s01["modal_dev_bp"] == "0" and s01["pct_q20"] == "1.0000",
          "S01 metrics exact", f"S01 row: {s01}")
    check(t["BP_long"]["core_frac"] == "1.0000",
          "core_frac counts only reads inside 1000-2000 bp",
          f"BP_long core_frac: {t['BP_long']['core_frac']}")
    check(t["Tiny"]["core_frac"] == "NA (n<200)"
          and t["Tiny"]["primer_found"] == "NA (n<200)",
          "n<200: stats NA, delivered_reads still flagged",
          f"Tiny row: {t['Tiny']}")
    check(t["NTC1"]["reasons"] == "" and t["NTC1"]["flag"] == "CTRL",
          "control flagged CTRL with no reasons", f"NTC1 row: {t['NTC1']}")
    head = (qd / "per_sample_qc.tsv").read_text()
    check("# config_version\t2.0" in head and "# filter_stage\tpost" in head,
          "TSV header records config version and filter stage",
          "TSV header lacks config version / filter stage")
    check((qd / "per_client" / "ClientA_qc.tsv").is_file(),
          "per-client table written", "per_client/ClientA_qc.tsv missing")
    with gzip.open(qd / "qc_bins.tsv.gz", "rt") as fh:
        kinds = {ln.split("\t")[1] for ln in fh.read().splitlines()[1:]}
    check(kinds == {"len", "q", "q_all"}, f"qc_bins.tsv.gz kinds {sorted(kinds)}",
          f"qc_bins.tsv.gz kinds: {kinds}")
    run_qc = (qd / "run_qc.txt").read_text()
    check("Filter stage: post" in run_qc and "unk" in run_qc
          and "NOTE: pre-mode" not in run_qc,
          "run_qc.txt: post stage, catch-all bins",
          f"run_qc.txt:\n{run_qc}")

    try:
        import matplotlib  # noqa: F401
        pngs = {p.name for p in qd.glob("*.png")}
        flagged = {p.stem for p in (qd / "flagged").glob("*.png")}
        want_f = {s for s, f in want.items() if f in ("CHECK", "FAIL")}
        check(pngs == {"length_heatmap.png", "q_heatmap.png"} and flagged == want_f,
              f"plots: {sorted(pngs)} + flagged {sorted(flagged)}",
              f"plots: {pngs}, flagged {flagged}")
    except ImportError:
        print("SKIP  plots (matplotlib not installed)")

    # ---- 3. pre mode -------------------------------------------------------
    out_pre = build_run(tmp, "run_pre", "pre")
    c = qc.parse_read_counts(out_pre / "run_log_20261013" / "read_counts_summary.txt")
    check(c.filter_stage == "pre" and c.rows["S01"].raw_demuxed is None
          and c.rows["S01"].pct_retained is None
          and c.rows["S01"].post_chopper == 300,
          "parser: pre-mode table with blank raw_demuxed / pct_retained",
          f"parser pre-mode: {c.filter_stage} {c.rows.get('S01')}")
    rc = run_report(out_pre, sheet)
    tp = table(out_pre / "qc_report" / "per_sample_qc.tsv")
    rq = (out_pre / "qc_report" / "run_qc.txt").read_text()
    check(rc == 0 and tp["S01"]["long_frac"] == "NA (pre mode)"
          and tp["BP_long"]["long_frac"] == "NA (pre mode)"
          and tp["S01"]["core_frac"] == "1.0000" and "NOTE: pre-mode run" in rq,
          "pre mode: long_frac NA (pre mode), noted at top of run_qc.txt",
          f"pre-mode S01: {tp.get('S01')}")

    # ---- 4. counts parser --------------------------------------------------
    cp = qc.parse_read_counts(out / "run_log_20261013" / "read_counts_summary.txt")
    r = cp.rows["BP_long"]
    check(cp.filter_stage == "post" and r.raw_demuxed == 300
          and r.post_chopper == 270 and r.pct_retained == 90.0
          and r.pct_pass1_kept == 97.0,
          "parser: post-mode raw_demuxed / post_chopper / pct_retained",
          f"parser post-mode: {cp.filter_stage} {r}")
    check("TOTAL DEMULTIPLEXED" not in cp.rows and "unk" in cp.rows,
          "parser: TOTAL rows ignored", f"parser rows: {sorted(cp.rows)}")

    long_name = "X" * 55
    out_long = build_run(tmp, "run_long", "post", long_name=long_name)
    rl = qc.parse_read_counts(out_long / "run_log_20261013" / "read_counts_summary.txt")
    check(rl.rows[long_name].reads == 12345
          and rl.rows[long_name].pct_pass1_kept == 91.2
          and rl.rows["BP_long"].pct_retained == 90.0,
          "parser: over-long sample name shifts columns correctly",
          f"parser long name: {rl.rows.get(long_name)}")

    counts_txt = (out / "run_log_20261013" / "read_counts_summary.txt").read_text()
    broken = tmp / "broken.txt"
    broken.write_text(counts_txt.replace("pct_pass1_kept", "pct_kept"))
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

    # ---- 5. demux_qc without Q histogram (older runs) ----------------------
    for q in (out / "demux_qc").glob("*.q.tsv"):
        q.unlink()
    run_report(out, sheet)
    t2 = table(qd / "per_sample_qc.tsv")
    rq = (qd / "run_qc.txt").read_text()
    check(t2["S01"]["pct_q20"] == "NA" and t2["BP_q"]["flag"] == "PASS"
          and "pct_q20 is NA" in rq,
          "no .q.tsv: pct_q20 NA with a warning, other metrics unaffected",
          f"no-Q run: {t2['S01']} / {t2['BP_q']['flag']}")
    msg = expect_exit(run_report, out, sheet, "--demux-qc-dir", str(tmp / "nope"))
    check("demux QC dir not found" in msg, "missing demux_qc/ fails loudly",
          f"missing demux_qc message: {msg!r}")

    # ---- 6. run-level relative mode ----------------------------------------
    sheet_rows = [qc.SheetSample("C", f"R{i:02d}", order=i) for i in range(12)]

    def synthetic(q20_pct, reads=30000):
        rows = {s.sample_id: qc.CountsRow(s.sample_id, reads, 1000, 1000, 100.0,
                                          970, 97.0) for s in sheet_rows}
        profs = {}
        for i, s in enumerate(sheet_rows):
            p = qc.DemuxProfile(s.sample_id, reads=1000, filter_stage="post",
                                modal_length=1500, len_bins={1500: 1000},
                                has_q=True)
            hi = q20_pct(i) * 10
            p.q_win = {25: hi, 17: 1000 - hi}
            profs[s.sample_id] = p
        return qc.Counts(Path("x"), "post", rows), profs

    counts, profs = synthetic(lambda i: 5 if i == 0 else 33 + i % 5)
    ev = qc.evaluate(sheet_rows, counts, profs, cfg, [])
    q_run = next(m for m in ev.runs if m.metric == "pct_q20")
    flags = {r.sample.sample_id: r.flag for r in ev.results}
    check(q_run.run_level == 2 and q_run.relative
          and any(ln.startswith("RUN: LOW_Q20 (median ") and ", FAIL)" in ln
                  for ln in qc.run_flag_lines(ev.runs)),
          "run flag RUN: LOW_Q20 raised and metric switched to relative mode",
          f"run-level pct_q20: {q_run}")
    check(flags["R00"] == "FAIL"
          and all(f == "PASS" for s, f in flags.items() if s != "R00"),
          "relative mode flags only the outlier, not the whole run",
          f"relative flags: {flags}")
    check("z=" in ev.results[0].reasons, f"outlier reason: {ev.results[0].reasons}",
          f"outlier reason lacks z: {ev.results[0].reasons}")

    # delivered_reads stays absolute even when the whole run is low
    counts, profs = synthetic(lambda i: 90)
    for c in counts.rows.values():
        c.reads = 9000
    ev = qc.evaluate(sheet_rows, counts, profs, cfg, [])
    rd = next(m for m in ev.runs if m.metric == "delivered_reads")
    check(rd.run_level == 1 and not rd.relative
          and all(r.levels["delivered_reads"] == 1 for r in ev.results),
          "delivered_reads: run flag raised but samples stay absolute",
          f"delivered_reads run: {rd}")

    # controls are left out of run medians
    counts, profs = synthetic(lambda i: 90)
    ctrl = qc.SheetSample("C", "POS", sample_type="POS_CTRL", order=99)
    counts.rows["POS"] = qc.CountsRow("POS", 30000, 1000, 1000, 100.0, 970, 97.0)
    profs["POS"] = qc.DemuxProfile("POS", reads=1000, modal_length=3000,
                                   len_bins={1500: 1000}, has_q=True,
                                   q_win={17: 1000})
    ev = qc.evaluate([*sheet_rows, ctrl], counts, profs, cfg, [])
    pos = ev.results[-1]
    check(ev.modal_median == 1500 and pos.flag == "CTRL"
          and pos.values["modal_dev_bp"] == 1500,
          "controls: flag CTRL, metrics reported, excluded from run medians",
          f"control: median {ev.modal_median}, {pos}")

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
