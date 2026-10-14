"""Fixture test for qc_report.py (v2.2) and the demux_qc.py Q histogram.

Builds a small fake run: per-sample demultiplexed FASTQs profiled by the real
demux_qc.run_demux_qc, and read_counts_summary.txt written by the production
writer, in post and in pre mode. Then runs the report end to end and checks
the Q bins and median_q, the metrics, the flags (Q is never flagged), the
parser (incl. pre-mode NA handling), run-level relative mode, controls, the
old-config warning and the self-contained qc_report.html.

Needs PyYAML (the config is YAML). matplotlib is optional; without it the
plot checks are skipped.

Run:  python tests/test_qc_report.py
"""
import contextlib
import getpass
import gzip
import importlib.util
import io
import json
import math
import re
import shutil
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Must never appear in an HTML report: facility paths, project, username
# (the username only when long enough not to match by accident).
INTERNAL = ("/g/data", "/scratch", "/home", "vz35",
            *([getpass.getuser()] if len(getpass.getuser()) >= 4 else []))

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
    "BP_q": [(150, 1500, Q40), (150, 1500, Q17)],      # median Q17: not flagged
    "ITS": [(300, 1350, Q40)],                         # modal 150 bp off
    "BP_primer": [(300, 1500, Q40)],                   # 70% primer found
    "NTC1": [(20, 1500, Q40)],
    "Tiny": [(100, 1500, Q40)],                        # n < 200
    "POS_B": [(300, 1500, Q40)],                       # ClientB's own control
    "unk": [(50, 1500, Q10)],
}
SAMPLES = [s for s in DESIGN if s != "unk"]
# Sample sheet: ClientB owns BP_core, ITS and the control POS_B; NTC1 has no
# client, so it is a shared control; everything else is ClientA.
CLIENT = {s: "ClientA" for s in SAMPLES} | {
    "BP_core": "ClientB", "ITS": "ClientB", "POS_B": "ClientB", "NTC1": ""}
TYPE = {"NTC1": "NTC", "POS_B": "POS_CTRL"}
COMMENT = {"Tiny": "low DNA input", "ITS": "clientB private note"}
BARCODE = {s: f"TBC{i + 1:03d}" for i, s in enumerate(SAMPLES)}
PASS1 = {s: 97.0 for s in SAMPLES} | {"BP_primer": 70.0}
DELIVERED = {s: 20000 for s in SAMPLES} | {"NTC1": 40, "Tiny": 100}


def median_low(xs):
    xs = sorted(xs)
    return xs[(len(xs) - 1) // 2] if xs else None


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
    path.write_text("Client,Sample_ID,Location,Barcode,Sample_Type,Comment\n" + "".join(
        f"{CLIENT[s]},{s},A{i + 1},{BARCODE[s]},{TYPE.get(s, '')},"
        f"{COMMENT.get(s, '')}\n"
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


def run_report(out: Path, sheet: Path, *extra, config=None) -> int:
    with contextlib.redirect_stdout(io.StringIO()):
        return qc.main(["--output-dir", str(out), "--samplesheet", str(sheet),
                        "--config", str(config or ROOT / "qc_config.yaml"),
                        *extra])


def html_blob(path: Path) -> dict:
    m = re.search(r'<script id="qc-data" type="application/json">(.*?)</script>',
                  path.read_text(encoding="utf-8"), re.S)
    return json.loads(m.group(1)) if m else {}


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
        direct = median_low(int(mean_q(q * L) + 1e-6) for L, q in reads_of(sid)
                            if 1000 <= L <= 2000)
        if (sum(p.q_all.values()) != n or sum(p.q_win.values()) != n_win
                or qc.median_q(p) != direct):
            bad.append((sid, sum(p.q_all.values()), n, qc.median_q(p), direct))
    check(not bad, f"Q bins sum to reads and median_q matches a direct "
                   f"calculation for {len(DESIGN)} bins",
          f"Q bins / median_q mismatch: {bad}")
    summ = table(dq / "demux_qc_summary.tsv", key=0)
    check(summ["S10"]["median_q"] == "20" and summ["BP_q"]["median_q"] == "17"
          and "pct_q20" not in summ["S10"],
          "demux_qc_summary.tsv median_q column (uniform Q20 lands in bin 20)",
          f"summary median_q: S10={summ['S10'].get('median_q')} "
          f"BP_q={summ['BP_q'].get('median_q')}")

    # ---- 2. post-mode end to end -------------------------------------------
    rc = run_report(out, sheet)
    qd = out / "qc_report"
    check(rc == 0, "post-mode end-to-end run completed", f"qc.main returned {rc}")
    t = table(qd / "per_sample_qc.tsv")
    want = {"S01": "PASS", "S10": "PASS", "BP_core": "FAIL", "BP_long": "CHECK",
            "BP_q": "PASS", "ITS": "FAIL", "BP_primer": "FAIL",
            "NTC1": "CTRL", "Tiny": "FAIL"}
    got = {s: t[s]["flag"] for s in want}
    check(got == want, f"flags {got}", f"flags: want {want}, got {got}")
    expect_reasons = {
        "BP_core": "FAIL: LOW_CORE(47%)",
        "BP_long": "CHECK: LONG(10%)",
        "BP_q": "",
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
          and s01["modal_dev_bp"] == "0" and s01["median_q"] == "40"
          and "pct_q20" not in s01,
          "S01 metrics exact, median_q informational, no pct_q20 column",
          f"S01 row: {s01}")
    check(t["BP_q"]["median_q"] == "17" and t["Tiny"]["median_q"] == "NA (n<200)",
          "median_q reported but never flagged; NA below min reads",
          f"median_q: BP_q={t['BP_q']} Tiny={t['Tiny']['median_q']}")
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
    check("# config_version\t2.2" in head and "# filter_stage\tpost" in head,
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

    # ---- 2b. qc_report.html -----------------------------------------------
    page = (qd / "qc_report.html").read_text(encoding="utf-8")
    blob = html_blob(qd / "qc_report.html")
    ids = {s["id"] for s in blob.get("samples", [])}
    check(ids == set(SAMPLES), f"qc_report.html JSON parses, {len(ids)} samples",
          f"qc_report.html samples: {sorted(ids)}")
    check("http://" not in page and "https://" not in page
          and "<script src" not in page and "<link" not in page,
          "qc_report.html is self-contained (no network references)",
          "qc_report.html references something external")
    hs = {s["id"]: s for s in blob["samples"]}
    check(hs["BP_core"]["flag"] == "FAIL" and hs["BP_core"]["levels"]["core_frac"] == 2
          and sum(hs["S01"]["len"]) == 300 and sum(hs["S01"]["q"]) == 300
          and [m["key"] for m in blob["metrics"]] == list(qc.METRICS)
          and hs["BP_q"]["values"]["median_q"] == 17,
          "qc_report.html data: flags, per-metric levels, bins, 5 flagged metrics",
          f"qc_report.html data: {hs['BP_core']['levels']} {blob['metrics']}")
    check(page.count("__QC_DATA__") == 0 and page.count('id="qc-data"') == 1,
          "qc_report.html template placeholders filled once",
          "qc_report.html placeholder left in or duplicated")

    try:
        import matplotlib  # noqa: F401
        pngs = {p.name for p in qd.glob("*.png")}
        check(pngs == {"length_heatmap.png", "q_heatmap.png"}
              and not (qd / "flagged").exists(),
              f"plots: {sorted(pngs)}, no flagged/ by default",
              f"plots: {pngs}, flagged/ exists: {(qd / 'flagged').exists()}")
        run_report(out, sheet, "--png-drilldown")
        flagged = {p.stem for p in (qd / "flagged").glob("*.png")}
        want_f = {s for s, f in want.items() if f in ("CHECK", "FAIL")}
        check(flagged == want_f, f"--png-drilldown: flagged {sorted(flagged)}",
              f"--png-drilldown flagged: {flagged}")
    except ImportError:
        print("SKIP  plots (matplotlib not installed)")

    # HTML must not depend on matplotlib.
    saved = {k: sys.modules[k] for k in list(sys.modules)
             if k == "matplotlib" or k.startswith("matplotlib.")}
    sys.modules["matplotlib"] = None          # makes 'import matplotlib' fail
    try:
        (qd / "qc_report.html").unlink()
        run_report(out, sheet)
    finally:
        del sys.modules["matplotlib"]
        sys.modules.update(saved)
    check((qd / "qc_report.html").is_file()
          and "matplotlib not available" in (qd / "run_qc.txt").read_text(),
          "qc_report.html still written without matplotlib",
          "qc_report.html missing when matplotlib is unavailable")

    # ---- 2c. old config with pct_q20 ---------------------------------------
    old_cfg = tmp / "old_config.yaml"
    old_cfg.write_text((ROOT / "qc_config.yaml").read_text().replace(
        "metrics:\n", "metrics:\n  pct_q20: {direction: higher_better, "
        "check: 0.60, fail: 0.40}\n", 1))
    rc = run_report(out, sheet, config=old_cfg)
    rq = (qd / "run_qc.txt").read_text()
    t3 = table(qd / "per_sample_qc.tsv")
    check(rc == 0 and "metrics.pct_q20 ignored: Q is no longer flagged" in rq
          and t3["BP_q"]["flag"] == "PASS",
          "old config with pct_q20: warning, not a crash, Q still not flagged",
          f"old config: rc={rc} BP_q={t3['BP_q']['flag']}")

    # ---- 2d. qc_report.html is client-safe ---------------------------------
    # A run under a path holding every facility marker, with an old config
    # and a second counts file beside it so warnings carry absolute paths.
    leak = build_run(tmp, "g/data/vz35/scratch/home/ONT_16S_LEAK", "post")
    shutil.copy(leak / "run_log_20261013" / "read_counts_summary.txt", leak)
    leak_cfg = leak / "old_config.yaml"
    leak_cfg.write_text(old_cfg.read_text())
    run_report(leak, sheet, config=leak_cfg)
    lpage = (leak / "qc_report" / "qc_report.html").read_text(encoding="utf-8")
    flat = lpage.replace("\\\\", "/").replace("\\", "/")
    found = [s for s in INTERNAL if s in flat]
    lblob = html_blob(leak / "qc_report" / "qc_report.html")
    check(not found and "output_dir" not in lblob["run"]
          and lblob["run"]["name"] == "ONT_16S_LEAK"
          and any("old_config.yaml: metrics.pct_q20" in w
                  for w in lblob["run"]["warnings"]),
          "qc_report.html has no facility paths (/g/data, /scratch, /home, "
          "vz35); run name kept, warnings cut to file names",
          f"qc_report.html leaks {found}; run={lblob.get('run')}")
    ls = {s["id"]: s for s in lblob["samples"]}
    check(ls["Tiny"]["comment"] == "low DNA input" and ls["S01"]["comment"] == ""
          and ls["NTC1"]["type"] == "NTC" and ls["NTC1"]["flag"] == "CTRL",
          "qc_report.html data keeps Comment and Sample_Type (controls -> CTRL)",
          f"comment/type: {ls['Tiny']} {ls['NTC1']}")
    check('{ key: "comment", label: "Comment"' in lpage
          and '{ key: "type", label: "Type"' not in lpage,
          "sample table shows Comment instead of Type",
          "sample table columns not as expected")

    leak_pages = sorted((leak / "integrated_demultiplexing").glob("*/*.html"))
    check(len(leak_pages) == 2, "per-client HTML present in the leak run",
          f"leak run per-client HTML: {leak_pages}")
    for f in leak_pages:
        fp = f.read_text(encoding="utf-8").replace("\\\\", "/").replace("\\", "/")
        found = [s for s in INTERNAL if s in fp]
        check(not found, f"{f.name} has no facility paths",
              f"{f.name} leaks {found}")

    # No integrated_demultiplexing/: per-client HTML falls back to per_client/.
    integ_l = leak / "integrated_demultiplexing"
    integ_l.rename(leak / "moved_away")
    try:
        run_report(leak, sheet)
        fb = sorted(p.name for p in (leak / "qc_report" / "per_client").glob("*.html"))
        rq = (leak / "qc_report" / "run_qc.txt").read_text()
    finally:
        (leak / "moved_away").rename(integ_l)
    check(fb == ["ClientA_qc_report.html", "ClientB_qc_report.html"]
          and "written to qc_report/per_client/ instead" in rq,
          "no integrated_demultiplexing/: per-client HTML in per_client/, with a warning",
          f"fallback: {fb}")

    # ---- 2e. per-client reports --------------------------------------------
    pc = qd / "per_client"
    integ = out / "integrated_demultiplexing"
    run_flags = {s["id"]: s["flag"] for s in html_blob(qd / "qc_report.html")["samples"]}
    placed = sorted(p.relative_to(integ).as_posix() for p in integ.glob("*/*.html"))
    check(placed == ["ClientA/ClientA_qc_report.html", "ClientB/ClientB_qc_report.html"]
          and not list(pc.glob("*.html")) and not list(integ.glob("*.html")),
          "<Client>_qc_report.html written into integrated_demultiplexing/<Client>/, "
          "none for shared controls, none left in per_client/",
          f"per-client HTML: {placed}, per_client/: {sorted(p.name for p in pc.glob('*'))}")
    for me, other in (("ClientA", "ClientB"), ("ClientB", "ClientA")):
        path = integ / me / f"{me}_qc_report.html"
        page = path.read_text(encoding="utf-8")
        blob = html_blob(path)
        ids = {s["id"] for s in blob["samples"]}
        mine = {s for s in SAMPLES if CLIENT[s] == me} | {"NTC1"}
        theirs = [s for s in SAMPLES if CLIENT[s] == other]
        leaked = [x for x in (other, *(f'"{s}"' for s in theirs),
                              *(BARCODE[s] for s in theirs),
                              *(COMMENT[s] for s in theirs if s in COMMENT))
                  if x in page]
        check(ids == mine and not leaked,
              f"{me} report: all {me} samples + shared control, nothing of {other}",
              f"{me} report: ids {sorted(ids)} (want {sorted(mine)}), leaked {leaked}")
        diff = {s["id"]: (s["flag"], run_flags[s["id"]]) for s in blob["samples"]
                if s["flag"] != run_flags[s["id"]]}
        counts = blob["run"]["flag_counts"]
        check(not diff and sum(counts.values()) == len(mine)
              and counts["CTRL"] == len({"NTC1", "POS_B"} & mine)
              and blob["run"]["client"] == me and not blob["run"]["warnings"]
              and not blob["run"]["catch_all"],
              f"{me} report: flags identical to the run report, counts "
              f"recomputed ({counts}), no warnings / catch-all",
              f"{me} report: flag diffs {diff}, counts {counts}, run {blob['run']}")
    check((pc / "ClientA_qc.tsv").is_file() and "NTC1" in table(pc / "ClientA_qc.tsv")
          and "POS_B" not in table(pc / "ClientA_qc.tsv"),
          "per-client TSV: shared control included, other client's control not",
          f"ClientA_qc.tsv rows: {sorted(table(pc / 'ClientA_qc.tsv'))}")
    tp = (qd / "per_sample_qc.tsv")
    check(table(tp)["NTC1"]["client"] == "" and table(tp)["NTC1"]["flag"] == "CTRL",
          "control with empty Client kept in the run table as CTRL",
          f"NTC1 row: {table(tp).get('NTC1')}")

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
    check(t2["S01"]["median_q"] == "NA" and t2["BP_q"]["flag"] == "PASS"
          and "median_q is NA" in rq,
          "no .q.tsv: median_q NA with a warning, other metrics unaffected",
          f"no-Q run: {t2['S01']} / {t2['BP_q']['flag']}")
    msg = expect_exit(run_report, out, sheet, "--demux-qc-dir", str(tmp / "nope"))
    check("demux QC dir not found" in msg, "missing demux_qc/ fails loudly",
          f"missing demux_qc message: {msg!r}")

    # ---- 6. run-level relative mode ----------------------------------------
    sheet_rows = [qc.SheetSample("C", f"R{i:02d}", order=i) for i in range(12)]

    def synthetic(primer_pct, reads=30000):
        rows = {s.sample_id: qc.CountsRow(s.sample_id, reads, 1000, 1000, 100.0,
                                          970, primer_pct(i))
                for i, s in enumerate(sheet_rows)}
        profs = {s.sample_id: qc.DemuxProfile(
                     s.sample_id, reads=1000, filter_stage="post",
                     modal_length=1500, len_bins={1500: 1000}, has_q=True,
                     q_win={17: 1000})
                 for s in sheet_rows}
        return qc.Counts(Path("x"), "post", rows), profs

    # Whole run at ~70% primer found (FAIL), one sample far below the rest.
    counts, profs = synthetic(lambda i: 30.0 if i == 0 else 68.0 + i % 5)
    ev = qc.evaluate(sheet_rows, counts, profs, cfg, [])
    p_run = next(m for m in ev.runs if m.metric == "primer_found")
    flags = {r.sample.sample_id: r.flag for r in ev.results}
    check(p_run.run_level == 2 and p_run.relative
          and any(ln.startswith("RUN: LOW_PRIMER (median ") and ", FAIL)" in ln
                  for ln in qc.run_flag_lines(ev.runs))
          and not any("Q20" in ln for ln in qc.run_flag_lines(ev.runs)),
          "run flag RUN: LOW_PRIMER raised and metric switched to relative mode",
          f"run-level primer_found: {p_run}")
    check(flags["R00"] == "FAIL"
          and all(f == "PASS" for s, f in flags.items() if s != "R00"),
          "relative mode flags only the outlier, not the whole run",
          f"relative flags: {flags}")
    check("z=" in ev.results[0].reasons, f"outlier reason: {ev.results[0].reasons}",
          f"outlier reason lacks z: {ev.results[0].reasons}")

    # delivered_reads stays absolute even when the whole run is low
    counts, profs = synthetic(lambda i: 97.0)
    for c in counts.rows.values():
        c.reads = 9000
    ev = qc.evaluate(sheet_rows, counts, profs, cfg, [])
    rd = next(m for m in ev.runs if m.metric == "delivered_reads")
    check(rd.run_level == 1 and not rd.relative
          and all(r.levels["delivered_reads"] == 1 for r in ev.results),
          "delivered_reads: run flag raised but samples stay absolute",
          f"delivered_reads run: {rd}")

    # controls are left out of run medians
    counts, profs = synthetic(lambda i: 97.0)
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
