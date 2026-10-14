#!/usr/bin/env python3
"""
qc_report.py

Standalone per-sample QC report (PASS / CHECK / FAIL) for a finished 16S run.

Judges only what the facility is responsible for: did each sample amplify and
sequence into usable reads. Chimera %, classification and NTC read ratios are
the client's responsibility and are deliberately not reported.

A pure reader: it never opens a FASTQ and runs no Minibar or Chopper, so it
takes seconds. Inputs:

    <output_dir>/demux_qc/<SampleID>.tsv     length histogram, modal length
    <output_dir>/demux_qc/<SampleID>.q.tsv   mean-read-Q histogram
    <output_dir>/run_log_*/read_counts_summary.txt
                                             delivered reads, raw_demuxed,
                                             pct_retained, pct_pass1_kept
    sample sheet                             Client, Sample_ID, Barcode,
                                             optional Sample_Type / Comment

Flagged metrics, all exact:

    delivered_reads  Reads in read_counts_summary.txt
    primer_found     pct_pass1_kept
    long_frac        1 - pct_retained. MinKNOW already drops < 1 kb and < Q15,
                     so in post mode Chopper only removes > 2 kb reads.
                     NA in pre mode, where Chopper ran before demux.
    core_frac        reads in core_window / reads in 1,000-2,000 bp
    modal_dev_bp     |modal length - run median modal length|

Informational, never flagged:

    median_q         median mean-read-Q bin of the reads in 1,000-2,000 bp.
                     Every delivered read already passed MinKNOW Q15, so
                     read quality is shown but not judged.

Anything < 1 kb never reaches the pipeline (MinKNOW), so a short-product
failure can only show as low delivered reads or a low core_frac.

Run at the end of run_script.qsub, or by hand / run_qc_report.qsub on any
past run. Thresholds come from qc_config.yaml; change them and re-run.

PASS / CHECK / FAIL is triage, not sign-off: staff review qc_report.html (or
run_qc.txt and the heatmaps) before release.

Usage:
    python3 qc_report.py \\
        --output-dir  /path/to/run_output \\
        --samplesheet /path/to/16s_samplesheet.csv \\
        [--config qc_config.yaml] [--counts read_counts_summary.txt] \\
        [--png-drilldown]
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import json
import re
import shutil
import statistics
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

# Local imports: sibling modules in the same tools/ dir.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import filter_chopper_demux_minibar as core  # noqa: E402
from demux_qc import BIN_WIDTH, CATCH_ALL, median_q_bin  # noqa: E402

log = core.log
banner = core.banner

DEFAULT_CONFIG = Path(__file__).resolve().parent / "qc_config.yaml"

# Flagged metrics. Order is also the tie-break order for reason codes.
METRICS = ("delivered_reads", "primer_found", "long_frac", "core_frac",
           "modal_dev_bp")
# Reported but never flagged.
INFO_METRICS = ("median_q",)
# NA below min_reads_for_stats; delivered_reads is always judged.
STATS_METRICS = (*METRICS[1:], *INFO_METRICS)
REASON_CODE = {
    "delivered_reads": "READS", "primer_found": "PRIMER", "long_frac": "LONG",
    "core_frac": "LOW_CORE", "modal_dev_bp": "MODAL",
}
RUN_CODE = {
    "delivered_reads": "LOW_READS", "primer_found": "LOW_PRIMER",
    "long_frac": "LONG", "core_frac": "LOW_CORE",
    "modal_dev_bp": "MODAL_SPREAD",
}
# Metrics dropped from the flagged set, with the warning an old config
# that still lists them gets.
RETIRED_METRICS = {
    "pct_q20": "ignored: Q is no longer flagged (config v2.2)",
}
LEVEL_NAME = {0: "PASS", 1: "CHECK", 2: "FAIL"}
FLAG_ORDER = {"FAIL": 0, "CHECK": 1, "PASS": 2, "CTRL": 3}

# MAD -> sigma for normally distributed data, so k_check/k_fail read as
# "standard deviations from the run median".
MAD_TO_SIGMA = 1.4826

# Plot ranges. Reads outside are clamped into the end bins.
LEN_PLOT_MIN = 1000
LEN_PLOT_MAX = 4000
Q_PLOT_MIN = 10
Q_PLOT_MAX = 40
HEATMAP_MAX_ROWS = 384


# =============================================================================
# Config
# =============================================================================
def load_config(path: Path) -> dict[str, Any]:
    try:
        import yaml
    except ImportError:
        sys.exit("ERROR: PyYAML is required to read the QC config. On Gadi use "
                 "the QC venv's python:\n"
                 "  /g/data/vz35/ONT_16s_workflow/tools/qc_env/bin/python3")
    if not path.is_file():
        sys.exit(f"ERROR: QC config not found: {path}")
    raw = path.read_bytes()
    cfg = yaml.safe_load(raw)

    def need(section: str, *keys: str) -> None:
        if not isinstance(cfg.get(section), dict):
            sys.exit(f"ERROR: {path}: missing section '{section}'")
        missing = [k for k in keys if k not in cfg[section]]
        if missing:
            sys.exit(f"ERROR: {path}: section '{section}' is missing "
                     f"{', '.join(missing)}")

    for key in ("config_version", "min_reads_for_stats"):
        if key not in cfg:
            sys.exit(f"ERROR: {path}: missing {key}")
    need("filter", "minknow_min_len", "minknow_min_q", "chopper_max_len")
    need("core_window", "min", "max")
    need("metrics", *METRICS)
    need("run_level", "outlier_k_check", "outlier_k_fail", "mad_floor")
    need("controls", "types")

    # Retired / unknown metric entries are dropped with a warning, so an old
    # config still runs. The caller adds _warnings to the run's warnings.
    cfg_warnings: list[str] = []
    for name in [m for m in cfg["metrics"] if m not in METRICS]:
        del cfg["metrics"][name]
        why = RETIRED_METRICS.get(name, "ignored: not a known metric")
        cfg_warnings.append(f"{path}: metrics.{name} {why}")

    for name in METRICS:
        spec = cfg["metrics"][name]
        for k in ("direction", "check", "fail"):
            if k not in spec:
                sys.exit(f"ERROR: {path}: metrics.{name} is missing '{k}'")
        if spec["direction"] == "higher_better":
            ok = spec["fail"] <= spec["check"]
        elif spec["direction"] == "lower_better":
            ok = spec["fail"] >= spec["check"]
        else:
            sys.exit(f"ERROR: {path}: metrics.{name}.direction must be "
                     f"higher_better or lower_better")
        if not ok:
            sys.exit(f"ERROR: {path}: metrics.{name}: fail threshold is "
                     f"less strict than check")

    f, cw = cfg["filter"], cfg["core_window"]
    if cw["min"] % BIN_WIDTH or cw["max"] % BIN_WIDTH:
        sys.exit(f"ERROR: {path}: core_window ends must be multiples of the "
                 f"{BIN_WIDTH} bp demux_qc bin width, got {cw['min']}-{cw['max']}")
    if not f["minknow_min_len"] <= cw["min"] < cw["max"] <= f["chopper_max_len"]:
        sys.exit(f"ERROR: {path}: core_window {cw['min']}-{cw['max']} must lie "
                 f"inside {f['minknow_min_len']}-{f['chopper_max_len']}")
    if f["minknow_min_len"] % BIN_WIDTH or f["chopper_max_len"] % BIN_WIDTH:
        sys.exit(f"ERROR: {path}: filter lengths must be multiples of "
                 f"{BIN_WIDTH} bp")

    cfg["run_level"].setdefault("min_samples", 8)
    cfg["controls"]["types"] = [str(t).upper() for t in cfg["controls"]["types"]]
    cfg["_path"] = str(path)
    cfg["_sha256"] = hashlib.sha256(raw).hexdigest()[:12]
    cfg["_warnings"] = cfg_warnings
    return cfg


def config_label(cfg: dict[str, Any]) -> str:
    return f"{cfg['config_version']} (sha256 {cfg['_sha256']})"


# =============================================================================
# Sample sheet
# =============================================================================
@dataclass
class SheetSample:
    client: str
    sample_id: str
    barcode: str = ""
    location: str = ""
    comment: str = ""
    sample_type: str = "SAMPLE"
    order: int = 0                 # row order in the sheet (= plate order)


def load_samplesheet(path: Path, control_types: list[str]) -> list[SheetSample]:
    """Client, Sample_ID, Barcode, plus optional Location, Comment and
    Sample_Type (SAMPLE / NTC / BLANK / POS_CTRL; blank means SAMPLE).

    Client and Sample_ID are sanitised exactly as the pipeline does, so they
    match the per-client folders and the demux_qc/<SampleID>.tsv names.
    """
    if not path.is_file():
        sys.exit(f"ERROR: sample sheet not found: {path}")
    allowed = {"SAMPLE", *control_types}
    out: list[SheetSample] = []
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh)
        fields = reader.fieldnames or []
        missing = [c for c in ("Client", "Sample_ID", "Barcode") if c not in fields]
        if missing:
            sys.exit(f"ERROR: {path} is missing column(s): {', '.join(missing)}")
        for i, row in enumerate(reader, start=2):
            client = (row.get("Client") or "").strip()
            sid = (row.get("Sample_ID") or "").strip()
            stype = (row.get("Sample_Type") or "").strip().upper() or "SAMPLE"
            # A control may leave Client empty: it is shared by every client.
            if not sid or (not client and stype not in control_types):
                continue
            if stype not in allowed:
                sys.exit(f"ERROR: {path} line {i}: Sample_Type {stype!r} is not "
                         f"one of {', '.join(sorted(allowed))}")
            out.append(SheetSample(
                client=core.sanitise(client) if client else "",
                sample_id=core.sanitise_sample_id(sid),
                barcode=(row.get("Barcode") or "").strip(),
                location=(row.get("Location") or "").strip(),
                comment=(row.get("Comment") or "").strip(),
                sample_type=stype,
                order=len(out),
            ))
    if not out:
        sys.exit(f"ERROR: no samples found in {path}")
    return out


# Client values that mark a control as shared by every client.
SHARED_CLIENTS = {"", "ALL"}


def is_shared(s: SheetSample, control_types: list[str]) -> bool:
    """A control with an empty Client or Client = ALL goes into every
    per-client report; a control with a named client only into that one."""
    return s.sample_type in control_types and s.client.upper() in SHARED_CLIENTS


def report_clients(sheet: list[SheetSample], control_types: list[str]) -> list[str]:
    """Clients that get a per-client report: every client of a non-shared row."""
    return sorted({s.client for s in sheet if not is_shared(s, control_types)})


def client_file_stem(client: str) -> str:
    """Client name for file names; same rule as Sample_ID sanitising."""
    return core.sanitise_sample_id(client)


# =============================================================================
# read_counts_summary.txt
# =============================================================================
@dataclass
class CountsRow:
    sample: str
    reads: int
    raw_demuxed: int | None
    post_chopper: int | None
    pct_retained: float | None     # percent, as printed
    pass1: int | None
    pct_pass1_kept: float | None   # percent, as printed


@dataclass
class Counts:
    path: Path
    filter_stage: str | None       # from the ' Filter stage:' header line
    rows: dict[str, CountsRow]


# Columns the report depends on, as they appear in the header line.
_COUNT_COLUMNS = ("Sample", "Reads", "% of Input", "raw_demuxed",
                  "post_chopper", "pct_retained", "post_cutadapt_pass1",
                  "pct_pass1_kept")


def find_counts_file(output_dir: Path) -> tuple[Path, list[str]]:
    """Locate read_counts_summary.txt: run_log_<date>/ after Step 7, or the
    run dir itself if Step 7 never ran. Newest wins if there are several."""
    cands = list(output_dir.glob("run_log_*/read_counts_summary.txt"))
    top = output_dir / "read_counts_summary.txt"
    if top.is_file():
        cands.append(top)
    if not cands:
        sys.exit(f"ERROR: no read_counts_summary.txt under {output_dir} "
                 f"(looked in run_log_*/ and the run dir). Pass --counts.")
    cands.sort(key=lambda p: p.stat().st_mtime)
    warnings = []
    if len(cands) > 1:
        warnings.append(
            f"{len(cands)} read_counts_summary.txt files found; using the newest "
            f"({cands[-1]}). Pass --counts to choose another.")
    return cands[-1], warnings


def _num(cell: str, path: Path, lineno: int, column: str) -> float | None:
    cell = cell.strip()
    if not cell:
        return None
    try:
        return float(cell.replace(",", "").removesuffix("%"))
    except ValueError:
        sys.exit(f"ERROR: {path} line {lineno}: column '{column}' holds "
                 f"{cell!r}, not a number -- has the table layout changed?")


def _int(v: float | None) -> int | None:
    return None if v is None else int(v)


def parse_read_counts(path: Path) -> Counts:
    """
    Per-sample rows of read_counts_summary.txt, keyed by sample ID (no
    sample_ prefix), plus the filter stage from the header block. The TOTAL
    rows below the second separator are ignored (known bugs).

    The table is fixed-width with right-aligned numeric cells, and several
    cells are legitimately blank (raw_demuxed and pct_retained in pre mode),
    so whitespace splitting cannot tell which column a number belongs to.
    Each cell is instead cut at the end position of its header label, shifted
    by however far an over-long sample name pushed the row to the right.
    """
    if not path.is_file():
        sys.exit(f"ERROR: read counts summary not found: {path}")
    lines = path.read_text().splitlines()

    stage = None
    for ln in lines:
        if ln.strip().startswith("Filter stage:"):
            stage = ln.split(":", 1)[1].strip() or None
            break

    hdr_idx = next((i for i, ln in enumerate(lines)
                    if ln.startswith("Sample") and "Reads" in ln), None)
    if hdr_idx is None:
        sys.exit(f"ERROR: {path}: no header line starting with 'Sample' -- "
                 f"has the table layout changed?")
    header = lines[hdr_idx]
    missing = [c for c in _COUNT_COLUMNS if c not in header]
    if missing:
        sys.exit(f"ERROR: {path}: header is missing column(s): "
                 f"{', '.join(missing)}. qc_report.py needs the extended "
                 f"read_counts_summary.txt written by the current pipeline.")
    ends = {c: header.index(c) + len(c) for c in _COUNT_COLUMNS[1:]}
    tail = list(_COUNT_COLUMNS[3:])        # cut by position; may be blank

    rows: dict[str, CountsRow] = {}
    i = hdr_idx + 1
    if i < len(lines) and lines[i].startswith("-"):
        i += 1
    for lineno in range(i + 1, len(lines) + 1):
        line = lines[lineno - 1]
        if line.startswith("-") or line.startswith("="):
            break                          # TOTAL block follows; ignore it
        if not line.strip():
            continue
        tokens = line.split()
        if len(tokens) < 3 or not tokens[2].endswith("%"):
            sys.exit(f"ERROR: {path} line {lineno}: expected 'Sample Reads "
                     f"%_of_Input ...', got {line!r}")
        sample, reads_cell, pct_cell = tokens[:3]
        # Where '% of Input' really ends on this line vs in the header.
        pos = line.index(reads_cell, len(sample))
        pos = line.index(pct_cell, pos + len(reads_cell)) + len(pct_cell)
        shift = pos - ends["% of Input"]

        cells: dict[str, float | None] = {}
        prev = ends["% of Input"]
        for col in tail:
            cells[col] = _num(line[prev + shift:ends[col] + shift],
                              path, lineno, col)
            prev = ends[col]

        reads = _num(reads_cell, path, lineno, "Reads")
        sid = sample.removeprefix("sample_")
        rows[sid] = CountsRow(
            sample=sid,
            reads=int(reads or 0),
            raw_demuxed=_int(cells["raw_demuxed"]),
            post_chopper=_int(cells["post_chopper"]),
            pct_retained=cells["pct_retained"],
            pass1=_int(cells["post_cutadapt_pass1"]),
            pct_pass1_kept=cells["pct_pass1_kept"],
        )
    if not rows:
        sys.exit(f"ERROR: {path}: no per-sample rows found under the header")
    return Counts(path, stage, rows)


# =============================================================================
# demux_qc/
# =============================================================================
@dataclass
class DemuxProfile:
    sample_id: str
    reads: int = 0
    filter_stage: str = ""
    length_window: str = ""
    modal_length: int | None = None
    len_bins: dict[int, int] = field(default_factory=dict)  # 50 bp bin start
    q_all: dict[int, int] = field(default_factory=dict)     # floor(mean Q)
    q_win: dict[int, int] = field(default_factory=dict)     # in-window reads
    has_q: bool = False


def _read_hist_tsv(path: Path) -> tuple[dict[str, str], list[dict[str, str]]]:
    """'#'-prefixed key<TAB>value header block, then a TSV table."""
    meta: dict[str, str] = {}
    body: list[str] = []
    for ln in path.read_text().splitlines():
        if ln.startswith("#"):
            k, _, v = ln[1:].strip().partition("\t")
            meta[k] = v
        elif ln.strip():
            body.append(ln)
    return meta, list(csv.DictReader(body, delimiter="\t"))


def read_demux_profile(qc_dir: Path, sid: str) -> DemuxProfile | None:
    path = qc_dir / f"{sid}.tsv"
    if not path.is_file():
        return None
    meta, rows = _read_hist_tsv(path)
    for col in ("bin_start", "count"):
        if rows and col not in rows[0]:
            sys.exit(f"ERROR: {path}: missing column '{col}'")
    bw = meta.get("bin_width_bp")
    if bw is not None and int(bw) != BIN_WIDTH:
        sys.exit(f"ERROR: {path}: bin width {bw} bp, expected {BIN_WIDTH}")
    reads = int(meta.get("reads", 0))
    prof = DemuxProfile(
        sample_id=sid,
        reads=reads,
        filter_stage=meta.get("filter_stage", ""),
        length_window=meta.get("length_window", ""),
        modal_length=int(meta["modal_length"]) if reads and "modal_length" in meta
        else None,
        len_bins={int(r["bin_start"]): int(r["count"]) for r in rows
                  if int(r["count"])},
    )
    qpath = qc_dir / f"{sid}.q.tsv"
    if qpath.is_file():
        _, qrows = _read_hist_tsv(qpath)
        for col in ("q_bin", "count_all", "count_in_window"):
            if qrows and col not in qrows[0]:
                sys.exit(f"ERROR: {qpath}: missing column '{col}'")
        prof.q_all = {int(r["q_bin"]): int(r["count_all"]) for r in qrows
                      if int(r["count_all"])}
        prof.q_win = {int(r["q_bin"]): int(r["count_in_window"]) for r in qrows
                      if int(r["count_in_window"])}
        prof.has_q = True
    return prof


def load_demux_qc(qc_dir: Path, sample_ids: list[str],
                  warnings: list[str]) -> dict[str, DemuxProfile]:
    if not qc_dir.is_dir():
        sys.exit(f"ERROR: demux QC dir not found: {qc_dir} -- was the run "
                 f"made with --no-demux-qc? Pass --demux-qc-dir if it moved.")
    profiles: dict[str, DemuxProfile] = {}
    no_q = []
    for sid in [*sample_ids, *CATCH_ALL]:
        prof = read_demux_profile(qc_dir, sid)
        if prof is None:
            continue
        profiles[sid] = prof
        if not prof.has_q and sid not in CATCH_ALL:
            no_q.append(sid)
    if no_q:
        warnings.append(
            f"{len(no_q)} sample(s) have no <SampleID>.q.tsv (demux_qc.py "
            f"predates the Q histogram); median_q is NA for them")
    return profiles


# =============================================================================
# Metrics from histograms
# =============================================================================
def core_frac(prof: DemuxProfile, cfg: dict[str, Any]) -> float | None:
    """Reads in core_window / reads in the delivered length window, both
    half-open on 50 bp bin starts (a read of exactly chopper_max_len bp falls
    in the next bin and is left out of both)."""
    f, cw = cfg["filter"], cfg["core_window"]
    den = sum(n for b, n in prof.len_bins.items()
              if f["minknow_min_len"] <= b < f["chopper_max_len"])
    num = sum(n for b, n in prof.len_bins.items() if cw["min"] <= b < cw["max"])
    return num / den if den else None


def median_q(prof: DemuxProfile) -> int | None:
    """Median floor(mean read Q) of the in-window reads (lower median).
    Informational only: read quality is not flagged."""
    if not prof.has_q:
        return None
    return median_q_bin(prof.q_win)


# =============================================================================
# Flagging
# =============================================================================
@dataclass
class SampleResult:
    sample: SheetSample
    n: int                                               # raw demuxed reads
    raw_demuxed: int | None
    modal_length: int | None
    values: dict[str, float | None]
    na: dict[str, str] = field(default_factory=dict)     # metric -> NA text
    levels: dict[str, int] = field(default_factory=dict)
    z: dict[str, float] = field(default_factory=dict)    # relative mode only
    flag: str = "PASS"
    reasons: str = ""


@dataclass
class MetricRun:
    metric: str
    median: float | None
    mad: float | None
    n: int
    run_level: int = 0
    relative: bool = False


@dataclass
class Evaluation:
    results: list[SampleResult]
    runs: list[MetricRun]
    filter_stage: str
    modal_median: float | None


def abs_level(value: float, spec: dict[str, Any]) -> int:
    if spec["direction"] == "higher_better":
        return 2 if value < spec["fail"] else 1 if value < spec["check"] else 0
    return 2 if value > spec["fail"] else 1 if value > spec["check"] else 0


def fmt_value(metric: str, value: float) -> str:
    if metric == "delivered_reads":
        return f"{int(value):,}"
    if metric == "modal_dev_bp":
        return f"{value:.0f}bp"
    return f"{100 * value:.0f}%"


def resolve_filter_stage(counts: Counts, profiles: dict[str, DemuxProfile],
                         warnings: list[str]) -> str:
    """Counts header first, demux_qc headers as a fallback; warn on any
    disagreement. 'pre' or 'post'."""
    demux = {p.filter_stage for p in profiles.values() if p.filter_stage}
    stage = counts.filter_stage or (demux.pop() if len(demux) == 1 else None)
    if stage is None:
        # Pre mode is the only one that leaves pct_retained blank everywhere.
        any_ret = any(r.pct_retained is not None for r in counts.rows.values())
        stage = "post" if any_ret else "pre"
        warnings.append(f"filter stage not recorded; inferred '{stage}' from "
                        f"the pct_retained column")
    elif demux and demux != {stage}:
        warnings.append(f"read_counts_summary.txt says filter stage '{stage}' "
                        f"but demux_qc/ says {', '.join(sorted(demux))}")
    if stage not in ("pre", "post"):
        sys.exit(f"ERROR: unknown filter stage {stage!r} in {counts.path}")
    return stage


def evaluate(sheet: list[SheetSample], counts: Counts,
             profiles: dict[str, DemuxProfile], cfg: dict[str, Any],
             warnings: list[str]) -> Evaluation:
    min_n = cfg["min_reads_for_stats"]
    rl = cfg["run_level"]
    controls = set(cfg["controls"]["types"])
    stage = resolve_filter_stage(counts, profiles, warnings)

    results: list[SampleResult] = []
    for s in sheet:
        row = counts.rows.get(s.sample_id)
        prof = profiles.get(s.sample_id)
        values: dict[str, float | None] = dict.fromkeys((*METRICS,
                                                         *INFO_METRICS))
        na: dict[str, str] = {}

        if row is None:
            warnings.append(f"{s.sample_id}: not in read_counts_summary.txt; "
                            f"delivered_reads taken as 0")
        if prof is None:
            warnings.append(f"{s.sample_id}: no demux_qc/{s.sample_id}.tsv "
                            f"(no reads demultiplexed?)")
        values["delivered_reads"] = row.reads if row else 0
        # Post mode: raw_demuxed. Pre mode has none; demux_qc then saw the
        # Chopper-filtered reads, which is the best count available.
        raw = row.raw_demuxed if row else None
        n = raw if raw is not None else (prof.reads if prof else 0)

        if n < min_n:
            for m in STATS_METRICS:
                na[m] = f"NA (n<{min_n})"
        else:
            if row is not None and row.pct_pass1_kept is not None:
                values["primer_found"] = row.pct_pass1_kept / 100.0
            if stage == "pre":
                na["long_frac"] = "NA (pre mode)"
            elif row is not None and row.pct_retained is not None:
                values["long_frac"] = 1.0 - row.pct_retained / 100.0
            if prof is not None:
                values["core_frac"] = core_frac(prof, cfg)
                values["median_q"] = median_q(prof)
        modal = prof.modal_length if prof and n >= min_n else None
        results.append(SampleResult(s, n, raw, modal, values, na))

    real = [r for r in results if r.sample.sample_type not in controls]

    # modal_dev_bp: distance from the run's median modal length.
    modals = [r.modal_length for r in real if r.modal_length is not None]
    modal_median = statistics.median(modals) if modals else None
    if modal_median is not None:
        for r in results:
            if r.modal_length is not None:
                r.values["modal_dev_bp"] = abs(r.modal_length - modal_median)

    runs: list[MetricRun] = []
    for m in METRICS:
        spec = cfg["metrics"][m]
        vals = [r.values[m] for r in real if r.values[m] is not None]
        if not vals:
            runs.append(MetricRun(m, None, None, 0))
            continue
        med = statistics.median(vals)
        mad = statistics.median(abs(v - med) for v in vals)
        mr = MetricRun(m, med, mad, len(vals), run_level=abs_level(med, spec))
        if mr.run_level and not spec.get("absolute_only", False):
            if len(vals) >= rl["min_samples"]:
                mr.relative = True
            else:
                warnings.append(
                    f"{m}: run median is at {LEVEL_NAME[mr.run_level]} but only "
                    f"{len(vals)} real sample(s) have a value (< min_samples "
                    f"{rl['min_samples']}); kept absolute thresholds")
        runs.append(mr)

    by_metric = {mr.metric: mr for mr in runs}
    for r in results:
        for m in METRICS:
            v = r.values[m]
            if v is None:
                continue
            spec = cfg["metrics"][m]
            mr = by_metric[m]
            if mr.relative:
                scale = max(MAD_TO_SIGMA * mr.mad,
                            spec.get("mad_floor", rl["mad_floor"]))
                z = ((mr.median - v) if spec["direction"] == "higher_better"
                     else (v - mr.median)) / scale
                r.z[m] = z
                r.levels[m] = (2 if z >= rl["outlier_k_fail"]
                               else 1 if z >= rl["outlier_k_check"] else 0)
            else:
                r.levels[m] = abs_level(v, spec)

        if r.sample.sample_type in controls:
            r.flag = "CTRL"
            continue
        worst = max(r.levels.values(), default=0)
        r.flag = LEVEL_NAME[worst]
        bad = sorted((m for m in METRICS if r.levels.get(m, 0) > 0),
                     key=lambda m: (-r.levels[m], METRICS.index(m)))
        parts = []
        for m in bad:
            txt = fmt_value(m, r.values[m])
            if m in r.z:
                txt += f",z={r.z[m]:.1f}"
            parts.append(f"{REASON_CODE[m]}({txt})")
        if parts:
            r.reasons = f"{r.flag}: " + ";".join(parts)
        elif r.n < min_n:
            r.reasons = f"stats NA (n={r.n})"
    return Evaluation(results, runs, stage, modal_median)


def run_flag_lines(runs: list[MetricRun]) -> list[str]:
    out = []
    for mr in runs:
        if not mr.run_level:
            continue
        mode = ("samples judged relative to the run (robust z)" if mr.relative
                else "samples judged on absolute thresholds")
        out.append(f"RUN: {RUN_CODE[mr.metric]} (median "
                   f"{fmt_value(mr.metric, mr.median)}, "
                   f"{LEVEL_NAME[mr.run_level]}) -> {mode}")
    return out


# =============================================================================
# Outputs
# =============================================================================
_TSV_COLUMNS = ["client", "sample", "barcode", "location", "type",
                "raw_demuxed", "delivered_reads", "primer_found", "long_frac",
                "core_frac", "modal_length", "modal_dev_bp", "median_q",
                "flag", "reasons", "comment"]


def _cell(r: SampleResult, m: str) -> str:
    v = r.values[m]
    if v is None:
        return r.na.get(m, "NA")
    if m in ("delivered_reads", "modal_dev_bp", "median_q"):
        return str(int(round(v)))
    return f"{v:.4f}"


def write_sample_tsv(path: Path, results: list[SampleResult],
                     header: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        for h in header:
            fh.write(f"# {h}\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(_TSV_COLUMNS)
        for r in results:
            s = r.sample
            vals = {m: _cell(r, m) for m in (*METRICS, *INFO_METRICS)}
            w.writerow([s.client, s.sample_id, s.barcode, s.location,
                        s.sample_type,
                        "" if r.raw_demuxed is None else r.raw_demuxed,
                        vals["delivered_reads"], vals["primer_found"],
                        vals["long_frac"], vals["core_frac"],
                        "" if r.modal_length is None else r.modal_length,
                        vals["modal_dev_bp"], vals["median_q"],
                        r.flag, r.reasons, s.comment])


def _stage_note(stage: str) -> str:
    if stage == "post":
        return "post (demultiplexed before Chopper; > 2 kb reads visible)"
    return ("pre (Chopper ran before demultiplexing): long_frac is NA and "
            "raw_demuxed is blank")


def tsv_header(cfg: dict[str, Any], ev: Evaluation) -> list[str]:
    f, cw = cfg["filter"], cfg["core_window"]
    lines = [
        f"generated\t{datetime.now():%Y-%m-%d %H:%M:%S}",
        f"config_version\t{config_label(cfg)}",
        f"filter_stage\t{_stage_note(ev.filter_stage)}",
        f"delivered_window\t{f['minknow_min_len']}-{f['chopper_max_len']} bp, "
        f"mean Q >= {f['minknow_min_q']}",
        f"core_window\t{cw['min']}-{cw['max']} bp",
        f"run_median_modal_length\t"
        f"{'NA' if ev.modal_median is None else f'{ev.modal_median:.0f}'}",
        f"fractions\t0-1; NA (n<{cfg['min_reads_for_stats']}) = too few "
        f"raw demuxed reads",
        "median_q\tinformational, not flagged: median mean-read-Q bin of "
        "the reads in the delivered window",
    ]
    lines += [f"run_flag\t{ln}" for ln in run_flag_lines(ev.runs)]
    return lines


def write_bins(path: Path, profiles: dict[str, DemuxProfile]) -> None:
    """Long format for re-plotting / re-thresholding. kind: len (50 bp bin
    start, all reads), q (floor mean Q, in-window reads), q_all (all reads)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "kind", "bin", "count"])
        for sid in sorted(profiles):
            p = profiles[sid]
            for kind, bins in (("len", p.len_bins), ("q", p.q_win),
                               ("q_all", p.q_all)):
                for b in sorted(bins):
                    w.writerow([sid, kind, b, bins[b]])


def write_run_qc(path: Path, args: argparse.Namespace, cfg: dict[str, Any],
                 counts: Counts, ev: Evaluation,
                 profiles: dict[str, DemuxProfile], warnings: list[str]) -> None:
    f, cw = cfg["filter"], cfg["core_window"]
    results, runs = ev.results, ev.runs
    L = [
        "========================================",
        " Per-sample QC report -- run summary",
        f" Generated:    {datetime.now():%Y-%m-%d %H:%M:%S}",
        f" Output dir:   {args.output_dir}",
        f" Counts:       {counts.path}",
        f" Demux QC:     {args.demux_qc_dir}",
        f" Config:       {cfg['_path']}",
        f" Version:      {config_label(cfg)}",
        f" Filter stage: {_stage_note(ev.filter_stage)}",
        "========================================",
    ]
    if ev.filter_stage == "pre":
        L += ["",
              " NOTE: pre-mode run. Chopper filtered the reads before Minibar,",
              " so reads > 2 kb were never demultiplexed: long_frac is NA for",
              " every sample and the length heatmap shows no > 2 kb tail."]
    L += [
        "",
        " Windows",
        f"   delivered:  {f['minknow_min_len']}-{f['chopper_max_len']} bp "
        f"(MinKNOW >= {f['minknow_min_len']} bp, Q{f['minknow_min_q']}; "
        f"Chopper <= {f['chopper_max_len']} bp)",
        f"   core:       {cw['min']}-{cw['max']} bp",
        f"   run median modal length: "
        f"{'NA' if ev.modal_median is None else f'{ev.modal_median:.0f} bp'}",
    ]

    L += ["", " Run-level flags"]
    L += [f"   {ln}" for ln in run_flag_lines(runs)] or ["   none"]

    L += ["", " Run medians (real samples, controls excluded)",
          f"   {'metric':<17} {'median':>9} {'MAD':>9} {'n':>5}  mode"]
    for mr in runs:
        if mr.median is None:
            L.append(f"   {mr.metric:<17} {'NA':>9} {'NA':>9} {0:>5}  -")
            continue
        fmt = ((lambda v: f"{v:,.0f}")
               if mr.metric in ("delivered_reads", "modal_dev_bp")
               else (lambda v: f"{v:.4f}"))
        L.append(f"   {mr.metric:<17} {fmt(mr.median):>9} {fmt(mr.mad):>9} "
                 f"{mr.n:>5}  {'relative' if mr.relative else 'absolute'}")

    L += ["", " Flag counts"]
    for flag in ("PASS", "CHECK", "FAIL", "CTRL"):
        L.append(f"   {flag:<6} {sum(r.flag == flag for r in results):>5}")

    # Share of all demultiplexed reads, from demux_qc (raw in post mode).
    total = sum(p.reads for p in profiles.values())
    L += ["", " Catch-all bins (share of all demultiplexed reads)"]
    for sid in CATCH_ALL:
        p = profiles.get(sid)
        n = p.reads if p else 0
        share = 100.0 * n / total if total else 0.0
        L.append(f"   {sid:<17} {n:>10,} reads  {share:5.1f}%")

    L += ["", " Warnings"]
    L += [f"   - {w}" for w in warnings] or ["   none"]
    L += [
        "",
        " PASS / CHECK / FAIL is triage, not sign-off. Review qc_report.html",
        " (or this file and length_heatmap.png / q_heatmap.png) before release.",
        "========================================",
    ]
    path.write_text("\n".join(L) + "\n")


# =============================================================================
# Interactive HTML report (qc_report.html)
# =============================================================================
# One self-contained file: the page template (vanilla JS + <canvas>, no
# external libraries or network references) with the run's data embedded as
# a JSON blob. Standard library only, so it is written even without
# matplotlib.
HTML_TEMPLATE = Path(__file__).resolve().parent / "qc_report_template.html"

METRIC_LABEL = {
    "delivered_reads": "Delivered reads", "primer_found": "Primer found",
    "long_frac": "Long reads (> 2 kb)", "core_frac": "Core window fraction",
    "modal_dev_bp": "Modal length deviation",
}


def _plot_counts(prof: DemuxProfile | None, kind: str) -> list[int]:
    """Counts per plot bin (see _edges), end bins clamped."""
    edges = _edges(kind)
    step = edges[1] - edges[0]
    out = [0] * len(edges)
    for b, n in _plot_bins(prof, kind).items():
        out[(b - edges[0]) // step] += n
    return out


def _run_median_profile(rows: list[list[int]]) -> list[float] | None:
    """Per-bin median of the row-normalised profiles, as a reference line."""
    fracs = [[c / sum(r) for c in r] for r in rows if sum(r)]
    if not fracs:
        return None
    return [round(statistics.median(col), 6) for col in zip(*fracs)]


# An absolute POSIX or Windows path, not a relative fragment like demux_qc/x.
_ABS_PATH = re.compile(r"(?<![\w.])(?:[A-Za-z]:)?(?:[\\/][^\\/\s:;,()'\"]+)+[\\/]?")


def _client_safe(text: str) -> str:
    """Reduce every absolute path to its last component. qc_report.html may
    be sent to clients, so it must not reveal facility paths (/g/data/...)."""
    return _ABS_PATH.sub(
        lambda m: re.split(r"[\\/]", m.group(0).rstrip("\\/"))[-1], text)


def _json_num(v: float | None) -> float | int | None:
    if v is None:
        return None
    return v if isinstance(v, int) else round(float(v), 6)


def html_data(args: argparse.Namespace, cfg: dict[str, Any], ev: Evaluation,
              profiles: dict[str, DemuxProfile], warnings: list[str],
              client: str | None = None) -> dict[str, Any]:
    """Data for qc_report.html. With client, only that client's samples and
    the shared controls are embedded, and anything naming other samples
    (warnings, catch-all counts) is left out, so the page can be sent to the
    client. Flags, run medians and the run-median overlay always come from
    the whole run, so a sample has the same flag in every report."""
    f, cw, rl = cfg["filter"], cfg["core_window"], cfg["run_level"]
    control_types = cfg["controls"]["types"]
    controls = set(control_types)
    by_metric = {mr.metric: mr for mr in ev.runs}
    shown = [r for r in ev.results
             if client is None or r.sample.client == client
             or is_shared(r.sample, control_types)]
    shown_ids = {id(r) for r in shown}

    metrics = []
    for m in METRICS:
        spec, mr = cfg["metrics"][m], by_metric[m]
        metrics.append({
            "key": m, "label": METRIC_LABEL[m], "code": REASON_CODE[m],
            "direction": spec["direction"], "check": spec["check"],
            "fail": spec["fail"], "relative": mr.relative,
            "run_level": mr.run_level, "median": _json_num(mr.median),
        })

    samples, real_len, real_q = [], [], []
    for r in ev.results:
        s, prof = r.sample, profiles.get(r.sample.sample_id)
        lc, qc_ = _plot_counts(prof, "len"), _plot_counts(prof, "q")
        if s.sample_type not in controls:
            real_len.append(lc)        # whole run: an anonymous aggregate
            real_q.append(qc_)
        if id(r) not in shown_ids:
            continue
        samples.append({
            "id": s.sample_id, "client": s.client, "barcode": s.barcode,
            "shared": is_shared(s, control_types),
            "location": s.location, "type": s.sample_type,
            "comment": s.comment, "order": s.order,
            "flag": r.flag, "reasons": r.reasons, "n": r.n,
            "raw": r.raw_demuxed, "modal_length": r.modal_length,
            "values": {m: _json_num(v) for m, v in r.values.items()},
            "na": r.na, "levels": r.levels,
            "z": {m: round(z, 3) for m, z in r.z.items()},
            "len": lc, "q": qc_,
        })

    # Internal only: warnings name samples, catch-all bins are run-wide.
    catch_all = []
    if client is None:
        total = sum(p.reads for p in profiles.values())
        for sid in CATCH_ALL:
            n = profiles[sid].reads if sid in profiles else 0
            catch_all.append({"sample": sid, "reads": n,
                              "share": 100.0 * n / total if total else 0.0})

    len_edges, q_edges = _edges("len"), _edges("q")
    return {
        "run": {
            # Run name only: no output_dir or other facility paths, every
            # HTML page may end up with a client.
            "name": args.output_dir.resolve().name,
            "client": client,
            "generated": f"{datetime.now():%Y-%m-%d %H:%M:%S}",
            "config_version": config_label(cfg),
            "filter_stage": ev.filter_stage,
            "filter_stage_note": _stage_note(ev.filter_stage),
            "flag_counts": {fl: sum(r.flag == fl for r in shown)
                            for fl in ("PASS", "CHECK", "FAIL", "CTRL")},
            # Anonymous text, e.g. "RUN: LONG (median 18%, FAIL) -> ...".
            "run_flags": run_flag_lines(ev.runs),
            "warnings": ([_client_safe(w) for w in warnings]
                         if client is None else []),
            "catch_all": catch_all,
            "core_window": {"min": cw["min"], "max": cw["max"]},
            "delivered": {"min": f["minknow_min_len"],
                          "max": f["chopper_max_len"]},
            "k_check": rl["outlier_k_check"], "k_fail": rl["outlier_k_fail"],
            "windows_note": (
                f"Delivered window {f['minknow_min_len']}-"
                f"{f['chopper_max_len']} bp; core window {cw['min']}-"
                f"{cw['max']} bp; config {config_label(cfg)}."),
        },
        "axes": {
            "len": {"min": len_edges[0], "max": len_edges[-1],
                    "step": len_edges[1] - len_edges[0], "n": len(len_edges)},
            "q": {"min": q_edges[0], "max": q_edges[-1],
                  "step": q_edges[1] - q_edges[0], "n": len(q_edges)},
        },
        "run_median": {"len": _run_median_profile(real_len),
                       "q": _run_median_profile(real_q)},
        "metrics": metrics,
        "samples": samples,
    }


def write_html(path: Path, data: dict[str, Any]) -> None:
    if not HTML_TEMPLATE.is_file():
        sys.exit(f"ERROR: HTML template not found: {HTML_TEMPLATE} -- it must "
                 f"sit beside qc_report.py")
    # '<' escaped so no string in the data (a sample Comment, say) can close
    # the <script> element early; the JSON unicode escape is still valid JSON.
    blob = json.dumps(data, ensure_ascii=False, allow_nan=False,
                      separators=(",", ":")).replace("<", "\\u003c")
    client = data["run"]["client"]
    title = html.escape(f"16S QC report - {data['run']['name']}"
                        + (f" - {client}" if client else ""))
    page = (HTML_TEMPLATE.read_text(encoding="utf-8")
            .replace("__TITLE__", title)
            .replace("__QC_DATA__", blob))
    path.write_text(page, encoding="utf-8")


# =============================================================================
# Plots
# =============================================================================
# Sequential blue ramp, surface -> dark, for row-normalised fractions.
_SEQ = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
        "#256abf", "#184f95", "#0d366b"]
_FLAG_COLOUR = {"PASS": "#0ca30c", "CHECK": "#fab219", "FAIL": "#d03b3b",
                "CTRL": "#898781"}
_INK = "#0b0b0b"
_INK2 = "#52514e"
_MUTED = "#898781"
_EMPTY = "#e1e0d9"
_GUIDE = "#d03b3b"


def _plot_bins(prof: DemuxProfile | None, kind: str) -> dict[int, int]:
    """Histogram clamped into the plot range, keyed by bin start."""
    if prof is None:
        return {}
    if kind == "len":
        lo, hi, src = LEN_PLOT_MIN, LEN_PLOT_MAX, prof.len_bins
    else:
        lo, hi, src = Q_PLOT_MIN, Q_PLOT_MAX, prof.q_win
    out: dict[int, int] = {}
    for b, n in src.items():
        k = min(max(b, lo), hi)
        out[k] = out.get(k, 0) + n
    return out


def _edges(kind: str) -> list[int]:
    if kind == "len":
        return list(range(LEN_PLOT_MIN, LEN_PLOT_MAX + BIN_WIDTH, BIN_WIDTH))
    return list(range(Q_PLOT_MIN, Q_PLOT_MAX + 1))


def _matrix(rows: list[SampleResult], profiles: dict[str, DemuxProfile],
            kind: str) -> tuple[Any, list[int]]:
    import numpy as np
    edges = _edges(kind)
    step = edges[1] - edges[0]
    mat = np.full((len(rows), len(edges)), np.nan)
    for i, r in enumerate(rows):
        bins = _plot_bins(profiles.get(r.sample.sample_id), kind)
        tot = sum(bins.values())
        if not tot:
            continue
        mat[i, :] = 0.0
        for b, n in bins.items():
            mat[i, (b - edges[0]) // step] += n / tot
    return mat, edges


def _guides(ax: Any, cfg: dict[str, Any], kind: str) -> None:
    """Length: shaded core window, dashed lines at the delivered window ends.
    Q is informational, so it gets no threshold line."""
    f, cw = cfg["filter"], cfg["core_window"]
    if kind == "len":
        ax.axvspan(cw["min"], cw["max"], color=_INK, alpha=0.08, lw=0)
        for x in (f["minknow_min_len"], f["chopper_max_len"]):
            ax.axvline(x, color=_GUIDE, lw=1.0, ls="--")


def _xlabel(cfg: dict[str, Any], kind: str) -> str:
    f, cw = cfg["filter"], cfg["core_window"]
    if kind == "len":
        return (f"Read length (bp, {BIN_WIDTH} bp bins; end bins hold "
                f"< {LEN_PLOT_MIN} / >= {LEN_PLOT_MAX}). Band = core window "
                f"{cw['min']}-{cw['max']} bp; dashed = MinKNOW min "
                f"{f['minknow_min_len']} bp and Chopper max "
                f"{f['chopper_max_len']} bp")
    return (f"Mean read Q of reads in {f['minknow_min_len']}-"
            f"{f['chopper_max_len']} bp (integer bins; end bins clamped). "
            f"Informational only -- read quality is not flagged")


def plot_heatmap(path: Path, rows: list[SampleResult],
                 profiles: dict[str, DemuxProfile], cfg: dict[str, Any],
                 kind: str, title: str) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, ListedColormap
    from matplotlib.patches import Patch
    import numpy as np

    mat, edges = _matrix(rows, profiles, kind)
    step = edges[1] - edges[0]
    n = len(rows)
    row_h = 0.11
    height = max(3.5, 1.6 + n * row_h)
    fig = plt.figure(figsize=(11, height), facecolor="#fcfcfb")
    gs = fig.add_gridspec(1, 3, width_ratios=[0.25, 10, 0.25], wspace=0.03,
                          left=0.16, right=0.94,
                          top=1 - 0.9 / height, bottom=0.75 / height)
    ax_flag = fig.add_subplot(gs[0])
    ax = fig.add_subplot(gs[1], sharey=ax_flag)
    ax_cb = fig.add_subplot(gs[2])

    flag_names = list(_FLAG_COLOUR)
    strip = np.array([[flag_names.index(r.flag)] for r in rows])
    ax_flag.imshow(strip, aspect="auto", interpolation="nearest",
                   cmap=ListedColormap([_FLAG_COLOUR[f] for f in flag_names]),
                   vmin=0, vmax=len(flag_names) - 1)
    ax_flag.set_xticks([])
    labelled = [(i, r) for i, r in enumerate(rows) if r.flag in ("CHECK", "FAIL")]
    ax_flag.set_yticks([i for i, _ in labelled])
    ax_flag.set_yticklabels([r.sample.sample_id for _, r in labelled],
                            fontsize=6, color=_INK)
    ax_flag.tick_params(axis="y", length=0)
    for sp in ax_flag.spines.values():
        sp.set_visible(False)

    cmap = LinearSegmentedColormap.from_list("seq", _SEQ)
    cmap.set_bad(_EMPTY)
    extent = (edges[0], edges[-1] + step, n - 0.5, -0.5)
    im = ax.imshow(np.ma.masked_invalid(mat), aspect="auto",
                   interpolation="nearest", cmap=cmap, vmin=0,
                   vmax=max(0.05, float(np.nanmax(mat)) if np.isfinite(mat).any() else 1),
                   extent=extent)
    _guides(ax, cfg, kind)
    ax.set_xlabel(_xlabel(cfg, kind), color=_INK2, fontsize=8)
    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.tick_params(axis="x", colors=_MUTED, labelsize=7)
    for sp in ax.spines.values():
        sp.set_color("#c3c2b7")

    cb = fig.colorbar(im, cax=ax_cb)
    cb.set_label("Fraction of the sample's reads", color=_INK2, fontsize=7)
    cb.ax.tick_params(labelsize=6, colors=_MUTED)
    cb.outline.set_visible(False)

    handles = [Patch(color=_FLAG_COLOUR[k], label=k) for k in flag_names]
    handles.append(Patch(color=_EMPTY, label="no data"))
    fig.legend(handles=handles, loc="upper right", ncol=len(handles),
               fontsize=7, frameon=False, bbox_to_anchor=(0.94, 1 - 0.1 / height))
    fig.text(0.16, 1 - 0.35 / height, title, fontsize=10, color=_INK,
             ha="left", va="top")
    fig.savefig(path, dpi=120, facecolor=fig.get_facecolor())
    plt.close(fig)


def plot_flagged(path: Path, r: SampleResult, prof: DemuxProfile | None,
                 cfg: dict[str, Any]) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.4), facecolor="#fcfcfb")
    colour = "#2a78d6"
    for a, kind in zip(axes, ("len", "q")):
        bins = _plot_bins(prof, kind)
        tot = sum(bins.values())
        edges = _edges(kind)
        step = edges[1] - edges[0]
        if tot:
            a.bar(list(bins), [c / tot for c in bins.values()],
                  width=step * 0.9, align="edge", color=colour, lw=0)
        else:
            a.text(0.5, 0.5, "no data", transform=a.transAxes,
                   ha="center", color=_MUTED)
        _guides(a, cfg, kind)
        a.set_xlim(edges[0], edges[-1] + step)
        a.set_xlabel("Read length (bp)" if kind == "len"
                     else f"Mean read Q (reads in {cfg['filter']['minknow_min_len']}-"
                          f"{cfg['filter']['chopper_max_len']} bp)",
                     fontsize=8, color=_INK2)
        a.set_ylabel("Fraction of reads", fontsize=8, color=_INK2)
        a.tick_params(colors=_MUTED, labelsize=7)
        a.grid(axis="y", color="#e1e0d9", lw=0.6)
        a.set_axisbelow(True)
        for side in ("top", "right"):
            a.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            a.spines[side].set_color("#c3c2b7")
    fig.suptitle(f"{r.sample.sample_id} ({r.sample.client}), n={r.n:,}"
                 f"  |  {r.reasons}", fontsize=9, color=_INK, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor=fig.get_facecolor())
    plt.close(fig)


def write_plots(qc_dir: Path, ev: Evaluation,
                profiles: dict[str, DemuxProfile], cfg: dict[str, Any],
                warnings: list[str], png_drilldown: bool = False) -> None:
    """Static PNG copies for email / records; qc_report.html is the main
    review view. Per-sample PNGs only with --png-drilldown."""
    for old in qc_dir.glob("*_heatmap*.png"):
        old.unlink()
    # Stale drill-downs from an earlier run would no longer match the flags.
    flagged_dir = qc_dir / "flagged"
    if flagged_dir.is_dir():
        shutil.rmtree(flagged_dir)

    try:
        import matplotlib
        matplotlib.use("Agg")
    except ImportError:
        warnings.append("matplotlib not available; static PNG heatmaps were "
                        "not drawn (qc_report.html is unaffected)")
        return

    ordered = sorted(ev.results, key=lambda r: (FLAG_ORDER[r.flag],
                                                r.sample.order, r.sample.barcode))
    if len(ordered) <= HEATMAP_MAX_ROWS:
        groups = [("", ordered)]
    else:
        clients = sorted({r.sample.client for r in ordered})
        groups = [(c, [r for r in ordered if r.sample.client == c])
                  for c in clients]
    when = ("raw demultiplexed, before Chopper" if ev.filter_stage == "post"
            else "pre mode: after Chopper")
    for client, rows in groups:
        suffix = f"_{client}" if client else ""
        for kind, stem, label in (("len", "length_heatmap", "Read length"),
                                  ("q", "q_heatmap",
                                   "Mean read Q, informational")):
            title = (f"{label} by sample ({when})"
                     + (f" -- {client}" if client else ""))
            plot_heatmap(qc_dir / f"{stem}{suffix}.png", rows, profiles, cfg,
                         kind, title)

    flagged = [r for r in ev.results if r.flag in ("CHECK", "FAIL")]
    if png_drilldown and flagged:
        flagged_dir.mkdir(parents=True, exist_ok=True)
        for r in flagged:
            plot_flagged(flagged_dir / f"{r.sample.sample_id}.png", r,
                         profiles.get(r.sample.sample_id), cfg)


# =============================================================================
# CLI
# =============================================================================
def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--output-dir", type=Path, required=True,
                   help="Pipeline run output dir; the report goes to "
                        "<output-dir>/qc_report/.")
    p.add_argument("--samplesheet", type=Path, required=True,
                   help="Sample sheet CSV (Client, Sample_ID, Barcode, optional "
                        "Location, Comment, Sample_Type).")
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG,
                   help=f"Thresholds (default: {DEFAULT_CONFIG.name} beside "
                        f"this script).")
    p.add_argument("--counts", type=Path, default=None,
                   help="read_counts_summary.txt (default: newest under "
                        "<output-dir>/run_log_*/, else <output-dir>/).")
    p.add_argument("--demux-qc-dir", type=Path, default=None,
                   help="demux_qc.py output (default: <output-dir>/demux_qc).")
    p.add_argument("--png-drilldown", action="store_true",
                   help="Also write flagged/<Sample>.png (length + Q "
                        "histograms) for CHECK / FAIL samples. qc_report.html "
                        "has the same drill-down interactively.")
    args = p.parse_args(argv)
    if args.demux_qc_dir is None:
        args.demux_qc_dir = args.output_dir / "demux_qc"
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    cfg = load_config(args.config)
    qc_dir = args.output_dir / "qc_report"
    warnings: list[str] = list(cfg["_warnings"])

    if not args.output_dir.is_dir():
        sys.exit(f"ERROR: output dir not found: {args.output_dir}")
    banner("QC REPORT: per-sample PASS / CHECK / FAIL")
    log(f" QC config: {cfg['_path']}  version {config_label(cfg)}")
    for w in cfg["_warnings"]:
        log(f" WARNING: {w}")

    sheet = load_samplesheet(args.samplesheet, cfg["controls"]["types"])
    if args.counts is not None:
        counts_path = args.counts
    else:
        counts_path, w = find_counts_file(args.output_dir)
        warnings += w
    counts = parse_read_counts(counts_path)
    log(f" Counts:    {counts_path} ({len(counts.rows)} rows)")
    profiles = load_demux_qc(args.demux_qc_dir,
                             [s.sample_id for s in sheet], warnings)
    log(f" Demux QC:  {args.demux_qc_dir} ({len(profiles)} profiles)")

    sheet_ids = {s.sample_id for s in sheet}
    extra = sorted(set(counts.rows) - sheet_ids - set(CATCH_ALL))
    if extra:
        warnings.append(f"read_counts_summary.txt rows not in the sample "
                        f"sheet: {', '.join(extra)}")

    # A --diagnostic run uses a wide Chopper window, so pct_retained no
    # longer means "> 2 kb removed".
    f = cfg["filter"]
    want = f"{f['minknow_min_len']}-{f['chopper_max_len']}"
    windows = {p.length_window for p in profiles.values() if p.length_window}
    if windows and windows != {want}:
        warnings.append(
            f"demux_qc length window {', '.join(sorted(windows))} differs from "
            f"the config's {want} (--diagnostic run?); long_frac reflects the "
            f"run's own Chopper window")

    ev = evaluate(sheet, counts, profiles, cfg, warnings)

    qc_dir.mkdir(parents=True, exist_ok=True)
    header = tsv_header(cfg, ev)
    write_sample_tsv(qc_dir / "per_sample_qc.tsv", ev.results, header)
    client_dir = qc_dir / "per_client"
    if client_dir.is_dir():
        shutil.rmtree(client_dir)
    control_types = cfg["controls"]["types"]
    clients = report_clients(sheet, control_types)
    stems = {c: client_file_stem(c) for c in clients}
    clash = sorted(c for c in clients if list(stems.values()).count(stems[c]) > 1)
    if clash:
        warnings.append(f"clients {', '.join(clash)} share a per-client file "
                        f"name after sanitising; only the last one is kept")
    for client in clients:
        write_sample_tsv(client_dir / f"{stems[client]}_qc.tsv",
                         [r for r in ev.results if r.sample.client == client
                          or is_shared(r.sample, control_types)],
                         header)
    write_bins(qc_dir / "qc_bins.tsv.gz", profiles)
    write_plots(qc_dir, ev, profiles, cfg, warnings, args.png_drilldown)
    write_html(qc_dir / "qc_report.html",
               html_data(args, cfg, ev, profiles, warnings))
    # Per-client HTML goes into the client's delivery folder (the folder name
    # is the sanitised client, as organise_by_client made it), so it ships
    # with the reads. Without integrated_demultiplexing/ it stays in
    # per_client/.
    integ = args.output_dir / "integrated_demultiplexing"
    if clients and not integ.is_dir():
        warnings.append(f"no {integ.name}/ in the run dir; per-client HTML "
                        f"reports written to qc_report/per_client/ instead")
    for client in clients:
        dest = integ / client if integ.is_dir() else client_dir
        dest.mkdir(parents=True, exist_ok=True)
        write_html(dest / f"{stems[client]}_qc_report.html",
                   html_data(args, cfg, ev, profiles, warnings, client=client))
    write_run_qc(qc_dir / "run_qc.txt", args, cfg, counts, ev, profiles,
                 warnings)

    log(f" Filter stage: {ev.filter_stage}")
    for ln in run_flag_lines(ev.runs):
        log(f" {ln}")
    for flag in ("PASS", "CHECK", "FAIL", "CTRL"):
        log(f" {flag:<6} {sum(r.flag == flag for r in ev.results):>5}")
    if warnings:
        log(f" {len(warnings)} warning(s); see run_qc.txt")
    log(f" Report: {qc_dir / 'qc_report.html'}")
    if clients:
        where = ("integrated_demultiplexing/<Client>/" if integ.is_dir()
                 else "qc_report/per_client/")
        log(f" Per-client reports ({len(clients)}): {where}<Client>_qc_report.html")
    log(f" QC report complete: {datetime.now():%Y-%m-%d %H:%M:%S}")
    log("========================================")
    return 0


if __name__ == "__main__":
    sys.exit(main())
