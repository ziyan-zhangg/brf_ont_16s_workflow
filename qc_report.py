#!/usr/bin/env python3
"""
qc_report.py

Standalone per-sample QC report (PASS / CHECK / FAIL) for a finished 16S run.

Judges only what the facility is responsible for: did each sample amplify and
sequence into usable reads. Chimera %, classification and NTC read ratios are
the client's responsibility and are deliberately not reported.

Runs as its own PBS job after the production pipeline (run_qc_report.qsub,
chained with qsub -W depend=afterok:<jid>). The production pipeline and its
outputs are never modified.

Where the numbers come from:

    exact      delivered_reads, primer_found
               parsed from <output_dir>/run_log_*/read_counts_summary.txt
    estimated  retention, short_frac, long_frac, lowq_frac
               from a diagnostic slice: Minibar (production flags) on every
               Nth raw fastq_pass file, so the reads are seen *before*
               Chopper. Production runs in filter_stage=pre and never sees
               them unfiltered.

Steps:

    1. Select slice   every Nth raw fastq_pass file, natural name order
    2. Minibar        on the slice, in parallel, work dirs on $PBS_JOBFS
    3. Scan           per sample: 50 bp length histogram, integer mean-Q
                      histogram, short / long / lowQ / pass categories
    4. Cross-check    Chopper (production settings) on the slice raw files;
                      its pass fraction must agree with the slice categories
    5. Counts         read_counts_summary.txt (per-sample rows only)
    6. Flag + write   <output_dir>/qc_report/

--report-only skips 1-4 and rebuilds everything from the saved
qc_report/qc_bins.tsv.gz, so thresholds can be changed in qc_config.yaml and
re-applied without touching reads.

PASS / CHECK / FAIL is triage, not sign-off: staff review run_qc.txt and the
heatmaps before release.

Usage:
    python3 qc_report.py \\
        --raw-input-dir /path/to/fastq_pass \\
        --output-dir    /path/to/run_output \\
        --primer-file   /path/to/16S_primer_setup_YYYYMMDD.txt \\
        --samplesheet   /path/to/16s_samplesheet.csv \\
        [--config qc_config.yaml] [--work-dir $PBS_JOBFS/qc_work]

    python3 qc_report.py --report-only \\
        --output-dir /path/to/run_output --samplesheet /path/to/sheet.csv
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import math
import os
import re
import shutil
import statistics
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, BinaryIO

# Local imports: sibling modules in the same tools/ dir.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import filter_chopper_demux_minibar as core  # noqa: E402
from demux_qc import (  # noqa: E402
    BIN_WIDTH,
    CATCH_ALL,
    _ERR_PROB,
    _open_fastq,
    sample_id_from_path,
)

log = core.log
banner = core.banner

DEFAULT_CONFIG = Path(__file__).resolve().parent / "qc_config.yaml"

# Mutually exclusive slice categories, in the order they are tested.
CATEGORIES = ("short", "long", "lowq", "pass")

# Metric order is also the tie-break order for reason codes.
METRICS = ("delivered_reads", "primer_found", "retention",
           "short_frac", "long_frac", "lowq_frac")
ESTIMATED = ("retention", "short_frac", "long_frac", "lowq_frac")
REASON_CODE = {
    "delivered_reads": "READS", "primer_found": "PRIMER",
    "retention": "RETENTION", "short_frac": "SHORT",
    "long_frac": "LONG", "lowq_frac": "LOWQ",
}
RUN_CODE = {
    "delivered_reads": "LOW_READS", "primer_found": "LOW_PRIMER",
    "retention": "LOW_RETENTION", "short_frac": "SHORT",
    "long_frac": "LONG", "lowq_frac": "LOW_Q",
}
LEVEL_NAME = {0: "PASS", 1: "CHECK", 2: "FAIL"}
FLAG_ORDER = {"FAIL": 0, "CHECK": 1, "PASS": 2, "CTRL": 3}

# Pseudo-sample holding run metadata inside qc_bins.tsv.gz.
META_SAMPLE = "__run__"

# MAD -> sigma for normally distributed data, so k_check/k_fail read as
# "standard deviations from the run median".
MAD_TO_SIGMA = 1.4826

# Plot ranges. Reads beyond the top bin are clamped into it.
LEN_PLOT_MAX = 4000
Q_PLOT_MAX = 40
HEATMAP_MAX_ROWS = 384


# =============================================================================
# Config
# =============================================================================
def load_config(path: Path) -> dict[str, Any]:
    try:
        import yaml
    except ImportError:
        sys.exit("ERROR: PyYAML is required to read the QC config. On Gadi:\n"
                 "  module load python3/3.12.1 && pip install --user pyyaml")
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

    if "config_version" not in cfg:
        sys.exit(f"ERROR: {path}: missing config_version")
    need("filter", "min_len", "max_len", "min_q")
    need("diagnostic_subsample", "every_nth_file", "min_reads_for_stats",
         "crosscheck_tolerance")
    need("metrics", *METRICS)
    need("run_level", "outlier_k_check", "outlier_k_fail", "mad_floor")
    need("controls", "types")

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
    if cfg["diagnostic_subsample"]["every_nth_file"] < 1:
        sys.exit(f"ERROR: {path}: every_nth_file must be >= 1")

    cfg["run_level"].setdefault("min_samples", 8)
    cfg["controls"]["types"] = [str(t).upper() for t in cfg["controls"]["types"]]
    cfg["_path"] = str(path)
    cfg["_sha256"] = hashlib.sha256(raw).hexdigest()[:12]
    return cfg


def ordinal(n: Any) -> str:
    """10 -> '10th', 2 -> '2nd'; anything non-integer passes through."""
    if not isinstance(n, int):
        return str(n)
    suffix = "th" if 10 <= n % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th")
    return f"{n}{suffix}"


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
    match the per-client folders and the sample_<id>.fastq names Minibar wrote.
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
            if not client or not sid:
                continue
            stype = (row.get("Sample_Type") or "").strip().upper() or "SAMPLE"
            if stype not in allowed:
                sys.exit(f"ERROR: {path} line {i}: Sample_Type {stype!r} is not "
                         f"one of {', '.join(sorted(allowed))}")
            out.append(SheetSample(
                client=core.sanitise(client),
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


# =============================================================================
# read_counts_summary.txt
# =============================================================================
@dataclass
class CountsRow:
    sample: str
    reads: int
    post_chopper: int | None
    pass1: int | None
    pct_pass1_kept: float | None   # percent, as printed


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


def parse_read_counts(path: Path) -> dict[str, CountsRow]:
    """
    Per-sample rows of read_counts_summary.txt, keyed by sample ID (no
    sample_ prefix). The TOTAL rows below the second separator are ignored.

    The table is fixed-width with right-aligned numeric cells, and several
    cells are legitimately blank (raw_demuxed and pct_retained in pre mode),
    so whitespace splitting cannot tell which column a number belongs to.
    Each cell is instead cut at the end position of its header label, shifted
    by however far an over-long sample name pushed the row to the right.
    """
    if not path.is_file():
        sys.exit(f"ERROR: read counts summary not found: {path}")
    lines = path.read_text().splitlines()

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
            post_chopper=None if cells["post_chopper"] is None
            else int(cells["post_chopper"]),
            pass1=None if cells["post_cutadapt_pass1"] is None
            else int(cells["post_cutadapt_pass1"]),
            pct_pass1_kept=cells["pct_pass1_kept"],
        )
    if not rows:
        sys.exit(f"ERROR: {path}: no per-sample rows found under the header")
    return rows


# =============================================================================
# Slice profiling
# =============================================================================
@dataclass
class SliceProfile:
    sample_id: str
    len_bins: dict[int, int] = field(default_factory=dict)  # 50 bp bin start
    q_bins: dict[int, int] = field(default_factory=dict)    # floor(mean Q)
    cats: dict[str, int] = field(
        default_factory=lambda: {c: 0 for c in CATEGORIES})

    @property
    def n(self) -> int:
        return sum(self.cats.values())

    def frac(self, cat: str) -> float:
        return self.cats.get(cat, 0) / self.n if self.n else 0.0

    def merge(self, other: "SliceProfile") -> None:
        for mine, theirs in ((self.len_bins, other.len_bins),
                             (self.q_bins, other.q_bins),
                             (self.cats, other.cats)):
            for k, v in theirs.items():
                mine[k] = mine.get(k, 0) + v


def scan_slice_fastq(path: Path, sample_id: str, min_len: int, max_len: int,
                     min_q: float) -> SliceProfile:
    """
    One streaming pass: length histogram, mean-Q histogram, categories.

    Mean read Q is demux_qc's definition (mean per-base error probability ->
    Phred), which is also how Chopper computes the quality it filters on.
    Categories are exclusive and tested in order:
        length < min_len -> short; length > max_len -> long;
        mean Q < min_q   -> lowq;  else pass.
    """
    prof = SliceProfile(sample_id)
    lens, qs, cats = prof.len_bins, prof.q_bins, prof.cats
    err_prob = _ERR_PROB.__getitem__
    with _open_fastq(path) as fh:
        readline = fh.readline
        while True:
            if not readline():
                break
            seq = readline().rstrip()
            plus = readline()
            qual = readline()
            if not plus or not qual:
                break                      # truncated final record
            qual = qual.rstrip()

            length = len(seq)
            b = (length // BIN_WIDTH) * BIN_WIDTH
            lens[b] = lens.get(b, 0) + 1

            q = None
            if qual:
                mean_err = sum(map(err_prob, qual)) / len(qual)
                q = -10.0 * math.log10(mean_err) if mean_err > 0 else 93.0
                qb = int(q)
                qs[qb] = qs.get(qb, 0) + 1

            if length < min_len:
                cats["short"] += 1
            elif length > max_len:
                cats["long"] += 1
            elif q is None or q < min_q:
                cats["lowq"] += 1
            else:
                cats["pass"] += 1
    return prof


def natural_key(path: Path) -> list[Any]:
    """'x_10.fastq.gz' after 'x_9.fastq.gz', so every Nth file is spread
    evenly over the run instead of clustering on lexical runs of digits."""
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", path.name)]


def select_slice(raw_input_dir: Path, every_nth: int) -> tuple[list[Path], int]:
    inputs = sorted(raw_input_dir.glob("*.fastq.gz"), key=natural_key)
    if not inputs:
        sys.exit(f"ERROR: no *.fastq.gz files found in {raw_input_dir}")
    return inputs[::every_nth], len(inputs)


def _scan_chunk_dir(chunk_dir: Path, min_len: int, max_len: int,
                    min_q: float) -> dict[str, SliceProfile]:
    """Worker: profile every sample_*.fastq Minibar wrote for one chunk."""
    out = {}
    for fq in sorted(chunk_dir.glob("sample_*.fastq")):
        sid = sample_id_from_path(fq)
        out[sid] = scan_slice_fastq(fq, sid, min_len, max_len, min_q)
    return out


def profile_slice(chunk_dirs: list[Path], cfg: dict[str, Any],
                  workers: int) -> dict[str, SliceProfile]:
    """Scan every chunk's per-sample output and merge across chunks.

    Pure-Python per-base work, so it uses processes rather than threads.
    """
    banner("QC STEP 3: Profiling slice reads")
    f = cfg["filter"]
    profiles: dict[str, SliceProfile] = {}
    n = max(1, min(workers, len(chunk_dirs)))
    with concurrent.futures.ProcessPoolExecutor(max_workers=n) as ex:
        futs = [ex.submit(_scan_chunk_dir, d, f["min_len"], f["max_len"],
                          f["min_q"]) for d in chunk_dirs]
        for fut in futs:
            for sid, prof in fut.result().items():
                profiles.setdefault(sid, SliceProfile(sid)).merge(prof)
    total = sum(p.n for p in profiles.values())
    log(f" {len(profiles)} bin(s), {total:,} slice reads profiled")
    log("========================================")
    return profiles


# =============================================================================
# Chopper cross-check
# =============================================================================
class _CountingReader:
    """File-like wrapper that counts newlines as Chopper's feeder reads."""

    def __init__(self, fh: BinaryIO) -> None:
        self.fh = fh
        self.lines = 0

    def read(self, n: int = -1) -> bytes:
        data = self.fh.read(n)
        self.lines += data.count(b"\n")
        return data


class _CountingSink:
    """Discards Chopper's output, keeping only its line count."""

    def __init__(self) -> None:
        self.lines = 0

    def write(self, data: bytes) -> int:
        self.lines += data.count(b"\n")
        return len(data)


def _chopper_count(raw: Path, chopper: Path, f: dict[str, Any]) -> tuple[int, int]:
    with gzip.open(raw, "rb") as zin:
        src, dst = _CountingReader(zin), _CountingSink()
        core._chopper_stream(chopper, src, dst, f["min_q"], f["min_len"],
                             f["max_len"], what=str(raw))
    return src.lines // 4, dst.lines // 4


def run_crosscheck(slice_files: list[Path], chopper: Path,
                   cfg: dict[str, Any], workers: int) -> tuple[int, int]:
    """Production Chopper settings on the slice raw files -> (reads_in,
    reads_out). Chopper is itself multi-threaded, so fewer workers."""
    banner("QC STEP 4: Chopper cross-check on slice raw files")
    n = max(1, min(workers // 4 or 1, len(slice_files)))
    with concurrent.futures.ThreadPoolExecutor(max_workers=n) as ex:
        res = list(ex.map(lambda p: _chopper_count(p, chopper, cfg["filter"]),
                          slice_files))
    n_in = sum(a for a, _ in res)
    n_out = sum(b for _, b in res)
    log(f" Chopper kept {n_out:,} of {n_in:,} slice reads")
    log("========================================")
    return n_in, n_out


# =============================================================================
# qc_bins.tsv.gz
# =============================================================================
def write_bins(path: Path, profiles: dict[str, SliceProfile],
               meta: dict[str, int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["sample", "kind", "bin", "count"])
        for k, v in meta.items():
            w.writerow([META_SAMPLE, "meta", k, v])
        for sid in sorted(profiles):
            p = profiles[sid]
            for c in CATEGORIES:
                w.writerow([sid, "category", c, p.cats.get(c, 0)])
            for b in sorted(p.len_bins):
                w.writerow([sid, "len", b, p.len_bins[b]])
            for b in sorted(p.q_bins):
                w.writerow([sid, "q", b, p.q_bins[b]])


def read_bins(path: Path) -> tuple[dict[str, SliceProfile], dict[str, int]]:
    if not path.is_file():
        sys.exit(f"ERROR: {path} not found -- --report-only needs a previous "
                 f"full run of qc_report.py on this output dir")
    profiles: dict[str, SliceProfile] = {}
    meta: dict[str, int] = {}
    with gzip.open(path, "rt", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sid, kind, b, n = row["sample"], row["kind"], row["bin"], int(row["count"])
            if kind == "meta":
                meta[b] = n
                continue
            p = profiles.setdefault(sid, SliceProfile(sid))
            if kind == "category":
                p.cats[b] = n
            elif kind == "len":
                p.len_bins[int(b)] = n
            elif kind == "q":
                p.q_bins[int(b)] = n
            else:
                sys.exit(f"ERROR: {path}: unknown kind {kind!r}")
    return profiles, meta


# =============================================================================
# Flagging
# =============================================================================
@dataclass
class SampleResult:
    sample: SheetSample
    slice_n: int
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


def abs_level(value: float, spec: dict[str, Any]) -> int:
    if spec["direction"] == "higher_better":
        return 2 if value < spec["fail"] else 1 if value < spec["check"] else 0
    return 2 if value > spec["fail"] else 1 if value > spec["check"] else 0


def fmt_value(metric: str, value: float) -> str:
    if metric == "delivered_reads":
        return f"{int(value):,}"
    return f"{100 * value:.0f}%"


def evaluate(sheet: list[SheetSample], counts: dict[str, CountsRow],
             profiles: dict[str, SliceProfile], cfg: dict[str, Any],
             warnings: list[str]) -> tuple[list[SampleResult], list[MetricRun]]:
    min_n = cfg["diagnostic_subsample"]["min_reads_for_stats"]
    rl = cfg["run_level"]
    controls = set(cfg["controls"]["types"])

    results: list[SampleResult] = []
    for s in sheet:
        row = counts.get(s.sample_id)
        prof = profiles.get(s.sample_id)
        n = prof.n if prof else 0
        values: dict[str, float | None] = {}
        na: dict[str, str] = {}

        if row is None:
            warnings.append(f"{s.sample_id}: not in read_counts_summary.txt; "
                            f"delivered_reads taken as 0")
        values["delivered_reads"] = row.reads if row else 0
        if row is not None and row.pct_pass1_kept is not None:
            values["primer_found"] = row.pct_pass1_kept / 100.0
        else:
            values["primer_found"] = None
            na["primer_found"] = "NA"

        if n < min_n:
            for m in ESTIMATED:
                values[m] = None
                na[m] = f"NA (n<{min_n})"
        else:
            values["retention"] = prof.frac("pass")
            values["short_frac"] = prof.frac("short")
            values["long_frac"] = prof.frac("long")
            values["lowq_frac"] = prof.frac("lowq")
        results.append(SampleResult(s, n, values, na))

    real = [r for r in results if r.sample.sample_type not in controls]

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
        elif r.slice_n < min_n:
            r.reasons = f"estimates NA (slice n={r.slice_n})"
    return results, runs


def run_flag_lines(runs: list[MetricRun], cfg: dict[str, Any]) -> list[str]:
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
_TSV_COLUMNS = ["client", "sample", "barcode", "location", "type", "slice_n",
                *METRICS, "flag", "reasons", "comment"]


def _cell(r: SampleResult, m: str) -> str:
    v = r.values[m]
    if v is None:
        return r.na.get(m, "NA")
    return str(int(v)) if m == "delivered_reads" else f"{v:.4f}"


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
            w.writerow([s.client, s.sample_id, s.barcode, s.location,
                        s.sample_type, r.slice_n,
                        *(_cell(r, m) for m in METRICS),
                        r.flag, r.reasons, s.comment])


def tsv_header(cfg: dict[str, Any], meta: dict[str, int],
               runs: list[MetricRun]) -> list[str]:
    lines = [
        f"generated\t{datetime.now():%Y-%m-%d %H:%M:%S}",
        f"config_version\t{config_label(cfg)}",
        f"slice\t{meta.get('slice_files', 0)} of {meta.get('total_files', 0)} "
        f"raw files (every {ordinal(meta.get('every_nth_file', '?'))})",
        f"exact_columns\tdelivered_reads, primer_found "
        f"(from read_counts_summary.txt)",
        f"estimated_columns\t{', '.join(ESTIMATED)} "
        f"(from the diagnostic slice, before Chopper)",
        "fractions\t0-1; NA (n<N) = too few slice reads to estimate",
    ]
    lines += [f"run_flag\t{ln}" for ln in run_flag_lines(runs, cfg)]
    return lines


def write_run_qc(path: Path, args: argparse.Namespace, cfg: dict[str, Any],
                 counts_path: Path, results: list[SampleResult],
                 runs: list[MetricRun], profiles: dict[str, SliceProfile],
                 meta: dict[str, int], warnings: list[str]) -> None:
    f = cfg["filter"]
    tol = cfg["diagnostic_subsample"]["crosscheck_tolerance"]
    total_slice = sum(p.n for p in profiles.values())
    slice_pass = sum(p.cats.get("pass", 0) for p in profiles.values())

    L = [
        "========================================",
        " Per-sample QC report -- run summary",
        f" Generated:  {datetime.now():%Y-%m-%d %H:%M:%S}",
        f" Output dir: {args.output_dir}",
        f" Counts:     {counts_path}",
        f" Config:     {cfg['_path']}",
        f" Version:    {config_label(cfg)}",
        f" Mode:       {'report-only (rebuilt from qc_bins.tsv.gz)' if args.report_only else 'full'}",
        "========================================",
        "",
        " Diagnostic slice",
        f"   files:   {meta.get('slice_files', 0)} of {meta.get('total_files', 0)} "
        f"raw fastq_pass files (every {ordinal(meta.get('every_nth_file', '?'))}, "
        f"natural name order)"
        + (f", {meta['empty_files']} empty skipped" if meta.get("empty_files") else ""),
        f"   reads:   {total_slice:,} demultiplexed (all bins)",
        f"   filter:  {meta.get('min_len', f['min_len'])}-"
        f"{meta.get('max_len', f['max_len'])} bp, mean Q >= "
        f"{meta.get('min_q', f['min_q'])}",
        "",
        " Chopper cross-check",
    ]
    if "chopper_in" in meta and meta["chopper_in"]:
        c_frac = meta["chopper_out"] / meta["chopper_in"]
        s_frac = slice_pass / total_slice if total_slice else 0.0
        diff = abs(c_frac - s_frac)
        verdict = "OK" if diff <= tol else "WARNING"
        L += [
            f"   Chopper pass fraction (slice raw):  {c_frac:.4f}  "
            f"({meta['chopper_out']:,} / {meta['chopper_in']:,})",
            f"   slice category pass fraction:       {s_frac:.4f}  "
            f"({slice_pass:,} / {total_slice:,})",
            f"   difference: {diff:.4f} (tolerance {tol}) -> {verdict}",
        ]
        if verdict == "WARNING":
            warnings.append(
                f"Chopper cross-check differs by {diff:.4f} (> {tol}): the "
                f"slice categories may not mirror production Chopper; check "
                f"the filter section of the config against run_script.qsub")
        if meta["chopper_in"] != total_slice:
            warnings.append(
                f"slice raw files hold {meta['chopper_in']:,} reads but Minibar "
                f"binned {total_slice:,}; estimates may be biased")
    else:
        L.append("   not run")

    L += ["", " Run-level flags"]
    flags = run_flag_lines(runs, cfg)
    L += [f"   {ln}" for ln in flags] or ["   none"]

    L += ["", " Run medians (real samples, controls excluded)",
          f"   {'metric':<17} {'median':>9} {'MAD':>9} {'n':>5}  mode"]
    for mr in runs:
        if mr.median is None:
            L.append(f"   {mr.metric:<17} {'NA':>9} {'NA':>9} {0:>5}  -")
            continue
        fmt = (lambda v: f"{v:,.0f}") if mr.metric == "delivered_reads" \
            else (lambda v: f"{v:.4f}")
        L.append(f"   {mr.metric:<17} {fmt(mr.median):>9} {fmt(mr.mad):>9} "
                 f"{mr.n:>5}  {'relative' if mr.relative else 'absolute'}")

    L += ["", " Flag counts"]
    for flag in ("PASS", "CHECK", "FAIL", "CTRL"):
        L.append(f"   {flag:<6} {sum(r.flag == flag for r in results):>5}")

    L += ["", " Catch-all bins (slice)"]
    for sid in CATCH_ALL:
        p = profiles.get(sid)
        if not p or not p.n:
            L.append(f"   {sid:<17} 0 reads")
            continue
        share = 100.0 * p.n / total_slice if total_slice else 0.0
        L.append(f"   {sid:<17} {p.n:,} reads ({share:.1f}% of slice): "
                 + ", ".join(f"{c} {100 * p.frac(c):.1f}%" for c in CATEGORIES))

    known = {s.sample.sample_id for s in results} | set(CATCH_ALL)
    stray = sorted(sid for sid, p in profiles.items() if sid not in known and p.n)
    if stray:
        warnings.append(f"slice bins not in the sample sheet: {', '.join(stray)}")

    L += ["", " Warnings"]
    L += [f"   - {w}" for w in warnings] or ["   none"]
    L += [
        "",
        " PASS / CHECK / FAIL is triage, not sign-off. Review this file and",
        " length_heatmap.png / q_heatmap.png before release.",
        "========================================",
    ]
    path.write_text("\n".join(L) + "\n")


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


def _matrix(rows: list[SampleResult], profiles: dict[str, SliceProfile],
            kind: str) -> tuple[Any, list[int]]:
    import numpy as np
    if kind == "len":
        edges = list(range(0, LEN_PLOT_MAX + BIN_WIDTH, BIN_WIDTH))
    else:
        edges = list(range(0, Q_PLOT_MAX + 1))
    top = edges[-1]
    step = edges[1] - edges[0]
    mat = np.full((len(rows), len(edges)), np.nan)
    for i, r in enumerate(rows):
        p = profiles.get(r.sample.sample_id)
        bins = (p.len_bins if kind == "len" else p.q_bins) if p else {}
        tot = sum(bins.values())
        if not tot:
            continue
        mat[i, :] = 0.0
        for b, n in bins.items():
            mat[i, min(b, top) // step] += n / tot
    return mat, edges


def plot_heatmap(path: Path, rows: list[SampleResult],
                 profiles: dict[str, SliceProfile], cfg: dict[str, Any],
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
    f = cfg["filter"]
    if kind == "len":
        ax.axvspan(f["min_len"], f["max_len"] + 1, color=_INK, alpha=0.06, lw=0)
        for x in (f["min_len"], f["max_len"] + 1):
            ax.axvline(x, color=_INK2, lw=0.8)
        ax.set_xlabel(f"Read length (bp, {BIN_WIDTH} bp bins; "
                      f">= {LEN_PLOT_MAX} in last bin). "
                      f"Band = Chopper window {f['min_len']}-{f['max_len']} bp",
                      color=_INK2, fontsize=8)
    else:
        ax.axvline(f["min_q"], color=_INK2, lw=1.0, ls="--")
        ax.set_xlabel(f"Mean read Q (integer bins; >= {Q_PLOT_MAX} in last bin). "
                      f"Dashed = Chopper min Q{f['min_q']}",
                      color=_INK2, fontsize=8)
    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.tick_params(axis="x", colors=_MUTED, labelsize=7)
    for sp in ax.spines.values():
        sp.set_color("#c3c2b7")

    cb = fig.colorbar(im, cax=ax_cb)
    cb.set_label("Fraction of the sample's slice reads", color=_INK2, fontsize=7)
    cb.ax.tick_params(labelsize=6, colors=_MUTED)
    cb.outline.set_visible(False)

    handles = [Patch(color=_FLAG_COLOUR[k], label=k) for k in flag_names]
    handles.append(Patch(color=_EMPTY, label="no slice reads"))
    fig.legend(handles=handles, loc="upper right", ncol=len(handles),
               fontsize=7, frameon=False, bbox_to_anchor=(0.94, 1 - 0.1 / height))
    fig.text(0.16, 1 - 0.35 / height, title, fontsize=10, color=_INK,
             ha="left", va="top")
    fig.savefig(path, dpi=120, facecolor=fig.get_facecolor())
    plt.close(fig)


def plot_flagged(path: Path, r: SampleResult, prof: SliceProfile | None,
                 cfg: dict[str, Any]) -> None:
    import matplotlib.pyplot as plt

    f = cfg["filter"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 3.4), facecolor="#fcfcfb")
    colour = "#2a78d6"
    if prof and prof.n:
        # Same ranges and top-bin clamping as the heatmaps.
        lens: dict[int, int] = {}
        for b, c in prof.len_bins.items():
            lens[min(b, LEN_PLOT_MAX)] = lens.get(min(b, LEN_PLOT_MAX), 0) + c
        qs: dict[int, int] = {}
        for b, c in prof.q_bins.items():
            qs[min(b, Q_PLOT_MAX)] = qs.get(min(b, Q_PLOT_MAX), 0) + c
        a1.bar(list(lens), [c / prof.n for c in lens.values()],
               width=BIN_WIDTH * 0.9, align="edge", color=colour, lw=0)
        a2.bar(list(qs), [c / prof.n for c in qs.values()], width=0.9,
               align="edge", color=colour, lw=0)
    else:
        for a in (a1, a2):
            a.text(0.5, 0.5, "no slice reads", transform=a.transAxes,
                   ha="center", color=_MUTED)
    a1.axvspan(f["min_len"], f["max_len"] + 1, color=_INK, alpha=0.06, lw=0)
    a1.set_xlabel(f"Read length (bp); band = {f['min_len']}-{f['max_len']} bp",
                  fontsize=8, color=_INK2)
    a1.set_xlim(0, LEN_PLOT_MAX + BIN_WIDTH)
    a2.axvline(f["min_q"], color=_INK2, lw=1.0, ls="--")
    a2.set_xlabel(f"Mean read Q; dashed = Q{f['min_q']}", fontsize=8, color=_INK2)
    a2.set_xlim(0, Q_PLOT_MAX + 1)
    for a in (a1, a2):
        a.set_ylabel("Fraction of slice reads", fontsize=8, color=_INK2)
        a.tick_params(colors=_MUTED, labelsize=7)
        a.grid(axis="y", color="#e1e0d9", lw=0.6)
        a.set_axisbelow(True)
        for side in ("top", "right"):
            a.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            a.spines[side].set_color("#c3c2b7")
    fig.suptitle(f"{r.sample.sample_id} ({r.sample.client}), slice n={r.slice_n:,}"
                 f"  |  {r.reasons}", fontsize=9, color=_INK, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(path, dpi=110, facecolor=fig.get_facecolor())
    plt.close(fig)


def write_plots(qc_dir: Path, results: list[SampleResult],
                profiles: dict[str, SliceProfile], cfg: dict[str, Any],
                warnings: list[str]) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
    except ImportError:
        warnings.append("matplotlib not available; heatmaps and flagged "
                        "plots were not drawn (pip install --user matplotlib)")
        return

    ordered = sorted(results, key=lambda r: (FLAG_ORDER[r.flag], r.sample.order,
                                             r.sample.barcode))
    if len(ordered) <= HEATMAP_MAX_ROWS:
        groups = [("", ordered)]
    else:
        clients = sorted({r.sample.client for r in ordered})
        groups = [(c, [r for r in ordered if r.sample.client == c])
                  for c in clients]
    for client, rows in groups:
        suffix = f"_{client}" if client else ""
        for kind, stem, label in (("len", "length_heatmap", "Read length"),
                                  ("q", "q_heatmap", "Mean read Q")):
            title = (f"{label} by sample (diagnostic slice, before Chopper)"
                     + (f" -- {client}" if client else ""))
            plot_heatmap(qc_dir / f"{stem}{suffix}.png", rows, profiles, cfg,
                         kind, title)

    flagged_dir = qc_dir / "flagged"
    if flagged_dir.is_dir():
        shutil.rmtree(flagged_dir)
    flagged = [r for r in results if r.flag in ("CHECK", "FAIL")]
    if flagged:
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
    p.add_argument("--raw-input-dir", type=Path,
                   help="fastq_pass directory (not needed with --report-only).")
    p.add_argument("--primer-file", type=Path,
                   help="Primer setup file (not needed with --report-only).")
    p.add_argument("--config", type=Path, default=DEFAULT_CONFIG,
                   help=f"Thresholds (default: {DEFAULT_CONFIG.name} beside "
                        f"this script).")
    p.add_argument("--counts", type=Path, default=None,
                   help="read_counts_summary.txt (default: newest under "
                        "<output-dir>/run_log_*/).")
    p.add_argument("--work-dir", type=Path, default=None,
                   help="Scratch for the slice Minibar run (default: "
                        "$PBS_JOBFS/qc_work, else $TMPDIR/qc_work).")
    p.add_argument("--workers", type=int, default=None,
                   help="Parallel workers (default: $PBS_NCPUS, else 1).")
    p.add_argument("--chopper", type=Path, default=core.DEFAULT_CHOPPER)
    p.add_argument("--minibar", type=Path, default=core.DEFAULT_MINIBAR)
    p.add_argument("--no-crosscheck", action="store_true",
                   help="Skip the Chopper cross-check (step 4).")
    p.add_argument("--report-only", action="store_true",
                   help="Skip steps 1-4; rebuild from qc_report/qc_bins.tsv.gz.")
    p.add_argument("--keep-work", action="store_true",
                   help="Do not delete the work dir afterwards.")
    args = p.parse_args(argv)
    if not args.report_only:
        for name in ("raw_input_dir", "primer_file"):
            if getattr(args, name) is None:
                p.error(f"--{name.replace('_', '-')} is required unless "
                        f"--report-only is given")
    return args


def resolve_work_dir(args: argparse.Namespace, qc_dir: Path,
                     warnings: list[str]) -> Path:
    if args.work_dir is not None:
        return args.work_dir
    for var in ("PBS_JOBFS", "TMPDIR"):
        if os.environ.get(var):
            return Path(os.environ[var]) / "qc_work"
    warnings.append(f"PBS_JOBFS and TMPDIR unset; slice work dir is "
                    f"{qc_dir / '_work'}")
    return qc_dir / "_work"


def build_slice(args: argparse.Namespace, cfg: dict[str, Any], qc_dir: Path,
                warnings: list[str]) -> tuple[dict[str, SliceProfile], dict[str, int]]:
    """Steps 1-4. Returns the merged per-sample profiles and run metadata."""
    for label, pth, is_dir in (("raw input dir", args.raw_input_dir, True),
                               ("primer file", args.primer_file, False),
                               ("minibar.py", args.minibar, False)):
        if not (pth.is_dir() if is_dir else pth.is_file()):
            sys.exit(f"ERROR: {label} not found: {pth}")
    if not args.no_crosscheck and not args.chopper.is_file():
        sys.exit(f"ERROR: chopper not found: {args.chopper}")

    workers = args.workers or int(os.environ.get("PBS_NCPUS", "1") or "1")
    sub = cfg["diagnostic_subsample"]
    f = cfg["filter"]

    banner("QC STEP 1: Selecting diagnostic slice")
    picked, total = select_slice(args.raw_input_dir, sub["every_nth_file"])
    slice_files = [p for p in picked if core.fastq_has_reads(p)]
    empty = len(picked) - len(slice_files)
    log(f" {len(picked)} of {total} raw files (every {ordinal(sub['every_nth_file'])})"
        + (f", {empty} empty skipped" if empty else ""))
    log("========================================")
    if not slice_files:
        sys.exit("ERROR: every file in the diagnostic slice is empty")

    work_dir = resolve_work_dir(args, qc_dir, warnings)
    try:
        # QC STEP 2: same flags as production, via the production function.
        chunk_dirs = core.run_minibar(slice_files, work_dir / "minibar",
                                      args.primer_file, args.minibar,
                                      qc_dir / "logs" / "minibar",
                                      max_workers=workers)
        profiles = profile_slice(chunk_dirs, cfg, workers)
    finally:
        if not args.keep_work:
            shutil.rmtree(work_dir, ignore_errors=True)

    meta = {
        "every_nth_file": sub["every_nth_file"], "total_files": total,
        "slice_files": len(slice_files), "empty_files": empty,
        "min_len": f["min_len"], "max_len": f["max_len"], "min_q": f["min_q"],
    }
    if args.no_crosscheck:
        warnings.append("Chopper cross-check skipped (--no-crosscheck)")
    else:
        meta["chopper_in"], meta["chopper_out"] = run_crosscheck(
            slice_files, args.chopper, cfg, workers)
    return profiles, meta


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    cfg = load_config(args.config)
    qc_dir = args.output_dir / "qc_report"
    bins_path = qc_dir / "qc_bins.tsv.gz"
    warnings: list[str] = []

    if not args.output_dir.is_dir():
        sys.exit(f"ERROR: output dir not found: {args.output_dir}")
    log(f" QC config: {cfg['_path']}  version {config_label(cfg)}")

    sheet = load_samplesheet(args.samplesheet, cfg["controls"]["types"])
    if args.counts is not None:
        counts_path = args.counts
    else:
        counts_path, w = find_counts_file(args.output_dir)
        warnings += w
    counts = parse_read_counts(counts_path)
    log(f" Counts:    {counts_path} ({len(counts)} rows)")

    if args.report_only:
        profiles, meta = read_bins(bins_path)
        f = cfg["filter"]
        built = (meta.get("min_len"), meta.get("max_len"), meta.get("min_q"))
        if built != (f["min_len"], f["max_len"], f["min_q"]):
            warnings.append(
                f"qc_bins.tsv.gz was built with filter {built} (min_len, max_len, "
                f"min_q) but the config now says "
                f"{(f['min_len'], f['max_len'], f['min_q'])}; the categories "
                f"still reflect the old filter -- rerun without --report-only")
    else:
        qc_dir.mkdir(parents=True, exist_ok=True)
        profiles, meta = build_slice(args, cfg, qc_dir, warnings)
        write_bins(bins_path, profiles, meta)

    sheet_ids = {s.sample_id for s in sheet}
    extra = sorted(set(counts) - sheet_ids - set(CATCH_ALL))
    if extra:
        warnings.append(f"read_counts_summary.txt rows not in the sample "
                        f"sheet: {', '.join(extra)}")

    banner("QC STEP 6: Flagging and writing outputs")
    results, runs = evaluate(sheet, counts, profiles, cfg, warnings)

    header = tsv_header(cfg, meta, runs)
    write_sample_tsv(qc_dir / "per_sample_qc.tsv", results, header)
    client_dir = qc_dir / "per_client"
    if client_dir.is_dir():
        shutil.rmtree(client_dir)
    for client in sorted({r.sample.client for r in results}):
        write_sample_tsv(client_dir / f"{client}_qc.tsv",
                         [r for r in results if r.sample.client == client], header)
    write_plots(qc_dir, results, profiles, cfg, warnings)
    write_run_qc(qc_dir / "run_qc.txt", args, cfg, counts_path, results, runs,
                 profiles, meta, warnings)

    for ln in run_flag_lines(runs, cfg):
        log(f" {ln}")
    for flag in ("PASS", "CHECK", "FAIL", "CTRL"):
        log(f" {flag:<6} {sum(r.flag == flag for r in results):>5}")
    log(f" Report: {qc_dir}")
    log(f" QC report complete: {datetime.now():%Y-%m-%d %H:%M:%S}")
    log("========================================")
    return 0


if __name__ == "__main__":
    sys.exit(main())
