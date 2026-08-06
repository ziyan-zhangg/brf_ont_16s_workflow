# ONT 16S Workflow - BRF @ Gadi

Pipeline for filtering, demultiplexing, orientation normalized and read-counting Oxford Nanopore 16S amplicon data on the NCI Gadi HPC.

---

## Overview

The workflow runs in a single PBS job (`run_script.qsub`) that calls a pre-step script and then the core pipeline:

| Stage | Script / Tool | What it does |
|-------|--------------|--------------|
| Pre-step | `generate_primer_setup.py` | Converts the sample sheet CSV into a per-sample primer setup file |
| Step 1 | `filter_chopper_demux_minibar.py` -> **Chopper** | Quality- and length-filters every raw `fastq.gz` — **`pre` mode only** |
| Step 2 | -> **Minibar** | Demultiplexes each input file into per-sample FASTQs, then merges across all `fastq_pass` files |
| Step 3 | -> `demux_qc.py` | Per-sample length/quality profile of the demultiplexed reads, before filtering |
| Step 4 | -> **Chopper** | Quality- and length-filters each merged per-sample FASTQ — **`post` mode only** |
| Step 5 | -> **Cutadapt** (two-pass) | Normalises read orientation and trims 5'/3' primers |
| Step 6 | -> organise + summarise | Groups reads by client, writes per-client and run-level read-count summaries |
| Step 7 | -> log collection | Moves all `*.log` / `*.txt` files into `run_log_<date>/` |

After running the pipeline on the control sample and confirming that demultiplexing counts and per-sample yields match expectations, the same `run_script.qsub` can be submitted unchanged for production samples.

### Filter stage: when Chopper runs

`filter_stage` in `run_script.qsub` decides whether Chopper runs **before** demultiplexing (step 1) or **after** it, per sample (step 4):

| Mode | Chopper runs | Minibar sees | Notes |
|------|--------------|--------------|-------|
| `pre` *(default)* | Step 1, on whole `fastq.gz` files | Filtered reads | Unchanged legacy behaviour |
| `post` | Step 4, on merged per-sample FASTQs | Raw `fastq_pass` reads | No pre-filter at all — `fastq_pass` is already MinKNOW min-qscore filtered |

Both modes share one step numbering, so step *N* means the same thing in either job log and the two can be read side by side. Each mode skips exactly one step, and the skip is announced in the log rather than silently absent:

```
pre                                   post
STEP 1: Chopper filtering             STEP 1: Chopper filtering SKIPPED
STEP 2: Minibar demultiplexing        STEP 2: Minibar demultiplexing
STEP 3: Per-sample demux QC           STEP 3: Per-sample demux QC
STEP 4: Per-sample Chopper SKIPPED    STEP 4: Chopper filtering (per sample)
STEP 5: Cutadapt two-pass             STEP 5: Cutadapt two-pass
STEP 6: Organising reads by client    STEP 6: Organising reads by client
STEP 7: Collecting logs               STEP 7: Collecting logs
```

`pre` is byte-for-byte identical to the pipeline as it stood before the two-stage option was added, so the two orderings can be run on the same input and compared directly.

**Why `post` matters.** Filtering before demultiplexing throws reads away before anyone can see what they looked like, so a low per-sample yield is ambiguous: too few reads, or the right reads outside the assumed 1–2 kb window? Running Minibar first makes the demux loss and the filter loss separately visible (steps 3 and 4), at the cost of Minibar processing 2–4× more reads.

---

## Repository contents

```
run_script.qsub                  # PBS job wrapper -- edit Variables section here
generate_primer_setup.py         # Pre-step: sample sheet -> primer setup file
filter_chopper_demux_minibar.py  # Core pipeline: Steps 1-7
demux_qc.py                      # Step 3: per-sample demux QC (imported by core script)
cutadapt_2pass.py                # Step 5: two-pass cutadapt (imported by core script)
```

`demux_qc.py` and `cutadapt_2pass.py` are imported as siblings of the core script, so all four `.py` files must sit in the same directory on Gadi.

---

## Prerequisites

All tools are expected to be present under `/g/data/vz35/ONT_16s_workflow/`:

| Tool | Default path |
|------|-------------|
| Python 3.12 | via `module load python3/3.12.1` |
| Chopper | `/g/data/vz35/ONT_16s_workflow/tools/chopper/chopper-linux-musl` |
| Minibar | `/g/data/vz35/ONT_16s_workflow/tools/minibar/minibar.py` |
| Cutadapt | `/g/data/vz35/ONT_16s_workflow/tools/cutadapt-env/bin/cutadapt` |
| `generate_primer_setup.py` | `/g/data/vz35/ONT_16s_workflow/tools/brf_ont_16s_workflow/` |
| `filter_chopper_demux_minibar.py` | `/g/data/vz35/ONT_16s_workflow/tools/brf_ont_16s_workflow/` |
| `demux_qc.py` | `/g/data/vz35/ONT_16s_workflow/tools/brf_ont_16s_workflow/` |
| `cutadapt_2pass.py` | `/g/data/vz35/ONT_16s_workflow/tools/brf_ont_16s_workflow/` |
| Twist 384 barcode reference | `/g/data/vz35/ONT_16s_workflow/tools/Twist_16S_384_barcode.txt` |

---

## Input files

### 1. Sample sheet (CSV)

Required columns:

| Column | Description |
|--------|-------------|
| `Client` | Client name -- used to organise output into subdirectories |
| `Sample_ID` | Sample identifier (alphanumeric, `-`, `_`; other characters are sanitised) |
| `Barcode` | Twist 384 barcode ID matching a row in the barcode reference |
| `Comment` | *(optional)* Shown in the per-client summary for samples below the low-read threshold |

### 2. Raw reads

PromethION `fastq_pass` directory located at:

```
/g/data/vz35/PromethION_data/sequencer_uploads/<run_name>/
```

The script locates the `fastq_pass` folder automatically from the `run_name` variable.

Only one `fastq_pass` folder is allowed per run.

---

## Configuration

Edit the **Variables** section near the top of `run_script.qsub`:

```bash
run_name=ONT_16S_20260422          # Subdirectory under sequencer_uploads/
samplesheet=/g/data/vz35/ONT_16s_workflow/sample_sheet/16s_samplesheet.csv
output_dir=/g/data/vz35/ONT_16s_workflow/minibar_output/ONT_16S_TBC_<date>
filter_stage=pre                   # pre | post -- see "Filter stage" above
```

Everything in the `DONT-CHANGE` section below resolves automatically from these four variables.

In `post` mode the wrapper additionally points Minibar's per-file working directories at `$PBS_JOBFS` (node-local disk) and passes `--minibar-work-dir` to the core script. Only the merged per-sample FASTQs are written to `output_dir`. Nothing outside `/g/data/vz35` and `$PBS_JOBFS` is touched, so `-l storage=gdata/vz35` remains correct in both modes.

---

## Usage

### 1. Prepare the sample sheet

Fill in `16s_samplesheet.csv` with `Client`, `Sample_ID`, and `Barcode` columns.
Include a control sample with a known expected yield to validate the run.

Optional `Comment` can be put based on the amplification tendency.

### 2. Update the Variables section

Open `run_script.qsub` and set `run_name`, `samplesheet`, `output_dir`, and `filter_stage`.

### 3. Submit the job

```bash
qsub run_script.qsub
```

PBS resources requested: 2 CPUs, 10 GB RAM, 200 GB jobfs, 20 h walltime, `storage=gdata/vz35`.

Notes on resources:

- **jobfs is 200 GB** (raised from 10 GB). In `post` mode Minibar's per-file working directories hold unfiltered reads and are several times larger than in `pre` mode. On queues that ration jobfs per core, 200 GB may need more than `ncpus=2` to be schedulable — raise `ncpus` before lowering `jobfs`.
- **Walltime is unchanged at 20 h**, but `post` mode feeds Minibar 2–4× more reads and Minibar is the dominant cost, so the margin is much thinner. Check the step 2 timestamps in the job log on the first `post` run before assuming 20 h holds.

### 4. Validate with the control run

After the job completes, open `read_counts_summary.txt` and the control sample's `summary.txt` inside `integrated_demultiplexing/<Client>/`.
Confirm that the control sample's read count and percentage match the expected values for that barcode.
Once validated, the same script can be resubmitted for any other sample set by updating the Variables section.

If a sample came out low, `read_counts_summary.txt` now says *where* the reads went — `pct_retained` for filter loss, `pct_pass1_kept` for non-specific amplification — and `demux_qc/demux_qc_summary.tsv` gives the length and quality profile behind those numbers.

### 5. (Optional) Generate the primer setup file standalone

The primer file is generated automatically inside the job. To create it independently:

```bash
python3 generate_primer_setup.py <samplesheet.csv> [-o OUTPUT_DIR] [-b BARCODE_FILE] [--date YYYYMMDD]
```

Writes `16S_primer_setup_<date>.txt` (tab-separated) with columns:
`SampleID`, `FwIndex`, `FwPrimer`, `RvIndex`, `RvPrimer`.

The default barcode reference is the Twist 384 file. Use `-b` to supply an alternative barcode file for runs with external barcodes.

---

## Pipeline steps in detail

### Pre-step: generate_primer_setup.py

Reads the sample sheet and looks up each barcode in the Twist 384 reference to produce a tab-separated primer setup file.
The file is written to the `sample_sheet/` directory and consumed by all downstream steps.

---

### Step 1: Chopper -- quality and length filtering (`pre` mode only)

Skipped entirely when `filter_stage=post`, where raw `fastq_pass/*.fastq.gz` go straight to Minibar and filtering happens at step 4 instead.

Each `*.fastq.gz` file in `fastq_pass/` is piped through Chopper, and reads that pass all three filters are retained:

| Filter | Parameter | Default |
|--------|-----------|---------|
| Quality score | `-q` / `--min-quality` | >= 15 |
| Minimum length | `--minlength` | 1,000 bp |
| Maximum length | `--maxlength` | 2,000 bp |

All three are now command-line options on the core script rather than hardcoded (`--min-length` / `--max-length` are accepted as aliases). The equivalent command run per file:

```
zcat <file>.fastq.gz | chopper -q 15 --minlength 1000 --maxlength 2000 | gzip > <file>_filtered.fastq.gz
```

Filtered files are written to `chopper_filtered/` and deleted automatically after demultiplexing.

---

### Step 2: Minibar -- demultiplexing

Minibar demultiplexes each input file into per-sample FASTQs using the primer setup file.
Each file is processed in its own subdirectory, and per-sample results are merged across all `fastq_pass` files at the end.

Flags are identical in both filter stages; only the input differs (Chopper-filtered files in `pre`, raw `fastq_pass` files in `post`).

| Flag | Value | Meaning |
|------|-------|---------|
| `-e 1` | 1 | Allowed mismatches in barcode |
| `-E 5` | 5 | Allowed mismatches in primer |
| `-l 200` | 200 | Search window at each read end (bp) |
| `-M 2` | 2 | Require barcode match on both ends |
| `-F` | -- | Write each sample to its own file |

Output files: `sample_<SampleID>.fastq`, `sample_unk.fastq`, `sample_Multiple_Matches.fastq`.
Per-file Minibar subdirectories are removed after merging to save storage. In `post` mode these subdirectories live on `$PBS_JOBFS` rather than in `output_dir`.

---

### Step 3: Per-sample demux QC

`demux_qc.py` makes a single streaming pass over each merged `sample_<SampleID>.fastq` **before** any filtering, and reports:

| Metric | Detail |
|--------|--------|
| Raw demuxed read count | Reads Minibar assigned to that barcode |
| Length histogram | 50 bp bins, plus exact N50, median, and modal length |
| Mean read Q | Per-read mean error probability converted to Phred, then averaged over reads |
| % within window | Share of reads inside `[minlength, maxlength]` |

`sample_unk` and `sample_Multiple_Matches` are profiled too — an inflated unassigned bin is the first sign of a demux problem, and its length profile says whether the cause is chemistry or barcode matching.

Outputs, written to `demux_qc/`:

| File | Contents |
|------|----------|
| `<SampleID>.tsv` | One row per 50 bp bin; scalar stats in a `#`-prefixed header block (read with `pandas.read_csv(..., sep='\t', comment='#')`) |
| `demux_qc_summary.tsv` | One row per sample, all scalar stats, real samples first |

The step runs in **both** filter stages, but means different things in each: in `post` mode it profiles genuinely raw demultiplexed reads; in `pre` mode those reads have already been through Chopper. Every QC output records which stage produced it, so the two are never confused. Skip with `--no-demux-qc`.

---

### Step 4: Chopper -- per-sample filtering (`post` mode only)

Skipped when `filter_stage=pre`, where Chopper already ran at step 1.

Applies the same `-q` / `--minlength` / `--maxlength` filter as step 1, but to each merged `sample_<SampleID>.fastq` individually, filtering in place. Pre-filter counts are carried over from step 3 so the demux loss and the filter loss stay separable in the run summary.

`sample_unk` and `sample_Multiple_Matches` are deliberately **left unfiltered** — their raw length profile is the whole reason they are diagnostic. They are therefore also excluded from the TOTAL `pct_retained` denominator in `read_counts_summary.txt`.

---

### Step 5: Cutadapt -- orientation normalisation and primer trimming

Run on each `sample_<SampleID>.fastq` after merging. Reads from `sample_unk` and `sample_Multiple_Matches` are left untouched.

**Pass 1 -- orient and trim 5' primer (strict)**

```
cutadapt -g FWD_PRIMER --revcomp --rename={header} --discard-untrimmed -e 0.2
```

- Searches for the forward primer at the 5' end `-g FWD`
- If the primer is found on the minus strand, the read is reverse-complemented `--revcomp` so all output reads face the same direction
- `--rename={header}` keep the name unchanged
- Reads where the forward primer cannot be found are discarded `--discard-untrimmed` — these are likely noise or off-target sequences

**Pass 2 -- trim 3' primer (tolerant)**

```
cutadapt -a REV_PRIMER_RC -e 0.2
```

- Trims the reverse-complemented reverse primer from the 3' end.
- Reads where the 3' construct is not found are **kept** (truncated reads are retained).

Output: `<SampleID>.fastq.gz` replaces `sample_<SampleID>.fastq`.
Per-sample logs are written to `cutadapt_logs/` and a run-level `cutadapt_summary.txt` is produced.

---

### Step 6: Organise by client and summarise

Demultiplexed files are moved into `integrated_demultiplexing/<Client>/` based on the sample sheet.
Two summary files are written:

| File | Contents |
|------|----------|
| `integrated_demultiplexing/<Client>/summary.txt` | Per-sample read counts for that client; flags samples below 15,000 reads |
| `read_counts_summary.txt` | Run-level totals plus the per-stage read-count breakdown below |

`read_counts_summary.txt` now tracks the read count through every stage, so the two-stage loss is separable:

| Column | Meaning |
|--------|---------|
| `Reads` | Final read count in the delivered FASTQ |
| `% of Input` | Share of total raw (`post`) or total filtered (`pre`) input |
| `raw_demuxed` | Reads Minibar assigned, before filtering |
| `post_chopper` | Reads surviving the Chopper length/Q filter |
| `pct_retained` | `post_chopper / raw_demuxed` |
| `post_cutadapt_pass1` | Reads with a findable forward primer |
| `pct_pass1_kept` | `post_cutadapt_pass1 / cutadapt input` |

**`pct_pass1_kept` is the non-specificity signal.** A barcode with a healthy `raw_demuxed` but a low `pct_pass1_kept` means Minibar found the right barcode on reads that are not the expected product — correct barcode, wrong amplicon. A `post_chopper` count alone hides this.

In `pre` mode no raw demuxed count exists (reads were already filtered when Minibar saw them), so `raw_demuxed` and `pct_retained` are left **blank** rather than being filled with a post-filter number.

### Step 7: log collection

All `*.log` and `*.txt` files (including `cutadapt_logs/`) are moved into `run_log_<date>/` to keep the output root clean. `demux_qc/` is a directory of `.tsv` files and stays where it is.

---

## Output structure

```
minibar_output/ONT_16S_TBC_<date>/
├── integrated_demultiplexing/
│   ├── <ClientA>/
│   │   ├── <SampleID>.fastq.gz
│   │   └── summary.txt
│   ├── <ClientB>/
│   │   └── ...
│   ├── sample_Multiple_Matches.fastq
│   └── sample_unk.fastq
├── demux_qc/
│   ├── <SampleID>.tsv
│   ├── unk.tsv
│   ├── Multiple_Matches.tsv
│   └── demux_qc_summary.tsv
└── run_log_<date>/
    ├── cutadapt_summary.txt
    ├── read_counts_summary.txt
    └── cutadapt_logs/
        ├── <SampleID>.pass1_orient.log
        └── <SampleID>.pass2_trim3.log
```

- `integrated_demultiplexing/` -- reads merged across all `fastq_pass` files and split by client
- `sample_unk.fastq` -- reads that did not match any barcode
- `sample_Multiple_Matches.fastq` -- reads that matched more than one barcode
- `summary.txt` -- per-client read counts; low-count samples (< 15,000 reads) show the `Comment` column from the sample sheet
- `demux_qc/` -- per-sample length histograms and quality stats taken before filtering
- `read_counts_summary.txt` -- overall run summary; ends up in `run_log_<date>/` after step 7

---

## Tool parameter reference

### Chopper

Configurable on the core script; the values below are the defaults.

| Parameter | Default | Effect |
|-----------|---------|--------|
| `--min-quality` / `-q` | 15 | Minimum mean quality score |
| `--minlength` | 1,000 bp | Discard reads shorter than this |
| `--maxlength` | 2,000 bp | Discard reads longer than this |

### Minibar

| Parameter | Value | Effect |
|-----------|-------|--------|
| `-e 1` | 1 | Barcode mismatch tolerance |
| `-E 5` | 5 | Primer mismatch tolerance |
| `-l 200` | 200 | End-window search length (bp) |
| `-M 2` | 2 | Match barcodes on both ends |
| `-F` | -- | Write separate file per sample |

### Cutadapt

| Parameter | Value | Effect |
|-----------|-------|--------|
| `-e 0.2` | 0.2 | Error rate for adapter matching (both passes) |
| `-j 4` | 4 | Threads per sample |
| `--revcomp` | -- | Pass 1: reverse-complement reads where primer found on minus strand |
| `--discard-untrimmed` | -- | Pass 1: drop reads with no forward primer match |
| `-g FWD` | forward primer | Pass 1: 5'-anchored adapter |
| `-a REV_RC` | RC of reverse primer | Pass 2: 3' adapter |

---

## Core script options

`run_script.qsub` covers the routine case. Run the core script directly for anything else:

```bash
python3 filter_chopper_demux_minibar.py \
    --raw-input-dir /path/to/fastq_pass \
    --output-dir    /path/to/run_output \
    --primer-file   /path/to/16S_primer_setup_YYYYMMDD.txt \
    --samplesheet   /path/to/16s_samplesheet.csv \
    [--filter-stage pre|post] [--diagnostic] \
    [--min-quality 15] [--minlength 1000] [--maxlength 2000] \
    [--minibar-work-dir $PBS_JOBFS] [--no-demux-qc]
```

| Option | Effect |
|--------|--------|
| `--filter-stage {pre,post}` | When Chopper runs. Default `pre`. |
| `--diagnostic` | Shorthand for `--filter-stage post --minlength 100 --maxlength 5000`. Explicitly given flags still win. |
| `--min-quality`, `--minlength`, `--maxlength` | Chopper filter parameters (`--min-length` / `--max-length` accepted as aliases). |
| `--minibar-work-dir` | Where Minibar's per-file working directories are created. Default: the output dir. |
| `--no-demux-qc` | Skip step 3. |
| `--chopper`, `--minibar`, `--cutadapt` | Override tool paths. |
| `--cutadapt-threads`, `--cutadapt-error-rate` | Cutadapt tuning (defaults 4 and 0.2). |

### Diagnostic mode

For a first look at an unfamiliar run, `--diagnostic` forces `post` mode with a 100–5,000 bp window so step 3 observes the **real** length distribution instead of only what survives an assumed 1–2 kb window. Read `demux_qc/demux_qc_summary.tsv`, pick a window that matches the observed modal length and N50, then do the production run with those values.

### Running demux QC standalone

`demux_qc.py` also works on an existing run directory without rerunning the pipeline:

```bash
python3 demux_qc.py \
    --integrated-dir /path/to/run_output/integrated_demultiplexing \
    [--qc-dir /path/to/demux_qc] \
    [--minlength 1000] [--maxlength 2000] [--filter-stage pre|post]
```

It only reads the FASTQs — it never rewrites reads — so it is safe to run against a completed run. Note that after step 6 the per-sample files have been renamed to `<SampleID>.fastq.gz` and moved into client subdirectories, so a standalone run against the top level will only find the `sample_unk` / `sample_Multiple_Matches` bins.
