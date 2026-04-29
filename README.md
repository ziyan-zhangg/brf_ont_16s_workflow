# ONT 16S Workflow — BRF @ Gadi

Pipeline for filtering, demultiplexing, and read-counting Oxford Nanopore 16S amplicon data on the NCI Gadi HPC.

---

## Overview

The workflow runs in three stages inside a single PBS job:

| Stage | Tool | What it does |
|-------|------|--------------|
| Pre-step | `generate_primer_setup.py` | Converts a sample sheet CSV into a Minibar primer file |
| Step 1 | Chopper | Quality- and length-filters raw reads (Q ≥ 15, 1–2 kb) |
| Step 2 | Minibar | Demultiplexes filtered reads per sample, merges across flow-cell files |
| Step 3 | Sort and count | Organises output by client, generates read-count summaries |

---

## Repository contents

```
Filter_Chopper_Demux_Minibar.qsub   # Main PBS job script
generate_primer_setup.py            # Primer setup file generator
changelog.txt                       # Version history
```

---

## Prerequisites

All tools are expected to be present under `/g/data/vz35/ONT_16s_workflow/`:

| Tool | Default path |
|------|-------------|
| Python 3.12 | via `module load python3/3.12.1` |
| Chopper | `/g/data/vz35/zpfeng/tools/chopper/chopper-linux-musl` |
| Minibar | `/g/data/vz35/ONT_16s_workflow/tools/minibar/minibar.py` |
| `generate_primer_setup.py` | `/g/data/vz35/ONT_16s_workflow/tools/generate_primer_setup.py` |
| Twist 384 barcode reference | `/g/data/vz35/ONT_16s_workflow/tools/Twist_16S_384_barcode.txt` |

---

## Input files

### 1. Sample sheet (CSV)

Required columns:

| Column | Description |
|--------|-------------|
| `Client` | Client name (used to organise output into subdirectories) |
| `Sample_ID` | Sample identifier (alphanumeric, `-`, `_`; other characters are sanitised) |
| `Barcode` | Twist 384 barcode ID matching a row in the barcode reference |

Sample sheet corresponds barcode info

---

### 2. Raw reads

PromethION `fastq_pass` directory located at:
```
/g/data/vz35/PromethION_data/sequencer_uploads/<run_name>/
```

The script finds the `fastq_pass` folder automatically.

---

## Configuration

Edit the **Variables** section near the top of `Filter_Chopper_Demux_Minibar.qsub`:

```bash
run_name=ONT_16S_20260422          # Subdirectory under sequencer_uploads/
samplesheet=/g/data/vz35/ONT_16s_workflow/sample_sheet/16s_samplesheet.csv
output_dir=/g/data/vz35/ONT_16s_workflow/minibar_output/ONT_16S_TBC_<date>
```

All other paths in the `DONT-CHANGE` section resolve automatically.

---

## Usage

### 1. Generate the primer setup file (standalone)

```bash
python3 generate_primer_setup.py <samplesheet.csv> [-o OUTPUT_DIR] [-b BARCODE_FILE] [--date YYYYMMDD]
```

Writes `16S_primer_setup_<date>.txt` (tab-separated) with columns:
`SampleID`, `FwIndex`, `FwPrimer`, `RvIndex`, `RvPrimer`

### 2. Submit the full pipeline

```bash
qsub Filter_Chopper_Demux_Minibar.qsub
```

PBS resources requested: 2 CPUs, 10 GB RAM, 10 GB jobfs, 20 h walltime, queue `biodev`.

---

## Output structure

```
minibar_output/ONT_16S_TBC_<date>/
├── integrated_demultiplexing/
│   ├── <ClientA>/
│   │   ├── sample_<SampleID>.fastq
│   │   └── summary.txt
│   ├── <ClientB>/
│   │   └── ...
│   └── sample_Multiple_Matches.fastq
│   └── sample_unk.fastq
└── read_counts_summary.txt
```

- `integrated_demultiplexing/` — demultiplexed reads merged across all fastq.gz files under fastq_pass
- `summary.txt` — per-client read counts
- `read_counts_summary.txt` — overall summary: total filtered input, total demultiplexed, successfully demultiplexed, and per-sample percentages
- Chopper filtered reads and per-file Minibar subdirectories are deleted automatically after merging to save storage

---

## Minibar parameters

| Flag | Value | Meaning |
|------|-------|---------|
| `-e 1` | 1 | Allowed errors in barcode |
| `-E 5` | 5 | Allowed errors in primer |
| `-l 200` | 200 | Search window length (bp) |
| `-M 2` | 2 | Match barcodes on both ends |
| `-T` | — | Trim barcode and primer from output |
| `-F` | — | Write each sample to its own file |

---

## Chopper parameters

Reads passing all three filters are kept:

- Quality score ≥ 15
- Length ≥ 1,000 bp
- Length ≤ 2,000 bp

---
