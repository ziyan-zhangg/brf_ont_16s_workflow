# ONT Amplicon Assembly Pipeline on Gadi

This toolkit prepares and submits ONT amplicon assembly jobs on the NCI Gadi HPC system using the [epi2me-labs/wf-amplicon](https://github.com/epi2me-labs/wf-amplicon) Nextflow pipeline.

---

## Files

| File | Description |
|---|---|
| `amplicon_prep_gadi.py` | Main Python script that generates all PBS job scripts and directory structure |
| `amplicon_setup.qsub` | PBS launcher script that runs `amplicon_prep_gadi.py` on Gadi |

---

## Prerequisites

### One-time setup (login node only)

These steps only need to be done once. Gadi compute nodes have no internet access, so all tools and containers must be pre-cached.

**1. Install Nextflow binary**
```bash
module load java/jdk-17.0.2
cd /g/data/vz35/amplicon_gadi
curl https://get.nextflow.io | bash
```

**2. Set environment and pull the pipeline**
```bash
module load java/jdk-17.0.2
export PATH=/g/data/vz35/amplicon_gadi:$PATH
export NXF_HOME=/g/data/vz35/amplicon_gadi
export NXF_VER=23.10.1
export NXF_DISABLE_CHECK_LATEST=true
export SINGULARITY_CACHEDIR=/g/data/vz35/amplicon_gadi/singularity_cache
export NXF_SINGULARITY_CACHEDIR=$SINGULARITY_CACHEDIR

nextflow pull epi2me-labs/wf-amplicon -r v1.2.2
```

**3. Pre-pull Singularity containers**

Run the pipeline once on the login node (it will fail due to memory limits but will cache the containers):
```bash
module load singularity
nextflow run epi2me-labs/wf-amplicon -r v1.2.2 \
  --fastq <your_fastq_dir> \
  --out_dir <output_dir> \
  --sample_sheet <sample_sheet.csv> \
  -profile singularity
```

Any containers not pulled automatically can be pulled manually:
```bash
singularity pull \
  --name ontresearch-medaka-sha<hash>.img \
  docker://ontresearch/medaka:sha<hash>

# Move to cache directory
mv ontresearch-medaka-sha<hash>.img /g/data/vz35/amplicon_gadi/singularity_cache/
```

---

## Input: Sample Sheet

A CSV file with 3 or 4 columns (reference is optional):

```
client,alias,barcode,reference
Sarah_Kaines,SK48,barcode44,
Sarah_Kaines,SK49,barcode45,
AnotherClient,AC01,barcode01,/path/to/reference.fa
```

| Column | Description |
|---|---|
| `client` | Client/project name — sets the output directory name |
| `alias` | Sample alias used in output filenames (e.g. SK48) |
| `barcode` | Barcode directory name in the PromethION data |
| `reference` | Optional path to a reference FASTA file |

Samples with a reference and samples without are handled separately — two sample sheets are generated automatically.

---

## Usage

### Step 1: Edit the launcher script

Open `amplicon_setup.qsub` and set the variables at the top:

```bash
Email="your.email@anu.edu.au"
PromData="ONT_PlasmidSeq_20260218/plasmidpool/20260218_1629_3F_PBE83525_e4822fd3"
AmpDir="amplicon_run_20260218"
SampleSheet="ONT_PlasmidSeq_20260218/plasmid_samplesheet_20260218.csv"
```

### Step 2: Submit the launcher

```bash
cd /path/to/your/working/directory
qsub amplicon_setup.qsub
```

This runs `amplicon_prep_gadi.py` which will:
1. Scan the PromethION directory for the barcodes listed in the sample sheet
2. Create a structured output directory under `AmpDir/`
3. Collapse all FASTQ files per barcode into a single `.fq.gz` file
4. Generate per-client sample sheets (with and without reference)
5. Generate per-client PBS job scripts (`run_<client>.qsub`)
6. Generate a top-level launcher script (`run_amplicons.sh`)

### Step 3: Submit the assembly jobs

```bash
cd <AmpDir>
./run_amplicons.sh
```

This submits a PBS job for each client.

---

## Output Directory Structure

```
AmpDir/
├── run_amplicons.sh               # Top-level launcher — runs all client jobs
├── ClientA_sample_sheet_noref.csv
├── ClientA_sample_sheet_ref.csv
├── run_ClientA.qsub               # PBS job script for ClientA
├── ClientA/
│   ├── barcode01/
│   │   └── barcode01.fq.gz        # Collapsed FASTQ
│   ├── barcode02/
│   │   ├── barcode02.fq.gz
│   │   └── reference/
│   │       └── reference.fa
│   └── output/                    # wf-amplicon results
│       ├── SK48.final.fasta
│       └── SK49.final.fasta
└── ClientB/
    └── ...
```

---

## Advanced Options

These can be passed to `amplicon_prep_gadi.py` directly or added to `amplicon_setup.qsub`:

| Option | Default | Description |
|---|---|---|
| `--pipeline_version` | `v1.2.2` | wf-amplicon version to use |
| `--basecaller_cfg` | None | Override basecaller e.g. `dna_r10.4.1_e8.2_400bps_sup@v5.0.0` |
| `--no_collapse` | False | Disable collapsing FASTQs into a single file per barcode |
| `--overwrite` | False | Overwrite existing output directory |
| `--nodata` | False | Dry run — generate scripts only, don't copy files |
| `-v` / `--verbose` | False | Print more detail during setup |

Example with basecaller override:
```bash
python3 amplicon_prep_gadi.py \
  -s samplesheet.csv \
  -p amplicon_run_20260218 \
  -e your.email@anu.edu.au \
  --basecaller_cfg dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
  ONT_PlasmidSeq_20260218/plasmidpool/20260218_1629_3F_PBE83525_e4822fd3
```

---

## Troubleshooting

**`curl: (7) Failed to connect to www.nextflow.io`**
The system `nextflow` module tries to download from the internet on compute nodes. Make sure you are using the local binary at `/g/data/vz35/amplicon_gadi/nextflow` and not loading the `nextflow` module.

**`Cannot find epi2me-labs/wf-amplicon`**
The pipeline cache was not found. Re-run `nextflow pull epi2me-labs/wf-amplicon -r v1.2.2` on the login node with `NXF_HOME=/g/data/vz35/amplicon_gadi` set.

**`Failed to pull singularity image ... network is unreachable`**
A Singularity container is missing from the cache. Pull it manually on the login node (see Prerequisites step 3).

**Job name in error files hard to read**
PBS error files are named `ampln_asm_<ClientName>.e<jobid>` — the client name is included automatically.