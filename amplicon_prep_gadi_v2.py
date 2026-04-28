from argparse import ArgumentParser as AP
from pathlib import Path
import os
from datetime import datetime
from shutil import copy2, rmtree
import gzip

"""
    amplicon_prep_gadi.py v1.2
    This script generates run scripts for assembling amplicon data on Gadi.
    It creates client-specific run scripts for the ONT amplicon pipeline.
"""

def generate_complete_run_script(top_dir_path, client_script_paths):
    """
    Make a top-level script that launches all client scripts via PBS
    top_dir_path = client_dir_path.parent
    """
    run_path = Path(top_dir_path) / 'run_amplicons.sh'
    with open(run_path,'wt') as fout:
        print('#!/bin/bash', file=fout)
        print('', file=fout)
        for csp in client_script_paths:
            print(f'qsub ./{csp.name}', file=fout)
    os.chmod(run_path, 0o755)
    print(f'Generated top-level script {run_path}')


def generate_client_run_script(client_sample_sheet_ref_path, client_sample_sheet_noref_path, client_info, client_sheet,
        client_path, pipeline_path, pipeline_version, minimap2_path, samtools_path, email, basecaller_cfg=None):
    """
    Inputs:
        client_sample_sheet_ref_path - path to client sample sheet with references
        client_sample_sheet_noref_path - path to client sample sheet without references
        client_info - client_info dictionary
        client_path - full path to client directory
        pipeline_path - singularity container path
        pipeline_version - e.g. v1.2.2
        minimap2_path - using module, so just name of executable
        samtools_path - using module, so just name of executable
        email - email address for PBS notifications
        basecaller_cfg - optional basecaller config override e.g. dna_r10.4.1_e8.2_400bps_sup@v5.0.0

    Runs:
    1) ONT Amplicon Pipeline wf-amplicon
    2) minimap2 of reads back to final assembly
    3) samtools index on the mapped BAM file

    """
    client_script_path = client_path.parent/f'run_{client_path.name}.qsub'
    client_name = client_path.name
    out_dn = client_name + "/output"
    nxf_base = '/g/data/vz35/amplicon_gadi'
    singularity_cache = f'{nxf_base}/singularity_cache'
    nextflow_path = 'nextflow'  # local binary, added to PATH below
    with open(client_script_path, 'wt') as fout:
        print('#!/bin/bash', file=fout)
        print('', file=fout)
        print(f'#PBS -N ampln_asm_{client_name}', file=fout)
        print('#PBS -P vz35', file=fout)
        print('#PBS -l mem=12GB,ncpus=4,walltime=4:00:00', file=fout)
        print('#PBS -q biodev', file=fout)
        print('#PBS -l storage=gdata/vz35', file=fout)
        print('#PBS -m abe', file=fout)
        print(f'#PBS -M {email}', file=fout)
        print('#PBS -l wd', file=fout)
        print('', file=fout)
        print('module load java/jdk-17.0.2', file=fout)
        print('module load singularity', file=fout)
        print('', file=fout)
        print(f'export PATH={nxf_base}:$PATH', file=fout)
        print(f'export NXF_HOME={nxf_base}', file=fout)
        print('export NXF_VER=23.10.1', file=fout)
        print(f'export SINGULARITY_CACHEDIR={singularity_cache}', file=fout)
        print('export NXF_SINGULARITY_CACHEDIR=$SINGULARITY_CACHEDIR', file=fout)
        print('export NXF_DISABLE_CHECK_LATEST=true', file=fout)
        print('', file=fout)
        if client_sample_sheet_ref_path:
            print('# ONT wf-amplicon pipeline with reference', file=fout)
            print(f'{nextflow_path} run {pipeline_path} -r {pipeline_version} \\', file=fout)
            print(f'  --fastq {client_name} \\', file=fout)
            print(f'  --out_dir {out_dn} \\', file=fout)
            print(f'  --sample_sheet ./{client_sample_sheet_ref_path.name} \\', file=fout)
            if basecaller_cfg:
                print(f'  --override_basecaller_cfg {basecaller_cfg} \\', file=fout)
            print(f'  -profile singularity \\', file=fout)
            print(f'  -offline', file=fout)
            print('', file=fout)
        if client_sample_sheet_noref_path:
            print('# ONT wf-amplicon pipeline without reference', file=fout)
            print(f'{nextflow_path} run {pipeline_path} -r {pipeline_version} \\', file=fout)
            print(f'  --fastq {client_name} \\', file=fout)
            print(f'  --out_dir {out_dn} \\', file=fout)
            print(f'  --sample_sheet ./{client_sample_sheet_noref_path.name} \\', file=fout)
            if basecaller_cfg:
                print(f'  --override_basecaller_cfg {basecaller_cfg} \\', file=fout)
            print(f'  -profile singularity \\', file=fout)
            print(f'  -offline', file=fout)
            print('', file=fout)

    os.chmod(client_script_path, 0o755)
    return client_script_path


def generate_sample_sheets(client_info: dict, client_path: Path, client_sheet: dict):
    """
    Generate two sample sheets for wf-amplicon, one with reference and one without.
    Comma separated sample sheet file covering all client samples.

    Inputs:
        client_info - dictionary of clients and samples
        client_path - full path to client directory
        client_sheet - user provided client/sample info (client,alias,barcode,reference)
    Returns:
        client_sample_sheet_noref_path, client_sample_sheet_ref_path

    note that 'barcode' is the name of the sample directory which contains fastq files for that sample
    type defaults to 'test_sample' but could also be 'positive_control','negative_control','no_template_control'
    headers for ref sheet: alias, barcode, type, reference
    headers for noref sheet: alias, barcode, type
    """
    client_sample_sheet_noref_path = None
    client_sample_sheet_ref_path = None
    samples_with_references = []
    samples_without_references = []
    for sample_name in client_info[client_path.name]:
        reference = str(client_info[client_path.name][sample_name].get('reference',''))
        if reference:
            samples_with_references.append(sample_name)
        else:
            samples_without_references.append(sample_name)

    if samples_with_references:
        client_sample_sheet_ref_path = client_path.parent / (str(client_path.name) + '_sample_sheet_ref.csv')
        with open(client_sample_sheet_ref_path, 'wt') as fout:
            print(','.join(['alias','barcode','type','reference']), file=fout)
            for sample_name in samples_with_references:
                alias = client_sheet[client_path.name][sample_name]['alias']
                barcode = sample_name
                reference = str(client_info[client_path.name][sample_name].get('reference',''))
                print(','.join([alias,barcode,'test_sample',reference]), file=fout)

    if samples_without_references:
        client_sample_sheet_noref_path = client_path.parent / (str(client_path.name) + '_sample_sheet_noref.csv')
        with open(client_sample_sheet_noref_path, 'wt') as fout:
            print(','.join(['alias','barcode','type']), file=fout)
            for sample_name in samples_without_references:
                alias = client_sheet[client_path.name][sample_name]['alias']
                barcode = sample_name
                print(','.join([alias,barcode,'test_sample']), file=fout)

    return client_sample_sheet_noref_path, client_sample_sheet_ref_path


def check_fastq_name(fn:str) -> bool:
    """
    Check that the FASTQ file name ends with the expected suffix. Ignore case.
    Returns True if name is good, otherwise False
    """
    suffix = ['.fq','.fq.gz','.fastq','.fastq.gz']
    for s in suffix:
        if str(fn).lower().endswith(s):
            return True
    return False


def check_fasta_name(fn:str) -> bool:
    """
    Check that the FASTA file name ends with the expected suffix. Ignore case
    Returns True if name is good, otherwise False
    """
    suffix = ['.fa','.fa.gz','.fasta','.fasta.gz']
    for s in suffix:
        if str(fn).lower().endswith(s):
            return True
    return False


def rename_fastq_to_bam(fp: str) -> Path|None:
    """
    Replace extension with .bam
    """
    suffix = ['.fq','.fq.gz','.fastq','.fastq.gz']
    for s in suffix:
        if str(fp).lower().endswith(s):
            return Path(str(fp).replace(s,'.bam'))


def parse_samplesheet(samplesheet: str) -> dict:
    """
    Reads an ONT wf-amplicon sample sheet: e.g.
    client,alias,barcode,reference
    A,amplicon1,barcode21,/path/to/reference.fa
    A,amplicon2,barcode22,ref.fasta
    B,amplicon3,barcode23,
    C,amplicon4,barcode24,

    Returns a dict [client]={barcode:{'alias':'','ref':'','fastqs':[]}}
    """
    if not Path(samplesheet).exists():
        print(f'Samplesheet {samplesheet} does not exist')
        exit(1)
    if not Path(samplesheet).is_file():
        print(f'Samplesheet {samplesheet} is not a file. Did you swap the parameters by mistake?')
        exit(1)

    client_info = {}
    client_barcode_aliases = {}  # per client set of (barcode,alias) to ensure uniqueness
    with open(samplesheet, 'rt') as f:
        for i,line in enumerate(f):
            cols = line.split(',')
            if i==0 and cols[0].lower().startswith('client'):
                continue  # header
            if len(cols) < 3:
                continue  # we must have at least 3 columns: client,alias,barcode
            client = cols[0].strip().replace(' ','_')
            alias = cols[1].strip().replace(' ','_')
            barcode = cols[2].strip()
            ref = ''
            if len(cols) > 3:  # optional reference in 4th column, ignore later columns
                ref = cols[3].strip()
            if client not in client_info:
                client_info[client] = {}
                client_barcode_aliases[client] = set()
            if (barcode,alias) not in client_barcode_aliases[client]:
                client_barcode_aliases[client].add((barcode,alias))
            else:
                print(f'barcode {barcode} and alias {alias} are not a unique combination in client {client}')
                exit(1)
            if barcode not in client_info[client]:
                client_info[client][barcode] = {'alias':alias,'ref':ref,'fastqs':[]}
            else:
                print(f'barcode {barcode} must be unique for client {client}')
                exit(1)
    return client_info


# now iterate through the PromethION directory tree looking for barcodes
def get_barcode_dirs(p: Path, all_barcodes: set, chosen_dirs: list) -> list:
    """
    iterate through directories finding all the barcode dirs we want
    """
    dirs = [x for x in p.iterdir() if x.is_dir()]
    for d in dirs:
        if d.name in ('fastq_pass', 'fastq_fail'):
            bc_dirs = [x for x in d.iterdir() if x.is_dir() if x.name in all_barcodes]
            if bc_dirs:
                chosen_dirs.extend(bc_dirs)
        else:
            get_barcode_dirs(d, all_barcodes, chosen_dirs)
    return chosen_dirs


def parse_input_dirs(prom_dir: str, client_sheet: dict) -> dict:
    """
    Scans a PromethION directory structure:
        Mla7_45_pool/
            -> 20241121_1136_3C_PAW74316_2656d858/
                -> fastq_pass/
                    -> barcode21/ (fastqs)
    Should be able to find everything listed in client_sheet
    returns a dict source_dirs[client] = {barcode:[list of barcode dirs from fastq_pass and/or fastq_fail]}
    """
    pdp = Path(prom_dir)
    if not pdp.exists():
        print(f"PromethION directory {pdp} does not exist")
        exit(1)
    if not pdp.is_dir():
        print(f"PromethION directory {pdp} is not a directory")
        exit(1)

    source_dirs = {}
    all_barcodes = set()
    for client in client_sheet:
        source_dirs[client] = {}
        for barcode in client_sheet[client]:
            source_dirs[client][barcode] = []  # list of src dirs (pass + fail)
            all_barcodes.add(barcode)

    barcode_dirs = get_barcode_dirs(pdp, all_barcodes, [])
    bcds = {}
    for bcd in barcode_dirs:
        bcds.setdefault(bcd.name, []).append(bcd)
    bcd_names = set(bcds.keys())
    if all_barcodes.difference(bcd_names):
        print(f"Barcodes not found {all_barcodes.difference(bcd_names)}")
        exit(2)
    if bcd_names.difference(all_barcodes):
        print(f"Extra barcodes found {bcd_names.difference(all_barcodes)}")
        exit(2)

    for client in source_dirs:
        for barcode in source_dirs[client]:
            source_dirs[client][barcode] = bcds[barcode]

    return source_dirs


def create_new_structure(plasmid_dir: Path, client_sheet: dict, source_dirs: dict, collapse=True, nodata=False, verbose=False) -> bool:
    """
    Create new plasmid directory tree

    args:
    plasmid_dir - Path to new plasmid directory
    client_sheet - dict of client and barcode info provided by the user
    source_dirs - dict of provided client and barcode directories
    collapse - bool, whether to collapse FASTQs into a single file
    nodata - bool, whether to skip copying data (for testing)
    verbose - bool, whether to display more information about the process

    returns:
    True if successful, otherwise False

    amplicon_run_20241217/
        -> clientA/
            -> barcode01/ (fastqs)
                -> reference/ (fasta)  optional
            -> barcode02/ (fastqs)
        -> clientB/
            -> barcode03/ (fastqs)
    By default the new barcode directories contain only the collapsed FASTQ file
    """
    #try:
        # make directories and copy files
    if True:
        if not plasmid_dir.exists():
            plasmid_dir.mkdir()
        for client in client_sheet:
            p = plasmid_dir / client
            if not p.exists():
                p.mkdir()
            for barcode in client_sheet[client]:
                bp = p/barcode
                if not bp.exists():
                    bp.mkdir()
                fps = [src_dir/f for src_dir in source_dirs[client][barcode]
                       for f in os.listdir(src_dir) if check_fastq_name(f)]
                if collapse:
                    collapse_fp = plasmid_dir/client/barcode/f'{barcode}.fq.gz'
                    if not nodata:
                        with gzip.open(collapse_fp,'wt') as fout:
                            for fp in fps:
                                if fp.name.lower().endswith('.gz'):
                                    if verbose:
                                        print(f'Collapsing {fp} to {collapse_fp}')
                                    with gzip.open(fp, 'rt') as f:
                                        for line in f:
                                            if line.strip():
                                                fout.write(line)
                                else:
                                    if verbose:
                                        print(f'Collapsing {fp} to {collapse_fp}')
                                    with open(fp, 'rt') as f:
                                        for line in f:
                                            if line.strip():
                                                fout.write(line)
                else:
                    for fp in fps:
                        if verbose:
                            print(f'Copying {fp} to {plasmid_dir/client/barcode}')
                        if not nodata:
                            dest_name = f'{fp.parent.parent.name}_{fp.name}'
                            copy2(fp, plasmid_dir/client/barcode/dest_name)

                ref = client_sheet[client][barcode]['ref']
                if ref:
                    ref_dp = bp/'reference'
                    if not ref_dp.exists():
                        ref_dp.mkdir()
                    if not nodata:
                        copy2(ref, ref_dp/Path(ref).name)
    # except Exception as exc:
    #     print(f'Failed to create new amplicon experiment directories {exc}')
    #     exit(3)
    return True


def main():
    """
    Generates run scripts for the ONT wf-amplicon pipeline on Gadi.

    Reads a (user provided) amplicon sample sheet of 3 or 4 columns (reference is optional): e.g.
        client,alias,barcode,reference
        A,amplicon1,barcode21,/path/to/reference.fa
        A,amplicon2,barcode22,ref.fasta
        B,amplicon3,barcode23,
        C,amplicon4,barcode24,

    Scans a PromethION directory structure:
        Mla7_45_pool/
            -> 20241121_1136_3C_PAW74316_2656d858/
                -> fastq_pass/
                    -> barcode21/ (fastqs)

    And creates the required directory structure and all scripts, etc, for the ONT
    wf-amplicon pipeline:
    amplicon_run_20241217/
        -> clientA/
            -> barcode01/ (fastqs)
                -> reference/ (fasta)  optional
            -> barcode02/ (fastqs)
        -> clientB/
            -> barcode03/ (fastqs)

    Requires one top-level directory for all the client amplicons you want to run in one go.
    Inside, you'd have one directory for each client. In each client directory you'd have
    an experiment directory for each separate amplicon (the barcode name, or sample name).
    Within that are all the FASTQ sequences files for that amplicon, and an optional
    directory called "reference" - if you have a reference.

    """
    dt = datetime.today().strftime('%Y%m%d')
    amplicon_dn = f'amplicon_run_{dt}'

    parser = AP()
    parser.add_argument('prom_dir', help='Path to input PromethION sequencing')
    parser.add_argument('-s','--samplesheet', required=True, help='Path to 3 or 4 column samplesheet to set up experiment')
    parser.add_argument('-p', '--amplicon_dir', default=amplicon_dn, help='Path to output folder containing all client amplicon data')
    parser.add_argument('-v', '--verbose', action='store_true', help='Display more information about the prep process')
    parser.add_argument('-o','--overwrite', action='store_true', help='Overwrite existing amplicon directory')
    parser.add_argument('--pipeline_path', default='epi2me-labs/wf-amplicon', help='Path to ONT wf-amplicon pipeline')
    parser.add_argument('--pipeline_version', default='v1.2.2', help='wf-amplicon pipeline version')
    parser.add_argument('--basecaller_cfg', default=None, help='Override basecaller config e.g. dna_r10.4.1_e8.2_400bps_sup@v5.0.0')
    parser.add_argument('--no_collapse', action='store_true', help='Disable collapsing FASTQs to a single file for each barcode')
    parser.add_argument('--minimap2', default='minimap2', help='Path to minimap2 executable (using module, so just name of executable)')
    parser.add_argument('--samtools', default='samtools', help='Path to samtools executable (using module, so just name of executable)')
    parser.add_argument('--nodata', action='store_true', help='Run the script without creating any files, for testing purposes')
    parser.add_argument('-e','--email', required=True, help='Email address for PBS notifications')

    args = parser.parse_args()

    prom_dir = Path(args.prom_dir)
    if not prom_dir.exists():
        print(f'PromethION sequencing directory {prom_dir} does not exist')
        exit(1)

    amplicon_dir = Path(args.amplicon_dir)
    if amplicon_dir.exists() and not args.nodata:
        if not args.overwrite:
            print(f'Amplicon run directory {amplicon_dir} already exists. '+\
                    f'Please delete it, choose to overwrite it, or name a different output directory')
            exit(1)
        else:
            rmtree(amplicon_dir)
    try:
        amplicon_dir.mkdir()
    except FileExistsError:
        if args.nodata:
            pass
        else:
            print(f'Amplicon run directory {amplicon_dir} already exists. Please delete it or name a different output directory')
            exit(1)

    # client sheet is the user input about each client and sample
    client_sheet = parse_samplesheet(args.samplesheet)

    # create a dictionary of provided barcode directories for each client
    source_dirs = parse_input_dirs(args.prom_dir, client_sheet)

    collapse_fastqs = True
    if args.no_collapse:
        collapse_fastqs = False
    success = create_new_structure(amplicon_dir, client_sheet, source_dirs, collapse=collapse_fastqs, nodata=args.nodata, verbose=args.verbose)

    if success:
        print(f'Successfully created amplicon directory {amplicon_dir}')

    # each sample/alias has its own set of records, alias is the barcode name and is the directory that holds the FASTQ files
    client_info = {}
    client_script_paths = []
    for client in client_sheet:
        cdir = amplicon_dir/client
        client_info[cdir.name] = {}
        sample_dirs = [d for d in cdir.glob('*') if d.is_dir() and str(d.name).startswith('barcode')]
        if not sample_dirs:
            print(f'Skipping client {cdir}, no sample directories found')
            continue
        if args.verbose:
            print(f'{cdir} contains samples: {sample_dirs}')

        # now create a sample_sheet.csv for each client so we run all their jobs together
        for sd in sample_dirs:
            client_info[cdir.name][sd.name] = {}

            if not args.nodata:
                seq_fns = [fp for fp in sd.glob('*') if fp.is_file() and check_fastq_name(fp.name)]
                if not seq_fns:
                    print(f'No FASTQ (.fq/.fastq/.fq.gz/.fastq.gz files found for client {cdir.name} sample {sd.name}')
                    exit(1)
                client_info[cdir.name][sd.name]['fastq_files'] = seq_fns
            else:
                client_info[cdir.name][sd.name]['fastq_files'] = []

            ref_dir = sd.joinpath('reference')  # optional

            if not ref_dir.exists():
                if args.verbose:
                    print(f'Reference directory {ref_dir} not found. Client {cdir.name} sample {sd.name}')
            else:
                if not ref_dir.is_dir():
                    print(f'Reference directory {ref_dir} is not a directory! Client {cdir.name} sample {sd.name}')
                    exit(1)
                ref_fp = [f for f in ref_dir.glob('*') if f.is_file() and check_fasta_name(f.name)]
                if len(ref_fp) != 1:
                    print(f'Reference files {ref_fp} found. There should be exactly one reference file')
                    exit(1)
                client_info[cdir.name][sd.name]['reference'] = Path(cdir.name)/sd.name/'reference'/ref_fp[0].name

        # generate client sample sheets without, and with, references
        client_sample_sheet_noref_path, client_sample_sheet_ref_path = generate_sample_sheets(client_info, cdir, client_sheet)
        if client_sample_sheet_noref_path:
            print(f'Created sample sheet without reference sequences {client_sample_sheet_noref_path} for client {cdir.name}')
        if client_sample_sheet_ref_path:
            print(f'Created sample sheet with reference sequences {client_sample_sheet_ref_path} for client {cdir.name}')

        client_run_script_path = generate_client_run_script(client_sample_sheet_ref_path,
                client_sample_sheet_noref_path, client_info, client_sheet, cdir,
                args.pipeline_path, args.pipeline_version, args.minimap2, args.samtools, args.email,
                args.basecaller_cfg)
        print(f'Created script {client_run_script_path} for client {cdir.name}')
        client_script_paths.append(client_run_script_path)

    # generate an overall run script that launches everything else
    generate_complete_run_script(amplicon_dir, client_script_paths)


if __name__ == '__main__':
    main()
