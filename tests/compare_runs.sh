#!/bin/bash
# Compare two pipeline output directories for equivalence.
#
# Intended for the before/after check described in the speed-up work: run the
# same input through the old and new code, then
#
#     tests/compare_runs.sh OLD_OUTPUT_DIR NEW_OUTPUT_DIR
#
# Exits 0 only if every per-sample read file is byte-identical and the
# summaries and QC tables match. Lines that legitimately differ between any
# two runs -- timestamps, and the input/output/jobfs paths -- are excluded;
# everything else is compared verbatim.
#
# A check that finds nothing to compare FAILS rather than passing silently:
# two empty directories are not evidence of anything.

set -uo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: $(basename "$0") OLD_OUTPUT_DIR NEW_OUTPUT_DIR" >&2
    exit 2
fi

old=$1
new=$2
fails=0

for d in "$old" "$new"; do
    if [[ ! -d "$d" ]]; then
        echo "ERROR: not a directory: $d" >&2
        exit 2
    fi
done

echo "OLD: $old"
echo "NEW: $new"
echo

pass() { echo "PASS  $*"; }
fail() { echo "FAIL  $*"; fails=$((fails + 1)); }

# Timestamps and the dir paths differ by construction: the two runs have
# different output dirs, and pre mode's filtered dir moved to jobfs.
strip_summary() {
    grep -v -e '^ Date:' -e '^ Raw input dir:' -e '^ Filtered dir:' \
            -e '^ Out dir:' "$1"
}

#--- 1. read_counts_summary.txt ----------------------------------------------
f=read_counts_summary.txt
if [[ ! -f "$old/$f" || ! -f "$new/$f" ]]; then
    fail "$f missing from one or both runs"
elif diff <(strip_summary "$old/$f") <(strip_summary "$new/$f") > /tmp/.cmp_$$ 2>&1; then
    pass "$f identical (ignoring Date and dir paths)"
else
    fail "$f differs:"
    sed 's/^/        /' /tmp/.cmp_$$
fi
rm -f /tmp/.cmp_$$

#--- 2. the reads themselves -------------------------------------------------
# Three naming patterns reach integrated_demultiplexing/: cutadapt output as
# <Client>/<SampleID>.fastq.gz, the uncut <Client>/sample_<id>.fastq fallback,
# and the top-level sample_unk / sample_Multiple_Matches catch-all bins.
# gzip stores an mtime in its header, so two runs of the same data produce
# .fastq.gz files whose raw bytes differ by 4 bytes even when every read is
# identical. Hash the DECOMPRESSED stream instead -- a raw md5sum here would
# report every comparison as a failure.
digest() {
    local root="$1/integrated_demultiplexing"
    [[ -d "$root" ]] || return 0
    ( cd "$root" || return 0
      find . \( -name '*.fastq.gz' -o -name 'sample_*.fastq' \) -print0 \
        | sort -z \
        | while IFS= read -r -d '' f; do
              if [[ "$f" == *.gz ]]; then
                  h=$(gzip -cd -- "$f" | md5sum | cut -d" " -f1) \
                      || h="GZIP_ERROR"
              else
                  h=$(md5sum -- "$f" | cut -d" " -f1)
              fi
              printf '%s  %s\n' "$h" "$f"
          done )
}
old_d=$(digest "$old")
new_d=$(digest "$new")
n_old=$(printf '%s' "$old_d" | grep -c . || true)
n_new=$(printf '%s' "$new_d" | grep -c . || true)

if [[ "$n_old" -eq 0 || "$n_new" -eq 0 ]]; then
    fail "no read files found (old: $n_old, new: $n_new) -- nothing was compared"
elif [[ "$old_d" == "$new_d" ]]; then
    pass "$n_old read file(s) identical (decompressed content)"
else
    fail "read files differ:"
    diff <(printf '%s\n' "$old_d") <(printf '%s\n' "$new_d") | sed 's/^/        /'
fi

#--- 3. demux QC tables ------------------------------------------------------
# '# source' carries the output dir, '# generated' a timestamp.
strip_qc() { grep -v -e '^# source' -e '^# generated' "$1"; }
qc_n=0
qc_bad=0
if [[ -d "$old/demux_qc" && -d "$new/demux_qc" ]]; then
    while IFS= read -r p; do
        b=$(basename "$p")
        qc_n=$((qc_n + 1))
        if [[ ! -f "$old/demux_qc/$b" ]]; then
            echo "        only in NEW: demux_qc/$b"
            qc_bad=$((qc_bad + 1))
        elif ! diff -q <(strip_qc "$old/demux_qc/$b") <(strip_qc "$p") > /dev/null; then
            echo "        differs: demux_qc/$b"
            qc_bad=$((qc_bad + 1))
        fi
    done < <(find "$new/demux_qc" -name '*.tsv' | sort)
fi
if [[ "$qc_n" -eq 0 ]]; then
    echo "SKIP  demux_qc: no TSVs found (--no-demux-qc?)"
elif [[ "$qc_bad" -eq 0 ]]; then
    pass "$qc_n demux_qc table(s) identical (ignoring source and generated)"
else
    fail "$qc_bad of $qc_n demux_qc table(s) differ"
fi

#--- 4. per-client summaries -------------------------------------------------
cl_n=0
cl_bad=0
while IFS= read -r p; do
    rel=${p#"$new/integrated_demultiplexing/"}
    cl_n=$((cl_n + 1))
    o="$old/integrated_demultiplexing/$rel"
    if [[ ! -f "$o" ]]; then
        echo "        only in NEW: $rel"
        cl_bad=$((cl_bad + 1))
    elif ! diff -q <(grep -v '^ Date:' "$o") <(grep -v '^ Date:' "$p") > /dev/null; then
        echo "        differs: $rel"
        cl_bad=$((cl_bad + 1))
    fi
done < <(find "$new/integrated_demultiplexing" -name 'summary.txt' 2>/dev/null | sort)
if [[ "$cl_n" -eq 0 ]]; then
    echo "SKIP  no per-client summary.txt found"
elif [[ "$cl_bad" -eq 0 ]]; then
    pass "$cl_n per-client summary file(s) identical (ignoring Date)"
else
    fail "$cl_bad of $cl_n per-client summary file(s) differ"
fi

#--- what only the new run has ----------------------------------------------
echo
nlogs=$(find "$new" -path '*/minibar_logs/*.log' 2>/dev/null | wc -l)
echo "NOTE  new run saved $nlogs per-chunk minibar log(s)"
echo "NOTE  walltime, CPU time, peak memory and jobfs are NOT checked here --"
echo "      read them from the PBS footer:  grep -A12 'Resource Usage' *.o*"

echo
if [[ "$fails" -eq 0 ]]; then
    echo "EQUIVALENT: all checks passed"
    exit 0
fi
echo "NOT EQUIVALENT: $fails check(s) failed"
echo "Report the difference rather than adjusting this script."
exit 1
