#!/bin/bash
# =============================================================================
# amplicon_cleanup.sh
# Cleans up and reorganises an amplicon run folder on Gadi.
#
# Usage: bash amplicon_cleanup.sh [--dry-run] <run_folder_name>
# Example: bash amplicon_cleanup.sh amplicon_run_20260416
#          bash amplicon_cleanup.sh --dry-run amplicon_run_20260416
#
# Working directory: /g/data/vz35/amplicon_gadi
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
BASE_DIR="/g/data/vz35/amplicon_gadi"
DRY_RUN=false

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Prints a message prefixed with [DRY-RUN] when in dry-run mode
log_action() {
    if [[ "$DRY_RUN" == true ]]; then
        echo "  [DRY-RUN] $*"
    else
        echo "  $*"
    fi
}

# Wrapper: runs a command for real, or just prints it in dry-run mode
run() {
    if [[ "$DRY_RUN" == true ]]; then
        echo "  [DRY-RUN] would run: $*"
    else
        "$@"
    fi
}

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
if [[ $# -eq 0 || $# -gt 2 ]]; then
    echo "Usage: bash amplicon_cleanup.sh [--dry-run] <run_folder_name>"
    echo "Example: bash amplicon_cleanup.sh amplicon_run_20260416"
    echo "         bash amplicon_cleanup.sh --dry-run amplicon_run_20260416"
    exit 1
fi

if [[ "$1" == "--dry-run" ]]; then
    DRY_RUN=true
    if [[ $# -ne 2 ]]; then
        echo "ERROR: --dry-run requires a run folder name."
        echo "Usage: bash amplicon_cleanup.sh --dry-run <run_folder_name>"
        exit 1
    fi
    RUN_FOLDER="$2"
else
    RUN_FOLDER="$1"
fi

RUN_DIR="${BASE_DIR}/${RUN_FOLDER}"

# Extract the date index from the folder name (e.g. "20260416" from "amplicon_run_20260416")
DATE_INDEX="${RUN_FOLDER##*_}"

# Validate the run directory exists
if [[ ! -d "$RUN_DIR" ]]; then
    echo "ERROR: Run directory not found: ${RUN_DIR}"
    exit 1
fi

# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
echo "============================================================"
echo "  Amplicon cleanup script"
echo "  Base dir : ${BASE_DIR}"
echo "  Run dir  : ${RUN_DIR}"
echo "  Date idx : ${DATE_INDEX}"
if [[ "$DRY_RUN" == true ]]; then
    echo "  Mode     : DRY-RUN — no changes will be made"
else
    echo "  Mode     : LIVE — changes will be applied"
fi
echo "============================================================"
echo ""

# ---------------------------------------------------------------------------
# Step 1 — Remove amplicon_setup scheduler logs from the base directory
# ---------------------------------------------------------------------------
echo "[1/4] Removing amplicon_setup scheduler logs from base directory..."

shopt -s nullglob
setup_logs=("${BASE_DIR}/amplicon_setup.e"* "${BASE_DIR}/amplicon_setup.o"*)
if [[ ${#setup_logs[@]} -gt 0 ]]; then
    for f in "${setup_logs[@]}"; do
        log_action "Removing: $f"
        run rm -f "$f"
    done
else
    echo "  No amplicon_setup.e* / amplicon_setup.o* files found."
fi
shopt -u nullglob

echo ""

# ---------------------------------------------------------------------------
# Step 2 — Remove run_* scripts and ampln_asm_cilent_name scheduler logs
#           from inside the run folder
# ---------------------------------------------------------------------------
echo "[2/4] Removing run_* scripts and ampln_asm_cilent_name logs from ${RUN_DIR}..."

shopt -s nullglob

# run_* scripts
run_scripts=("${RUN_DIR}/run_"*)
if [[ ${#run_scripts[@]} -gt 0 ]]; then
    for f in "${run_scripts[@]}"; do
        log_action "Removing: $f"
        run rm -f "$f"
    done
else
    echo "  No run_* scripts found."
fi

# scheduler logs for ampln_asm_cilent_name
asm_logs=("${RUN_DIR}/ampln_asm_cilent_name.e"* "${RUN_DIR}/ampln_asm_cilent_name.o"*)
if [[ ${#asm_logs[@]} -gt 0 ]]; then
    for f in "${asm_logs[@]}"; do
        log_action "Removing: $f"
        run rm -f "$f"
    done
else
    echo "  No ampln_asm_cilent_name.e* / ampln_asm_cilent_name.o* files found."
fi

shopt -u nullglob
echo ""

# ---------------------------------------------------------------------------
# Step 3 — Inside each client folder: clean up output* subfolders
#           Remove: execution/ folder, params.json, versions.txt
# ---------------------------------------------------------------------------
echo "[3/4] Cleaning up output* folders inside each client directory..."

for client_dir in "${RUN_DIR}"/*/; do
    [[ -d "$client_dir" ]] || continue
    client_name="$(basename "$client_dir")"

    # Skip folders that have already been renamed (start with the date index)
    if [[ "$client_name" == "${DATE_INDEX}_"* ]]; then
        continue
    fi

    echo "  Client: ${client_name}"

    shopt -s nullglob
    output_dirs=("${client_dir}"output*/)
    shopt -u nullglob

    if [[ ${#output_dirs[@]} -eq 0 ]]; then
        echo "    No output* folders found."
        continue
    fi

    for output_dir in "${output_dirs[@]}"; do
        [[ -d "$output_dir" ]] || continue
        echo "    Output folder: $(basename "$output_dir")"

        # Remove execution/ subdirectory
        if [[ -d "${output_dir}execution" ]]; then
            log_action "Removing: execution/"
            run rm -rf "${output_dir}execution"
        fi

        # Remove params.json
        if [[ -f "${output_dir}params.json" ]]; then
            log_action "Removing: params.json"
            run rm -f "${output_dir}params.json"
        fi

        # Remove versions.txt
        if [[ -f "${output_dir}versions.txt" ]]; then
            log_action "Removing: versions.txt"
            run rm -f "${output_dir}versions.txt"
        fi
    done
done

echo ""

# ---------------------------------------------------------------------------
# Step 4 — Rename client folders inside the run directory
#           Pattern: <ClientName> -> <DATE_INDEX>_<ClientName>
# ---------------------------------------------------------------------------
echo "[4/4] Renaming client folders to ${DATE_INDEX}_<ClientName>..."

for client_dir in "${RUN_DIR}"/*/; do
    [[ -d "$client_dir" ]] || continue
    client_name="$(basename "$client_dir")"

    # Skip if already prefixed with the date index
    if [[ "$client_name" == "${DATE_INDEX}_"* ]]; then
        echo "  Already renamed, skipping: ${client_name}"
        continue
    fi

    new_name="${DATE_INDEX}_${client_name}"
    new_path="${RUN_DIR}/${new_name}"

    if [[ -d "$new_path" ]]; then
        echo "  WARNING: Target already exists, skipping: ${new_name}"
        continue
    fi

    log_action "Renaming: ${client_name} -> ${new_name}"
    run mv "${client_dir}" "${new_path}"
done

echo ""
echo "============================================================"
if [[ "$DRY_RUN" == true ]]; then
    echo "  Dry-run complete. No files were modified."
    echo "  Re-run without --dry-run to apply changes."
else
    echo "  Cleanup complete."
fi
echo "============================================================"
