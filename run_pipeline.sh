#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh - LungScape Python Pipeline Orchestrator
#
# Runs the full Python preprocessing pipeline in sequence:
#   Step 0  - Dataset preparation       (prepare_dataset.py)
#   Step 1  - TotalSegmentator          (run_TotalSegmentator.py)
#   Step 2a - nnU-Net Model201 setup    (setup_model201.py)   [skipped if already done]
#   Step 2b - nnU-Net Model201 airways  (run_model201.py)     [skipped with --use-ts-airways]
#   Step 3a - nnU-Net Model191 setup    (setup_model191.py)   [skipped if already done]
#   Step 3b - nnU-Net Model191 pathology (run_model191.py)    [skipped with --skip-model191]
#
# Usage:
#   ./run_pipeline.sh --input <input_dir> --output <output_dir> \
#                     --model201-dir <path> --model191-dir <path> [options]
#
# Required:
#   --input,       -i  <path>   Input directory (DICOM folders or NIfTI files)
#   --output,      -o  <path>   Working directory (DATA/ and nnUNet_INPUT/ created here)
#   --model201-dir     <path>   Directory for Model201 installation/reuse
#                               (not required with --use-ts-airways)
#   --model191-dir     <path>   Directory for Model191 installation/reuse
#                               (not required with --skip-model191)
#
# Optional:
#   --skip-model191             Skip Model191 high attenuation segmentation
#   --nifti                     Input is NIfTI instead of DICOM
#   --skip-resampling           Skip volume resampling in Step 0
#   --use-ts-airways            Use TotalSegmentator airways (skips Model201 steps)
#   --no-cleanup                Keep nnU-Net folders after inference (for debugging/reuse)
#   --help,        -h           Show this help message
#
# Notes:
#   - The script uses whichever Python is active in the current environment.
#     Activate your virtualenv/conda env before running.
#   - Model setup steps (2a, 3a) are skipped automatically if the model is
#     already installed (nnUNet_results/ folder exists and is non-empty).
#   - Full log is saved to <output_dir>/pipeline.log.
#   - After completion, set DATA_DIRECTORY in LungScape.m and run the MATLAB pipeline.
#
# Example:
#   conda activate lungscape
#   ./run_pipeline.sh \
#     --input  /data/dicom \
#     --output /data/lungscape_output \
#     --model201-dir /models/model201 \
#     --model191-dir /models/model191
#
#   # Without Model191 (patients without significant high attenuation pathology):
#   ./run_pipeline.sh \
#     --input  /data/dicom \
#     --output /data/lungscape_output \
#     --model201-dir /models/model201 \
#     --model191-dir /models/model191 \
#     --skip-model191
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Defaults ──────────────────────────────────────────────────────────────────
INPUT_DIR=""
OUTPUT_DIR=""
MODEL201_DIR=""
MODEL191_DIR=""
SKIP_MODEL191=false
NIFTI_FLAG=""
SKIP_RESAMPLING_FLAG=""
USE_TS_AIRWAYS=false
NO_CLEANUP_FLAG=""

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
    sed -n '/^# Usage:/,/^# ====/p' "$0" | sed 's/^# \?//'
}

# ── Argument parsing ──────────────────────────────────────────────────────────
if [[ $# -eq 0 ]]; then
    usage
    exit 1
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input|-i)          INPUT_DIR="$2";       shift 2 ;;
        --output|-o)         OUTPUT_DIR="$2";      shift 2 ;;
        --model201-dir)      MODEL201_DIR="$2";    shift 2 ;;
        --model191-dir)      MODEL191_DIR="$2";    shift 2 ;;
        --skip-model191)     SKIP_MODEL191=true;   shift ;;
        --nifti)             NIFTI_FLAG="--nifti"; shift ;;
        --skip-resampling)   SKIP_RESAMPLING_FLAG="--skip-resampling"; shift ;;
        --use-ts-airways)    USE_TS_AIRWAYS=true;  shift ;;
        --no-cleanup)        NO_CLEANUP_FLAG="--no-cleanup"; shift ;;
        --help|-h)           usage; exit 0 ;;
        *) echo "[ERROR] Unknown argument: $1"; echo; usage; exit 1 ;;
    esac
done

# ── Validation ────────────────────────────────────────────────────────────────
errors=()

[[ -z "$INPUT_DIR" ]]  && errors+=("--input is required.")
[[ -z "$OUTPUT_DIR" ]] && errors+=("--output is required.")

if ! $USE_TS_AIRWAYS && [[ -z "$MODEL201_DIR" ]]; then
    errors+=("--model201-dir is required (or use --use-ts-airways to skip Model201).")
fi

if ! $SKIP_MODEL191 && [[ -z "$MODEL191_DIR" ]]; then
    errors+=("--model191-dir is required (or use --skip-model191 to exclude high attenuation segmentation).")
fi

if [[ ${#errors[@]} -gt 0 ]]; then
    for e in "${errors[@]}"; do echo "[ERROR] $e"; done
    echo
    usage
    exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "[ERROR] Input directory not found: $INPUT_DIR"
    exit 1
fi

# ── Python detection ──────────────────────────────────────────────────────────
if command -v python3 &>/dev/null; then
    PYTHON="python3"
elif command -v python &>/dev/null; then
    PYTHON="python"
else
    echo "[ERROR] Python not found in PATH. Activate your environment and retry."
    exit 1
fi

# ── Logging setup ─────────────────────────────────────────────────────────────
mkdir -p "$OUTPUT_DIR"
LOG_FILE="$OUTPUT_DIR/pipeline.log"
exec > >(tee -a "$LOG_FILE") 2>&1

# ── Helpers ───────────────────────────────────────────────────────────────────
step_header() {
    local step="$1" desc="$2"
    echo ""
    echo "================================================================="
    echo "  STEP ${step}: ${desc}"
    echo "  $(date '+%Y-%m-%d %H:%M:%S')"
    echo "================================================================="
}

step_skip() {
    local step="$1" reason="$2"
    echo ""
    echo "  [SKIP] Step ${step}: ${reason}"
}

model_is_installed() {
    local model_dir="$1"
    local results_dir="${model_dir}/nnUNet_results"
    [[ -d "$results_dir" ]] && [[ -n "$(ls -A "$results_dir" 2>/dev/null)" ]]
}

# ── Banner ────────────────────────────────────────────────────────────────────
echo "================================================================="
echo "  LungScape - Python Pipeline"
echo "  $(date '+%Y-%m-%d %H:%M:%S')"
echo "================================================================="
echo ""
echo "  LungSCAPE was developed by Alberto Arrigoni"
echo "  Istituto di Ricerche Farmacologiche Mario Negri IRCCS."
echo ""
echo "  If you use LungScape in your work, find it useful for your"
echo "  research, or build upon it for further development, please cite"
echo "  the corresponding publication:"
echo "    Arrigoni, A. et al. Radiol med (2025)."
echo "    https://doi.org/10.1007/s11547-025-02166-w"
echo "  and the GitHub repository."
echo "================================================================="
echo "  Input:          $INPUT_DIR"
echo "  Output:         $OUTPUT_DIR"
[[ -n "$MODEL201_DIR" ]] && echo "  Model201 dir:   $MODEL201_DIR"
$SKIP_MODEL191             && echo "  Model191:       SKIPPED (--skip-model191)" \
                           || echo "  Model191 dir:   $MODEL191_DIR"
$USE_TS_AIRWAYS            && echo "  Airways:        TotalSegmentator (--use-ts-airways)"
[[ -n "$NIFTI_FLAG" ]]     && echo "  Input format:   NIfTI"
echo "  Log file:       $LOG_FILE"
echo "================================================================="

# ── Step 0: Dataset Preparation ───────────────────────────────────────────────
step_header "0" "Dataset Preparation (prepare_dataset.py)"

"$PYTHON" "$SCRIPT_DIR/prepare_dataset.py" \
    "$INPUT_DIR" "$OUTPUT_DIR" \
    $NIFTI_FLAG $SKIP_RESAMPLING_FLAG

# ── Step 1: TotalSegmentator ──────────────────────────────────────────────────
step_header "1" "TotalSegmentator (run_TotalSegmentator.py)"

TS_CMD=("$PYTHON" "$SCRIPT_DIR/run_TotalSegmentator.py" "--data-dir" "$OUTPUT_DIR/DATA")
$USE_TS_AIRWAYS && TS_CMD+=("--use-ts-airways")

"${TS_CMD[@]}"

# ── Step 2a: Model201 Setup ───────────────────────────────────────────────────
if $USE_TS_AIRWAYS; then
    step_skip "2a+2b" "nnU-Net Model201 skipped (--use-ts-airways active)"
else
    step_header "2a" "nnU-Net Model201 Setup (setup_model201.py)"

    if model_is_installed "$MODEL201_DIR"; then
        echo "  Model201 already installed in $MODEL201_DIR — skipping download."
    else
        "$PYTHON" "$SCRIPT_DIR/setup_model201.py" "$MODEL201_DIR"
    fi

    # ── Step 2b: Model201 Inference ───────────────────────────────────────────
    step_header "2b" "nnU-Net Model201 Airways Segmentation (run_model201.py)"

    MODEL201_CMD=("$PYTHON" "$SCRIPT_DIR/run_model201.py" "$MODEL201_DIR"
                  "--project-dir" "$OUTPUT_DIR")
    [[ -n "$NO_CLEANUP_FLAG" ]] && MODEL201_CMD+=("--no-cleanup")

    "${MODEL201_CMD[@]}"
fi

# ── Step 3a+3b: Model191 ─────────────────────────────────────────────────────
if $SKIP_MODEL191; then
    step_skip "3a+3b" "nnU-Net Model191 skipped (--skip-model191)"
else
    step_header "3a" "nnU-Net Model191 Setup (setup_model191.py)"

    if model_is_installed "$MODEL191_DIR"; then
        echo "  Model191 already installed in $MODEL191_DIR — skipping download."
    else
        "$PYTHON" "$SCRIPT_DIR/setup_model191.py" "$MODEL191_DIR"
    fi

    step_header "3b" "nnU-Net Model191 High Attenuation Segmentation (run_model191.py)"

    MODEL191_CMD=("$PYTHON" "$SCRIPT_DIR/run_model191.py" "$MODEL191_DIR"
                  "--project-dir" "$OUTPUT_DIR")
    [[ -n "$NO_CLEANUP_FLAG" ]] && MODEL191_CMD+=("--no-cleanup")

    "${MODEL191_CMD[@]}"
fi

# ── Data validation ───────────────────────────────────────────────────────────
step_header "✓" "Pre-MATLAB Data Validation (check_data_ready.py)"
"$PYTHON" "$SCRIPT_DIR/check_data_ready.py" "$OUTPUT_DIR/DATA" || {
    echo ""
    echo "[WARNING] Some patients are not ready for MATLAB processing."
    echo "          Check the report above, resolve the issues, and re-run"
    echo "          check_data_ready.py before launching LungScape.m."
    echo ""
}

# ── Done ──────────────────────────────────────────────────────────────────────
echo ""
echo "================================================================="
echo "  Python pipeline completed successfully."
echo "  $(date '+%Y-%m-%d %H:%M:%S')"
echo ""
echo "  DATA directory:  $OUTPUT_DIR/DATA"
echo "  Full log:        $LOG_FILE"
echo ""
echo "  Next step: open LungScape.m, set"
echo "    DATA_DIRECTORY = '$OUTPUT_DIR/DATA'"
echo "  and run the MATLAB pipeline."
echo "================================================================="
