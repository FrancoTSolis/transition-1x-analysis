#!/bin/bash
# Batch runner for LUCJ computation via ffsim.
#
# Runs compute_lucj.py in parallel across all completed CCSD job directories.
# Each job reads/writes only within its own subdirectory (no race conditions).
#
# IMPORTANT: ffsim is installed via pip but the project root contains a local
# ffsim/ source tree that would shadow it. This script runs from /tmp to avoid
# the issue.
#
# Usage:
#   bash run_lucj_batch.sh [options]
#
# Options:
#   --jobs-dir DIR       Directory containing job subdirs (default: jobs)
#   --parallel-jobs N    Number of simultaneous jobs (default: 6)
#   --n-reps N           LUCJ circuit layer repetitions (default: 2)
#   --maxiter N          Max optimization iterations (default: 50)
#   --threads-per-job N  CPU threads per job (default: 6)
#   --log FILE           Log file (default: run_lucj_batch.log)
#   --max-jobs N         Limit jobs to process (0 = all, default: 0)
#   --dry-run            Preview without executing
#
# Background execution:
#   nohup bash run_lucj_batch.sh --parallel-jobs 6 > lucj_batch.log 2>&1 &
#   disown

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
COMPUTE_SCRIPT="$SCRIPT_DIR/compute_lucj.py"

JOBS_DIR="$SCRIPT_DIR/jobs"
PARALLEL_JOBS=6
N_REPS=2
MAXITER=50
THREADS_PER_JOB=6
LOG_FILE="$SCRIPT_DIR/run_lucj_batch.log"
MAX_JOBS=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --jobs-dir)        JOBS_DIR="$2";        shift 2;;
        --parallel-jobs)   PARALLEL_JOBS="$2";   shift 2;;
        --n-reps)          N_REPS="$2";          shift 2;;
        --maxiter)         MAXITER="$2";         shift 2;;
        --threads-per-job) THREADS_PER_JOB="$2"; shift 2;;
        --log)             LOG_FILE="$2";        shift 2;;
        --max-jobs)        MAX_JOBS="$2";        shift 2;;
        --dry-run)         DRY_RUN=1;            shift;;
        *) echo "Unknown option: $1"; exit 1;;
    esac
done

TOTAL_CORES=$(nproc)
USED_CORES=$((THREADS_PER_JOB * PARALLEL_JOBS))

echo "=== LUCJ Batch Runner (ffsim) ==="
echo "Jobs directory:     $JOBS_DIR"
echo "Parallel jobs:      $PARALLEL_JOBS"
echo "Threads per job:    $THREADS_PER_JOB"
echo "Total cores used:   $USED_CORES / $TOTAL_CORES"
echo "n_reps:             $N_REPS"
echo "maxiter:            $MAXITER"
echo "Log file:           $LOG_FILE"
echo ""

if [ "$USED_CORES" -gt "$TOTAL_CORES" ]; then
    echo "WARNING: Requesting $USED_CORES cores but only $TOTAL_CORES available!"
    echo ""
fi

PENDING_JOBS=()
DONE_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0

for JOB in "$JOBS_DIR"/*/; do
    [ -d "$JOB" ] || continue
    if [ -f "$JOB/.lucj_done" ]; then
        DONE_COUNT=$((DONE_COUNT + 1))
    elif [ -f "$JOB/.lucj_failed" ]; then
        FAIL_COUNT=$((FAIL_COUNT + 1))
        PENDING_JOBS+=("$JOB")
    elif [ ! -f "$JOB/ccsd_t1.dat" ] || [ ! -f "$JOB/ccsd_t2.dat" ]; then
        SKIP_COUNT=$((SKIP_COUNT + 1))
    else
        PENDING_JOBS+=("$JOB")
    fi
done

TOTAL_PENDING=${#PENDING_JOBS[@]}
echo "Total job directories: $((DONE_COUNT + FAIL_COUNT + SKIP_COUNT + TOTAL_PENDING))"
echo "Already completed:     $DONE_COUNT"
echo "Previously failed:     $FAIL_COUNT (will retry)"
echo "No amplitude files:    $SKIP_COUNT (skipped)"
echo "Pending:               $TOTAL_PENDING"
echo ""

if [ "$MAX_JOBS" -gt 0 ] && [ "$MAX_JOBS" -lt "$TOTAL_PENDING" ]; then
    PENDING_JOBS=("${PENDING_JOBS[@]:0:$MAX_JOBS}")
    echo "Limiting to first $MAX_JOBS jobs"
fi

if [ "$DRY_RUN" -eq 1 ]; then
    echo "--- DRY RUN ---"
    for JOB in "${PENDING_JOBS[@]}"; do
        echo "  $(basename "$JOB")"
    done
    exit 0
fi

if [ "$TOTAL_PENDING" -eq 0 ]; then
    echo "Nothing to do!"
    exit 0
fi

echo "Starting batch at $(date)" | tee -a "$LOG_FILE"
echo "Config: $PARALLEL_JOBS parallel x $THREADS_PER_JOB threads, n_reps=$N_REPS, maxiter=$MAXITER" >> "$LOG_FILE"

run_lucj_job() {
    local job_dir="$1"
    local job_name
    job_name=$(basename "$job_dir")

    rm -f "$job_dir/.lucj_failed"

    export OMP_NUM_THREADS="$THREADS_PER_JOB"
    export OPENBLAS_NUM_THREADS="$THREADS_PER_JOB"
    export MKL_NUM_THREADS="$THREADS_PER_JOB"
    export JAX_NUM_THREADS="$THREADS_PER_JOB"

    cd /tmp
    python3 "$COMPUTE_SCRIPT" "$job_dir" \
        --n-reps "$N_REPS" \
        --maxiter "$MAXITER" \
        2>&1
}
export -f run_lucj_job
export COMPUTE_SCRIPT THREADS_PER_JOB N_REPS MAXITER

printf '%s\n' "${PENDING_JOBS[@]}" | \
    xargs -P "$PARALLEL_JOBS" -I {} bash -c 'run_lucj_job "$@"' _ {} \
    2>&1 | tee -a "$LOG_FILE"

FINAL_DONE=$(find "$JOBS_DIR" -name ".lucj_done" | wc -l)
FINAL_FAIL=$(find "$JOBS_DIR" -name ".lucj_failed" | wc -l)
TOTAL_DIRS=$(find "$JOBS_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l)

echo ""
echo "=== Batch Complete ==="
echo "Completed: $FINAL_DONE / $TOTAL_DIRS"
echo "Failed:    $FINAL_FAIL / $TOTAL_DIRS"
echo "Finished at $(date)"
