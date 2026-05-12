#!/bin/bash
# Run Q-Chem CCSD amplitude jobs in parallel batches.
#
# Usage:
#   bash run_batch.sh [options]
#
# Options:
#   --jobs-dir DIR       Directory containing job subdirs (default: jobs)
#   --threads-per-job N  OpenMP threads per Q-Chem process (default: 4)
#   --parallel-jobs N    Number of simultaneous jobs (default: 6)
#   --log FILE           Log file path (default: run_batch.log)
#   --max-jobs N         Limit jobs to run (0 = all, default: 0)
#   --dry-run            Print job list without executing
#
# Recommended parallelism (48 cores, 220 GB RAM):
#   4 threads x 12 parallel = 48 cores fully utilized
#   4 threads x 6 parallel  = 24 cores (leaves room for other work)
#   8 threads x 6 parallel  = 48 cores (better per-job performance)
#
# Background execution:
#   nohup bash run_batch.sh --parallel-jobs 12 --threads-per-job 4 > batch.log 2>&1 &
#   disown

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

JOBS_DIR="$SCRIPT_DIR/jobs"
THREADS_PER_JOB=4
PARALLEL_JOBS=6
LOG_FILE="$SCRIPT_DIR/run_batch.log"
MAX_JOBS=0
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --jobs-dir)        JOBS_DIR="$2";          shift 2;;
        --threads-per-job) THREADS_PER_JOB="$2";   shift 2;;
        --parallel-jobs)   PARALLEL_JOBS="$2";      shift 2;;
        --log)             LOG_FILE="$2";           shift 2;;
        --max-jobs)        MAX_JOBS="$2";           shift 2;;
        --dry-run)         DRY_RUN=1;               shift;;
        *) echo "Unknown option: $1"; exit 1;;
    esac
done

TOTAL_CORES=$(nproc)
USED_CORES=$((THREADS_PER_JOB * PARALLEL_JOBS))

echo "=== Q-Chem CCSD Amplitude Batch Runner ==="
echo "Jobs directory:     $JOBS_DIR"
echo "Threads per job:    $THREADS_PER_JOB"
echo "Parallel jobs:      $PARALLEL_JOBS"
echo "Total cores used:   $USED_CORES / $TOTAL_CORES"
echo "Log file:           $LOG_FILE"
echo ""

if [ "$USED_CORES" -gt "$TOTAL_CORES" ]; then
    echo "WARNING: Requesting $USED_CORES cores but only $TOTAL_CORES available!"
    echo ""
fi

PENDING_JOBS=()
DONE_COUNT=0
FAIL_COUNT=0

for JOB in "$JOBS_DIR"/*/; do
    [ -d "$JOB" ] || continue
    if [ -f "$JOB/.done" ]; then
        DONE_COUNT=$((DONE_COUNT + 1))
    elif [ -f "$JOB/.failed" ]; then
        FAIL_COUNT=$((FAIL_COUNT + 1))
        PENDING_JOBS+=("$JOB")
    else
        PENDING_JOBS+=("$JOB")
    fi
done

TOTAL_PENDING=${#PENDING_JOBS[@]}
echo "Total job directories: $((DONE_COUNT + FAIL_COUNT + TOTAL_PENDING))"
echo "Already completed:     $DONE_COUNT"
echo "Previously failed:     $FAIL_COUNT (will retry)"
echo "Pending:               $TOTAL_PENDING"
echo ""

if [ "$MAX_JOBS" -gt 0 ] && [ "$MAX_JOBS" -lt "$TOTAL_PENDING" ]; then
    PENDING_JOBS=("${PENDING_JOBS[@]:0:$MAX_JOBS}")
    echo "Limiting to first $MAX_JOBS jobs"
fi

if [ "$DRY_RUN" -eq 1 ]; then
    echo "--- DRY RUN: would process these jobs ---"
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
echo "Config: $PARALLEL_JOBS parallel x $THREADS_PER_JOB threads" >> "$LOG_FILE"

run_job() {
    local job_dir="$1"
    bash "$SCRIPT_DIR/run_single.sh" "$job_dir" "$THREADS_PER_JOB" 2>&1
}
export -f run_job
export SCRIPT_DIR THREADS_PER_JOB

printf '%s\n' "${PENDING_JOBS[@]}" | \
    xargs -P "$PARALLEL_JOBS" -I {} bash -c 'run_job "$@"' _ {} \
    2>&1 | tee -a "$LOG_FILE"

FINAL_DONE=$(find "$JOBS_DIR" -name ".done" | wc -l)
FINAL_FAIL=$(find "$JOBS_DIR" -name ".failed" | wc -l)
TOTAL_DIRS=$(find "$JOBS_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l)

echo ""
echo "=== Batch Complete ==="
echo "Completed: $FINAL_DONE / $TOTAL_DIRS"
echo "Failed:    $FINAL_FAIL / $TOTAL_DIRS"
echo "Finished at $(date)"
