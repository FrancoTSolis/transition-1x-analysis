#!/bin/bash
# Run a single Q-Chem CCSD job with PRINT_QIS and collect amplitude files.
#
# Usage: bash run_single.sh <job_dir> [threads]
# Example: bash run_single.sh jobs/C2H2N2O_rxn2091_TS 1
#
# Uses local scratch (/tmp) instead of NFS for I/O performance.
# Runs Q-Chem from inside the job directory so .fchk lands there.

set -euo pipefail

JOB_DIR="$(cd "$1" && pwd)"
THREADS="${2:-1}"

if [ ! -d "$JOB_DIR" ]; then
    echo "ERROR: Job directory not found: $JOB_DIR"
    exit 1
fi

INPUT_FILE=$(ls "$JOB_DIR"/*.in 2>/dev/null | head -1)
if [ -z "$INPUT_FILE" ]; then
    echo "ERROR: No .in file found in $JOB_DIR"
    exit 1
fi

JOB_NAME=$(basename "$JOB_DIR")
INPUT_BASENAME=$(basename "$INPUT_FILE")
OUTPUT_FILE="$JOB_DIR/${JOB_NAME}.out"
SCRATCH_NAME="ccsd_amp_${JOB_NAME}"

source /xuanwu-tank/east/fts/qchem_compile/qcenv_dev.sh 2>/dev/null

# Use local disk for scratch -- much faster than NFS
export QCSCRATCH="/tmp/qchem_scratch"
mkdir -p "$QCSCRATCH"

DONE_MARKER="$JOB_DIR/.done"
FAIL_MARKER="$JOB_DIR/.failed"

if [ -f "$DONE_MARKER" ]; then
    echo "SKIP: $JOB_NAME already completed"
    exit 0
fi

rm -f "$FAIL_MARKER"
rm -rf "$QCSCRATCH/$SCRATCH_NAME" 2>/dev/null || true

# cd into job dir so Q-Chem writes .fchk there
cd "$JOB_DIR"

echo "$(date +%Y-%m-%dT%H:%M:%S) START $JOB_NAME (threads=$THREADS)"

START_TIME=$(date +%s)
if qchem -save -nt "$THREADS" "$INPUT_BASENAME" "${JOB_NAME}.out" "$SCRATCH_NAME" 2>/dev/null; then
    END_TIME=$(date +%s)
    WALL=$(( END_TIME - START_TIME ))

    CCSD_E=$(grep "CCSD total energy" "${JOB_NAME}.out" 2>/dev/null | tail -1 | awk '{print $NF}')
    CCSD_E="${CCSD_E:-N/A}"

    for DATNAME in ccsd_t1.dat ccsd_t2.dat mp2_t2.dat mp2_t1.dat; do
        DATPATH="$QCSCRATCH/$SCRATCH_NAME/$DATNAME"
        if [ -f "$DATPATH" ]; then
            cp "$DATPATH" .
        fi
    done

    echo "$CCSD_E" > ccsd_energy.txt
    echo "$WALL" > wall_time.txt
    touch .done

    echo "$(date +%Y-%m-%dT%H:%M:%S) DONE  $JOB_NAME  E=$CCSD_E  wall=${WALL}s"
else
    END_TIME=$(date +%s)
    WALL=$(( END_TIME - START_TIME ))
    echo "QCHEM_FAILURE" > .failed
    echo "$(date +%Y-%m-%dT%H:%M:%S) FAIL  $JOB_NAME  wall=${WALL}s"
fi

rm -rf "$QCSCRATCH/$SCRATCH_NAME" 2>/dev/null || true
