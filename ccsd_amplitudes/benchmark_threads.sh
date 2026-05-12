#!/bin/bash
# Benchmark Q-Chem CCSD with PRINT_QIS at different thread counts.
# Uses a small molecule (7 atoms, STO-3G) for quick iteration.
#
# Usage: bash benchmark_threads.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BENCH_DIR="$SCRIPT_DIR/benchmark"
mkdir -p "$BENCH_DIR"

source /xuanwu-tank/east/fts/qchem_compile/qcenv_dev.sh

cat > "$BENCH_DIR/bench_small.in" << 'EOF'
$molecule
0 1
O    0.23947457   -1.34716795   -0.12501468
C    0.97146208   -0.34981609   -0.13612089
N    0.58929970    0.93141869   -0.08536160
N   -0.90628544    1.13570298   -0.63554299
C   -1.21230917    0.26282494    0.11684332
H    2.06966965   -0.41696990   -0.14543590
H   -1.78569220   -0.31935637    0.81484764
$end

$rem
    METHOD          ccsd
    BASIS           STO-3G
    UNRESTRICTED    TRUE
    INTERNAL_STABILITY TRUE
    SCF_ALGORITHM   DIIS_GDM
    MAX_SCF_CYCLES  200
    THRESH          15
    SCF_CONVERGENCE 10
    CC_CONVERGENCE  7
    N_FROZEN_CORE   FC
    GUI             2
    PRINT_QIS       TRUE
$end
EOF

echo "=== Benchmarking Q-Chem CCSD/STO-3G (7 atoms, 27 BFs) ==="
echo "Machine: $(nproc) cores, $(free -h | awk '/Mem:/{print $2}') RAM"
echo ""

THREAD_COUNTS="1 2 4 8 12 16 24"
echo "Thread count | Wall time (s) | Speedup"
echo "-------------|---------------|--------"

BASELINE=""
for NT in $THREAD_COUNTS; do
    SCRATCH_NAME="bench_${NT}t"
    rm -rf "$QCSCRATCH/$SCRATCH_NAME" 2>/dev/null || true

    START=$(date +%s%N)
    qchem -save -nt "$NT" "$BENCH_DIR/bench_small.in" "$BENCH_DIR/bench_${NT}t.out" "$SCRATCH_NAME" 2>/dev/null
    END=$(date +%s%N)
    ELAPSED=$(( (END - START) / 1000000 ))

    if [ -z "$BASELINE" ]; then
        BASELINE=$ELAPSED
    fi

    CCSD_E=$(grep "CCSD total energy" "$BENCH_DIR/bench_${NT}t.out" | tail -1 | awk '{print $NF}')
    AMP_OK=$(grep -c "Full.*amplitudes" "$BENCH_DIR/bench_${NT}t.out" || echo 0)

    python3 -c "
elapsed_ms = $ELAPSED
baseline_ms = $BASELINE
s = elapsed_ms / 1000
sp = baseline_ms / elapsed_ms if elapsed_ms > 0 else 0
print(f'${NT:>13} | {s:>13.1f} | {sp:.2f}x  (E=${CCSD_E}, ${AMP_OK} amp files)')
"

    rm -rf "$QCSCRATCH/$SCRATCH_NAME" 2>/dev/null || true
done

echo ""
echo "=== Scratch size after a single CCSD/STO-3G (7 atoms) ==="
SCRATCH_NAME="bench_size_test"
rm -rf "$QCSCRATCH/$SCRATCH_NAME" 2>/dev/null || true
qchem -save -nt 4 "$BENCH_DIR/bench_small.in" "$BENCH_DIR/bench_size.out" "$SCRATCH_NAME" 2>/dev/null
du -sh "$QCSCRATCH/$SCRATCH_NAME"
echo "Breakdown of large files:"
du -sh "$QCSCRATCH/$SCRATCH_NAME"/* 2>/dev/null | sort -rh | head -10
echo ""
echo "Amplitude files:"
ls -lh "$QCSCRATCH/$SCRATCH_NAME"/*t*.dat 2>/dev/null

echo ""
echo "=== Done ==="
