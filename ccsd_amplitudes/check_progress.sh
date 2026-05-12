#!/bin/bash
# Quick progress check for CCSD amplitude batch runs.
# Usage: bash check_progress.sh [jobs_dir]

JOBS_DIR="${1:-$(dirname "$0")/jobs}"

TOTAL=$(find "$JOBS_DIR" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | wc -l)
DONE=$(find "$JOBS_DIR" -name ".done" 2>/dev/null | wc -l)
FAILED=$(find "$JOBS_DIR" -name ".failed" 2>/dev/null | wc -l)
PENDING=$((TOTAL - DONE - FAILED))

echo "=== CCSD Amplitude Job Progress ==="
echo "Total:     $TOTAL"
echo "Completed: $DONE"
echo "Failed:    $FAILED"
echo "Pending:   $PENDING"

if [ "$TOTAL" -gt 0 ]; then
    PCT=$((DONE * 100 / TOTAL))
    echo "Progress:  ${PCT}%"
fi

echo ""

if [ "$DONE" -gt 0 ]; then
    echo "--- Last 5 completed ---"
    find "$JOBS_DIR" -name ".done" -printf '%T@ %h
' 2>/dev/null |         sort -rn | head -5 | while read ts dir; do
            name=$(basename "$dir")
            wall=$(cat "$dir/wall_time.txt" 2>/dev/null || echo "?")
            energy=$(cat "$dir/ccsd_energy.txt" 2>/dev/null || echo "?")
            echo "  $name  E=$energy  wall=${wall}s"
        done
fi

if [ "$FAILED" -gt 0 ]; then
    echo ""
    echo "--- Failed jobs ---"
    find "$JOBS_DIR" -name ".failed" 2>/dev/null | while read f; do
        echo "  $(basename "$(dirname "$f")")"
    done
fi

echo ""
AMP_SIZE=$(find "$JOBS_DIR" -name "*.dat" -exec du -ch {} + 2>/dev/null | tail -1 | awk '{print $1}')
TOTAL_SIZE=$(du -sh "$JOBS_DIR" 2>/dev/null | awk '{print $1}')
echo "Disk usage: ${TOTAL_SIZE:-0} total, ${AMP_SIZE:-0} amplitude files"
