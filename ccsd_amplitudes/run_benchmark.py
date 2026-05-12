#!/usr/bin/env python3
"""Benchmark Q-Chem CCSD at different thread counts."""

import os, subprocess, time

BENCH_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "benchmark")
os.makedirs(BENCH_DIR, exist_ok=True)

ENV_SCRIPT = "/xuanwu-tank/east/fts/qchem_compile/qcenv_dev.sh"
QCSCRATCH = "/xuanwu-tank/east/fts/qchem/scratch"
INPUT_FILE = os.path.join(BENCH_DIR, "bench_small.in")

THREAD_COUNTS = [1, 2, 4, 8, 12, 16, 24]

print("=== Q-Chem CCSD/STO-3G Thread Benchmark (7 atoms, 27 BFs) ===")
print(f"{'Threads':>8} {'Wall (s)':>10} {'Speedup':>8} {'CCSD Energy':>20} {'Amps':>6}")
print("-" * 60)

baseline = None
results = []

for nt in THREAD_COUNTS:
    scratch_name = f"bench_{nt}t"
    scratch_path = os.path.join(QCSCRATCH, scratch_name)
    out_file = os.path.join(BENCH_DIR, f"bench_{nt}t.out")

    subprocess.run(["rm", "-rf", scratch_path], capture_output=True)

    cmd = f"source {ENV_SCRIPT} 2>/dev/null && qchem -save -nt {nt} {INPUT_FILE} {out_file} {scratch_name}"

    start = time.time()
    subprocess.run(["bash", "-c", cmd], capture_output=True, text=True)
    elapsed = time.time() - start

    if baseline is None:
        baseline = elapsed
    speedup = baseline / elapsed if elapsed > 0 else 0

    ccsd_e = "N/A"
    amp_count = 0
    try:
        with open(out_file) as f:
            for line in f:
                if "CCSD total energy" in line:
                    ccsd_e = line.strip().split()[-1]
                if "Full" in line and "amplitudes" in line:
                    amp_count += 1
    except FileNotFoundError:
        pass

    print(f"{nt:>8} {elapsed:>10.1f} {speedup:>7.2f}x {ccsd_e:>20} {amp_count:>6}")
    results.append((nt, elapsed, speedup, ccsd_e, amp_count))

    subprocess.run(["rm", "-rf", scratch_path], capture_output=True)

print()
print("=== Throughput Analysis (48 cores total) ===")

for nt, elapsed, speedup, _, _ in results:
    efficiency = speedup / nt * 100
    parallel_jobs = 48 // nt
    throughput = parallel_jobs / elapsed * 3600
    print(f"  {nt:>2}t x {parallel_jobs:>2} parallel = {nt*parallel_jobs:>2} cores, "
          f"eff={efficiency:.0f}%, "
          f"~{throughput:.0f} jobs/hr")
