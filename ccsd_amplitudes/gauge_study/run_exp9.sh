#!/bin/bash
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
exec taskset -c 0-27 .lucj_venv/bin/python3 -u -m gauge_study.exp9_reg_energy_sweep --n-procs 24 --connectivities square all-to-all
