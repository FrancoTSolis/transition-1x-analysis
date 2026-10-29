#!/bin/bash
# usage: chain_worker.sh <pid to wait for> <jobs file> <gpu>   -- start a worker once the given process has exited
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
while kill -0 "$1" 2>/dev/null; do sleep 60; done
exec ucj_init/worker.sh "$2" "$3"
