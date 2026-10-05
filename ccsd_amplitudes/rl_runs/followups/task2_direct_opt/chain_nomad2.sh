#!/bin/bash
# Resume (5 Oct, successor agent): the objective-A NOMAD jobs of jobs_A_ext.txt (lines 16-27), reordered so that the
# label starts of all molecules come first, run one at a time on scai2 GPU 3 after the job the old queue was running
# (C2H4O_rxn2507_TS rl4n29f NOMAD) has finished.  No job is started after the deadline (finish within the ~4 h window).
cd /xuanwu-tank/east/fts/projects/transition-1x-analysis/ccsd_amplitudes
D=rl_runs/followups/task2_direct_opt
DEADLINE=$(date -d "2026-10-05 13:15" +%s)
while pgrep -f "task2_opt --name C2H4O_rxn2507_TS --start rl4n29f --objective lucj --optimizer nomad" > /dev/null; do sleep 20; done
n=0
while IFS= read -r line; do
  [[ -z "$line" || "$line" == \#* ]] && continue
  if [ "$(date +%s)" -gt "$DEADLINE" ]; then echo "deadline passed $(date '+%F %T'); not started: $line"; continue; fi
  n=$((n + 1)); f=$D/.job_nomad2_$n.txt; echo "$line" > $f
  bash pretrain/followups/task2_queue2.sh $f < /dev/null
  rm -f $f
done < $D/jobs_A_nomad2.txt
echo "nomad2 queue done $(date '+%F %T')"
