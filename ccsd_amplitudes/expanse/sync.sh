#!/bin/bash
# Mirror the code + small-molecule data needed for energy/RL jobs to Expanse Lustre.
#   bash expanse/sync.sh [extra local paths relative to ccsd_amplitudes/ ...]
# Remote layout: /expanse/lustre/scratch/fsun3/temp_project/ml_lucj/code/ccsd_amplitudes/
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
CA=$(cd "$HERE/.." && pwd)
REMOTE=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj/code/ccsd_amplitudes
cd "$CA"
LIST=$(mktemp)
{
  ls pretrain/*.py pretrain/rl/*.py pretrain/opt_true/*.py pretrain/opt_true/*.txt pretrain/opt_true/*.npz 2>/dev/null
  ls gauge_study/__init__.py gauge_study/common.py gauge_study/compressed_canonical.py gauge_study/names_*.txt
  ls expanse/*.py expanse/*.sh expanse/*.slurm 2>/dev/null || true
  # small molecules (norb <= 20): inputs + cached active-space Hamiltonians
  python3 - <<'EOF'
import json
idx = json.load(open("rhf_dataset/_index.json"))
for k, v in sorted(idx.items()):
    if v[0] <= 20:
        print(f"rhf_dataset/{k}.npz")
EOF
  ls rhf_hamiltonians/*.npz
  for x in "$@"; do echo "$x"; done
} | sort -u > "$LIST"
python3 - "$LIST" <<'EOF'
import json, sys
# remote index restricted to the synced molecules (CCSDAmplitudeDataset reads _index.json)
idx = json.load(open("rhf_dataset/_index.json"))
names = {l.strip()[len("rhf_dataset/"):-4] for l in open(sys.argv[1]) if l.startswith("rhf_dataset/") and l.strip().endswith(".npz")}
json.dump({k: v for k, v in idx.items() if k in names}, open("/tmp/_index_small.json", "w"))
EOF
python3 "$HERE/expanse.py" run "mkdir -p $REMOTE/rhf_dataset $REMOTE/rhf_hamiltonians"
python3 "$HERE/expanse.py" put "$CA/" "$REMOTE/" --files-from="$LIST"
python3 "$HERE/expanse.py" put /tmp/_index_small.json "$REMOTE/rhf_dataset/_index.json"
echo "synced $(wc -l < "$LIST") files -> $REMOTE"
rm -f "$LIST"
