#!/bin/bash
# One-time Expanse setup for the ML-LUCJ project (run ON the login node, in the background):
#   python3 expanse.py put setup_env.sh /home/fsun3/ml_lucj_setup_env.sh
#   python3 expanse.py run "nohup bash /home/fsun3/ml_lucj_setup_env.sh > /home/fsun3/ml_lucj_setup_env.log 2>&1 &"
# Creates our own directories and a separate conda env prefix; does not touch the Fermi-arc
# project's directories or its `negf` env.
set -euo pipefail
ROOT=/expanse/lustre/scratch/fsun3/temp_project/ml_lucj
ENV=/home/fsun3/envs/ml_lucj
mkdir -p "$ROOT"/{code,data,runs,logs,results} /home/fsun3/envs

source /etc/profile.d/modules.sh
module purge
module load cpu/0.17.3b gcc/10.2.0/npcyll4 anaconda3/2021.05/q4munrg
source "$(conda info --base)/etc/profile.d/conda.sh"

if [ ! -x "$ENV/bin/python" ]; then
  echo "[setup] creating env $ENV $(date)"
  conda create -y -p "$ENV" -c conda-forge python=3.11 pip
fi
conda activate "$ENV"
python -m pip install -q --upgrade pip
python -m pip install -q numpy scipy
python -m pip install -q torch --index-url https://download.pytorch.org/whl/cpu
python -m pip install -q ffsim pyscf quimb cotengra autoray qiskit qiskit-quimb qiskit-addon-sqd jax jaxlib opt_einsum
python - <<'EOF'
import ffsim, pyscf, quimb, torch, jax, qiskit_addon_sqd, qiskit_quimb
print("[setup] OK torch", torch.__version__, "jax", jax.__version__, "quimb", quimb.__version__, "pyscf", pyscf.__version__)
EOF
echo "[setup] done $(date)"
