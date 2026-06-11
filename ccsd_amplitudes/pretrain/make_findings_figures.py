#!/usr/bin/env python3
"""Generate the figures backing docs/FINDINGS.md.

Uses the RHF diagnostic cache (rhf_cache/) and the training log
(train_pretrain.log). Outputs PNGs to docs/.
"""
import json
import re
from itertools import combinations
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score

OUT = Path("docs")
OUT.mkdir(exist_ok=True)
CACHE = Path("rhf_cache")


def load_cache():
    S = [np.load(f) for f in sorted(CACHE.glob("*.npz"))]
    # keep the dominant size for clean fixed-size analysis
    from collections import Counter
    norb = Counter(int(s["norb"]) for s in S).most_common(1)[0][0]
    S = [s for s in S if int(s["norb"]) == norb]
    return S, norb


def fig_learnability(S, norb):
    nocc = int(S[0]["nocc"]); nvirt = int(S[0]["nvirt"])
    n = len(S); ntr = int(0.8 * n)
    r2_J, r2_k = [], []
    for i in range(nocc):
        for j in range(i + 1, nocc):
            X = np.array([s["t2"][i, j].reshape(-1) for s in S])
            yJ = np.array([s["J"][0, 0, i, j] for s in S])
            yk = np.array([s["kappa_real"][0, i, j] for s in S])
            rJ = Ridge(0.1).fit(X[:ntr], yJ[:ntr])
            rk = Ridge(0.1).fit(X[:ntr], yk[:ntr])
            r2_J.append(r2_score(yJ[ntr:], rJ.predict(X[ntr:])))
            r2_k.append(r2_score(yk[ntr:], rk.predict(X[ntr:])))
    r2_J = np.clip(r2_J, -1, 1); r2_k = np.clip(r2_k, -1, 1)

    fig, ax = plt.subplots(figsize=(6, 4.2))
    ax.bar(["J / Z\n(diag Coulomb)", "kappa\n(log U)"],
           [np.mean(r2_J), np.mean(r2_k)],
           yerr=[np.std(r2_J), np.std(r2_k)], capsize=6,
           color=["#2a9d8f", "#e76f51"])
    ax.axhline(0, color="k", lw=0.8)
    ax.set_ylabel("local-linear test $R^2$\n(t2 slice → target entry)")
    ax.set_title(f"Diagonal-Coulomb J is learnable; orbital-rotation κ is not\n"
                 f"(RHF cache, {n} molecules, norb={norb})")
    ax.text(0, np.mean(r2_J) + 0.05, f"{np.mean(r2_J):+.2f}", ha="center", fontweight="bold")
    ax.text(1, 0.05, f"{np.mean(r2_k):+.2f}", ha="center", fontweight="bold")
    fig.tight_layout(); fig.savefig(OUT / "fig1_learnability.png", dpi=130); plt.close(fig)
    return np.mean(r2_J), np.mean(r2_k)


def fig_smoothness(S):
    nocc = int(S[0]["nocc"])
    t2 = [s["t2"].ravel() for s in S]
    kr = [s["kappa_real"][0].ravel() for s in S]
    J = [s["J"][0, 0].ravel() for s in S]
    dt, dk, dj = [], [], []
    for a, b in combinations(range(len(S)), 2):
        dt.append(np.linalg.norm(t2[a] - t2[b]))
        dk.append(np.linalg.norm(kr[a] - kr[b]))
        dj.append(np.linalg.norm(J[a] - J[b]))
    dt, dk, dj = map(np.array, (dt, dk, dj))
    ck = np.corrcoef(dt, dk)[0, 1]; cj = np.corrcoef(dt, dj)[0, 1]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, dy, c, name, col in [
        (axes[0], dj, cj, "J / Z (diag Coulomb)", "#2a9d8f"),
        (axes[1], dk, ck, "kappa (log U)", "#e76f51")]:
        ax.scatter(dt, dy, s=6, alpha=0.35, color=col)
        ax.set_xlabel(r"$\|\Delta t_2\|$ between molecules")
        ax.set_ylabel(rf"$\|\Delta\,${name}$\|$")
        ax.set_title(f"{name}\ncorr = {c:+.2f}  ({'smooth' if c > 0.3 else 'NOT smooth'})")
    fig.suptitle("Is the target a smooth function of t2? "
                 "Nearby t2 → nearby J, but NOT nearby κ", fontweight="bold")
    fig.tight_layout(); fig.savefig(OUT / "fig2_smoothness.png", dpi=130); plt.close(fig)
    return cj, ck


def fig_training():
    log = Path("train_pretrain.log")
    if not log.exists():
        return
    eps, vz, vk = [], [], []
    for line in log.read_text().splitlines():
        m = re.search(r"Epoch (\d+)/\d+.*val_z=([\d.eE+-]+).*val_kappa=([\d.eE+-]+)", line)
        if m:
            eps.append(int(m.group(1))); vz.append(float(m.group(2))); vk.append(float(m.group(3)))
    if not eps:
        return
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.plot(eps, vz, "o-", color="#2a9d8f", label="val Z/J loss (MSE)")
    ax.plot(eps, vk, "s-", color="#e76f51", label="val kappa loss (MSE)")
    ax.set_yscale("log"); ax.set_xlabel("epoch"); ax.set_ylabel("validation MSE (log)")
    ax.set_title("Supervised pretraining: Z/J converges low; κ plateaus (gauge wall)")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(OUT / "fig3_training.png", dpi=130); plt.close(fig)
    return eps[-1], vz[-1], vk[-1]


def main():
    S, norb = load_cache()
    rJ, rk = fig_learnability(S, norb)
    cj, ck = fig_smoothness(S)
    tr = fig_training()
    print(json.dumps({
        "n_cache": len(S), "norb": norb,
        "learnability_J_R2": round(rJ, 3), "learnability_kappa_R2": round(rk, 3),
        "smooth_corr_J": round(cj, 3), "smooth_corr_kappa": round(ck, 3),
        "last_epoch": tr,
    }, indent=2))


if __name__ == "__main__":
    main()
