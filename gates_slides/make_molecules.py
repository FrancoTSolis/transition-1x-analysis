#!/usr/bin/env python3
"""Render real molecules from the training set as clean deck illustrations.

Everything here is drawn from the actual geometries in ccsd_amplitudes/jobs,
so the pictures on the slides are the molecules the model was trained and
tested on -- not stock art.

Outputs (art/):
  one_molecule.png   a single molecule, drawn in full detail  -> "quantum"
  many_molecules.png a grid of different molecules             -> "AI at scale"
  reaction.png       reactant -> transition state -> product,
                     with the bond that breaks marked          -> the case

Run from gates_slides/:  python3 make_molecules.py
"""
import random
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

JOBS = Path(__file__).resolve().parents[1] / "ccsd_amplitudes" / "jobs"
OUT = Path(__file__).resolve().parent / "art"

COV = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "S": 1.05}
COLOR = {"H": "#ffffff", "C": "#44545f", "N": "#2f6fb3", "O": "#d1503c",
         "F": "#5aa469", "S": "#c8a020"}
EDGE = {"H": "#8b9aa6", "C": "#2d3a44", "N": "#1d4f86", "O": "#9c3527",
        "F": "#3d7a4b", "S": "#8f7317"}
SIZE = {"H": 0.30, "C": 0.52, "N": 0.50, "O": 0.49, "F": 0.47, "S": 0.62}
BOND = "#9aa9b5"


def read_xyz(path: Path):
    lines = path.read_text().splitlines()
    n = int(lines[0].split()[0])
    syms, xyz = [], []
    for line in lines[2:2 + n]:
        p = line.split()
        syms.append(p[0])
        xyz.append([float(v) for v in p[1:4]])
    return syms, np.asarray(xyz)


def bond_list(syms, xyz, slack=1.30):
    out = []
    for i in range(len(syms)):
        for j in range(i + 1, len(syms)):
            d = np.linalg.norm(xyz[i] - xyz[j])
            if d < slack * (COV.get(syms[i], .7) + COV.get(syms[j], .7)):
                out.append((i, j, d))
    return out


def frame(xyz):
    """Principal-axis view: widest spread in the plane, depth along the third."""
    c = xyz - xyz.mean(0)
    _, _, vt = np.linalg.svd(c, full_matrices=False)
    return vt


def kabsch(mobile, target):
    """Rotation aligning mobile onto target (same atom order)."""
    a = mobile - mobile.mean(0)
    b = target - target.mean(0)
    u, _, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    return u @ np.diag([1, 1, d]) @ vt


def draw(ax, syms, xyz, basis, scale=1.0, mark=None, lw=6.0):
    p = (xyz - xyz.mean(0)) @ basis.T
    order = np.argsort(p[:, 2])
    for i, j, _ in bond_list(syms, xyz):
        ax.plot([p[i, 0], p[j, 0]], [p[i, 1], p[j, 1]], color=BOND,
                lw=lw * scale, solid_capstyle="round", zorder=1)
    if mark is not None:
        i, j = mark
        ax.plot([p[i, 0], p[j, 0]], [p[i, 1], p[j, 1]], color="#e0863a",
                lw=lw * scale * 1.05, ls=(0, (2.2, 1.6)),
                solid_capstyle="round", zorder=2)
    for k in order:
        s = syms[k]
        ax.add_patch(Circle((p[k, 0], p[k, 1]), SIZE.get(s, .5) * scale,
                            facecolor=COLOR.get(s, "#888"),
                            edgecolor=EDGE.get(s, "#555"),
                            lw=1.6 * scale, zorder=3 + k * 1e-3))
    ax.set_aspect("equal")
    ax.axis("off")
    return p


def pick(name):
    return read_xyz(JOBS / name / f"{name}.xyz")


def one_molecule(name):
    syms, xyz = pick(name)
    fig, ax = plt.subplots(figsize=(4.4, 3.3))
    p = draw(ax, syms, xyz, frame(xyz), scale=1.0, lw=7.0)
    pad = 1.0
    ax.set_xlim(p[:, 0].min() - pad, p[:, 0].max() + pad)
    ax.set_ylim(p[:, 1].min() - pad, p[:, 1].max() + pad)
    fig.tight_layout(pad=0.1)
    fig.savefig(OUT / "one_molecule.png", dpi=220, transparent=True)
    plt.close(fig)
    print("one_molecule:", name, len(syms), "atoms")


def many_molecules(names, cols=5, rows=3):
    fig, axes = plt.subplots(rows, cols, figsize=(4.9, 3.3))
    for ax, name in zip(axes.ravel(), names):
        syms, xyz = pick(name)
        p = draw(ax, syms, xyz, frame(xyz), scale=0.72, lw=4.4)
        pad = 0.9
        ax.set_xlim(p[:, 0].min() - pad, p[:, 0].max() + pad)
        ax.set_ylim(p[:, 1].min() - pad, p[:, 1].max() + pad)
    for ax in axes.ravel()[len(names):]:
        ax.axis("off")
    fig.subplots_adjust(0, 0, 1, 1, 0.04, 0.04)
    fig.savefig(OUT / "many_molecules.png", dpi=220, transparent=True)
    plt.close(fig)
    print("many_molecules:", len(names))


def reaction(stem):
    syms, ts = pick(f"{stem}_TS")
    _, r = pick(f"{stem}_R")
    _, p_ = pick(f"{stem}_P")
    basis = frame(ts)
    r = r @ kabsch(r, ts)
    p_ = p_ @ kabsch(p_, ts)

    br = {(i, j) for i, j, _ in bond_list(syms, r)}
    bp = {(i, j) for i, j, _ in bond_list(syms, p_)}
    changed = sorted(br ^ bp, key=lambda e: -abs(
        np.linalg.norm(r[e[0]] - r[e[1]]) - np.linalg.norm(p_[e[0]] - p_[e[1]])))
    mark = changed[0] if changed else None

    # one shared scale so the three frames are visually comparable
    half = 1.15 + max(
        np.abs((g - g.mean(0)) @ basis.T)[:, :2].max() for g in (r, ts, p_))

    fig, axes = plt.subplots(1, 3, figsize=(9.4, 3.1))
    titles = ["reactant", "transition state", "product"]
    for ax, geom, title, m in zip(axes, [r, ts, p_], titles,
                                  [None, mark, None]):
        draw(ax, syms, geom, basis, scale=0.92, mark=m, lw=6.2)
        ax.set_xlim(-half, half)
        ax.set_ylim(-half * 0.72, half * 0.72)
        ax.set_title(title, fontsize=14.5, color="#10344d",
                     fontweight="bold" if title == "transition state" else "normal",
                     pad=1)
    for x in (0.352, 0.672):
        fig.text(x, 0.45, "→", fontsize=26, color="#9aa9b5",
                 ha="center", va="center")
    fig.subplots_adjust(0.005, 0.02, 0.995, 0.91, 0.10, 0)
    fig.savefig(OUT / "reaction.png", dpi=220, transparent=True)
    plt.close(fig)
    print("reaction:", stem, "marked bond:", mark)


GRID = [
    "C2H3NO_rxn3841_R", "C4HNO_rxn3724_R", "C2H2N2O_rxn2091_R",
    "C3H2N2_rxn1185_R", "C2H2N2O2_rxn3534_R", "C2H2N2O_rxn2545_R",
    "C4HNO_rxn0475_R", "C2H3NO_rxn3839_P", "C2H2N2O_rxn3121_R",
    "C2H2N2O2_rxn3878_R", "C4HNO_rxn8987_R", "C2H2N2O_rxn3125_P",
    "C2H3NO_rxn3842_R", "C2H2N2O_rxn9455_R", "C4HNO_rxn3727_P",
]


def main():
    OUT.mkdir(exist_ok=True)
    one_molecule("C2H2N2O_rxn2091_R")
    reaction("C2H3NO_rxn3841")
    many_molecules(GRID)


if __name__ == "__main__":
    main()
