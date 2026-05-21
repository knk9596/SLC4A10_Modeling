#!/usr/bin/env python3
"""Save Na analysis as CSV and generate plots.

Na1 is assigned from chain B, Na2 from chain C in AF3 predictions.
"""

import glob
import os
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

ALIGNED_DIRS = [
    "/sc/arion/scratch/zhaoy14/ALineMol/AF3_ion/aligned",
    "/sc/arion/scratch/zhaoy14/ALineMol/AF3_ion/randomseed_1100/aligned",
]
REF_PDB = "/sc/arion/scratch/zhaoy14/4A10_revision/0316/exp_chainB.pdb"
OUT_DIR = "/sc/arion/scratch/zhaoy14/ALineMol/AF3_ion/randomseed_1100/results"

os.makedirs(OUT_DIR, exist_ok=True)


def parse_na_atoms(pdb_path):
    na_atoms = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("HETATM") and " NA " in line and "NA" in line[17:20].strip():
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                chain = line[21].strip()
                resseq = int(line[22:26])
                bfactor = float(line[60:66])
                na_atoms.append({
                    "chain": chain, "resseq": resseq,
                    "coords": np.array([x, y, z]),
                    "bfactor": bfactor
                })
    return na_atoms


def parse_ca_atoms(pdb_path):
    """Extract CA atom coordinates keyed by residue number from a PDB file."""
    cas = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                resseq = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                cas[resseq] = np.array([x, y, z])
    return cas


def compute_ca_rmsd(ref_cas, mov_cas):
    """Compute RMSD between two dicts of CA coordinates matched by residue number."""
    common = sorted(set(ref_cas) & set(mov_cas))
    if len(common) == 0:
        return float("inf")
    ref = np.array([ref_cas[r] for r in common])
    mov = np.array([mov_cas[r] for r in common])
    return np.sqrt(np.mean(np.sum((ref - mov) ** 2, axis=1)))


# Reference — Na1: chain C res 6, Na2: chain C res 7
ref_nas = parse_na_atoms(REF_PDB)
_ref_by_key = {(n["chain"], n["resseq"]): n for n in ref_nas}
ref_na1 = _ref_by_key[("C", 6)]["coords"]
ref_na2 = _ref_by_key[("C", 7)]["coords"]
ref_cas = parse_ca_atoms(REF_PDB)

# Collect all data from both directories
aligned_files = sorted(
    f for d in ALIGNED_DIRS for f in glob.glob(os.path.join(d, "*_aligned.pdb"))
)
rows = []

for pdb_file in aligned_files:
    name = os.path.basename(pdb_file).replace("_aligned.pdb", "")
    nas = parse_na_atoms(pdb_file)
    if len(nas) < 2:
        continue

    # Compute backbone RMSD
    mov_cas = parse_ca_atoms(pdb_file)
    backbone_rmsd = compute_ca_rmsd(ref_cas, mov_cas)

    # Assign Na1/Na2 by chain (B = Na1, C = Na2), then fall back to proximity
    na_b = [n for n in nas if n["chain"] == "B"]
    na_c = [n for n in nas if n["chain"] == "C"]

    if na_b and na_c:
        pred_na1 = min(na_b, key=lambda n: np.linalg.norm(n["coords"] - ref_na1))
        pred_na2 = min(na_c, key=lambda n: np.linalg.norm(n["coords"] - ref_na2))
    else:
        # Fall back to proximity-based assignment
        d1_to_ref1 = np.linalg.norm(nas[0]["coords"] - ref_na1)
        d1_to_ref2 = np.linalg.norm(nas[0]["coords"] - ref_na2)
        d2_to_ref1 = np.linalg.norm(nas[1]["coords"] - ref_na1)
        d2_to_ref2 = np.linalg.norm(nas[1]["coords"] - ref_na2)
        if d1_to_ref1 + d2_to_ref2 < d1_to_ref2 + d2_to_ref1:
            pred_na1, pred_na2 = nas[0], nas[1]
        else:
            pred_na1, pred_na2 = nas[1], nas[0]

    dist_na1 = np.linalg.norm(pred_na1["coords"] - ref_na1)
    dist_na2 = np.linalg.norm(pred_na2["coords"] - ref_na2)

    rows.append({
        "model": name,
        "backbone_rmsd": backbone_rmsd,
        "na1_x": pred_na1["coords"][0],
        "na1_y": pred_na1["coords"][1],
        "na1_z": pred_na1["coords"][2],
        "na1_dist_to_exp": dist_na1,
        "na1_plddt": pred_na1["bfactor"],
        "na2_x": pred_na2["coords"][0],
        "na2_y": pred_na2["coords"][1],
        "na2_z": pred_na2["coords"][2],
        "na2_dist_to_exp": dist_na2,
        "na2_plddt": pred_na2["bfactor"],
    })

# ── Save CSV ──
csv_path = os.path.join(OUT_DIR, "na_analysis.csv")
fields = list(rows[0].keys())
with open(csv_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)
print(f"Saved CSV: {csv_path}")

# ── Prepare arrays ──
na1_coords = np.array([[r["na1_x"], r["na1_y"], r["na1_z"]] for r in rows])
na2_coords = np.array([[r["na2_x"], r["na2_y"], r["na2_z"]] for r in rows])
na1_dists = np.array([r["na1_dist_to_exp"] for r in rows])
na2_dists = np.array([r["na2_dist_to_exp"] for r in rows])

# ── Plot 1: 3D scatter of Na1 and Na2 coordinates ──
fig = plt.figure(figsize=(14, 6))

ax1 = fig.add_subplot(121, projection="3d")
ax1.scatter(*na1_coords.T, c="dodgerblue", s=50, alpha=0.7, label="AF3 Na1")
ax1.scatter(*ref_na1, c="red", s=200, marker="*", label="Exp Na1", zorder=5)
ax1.set_title("Na1 positions (tight cluster)", fontsize=12, fontweight="bold")
ax1.set_xlabel("X (Å)")
ax1.set_ylabel("Y (Å)")
ax1.set_zlabel("Z (Å)")
ax1.legend()

ax2 = fig.add_subplot(122, projection="3d")
ax2.scatter(*na2_coords.T, c="orange", s=50, alpha=0.7, label="AF3 Na2")
ax2.scatter(*ref_na2, c="red", s=200, marker="*", label="Exp Na2", zorder=5)
ax2.set_title("Na2 positions (dispersed)", fontsize=12, fontweight="bold")
ax2.set_xlabel("X (Å)")
ax2.set_ylabel("Y (Å)")
ax2.set_zlabel("Z (Å)")
ax2.legend()

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "na_3d_scatter.png"), dpi=300, bbox_inches="tight")
print(f"Saved: na_3d_scatter.png")
plt.close()

# ── Plot 2: Distance from experimental site (bar chart) ──
fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=False)
models = [r["model"].replace("seed-", "s").replace("_sample-", "_") for r in rows]
x = np.arange(len(models))

# Na1
bars1 = axes[0].bar(x, na1_dists, color="dodgerblue", alpha=0.8, edgecolor="navy", linewidth=0.5)
axes[0].axhline(y=1.5, color="green", linestyle="--", alpha=0.6, label="1.5 Å threshold")
axes[0].set_title("Na1: Distance from Experimental Site", fontsize=12, fontweight="bold")
axes[0].set_ylabel("Distance (Å)")
axes[0].set_xticks(x)
axes[0].set_xticklabels(models, rotation=90, fontsize=6)
axes[0].legend()
axes[0].set_ylim(0, max(na1_dists) * 1.2)

# Na2
bars2 = axes[1].bar(x, na2_dists, color="orange", alpha=0.8, edgecolor="darkorange", linewidth=0.5)
axes[1].axhline(y=3.0, color="green", linestyle="--", alpha=0.6, label="3.0 Å threshold")
axes[1].set_title("Na2: Distance from Experimental Site", fontsize=12, fontweight="bold")
axes[1].set_ylabel("Distance (Å)")
axes[1].set_xticks(x)
axes[1].set_xticklabels(models, rotation=90, fontsize=6)
axes[1].legend()

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "na_displacement_bars.png"), dpi=300, bbox_inches="tight")
print(f"Saved: na_displacement_bars.png")
plt.close()

# ── Plot 3: Combined violin + strip plot (all models) ──
n_all = len(rows)

fig, ax = plt.subplots(figsize=(7, 6))
data = [na1_dists, na2_dists]
parts = ax.violinplot(data, positions=[1, 2], showmedians=True, showextrema=False)
parts["bodies"][0].set_facecolor("dodgerblue")
parts["bodies"][0].set_alpha(0.4)
parts["bodies"][1].set_facecolor("orange")
parts["bodies"][1].set_alpha(0.4)

# Overlay individual points
np.random.seed(42)
jitter1 = 1 + np.random.normal(0, 0.04, len(na1_dists))
jitter2 = 2 + np.random.normal(0, 0.04, len(na2_dists))
ax.scatter(jitter1, na1_dists, c="dodgerblue", s=30, alpha=0.7, edgecolors="navy", linewidths=0.5, zorder=3)
ax.scatter(jitter2, na2_dists, c="orange", s=30, alpha=0.7, edgecolors="darkorange", linewidths=0.5, zorder=3)

ax.set_xticks([1, 2])
ax.set_xticklabels(["Na1", "Na2"], fontsize=13)
ax.set_ylabel("Distance from Experimental Position (Å)", fontsize=11)
ax.set_title("AF3 Sodium Ion Prediction Without Anion", fontsize=13, fontweight="bold")
ax.axhline(y=1.5, color="gray", linestyle=":", alpha=0.5)

# Annotations with stats over all models
na1_txt = (f"n={n_all}\n"
           f"{(na1_dists < 1.5).sum()}/{n_all} within 1.5 Å\n"
           f"{na1_dists.mean():.2f} ± {na1_dists.std():.2f} Å")
na2_txt = (f"n={n_all}\n"
           f"{(na2_dists < 3.0).sum()}/{n_all} within 3.0 Å\n"
           f"{(na2_dists > 10).sum()}/{n_all} displaced >10 Å\n"
           f"{na2_dists.mean():.1f} ± {na2_dists.std():.1f} Å")
y_top = max(na1_dists.max(), na2_dists.max())
ax.annotate(na1_txt, xy=(1, na1_dists.mean()), xytext=(1, y_top * 0.72),
            fontsize=8.5, color="navy", ha="center",
            bbox=dict(boxstyle="round,pad=0.3", fc="aliceblue", ec="navy", alpha=0.8))
ax.annotate(na2_txt, xy=(2, na2_dists.mean()), xytext=(2.6, y_top * 0.72),
            fontsize=8.5, color="darkorange", ha="center",
            bbox=dict(boxstyle="round,pad=0.3", fc="lemonchiffon", ec="darkorange", alpha=0.8))

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "na_violin.png"), dpi=300, bbox_inches="tight")
print(f"Saved: na_violin.png")
plt.close()

