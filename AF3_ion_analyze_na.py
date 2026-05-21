#!/usr/bin/env python3
"""Analyze Na1 and Na2 positions across all aligned AF3 models vs experimental reference.

Na1 is assigned from chain B, Na2 from chain C in AF3 predictions.
"""

import glob
import os
import numpy as np

ALIGNED_DIRS = [
    "/sc/arion/scratch/zhaoy14/ALineMol/AF3_ion/aligned",
    "/sc/arion/scratch/zhaoy14/ALineMol/AF3_ion/randomseed_1100/aligned",
]
REF_PDB = "/sc/arion/scratch/zhaoy14/4A10_revision/0316/exp_chainB.pdb"


def parse_na_atoms(pdb_path):
    """Extract Na atom coordinates from a PDB file."""
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


# Reference Na positions — Na1: chain C res 6, Na2: chain C res 7
ref_nas = parse_na_atoms(REF_PDB)
_ref_by_key = {(n["chain"], n["resseq"]): n for n in ref_nas}
ref_na1_atom = _ref_by_key[("C", 6)]
ref_na2_atom = _ref_by_key[("C", 7)]

print("=" * 70)
print("EXPERIMENTAL REFERENCE Na POSITIONS")
print("=" * 70)
for label, na in [("Na1", ref_na1_atom), ("Na2", ref_na2_atom)]:
    print(f"  {label} (chain {na['chain']}, res {na['resseq']}): "
          f"({na['coords'][0]:.3f}, {na['coords'][1]:.3f}, {na['coords'][2]:.3f})  "
          f"B-factor: {na['bfactor']:.2f}")

ref_na1 = ref_na1_atom["coords"]
ref_na2 = ref_na2_atom["coords"]

print(f"\n  Distance between Na1 and Na2: {np.linalg.norm(ref_na1 - ref_na2):.2f} Å")

# Analyze all aligned models from both directories
aligned_files = sorted(
    f for d in ALIGNED_DIRS for f in glob.glob(os.path.join(d, "*_aligned.pdb"))
)

na1_dists = []
na2_dists = []
na1_coords_all = []
na2_coords_all = []
results = []

print("\n" + "=" * 70)
print("AF3 PREDICTED Na POSITIONS (aligned to experimental)")
print("=" * 70)
print(f"{'Model':<28} {'Na1 dist (Å)':>12} {'Na2 dist (Å)':>12} {'Na1 B-fac':>10} {'Na2 B-fac':>10}")
print("-" * 70)

for pdb_file in aligned_files:
    name = os.path.basename(pdb_file).replace("_aligned.pdb", "")
    nas = parse_na_atoms(pdb_file)

    if len(nas) < 2:
        print(f"{name:<28} {'< 2 Na atoms found':>40}")
        continue

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

    na1_dists.append(dist_na1)
    na2_dists.append(dist_na2)
    na1_coords_all.append(pred_na1["coords"])
    na2_coords_all.append(pred_na2["coords"])

    results.append((name, dist_na1, dist_na2, pred_na1["bfactor"], pred_na2["bfactor"]))
    print(f"{name:<28} {dist_na1:>12.3f} {dist_na2:>12.3f} {pred_na1['bfactor']:>10.2f} {pred_na2['bfactor']:>10.2f}")

# Summary statistics
print("\n" + "=" * 70)
print("SUMMARY STATISTICS")
print("=" * 70)
na1_dists = np.array(na1_dists)
na2_dists = np.array(na2_dists)

print(f"\nNa1 displacement from experimental site:")
print(f"  Mean: {na1_dists.mean():.3f} ± {na1_dists.std():.3f} Å")
print(f"  Min:  {na1_dists.min():.3f} Å")
print(f"  Max:  {na1_dists.max():.3f} Å")

print(f"\nNa2 displacement from experimental site:")
print(f"  Mean: {na2_dists.mean():.3f} ± {na2_dists.std():.3f} Å")
print(f"  Min:  {na2_dists.min():.3f} Å")
print(f"  Max:  {na2_dists.max():.3f} Å")

print(f"\n{'→ Na1 is well-predicted' if na1_dists.mean() < 2.0 else '→ Na1 is NOT well-predicted'} "
      f"(mean displacement {na1_dists.mean():.2f} Å)")
print(f"{'→ Na2 is well-predicted' if na2_dists.mean() < 2.0 else '→ Na2 is NOT well-predicted'} "
      f"(mean displacement {na2_dists.mean():.2f} Å)")

# B-factor comparison
na1_bfacs = [r[3] for r in results]
na2_bfacs = [r[4] for r in results]
print(f"\nNa1 mean pLDDT/B-factor: {np.mean(na1_bfacs):.2f}")
print(f"Na2 mean pLDDT/B-factor: {np.mean(na2_bfacs):.2f}")

# Coordinate spread
na1_coords = np.array(na1_coords_all)
na2_coords = np.array(na2_coords_all)
print(f"\nNa1 coordinate spread (std): x={na1_coords[:,0].std():.2f}, y={na1_coords[:,1].std():.2f}, z={na1_coords[:,2].std():.2f}")
print(f"Na2 coordinate spread (std): x={na2_coords[:,0].std():.2f}, y={na2_coords[:,1].std():.2f}, z={na2_coords[:,2].std():.2f}")
