#!/usr/bin/env python3
"""Align all AF3 seed/sample CIF models to a reference PDB using CA atoms."""

import glob
import os
import warnings
from Bio.PDB import PDBParser, Superimposer, PDBIO, Select
from Bio.PDB.MMCIFParser import MMCIFParser

warnings.filterwarnings("ignore")

REF_PDB = "/sc/arion/scratch/zhaoy14/4A10_revision/0316/exp_chainB.pdb"
AF3_DIR = "/sc/arion/scratch/zhaoy14/ALineMol/AF3_ion/output/4a10"
OUT_DIR = "/sc/arion/scratch/zhaoy14/ALineMol/AF3_ion/aligned"

os.makedirs(OUT_DIR, exist_ok=True)

# Parse reference
pdb_parser = PDBParser(QUIET=True)
ref_struct = pdb_parser.get_structure("ref", REF_PDB)
ref_model = ref_struct[0]

# Get reference CA atoms (first chain)
ref_chain = list(ref_model.get_chains())[0]
ref_cas = []
ref_resids = set()
for res in ref_chain:
    if res.id[0] == " " and "CA" in res:
        ref_cas.append(res["CA"])
        ref_resids.add(res.id[1])

print(f"Reference: {REF_PDB}")
print(f"Reference chain {ref_chain.id}: {len(ref_cas)} CA atoms")
print(f"Output dir: {OUT_DIR}\n")

# Find all seed folders
cif_parser = MMCIFParser(QUIET=True)
seed_dirs = sorted(glob.glob(os.path.join(AF3_DIR, "seed-*_sample-*")))

for seed_dir in seed_dirs:
    cif_files = glob.glob(os.path.join(seed_dir, "*.cif"))
    if not cif_files:
        continue

    cif_path = cif_files[0]
    name = os.path.basename(seed_dir)
    print(f"Processing {name}...", end=" ")

    # Parse the model CIF
    mov_struct = cif_parser.get_structure(name, cif_path)
    mov_model = mov_struct[0]

    # Get first chain (protein) CA atoms matching reference residues
    mov_chain = list(mov_model.get_chains())[0]
    mov_cas = []
    matched_ref_cas = []

    for res in mov_chain:
        if res.id[0] == " " and "CA" in res and res.id[1] in ref_resids:
            mov_cas.append(res["CA"])

    # Match by sequential order (both should be same protein)
    # Use minimum overlapping length
    n = min(len(ref_cas), len(mov_cas))
    if n < 10:
        print(f"SKIPPED (only {n} matching CA atoms)")
        continue

    sup = Superimposer()
    sup.set_atoms(ref_cas[:n], mov_cas[:n])
    sup.apply(mov_struct.get_atoms())

    rmsd = sup.rms
    print(f"RMSD = {rmsd:.3f} A ({n} CA atoms)")

    # Save aligned structure as PDB
    out_path = os.path.join(OUT_DIR, f"{name}_aligned.pdb")
    io = PDBIO()
    io.set_structure(mov_struct)
    io.save(out_path)

print(f"\nDone! Aligned structures saved to {OUT_DIR}")
