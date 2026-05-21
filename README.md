MD simulations: plot ion distance to key binding-site residue(s) over simulation time.
AF3 modeling: align AF3 models to the experimental structure and compare predicted vs experimental Na binding sites.

AF3_ion_align_all.py
Aligns all AF3 models to the experimental reference structure using CA atoms, then saves the aligned models as PDBs.

AF3_ion_analyze_na.py
Checks the predicted Na1/Na2 positions in the aligned AF3 models against the experimental Na sites, then prints per-model distances and summary stats.

AF3_ion_save_results.py
Saves the results to a CSV and generates plots.

MD_distance_ion.py
Analyzes the MD trajectory and plots distance vs time between the ion and key binding-site residue(s): CO3–D831, SOD1–D831, and SOD2–D831, plus a combined plot.
