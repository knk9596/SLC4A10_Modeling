# NBCn2 Structural Modeling and Simulation Analysis

Analysis scripts for MD simulations and AF3 modeling of the neuronal sodium-dependent bicarbonate transporter
**NBCn2 (SLC4A10)**.

These scripts support this study:

> **Structural insights enable drug discovery for the neuronal NBCn2
> carbonate transporter**\
> *Nature Communications* (2026)\
> DOI: 10.1038/s41467-026-75444-4

## Overview

This repository contains scripts used for post-processing and analysis
of **AlphaFold3 models** and **molecular dynamics simulations** of
NBCn2.

The analyses include:

-   alignment of AF3 models to the experimental reference structure;
-   comparison of predicted Na+ positions with experimentally observed
    ion sites;
-   calculation and visualization of AF3 ion-positioning results; and
-   analysis of distances between ions and key binding-site residues
    during MD simulations.

## Repository structure

``` text
.
├── AF3_ion_align_all.py
├── AF3_ion_analyze_na.py
├── AF3_ion_save_results.py
├── MD_distance_ion.py
└── README.md
```

## Scripts

### `AF3_ion_align_all.py`

Aligns AF3 models to the experimental reference structure using C-alpha
atoms and saves the aligned models as PDB files.

### `AF3_ion_analyze_na.py`

Checks predicted Na1/Na2 positions in the aligned AF3 models against the
experimental Na+ sites and reports per-model distances and summary
statistics.

### `AF3_ion_save_results.py`

Saves the AF3 ion-position analysis results to CSV and generates plots.

### `MD_distance_ion.py`

Analyzes an MD trajectory and plots distance versus time between the ion
and selected binding-site residues, including:

-   CO3-D831
-   SOD1-D831
-   SOD2-D831

A combined distance plot is also generated.

## Installation

A Python environment is recommended for running the analysis scripts.

For example:

``` bash
conda create -n nbcn2-analysis python=3.10
conda activate nbcn2-analysis
```

> **Note:** The exact Python dependencies should be determined from the
> imports in the analysis scripts. The associated study used MDAnalysis
> and ProLIF for MD trajectory analysis.

> Input paths, topology/trajectory filenames, and analysis parameters
> may need to be adjusted in the scripts for your local data
> organization.

## AlphaFold3 analysis

In the associated study, AlphaFold3 was used to predict NBCn2 under two
ion conditions:

1.  NBCn2 with **two sodium ions, one carbonate ion, and three water
    molecules**.
2.  NBCn2 with **two sodium ions only**.

A total of **1,100 models** were generated using **220 random seeds**,
with five models per seed.

Predicted structures were aligned to the experimentally determined
substrate-bound NBCn2 structure, and models with RMSD values greater
than 3 Å relative to the reference were excluded from downstream
analysis.

The analysis was used to evaluate how accurately AF3 reproduced the
experimentally observed ion-binding configuration, particularly the
Na+(1) and Na+(2) sites.

## Molecular dynamics analysis

The associated MD simulations were performed using **GROMACS** with the
**CHARMM36m** all-atom force field and **TIP3P** water.

Four independent replicas were generated using the same initial
coordinates with different initial velocities.

The combined simulation times reported in the study were:

-   **1.0 μs** for Cmpd38J-bound NBCn2
-   **1.2 μs** for substrate-bound NBCn2

MDAnalysis and ProLIF were used for trajectory and protein-ligand
interaction analysis in the study.

## Data availability

Structural coordinates and density maps are available through the PDB
and EMDB:

-   NBCn2 apo: **PDB 9MIX / EMD-48304**
-   NBCn2-NaHCO3: **PDB 9MJO / EMD-48318**
-   NBCn2-Cmpd38J: **PDB 9MK6 / EMD-48320**

MD simulation and modeling data associated with the study are available
through Zenodo (https://doi.org/10.5281/zenodo.20278709).

## Citation

If you use these scripts or the associated data, please cite:

**Yang S, Zhao Y, Vatansever S, et al. Structural insights enable drug
discovery for the neuronal NBCn2 carbonate transporter. Nature
Communications. 2026;17:8923.**

DOI: **10.1038/s41467-026-75444-4**

## License

Please refer to the repository license for terms of use.

