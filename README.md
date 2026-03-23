<img width="1536" height="1024" alt="peptifold_logo" src="https://github.com/user-attachments/assets/daa93cb6-9659-4da0-a3da-386f47ce1d1e" />

# PeptiFold3D

Automated pipeline that reads antimicrobial peptides from an APD-format FASTA file,
generates predicted 3D structures using **ESMFold**, and produces per-residue quality
(pLDDT) plots for each peptide.

---

## Overview

```
Input FASTA  →  ESMFold V1 (ESM Atlas API)  →  PDB files  →  pLDDT quality plots
                       ↓ (if API unavailable)
               OpenMM AMBER14 energy minimization  →  PDB files
```

**Primary method – ESMFold (ESM Atlas public API)**
Calls Meta's free [ESM Atlas](https://esmatlas.com) inference endpoint. No API key or
local GPU required. Returns full-atom structures with per-residue confidence scores
(pLDDT).

**Fallback method – OpenMM energy minimization**
Builds an extended backbone with PeptideBuilder, adds all heavy atoms and hydrogens
with PDBFixer, then minimises with the AMBER14 force field (implicit OBC2 solvent).
Used automatically when the ESM Atlas API is unreachable.

---

## Installation

```bash
conda env create -f environment.yml
conda activate peptifold3D
```

> **Note:** `openmm` is installed via conda-forge because it ships compiled platform
> plugins. All other packages are installed via pip inside the same environment.

---

## Input format

Standard FASTA with APD codes as sequence IDs:

```
>Pep_1
ACYCRIPACIAGERRYGTCIYQGRLWAFCC

>Pep_2
CYCRIPACIAGERRYGTCIYQGRLWAFCC
```

---

## Usage

```bash
# Basic – outputs to ./pdb_output/
python peptifold3D.py peptides_example.txt

# Custom output directory
python peptifold3D.py peptides_example.txt -o results/

# Increase OpenMM minimisation steps (fallback only)
python peptifold3D.py peptides_example.txt --steps 5000
```

### Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `input` | *(required)* | FASTA file with APD codes and sequences |
| `-o / --output` | `pdb_output` | Output directory for PDB files and plots |
| `--steps` | `1000` | Max OpenMM minimisation iterations (fallback only) |

---

## Output

For each peptide the pipeline writes:

| File | Description |
|------|-------------|
| `<APD_CODE>.pdb` | 3D structure (ESMFold or OpenMM) |
| `quality_plddt_<APD_CODE>.png` | Per-residue pLDDT plot |
| `quality_plddt_summary.png` | Mean pLDDT bar chart for all peptides |

### Example output (pLDDT per residue)

pLDDT (predicted Local Distance Difference Test) is ESMFold's per-residue confidence
score. Higher is better.

| Range | Colour | Interpretation |
|-------|--------|----------------|
| 90 – 100 | Blue | Very high confidence |
| 70 – 90 | Cyan | Confident |
| 50 – 70 | Yellow | Low confidence |
| 0 – 50 | Orange | Very low confidence |

---

## Requirements

| Package | Version | Purpose |
|---------|---------|---------|
| Python | 3.10 | |
| biopython | 1.79 | PDB parsing & I/O |
| PeptideBuilder | 1.1.0 | Extended-chain backbone generation |
| pdbfixer | 1.12.0 | Missing atom / hydrogen placement |
| openmm | 8.4.0 | AMBER14 energy minimisation |
| numpy | 1.26.0 | Numerical operations |
| matplotlib | 3.5.1 | Quality plots |
| requests | 2.32.3 | ESM Atlas API calls |

---

## Notes

- The ESM Atlas API accepts sequences up to ~400 residues.
- pLDDT scores are only available for **ESMFold** outputs. OpenMM-minimised
  structures do not contain meaningful B-factor quality scores.
- For large batches, be mindful of ESM Atlas rate limits; the pipeline retries
  automatically on `503` / `429` responses.
