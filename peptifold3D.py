#!/usr/bin/env python3
"""
Peptide 3D Structure Pipeline
Reads peptides from a FASTA-like APD-format file and generates
energy-minimized PDB files named after each peptide's APD code.

Pipeline per peptide:
  1. Build extended-chain 3D backbone with PeptideBuilder
  2. Add missing heavy atoms & hydrogens with PDBFixer
  3. Energy-minimize with OpenMM (AMBER14 force field, implicit solvent)
  4. Strip solvent → save final PDB

Usage:
    python peptifold3D.py peptides_example.txt
    python peptifold3D.py peptides_example.txt -o my_pdbs/ --steps 2000
"""

import io
import os
import sys
import time
import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import requests
import PeptideBuilder
import Bio.PDB
import pdbfixer
import openmm as mm
from openmm import app, unit


SUPPORTED_AA = set("ACDEFGHIKLMNPQRSTVWY")
ESMFOLD_API = "https://api.esmatlas.com/foldSequence/v1/pdb/"


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------
def parse_apd_fasta(filepath: str) -> dict[str, str]:
    peptides: dict[str, str] = {}
    current_id = None
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_id = line[1:].strip()
            elif current_id:
                peptides[current_id] = line.upper()
                current_id = None
    return peptides


def validate_sequence(seq: str) -> tuple[bool, str]:
    bad = set(seq) - SUPPORTED_AA
    if bad:
        return False, f"Unsupported residues: {', '.join(sorted(bad))}"
    return True, ""


# ---------------------------------------------------------------------------
# Method 1 – ESMFold via ESM Atlas public API (no auth required)
# ---------------------------------------------------------------------------
def predict_esmfold(sequence: str, retries: int = 3) -> str:
    """Call ESM Atlas fold endpoint. Returns PDB string."""
    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    for attempt in range(1, retries + 1):
        resp = requests.post(ESMFOLD_API, headers=headers, data=sequence, timeout=120)
        if resp.status_code == 200 and resp.text.startswith("HEADER"):
            return resp.text
        if resp.status_code == 503 or resp.status_code == 429:
            wait = 15 * attempt
            print(f"\n      (sunucu meşgul, {wait}s bekleniyor...)", flush=True)
            time.sleep(wait)
            continue
        raise RuntimeError(f"ESMFold API {resp.status_code}: {resp.text[:200]}")
    raise RuntimeError("ESMFold API erişilemedi (retry sonrası)")


# ---------------------------------------------------------------------------
# Method 2 – build extended backbone (fallback)
# ---------------------------------------------------------------------------
def build_extended_pdb(sequence: str) -> str:
    """Return PDB string of extended-chain structure (backbone + side-chain atoms)."""
    structure = PeptideBuilder.make_extended_structure(sequence)
    PeptideBuilder.add_terminal_OXT(structure)
    buf = io.StringIO()
    io_obj = Bio.PDB.PDBIO()
    io_obj.set_structure(structure)
    io_obj.save(buf)
    return buf.getvalue()


# ---------------------------------------------------------------------------
# Step 2+3 – fix missing atoms, add H, minimize
# ---------------------------------------------------------------------------
def minimize(pdb_string: str, max_iterations: int) -> str:
    """
    Run PDBFixer + OpenMM energy minimization.
    Returns PDB string of minimized structure (no solvent).
    """
    buf = io.StringIO(pdb_string)
    fixer = pdbfixer.PDBFixer(pdbfile=buf)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    # AMBER14 + implicit solvent (OBC2) – no explicit water in output
    forcefield = app.ForceField("amber14-all.xml", "implicit/obc2.xml")

    modeller = app.Modeller(fixer.topology, fixer.positions)

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds,
    )

    integrator = mm.LangevinMiddleIntegrator(
        300 * unit.kelvin,
        1 / unit.picosecond,
        0.002 * unit.picoseconds,
    )

    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(maxIterations=max_iterations)

    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions()

    # Write to PDB string
    out = io.StringIO()
    app.PDBFile.writeFile(modeller.topology, positions, out)
    return out.getvalue()


# ---------------------------------------------------------------------------
# Per-peptide entry point
# ---------------------------------------------------------------------------
def generate_pdb(
    apd_code: str,
    sequence: str,
    output_dir: str,
    max_iterations: int,
) -> tuple[bool, str, str]:
    """Generate PDB. Returns (success, path_or_error, method_used)."""
    valid, msg = validate_sequence(sequence)
    if not valid:
        return False, msg, "-"

    output_path = os.path.join(output_dir, f"{apd_code}.pdb")

    # --- ESMFold (AI prediction via ESM Atlas) ---
    try:
        pdb_text = predict_esmfold(sequence)
        with open(output_path, "w") as f:
            f.write(pdb_text)
        return True, output_path, "ESMFold"
    except Exception as e:
        print(f"\n      [ESMFold hata: {e} → OpenMM fallback]", flush=True)

    # --- OpenMM minimization (fallback) ---
    try:
        pdb_raw = build_extended_pdb(sequence)
        pdb_min = minimize(pdb_raw, max_iterations)
        with open(output_path, "w") as f:
            f.write(pdb_min)
        return True, output_path, "OpenMM"
    except Exception as e:
        return False, str(e), "-"


# ---------------------------------------------------------------------------
# Quality analysis
# ---------------------------------------------------------------------------
# pLDDT confidence bands (ESMFold stores 0-1, AlphaFold uses 0-100;
# ESM Atlas returns 0-1 scale → multiply by 100 for standard display)
PLDDT_BANDS = [
    (90, 100, "#0053D6", "Very high (>90)"),
    (70,  90, "#65CBF3", "Confident (70–90)"),
    (50,  70, "#FFDB13", "Low (50–70)"),
    ( 0,  50, "#FF7D45", "Very low (<50)"),
]


def parse_plddt(pdb_path: str) -> tuple[list[int], list[float]]:
    """Extract per-residue pLDDT from CA B-factors. Returns (residue_ids, plddt_100)."""
    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_path)
    residues, scores = [], []
    for model in structure:
        for chain in model:
            for res in chain:
                for atom in res:
                    if atom.get_name() == "CA":
                        val = atom.get_bfactor()
                        # ESM Atlas: values in [0,1] → scale to [0,100]
                        if val <= 1.0:
                            val *= 100
                        residues.append(res.get_id()[1])
                        scores.append(round(val, 2))
                        break
    return residues, scores


def band_color(score: float) -> str:
    for lo, hi, color, _ in PLDDT_BANDS:
        if lo <= score <= hi:
            return color
    return "#FF7D45"


def quality_report(results: list[tuple[str, str, str]], output_dir: str) -> None:
    """
    Parse pLDDT from generated PDB files, print summary table,
    and save two figures:
      quality_plddt_perresidue.png  – per-residue pLDDT line plots
      quality_plddt_summary.png     – mean pLDDT bar chart
    """
    records = []  # (apd_code, method, residues, scores)
    for apd_code, pdb_path, method in results:
        if not os.path.isfile(pdb_path):
            continue
        res_ids, scores = parse_plddt(pdb_path)
        if scores:
            records.append((apd_code, method, res_ids, scores))

    if not records:
        print("  (no ESMFold PDB files found for quality analysis)")
        return

    # ── summary table ────────────────────────────────────────────────────────
    print("\n" + "=" * 62)
    print(f"  {'APD Code':<12} {'Method':<10} {'Length':>6} {'Mean pLDDT':>11}  Confidence")
    print("  " + "-" * 58)
    for apd_code, method, _, scores in records:
        mean_p = np.mean(scores)
        band = next(label for lo, hi, _, label in PLDDT_BANDS if lo <= mean_p <= hi)
        print(f"  {apd_code:<12} {method:<10} {len(scores):>6} {mean_p:>10.1f}  {band}")
    print("=" * 62)

    # ── figure 1: per-residue pLDDT – one file per peptide ───────────────────
    from matplotlib.patches import Patch
    legend_handles = [Patch(color=c, alpha=0.7, label=lbl) for _, _, c, lbl in PLDDT_BANDS]

    print()
    for apd_code, method, res_ids, scores in records:
        arr = np.array(scores)
        fig, ax = plt.subplots(figsize=(8, 3.5))

        for lo, hi, color, _ in reversed(PLDDT_BANDS):
            ax.axhspan(lo, hi, alpha=0.12, color=color, zorder=0)

        ax.plot(res_ids, arr, color="#444", lw=1.4, zorder=2)
        sc_colors = [band_color(v) for v in arr]
        ax.scatter(res_ids, arr, c=sc_colors, s=35, zorder=3, edgecolors="none")

        ax.set_xlim(res_ids[0] - 0.5, res_ids[-1] + 0.5)
        ax.set_ylim(0, 100)
        ax.set_title(f"{apd_code}  –  Per-residue pLDDT  ({method})",
                     fontsize=11, fontweight="bold")
        ax.set_xlabel("Residue index", fontsize=9)
        ax.set_ylabel("pLDDT", fontsize=9)
        ax.axhline(70, color="#65CBF3", lw=0.9, ls="--", zorder=1)
        ax.axhline(90, color="#0053D6", lw=0.9, ls="--", zorder=1)
        ax.legend(handles=legend_handles, fontsize=7, ncol=4,
                  loc="lower center", bbox_to_anchor=(0.5, -0.38),
                  title="pLDDT confidence", title_fontsize=7)
        plt.tight_layout()

        path1 = os.path.join(output_dir, f"quality_plddt_{apd_code}.png")
        fig.savefig(path1, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {path1}")

    # ── figure 2: mean pLDDT bar chart ───────────────────────────────────────
    labels = [r[0] for r in records]
    means  = [np.mean(r[3]) for r in records]
    colors = [band_color(m) for m in means]

    n = len(records)
    fig2, ax2 = plt.subplots(figsize=(max(7, n * 1.1), 4))
    bars = ax2.bar(labels, means, color=colors, edgecolor="white", linewidth=0.8)
    ax2.bar_label(bars, fmt="%.1f", padding=3, fontsize=8)

    for lo, hi, color, _ in PLDDT_BANDS:
        ax2.axhspan(lo, hi, alpha=0.08, color=color)
    ax2.axhline(70, color="#65CBF3", lw=1, ls="--", label="Confident (70)")
    ax2.axhline(90, color="#0053D6", lw=1, ls="--", label="Very high (90)")

    ax2.set_ylim(0, 105)
    ax2.set_ylabel("Mean pLDDT", fontsize=10)
    ax2.set_title("Mean pLDDT per Peptide (ESMFold)", fontsize=11)
    ax2.tick_params(axis="x", labelsize=8, rotation=30)
    ax2.legend(fontsize=8, loc="lower right")
    plt.tight_layout()
    path2 = os.path.join(output_dir, "quality_plddt_summary.png")
    fig2.savefig(path2, dpi=150, bbox_inches="tight")
    plt.close(fig2)
    print(f"  Saved: {path2}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Generate energy-minimized 3D PDB structures from APD peptide sequences."
    )
    parser.add_argument("input", help="Input FASTA file with APD codes and sequences")
    parser.add_argument(
        "-o", "--output", default="pdb_output",
        help="Output directory (default: pdb_output)",
    )
    parser.add_argument(
        "--steps", type=int, default=1000,
        help="Max OpenMM minimization iterations (fallback only, default: 1000)",
    )
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    method_note = "ESMFold (ESM Atlas API) → OpenMM AMBER14 fallback"
    print(f"Input       : {args.input}")
    print(f"Output dir  : {output_dir.resolve()}")
    print(f"Method      : {method_note}")
    print("-" * 60)

    peptides = parse_apd_fasta(args.input)
    if not peptides:
        print("ERROR: No peptides found in input file.")
        sys.exit(1)

    print(f"Found {len(peptides)} peptide(s)\n")

    ok_count = fail_count = 0
    completed: list[tuple[str, str, str]] = []  # (apd_code, pdb_path, method)

    for apd_code, sequence in peptides.items():
        label = f"{apd_code} ({len(sequence)} aa)"
        print(f"  {label} ... ", end="", flush=True)
        ok, result, method = generate_pdb(apd_code, sequence, str(output_dir), args.steps)
        if ok:
            print(f"OK [{method}]  →  {result}")
            ok_count += 1
            completed.append((apd_code, result, method))
        else:
            print(f"FAILED: {result}")
            fail_count += 1

    print("\n" + "-" * 60)
    print(f"Done: {ok_count} succeeded, {fail_count} failed")
    print(f"PDB files saved to: {output_dir.resolve()}")

    # Quality report (only for ESMFold outputs which have real pLDDT)
    esmfold_results = [(c, p, m) for c, p, m in completed if m == "ESMFold"]
    if esmfold_results:
        print("\nGenerating quality report...")
        quality_report(esmfold_results, str(output_dir))


if __name__ == "__main__":
    main()
