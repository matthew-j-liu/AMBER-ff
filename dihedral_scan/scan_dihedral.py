"""
Runs the C++ `dihedral_scan` binary file (which does the rotation + energy eval
using src/ routines) and plots each energy term and the total vs. dihedral angle.
"""

import io
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent
BUILD_BIN = REPO_ROOT / "build" / "dihedral_scan"
OUT_DIR   = REPO_ROOT / "dihedral_scan" / "results"

def main():
    if not BUILD_BIN.exists():
        sys.exit(f"binary not found: {BUILD_BIN} - build the project first")

    res = subprocess.run([str(BUILD_BIN)], cwd=REPO_ROOT,
                         capture_output=True, text=True)
    if res.returncode != 0:
        sys.stderr.write(res.stderr)
        sys.exit("dihedral_scan binary failed")

    data = np.genfromtxt(io.StringIO(res.stdout), delimiter=",", names=True)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_DIR / f"But-1-ene_scan.csv"
    with csv_path.open("w") as f:
        f.write("angle_deg,bonds,angles,dihedrals,vdw,total\n")
        for r in data:
            f.write(f"{r['angle_deg']:.4f},{r['bonds']:.6f},{r['angles']:.6f},"
                    f"{r['dihedrals']:.6f},{r['vdw']:.6f},{r['total']:.6f}\n")
    print(f"wrote {csv_path.relative_to(REPO_ROOT)}")

    i, j, k, l = 3, 2, 1, 0

    fig, ax = plt.subplots(figsize=(7, 4.5))

    plot_labels = {
    "bonds":     "Bond stretching energy",
    "angles":    "Bond Rotation energy",
    "dihedrals": "Torsion energy",
    "vdw":       "Van der waals energy",
    "total":     "Total energy"}

    for key, label in plot_labels.items():
        ax.plot(data["angle_deg"], data[key], label=label, lw=1.5, alpha=0.6 if key != "total" else 1.0)
    ax.set_xlabel(f"dihedral {i}-{j}-{k}-{l} (deg)")
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_ylabel("energy (kcal/mol)")
    ax.set_title(f"But-1-ene torsion scan")
    ax.axhline(0, color="k", lw=0.5, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig_path = OUT_DIR / f"But-1-ene_scan.png"
    fig.savefig(fig_path, dpi=150)
    print(f"wrote {fig_path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
