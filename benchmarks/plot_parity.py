"""
Make parity plots for each energy term 
"""
from pathlib import Path
from typing import Iterable
import h5py
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "results"
FIGURES_DIR = REPO_ROOT / "benchmarks" / "figures"

TERMS = {
    "bonds_energy": "Bond stretching",
    "angles_energy": "Bond angle bending",
    "dihedrals_energy": "Dihedral torsion",
    "vdw_energy": "Van der Waals",
    "total_energy": "Total",
}


def load(mol: str, source: str) -> dict[str, float]:
    with h5py.File(RESULTS_DIR / mol / f"{source}.h5", "r") as f:
        return {k: float(f[k][()]) for k in TERMS}


def plot_all(molecules: Iterable[str]) -> list[Path]:
    molecules = list(molecules)
    custom = {m: load(m, "custom") for m in molecules}
    amber = {m: load(m, "amber") for m in molecules}

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    palette = plt.cm.viridis(np.linspace(0.1, 0.85, len(molecules)))

    saved: list[Path] = []
    for term_key, term_label in TERMS.items():
        x = np.array([custom[m][term_key] for m in molecules])
        y = np.array([amber[m][term_key] for m in molecules])

        fig, ax = plt.subplots(figsize=(7.5, 7.5))

        for i, mol in enumerate(molecules):
            ax.scatter(x[i], y[i], s=180, color=palette[i],
                       edgecolor="black", linewidth=1.0, label=mol, zorder=3)

        lo = min(x.min(), y.min())
        hi = max(x.max(), y.max())
        pad = 0.05 * max(hi - lo, 1e-6)
        line = np.array([lo - pad, hi + pad])
        ax.plot(line, line, "k:", linewidth=2.0, label="y = x", zorder=2)
        ax.set_xlim(line[0], line[1])
        ax.set_ylim(line[0], line[1])

        ax.set_xlabel("Custom implementation (kcal/mol)", fontsize=18)
        ax.set_ylabel("AMBER reference (kcal/mol)", fontsize=18)
        ax.set_title(f"{term_label} energy parity", fontsize=20, pad=14)
        ax.tick_params(axis="both", which="major", labelsize=14)
        ax.grid(True, linestyle="--", alpha=0.4, zorder=1)
        ax.set_aspect("equal", adjustable="box")
        ax.legend(fontsize=14, loc="best", frameon=True, framealpha=0.9)

        out = FIGURES_DIR / f"parity_{term_key.replace('_energy', '')}.png"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(fig)
        saved.append(out)

    return saved


if __name__ == "__main__":
    import sys
    mols = sys.argv[1:] or sorted(p.name for p in RESULTS_DIR.iterdir() if p.is_dir())
    if not mols:
        sys.exit("no molecules found in results/ -- run run_parity.py first")
    for q in plot_all(mols):
        print(f"wrote {q.relative_to(REPO_ROOT)}")
