"""
Parity plot colored by broad chemical class (shareable summary figure).

Classes shown: alkanes, alkenes, alkynes, carboxylic acids, halogenated.

Outputs:
  benchmarks/figures/by_class/parity_<term>.png  — one per energy term

Run from repo root:
  python benchmarks/plot_by_class.py
"""

import sys
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))

REPO_ROOT   = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "results"
FIGURES_DIR = REPO_ROOT / "benchmarks" / "figures" / "by_class"

TERMS = {
    "bonds_energy":     "Bond stretching",
    "angles_energy":    "Bond angle bending",
    "dihedrals_energy": "Dihedral torsion",
    "vdw_energy":       "Van der Waals",
    "total_energy":     "Total",
}

# Map functional-group folder names → broad chemical class label
CLASS_MAP: dict[str, str] = {
    "straight_chain_alkanes": "Straight-chain alkanes",
    "branched_alkanes":       "Branched alkanes",
    "alkenes":                "Alkenes",
    "alkynes":                "Alkynes",
    "alcohols":               "Alcohols",
    "amines":                 "Amines",
    "carboxylic_acids":       "Carboxylic acids",
    "halogenated":            "Haloalkanes",
    "2+_functional_groups":   "2+ functional groups",
}

CLASS_COLORS: dict[str, str] = {
    "Straight-chain alkanes": "#1f77b4",
    "Branched alkanes":       "#17becf",
    "Alkenes":                "#29e429",
    "Alkynes":                "#46853A",
    "Alcohols":               "#e11d24",
    "Amines":                 "#a52880",
    "Carboxylic acids":       "#f28220",
    "Haloalkanes":            "#7B4BC8",
    "2+ functional groups":   "#eaf111",
}

CLASS_MARKERS: dict[str, str] = {cls: "o" for cls in CLASS_COLORS}


def _load(group: str, mol: str, source: str) -> dict[str, float]:
    path = RESULTS_DIR / group / mol / f"{source}.h5"
    with h5py.File(path, "r") as f:
        return {k: float(f[k][()]) for k in TERMS}


def _collect() -> dict[str, list[tuple[str, str, dict, dict]]]:
    """Return {class_label: [(group, mol, amber_vals, custom_vals), ...]}."""
    from run_parity import MOLECULE_GROUPS

    data: dict[str, list] = {cls: [] for cls in CLASS_COLORS}
    missing: list[str] = []
    unmapped: list[str] = []
    for group, molecules in MOLECULE_GROUPS.items():
        cls = CLASS_MAP.get(group)
        if cls is None:
            unmapped.append(group)
            continue
        for mol, _ in molecules:
            a_h5 = RESULTS_DIR / group / mol / "amber.h5"
            c_h5 = RESULTS_DIR / group / mol / "custom.h5"
            if not (a_h5.exists() and c_h5.exists()):
                missing.append(f"{group}/{mol}")
                continue
            data[cls].append((group, mol, _load(group, mol, "amber"),
                               _load(group, mol, "custom")))

    if unmapped:
        print(f"[warn] groups in input_molecules_copy with no class mapping: {unmapped}")
    if missing:
        print(f"[warn] {len(missing)} molecules in input_molecules_copy missing h5 results "
              f"(run run_parity.py to generate):")
        for m in missing:
            print(f"   - {m}")
    return data


def _parity_line(ax, x_all: list, y_all: list, symlog: bool = False) -> None:
    if symlog:
        vals = list(x_all) + list(y_all)
        if not vals:
            return
        lo, hi = min(vals), max(vals)
        span = max(abs(lo), abs(hi)) * 1.5
        line = np.array([-span, span]) if lo < 0 else np.array([lo / 1.5, hi * 1.5])
    else:
        pos = [v for v in (*x_all, *y_all) if v > 0]
        if not pos:
            return
        lo = min(pos)
        hi = max(pos)
        line = np.array([lo / 1.5, hi * 1.5])
    ax.plot(line, line, "k:", linewidth=1.8, label="y = x", zorder=2)
    ax.set_xlim(line[0], line[1])
    ax.set_ylim(line[0], line[1])
    ax.set_aspect("equal", adjustable="box")


def plot_by_class() -> list[Path]:
    data = _collect()

    total_mols = sum(len(v) for v in data.values())
    if total_mols == 0:
        sys.exit("No H5 data found — run run_parity.py first.")

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    saved: list[Path] = []

    for term_key, term_label in TERMS.items():
        fig, ax = plt.subplots(figsize=(8, 8))
        all_x, all_y = [], []
        legend_added: set[str] = set()

        for cls in CLASS_COLORS:
            entries = data[cls]
            if not entries:
                continue
            color  = CLASS_COLORS[cls]
            marker = CLASS_MARKERS[cls]
            for group, mol, amber, custom in entries:
                xi, yi = custom[term_key], amber[term_key]
                label = cls if cls not in legend_added else ""
                ax.scatter(
                    xi, yi,
                    s=120, color=color, marker=marker,
                    edgecolor="black", linewidth=0.7,
                    alpha=0.7,
                    label=label, zorder=3,
                )
                legend_added.add(cls)
                all_x.append(xi)
                all_y.append(yi)

        if not all_x:
            plt.close(fig)
            continue

        _parity_line(ax, all_x, all_y, symlog=(term_key == "vdw_energy"))
        if term_key == "vdw_energy":
            ax.set_xscale("symlog", linthresh=1e-2)
            ax.set_yscale("symlog", linthresh=1e-2)
        else:
            ax.set_xscale("log")
            ax.set_yscale("log")
        if term_key == "dihedrals_energy":
            xlo, xhi = ax.get_xlim()
            ylo, yhi = ax.get_ylim()
            ax.set_xlim(max(xlo, 1e-2), xhi)
            ax.set_ylim(max(ylo, 1e-2), yhi)
        ax.set_xlabel("Custom implementation (kcal/mol)", fontsize=14)
        ax.set_ylabel("AMBER reference (kcal/mol)", fontsize=14)
        titles = {
            "total_energy":     "Total Energy: by functional group",
            "vdw_energy":       "Van der Waals forces: by functional groups",
            "angles_energy":    "Bond Rotation Energy: by functional groups",
            "dihedrals_energy": "Torsion Energy: by functional groups",
            "bonds_energy":     "Bond Stretching Energy: by functional groups",
        }
        title = titles.get(term_key, f"{term_label} — parity by chemical class")
        ax.set_title(title, fontsize=14, pad=12)
        ax.tick_params(axis="both", labelsize=12)
        ax.grid(True, linestyle="--", alpha=0.35, zorder=1)

        handles, labels = ax.get_legend_handles_labels()
        non_empty = [(h, l) for h, l in zip(handles, labels) if l]
        if non_empty:
            ax.legend(*zip(*non_empty), fontsize=12, loc="best",
                      frameon=True, framealpha=0.9)

        if term_key == "total_energy":
            out = FIGURES_DIR / "Total Energy : by functional groups.png"
        elif term_key == "vdw_energy":
            out = FIGURES_DIR / "vdw_symlog.png"
        else:
            out = FIGURES_DIR / f"parity_{term_key.replace('_energy', '')}.png"
        fig.savefig(out, dpi=300, bbox_inches="tight")
        plt.close(fig)
        saved.append(out)
        print(f"wrote {out.relative_to(REPO_ROOT)}")

    return saved


if __name__ == "__main__":
    plot_by_class()
