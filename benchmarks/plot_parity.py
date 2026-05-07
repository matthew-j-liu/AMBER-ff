"""
Generate parity plots grouped by functional group.

Outputs:
  benchmarks/figures/<group>/parity_<term>.png  — one plot per term per group
  benchmarks/figures/parity_<term>_all.png      — combined overview (all groups)

Run standalone (after run_parity.py has populated results/):
  python benchmarks/plot_parity.py
  python benchmarks/plot_parity.py straight_chain_alkanes alcohols
"""

from pathlib import Path
from typing import Iterable
import h5py
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT   = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "results"
FIGURES_DIR = REPO_ROOT / "benchmarks" / "figures"

TERMS = {
    "bonds_energy":     "Bond stretching",
    "angles_energy":    "Bond angle bending",
    "dihedrals_energy": "Dihedral torsion",
    "vdw_energy":       "Van der Waals",
    "total_energy":     "Total",
}

GROUP_COLORS = {
    "straight_chain_alkanes": "#1f77b4",
    "branched_alkanes":       "#ff7f0e",
    "alkenes":                "#2ca02c",
    "alkynes":                "#d62728",
    "alcohols":               "#9467bd",
    "amines":                 "#8c564b",
    "aldehydes":              "#e377c2",
    "carboxylic_acids":       "#7f7f7f",
    "halogenated":            "#bcbd22",
    "amino_acids":            "#17becf",
}


def load(group: str, mol: str, source: str) -> dict[str, float]:
    with h5py.File(RESULTS_DIR / group / mol / f"{source}.h5", "r") as f:
        return {k: float(f[k][()]) for k in TERMS}


def _collect(
    molecule_groups: dict[str, list[tuple[str, str]]],
) -> dict[str, dict[str, dict]]:
    """Return {group: {mol: {amber: {...}, custom: {...}}}} for existing H5 pairs."""
    data: dict[str, dict[str, dict]] = {}
    for group, molecules in molecule_groups.items():
        data[group] = {}
        for mol, _ in molecules:
            a_h5 = RESULTS_DIR / group / mol / "amber.h5"
            c_h5 = RESULTS_DIR / group / mol / "custom.h5"
            if a_h5.exists() and c_h5.exists():
                data[group][mol] = {
                    "amber":  load(group, mol, "amber"),
                    "custom": load(group, mol, "custom"),
                }
    return data


def _parity_axes(ax, x_all, y_all):
    """Draw y=x reference line sized to the data."""
    lo = min(min(x_all), min(y_all))
    hi = max(max(x_all), max(y_all))
    pad = 0.05 * max(hi - lo, 1e-6)
    line = np.array([lo - pad, hi + pad])
    ax.plot(line, line, "k:", linewidth=2.0, label="y = x", zorder=2)
    ax.set_xlim(line[0], line[1])
    ax.set_ylim(line[0], line[1])
    ax.set_aspect("equal", adjustable="box")


def _save(fig, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_per_group(
    molecule_groups: dict[str, list[tuple[str, str]]],
    all_data: dict,
) -> list[Path]:
    """One plot per energy term per group."""
    saved: list[Path] = []
    for group, mol_data in all_data.items():
        if not mol_data:
            continue
        mols = list(mol_data)
        palette = plt.cm.viridis(np.linspace(0.1, 0.85, len(mols)))
        group_dir = FIGURES_DIR / group

        for term_key, term_label in TERMS.items():
            x = [mol_data[m]["custom"][term_key] for m in mols]
            y = [mol_data[m]["amber"][term_key]  for m in mols]

            fig, ax = plt.subplots(figsize=(7, 7))
            for i, mol in enumerate(mols):
                ax.scatter(x[i], y[i], s=160, color=palette[i],
                           edgecolor="black", linewidth=0.8, label=mol, zorder=3)

            _parity_axes(ax, x, y)
            ax.set_xlabel("Custom implementation (kcal/mol)", fontsize=15)
            ax.set_ylabel("AMBER reference (kcal/mol)", fontsize=15)
            ax.set_title(
                f"{term_label} — {group.replace('_', ' ')}",
                fontsize=15, pad=12,
            )
            ax.tick_params(axis="both", labelsize=12)
            ax.grid(True, linestyle="--", alpha=0.4, zorder=1)
            ax.legend(fontsize=11, loc="best", frameon=True, framealpha=0.9)

            out = group_dir / f"parity_{term_key.replace('_energy', '')}.png"
            saved.append(_save(fig, out))

    return saved


def plot_combined(
    molecule_groups: dict[str, list[tuple[str, str]]],
    all_data: dict,
) -> list[Path]:
    """One overview plot per energy term, all groups coloured differently."""
    groups_present = [g for g in molecule_groups if all_data.get(g)]
    if not groups_present:
        return []

    saved: list[Path] = []
    for term_key, term_label in TERMS.items():
        fig, ax = plt.subplots(figsize=(9, 9))
        all_x, all_y = [], []

        for group in groups_present:
            mol_data = all_data[group]
            if not mol_data:
                continue
            color = GROUP_COLORS.get(group, "#333333")
            first = True
            for mol, vals in mol_data.items():
                xi = vals["custom"][term_key]
                yi = vals["amber"][term_key]
                ax.scatter(
                    xi, yi, s=160, color=color,
                    edgecolor="black", linewidth=0.8,
                    label=group.replace("_", " ") if first else "",
                    zorder=3,
                )
                all_x.append(xi)
                all_y.append(yi)
                first = False

        if not all_x:
            plt.close(fig)
            continue

        _parity_axes(ax, all_x, all_y)
        ax.set_xlabel("Custom implementation (kcal/mol)", fontsize=16)
        ax.set_ylabel("AMBER reference (kcal/mol)", fontsize=16)
        ax.set_title(f"{term_label} energy parity — all groups", fontsize=16, pad=12)
        ax.tick_params(axis="both", labelsize=13)
        ax.grid(True, linestyle="--", alpha=0.4, zorder=1)

        handles, labels = ax.get_legend_handles_labels()
        non_empty = [(h, l) for h, l in zip(handles, labels) if l]
        if non_empty:
            ax.legend(*zip(*non_empty), fontsize=12, loc="best",
                      frameon=True, framealpha=0.9)

        out = FIGURES_DIR / f"parity_{term_key.replace('_energy', '')}_all.png"
        saved.append(_save(fig, out))

    return saved


def plot_all(
    molecule_groups: dict[str, list[tuple[str, str]]],
) -> list[Path]:
    all_data = _collect(molecule_groups)
    saved  = plot_per_group(molecule_groups, all_data)
    saved += plot_combined(molecule_groups, all_data)
    return saved


if __name__ == "__main__":
    import sys
    from run_parity import MOLECULE_GROUPS

    requested_groups = sys.argv[1:]
    if requested_groups:
        groups = {g: MOLECULE_GROUPS[g] for g in requested_groups if g in MOLECULE_GROUPS}
    else:
        groups = MOLECULE_GROUPS

    if not groups:
        sys.exit("No valid groups specified.")

    for q in plot_all(groups):
        print(f"wrote {q.relative_to(REPO_ROOT)}")
