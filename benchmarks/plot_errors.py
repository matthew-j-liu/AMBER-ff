
from pathlib import Path
import h5py
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT   = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "results"
INPUT_DIR   = REPO_ROOT / "input_molecules_copy"
FIGURES_DIR = REPO_ROOT / "benchmarks" / "figures"

TERMS = {
    "bonds_energy":     "Bond",
    "angles_energy":    "Angle",
    "dihedrals_energy": "Dihedral",
    "vdw_energy":       "vdW",
    "total_energy":     "Total",
}

GROUP_COLORS = {
    "straight_chain_alkanes": "#1f77b4",
    "branched_alkanes":       "#17becf",
    "alkenes":                "#29e429",
    "alkynes":                "#46853A",
    "alcohols":               "#e11d24",
    "amines":                 "#a52880",
    "carboxylic_acids":       "#f28220",
    "halogenated":            "#7B4BC8",
    "amino_acids":            "#8B4513",
    "aldehydes":              "#FFB347",
    "2+_functional_groups":   "#eaf111",
}


def _load(group: str, mol: str, source: str) -> dict[str, float]:
    with h5py.File(RESULTS_DIR / group / mol / f"{source}.h5", "r") as f:
        return {k: float(f[k][()]) for k in TERMS}


def _atom_count(group: str, mol: str) -> int | None:
    for xyz in (INPUT_DIR / group).glob("*.xyz"):
        stem = xyz.stem
        for suffix in ("_PCS2", "_SE2", "_SE"):
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
        if stem == mol:
            try:
                return int(xyz.read_text().splitlines()[0].strip())
            except (ValueError, IndexError):
                return None
    return None


def _collect(groups: list[str]) -> dict:
    data: dict[str, list] = {}
    for group in groups:
        gdir = RESULTS_DIR / group
        if not gdir.is_dir():
            continue
        entries = []
        for mol_dir in sorted(gdir.iterdir()):
            if not mol_dir.is_dir():
                continue
            a_h5 = mol_dir / "amber.h5"
            c_h5 = mol_dir / "custom.h5"
            if not (a_h5.exists() and c_h5.exists()):
                continue
            amber  = _load(group, mol_dir.name, "amber")
            custom = _load(group, mol_dir.name, "custom")
            errors = {k: custom[k] - amber[k] for k in TERMS}
            n = _atom_count(group, mol_dir.name)
            entries.append((mol_dir.name, n, errors))
        if entries:
            data[group] = entries
    return data


def _save(fig, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return path


RMSE_GROUP_ORDER = [
    "straight_chain_alkanes",
    "branched_alkanes",
    "alkenes",
    "amines",
    "halogenated",
    "alcohols",
    "carboxylic_acids",
]

def _rmse_matrix(data: dict):
    groups = [g for g in RMSE_GROUP_ORDER if g in data]
    term_keys = ["bonds_energy", "angles_energy", "dihedrals_energy", "vdw_energy"]
    term_labels = [TERMS[t] for t in term_keys]
    rmse = np.zeros((len(groups), len(term_keys)))
    for i, g in enumerate(groups):
        for j, t in enumerate(term_keys):
            errs = np.array([e[t] for _, _, e in data[g]])
            if len(errs):
                rmse[i, j] = np.sqrt(np.mean(errs ** 2))
    return groups, term_keys, term_labels, rmse


def plot_rmse(data: dict) -> Path:
    groups, term_keys, term_labels, rmse = _rmse_matrix(data)

    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    axes = axes.flatten()
    x = np.arange(len(groups))

    for j, (key, label) in enumerate(zip(term_keys, term_labels)):
        ax = axes[j]
        colors = [GROUP_COLORS.get(g, "#333333") for g in groups]
        ax.bar(x, rmse[:, j], color=colors, edgecolor="black", linewidth=0.5)
        ax.set_title(f"{label} RMSE", fontsize=13, pad=6)
        ax.set_ylabel("RMSE (kcal/mol)", fontsize=11)
        ax.set_xticks(x)
        ax.set_xticklabels([g.replace("_", " ") for g in groups],
                           rotation=25, ha="right", fontsize=10)
        ax.grid(True, axis="y", linestyle="--", alpha=0.4)

    fig.suptitle("RMSE (Custom vs Benchmark)", fontsize=15, y=1.00)
    fig.tight_layout()
    return _save(fig, FIGURES_DIR / "error_rmse_small_multiples.png")


if __name__ == "__main__":
    groups = RMSE_GROUP_ORDER 
    data = _collect(groups)
    print(f"wrote {plot_rmse(data).relative_to(REPO_ROOT)}")
