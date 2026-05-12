"""
Run main.cpp and produce per-group parity plots.

Usage (from repo root):
    python benchmarks/run_parity.py                  # all groups
    python benchmarks/run_parity.py straight_chain_alkanes alcohols
"""

import re
import sys
import shutil
import subprocess
from pathlib import Path
import h5py

sys.path.insert(0, str(Path(__file__).parent))

REPO_ROOT = Path(__file__).resolve().parent.parent
INPUT_DIR = REPO_ROOT / "input_molecules_copy"
PARAMS    = REPO_ROOT / "params" / "selected_atoms.dat"
BUILD_BIN = shutil.which("amber_ff") or str(REPO_ROOT / "build" / "amber_ff")

def _mol_name(stem: str) -> str:
    """Strip source-file artifact suffixes so names match existing results dirs."""
    for suffix in ("_PCS2", "_SE2", "_SE"):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return stem


def discover_molecule_groups() -> dict[str, list[tuple[str, str]]]:
    """Scan INPUT_DIR for subdirectories and collect all .xyz files within each."""
    groups: dict[str, list[tuple[str, str]]] = {}
    for group_dir in sorted(INPUT_DIR.iterdir()):
        if not group_dir.is_dir():
            continue
        molecules = [
            (_mol_name(xyz.stem), f"{group_dir.name}/{xyz.name}")
            for xyz in sorted(group_dir.glob("*.xyz"))
        ]
        if molecules:
            groups[group_dir.name] = molecules
    return groups


# Molecules organised by functional group, auto-discovered from INPUT_DIR.
# Each entry: (display_name, path_relative_to_INPUT_DIR)
MOLECULE_GROUPS: dict[str, list[tuple[str, str]]] = discover_molecule_groups()


def parse_amber_energies(sp_out: Path) -> dict[str, float]:
    """Extract bond/angle/dihed/vdwaals/1-4 VDW from a sander .out file."""
    text = sp_out.read_text()

    def grab(label):
        m = re.search(rf"{label}\s*=\s*(-?\d+\.\d+)", text)
        if m is None:
            raise RuntimeError(f"missing '{label}' in {sp_out}")
        return float(m.group(1))

    return {
        "bond":    grab("BOND"),
        "angle":   grab("ANGLE"),
        "dihed":   grab("DIHED"),
        "vdwaals": grab("VDWAALS"),
        "nb14":    grab("1-4 VDW"),
    }


def parse_custom_output(stdout: str) -> dict[str, float]:
    """Parse energy lines from the custom C++ binary stdout."""
    def grab(label):
        m = re.search(
            rf"{re.escape(label)}\s*=\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)",
            stdout,
        )
        if m is None:
            raise RuntimeError(f"missing '{label}' in custom binary output")
        return float(m.group(1))

    return {
        "bonds_energy":     grab("Bond stretching energy"),
        "angles_energy":    grab("Bond Rotation energy"),
        "dihedrals_energy": grab("Torsion energy"),
        "vdw_energy":       grab("Van der waals energy"),
        "total_energy":     grab("Total energy"),
    }


def write_amber_h5(path: Path, terms: dict[str, float]) -> None:
    vdw_total = terms["vdwaals"] + terms["nb14"]
    total     = terms["bond"] + terms["angle"] + terms["dihed"] + vdw_total
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f.create_dataset("bonds_energy",     data=terms["bond"])
        f.create_dataset("angles_energy",    data=terms["angle"])
        f.create_dataset("dihedrals_energy", data=terms["dihed"])
        f.create_dataset("vdw_energy",       data=vdw_total)
        f.create_dataset("total_energy",     data=total)


def write_custom_h5(path: Path, energies: dict[str, float]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        for k, v in energies.items():
            f.create_dataset(k, data=v)


def run_custom(mol: str, xyz_path: Path) -> dict[str, float]:
    print(f"  [custom] {mol}")
    result = subprocess.run(
        [str(BUILD_BIN), str(xyz_path), str(PARAMS)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(f"custom run failed for {mol}")
    return parse_custom_output(result.stdout)


def main():
    # Allow filtering to specific groups via CLI args
    requested = set(sys.argv[1:]) if len(sys.argv) > 1 else set(MOLECULE_GROUPS)
    unknown = requested - set(MOLECULE_GROUPS)
    if unknown:
        sys.exit(f"Unknown groups: {unknown}\nAvailable: {list(MOLECULE_GROUPS)}")

    ran_any = False
    for group, molecules in MOLECULE_GROUPS.items():
        if group not in requested:
            continue
        print(f"\n=== {group} ===")
        for mol, rel_path in molecules:
            xyz_path = INPUT_DIR / rel_path
            if not xyz_path.exists():
                print(f"  [skip] {mol} — xyz not found: {xyz_path}")
                continue

            # sp.out is under amber_runs/<group>/<stem>/ — matches run_AMBER_benchmarking.sh
            safe_rel = rel_path.replace(" ", "_")
            sp_out    = REPO_ROOT / "amber_runs" / Path(safe_rel).with_suffix("") / "sp.out"
            amber_h5  = REPO_ROOT / "results" / group / mol / "amber.h5"
            custom_h5 = REPO_ROOT / "results" / group / mol / "custom.h5"

            if not sp_out.exists():
                print(f"  [skip] {mol} — sp.out not found: {sp_out}")
                continue

            custom_energies = run_custom(mol, xyz_path)
            write_amber_h5(amber_h5, parse_amber_energies(sp_out))
            write_custom_h5(custom_h5, custom_energies)
            print(f"    wrote {amber_h5.relative_to(REPO_ROOT)}")
            print(f"    wrote {custom_h5.relative_to(REPO_ROOT)}")
            ran_any = True

    if not ran_any:
        print("Nothing to run — check group names or xyz files.")
        return

    print("\nplotting parity figures...")
    from plot_parity import plot_all
    for q in plot_all(MOLECULE_GROUPS):
        print(f"  wrote {q.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
