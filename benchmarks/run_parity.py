"""
Run main.cpp and produce per-group parity plots.

Usage (from repo root):
    python benchmarks/run_parity.py                  # all groups
    python benchmarks/run_parity.py straight_chain_alkanes alcohols
"""

import os
import re
import sys
import subprocess
from pathlib import Path
import h5py

sys.path.insert(0, str(Path(__file__).parent))

REPO_ROOT = Path(__file__).resolve().parent.parent
INPUT_DIR = REPO_ROOT / "input_molecules_copy"
PARAMS    = REPO_ROOT / "params" / "selected_atoms.dat"
BUILD_BIN = REPO_ROOT / "build" / "amber_ff"

# Molecules organised by functional group.
# Each entry: (display_name, path_relative_to_INPUT_DIR)
MOLECULE_GROUPS: dict[str, list[tuple[str, str]]] = {
    "straight_chain_alkanes": [
        ("methane",   "straight_chain_alkanes/methane.xyz"),
        ("ethane",    "straight_chain_alkanes/ethane.xyz"),
        ("propane",   "straight_chain_alkanes/propane.xyz"),
        ("n-butane",  "straight_chain_alkanes/n-butane.xyz"),
        ("n-pentane", "straight_chain_alkanes/n-pentane.xyz"),
        ("n-hexane",  "straight_chain_alkanes/n-hexane.xyz"),
    ],
    "branched_alkanes": [
        ("isobutane",          "branched_alkanes/isobutane.xyz"),
        ("isopentane",         "branched_alkanes/isopentane.xyz"),
        ("neopentane",         "branched_alkanes/neopentane.xyz"),
        ("2-methylpropane",    "branched_alkanes/2-Methylpropane.xyz"),
        ("2_3-dimethylbutane", "branched_alkanes/2_3-dimethylbutane.xyz"),
    ],
    "alkenes": [
        ("ethene",             "alkenes/ethene.xyz"),
        ("propene",            "alkenes/propene.xyz"),
        ("isobutylene",        "alkenes/isobutylene.xyz"),
        ("2-methyl-1-propene", "alkenes/2-methyl-1-propene_PCS2.xyz"),
        ("but-1-ene",          "alkenes/But-1-ene.xyz"),
    ],
    "alkynes": [
        ("ethyne",                  "alkynes/ethyne.xyz"),
        ("propyne",                 "alkynes/propyne.xyz"),
        ("1-butyne",                "alkynes/1-butyne.xyz"),
        ("3-methyl-but-3-en-1-yne", "alkynes/3-methyl-but-3-en-1-yne_PCS2.xyz"),
    ],
    "alcohols": [
        ("methanol",        "alcohols/methanol.xyz"),
        ("ethanol",         "alcohols/ethanol.xyz"),
        ("isopropanol",     "alcohols/isopropanol.xyz"),
        ("1-propanol",      "alcohols/1-propanol.xyz"),
        ("ethylene_glycol", "alcohols/ethylene_glycol.xyz"),
        ("glycerol",        "alcohols/glycerol.xyz"),
    ],
    "aldehydes": [
        ("propanal", "aldehydes/propanal.xyz"),
    ],
    "amines": [
        ("methylamine", "amines/methylamine.xyz"),
        ("ethylamine",  "amines/ethylamine.xyz"),
        ("azepane",     "amines/azepane.xyz"),
    ],
    "carboxylic_acids": [
        ("acetic_acid",        "carboxylic_acids/acetic_acid.xyz"),
        ("glycolic_acid",      "carboxylic_acids/glycolic_acid.xyz"),
        ("cis_acrylic_acid",   "carboxylic_acids/cis_acrylic_acid.xyz"),
        ("trans_acrylic_acid", "carboxylic_acids/trans_acrylic_acid.xyz"),
    ],
    "halogenated": [
        ("methylene_chloride", "halogenated/methylene_chloride.xyz"),
        ("ethyl_chloride",     "halogenated/ethyl_chloride.xyz"),
    ],
    "amino_acids": [
        ("glycine",  "amino_acids/glycine.xyz"),
        ("alanine",  "amino_acids/alanine.xyz"),
    ],
}


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


def run_amber(mol: str, xyz_path: Path) -> None:
    print(f"  [AMBER]  {mol}")
    script = REPO_ROOT / "run_AMBER_benchmarking_simple.sh"
    ret = os.system(f"'{script}' '{xyz_path}'")
    if ret != 0:
        sys.exit(f"AMBER run failed for {mol}")


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

            run_amber(mol, xyz_path)
            custom_energies = run_custom(mol, xyz_path)

            amber_mol_name = xyz_path.stem  # matches workdir created by the bash script
            sp_out    = REPO_ROOT / "amber_runs" / amber_mol_name / "sp.out"
            amber_h5  = REPO_ROOT / "results" / group / mol / "amber.h5"
            custom_h5 = REPO_ROOT / "results" / group / mol / "custom.h5"

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
