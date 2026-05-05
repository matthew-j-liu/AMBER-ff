"""
Run main.cpp and produce parity plots
"""

import os
import re # Python's built in regular expressions module 
import sys
from pathlib import Path
import h5py

# default list of molecules to run
ALKANES = ["ethane", "propane", "n-butane", "n-pentane", "n-hexane"]


def parse_energies(sp_out: Path) -> dict[str, float]:
    """
    Extract bond/angle/dihed/vdw/1-4 NB from a single .out file
    """
    text = sp_out.read_text()
    def grab(label: str): 
        """
        specify `label` to be the energy term that we want to extract 
        """
        m = re.search(rf"{label}\s*=\s*(-?\d+\.\d+)", text) # derived the regex by examining the .out file 
        if m is None:
            raise RuntimeError(f"missing '{label}' in {sp_out}")
        return float(m.group(1))
    return {
        "bond": grab("BOND"),
        "angle": grab("ANGLE"),
        "dihed": grab("DIHED"),
        "vdw": grab("VDWAALS"),
        "nb14": grab("1-4 NB"),
    }


def write_amber_h5(path: Path, terms: dict[str, float]) -> None:
    """
    Write AMBER reference energies to HDF5 with dataset names matching our C++ implementation
    """
    vdw_total = terms["vdw"] + terms["nb14"]
    total = terms["bond"] + terms["angle"] + terms["dihed"] + vdw_total
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f.create_dataset("bonds_energy", data=terms["bond"])
        f.create_dataset("angles_energy", data=terms["angle"])
        f.create_dataset("dihedrals_energy", data=terms["dihed"])
        f.create_dataset("vdw_energy", data=vdw_total)
        f.create_dataset("total_energy", data=total)


def run_amber(mol: str) -> None:
    print(f"[AMBER] {mol}")
    if os.system(f"./run_AMBER_benchmarking_simple.sh {mol}") != 0:
        sys.exit(f"AMBER run failed for {mol}")


def run_custom(mol: str) -> None:
    print(f"[custom] {mol}")
    # exe writes results/{mol}/custom.h5 relative to cwd
    if os.system(f"./build/amber_ff input_molecules/{mol}.xyz params/hydrocarbons.dat") != 0:
        sys.exit(f"custom run failed for {mol}")


def main():
    # pass molecule names on the command line (no arguments defaults to ALKANES)
    mols = sys.argv[1:] or ALKANES

    for mol in mols:
        print(mol + ":")
        run_amber(mol)
        run_custom(mol)
        amber_h5 = Path("results") / mol / "amber.h5"
        sp_out = Path("amber_runs") / mol / "sp.out"
        write_amber_h5(amber_h5, parse_energies(sp_out))
        print(f"wrote {amber_h5}")

    print("\nplotting parity figures...")
    from plot_parity import plot_all
    for q in plot_all(mols):
        print(f"  wrote {q}")


if __name__ == "__main__":
    main()
