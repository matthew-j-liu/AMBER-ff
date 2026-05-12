# AMBER-ff

This is a from-scratch C++ implementation of the AMBER molecular-mechanics force field. The input is the coordinates of a small organic molecule (XYZ file) and the force field parameter file, and the amber_ff binary computes the total potential energy as the sum of four terms: bond stretching, angle bending, dihedral torsion, and van der Waals. 

The repository also has a benchmarking script that compares this custom implementation against the reference AmberTools pipeline (obabel, antechamber, parmchk2, tleap, sander) across 40 small molecules covering alkanes, alkenes, alkynes, alcohols, amines, carboxylic acids, halogenated compounds, and amino acids. It also produces parity plots, error plots, and a torsion study from the same data.

The two pipelines compared throughout this document are:
- AMBER reference — the AmberTools 24.8 pipeline: obabel, antechamber, parmchk2, tleap, sander
- Custom implementation — the amber_ff C++ binary built from this repo

The complete list of molecules tested is included in input_molecules_copy/list_of_molecules.md

## Prerequisites

From the repo root, build the image once (takes about 5 min):
```bash
docker build -t amber-ff .
```
```bash
docker build -t amber-ff .
```
On Apple Silicon Macs, add `--platform=linux/amd64` to the command above.

Then open a shell inside the container. 
```bash
docker run --rm -it \
    -v "$PWD/amber_runs:/work/amber_runs" \
    -v "$PWD/results:/work/results" \
    -v "$PWD/benchmarks/figures:/work/benchmarks/figures" \
    amber-ff
```

Run all commands below from inside that container shell.

---

## Demonstrating one molecule at a time

### AMBER reference pipeline
By molecule name 

```bash
./run_AMBER_benchmarking.sh ethane
```

By full xyz path
```bash
./run_AMBER_benchmarking.sh input_molecules_copy/straight_chain_alkanes/ethane.xyz
```
The script prints the BOND / ANGLE / DIHED / VDWAALS energy block from sander at the end.

### Custom implementation

By molecule name (searches input_molecules_copy/ subdirectories automatically)
```bash
./run_AMBER_custom.sh ethane
```

By full xyz path with an explicit params file
```bash
./run_AMBER_custom.sh input_molecules_copy/straight_chain_alkanes/ethane.xyz params/selected_atoms.dat
```

The script prints the BOND / ANGLE / DIHED / VDWAALS energy at the end.

---

## Running both pipelines as a loop

### Loop over all 40 molecules

**AMBER reference loop** 
```bash
./run_AMBER_benchmarking.sh
```

**Custom implementation loop**
```bash
./run_AMBER_custom.sh
```

### Full parity suite across all functional groups

`benchmarks/run_parity.py` runs both pipelines for every molecule in every group, writes
results to HDF5 under `results/`, and generates parity plots under `benchmarks/figures/`.

**All groups at once:**
```bash
python benchmarks/run_parity.py
```

**One or more specific groups:**
```bash
python benchmarks/run_parity.py straight_chain_alkanes
python benchmarks/run_parity.py straight_chain_alkanes alcohols alkynes
```

Available groups:

| Group | Molecules |
|---|---|
| `straight_chain_alkanes` | methane, ethane, propane, n-butane, n-pentane, n-hexane |
| `branched_alkanes` | isobutane, isopentane, neopentane, 2-methylpropane, 2,3-dimethylbutane |
| `alkenes` | ethene, propene, isobutylene, 2-methyl-1-propene, but-1-ene |
| `alkynes` | ethyne, propyne, 1-butyne, 3-methyl-but-3-en-1-yne |
| `alcohols` | methanol, ethanol, isopropanol, 1-propanol, ethylene glycol, glycerol |
| `aldehydes` | propanal |
| `amines` | methylamine, ethylamine, azepane |
| `carboxylic_acids` | acetic acid, glycolic acid, cis/trans acrylic acid |
| `halogenated` | methylene chloride, ethyl chloride |
| `amino_acids` | glycine, alanine |

Outputs per molecule:
- `results/<group>/<mol>/amber.h5` — AMBER reference energies (bond, angle, dihed, vdw, total)
- `results/<group>/<mol>/custom.h5` — custom binary energies (same terms)

Outputs per group:
- `benchmarks/figures/<group>_parity.png` — parity plot (custom vs. AMBER, one panel per energy term)

---

## Dihedral scan

Runs the C++ `dihedral_scan` binary on but-1-ene, sweeping the dihedral -180 to 180° and reporting all energy terms vs. angle.

```bash
python dihedral_scan/scan_dihedral.py
```

Outputs:
- `dihedral_scan/results/But-1-ene_scan.csv` — angle, bonds, angles, dihedrals, vdw, total
- `dihedral_scan/results/But-1-ene_scan.png` — all terms plotted vs. dihedral angle

Requires `build/dihedral_scan` 

---

## Plotting from existing results

Both scripts below read the HDF5 files written by `benchmarks/run_parity.py` under `results/<group>/<mol>/{amber,custom}.h5` and write PNGs under `benchmarks/figures/`. Run `run_parity.py` first.

### Per-functional-group error plots

```bash
python benchmarks/plot_errors.py
```

Outputs to `benchmarks/figures/error_rmse_small_multiples.png`

### Parity plot colored by chemical class

```bash
python benchmarks/plot_by_class.py
```

Outputs to `benchmarks/figures/by_class/`:
- `parity_<term>.png` — one per energy term (bond, angle, dihedral, vdw, total)

---

## Script summary

| Script | What it runs | Single molecule | Full loop |
|---|---|---|---|
| `run_AMBER_benchmarking.sh [mol]` | AMBER pipeline | `./script ethane` | `./script` |
| `run_AMBER_custom.sh [mol] [params]` | custom `amber_ff` | `./script ethane` | `./script` |
| `python benchmarks/run_parity.py [groups]` | both + HDF5 + plots | `… straight_chain_alkanes` | `python …` |
| `python dihedral_scan/scan_dihedral.py` | `dihedral_scan` binary on but-1-ene + plot | — | `python …` |
| `python benchmarks/plot_errors.py` | per-group error plots from HDF5 | — | `python …` |
| `python benchmarks/plot_by_class.py` | parity plots colored by class | — | `python …` |
