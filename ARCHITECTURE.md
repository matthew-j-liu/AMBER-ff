# Codebase Architecture

## What it does

Reads a molecule from an XYZ file and a set of AMBER force field parameters, then computes the total potential energy as a sum of four terms: bond stretching, angle bending, dihedral torsion, and van der Waals.

---

## High-level pipeline

```
main()
  │
  ├── load_params(params_file)   →  ForceField
  │
  ├── read_xyz(mol_file)
  │     ├── assign_bonds(mol)          ← distance-based bond detection
  │     └── assign_amber_types(mol)   ← classify each atom into an AMBER type
  │
  ├── calculate_bonds_term(mol, ff)
  ├── calculate_angles_term(mol, ff)
  ├── calculate_dihedrals_term(mol, ff)
  └── calculate_vdw_term(mol, ff)
        └── total = sum of all four terms
```

---

## Source files

| File | Role |
|------|------|
| [main.cpp](src/main.cpp) | Entry point. Wires the pipeline together; prints per-term and total energies. |
| [molecule.cpp/.hpp](src/molecule.hpp) | `MoleculeGraph` data structure. Stores atoms + adjacency list; provides `find_all_bonds()`, `find_all_bond_angle_triplets()`, `find_all_torsion_quadruplets()`. |
| [ff_params.hpp](src/ff_params.hpp) | POD structs for parameters (`BondParams`, `AngleParams`, `DihedralParams`, `VdwParams`) and the `ForceField` container that maps atom-type tuples to those structs. |
| [read_params.cpp](src/read_params.cpp) | `load_params()` — parses `selected_atoms.dat` (INI-style sections) into a `ForceField`. One helper per section: `parse_bond_row`, `parse_angle_row`, `parse_dihedral_row`, `parse_vdw_row`. |
| [parsing.hpp](src/parsing.hpp) | Declares `load_params()` and `read_xyz()`; holds the hard-coded `BOND_LENGTHS` reference table used for bond detection. |
| [read_xyz.cpp](src/read_xyz.cpp) | `read_xyz()` — reads atom coordinates from an XYZ file, then calls `assign_bonds` and `assign_amber_types`. |
| [mol_utils.cpp/.hpp](src/mol_utils.hpp) | AMBER type assignment. `assign_amber_types()` calls four sub-routines in order: `type_halogens` → `type_oxygen_and_nitrogen` → `type_carbon_backbones` → `type_hydrogens`. Order matters because hydrogen typing depends on what its neighbor was already assigned. |
| [geom_util.cpp/.hpp](src/geom_util.hpp) | `bond_length()`, `angle_deg()`, `dihedral_angle()` — pure geometry using Armadillo vectors. Each has an overload that accepts `Atom` objects and delegates to the vector version. |
| [ff.cpp/.hpp](src/ff.hpp) | Energy calculations. Low-level functions (`bond_streching_energy`, `bond_rotation_energy`, `bond_torsion_energy`, `vdw_energy`) implement the AMBER functional forms; high-level functions (`calculate_*_term`) enumerate the relevant atom groups from the graph and sum the contributions. |

---

## Key data structures

```
Atom
  element       string       "C", "H", "O" ...
  position      arma::vec    x, y, z in Angstroms
  amber_type    string       "c3", "hc", "oh" ...

MoleculeGraph
  atoms         vector<Atom>
  bonds         vector<vector<size_t>>   adjacency list
  special_notes string       "aromatic" / "conjugated" / ""

ForceField
  bonds         map< (type_i, type_j),               BondParams     >
  angles        map< (type_i, type_j, type_k),        AngleParams    >
  dihedrals     map< (type_i, type_j, type_k, type_l), [DihedralParams] >
  vdw           map< type,                             VdwParams      >
```

---

## Energy term details

| Term | Functional form | Parameters from ff |
|------|----------------|-------------------|
| Bond stretch | `K_r * (r − r_eq)²` | `BondParams.K_r`, `.r_eq` |
| Angle bend | `K_θ * (θ − θ_eq)²` | `AngleParams.K_theta`, `.theta_eq` |
| Dihedral torsion | `Σ (Vn/divider) * (1 + cos(n·φ − γ))` | `DihedralParams.Vn`, `.divider`, `.gamma`, `.n`; `X` wildcards allowed |
| van der Waals | `ε·[(R*/r)¹² − 2·(R*/r)⁶]`; 1-2 and 1-3 pairs excluded; 1-4 pairs scaled by 0.5 | `VdwParams.R_star`, `.epsilon`; combined with additive R* and geometric-mean ε |

---

## AMBER type assignment order

`assign_amber_types` must run its sub-routines in this fixed order because later steps query types set by earlier ones:

1. **`type_halogens`** — `cl`, `br`, `f` (unambiguous: just check element)
2. **`type_oxygen_and_nitrogen`** — `oh`/`o` by bond count; `n4` by bond count
3. **`type_carbon_backbones`** — `c3`/`c2`/`c1`/`c` by neighbor bond distances; uses `special_notes` for aromatic/conjugated overrides
4. **`type_hydrogens`** — `ho`/`hn`/`hc`/`ha` by the already-assigned type of the atom H is bonded to

Bond detection (`assign_bonds`) runs before typing; it uses the `BOND_LENGTHS` table in `parsing.hpp` to accept a distance if it falls within ±10 % of any known bond length for that element pair.
