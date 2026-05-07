# Parameter File Formats

## gaff2.dat (AmberTools, Version 2.2.20)

The full GAFF2 force field parameter file. Sections are separated by blank lines; there are no explicit section headers except for `MOD4` before the vdW block.

| Section | Format | Example |
|---------|--------|---------|
| **Atom types** (lines 2–100) | `type  mass [amu]  polarizability [Å³]  description` | `c3  12.01  0.878  Sp3 C` |
| **Bonds** (lines 102–1437) | `type_i-type_j  K_r [kcal/mol/Å²]  r_eq [Å]  source` | `c3-c3  228.89  1.5354  SOURCE3` |
| **Angles** (lines 1439–11150) | `type_i-type_j-type_k  K_theta [kcal/mol/rad²]  theta_eq [deg]  source` | `hc-c3-hc  35.80  107.73  SOURCE3` |
| **Proper dihedrals** (lines 11152–12968) | `type_i-type_j-type_k-type_l  divider  Vn [kcal/mol]  gamma [deg]  n`; `X` is a wildcard | `X-c3-c3-X  9  1.400  0.000  3` |
| **Improper dihedrals** (lines 12970–13009) | Same format as proper dihedrals; central atom is third position | `X-X-c-o  1  10.500  180.000  2` |
| **vdW** (`MOD4 RE`, lines 13013–13110) | `type  R* [Å]  epsilon [kcal/mol]`; `R*` is the per-atom radius, combined additively for pairs | `c3  1.9069  0.1078` |

---

## selected_atoms.dat (project-local subset)

A hand-curated subset of gaff2.dat covering the 14 atom types used in this project (c, c1, c2, c3, hc, ha, hn, ho, oh, o, n4, f, cl, br). Uses bracket headers so `load_params()` in [read_params.cpp](../src/read_params.cpp) can parse it with a simple section-dispatch loop.

| Section | Format | Example |
|---------|--------|---------|
| `[bonds]` | `type_i  type_j  K_r [kcal/mol/Å²]  r_eq [Å]` | `c3 c3  228.89  1.5354` |
| `[angles]` | `type_i  type_j  type_k  K_theta [kcal/mol/rad²]  theta_eq [deg]` | `hc c3 hc  35.80  107.73` |
| `[dihedrals]` | `type_i  type_j  type_k  type_l  divider  Vn [kcal/mol]  gamma [deg]  n`; `X` = wildcard | `X c3 c3 X  9  1.400  0.000  3` |
| `[vdw]` | `type  R* [Å]  epsilon [kcal/mol]` | `c3  1.9069  0.1078` |

Lines beginning with `#` are comments. Values are taken verbatim from the corresponding sections of gaff2.dat; the only difference is the bracket headers and the omission of source/quality annotations.
