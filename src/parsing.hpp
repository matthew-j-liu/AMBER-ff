#pragma once
#include "ff_params.hpp"
#include <map>
#include <string>
#include <utility>
#include <vector>
#include "molecule.hpp"


// Parses parameter file and returns a populated ForceField
ForceField load_params(const std::string& filepath);

/*
Reference bond lengths (Angstroms) from gaff2.dat, keyed by sorted element pair.
Each entry lists all bond orders for that pair (e.g. C-C: triple, double, single).
Keys are sorted alphabetically so lookup is order-independent.
*/

using ElemPair = std::pair<std::string, std::string>;
static const std::map<ElemPair, std::vector<double>> BOND_LENGTHS = {
    {{"C",  "C"}, {1.2095, 1.3185, 1.3379, 1.5088, 1.5354}}, // triple, c1c2, double, c2c3, single
    {{"C",  "H"}, {1.0962}}, // C-H
    {{"C",  "N"}, {1.1619, 1.2809, 1.3967, 1.4626}}, // triple, double, single sp2, single sp3
    {{"C",  "O"}, {1.2190, 1.2249, 1.3206, 1.3423, 1.3556, 1.4242}}, // C=O to C-OH
    {{"H",  "N"}, {1.0095, 1.0120, 1.0130, 1.0190}}, // N-H variants
    {{"H",  "O"}, {0.9725}}, // O-H
    {{"Br", "C"}, {1.9054, 1.9805}}, // C-Br: sp2, sp3
    {{"C",  "Cl"},{1.7415, 1.8092}}, // C-Cl: sp2, sp3
};


/*
Reads an .xyz file and returns a MoleculeGraph.

Expected format
----------------
    -Line 1: num_atoms
    -Line 2: comment
    -Lines 3-A: section for atomic positions
        -Each line: element x y z
    -Lines A-Z: section for bonds 
        -Each line: atom_i atom_j 
*/
MoleculeGraph read_xyz(const std::string& filepath);