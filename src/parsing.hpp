#pragma once
#include "ff_params.hpp"
#include <map>
#include <string>
#include <utility>
#include <vector>
#include "molecule.hpp"



// Parses parameter file and returns a populated ForceField
ForceField load_params(const std::string& filepath);

// Reference bond lengths (Angstroms) from gaff2.dat, keyed by sorted element pair.
// Each entry lists all bond orders for that pair (e.g. C-C: triple, double, single).
// Keys are sorted alphabetically so lookup is order-independent.
// Sources: c1-c1 1.2095, c1-c2 1.3185, c2-c2 1.3379, c2-c3 1.5088, c3-c3 1.5354,
//          c3-hc 1.0962, c1-n1 1.1619, c2-n2 1.2809, c-n 1.3967, c3-n3 1.4626,
//          c-o 1.2190, c2-o 1.2249, c3-o 1.3206, c2-oh 1.3423, c-oh 1.3556, c3-oh 1.4242,
//          c2-cl 1.7415, c3-cl 1.8092, br-c2 1.9054, br-c3 1.9805,
//          hn-n* ~1.01-1.02, ho-oh 0.9725
using ElemPair = std::pair<std::string, std::string>;
static const std::map<ElemPair, std::vector<double>> BOND_LENGTHS = {
    {{"C",  "C"}, {1.2095, 1.3185, 1.3379, 1.5088, 1.5354}}, // triple, c1c2, double, c2c3, single
    {{"C",  "H"}, {1.0962}},                    // C-H
    {{"C",  "N"}, {1.1619, 1.2809, 1.3967, 1.4626}}, // triple, double, single sp2, single sp3
    {{"C",  "O"}, {1.2190, 1.2249, 1.3206, 1.3423, 1.3556, 1.4242}}, // C=O to C-OH
    {{"H",  "N"}, {1.0095, 1.0120, 1.0130, 1.0190}}, // N-H variants
    {{"H",  "O"}, {0.9725}},                    // O-H
    {{"Br", "C"}, {1.9054, 1.9805}},           // C-Br: sp2, sp3
    {{"C",  "Cl"},{1.7415, 1.8092}},            // C-Cl: sp2, sp3
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