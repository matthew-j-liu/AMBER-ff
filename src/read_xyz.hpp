#pragma once
#include "molecule.hpp"
#include <string>

/*
Reads an .xyz file and returns a Molecule.

Expected format
----------------
    -Line 1: num_atos
    -Line 2: comment
    -Lines 3-A: section for atomic positions
        -Each line: element x y z
    -Lines A-Z: section for bonds 
        -Each line: atom_i atom_j 
*/
Molecule read_xyz(const std::string& filepath);