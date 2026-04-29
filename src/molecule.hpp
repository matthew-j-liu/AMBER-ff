#pragma once
#include <vector>
#include <array>
#include <string>
#include <armadillo>


struct Atom
{
    std::string element; // e.g., "H", "C"
    arma::vec position; // in units Angstroms 
    std::string amber_type; // AMBER code. For example, "c3" is sp3 carbon 
};


/*
We will use composition such that many Atom object comprise a Molecule object.
The molecule will have a vector of properties (i.e., bonds, angles, and dihedrals). 
We mirror the logic in `ff_params.hpp` by using std::array<int, N>, where N is 2 for bonds,
3 for angles, and 4 for dihedrals. The int denotes the atom #. For example, if we had an array
{0, 1}, this means a single bond between atoms 0 and 1 in the molecule. 
*/
struct Molecule
{
    std::vector<Atom> atoms; 
    std::string comment; 
    std::vector<std::array<int, 2>> bonds; // read directly from .xyz file 
    std::vector<std::array<int, 3>> angles; // derived from bonds (two bonds sharing a central atom)
    std::vector<std::array<int, 4>> dihedrals; // derived from bonds (chain of 3 consecutive bonds)
};