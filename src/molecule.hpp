#pragma once
#include <vector>
#include <string>
#include <array>
#include <armadillo>


// Node in the molecular graph
struct Atom
{
    std::string element;    // e.g., "H", "C"
    arma::vec   position;   // Angstroms, 3D
    std::string amber_type; // AMBER code. For example, "c3" is sp3 carbon
};


// Representing the molecule as an undirected graph, using adjacency list style to represent bonded atoms:
// atoms are nodes, bonds are edges.
// bonds[i] holds the indices of all atoms bonded to atom i.
struct MoleculeGraph
{
    std::vector<Atom> atoms;
    std::vector<std::vector<size_t>> bonds;

    void add_atom(const Atom& a);
    void add_bond(size_t i, size_t j);

    const std::vector<size_t>& get_bonds(size_t i) const;

    size_t num_atoms() const;
    size_t num_bonds() const;

    std::vector<std::array<size_t, 2>> find_all_bonds() const;
    std::vector<std::array<size_t, 3>> find_all_bond_angle_triplets() const;
    std::vector<std::array<size_t, 4>> find_all_torsion_quadruplets() const;

    void print_all() const;
};
