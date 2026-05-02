#pragma once
#include <vector>
#include <string>
#include <armadillo>


// Node in the molecular graph
struct Atom
{
    std::string element;    // e.g., "H", "C"
    arma::vec   position;   // Angstroms, 3D
    std::string amber_type; // AMBER code. For example, "c3" is sp3 carbon 
};


// Representing the molecule as an undirected graph, using adjacency lists: 
// atoms are nodes, bonds are edges.
// adjacency[i] holds the indices of all atoms bonded to atom i.
struct MoleculeGraph
{
    std::vector<Atom>              atoms;      // nodes
    std::vector<std::vector<int>>  adjacency;  // edges

    int  add_atom(const Atom& a);
    void add_bond(int i, int j);
    const std::vector<int>& neighbors(int i) const;
    int  num_atoms() const;
    int  num_bonds() const;
};
