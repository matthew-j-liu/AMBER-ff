#include "molecule.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

// in the constructor add num_atoms = num_bonds = 0 (default zero initialize)
MoleculeGraph::MoleculeGraph(){
    
}

int MoleculeGraph::add_atom(const Atom& a)
{
    atoms.push_back(a);
    adjacency.emplace_back();
    return static_cast<int>(atoms.size()) - 1;
}

void MoleculeGraph::add_bond(int i, int j)
{
    adjacency[i].push_back(j);
    adjacency[j].push_back(i);
}

const std::vector<int>& MoleculeGraph::neighbors(int i) const
{
    return adjacency[i];
}

int MoleculeGraph::num_atoms() const
{
    return static_cast<int>(atoms.size());
}

int MoleculeGraph::num_bonds() const
{
    int total = 0;
    for (const auto& nbrs : adjacency) total += static_cast<int>(nbrs.size());
    return total / 2;
}
