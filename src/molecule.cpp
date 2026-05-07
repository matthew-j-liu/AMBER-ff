#include "molecule.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

void MoleculeGraph::add_atom(const Atom& a)
{
    atoms.push_back(a);
    bonds.push_back({});
}

// add a new bond to an atom
void MoleculeGraph::add_bond(size_t i, size_t j)
{
    bonds[i].push_back(j);
    bonds[j].push_back(i);
}

const std::vector<size_t>& MoleculeGraph::get_bonds(size_t i) const
{
    return bonds.at(i);
}

size_t MoleculeGraph::num_atoms() const
{
    return atoms.size();
}

size_t MoleculeGraph::num_bonds() const
{
    size_t total = 0;
    for (size_t i = 0; i < num_atoms(); i++)
        total += bonds[i].size();
    return total / 2;
}

// Each bond is emitted once by only including pairs where i < j.
std::vector<std::array<size_t, 2>> MoleculeGraph::find_all_bonds() const
{
    std::vector<std::array<size_t, 2>> result;
    for (size_t i = 0; i < num_atoms(); i++)
        for (size_t j : bonds[i])
            if (i < j)
                result.push_back({i, j});
    return result;
}

std::vector<std::array<size_t, 3>> MoleculeGraph::find_all_bond_angle_triplets() const
{
    std::vector<std::array<size_t, 3>> result;
    for (size_t j = 0; j < num_atoms(); j++) {
        const auto& nbrs = bonds[j];
        for (size_t i = 0; i < nbrs.size(); i++)
            for (size_t k = i + 1; k < nbrs.size(); k++)
                result.push_back({nbrs[i], j, nbrs[k]});
    }
    return result;
}

std::vector<std::array<size_t, 4>> MoleculeGraph::find_all_torsion_quadruplets() const
{
    std::vector<std::array<size_t, 4>> result;
    for (size_t j = 0; j < num_atoms(); j++) {
        for (size_t k : bonds[j]) {
            if (k <= j) continue;
            for (size_t i : bonds[j]) {
                if (i == k) continue;
                for (size_t l : bonds[k]) {
                    if (l == j || l == i) continue;
                    result.push_back({i, j, k, l});
                }
            }
        }
    }
    return result;
}


void MoleculeGraph::print_all() const
{
    std::cout << "=== Bonds (" << find_all_bonds().size() << ") ===\n";
    for (const auto& b : find_all_bonds())
        std::cout << b[0] << " -- " << b[1] << "\n";

    std::cout << "\n=== Bond Angle Triplets (" << find_all_bond_angle_triplets().size() << ") ===\n";
    for (const auto& t : find_all_bond_angle_triplets())
        std::cout << t[0] << " -- " << t[1] << " -- " << t[2] << "\n";

    std::cout << "\n=== Torsion Quadruplets (" << find_all_torsion_quadruplets().size() << ") ===\n";
    for (const auto& q : find_all_torsion_quadruplets())
        std::cout << q[0] << " -- " << q[1] << " -- " << q[2] << " -- " << q[3] << "\n";
}