#include "atom_typing.hpp"
#include <stdexcept>
#include <string>

/*
classify the carbon backbone of hydrocarbon, with
plans to generalize to other elements later 
*/
void type_heavy_atoms(MoleculeGraph& mol)
{
    for (int i = 0; i < mol.num_atoms(); i++)
    {
        const std::string& elem = mol.atoms[i].element;
        if (elem == "H") { continue; }

        else if (elem == "C")
        {
            size_t n_bonds = mol.neighbors(i).size();
            if (n_bonds != 4) // only supporting sp3 carbons for now 
            {
                throw std::runtime_error("non-sp3 carbon currently not supported!"); 
            }
            mol.atoms[i].amber_type = "c3";
        }

        else
        {
            throw std::runtime_error("Atom " + std::to_string(i) + " has unsupported element: " + elem);
        }
    }
}


/*
Each H can only be bonded to one (heavy) atom.
The identity of that atom dictates the identity of the H 
For now, only support "hc", look into expanding later 
*/
void type_hydrogens(MoleculeGraph& mol)
{
    for (int i = 0; i < mol.num_atoms(); i++)
    {
        if (mol.atoms[i].element != "H") { continue; }

        std::vector<int> neighbors = mol.neighbors(i);
        std::string neighbor_type = mol.atoms[neighbors[0]].amber_type;
        
        if (neighbor_type != "c3")
        {
            throw std::runtime_error("non-hc hydrogen currently not supported!"); 
        }

        mol.atoms[i].amber_type = "hc";
    }
}


/*
leverage both functions, first identifying the backbone and then the hydrogens
*/
void assign_amber_types(MoleculeGraph& mol)
{
    type_heavy_atoms(mol);
    type_hydrogens(mol);
}