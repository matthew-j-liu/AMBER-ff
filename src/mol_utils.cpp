#include "mol_utils.hpp"
#include "read_params.hpp"
#include <algorithm>
#include <stdexcept>
#include <string>

// check this function 
bool is_bond(const Atom& a, const Atom& b)
{
    std::string e1 = a.element, e2 = b.element;
    if (e1 > e2) std::swap(e1, e2);

    auto it = BOND_LENGTHS.find({e1, e2});
    if (it == BOND_LENGTHS.end()) return false;

    const auto& refs = it->second;
    double lo = *std::min_element(refs.begin(), refs.end()) * 0.9;
    double hi = *std::max_element(refs.begin(), refs.end()) * 1.1;

    double d = dist(a.position, b.position);
    return d >= lo && d <= hi;
}


void assign_bonds(MoleculeGraph& mol)
{
    for (int i = 0; i < mol.num_atoms(); i++) {
        for (int j = i + 1; j < mol.num_atoms(); j++) {
            if (is_bond(mol.atoms[i], mol.atoms[j])) {
                mol.add_bond(i, j);
            }
        }
    }
}


// only H, Cl and Br 
// make this code better, no if statements like this
void type_hydrogen_and_halogens(MoleculeGraph& mol)
{
    for (int i = 0; i < mol.num_atoms(); i++)
    {
        if (mol.atoms[i].element == "H") 
        { 
            mol.atoms[i].amber_type = "H";
        }   
        if (mol.atoms[i].element == "Cl") 
        { 
            mol.atoms[i].amber_type = "Cl";
        }    
        if (mol.atoms[i].element == "Br") 
        { 
            mol.atoms[i].amber_type = "Br";
        }      
    }
}

void type_oxygen_and_nitrogen(MoleculeGraph& mol)
{
    for (int i = 0; i < mol.num_atoms(); i++)
    {
        // count bonds
        std::vector<int> bonds = mol.bonds(i);
        int num_bonds = bonds.size(); 

        if (mol.atoms[i].element == "O") 
        { 
            // OH hydroxyl
            if (num_bonds == 2){
                mol.atoms[i].amber_type = "oh"; 
            }
            // =O in carbonyl
            else {
                mol.atoms[i].amber_type = "o"; 
            }
        }   
        // check logic here
        if (mol.atoms[i].element == "N") 
        { 
            // double bond amide
            if (num_bonds == 2){
                mol.atoms[i].amber_type = "n"; 
            }
            // single bond amine or amide
            else {
                mol.atoms[i].amber_type = "n3"; 
            }
        }    
    }
}


/*
Classify carbons based on the neighbors
*/
void type_carbon_backbone(MoleculeGraph& mol)
{
    for (int i = 0; i < mol.num_atoms(); i++)
    {

    }
}

void assign_amber_types(MoleculeGraph& mol)
{
    type_hydrogen_and_halogens(mol); // only 1 option
    type_oxygen_and_nitrogen(mol);
    type_carbon_backbone(mol);
}