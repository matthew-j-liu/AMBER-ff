// functions to help identify the geometry of the molecule from the coordinates

/* 
Ideas for future improvement:
- do not loop over all atoms for each element type
- start with backbone and stick to that only
- use Boost graph functions for traversal etc (MoleculeGraph should inherit from Boost graphs)
- make code more generalizable for bond orders
*/ 

#include "mol_utils.hpp"
#include "molecule.hpp"
#include "parsing.hpp"
#include "geom_util.hpp"
#include <algorithm>
#include <stdexcept>
#include <string>

bool is_bond(const Atom& a, const Atom& b)
{
    std::string e1 = a.element, e2 = b.element;
    // arrange so lower index is always first
    if (e1 > e2) std::swap(e1, e2);

    // std::map.find(), returns iterator to std::pair<key, value>, accessed by first and second
    // use hardcoded element wise lookup, since amber type has not been assigned yet
    auto actual_bond_length = BOND_LENGTHS.find({e1, e2});
    if (actual_bond_length == BOND_LENGTHS.end()) return false;

    const auto& refs = actual_bond_length->second;
    double lo = *std::min_element(refs.begin(), refs.end()) * 0.9;
    double hi = *std::max_element(refs.begin(), refs.end()) * 1.1;

    double d = bond_length(a, b);
    // return true only if d is within range of bond limits
    return (lo <= d && hi >= d);
}

// loop through all atom pairs and consider it a bond if it lies within 10% of the expected bond lengths
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

/*
Assign amber_types to each atom based on the bonds information we have now in the molecule.
The element wise order of assignment is important.
Smarter way is to traverse the backbone and then assign the neighbors as you go, instead of looping through all atoms each time

Types Considered
------------------
c   = Sp2 C carbonyl group
c1  = Sp C (alkyne, nitrile)
ha  = H bonded to sp2/sp C
c2  = Sp2 C (non-aromatic)
c3  = Sp3 C
hc  = H bonded to aliphatic C without electronegative groups
hn  = H bonded to N
ho  = H in hydroxyl group
n4  = Sp3 N with four connected atoms
o   = O with one connected atom (carbonyl)
oh  = O in hydroxyl group
cl  = Chlorine
br  = Bromine
f   = Fluorine
*/
void assign_amber_types(MoleculeGraph& mol)
{
    type_halogens(mol);
    type_oxygen_and_nitrogen(mol);
    type_carbon_backbones(mol);
    type_hydrogens(mol);
}

// Step 1: we can assign halogens right away 
void type_halogens(MoleculeGraph& mol)
{
    for (int i = 0; i < mol.num_atoms(); i++)
    {
        if (mol.atoms[i].element == "Cl")
            mol.atoms[i].amber_type = "cl";
        else if (mol.atoms[i].element == "Br")
            mol.atoms[i].amber_type = "br";
        else if (mol.atoms[i].element == "F")
            mol.atoms[i].amber_type = "f";
    }
}


// Oxygen: single bond means -OH (oh) (no ethers), double bond means C=O (c) 
// Nitrogen: assume 4 bonds
void type_oxygen_and_nitrogen(MoleculeGraph& mol)
{
    for (int i = 0; i < mol.num_atoms(); i++)
    {
        const auto& bonds = mol.get_bonds(i);
        int num_bonds = bonds.size();

        if (mol.atoms[i].element == "O")
        {
            // OH hydroxyl, no ethers. Otherwise check if the neighbor is H or C. 
            if (num_bonds == 2)
                mol.atoms[i].amber_type = "oh";
            // =O in carbonyl
            else
                mol.atoms[i].amber_type = "o";
        }
        if (mol.atoms[i].element == "N")
        {
            // sp3 amine N with lone pair
            // amides etc not considered here
            if (num_bonds == 3)
                mol.atoms[i].amber_type = "n8";
        }
    }
}


/*
Classify carbons based on bond distances to neighbors, and number of atoms.
Priority: c1 (sp, triple) > c (sp2, C=O) > c2 (sp2, C=C) > c3 (sp3)

Distance thresholds from selected_atoms.dat/ gaff2.dat:
(Hardcoded here)

C-C triple  (c1-c1): ~1.21 A  -> d < 1.28
C-C double  (c2-c2): ~1.34 A  -> 1.28 <= d < 1.45
C-C single  (c3-c3): ~1.54 A  -> d >= 1.45
C=O double  : ~1.22 A  -> d < 1.28
C-O single  (c3-oh): ~1.43 A  -> d >= 1.28  (c-oh in carboxylate: ~1.36 A)
C-N single  (c3-n4): ~1.52 A  -> any distance (n4 is sp3, no C=N in scope)
*/
void type_carbon_backbones(MoleculeGraph& mol)
{
    for (int i = 0; i < mol.num_atoms(); i++)
    {
        if (mol.atoms[i].element != "C") continue;
        
        // assume all carbons are aromatic.
        if (mol.special_notes == "aromatic")
        {
            mol.atoms[i].amber_type = "ca";
        }

        // assume all carbons are in conjugated double bonds
        if (mol.special_notes == "conjugated")
        {
            mol.atoms[i].amber_type = "ce";
        }
        else
        {
            bool has_triple   = false;
            bool has_C_double = false;
            bool has_O_double = false;
            bool has_C_single = false;
            bool has_O_single = false;
            bool has_N_single = false;

            for (size_t j : mol.get_bonds(i))
            {
                const Atom& c_neighbor = mol.atoms[j];
                double d = bond_length(mol.atoms[i], c_neighbor);
                // find c_neighbors here and if halogen or H, assign its type. 
                // Prevents 2 additional loops through the molecule

                if (c_neighbor.element == "C") {
                    if (d < 1.28)       has_triple   = true;
                    else if (d < 1.45)  has_C_double = true;
                    else                has_C_single = true;   
                }
                else if (c_neighbor.element == "O") {
                    if (d < 1.28)       has_O_double = true;
                    // O is now carbonyl (co) 
                    else                has_O_single = true;   
                    // O is now hydroxyl or ether
                }
                else if (c_neighbor.element == "N") {
                    has_N_single = true;   // n4 is sp3 only — no C=N double/triple in scope 
                }
            }

            if (has_triple)
                mol.atoms[i].amber_type = "c1";

            else if (has_O_double)
                mol.atoms[i].amber_type = "c";

            else if (has_C_double)
                mol.atoms[i].amber_type = "c2";

            else
                mol.atoms[i].amber_type = "c3";
        }
    }
}

/*
Lastly, type hydrogens based on what they are bonded to.

hc — H on sp3 C (c3)
ha — H on sp2/sp C (c, c2, c1)
ho — H on O
hn — H on N
*/
void type_hydrogens(MoleculeGraph& mol)
{
    for (int i = 0; i < mol.num_atoms(); i++)
    {
        if (mol.atoms[i].element != "H") continue;

        const auto& h_neighbors = mol.get_bonds(i);
        if (h_neighbors.empty()) { mol.atoms[i].amber_type = "H"; continue; }

        const Atom& neighbor = mol.atoms[h_neighbors[0]];  // H always has exactly 1 bond

        if (neighbor.element == "O")
            mol.atoms[i].amber_type = "ho";
        else if (neighbor.element == "N")
            mol.atoms[i].amber_type = "hn";
        else if (neighbor.element == "C") {
            if (neighbor.amber_type == "c3")
                mol.atoms[i].amber_type = "hc";
            else  // c, c2, c1, ca
                mol.atoms[i].amber_type = "ha";
        }
        else
            mol.atoms[i].amber_type = "H";
    }
}
