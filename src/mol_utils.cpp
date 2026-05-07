#include "mol_utils.hpp"
#include "molecule.hpp"
#include "read_params.hpp"
#include "geom_util.hpp"
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

    double d = dist(a, b);
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


static void type_halogens(MoleculeGraph& mol)
{
    for (int i = 0; i < mol.num_atoms(); i++)
    {
        if (mol.atoms[i].element == "Cl")
            mol.atoms[i].amber_type = "cl";
        else if (mol.atoms[i].element == "Br")
            mol.atoms[i].amber_type = "br";
    }
}

static void type_oxygen_and_nitrogen(MoleculeGraph& mol)
{
    for (int i = 0; i < mol.num_atoms(); i++)
    {
        const auto& bonds = mol.get_bonds(i);
        int num_bonds = bonds.size();

        if (mol.atoms[i].element == "O")
        {
            // OH hydroxyl
            if (num_bonds == 2)
                mol.atoms[i].amber_type = "oh";
            // =O in carbonyl
            else
                mol.atoms[i].amber_type = "o";
        }
        if (mol.atoms[i].element == "N")
        {
            // sp2 amide N
            if (num_bonds == 2)
                mol.atoms[i].amber_type = "n";
            // sp3 amine/amide N
            else
                mol.atoms[i].amber_type = "n3";
        }
    }
}


/*
Classify carbons based on bond distances to neighbors.
Priority: c1 (sp, triple) > c (sp2, C=O/C=N) > c2 (sp2, C=C) > c3 (sp3)

Distance thresholds from gaff2.dat:
  C-C triple  (c1-c1): ~1.21 A  -> d < 1.28
  C-C double  (c2-c2): ~1.34 A  -> 1.28 <= d < 1.45
  C=O                : ~1.22 A  -> d < 1.28
  C-N triple  (c1-n1): ~1.16 A  -> d < 1.22
  C=N double  (c=n)  : ~1.28 A  -> 1.22 <= d < 1.38
*/
void type_carbon_backbones(MoleculeGraph& mol)
{
    for (int i = 0; i < mol.num_atoms(); i++)
    {
        if (mol.atoms[i].element != "C") continue;

        bool has_triple   = false;
        bool has_C_double = false;
        bool has_O_double = false;
        bool has_N_double = false;

        for (size_t j : mol.get_bonds(i))
        {
            const Atom& nb = mol.atoms[j];
            double d = dist(mol.atoms[i], nb);

            if (nb.element == "C") {
                if (d < 1.28)       has_triple   = true;
                else if (d < 1.45)  has_C_double = true;
            }
            else if (nb.element == "O") {
                if (d < 1.28)       has_O_double = true;
            }
            else if (nb.element == "N") {
                if (d < 1.22)       has_triple   = true;
                else if (d < 1.38)  has_N_double = true;
            }
        }

        if (has_triple)
            mol.atoms[i].amber_type = "c1";
        else if (has_O_double || has_N_double)
            mol.atoms[i].amber_type = "c";
        else if (has_C_double)
            mol.atoms[i].amber_type = "c2";
        else
            mol.atoms[i].amber_type = "c3";
    }
}

/*
Type hydrogens based on what they are bonded to.
Must run after carbons, oxygens, and nitrogens are typed.
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

        const auto& nbrs = mol.get_bonds(i);
        if (nbrs.empty()) { mol.atoms[i].amber_type = "H"; continue; }

        const Atom& nb = mol.atoms[nbrs[0]];  // H always has exactly 1 bond

        if (nb.element == "O")
            mol.atoms[i].amber_type = "ho";
        else if (nb.element == "N")
            mol.atoms[i].amber_type = "hn";
        else if (nb.element == "C") {
            if (nb.amber_type == "c3")
                mol.atoms[i].amber_type = "hc";
            else  // c, c2, c1
                mol.atoms[i].amber_type = "ha";
        }
        else
            mol.atoms[i].amber_type = "H";
    }
}

void assign_amber_types(MoleculeGraph& mol)
{
    type_carbon_backbones(mol);      
    type_oxygen_and_nitrogen(mol);
    type_halogens(mol);
    type_hydrogens(mol);            
}