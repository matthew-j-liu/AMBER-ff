#pragma once
#include "molecule.hpp"

/*
Given a MoleculeGraph, assign each atom to an AMBER type 
based on the element and bonding environment 

Given we will be performing calculations on organic molecules, we
can divide up the atom typing into two functions:

1. The backbone 
2. The hydrogens 

Where the type of the latter depends on the former being classified first

*/

bool is_bond(const Atom& a, const Atom& b);
void assign_bonds(MoleculeGraph& mol);
void type_carbon_backbones(MoleculeGraph& mol);
void type_hydrogens(MoleculeGraph& mol);
void type_halogens(MoleculeGraph& mol);
void type_oxygen_and_nitrogen(MoleculeGraph& mol);
void assign_amber_types(MoleculeGraph& mol);







