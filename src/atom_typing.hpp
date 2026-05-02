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

// functions in atom_typing.cpp 
void type_carbon_backbones(MoleculeGraph& mol);
void type_hydrogens(MoleculeGraph& mol);
void assign_amber_types(MoleculeGraph& mol);
