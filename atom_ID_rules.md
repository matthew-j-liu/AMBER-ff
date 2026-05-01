# Atom Type ID rules for FF 

## Atom connections to consider: 

C3-C3 
C2-C2 (	non-aromatic)
C3-C2 
C3-Cl
C2-Cl
C2-Br 
C3-Br 
C3-O (sp3)
C2-O 
C2=O
C-H
C-N
C=N
C=-N (triple bond)
N-H (consider only terminal)
O-H 

## Logic: 
0. Create a lookup for bondlengths to be considered (from gaff.dat in params). std::map, in util.hpp 
In constructor of MoleculeGraph: 
1. Create a backbone of the molecules (all Cs) from the molecule. Vector called backbone, if C, add to backbone
2. Use this backbone and find distances between all other atoms and Cs 
3. For each C, and all other atoms (including other C), check if it meets the bond length requirement for being a bond. Also check bond order. Add to the edges (adjacency list). Bond order can be added to the graph. 
4. Loop through all the nodes and check adjacency list. Check if it follows classical chemistry rules (valency), throw an error if not. Assign atom type to each atom.  
5. How many molecules in test cases does this satisfy?

This logic ignores all molecules with N or O in the backbone (no esters). Will hydroxyls work with this?  
This logic also ignores aromatics for now. 


