## Chem 279 | MSSE | Spring 2026
### Matthew Liu and Jahnavi Gandhi

### ML 2026-04-28 log

* I downloaded a database of .xyz files from: https://www.skies-village.it/webtools/databases/LCB25/
    * Then I pasted in the following subset of .xyz files into the repo:
        * ethane.xyz
        * propane.xyz
        * n-butane.xyz
        * n-pentane.xyz
        * n-hexane.xyz
* For each of the .xyz files, I added additional lines after the coordinates that encode which atoms are bonded to one another. For example, here is the ethane file with my annotations: 
```
8 // the number of atoms in the molecule 

C -0.000000 0.761550 -0.000000 // xyz coordinates of atom 0 
C -0.000000 -0.761550 -0.000000 // xyz coordinates of atom 1
H -1.015397 1.155664 -0.000000 // xyz coordinates of atom 2 
H 0.507699 1.155662 0.879360 // etc. 
H 0.507699 1.155662 -0.879360
H 1.015397 -1.155664 0.000000
H -0.507699 -1.155662 0.879360
H -0.507699 -1.155662 -0.879360
* Bonds // note that we can eventually add a column for bond order. for now, all implicitly single bonds 
0 1 // means we have a bond between atom 0 and atom 1, i.e., carbon to carbon
0 2 // bond between atom 0 and atom 2, i.e., C-H bond 
0 3
0 4
1 5
1 6
1 7
```

* From here, I made two independent header files. I anticipate using both in downstream code. Please see the source files for each, which have docstrings detailing the logic.
    * `ff_params.hpp`: ready to store AMBER force field parameters.
    * `molecule.hpp`: used composition to define a `Molecule` struct based on many `Atom` structs. 


### ML 2026-04-29 log 
* Given that `gaff2.dat` is this monolith of parameters, of which we only need a small subset, I extracted just the rows relevant to saturated hydrocarbons. The smaller parameter file is `hydrocarbons.dat`, where parameter sections are delimited by bracketed headers to help with parsing.
* Wrote a parser for `hydrocarbons.dat` in `read_params.cpp`. A single function, `load_params(filepath)`, uses 4 helper functions (one per parameter struct). 


### JG 2026-04-29 log 
* implemented the basic version of a molecule as a graph structure - nodes as the atoms and egdes as the bonds associated with it. 
    * on initialization, it needs functions that identify which atoms are bonded (is there a threshold to use?)
* implemented the functions and associated helper functions to calculate each term in the FF. These can all be found in ff.cpp
    * Imp Question: where do the A, B, C, D parameters come from? I don't think we read them from the .dat
* implemented functions to calculate bond lengths, bond angles, and dihedral angles from coordinates file and molecule as a graph

### ML 2026-04-30 log 
* Fixed a few bugs in `ff.cpp` and `util.cpp` that would prevent the energy code from compliling and linking. For example, some methods were marked `static` that shouldn't have been. Some additional small changes
    * Moved `DEG2RAD` INTO `util.hpp` as a constant 
    * Fixed variable name mismatch in `calculate_bonds_term` 
    * Parameter lookup happens through the `ForceField` map, so got rid of the 4 placeholder functions in `ff.cpp` 
* To your question about how A and B terms are derived - the params file assumes one form of the LJ equation using $R^*$ and $\epsilon$. The equation we've been looking at uses another form with $A$ and $B$. The two are related as:

    ``` 
    double A = eps_ij * std::pow(R_star_ij, 12);
    double B = 2.0 * eps_ij * std::pow(R_star_ij, 6); 
    ```

    Note I have added this conversion in `ff.cpp`. 

* I wrote `read_xyz.cpp` leveraging the `MoleculeGraph` and got rid of the original `Molecule` struct I had in `molecule.hpp`.
* To verify both parameter and xyz parsing, I filled in `main` with a quick test case on ethane. The correct parameters, atom positions, and bond information are printed out. I have yet to test on the other hydrocarbons, but this was encouraging! Note that I compiled and ran as follows:

    ```
    g++ -std=c++17 src/main.cpp src/read_xyz.cpp src/read_params.cpp src/molecule.cpp src/util.cpp src/ff.cpp -o main -larmadillo 

    ./main input_molecules/ethane.xyz params/hydrocarbons.dat
    ```

    For your convenience, here was the output:

    ```
    Molecule loaded!
    atoms: 8
    bonds: 7

    Forcefield parameters loaded!
    # of bond params: 2
    # of angle params: 3
    # of dihedral params: 3
    # of vdw params: 2
    0  C  (-0, 0.76155, -0)  neighbors={1,2,3,4}
    1  C  (-0, -0.76155, -0)  neighbors={0,5,6,7}
    2  H  (-1.0154, 1.15566, -0)  neighbors={0}
    3  H  (0.507699, 1.15566, 0.87936)  neighbors={0}
    4  H  (0.507699, 1.15566, -0.87936)  neighbors={0}
    5  H  (1.0154, -1.15566, 0)  neighbors={1}
    6  H  (-0.507699, -1.15566, 0.87936)  neighbors={1}
    7  H  (-0.507699, -1.15566, -0.87936)  neighbors={1}
    ```
### ML 2026-04-30 log (part II)
* Added `atom_typing.cpp` to assign sp3 carbon and `hc` hydrogen in saturated hydrocarbons. Throws an error if there are other types of atoms, we can augment the code to identiy other atom types later. The logic is quite simple, just using the element identity and number of bonds associated to that element
* I think from here, we have all the pieces in place to calculate energy of our hydrocarbons! I have ethane up to hexane in the `input_molecules` directory.

### JG 2026-04-30 log

1. on branch jahnavi-misc, I have started with some code for making and building the repo + some bash scripts for easy build_and_run. Also added some logic for to implement the atom type ID.  
2. Tested all 5 hydrocarbons currently considered. Edited code in main for this. 
Note to myself: Compilation code for me
    ```
    g++ -std=c++17 
    src/main.cpp src/read_xyz.cpp src/read_params.cpp src/molecule.cpp src/util.cpp src/ff.cpp 
    -I$CONDA_PREFIX/include 
    -L$CONDA_PREFIX/lib -Wl, -rpath, $CONDA_PREFIX/lib
    -o main 
    -larmadillo
    ```

### JG 2026-05-01 log
1. added some starting code for atom typing 
2. adding some starting code for unit tests using google test suite and containerization 
3. merged branches jahnavi-misc and jahnavi_updates (no reason why separate)



### JG 2026-05-04 log
1. Over the weekend, I worked on installing the actual AMBER FF package, and got it running to calculate energies. I have a bash file that you can run that reads the output file and prints the final energy readout. Important thing is to compare term by term and not the total energy as electrostatic interactions are not accounted for in our implementation (those params not in .gaff file)
2. updated gaff.dat file (and hydrocarbons.dat) to the version used by the amber package that I installed (so that results match exactly). 
3. updated these results in slides shared with you!
4. updated van der waals term calculation to fix an error

### Next steps
* Complete atom typing so we can do more molecules 
* think about what next/ how to present results now that we can calculate energies

### Notes/ ongoing questions
* Modern FF do not use the explicit H bond term. 
* where do the A, B, C, D parameters come from? I don't think we read them from the .dat
    * Matthew's note (2026-04-30): I addressed $A$ and $B$ in my log, good question for C and D. I think as you noted above, various versions of AMBER do not have an explicit hydrogen bonding term. Rather, many parameters have the H-bonding folded into the parameter itself. I forget which is the case for our .dat file, but for now I propose we don't have the H-bonding term for simplicity. 



