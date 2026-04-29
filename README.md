## Chem 279 | MSSE | Spring 2026
### Matthew Liu and Jahnavi Gandhi

### 2026-04-28 log

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


### 2026-04-29 log
* Given that `gaff2.dat` is this monolith of parameters, of which we only need a small subset, I extracted just the rows relevant to saturated hydrocarbons. The smaller parameter file is `hydrocarbons.dat`, where parameter sections are delimited by bracketed headers to help with parsing.
* Wrote a parser for `hydrocarbons.dat` in `read_params.cpp`. A single function, `load_params(filepath)`, uses 4 helper functions (one per parameter struct). 


### Next steps
* Parsing the .xyz files into `Molecule` objects 



