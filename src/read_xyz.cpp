#include "read_xyz.hpp"
#include "read_params.hpp"
#include "mol_utils.hpp"
#include "util.hpp"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>


MoleculeGraph read_xyz(const std::string& filepath)
{
    std::ifstream f(filepath);

    if (!f.is_open())
    {
        throw std::runtime_error("Cannot open: " + filepath);
    }

    MoleculeGraph mol;

    std::string line;

    // line 1: atom count
    std::getline(f, line);
    std::istringstream count_iss(line);
    int num_atoms = 0;

    if (!(count_iss >> num_atoms))
    {
        throw std::runtime_error("Atom count needed on line 1!");
    }

    // line 2: comment (discard, like in the homeworks)
    std::getline(f, line);

    // Lines 3-A: atom rows (element x y z)
    // future: read till end of file and make sure num_atoms matches num lines
    for (int i = 0; i < num_atoms; i++)
    {
        std::getline(f, line);
        std::istringstream iss(line);
        Atom a;
        double x, y, z;

        // assign coordinates to position vec  
        if (!(iss >> a.element >> x >> y >> z))
        {
            throw std::runtime_error("Malformed atom row: " + line);
        }
        a.position = arma::vec({x, y, z});

        // assign default to amber_type
        a.amber_type = "n/a";

        mol.add_atom(a);
    }

    // assign bonds based on the coordinates
    assign_bonds(mol); 

    // assign amber_type based on bonds 
    assign_amber_types(mol); 
    

    // OLD CODE (HARD CODED BONDS)
    // Lines A-Z: bond rows (atom_i atom_j)
    //std::getline(f, line); // section header

    //while (std::getline(f, line))
    //{
    //   std::istringstream iss(line);
    //    int i, j;

    //    if (!(iss >> i >> j))
    //    {
    //        throw std::runtime_error("Malformed bond row: " + line);
    //    }

    //    mol.add_bond(i, j);
    //}

    return mol;
}
