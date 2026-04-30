#include "molecule.hpp"
#include "read_xyz.hpp"
#include "ff_params.hpp"
#include "read_params.hpp"
#include "atom_typing.hpp"
#include <iostream>
#include <stdexcept>

/*
Test parameter parsing, molecule parsing, and
molecule graph construction
*/
int main(int argc, char** argv)
{
    // adapting from the homework
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << "path/to/molecule" << " path/to/params" << std::endl;
        return EXIT_FAILURE;
    }

    // parse molecule (first argument in CLI)
    MoleculeGraph mol = read_xyz(argv[1]);
    assign_amber_types(mol);
    std::cout << "Molecule loaded!" << std::endl << std::endl;
    std::cout << "atoms: " << mol.num_atoms() << std::endl;
    std::cout << "bonds: " << mol.num_bonds() << std::endl << std::endl;

    // parse parameters (second of two arguments in CLI)]
    ForceField ff = load_params(argv[2]);
    std::cout << "Forcefield parameters loaded!" << std::endl;
    std::cout << "# of bond params: " << ff.bonds.size() << std::endl;
    std::cout << "# of angle params: " << ff.angles.size() << std::endl;
    std::cout << "# of dihedral params: " << ff.dihedrals.size() << std::endl;
    std::cout << "# of vdw params: " << ff.vdw.size() << std::endl << std::endl;

    for (int i = 0; i < mol.num_atoms(); i++)
    {
        const Atom& a = mol.atoms[i];

        std::cout << i << "  " << a.element << "  (" << a.position(0)
                  << ", " << a.position(1) << ", " << a.position(2) << ")" 
                  << " | AMBER type: " << a.amber_type << " | "
                  << "neighbors={";

        std::vector<int> neighbors = mol.neighbors(i);

        for (size_t k = 0; k < neighbors.size(); k++)
        {
            std::cout << neighbors[k];
            if (k + 1 < neighbors.size()) std::cout << ",";
        }
        std::cout << "}" << std::endl;
    }
    return 0;
}