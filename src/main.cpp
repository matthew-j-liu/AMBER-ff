#include "molecule.hpp"
#include "read_xyz.hpp"
#include "ff_params.hpp"
#include "read_params.hpp"
#include "ff.hpp"
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
    std::cout << "Molecule loaded!" << std::endl;
    std::cout << "atoms: " << mol.num_atoms() << std::endl;
    std::cout << "bonds: " << mol.num_bonds() << std::endl << std::endl;

    // parse parameters (second of two arguments in CLI)]
    ForceField ff = load_params(argv[2]);
    std::cout << "Forcefield parameters loaded!" << std::endl;
    std::cout << "# of bond params: " << ff.bonds.size() << std::endl;
    std::cout << "# of angle params: " << ff.angles.size() << std::endl;
    std::cout << "# of dihedral params: " << ff.dihedrals.size() << std::endl;
    std::cout << "# of vdw params: " << ff.vdw.size() << std::endl;

    for (int i = 0; i < mol.num_atoms(); i++)
    {
        const Atom& a = mol.atoms[i];

        std::cout << i << "  " << a.element << "  (" << a.position(0)
                  << ", " << a.position(1) << ", " << a.position(2) << ") | "
                  << "  neighbors={";

        std::vector<int> neighbors = mol.neighbors(i);

        for (size_t k = 0; k < neighbors.size(); k++)
        {
            std::cout << neighbors[k];
            if (k + 1 < neighbors.size()) std::cout << ",";
        }
        std::cout << "}" << std::endl;
    }

    // calculate energy 
    double bonds_term = calculate_bonds_term(mol, ff);
    double angles_term = calculate_angles_term(mol, ff); 
    double dihedrals_term = calculate_dihedrals_term(mol, ff); 
    double electrostatic_vdw_term = calculate_electrostatic_vdw_term(mol, ff); 

    double total_energy = bonds_term + angles_term + dihedrals_term + electrostatic_vdw_term; 

    std::cout << "Bond stretching energy = " << bonds_term << std::endl;
    std::cout << "Bond Rotation energy = " << angles_term << std::endl;
    std::cout << "Torsion energy = " << dihedrals_term << std::endl;
    std::cout << "Van der waals energy = " << electrostatic_vdw_term << std::endl;

    std::cout << "Total energy = " << total_energy << std::endl;
    return 0;
}