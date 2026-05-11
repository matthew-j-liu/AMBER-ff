#include "molecule.hpp"
#include "ff_params.hpp"
#include "parsing.hpp"
#include "ff.hpp"
#include "geom_util.hpp"
#include <armadillo>
#include <cmath>
#include <iostream>


const char* xyz_path    = "input_molecules_copy/alkenes/But-1-ene.xyz";
const char* params_path = "params/selected_atoms.dat";
const int i = 3, j = 2, k = 1, l = 0;
const double angle_min = -180.0, angle_max = 180.0, step = 1.0;

// Atoms being rotated 
const int moving_atoms[] = {0, 1, 5, 6, 7, 8, 11};


int main()
{
    ForceField ff = load_params(params_path);
    MoleculeGraph mol = read_xyz(xyz_path);

    std::vector<arma::vec> ref(mol.num_atoms());
    for (size_t a = 0; a < mol.num_atoms(); a++) 
        ref[a] = mol.atoms[a].position;

    const double phi0 = dihedral_angle(ref[i], ref[j], ref[k], ref[l]);
    const arma::vec origin = ref[j];
    const arma::vec axis   = arma::normalise(ref[k] - ref[j]);

    std::cout << "angle_deg,bonds,angles,dihedrals,vdw,total\n";
    std::cout.setf(std::ios::fixed);
    std::cout.precision(6);

    for (double angle = angle_min; angle <= angle_max + step * 0.5; angle += step)
    {
        const double delta = (angle - phi0) * DEG_to_RAD;
        const double c = std::cos(delta), s = std::sin(delta);

        // rotation of each moving atom about (origin, axis).
        for (int idx : moving_atoms)
        {
            const arma::vec v = ref[idx] - origin;
            mol.atoms[idx].position = origin + v * c
                                    + arma::cross(axis, v) * s
                                    + axis * (arma::dot(axis, v) * (1.0 - c));
        }

        const double e_bonds = calculate_bonds_term(mol, ff);
        const double e_angle = calculate_angles_term(mol, ff);
        const double e_dihedral = calculate_dihedrals_term(mol, ff);
        const double e_vdw = calculate_vdw_term(mol, ff);

        std::cout << angle << ',' << e_bonds << ',' << e_angle << ',' << e_dihedral << ','
                  << e_vdw << ',' << (e_bonds + e_angle + e_dihedral + e_vdw) << '\n';
    }

    return 0;
}
