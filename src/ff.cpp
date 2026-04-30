
/*
Loop through atoms -
1. for bonds term, loop through all bonded atoms, from the MoleculeGraph struct
2. for the angles term,
3. for the dihedrals
4. for the vdw term, all non bonded atoms
5. for the h-bonds term, all h-bonded atoms, identified by atomic number of the atom

within the loops, calculate the individual terms,
use the read_params files for getting values of the parameters corresponding to each atom type and then structural info from the ff_params.hpp
formula used are the original amber force field terms


Imp Question: where do A, B, C, D terms come from for the VDW and HBond terms?? 
*/

#include "ff.hpp"
#include <cmath>
#include <armadillo>


// later move all constants to a constants namespace in a header file
static const double DEG2RAD = M_PI / 180.0;

// E = K_r * (r - r_eq)^2
double bond_streching_energy(double bond_length, double K_r, double r_eq)
{
    double dr = bond_length - r_eq;
    return K_r * dr * dr;
}

// include only bonded atoms
double calculate_bonds_term(const MoleculeGraph& mol, const ForceField& ff)
{
    double total_be = 0.0;
    for (int i = 0; i < mol.num_atoms(); ++i)
    {
        for (int j : mol.neighbors(i))
        {
            if (j <= i) continue; // each bond once
            const std::string& ti = mol.atoms[i].amber_type;
            const std::string& tj = mol.atoms[j].amber_type;
            double r = dist(mol.atoms[i].position, mol.atoms[j].position);
            auto it = ff.bonds.find({ti, tj});
            if (it == ff.bonds.end()) it = ff.bonds.find({tj, ti});
            if (it == ff.bonds.end()) continue;
            total += bond_streching_energy(r, it->second.K_r, it->second.r_eq);
        }
    }
    return total_be;
}

// E = K_theta * (theta_rad - theta_eq_rad)^2
double bond_rotation_energy(double angle_deg, double K_theta, double theta_eq_deg)
{
    double dtheta = (angle_deg - theta_eq_deg) * DEG2RAD;
    return K_theta * dtheta * dtheta;
}

// only adjacent bonds
double calculate_angles_term(const MoleculeGraph& mol, const ForceField& ff)
{
    double total = 0.0;
    for (int j = 0; j < mol.num_atoms(); ++j)
    {
        const auto& nbrs = mol.neighbors(j);
        for (int a = 0; a < (int)nbrs.size(); ++a)
        {
            for (int b = a + 1; b < (int)nbrs.size(); ++b)
            {
                int i = nbrs[a], k = nbrs[b];
                const std::string& ti = mol.atoms[i].amber_type;
                const std::string& tj = mol.atoms[j].amber_type;
                const std::string& tk = mol.atoms[k].amber_type;
                double theta = angle_deg(mol.atoms[i].position,
                                         mol.atoms[j].position,
                                         mol.atoms[k].position);
                auto it = ff.angles.find({ti, tj, tk});
                if (it == ff.angles.end()) it = ff.angles.find({tk, tj, ti});
                if (it == ff.angles.end()) continue;
                total += bond_rotation_energy(theta, it->second.K_theta, it->second.theta_eq);
            }
        }
    }
    return total;
}

// E = sum_k (Vn/divider) * (1 + cos(n * phi - gamma)), angles in degrees
double bond_torsion_energy(double phi_deg, const std::vector<DihedralParams>& terms)
{
    double e = 0.0;
    double phi_rad = phi_deg * DEG2RAD;
    for (const auto& t : terms)
        e += (t.Vn / t.divider) * (1.0 + std::cos(t.n * phi_rad - t.gamma * DEG2RAD));
    return e;
}

double calculate_dihedrals_term(const MoleculeGraph& mol, const ForceField& ff)
{
    double total = 0.0;
    for (int j = 0; j < mol.num_atoms(); ++j)
    {
        for (int k : mol.neighbors(j))
        {
            if (k <= j) continue; // each central bond once
            for (int i : mol.neighbors(j))
            {
                if (i == k) continue;
                for (int l : mol.neighbors(k))
                {
                    if (l == j || l == i) continue;
                    const std::string& ti = mol.atoms[i].amber_type;
                    const std::string& tj = mol.atoms[j].amber_type;
                    const std::string& tk = mol.atoms[k].amber_type;
                    const std::string& tl = mol.atoms[l].amber_type;
                    double phi = dihedral_deg(mol.atoms[i].position, mol.atoms[j].position,
                                              mol.atoms[k].position, mol.atoms[l].position);
                    auto it = ff.dihedrals.find({ti, tj, tk, tl});
                    if (it == ff.dihedrals.end()) it = ff.dihedrals.find({tl, tk, tj, ti});
                    if (it == ff.dihedrals.end()) continue;
                    total += bond_torsion_energy(phi, it->second);
                }
            }
        }
    }
    return total;
}

// E = A/r^12 - B/r^6
// A = epsilon_ij * R_ij^12, B = epsilon_ij * R_ij^6
double vdw_energy(double r, double A, double B)
{
    double r6 = std::pow(r, 6);
    return A / (r6 * r6) - B / r6;
}

double calculate_electrostatic_vdw_term(const MoleculeGraph& mol, const ForceField& ff)
{
    double total = 0.0;
    int n = mol.num_atoms();
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            bool bonded = false;
            for (int nb : mol.neighbors(i))
                if (nb == j) { bonded = true; break; }
            if (bonded) continue;

            const std::string& ti = mol.atoms[i].amber_type;
            const std::string& tj = mol.atoms[j].amber_type;
            auto it_i = ff.vdw.find(ti);
            auto it_j = ff.vdw.find(tj);
            if (it_i == ff.vdw.end() || it_j == ff.vdw.end()) continue;

            double R_star_ij = it_i->second.R_star + it_j->second.R_star;
            double eps_ij    = std::sqrt(it_i->second.epsilon * it_j->second.epsilon);
            double r         = dist(mol.atoms[i].position, mol.atoms[j].position);
            total += vdw_energy(r, R_star_ij, eps_ij);
        }
    }
    return total;
}

// Weak explicit H-bond term
// E = C/r^12 - D/r^10
double h_bond_energy(double r, double C, double D)
{
    double r10 = std::pow(r, 10);
    return C / (r10 * r * r) - D / r10;
}


double calculate_hbond_term(const MoleculeGraph& mol, const ForceField& ff)
{
    double total = 0.0;
    int n = mol.num_atoms();
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            bool i_is_H = (mol.atoms[i].element == "H");
            bool j_is_H = (mol.atoms[j].element == "H");
            if (!i_is_H && !j_is_H) continue;

            bool bonded = false;
            for (int nb : mol.neighbors(i))
                if (nb == j) { bonded = true; break; }
            if (bonded) continue;

            const std::string& ti = mol.atoms[i].amber_type;
            const std::string& tj = mol.atoms[j].amber_type;
            auto it_i = ff.vdw.find(ti);
            auto it_j = ff.vdw.find(tj);
            if (it_i == ff.vdw.end() || it_j == ff.vdw.end()) continue;

            double R_star_ij = it_i->second.R_star + it_j->second.R_star;
            double eps_ij    = std::sqrt(it_i->second.epsilon * it_j->second.epsilon);
            double C = eps_ij * std::pow(R_star_ij, 12);
            double D = eps_ij * std::pow(R_star_ij, 10);
            double r = dist(mol.atoms[i].position, mol.atoms[j].position);
            total += h_bond_energy(r, C, D);
        }
    }
    return total;
}
