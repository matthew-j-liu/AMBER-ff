#include "geom_util.hpp"
#include "ff.hpp"
#include <cmath>
#include <set>

// E = K_r * (r - r_eq)^2
double bond_streching_energy(double bond_length, double K_r, double r_eq)
{
    double dr = bond_length - r_eq;
    return K_r * dr * dr;
}

double calculate_bonds_term(const MoleculeGraph& mol, const ForceField& ff)
{
    double total = 0.0;
    for (const auto& bond : mol.find_all_bonds())
    {
        size_t i = bond[0], j = bond[1];
        const std::string& ti = mol.atoms[i].amber_type;
        const std::string& tj = mol.atoms[j].amber_type;
        double r = dist(mol.atoms[i], mol.atoms[j]);
        auto it = ff.bonds.find({ti, tj});
        if (it == ff.bonds.end()) it = ff.bonds.find({tj, ti});
        if (it == ff.bonds.end()) continue;
        total += bond_streching_energy(r, it->second.K_r, it->second.r_eq);
    }
    return total;
}

// E = K_theta * (theta_rad - theta_eq_rad)^2
double bond_rotation_energy(double angle_deg, double K_theta, double theta_eq_deg)
{
    double dtheta = (angle_deg - theta_eq_deg) * DEG2RAD;
    return K_theta * dtheta * dtheta;
}

double calculate_angles_term(const MoleculeGraph& mol, const ForceField& ff)
{
    double total = 0.0;
    for (const auto& triplet : mol.find_all_bond_angle_triplets())
    {
        size_t i = triplet[0], j = triplet[1], k = triplet[2];
        const std::string& ti = mol.atoms[i].amber_type;
        const std::string& tj = mol.atoms[j].amber_type;
        const std::string& tk = mol.atoms[k].amber_type;
        double theta = angle_deg(mol.atoms[i], mol.atoms[j], mol.atoms[k]);
        auto it = ff.angles.find({ti, tj, tk});
        if (it == ff.angles.end()) it = ff.angles.find({tk, tj, ti});
        if (it == ff.angles.end()) continue;
        total += bond_rotation_energy(theta, it->second.K_theta, it->second.theta_eq);
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
    for (const auto& quad : mol.find_all_torsion_quadruplets())
    {
        size_t i = quad[0], j = quad[1], k = quad[2], l = quad[3];
        const std::string& ti = mol.atoms[i].amber_type;
        const std::string& tj = mol.atoms[j].amber_type;
        const std::string& tk = mol.atoms[k].amber_type;
        const std::string& tl = mol.atoms[l].amber_type;
        double phi = dihedral_deg(mol.atoms[i], mol.atoms[j], mol.atoms[k], mol.atoms[l]);
        auto it = ff.dihedrals.find({ti, tj, tk, tl});
        if (it == ff.dihedrals.end()) it = ff.dihedrals.find({tl, tk, tj, ti});
        if (it == ff.dihedrals.end()) continue;
        total += bond_torsion_energy(phi, it->second);
    }
    return total;
}

// E = A/r^12 - B/r^6
double vdw_energy(double r, double A, double B)
{
    double r6 = std::pow(r, 6);
    return A / (r6 * r6) - B / r6;
}

double calculate_electrostatic_vdw_term(const MoleculeGraph& mol, const ForceField& ff)
{
    double total = 0.0;
    size_t n = mol.num_atoms();

    // AMBER excludes 1-2 and 1-3 pairs; 1-4 pairs are scaled by 0.5 (GAFF scee/scnb)
    std::set<std::pair<size_t, size_t>> excluded;
    std::set<std::pair<size_t, size_t>> scaled_14;

    for (const auto& bond : mol.find_all_bonds())
        excluded.insert({bond[0], bond[1]});

    for (const auto& triplet : mol.find_all_bond_angle_triplets())
    {
        size_t i = triplet[0], k = triplet[2];
        excluded.insert({std::min(i, k), std::max(i, k)});
    }

    for (const auto& quad : mol.find_all_torsion_quadruplets())
    {
        size_t i = quad[0], l = quad[3];
        auto key = std::make_pair(std::min(i, l), std::max(i, l));
        if (!excluded.count(key)) scaled_14.insert(key);
    }

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = i + 1; j < n; ++j)
        {
            auto key = std::make_pair(i, j);
            if (excluded.count(key)) continue;

            const std::string& ti = mol.atoms[i].amber_type;
            const std::string& tj = mol.atoms[j].amber_type;
            auto it_i = ff.vdw.find(ti);
            auto it_j = ff.vdw.find(tj);
            if (it_i == ff.vdw.end() || it_j == ff.vdw.end()) continue;

            double R_star_ij = it_i->second.R_star + it_j->second.R_star;
            double eps_ij    = std::sqrt(it_i->second.epsilon * it_j->second.epsilon);
            double A = eps_ij * std::pow(R_star_ij, 12);
            double B = 2.0 * eps_ij * std::pow(R_star_ij, 6);
            double r = dist(mol.atoms[i], mol.atoms[j]);

            double e = vdw_energy(r, A, B);
            if (scaled_14.count(key)) e *= 0.5; // GAFF scnb = 2 → factor 1/2
            total += e;
        }
    }
    return total;
}
