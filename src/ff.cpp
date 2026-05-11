// Calculates the energy terms of the Amber Force Field

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

// E = K_theta * (theta_rad - theta_eq_rad)^2
double bond_rotation_energy(double angle_deg, double K_theta, double theta_eq_deg)
{
    double dtheta = (angle_deg - theta_eq_deg) * DEG_to_RAD;
    return K_theta * dtheta * dtheta;
}

// E = sum_k (Vn/divider) * (1 + cos(n * phi - gamma)), angles in degrees
double bond_torsion_energy(double phi_deg, const std::vector<DihedralParams>& terms)
{
    double e = 0.0;
    double phi_rad = phi_deg * DEG_to_RAD;
    for (const auto& t : terms)
        e += (t.Vn / t.divider) * (1.0 + std::cos(t.n * phi_rad - t.gamma * DEG_to_RAD));
    return e;
}

// E = A/r^12 - B/r^6
double vdw_energy(double r, double A, double B)
{
    double r6 = std::pow(r, 6);
    return A / (r6 * r6) - B / r6;
}

// summation over the full molecule 

// begin with direct bonds 
double calculate_bonds_term(const MoleculeGraph& mol, const ForceField& ff)
{
    double total = 0.0;
    // find_all_bonds() returns as std::vector<std::array<size_t, 2>> 
    for (const auto& bond : mol.find_all_bonds())
    {
        // get indices and amber types of all atoms
        size_t i = bond[0], j = bond[1];  
        const std::string& type_i = mol.atoms[i].amber_type;
        const std::string& type_j = mol.atoms[j].amber_type;

        // find bond length 
        double r = bond_length(mol.atoms[i], mol.atoms[j]);

        // read params using BondKey = std::tuple<std::string, std::string> 
        // std::map.find(key) returns an iterator to the element as std::pair<const KeyType, ParamsType>
        auto bond_params = ff.bonds.find({type_i, type_j});
        if (bond_params == ff.bonds.end()) 
        {
            // reverse search 
            bond_params = ff.bonds.find({type_j, type_i});
        }
        // give up and continue if still not found (i.e. +energy = 0)
        if (bond_params == ff.bonds.end()) continue;

        total += bond_streching_energy(r, bond_params->second.K_r, bond_params->second.r_eq);
    }
    return total;
}

// loop over all triplets 
double calculate_angles_term(const MoleculeGraph& mol, const ForceField& ff)
{
    double total = 0.0;
    // find_all_bond_angle_triplets() returns as std::vector<std::array<size_t, 3>> 
    for (const auto& triplet : mol.find_all_bond_angle_triplets())
    {
        // get indices and amber types of all atoms
        size_t i = triplet[0], j = triplet[1], k = triplet[2];
        const std::string& type_i = mol.atoms[i].amber_type;
        const std::string& type_j = mol.atoms[j].amber_type;
        const std::string& type_k = mol.atoms[k].amber_type;

        // find bond angle 
        double theta = angle_deg(mol.atoms[i], mol.atoms[j], mol.atoms[k]);

        // read params using AngleKey = std::tuple<std::string, std::string, std::string>
        auto angle_params = ff.angles.find({type_i, type_j, type_k});
        if (angle_params == ff.angles.end()) 
        {
            // reverse search 
            angle_params = ff.angles.find({type_k, type_j, type_i});
        }
        // give up and continue if still not found (i.e. +energy = 0)
        if (angle_params == ff.angles.end()) continue;

        total += bond_rotation_energy(theta, angle_params->second.K_theta, angle_params->second.theta_eq);
    }
    return total;
}

// loops over all quadruplets
double calculate_dihedrals_term(const MoleculeGraph& mol, const ForceField& ff)
{
    double total = 0.0;
    // find_all_torsion_quadruplets() returns as std::vector<std::array<size_t, 4>> 
    for (const auto& quad : mol.find_all_torsion_quadruplets())
    {
        // get indices and amber types of all atoms
        size_t i = quad[0], j = quad[1], k = quad[2], l = quad[3];
        const std::string& type_i = mol.atoms[i].amber_type;
        const std::string& type_j = mol.atoms[j].amber_type;
        const std::string& type_k = mol.atoms[k].amber_type;
        const std::string& type_l = mol.atoms[l].amber_type;
        const std::string& type_wildcard = "X"; 

        // find dihedral angle 
        double phi = dihedral_angle(mol.atoms[i], mol.atoms[j], mol.atoms[k], mol.atoms[l]);

        // read params using DihedralKey = std::tuple<std::string, std::string, std::string, std::string>
        // note that many dihedral params use X as wild-card

        auto dihedral_params = ff.dihedrals.find({type_i, type_j, type_k, type_l});
        // try forward search 
        if (dihedral_params == ff.dihedrals.end()) 
        {
            // try reverse search
            dihedral_params = ff.dihedrals.find({type_l, type_k, type_j, type_i});
            if (dihedral_params == ff.dihedrals.end()) 
            {
                // try forward with wildcard
                dihedral_params = ff.dihedrals.find({type_wildcard, type_j, type_k, type_wildcard});
                if (dihedral_params == ff.dihedrals.end()) 
                {
                    // try reverse with wildcard
                    dihedral_params = ff.dihedrals.find({type_wildcard, type_k, type_j, type_wildcard});
                }
            }
            // give up and continue if still not found (i.e. +energy = 0)
            if (dihedral_params == ff.dihedrals.end()) continue; 
        }
        total += bond_torsion_energy(phi, dihedral_params->second);
    }
    return total;
}

// loops over all atom pairs first, and excludes those counted by other terms
double calculate_vdw_term(const MoleculeGraph& mol, const ForceField& ff)
{
    double total = 0.0;
    size_t n = mol.num_atoms();

    // Only 1-4 pairs considered and scaled by 0.5 
    std::set<std::pair<size_t, size_t>> excluded;
    std::set<std::pair<size_t, size_t>> scaled_14;

    // 1-2 pairs already accounted by bond stretching, add to excluded (to prevent double counting)
    for (const auto& bond : mol.find_all_bonds())
        excluded.insert({bond[0], bond[1]});

    // 1-3 pairs already accounted by bond angles, add to excluded (to prevent double counting)
    for (const auto& triplet : mol.find_all_bond_angle_triplets())
    {
        // indices of i-j-k in triplet
        size_t i = triplet[0], k = triplet[2];
        excluded.insert({std::min(i, k), std::max(i, k)});
    }

    // 1-4 pairs to count
    for (const auto& quad : mol.find_all_torsion_quadruplets())
    {
        // indices of i-j-k-l in quad
        size_t i = quad[0], l = quad[3];
        auto vdw_atom_pair = std::make_pair(std::min(i, l), std::max(i, l));
        // pair addded to 
        if (!excluded.count(vdw_atom_pair)) scaled_14.insert(vdw_atom_pair);
    }

    // loop through all atoms
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = i + 1; j < n; ++j)
        {
            auto atom_pair = std::make_pair(i, j);
            if (excluded.count(atom_pair)) continue;

            const std::string& type_i = mol.atoms[i].amber_type;
            const std::string& type_j = mol.atoms[j].amber_type;

            auto vdw_params_i = ff.vdw.find(type_i);
            auto vdw_params_j = ff.vdw.find(type_j);
            if (vdw_params_i == ff.vdw.end() || vdw_params_j == ff.vdw.end()) continue;

            // R*_ij = R*_i + R*_j, E_ij = sqrt(E_i * E_j)
            double R_star_ij = vdw_params_i->second.R_star + vdw_params_j->second.R_star;
            double eps_ij = std::sqrt(vdw_params_i->second.epsilon * vdw_params_j->second.epsilon);

            double A = eps_ij * std::pow(R_star_ij, 12);
            double B = 2.0 * eps_ij * std::pow(R_star_ij, 6);

            double r = bond_length(mol.atoms[i], mol.atoms[j]);

            double e = vdw_energy(r, A, B);
            if (scaled_14.count(atom_pair)) e *= 0.5; 
            total += e;
        }
    }
    return total;
}
