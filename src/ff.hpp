#pragma once
#include "molecule.hpp"
#include "ff_params.hpp"
#include <vector>

double calculate_bonds_term(const MoleculeGraph& mol, const ForceField& ff);
double calculate_angles_term(const MoleculeGraph& mol, const ForceField& ff);
double calculate_dihedrals_term(const MoleculeGraph& mol, const ForceField& ff);
double calculate_electrostatic_vdw_term(const MoleculeGraph& mol, const ForceField& ff);

double bond_streching_energy(double bond_length, double K_r, double r_eq);
double bond_rotation_energy(double angle_deg, double K_theta, double theta_eq_deg);
double bond_torsion_energy(double phi_deg, const std::vector<DihedralParams>& terms);
double vdw_energy(double r, double A, double B);
