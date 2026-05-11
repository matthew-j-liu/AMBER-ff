#pragma once
#include "molecule.hpp"
#include <armadillo>
#include <cmath>

inline const double DEG_to_RAD = M_PI / 180.0;

// arma::vec overloads 
double bond_length(const arma::vec& a, const arma::vec& b);
double angle_deg(const arma::vec& pi, const arma::vec& pj, const arma::vec& pk);
double dihedral_angle(const arma::vec& pi, const arma::vec& pj, const arma::vec& pk, const arma::vec& pl);

// Atom overloads 
double bond_length(const Atom& a, const Atom& b);
double angle_deg(const Atom& ai, const Atom& aj, const Atom& ak);
double dihedral_angle(const Atom& ai, const Atom& aj, const Atom& ak, const Atom& al);