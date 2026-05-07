#pragma once
#include "molecule.hpp"
#include <armadillo>
#include <cmath>

inline const double DEG2RAD = M_PI / 180.0;

// arma::vec overloads (core implementations)
double dist(const arma::vec& a, const arma::vec& b);
double angle_deg(const arma::vec& pi, const arma::vec& pj, const arma::vec& pk);
double dihedral_deg(const arma::vec& pi, const arma::vec& pj, const arma::vec& pk, const arma::vec& pl);

// Atom overloads (convenience wrappers)
double dist(const Atom& a, const Atom& b);
double angle_deg(const Atom& ai, const Atom& aj, const Atom& ak);
double dihedral_deg(const Atom& ai, const Atom& aj, const Atom& ak, const Atom& al);