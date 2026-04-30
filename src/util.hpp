#pragma once
#include <armadillo>
#include <cmath>

// degrees to radians
inline const double DEG2RAD = M_PI / 180.0; 

// functions in util.cpp
double dist(const arma::vec& a, const arma::vec& b);
double angle_deg(const arma::vec& pi, const arma::vec& pj, const arma::vec& pk);
double dihedral_deg(const arma::vec& pi, const arma::vec& pj, const arma::vec& pk, const arma::vec& pl);