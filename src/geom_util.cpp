// Utility functions related to the structure/ geometry of the molecule

#include "geom_util.hpp"
#include "ff_params.hpp"
#include <algorithm>

// Calculate the bond distance, bond angle and dihedral angles
// euclidean distance between two atom positions
double bond_length(const arma::vec& a, const arma::vec& b)
{
    return arma::norm(b - a);
}

// bond angle at j in the triplet i-j-k, returned in degrees
double angle_deg(const arma::vec& pi, const arma::vec& pj, const arma::vec& pk)
{
    arma::vec v1 = pi - pj;
    arma::vec v2 = pk - pj;
    double c = arma::dot(v1, v2) / (arma::norm(v1) * arma::norm(v2));
    c = std::max(-1.0, std::min(1.0, c));
    return std::acos(c) / DEG_to_RAD;
}

// dihedral angle for i-j-k-l, returned in degrees
double dihedral_angle(const arma::vec& pi, const arma::vec& pj,
                            const arma::vec& pk, const arma::vec& pl)
{
    arma::vec b1 = pj - pi;
    arma::vec b2 = pk - pj;
    arma::vec b3 = pl - pk;
    arma::vec n1 = arma::cross(b1, b2);
    arma::vec n2 = arma::cross(b2, b3);
    double c = arma::dot(n1, n2) / (arma::norm(n1) * arma::norm(n2));
    c = std::max(-1.0, std::min(1.0, c));
    double phi = std::acos(c) / DEG_to_RAD;
    if (arma::dot(arma::cross(n1, n2), b2) < 0.0) phi = -phi;
    return phi;
}

// Atom overloads — delegate to the arma::vec versions above
double bond_length(const Atom& a, const Atom& b)
{
    return bond_length(a.position, b.position);
}

double angle_deg(const Atom& ai, const Atom& aj, const Atom& ak)
{
    return angle_deg(ai.position, aj.position, ak.position);
}

double dihedral_angle(const Atom& ai, const Atom& aj, const Atom& ak, const Atom& al)
{
    return dihedral_angle(ai.position, aj.position, ak.position, al.position);
}

