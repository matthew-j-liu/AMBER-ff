#include <sstream> 
#include <fstream> 

/* 
Parameter lookups for all 4 functions from the dat file
*/ 

double bond_term_params(std::string atom_type)
{
    
}


double angles_term_params(std::string atom_type1, std::string atom_type2)
{
    
}


double dihedrals_term_params(std::string atom_type1, std::string atom_type2, std::string atom_type3)
{
    
}


double electrostatic_vdw_term_params(std::string atom_type1, std::string atom_type2)
{
    
}

// Calculate the bond distance, bond angle and dihedral angles
// euclidean distance between two atom positions
static double dist(const arma::vec& a, const arma::vec& b)
{
    return arma::norm(b - a);
}

// bond angle at j in the triplet i-j-k, returned in degrees
static double angle_deg(const arma::vec& pi, const arma::vec& pj, const arma::vec& pk)
{
    arma::vec v1 = pi - pj;
    arma::vec v2 = pk - pj;
    double c = arma::dot(v1, v2) / (arma::norm(v1) * arma::norm(v2));
    c = std::max(-1.0, std::min(1.0, c));
    return std::acos(c) / DEG2RAD;
}

// dihedral angle for i-j-k-l, returned in degrees
static double dihedral_deg(const arma::vec& pi, const arma::vec& pj,
                            const arma::vec& pk, const arma::vec& pl)
{
    arma::vec b1 = pj - pi;
    arma::vec b2 = pk - pj;
    arma::vec b3 = pl - pk;
    arma::vec n1 = arma::cross(b1, b2);
    arma::vec n2 = arma::cross(b2, b3);
    double c = arma::dot(n1, n2) / (arma::norm(n1) * arma::norm(n2));
    c = std::max(-1.0, std::min(1.0, c));
    double phi = std::acos(c) / DEG2RAD;
    if (arma::dot(arma::cross(n1, n2), b2) < 0.0) phi = -phi;
    return phi;
}

