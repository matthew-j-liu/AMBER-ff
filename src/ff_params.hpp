#pragma once 
#include <vector>
#include <map>
#include <tuple>
#include <string>

struct BondParams
{ double K_r; double r_eq; };


struct AngleParams
{ double K_theta; double theta_eq; };


struct DihedralParams
{ int divider; double Vn; double gamma; int n; };


struct VdwParams
{ double R_star; double epsilon; };

/*
Let's reason through how many atoms are needed to encode a bond, angle, and dihedral 

bond: 2 atoms (e.g., C-H bond)
angle: 3 atoms (e.g., H-O-H angle is 104.5 degrees)
dihedral: 4 atoms (2 atoms form axis of rotation. Atoms attached to each define the full rotation)"

A logical structure would therefore use the atoms involved in each of these properties (e.g., a bond)
as a key to the parameter values. 
*/

using BondKey = std::tuple<std::string, std::string>; // 2 atoms define a bond 
using AngleKey = std::tuple<std::string, std::string, std::string>; // 3 atoms define a bond angle 
using DihedralKey = std::tuple<std::string, std::string, std::string, std::string>; // 4 atoms define a dihedral angle 

/*
We now construct a dictionary by creating key-value pairs. Note that Dihedral Key will
need to map to a vector of DihedralParams because the dihedral contribution to energy 
is a Fourier series that may contain multiple terms. For example, see Table III 
in the following paper, which tabulates values for V1 through V6 for n-butane:
https://pubs.acs.org/doi/10.1021/j100463a018                 
*/

struct ForceField
{
    std::map<BondKey, BondParams> bonds;
    std::map<AngleKey, AngleParams> angles;
    std::map<DihedralKey, std::vector<DihedralParams>> dihedrals; 
    std::map<std::string, VdwParams> vdw;  
};
