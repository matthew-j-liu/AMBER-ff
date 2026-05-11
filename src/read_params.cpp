#include "parsing.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>


/*
define 4 helper functions for parsing - one for each parameter struct

the paradigm for each parsing function is as follows

inputs
-------
    -iss: std::istringstream with one line of data
    -ff: ForceField to populate

logic
------
    -declare empty variables. Strings for atom-types and our custom `struct` for parameters
    -parse the text, with whitespace as separator. note: hard-coded the number of >> based on data file
    -store key-value pair, where the key is a tuple of strings (or just one, for vdw) and the value is the `struct`
*/

void parse_bond_row(std::istringstream& iss, ForceField& ff)
{
    std::string ti, tj;
    BondParams p;

    if (!(iss >> ti >> tj >> p.K_r >> p.r_eq))
    {
        throw std::runtime_error("Malformed bond row");
    }

    ff.bonds[{ti, tj}] = p;
}

void parse_angle_row(std::istringstream& iss, ForceField& ff)
{
    std::string ti, tj, tk;
    AngleParams p;

    if (!(iss >> ti >> tj >> tk >> p.K_theta >> p.theta_eq))
    {
        throw std::runtime_error("Malformed angle row");
    }

    ff.angles[{ti, tj, tk}] = p;
}

void parse_dihedral_row(std::istringstream& iss, ForceField& ff)
{
    std::string ti, tj, tk, tl;
    DihedralParams p;

    if (!(iss >> ti >> tj >> tk >> tl >> p.divider >> p.Vn >> p.gamma >> p.n))
    {
        throw std::runtime_error("Malformed dihedral row");
    }

    ff.dihedrals[{ti, tj, tk, tl}].push_back(p);
}

void parse_vdw_row(std::istringstream& iss, ForceField& ff)
{
    std::string type;
    VdwParams p;

    if (!(iss >> type >> p.R_star >> p.epsilon))
    {
        throw std::runtime_error("Malformed vdw row");
    }

    ff.vdw[type] = p;
}

// make use of all 4 parsers
ForceField load_params(const std::string& filepath)
{
    std::ifstream f(filepath);

    if (!f.is_open())
    {
        throw std::runtime_error("Cannot open: " + filepath);
    }

    ForceField ff;
    std::string section = "";
    std::string line;

    while (std::getline(f, line))
    {
        // skip blank lines and comments
        if (line.empty() || line[0] == '#') { continue; }

        // identify which parsing function to use
        if (line[0] == '[')
        {
            if (line == "[bonds]") { section = "bonds"; }
            else if (line == "[angles]") {section = "angles"; }
            else if (line == "[dihedrals]") { section = "dihedrals"; }
            else if (line == "[vdw]") { section = "vdw"; }
            else { throw std::runtime_error("Unknown section: " + line); }
            continue;
        }

        // deploy the parsing function
        std::istringstream iss(line);
        if (section == "bonds") { parse_bond_row(iss, ff); }
        else if (section == "angles") { parse_angle_row(iss, ff); }
        else if (section == "dihedrals") { parse_dihedral_row(iss, ff); }
        else if (section == "vdw") { parse_vdw_row(iss, ff); }
    }

    return ff;
}