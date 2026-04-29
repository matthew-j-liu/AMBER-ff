#include "read_xyz.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>

Molecule read_xyz(const std::string& filepath)
{
    std::ifstream f(filepath);  

    if (!f.is_open())     
    {
        throw std::runtime_error("Cannot open: " + filepath);
    } 

    Molecule molecule;

    // adapt xyz parsing function from homeworks 

    return molecule; 
}