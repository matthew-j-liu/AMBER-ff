#pragma once
#include "ff_params.hpp"
#include <string>


// Parses parameter file and returns a populated ForceField 
ForceField load_params(const std::string& filepath); 