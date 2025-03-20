#pragma once

#include <string>
#include <vector>

typedef std::vector<std::pair<unsigned int, unsigned int>> IgnitionCells;

IgnitionCells read_ignition_cells(std::string filename);
