#pragma once

#include <cstddef>
#include <vector>

#include "landscape.hpp"
#include "matrix.hpp"

struct Fire {
  unsigned int width;
  unsigned int height;

  Matrix<bool> burned_layer;

  std::vector<std::pair<unsigned int, unsigned int>> burned_ids;

  // Positions in burned_ids where a new step starts, empty if the fire was not simulated
  std::vector<unsigned int> burned_ids_steps;

  bool operator==(const Fire& other) const {
    return width == other.width && height == other.height &&
           burned_layer == other.burned_layer && burned_ids == other.burned_ids;
  }
};

Fire read_fire(unsigned int width, unsigned int height, std::string filename);

/* Type and function useful for comparing fires */

struct FireStats {
  unsigned int counts_veg_matorral;
  unsigned int counts_veg_subalpine;
  unsigned int counts_veg_wet;
  unsigned int counts_veg_dry;
};

FireStats get_fire_stats(const Fire& fire, const Landscape& landscape);
