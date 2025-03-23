#include "many_simulations.hpp"

#include <cmath>

Matrix<unsigned int> burned_amounts_per_cell(
    const Landscape& landscape, const std::vector<std::pair<unsigned int, unsigned int>>& ignition_cells,
    SimulationParams params, float distance, float elevation_mean, float elevation_sd,
    float upper_limit, unsigned int n_replicates
) {

  Matrix<unsigned int> burned_amounts(landscape.width, landscape.height);

  for (unsigned int i = 0; i < n_replicates; i++) {
    Fire fire = simulate_fire(
        landscape, ignition_cells, params, distance, elevation_mean, elevation_sd, upper_limit
    );

    for (unsigned int col = 0; col < landscape.width; col++) {
      for (unsigned int row = 0; row < landscape.height; row++) {
        if (fire.burned_layer[{col, row}]) {
          burned_amounts[{col, row}] += 1;
        }
      }
    }
  }

  return burned_amounts;
}
