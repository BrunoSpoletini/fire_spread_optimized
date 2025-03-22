#include "spread_functions.hpp"

#define _USE_MATH_DEFINES
#include <cmath>
#include <random>
#include <vector>
#include <omp.h>
#include <cstring>

#include "fires.hpp"
#include "landscape.hpp"

float spread_probability(
    const Cell& burning, const Cell& neighbour, SimulationParams params, float angle,
    float distance, float elevation_mean, float elevation_sd, float upper_limit = 1.0
) {

  float slope_term = sin(atan((neighbour.elevation - burning.elevation) / distance));
  float wind_term = cos(angle - burning.wind_direction);
  float elev_term = (neighbour.elevation - elevation_mean) / elevation_sd;

  float linpred = params.independent_pred;

  if (neighbour.vegetation_type == SUBALPINE) {
    linpred += params.subalpine_pred;
  } else if (neighbour.vegetation_type == WET) {
    linpred += params.wet_pred;
  } else if (neighbour.vegetation_type == DRY) {
    linpred += params.dry_pred;
  }

  linpred += params.fwi_pred * neighbour.fwi;
  linpred += params.aspect_pred * neighbour.aspect;

  linpred += wind_term * params.wind_pred + elev_term * params.elevation_pred +
             slope_term * params.slope_pred;

             float prob = upper_limit / (1 + exp(-linpred));

  return prob;
}

Fire simulate_fire(
    const Landscape& landscape, const std::vector<std::pair<unsigned int, unsigned int>>& ignition_cells,
    SimulationParams params, float distance, float elevation_mean, float elevation_sd,
    float upper_limit = 1.0
) {

  unsigned int n_row = landscape.height;
  unsigned int n_col = landscape.width;

  std::vector<std::pair<unsigned int, unsigned int>> burned_ids;

  unsigned int start = 0;
  unsigned int end = ignition_cells.size();

  for (unsigned int i = 0; i < end; i++) {
    burned_ids.push_back(ignition_cells[i]);
  }

  std::vector<unsigned int> burned_ids_steps;
  burned_ids_steps.push_back(end);

  unsigned int burning_size = end + 1;

  Matrix<bool> burned_bin = Matrix<bool>(n_col, n_row);

  for (unsigned int i = 0; i < end; i++) {
    unsigned int cell_0 = ignition_cells[i].first;
    unsigned int cell_1 = ignition_cells[i].second;
    burned_bin[{ cell_0, cell_1 }] = 1;
  }

  double t = omp_get_wtime();
  while (burning_size > 0) {
    unsigned int end_forward = end;

    // Loop over burning cells in the cycle

    // b is going to keep the position in burned_ids that have to be evaluated
    // in this burn cycle
    for (unsigned int b = start; b < end; b++) {
      unsigned int burning_cell_0 = burned_ids[b].first;
      unsigned int burning_cell_1 = burned_ids[b].second;

      const Cell& burning_cell = landscape[{ burning_cell_0, burning_cell_1 }];

      constexpr int moves[8][2] = { { -1, -1 }, { -1, 0 }, { -1, 1 }, { 0, -1 },
                                    { 0, 1 },   { 1, -1 }, { 1, 0 },  { 1, 1 } };

      int neighbors_coords[2][8];

      for (unsigned int i = 0; i < 8; i++) {
        neighbors_coords[0][i] = int(burning_cell_0) + moves[i][0];
        neighbors_coords[1][i] = int(burning_cell_1) + moves[i][1];
      }
      // Note that in the case 0 - 1 we will have size_t_MAX

      // Loop over neighbors_coords of the focal burning cell

      for (unsigned int n = 0; n < 8; n++) {

        int neighbour_cell_0 = neighbors_coords[0][n];
        int neighbour_cell_1 = neighbors_coords[1][n];

        // Is the cell in range?
        bool out_of_range = 0 > neighbour_cell_0 || neighbour_cell_0 >= int(n_col) ||
                            0 > neighbour_cell_1 || neighbour_cell_1 >= int(n_row);

        if (out_of_range)
          continue;

        const Cell& neighbour_cell = landscape[{ neighbour_cell_0, neighbour_cell_1 }];

        // Is the cell burnable?
        bool burnable_cell =
            !burned_bin[{ neighbour_cell_0, neighbour_cell_1 }] && neighbour_cell.burnable;

        if (!burnable_cell)
          continue;

        constexpr float angles[8] = { M_PI * 3 / 4, M_PI, M_PI * 5 / 4, M_PI / 2, M_PI * 3 / 2,
                                       M_PI / 4,     0,    M_PI * 7 / 4 };

        // simulate fire
        float prob = spread_probability(
            burning_cell, neighbour_cell, params, angles[n], distance, elevation_mean,
            elevation_sd, upper_limit
        );

        // Burn with probability prob (Bernoulli)
        bool burn = (float)rand() / (float)RAND_MAX < prob;

        if (burn == 0)
          continue;

        // If burned, store id of recently burned cell and set 1 in burned_bin
        end_forward += 1;
        burned_ids.push_back({ neighbour_cell_0, neighbour_cell_1 });
        burned_bin[{ neighbour_cell_0, neighbour_cell_1 }] = true;
      }
    }

    // update start and end
    start = end;
    end = end_forward;
    burning_size = end - start;

    burned_ids_steps.push_back(end);
  }
  //printf("celdas incendiadas por segundo : %lf\n", burned_ids.size() / (omp_get_wtime() - t));
  fprintf(stderr, "celdas incendiadas por microsegundo: %lf\n", burned_ids.size() / (1E06 * (omp_get_wtime()-t)));

  return { n_col, n_row, burned_bin, burned_ids, burned_ids_steps };
}
