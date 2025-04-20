#include "spread_functions.hpp"

#define _USE_MATH_DEFINES
#include <cmath>
#include <random>
#include <vector>
#include <omp.h>
#include <cstring>

#include "fires.hpp"
#include "landscape.hpp"

#pragma omp declare simd
float spread_probability(
    float burning_cell_elevation, float burning_cell_wind_direction, float neighbour_fwi, float neighbour_aspect, float neighbour_elevation, 
    float linpred,
    float fwi_pred, float aspect_pred, float wind_pred,
    float elevation_pred, float slope_pred, float angle, float distance, float elevation_mean,
    float elevation_sd, float upper_limit = 1.0
) {
  float slope_term = sin(atan((neighbour_elevation - burning_cell_elevation) / distance));
  float wind_term = cos(angle - burning_cell_wind_direction);
  float elev_term = (neighbour_elevation - elevation_mean) / elevation_sd;

  linpred += fwi_pred * neighbour_fwi;
  linpred += aspect_pred * neighbour_aspect;

  linpred += wind_term * wind_pred + elev_term * elevation_pred + slope_term * slope_pred;

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

  constexpr float angles[8] = { M_PI * 3 / 4, M_PI, M_PI * 5 / 4, M_PI / 2, M_PI * 3 / 2,
                                M_PI / 4,     0,    M_PI * 7 / 4 };

  constexpr int moves[8][2] = { { -1, -1 }, { -1, 0 }, { -1, 1 }, { 0, -1 },
                                { 0, 1 },   { 1, -1 }, { 1, 0 },  { 1, 1 } };

  double t = omp_get_wtime();
  while (burning_size > 0) {
    unsigned int end_forward = end;

    // Loop over burning cells in the cycle

    // b is going to keep the position in burned_ids that have to be evaluated
    // in this burn cycle

		// LOOP POR CELULAS QUEMADAS
    for (unsigned int b = start; b < end; b++) {
      unsigned int burning_cell_0 = burned_ids[b].first;
      unsigned int burning_cell_1 = burned_ids[b].second;

      const Cell& burning_cell = landscape[{ burning_cell_0, burning_cell_1 }];

      // Loop over neighbors_coords of the focal burning cell

      // Arreglos para almacenar los resultados del bucle vectorizado
      int neighbor_x[8], neighbor_y[8];
      bool in_range[8], is_burnable[8];
      bool should_burn[8];
			float prob[8];
			Cell neighbour_cell_list[8];

      // SOA for neighbour_cell_list
      float neighbour_fwi[8];
      float neighbour_aspect[8];
      float neighbour_elevation[8];
      float neighbour_vegetation_type[8];
      float linpred_list[8];
      float burning_cell_elevation = burning_cell.elevation;
      float burning_cell_wind_direction = burning_cell.wind_direction;

      #pragma omp simd
      for (unsigned int i = 0; i < 8; i++) {
        neighbor_x[i] = int(burning_cell_0) + moves[i][0];
        neighbor_y[i] = int(burning_cell_1) + moves[i][1];
			}
			
      #pragma omp simd
      for (unsigned int i = 0; i < 8; i++) {
        // Inicializamos la probabilidad de quema en 0
        prob[i] = 0.0;
        // Verificar si está en rango
        in_range[i] = !(0 > neighbor_x[i] || neighbor_x[i] >= int(n_col) ||
                       0 > neighbor_y[i] || neighbor_y[i] >= int(n_row));

				// Verificar si está quemado
				if (in_range[i]) {
          // Load SOA
          neighbour_cell_list[i] = landscape[{ neighbor_x[i], neighbor_y[i] }];
          neighbour_fwi[i] = neighbour_cell_list[i].fwi;
          neighbour_aspect[i] = neighbour_cell_list[i].aspect;
          neighbour_elevation[i] = neighbour_cell_list[i].elevation;
          neighbour_vegetation_type[i] = neighbour_cell_list[i].vegetation_type;

					// Verificar si es quemable
					is_burnable[i] = !burned_bin[{neighbor_x[i], neighbor_y[i]}] && neighbour_cell_list[i].burnable;
          if (is_burnable[i]) {
            // linpred
            linpred_list[i] = params.independent_pred;
            if (neighbour_vegetation_type[i] == SUBALPINE) {
              linpred_list[i] += params.subalpine_pred;
            } else if (neighbour_vegetation_type[i] == WET) {
                linpred_list[i] += params.wet_pred;
            } else if (neighbour_vegetation_type[i] == DRY) {
                linpred_list[i] += params.dry_pred;
            }
          }

				}
			}
			// Calcular la probabilidad de que se queme cada cel neighbor_cell
			#pragma omp simd
			for (int n = 0; n < 8; n++) {

        float mask = (in_range[n] && is_burnable[n]) ? 1.0f : 0.0f;

				//if (in_range[n] && is_burnable[n]){
        	// Calcular probabilidad

					float temp_prob = {spread_probability(
						burning_cell_elevation, burning_cell_wind_direction, neighbour_fwi[n], neighbour_aspect[n], neighbour_elevation[n], 
            linpred_list[n],
						params.fwi_pred, params.aspect_pred,
						params.wind_pred, params.elevation_pred, params.slope_pred, angles[n],
						distance, elevation_mean, elevation_sd, upper_limit
				)};

        prob[n] = temp_prob * mask;
			//}
		}
			for (unsigned int n = 0; n < 8; n++) {
					if (in_range[n] && is_burnable[n]){
					// Determinar si se quema
					should_burn[n] = (float)rand() / (float)RAND_MAX < prob[n];
					
					if (should_burn[n]) {
						// Si se quema, actualizar las estructuras de datos
						end_forward += 1;
						burned_ids.push_back({neighbor_x[n], neighbor_y[n]});
						burned_bin[{neighbor_x[n], neighbor_y[n]}] = true;
					}
				}
      }
    }

    // update start and end
    start = end;
    end = end_forward;
    burning_size = end - start;

    burned_ids_steps.push_back(end);
  }
  fprintf(stderr, "celdas incendiadas por microsegundo: %lf\n", burned_ids.size() / (1E06 * (omp_get_wtime()-t)));

  return { n_col, n_row, burned_bin, burned_ids, burned_ids_steps };
}
