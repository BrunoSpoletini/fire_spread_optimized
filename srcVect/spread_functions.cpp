#include "spread_functions.hpp"

#define _USE_MATH_DEFINES
#include <cmath>

#include <random>
#include <vector>
#include <omp.h>
#include <cstring>

#include "fires.hpp"
#include "landscape.hpp"

#include "spread_functions_vectorized.hpp"

Fire simulate_fire(
    const Landscape& landscape, const std::vector<std::pair<unsigned int, unsigned int>>& ignition_cells,
    SimulationParams params, float distance, float elevation_mean, float elevation_sd,
    float upper_limit = 1.0
  ) {
  double t = omp_get_wtime();

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


  // random number generation
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);


  // Arreglos para almacenar los resultados del bucle vectorizado
  int neighbor_x[8], neighbor_y[8];
  bool in_range[8], is_burnable[8], should_burn[8];
  Cell neighbour_cell_list[8];
  // SOA for neighbour_cell_list
  float prob[8],
        neighbour_fwi[8],
        neighbour_aspect[8],
        neighbour_elevation[8],
        neighbour_vegetation_type[8],
        linpred_list[8],
        burning_cell_elevation[8],
        burning_cell_wind_direction[8];
  int batch_size = 8;

  while (burning_size > 0) {
    unsigned int end_forward = end;

    // Loop over burning cells in the cycle

		// LOOP POR CELULAS QUEMADAS
    for (unsigned int b = start; b < end; b++) {
      unsigned int burning_cell_0 = burned_ids[b].first;
      unsigned int burning_cell_1 = burned_ids[b].second;

      const Cell& burning_cell = landscape[{ burning_cell_0, burning_cell_1 }];

      for (unsigned int i = 0; i < 8; i++) {
        neighbor_x[i] = int(burning_cell_0) + moves[i][0];
        neighbor_y[i] = int(burning_cell_1) + moves[i][1];

        // Inicializamos la probabilidad de quema en 0
        prob[i] = 0.0;
        // Verificar si está en rango
        in_range[i] = !(0 > neighbor_x[i] || neighbor_x[i] >= int(n_col) ||
                       0 > neighbor_y[i] || neighbor_y[i] >= int(n_row));

				// Verificar si está quemado
				if (in_range[i]) {
          
					// Verificar si es quemable
					is_burnable[i] = !burned_bin[{neighbor_x[i], neighbor_y[i]}] && landscape[{ neighbor_x[i], neighbor_y[i] }].burnable;
          if (is_burnable[i]) {
            // Load SOA
            neighbour_fwi[i] = landscape[{ neighbor_x[i], neighbor_y[i] }].fwi;
            neighbour_aspect[i] = landscape[{ neighbor_x[i], neighbor_y[i] }].aspect;
            neighbour_elevation[i] = landscape[{ neighbor_x[i], neighbor_y[i] }].elevation;
            neighbour_vegetation_type[i] = landscape[{ neighbor_x[i], neighbor_y[i] }].vegetation_type;
            burning_cell_elevation[i] = burning_cell.elevation;
            burning_cell_wind_direction[i] = burning_cell.wind_direction;

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
        if (!in_range[i] || !is_burnable[i]) {
          neighbour_fwi[i] = 0.0f;
          neighbour_aspect[i] = 0.0f;
          neighbour_elevation[i] = 0.0f;
          neighbour_vegetation_type[i] = 0.0f;
          linpred_list[i] = 0.0f;
          burning_cell_elevation[i] = 0.0f;
          burning_cell_wind_direction[i] = 0.0f;
        }
			}

      spread_probability_vectorized(
        burning_cell_elevation, 
        burning_cell_wind_direction, 
        neighbour_fwi,
        neighbour_aspect, 
        neighbour_elevation, 
        linpred_list,
        params.fwi_pred, 
        params.aspect_pred,
        params.wind_pred, 
        params.elevation_pred,
        params.slope_pred,
        angles,
        distance, 
        elevation_mean, 
        elevation_sd, 
        upper_limit,
        prob,
        batch_size
      );

			for (unsigned int n = 0; n < 8; n++) {
        float mask = (in_range[n] && is_burnable[n]) ? 1.0f : 0.0f;
        prob[n] = prob[n] * mask;
        prob[n] = prob[n] * 0.96f; // To counteract the effects of exp aproximation (tends to underestimate -> 1/(1+exp(-linpred)) yields hig)

					// Determinar si se quema
        if (prob[n] > 0.0f) {
          should_burn[n] = dist(gen) < prob[n];
					
          if (should_burn[n]) {
            // Si se quema, actualizar las estructuras de datos
            end_forward += 1;
            burned_ids.push_back({neighbor_x[n], neighbor_y[n]});
            burned_bin[{neighbor_x[n], neighbor_y[n]}] = true;
          }
        }
      }
        //}
    }

    // update start and end
    start = end;
    end = end_forward;
    burning_size = end - start;

    burned_ids_steps.push_back(end);
  }
  fprintf(stderr, "Celdas incendiadas: %ld\n", burned_ids.size());
  fprintf(stderr, "celdas incendiadas por microsegundo: %lf\n", burned_ids.size() / (1E06 * (omp_get_wtime()-t)));

  return { n_col, n_row, burned_bin, burned_ids, burned_ids_steps };
}

