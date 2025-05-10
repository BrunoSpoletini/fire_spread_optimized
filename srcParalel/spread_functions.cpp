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
  int max_threads = omp_get_max_threads();
  unsigned int base_seed = rd();
  std::vector<std::mt19937> thread_generators(max_threads);
  for (int i = 0; i < max_threads; ++i) {
      thread_generators[i].seed(base_seed + i);
  }

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
    
    // Create a thread-local buffer for new burned cells
    std::vector<std::vector<std::pair<unsigned int, unsigned int>>> thread_local_burned_ids(max_threads);
    
    // Parallel region for processing burning cells
    #pragma omp parallel 
    {
      int thread_id = omp_get_thread_num();
      std::mt19937& local_gen = thread_generators[thread_id];
      std::uniform_real_distribution<float> local_dist(0.0f, 1.0f);
      
      // Thread-local arrays
      int local_neighbor_x[8], local_neighbor_y[8];
      bool local_in_range[8], local_is_burnable[8], local_should_burn[8];
      float local_prob[8],
            local_neighbour_fwi[8],
            local_neighbour_aspect[8],
            local_neighbour_elevation[8],
            local_neighbour_vegetation_type[8],
            local_linpred_list[8],
            local_burning_cell_elevation[8],
            local_burning_cell_wind_direction[8];
      
      // Parallelize the loop over burning cells
      #pragma omp for 
      for (unsigned int b = start; b < end; b++) {
        unsigned int burning_cell_0 = burned_ids[b].first;
        unsigned int burning_cell_1 = burned_ids[b].second;

        const Cell& burning_cell = landscape[{ burning_cell_0, burning_cell_1 }];

        for (unsigned int i = 0; i < 8; i++) {
          local_neighbor_x[i] = int(burning_cell_0) + moves[i][0];
          local_neighbor_y[i] = int(burning_cell_1) + moves[i][1];

          // Inicializamos la probabilidad de quema en 0
          local_prob[i] = 0.0;
          // Verificar si está en rango
          local_in_range[i] = !(0 > local_neighbor_x[i] || local_neighbor_x[i] >= int(n_col) ||
                         0 > local_neighbor_y[i] || local_neighbor_y[i] >= int(n_row));

          // Verificar si está quemado
          if (local_in_range[i]) {
            // Use critical section to check burned_bin
            bool already_burned;
            //#pragma omp critical(burned_bin_read)
            //{
                already_burned = burned_bin[{local_neighbor_x[i], local_neighbor_y[i]}];
            //}
            
            // Verificar si es quemable
            local_is_burnable[i] = !already_burned && landscape[{ local_neighbor_x[i], local_neighbor_y[i] }].burnable;
            if (local_is_burnable[i]) {
              // Load SOA
              local_neighbour_fwi[i] = landscape[{ local_neighbor_x[i], local_neighbor_y[i] }].fwi;
              local_neighbour_aspect[i] = landscape[{ local_neighbor_x[i], local_neighbor_y[i] }].aspect;
              local_neighbour_elevation[i] = landscape[{ local_neighbor_x[i], local_neighbor_y[i] }].elevation;
              local_neighbour_vegetation_type[i] = landscape[{ local_neighbor_x[i], local_neighbor_y[i] }].vegetation_type;
              local_burning_cell_elevation[i] = burning_cell.elevation;
              local_burning_cell_wind_direction[i] = burning_cell.wind_direction;

              // linpred
              local_linpred_list[i] = params.independent_pred;
              if (local_neighbour_vegetation_type[i] == SUBALPINE) {
                local_linpred_list[i] += params.subalpine_pred;
              } else if (local_neighbour_vegetation_type[i] == WET) {
                  local_linpred_list[i] += params.wet_pred;
              } else if (local_neighbour_vegetation_type[i] == DRY) {
                  local_linpred_list[i] += params.dry_pred;
              }
            }
          }
          
          if (!local_in_range[i] || !local_is_burnable[i]) {
            local_neighbour_fwi[i] = 0.0f;
            local_neighbour_aspect[i] = 0.0f;
            local_neighbour_elevation[i] = 0.0f;
            local_neighbour_vegetation_type[i] = 0.0f;
            local_linpred_list[i] = 0.0f;
            local_burning_cell_elevation[i] = 0.0f;
            local_burning_cell_wind_direction[i] = 0.0f;
          }
        }

        spread_probability_vectorized(
          local_burning_cell_elevation, 
          local_burning_cell_wind_direction, 
          local_neighbour_fwi,
          local_neighbour_aspect, 
          local_neighbour_elevation, 
          local_linpred_list,
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
          local_prob,
          batch_size
        );

        for (unsigned int n = 0; n < 8; n++) {
          float mask = (local_in_range[n] && local_is_burnable[n]) ? 1.0f : 0.0f;
          local_prob[n] = local_prob[n] * mask;
          local_prob[n] = local_prob[n] * 0.96f; // To counteract the effects of exp aproximation

          // Determinar si se quema
          if (local_prob[n] > 0.0f) {
            local_should_burn[n] = local_dist(local_gen) < local_prob[n];
            
            if (local_should_burn[n]) {
              // Use critical sections to update the burned_bin matrix
              bool already_marked = false;
              //#pragma omp critical(burned_bin_read)
              //{
                  //already_marked = burned_bin[{local_neighbor_x[n], local_neighbor_y[n]}];
              //}
              //if (!already_marked) {
                  //#pragma omp critical(burned_bin_write)
                  //{
                      //burned_bin[{local_neighbor_x[n], local_neighbor_y[n]}] = true;
                  //}
                  // Add to thread-local buffer
                  //thread_local_burned_ids[thread_id].push_back({local_neighbor_x[n], local_neighbor_y[n]});
              if (!burned_bin[{local_neighbor_x[n], local_neighbor_y[n]}]) {
                thread_local_burned_ids[thread_id].emplace_back(local_neighbor_x[n], local_neighbor_y[n]);
              }


                // #pragma omp critical
                // {
                //   if (!burned_bin[{local_neighbor_x[n], local_neighbor_y[n]}]) {
                //     burned_bin[{local_neighbor_x[n], local_neighbor_y[n]}] = true;
                //     burned_ids.push_back({local_neighbor_x[n], local_neighbor_y[n]});
                //     end_forward++;
                //   }
                // }
              //}
            }
          }
        }
      }
    }
    
    // Merge thread-local burned cells into the main vector
    for (int t = 0; t < max_threads; t++) {
      for (const auto& cell : thread_local_burned_ids[t]) {
        auto [x, y] = cell;
        if (!burned_bin[{x, y}]) {
          burned_bin[{x, y}] = true;
          burned_ids.push_back(cell);
          end_forward++;
        }
      }
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

