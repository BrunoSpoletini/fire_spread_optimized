#include "spread_functions.hpp"

#define _USE_MATH_DEFINES
#include <cmath>
#include <random>
#include <vector>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include "fires.hpp"
#include "landscape.hpp"

// Helper function for CUDA error checking
#define CUDA_CHECK(call) \
    do { \
        cudaError_t error = call; \
        if (error != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(error)); \
            exit(1); \
        } \
    } while(0)

__device__ float spread_probability(
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

// CUDA kernel for spread probability calculation
__global__ void calculate_spread_probabilities(
    const Cell* __restrict__ landscape,
    const bool* __restrict__ burned_bin,
    const unsigned int* __restrict__ burning_cells,
    unsigned int n_burning_cells,
    unsigned int n_col,
    unsigned int n_row,
    const SimulationParams* params,
    float distance,
    float elevation_mean,
    float elevation_sd,
    float upper_limit,
    curandState* states,
    unsigned int* new_burned_cells,
    unsigned int* n_new_burned
) {
    const float angles[8] = { M_PI * 3 / 4, M_PI, M_PI * 5 / 4, M_PI / 2, M_PI * 3 / 2,
                              M_PI / 4,     0,    M_PI * 7 / 4 };
    const int moves[8][2] = { { -1, -1 }, { -1, 0 }, { -1, 1 }, { 0, -1 },
                              { 0, 1 },   { 1, -1 }, { 1, 0 },  { 1, 1 } };

    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_burning_cells) return;

    // Get burning cell coordinates
    unsigned int burning_cell_0 = burning_cells[idx * 2];
    unsigned int burning_cell_1 = burning_cells[idx * 2 + 1];
    const Cell& burning_cell = landscape[burning_cell_1 * n_col + burning_cell_0];

    // Get the random state for this cell
    unsigned int cell_idx = burning_cell_1 * n_col + burning_cell_0;
    curandState localState = states[cell_idx];

    // Process each neighbor
    for (int n = 0; n < 8; n++) {
        int neighbor_x = burning_cell_0 + moves[n][0];
        int neighbor_y = burning_cell_1 + moves[n][1];

        // Check if neighbor is in range
        if (neighbor_x < 0 || neighbor_x >= n_col || neighbor_y < 0 || neighbor_y >= n_row) {
            continue;
        }

        // Check if already burned
        if (burned_bin[neighbor_y * n_col + neighbor_x]) {
            continue;
        }

        const Cell& neighbor = landscape[neighbor_y * n_col + neighbor_x];
        if (!neighbor.burnable) {
            continue;
        }

        // Calculate spread probability using the helper function
        float prob = spread_probability(burning_cell, neighbor, *params, angles[n], 
                                      distance, elevation_mean, elevation_sd, upper_limit);

        // Random number generation and burn decision
        float random_value = curand_uniform(&localState);

        if (random_value < prob) {
            // Atomically add new burned cell
            unsigned int new_idx = atomicAdd(n_new_burned, 1);
            new_burned_cells[new_idx * 2] = neighbor_x;
            new_burned_cells[new_idx * 2 + 1] = neighbor_y;
        }
    }

    // Save the updated random state
    states[cell_idx] = localState;
}

// Initialize random states
__global__ void init_random_states(curandState* states, unsigned int seed) {
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    curand_init(seed + idx, 0, 0, &states[idx]);
}

Fire simulate_fire(
    const Landscape& landscape, const std::vector<std::pair<unsigned int, unsigned int>>& ignition_cells,
    SimulationParams params, float distance, float elevation_mean, float elevation_sd,
    float upper_limit = 1.0
) {
    // Use CUDA timing instead of OpenMP
    cudaEvent_t start_event, stop_event;
    CUDA_CHECK(cudaEventCreate(&start_event));
    CUDA_CHECK(cudaEventCreate(&stop_event));
    CUDA_CHECK(cudaEventRecord(start_event));

    unsigned int n_row = landscape.height;
    unsigned int n_col = landscape.width;

    // Initialize host arrays
    std::vector<std::pair<unsigned int, unsigned int>> burned_ids;
    std::vector<unsigned int> burned_ids_steps;
    Matrix<bool> burned_bin(n_col, n_row);

    // Add ignition cells
    for (const auto& cell : ignition_cells) {
        burned_ids.push_back(cell);
        burned_bin[{cell.first, cell.second}] = true;
    }
    burned_ids_steps.push_back(ignition_cells.size());

    // Allocate device memory with error checking
    Cell* d_landscape = nullptr;
    bool* d_burned_bin = nullptr;
    unsigned int* d_burning_cells = nullptr;
    unsigned int* d_new_burned_cells = nullptr;
    unsigned int* d_n_new_burned = nullptr;
    curandState* d_states = nullptr;
    SimulationParams* d_params = nullptr;

    size_t landscape_size = n_col * n_row * sizeof(Cell);
    size_t burned_bin_size = n_col * n_row * sizeof(bool);
    size_t burning_cells_size = 2 * n_col * n_row * sizeof(unsigned int);
    size_t new_burned_cells_size = 2 * n_col * n_row * sizeof(unsigned int);

    CUDA_CHECK(cudaMalloc(&d_landscape, landscape_size));
    CUDA_CHECK(cudaMalloc(&d_burned_bin, burned_bin_size));
    CUDA_CHECK(cudaMalloc(&d_burning_cells, burning_cells_size));
    CUDA_CHECK(cudaMalloc(&d_new_burned_cells, new_burned_cells_size));
    CUDA_CHECK(cudaMalloc(&d_n_new_burned, sizeof(unsigned int)));
    CUDA_CHECK(cudaMalloc(&d_states, n_col * n_row * sizeof(curandState)));
    CUDA_CHECK(cudaMalloc(&d_params, sizeof(SimulationParams)));

    // Initialize device memory to zero
    CUDA_CHECK(cudaMemset(d_burned_bin, 0, burned_bin_size));
    CUDA_CHECK(cudaMemset(d_burning_cells, 0, burning_cells_size));
    CUDA_CHECK(cudaMemset(d_new_burned_cells, 0, new_burned_cells_size));
    CUDA_CHECK(cudaMemset(d_n_new_burned, 0, sizeof(unsigned int)));

    // Copy params to device
    CUDA_CHECK(cudaMemcpy(d_params, &params, sizeof(SimulationParams), cudaMemcpyHostToDevice));

    // Copy landscape to device
    Cell* landscape_data = new Cell[n_col * n_row];
    if (!landscape_data) {
        fprintf(stderr, "Failed to allocate landscape_data\n");
        exit(1);
    }

    for (unsigned int i = 0; i < n_col; i++) {
        for (unsigned int j = 0; j < n_row; j++) {
            landscape_data[j * n_col + i] = landscape[{i, j}];
        }
    }
    CUDA_CHECK(cudaMemcpy(d_landscape, landscape_data, landscape_size, cudaMemcpyHostToDevice));
    delete[] landscape_data;

    // Copy burned_bin to device
    bool* burned_bin_data = new bool[n_col * n_row];
    if (!burned_bin_data) {
        fprintf(stderr, "Failed to allocate burned_bin_data\n");
        exit(1);
    }

    for (unsigned int i = 0; i < n_col; i++) {
        for (unsigned int j = 0; j < n_row; j++) {
            burned_bin_data[j * n_col + i] = burned_bin[{i, j}];
        }
    }
    CUDA_CHECK(cudaMemcpy(d_burned_bin, burned_bin_data, burned_bin_size, cudaMemcpyHostToDevice));
    delete[] burned_bin_data;

    // Initialize random states
    int block_size = 256;
    int num_blocks = (n_col * n_row + block_size - 1) / block_size;
    init_random_states<<<num_blocks, block_size>>>(d_states, time(NULL));
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    unsigned int current_start = 0;
    unsigned int current_end = ignition_cells.size();
    unsigned int burning_size = current_end;

    while (burning_size > 0) {
        // Verify indices are valid
        if (current_start > burned_ids.size()) {
            fprintf(stderr, "Error: current_start (%u) exceeds burned_ids size (%zu)\n", 
                    current_start, burned_ids.size());
            break;
        }

        // Calculate actual burning size based on available cells
        burning_size = std::min(burning_size, static_cast<unsigned int>(burned_ids.size() - current_start));
        if (burning_size == 0) break;

        // Copy current burning cells to device
        unsigned int* current_burning = new unsigned int[2 * burning_size];
        if (!current_burning) {
            fprintf(stderr, "Failed to allocate current_burning\n");
            exit(1);
        }

        // Copy burning cells with bounds checking
        for (unsigned int i = 0; i < burning_size; i++) {
            unsigned int idx = current_start + i;
            if (idx >= burned_ids.size()) {
                fprintf(stderr, "Error: Index %u out of bounds in burned_ids (size: %zu)\n", 
                        idx, burned_ids.size());
                delete[] current_burning;
                exit(1);
            }
            current_burning[i * 2] = burned_ids[idx].first;
            current_burning[i * 2 + 1] = burned_ids[idx].second;
        }

        CUDA_CHECK(cudaMemcpy(d_burning_cells, current_burning, 2 * burning_size * sizeof(unsigned int), cudaMemcpyHostToDevice));
        delete[] current_burning;

        // Reset counter for new burned cells
        unsigned int zero = 0;
        CUDA_CHECK(cudaMemcpy(d_n_new_burned, &zero, sizeof(unsigned int), cudaMemcpyHostToDevice));

        // Launch kernel
        num_blocks = (burning_size + block_size - 1) / block_size;
        calculate_spread_probabilities<<<num_blocks, block_size>>>(
            d_landscape, d_burned_bin, d_burning_cells, burning_size,
            n_col, n_row, d_params, distance, elevation_mean, elevation_sd,
            upper_limit, d_states, d_new_burned_cells, d_n_new_burned
        );
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());

        // Get number of new burned cells
        unsigned int n_new_burned = 0;
        CUDA_CHECK(cudaMemcpy(&n_new_burned, d_n_new_burned, sizeof(unsigned int), cudaMemcpyDeviceToHost));

        // Debug print
        //fprintf(stderr, "n_new_burned = %u\n", n_new_burned);
        if (n_new_burned > 2 * n_col * n_row) {
            fprintf(stderr, "ERROR: n_new_burned (%u) exceeds max possible cells (%u)!\n", n_new_burned, 2 * n_col * n_row);
            exit(1);
        }

        if (n_new_burned > 0) {
            // Copy new burned cells back to host
            unsigned int* new_burned = new unsigned int[2 * n_new_burned];
            if (!new_burned) {
                fprintf(stderr, "Failed to allocate new_burned\n");
                exit(1);
            }

            CUDA_CHECK(cudaMemcpy(new_burned, d_new_burned_cells, 2 * n_new_burned * sizeof(unsigned int), cudaMemcpyDeviceToHost));

            // Update burned_bin and burned_ids
            for (unsigned int i = 0; i < n_new_burned; i++) {
                unsigned int x = new_burned[i * 2];
                unsigned int y = new_burned[i * 2 + 1];
                
                // Strict bounds checking
                if (x >= n_col || y >= n_row) {
                    fprintf(stderr, "Warning: Out of bounds cell coordinates (%u, %u) with grid size (%u, %u)\n", 
                            x, y, n_col, n_row);
                    continue;
                }

                // Check if cell is already burned
                if (burned_bin[{x, y}]) {
                    continue;
                }

                burned_bin[{x, y}] = true;
                burned_ids.push_back({x, y});
            }
            delete[] new_burned;
        }

        // Update indices
        current_start = current_end;
        current_end = burned_ids.size();  // Use actual size instead of adding n_new_burned
        burning_size = n_new_burned;
        burned_ids_steps.push_back(current_end);
    }

    // Free device memory
    if (d_landscape) CUDA_CHECK(cudaFree(d_landscape));
    if (d_burned_bin) CUDA_CHECK(cudaFree(d_burned_bin));
    if (d_burning_cells) CUDA_CHECK(cudaFree(d_burning_cells));
    if (d_new_burned_cells) CUDA_CHECK(cudaFree(d_new_burned_cells));
    if (d_n_new_burned) CUDA_CHECK(cudaFree(d_n_new_burned));
    if (d_states) CUDA_CHECK(cudaFree(d_states));
    if (d_params) CUDA_CHECK(cudaFree(d_params));

    // Get elapsed time using CUDA events
    CUDA_CHECK(cudaEventRecord(stop_event));
    CUDA_CHECK(cudaEventSynchronize(stop_event));
    float milliseconds = 0;
    CUDA_CHECK(cudaEventElapsedTime(&milliseconds, start_event, stop_event));
    double seconds = milliseconds / 1000.0;

    fprintf(stderr, "Celdas incendiadas: %ld\n", burned_ids.size());
    fprintf(stderr, "celdas incendiadas por microsegundo: %lf\n", burned_ids.size() / (1E06 * seconds));

    // Clean up CUDA events
    CUDA_CHECK(cudaEventDestroy(start_event));
    CUDA_CHECK(cudaEventDestroy(stop_event));

    return { n_col, n_row, burned_bin, burned_ids, burned_ids_steps };
}
