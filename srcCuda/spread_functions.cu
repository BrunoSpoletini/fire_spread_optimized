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
// Each thread processes one cell, if some of its neighbors are burning,
// it calculates the probability of burning itself based on the parameters
__global__ void calculate_spread_probabilities(
    const Cell* __restrict__ landscape,
    const unsigned short* __restrict__ burned_ids,
    unsigned int burned_size,
    unsigned int n_col,
    unsigned int n_row,
    const SimulationParams* params,
    float distance,
    float elevation_mean,
    float elevation_sd,
    float upper_limit,
    curandState* states,
    bool* burned_bin,
    char* cell_states_initial,
    char* cell_states_final
) {
    const float angles[8] = { M_PI * 3 / 4, M_PI, M_PI * 5 / 4, M_PI / 2, M_PI * 3 / 2,
                              M_PI / 4,     0,    M_PI * 7 / 4 };
    const int moves[8][2] = { { -1, -1 }, { -1, 0 }, { -1, 1 }, { 0, -1 },
                              { 0, 1 },   { 1, -1 }, { 1, 0 },  { 1, 1 } };

    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Loop
    if (idx >= (n_col * n_row)) return; // Ensure we don't exceed the number of cells in the grid
    
    // Get the current cell coordinates
    unsigned int cell_x = idx % n_col;
    unsigned int cell_y = idx / n_col;
    const Cell& current_cell = landscape[cell_y * n_col + cell_x];


    // // Get burning cell coordinates
    // unsigned int burning_cell_0 = burning_cells[idx * 2];
    // unsigned int burning_cell_1 = burning_cells[idx * 2 + 1];
    // const Cell& burning_cell = landscape[burning_cell_1 * n_col + burning_cell_0];

    // // Get the random state for this cell
    // unsigned int cell_idx = burning_cell_1 * n_col + burning_cell_0;
    // curandState localState = states[cell_idx];

    // // Process each neighbor
    // for (int n = 0; n < 8; n++) {
    //     int neighbor_x = burning_cell_0 + moves[n][0];
    //     int neighbor_y = burning_cell_1 + moves[n][1];

    //     // Check if neighbor is in range
    //     if (neighbor_x < 0 || neighbor_x >= n_col || neighbor_y < 0 || neighbor_y >= n_row) {
    //         continue;
    //     }

    //     // Check if already burned
    //     if (burned_bin[neighbor_y * n_col + neighbor_x]) {
    //         continue;
    //     }

    //     const Cell& neighbor = landscape[neighbor_y * n_col + neighbor_x];
    //     if (!neighbor.burnable) {
    //         continue;
    //     }

    //     // Calculate spread probability using the helper function
    //     float prob = spread_probability(burning_cell, neighbor, *params, angles[n], 
    //                                   distance, elevation_mean, elevation_sd, upper_limit);

    //     // Random number generation and burn decision
    //     float random_value = curand_uniform(&localState);

    //     if (random_value < prob) {
    //         // Atomically add new burned cell
    //         unsigned int new_idx = atomicAdd(n_new_burned, 1);
    //         new_burned_cells[new_idx * 2] = neighbor_x;
    //         new_burned_cells[new_idx * 2 + 1] = neighbor_y;
    //     }
    // }

    // // Save the updated random state
    // states[cell_idx] = localState;
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
    // std::vector<std::pair<unsigned int, unsigned int>> burned_ids;
    // std::vector<unsigned int> burned_ids_steps;
    //Matrix<bool> burned_bin(n_col, n_row);

    // Host variables
    Matrix<bool> burned_bin(n_col, n_row);
    char* cell_states_initial_h = new char[n_col * n_row];
    char* cell_states_final_h = new char[n_col * n_row];
    bool h_burned = true;
    unsigned int h_burned_size = 0;
    std::vector<std::pair<unsigned int, unsigned int>> h_burned_ids;
    std::vector<unsigned int> h_burned_ids_steps;
    size_t landscape_size = n_col * n_row * sizeof(Cell);

    // Device variables
    bool burned_d;
    char* cell_states_initial_d = nullptr;
    char* cell_states_final_d = nullptr;
    unsigned int d_burned_size = nullptr;
    unsigned short* d_burned_ids = nullptr;
    unsigned int* d_burned_ids_steps = nullptr;
    Cell* d_landscape = nullptr;
    curandState* d_states = nullptr;
    SimulationParams* d_params = nullptr;

    // Variables sizes
    size_t d_burned_ids_size = 2 * n_col * n_row * sizeof(unsigned short);
    size_t d_burned_ids_steps_size = (n_col * n_row) * sizeof(unsigned int);

    // Allocate device memory with error checking
    CUDA_CHECK(cudaMalloc(&cell_states_initial_d, n_col * n_row * sizeof(char)));
    CUDA_CHECK(cudaMalloc(&cell_states_final_d, n_col * n_row * sizeof(char)));
    CUDA_CHECK(cudaMalloc(&d_burned_size, sizeof(int)));
    CUDA_CHECK(cudaMalloc(&burned_d, sizeof(bool)));
    CUDA_CHECK(cudaMalloc(&d_burned_ids, d_burned_ids_size));
    CUDA_CHECK(cudaMalloc(&d_burned_ids_steps, d_burned_ids_steps_size));
    CUDA_CHECK(cudaMalloc(&d_landscape, landscape_size));
    CUDA_CHECK(cudaMalloc(&d_states, n_col * n_row * sizeof(curandState)));
    CUDA_CHECK(cudaMalloc(&d_params, sizeof(SimulationParams)));

    // Initialize device memory to zero
    CUDA_CHECK(cudaMemset(d_burned_ids, 0, d_burned_ids_size));
    CUDA_CHECK(cudaMemset(d_burned_ids_steps, 0, d_burned_ids_steps_size));

    // Initialize cell states in HOST
    for (unsigned int i = 0; i < n_col; i++) {
        for (unsigned int j = 0; j < n_row; j++) {
            cell_states_initial_h[j * n_col + i] = 'U'; 
        }
    }

    // Add ignition cells
    for (const auto& cell : ignition_cells) {
        cell_states_initial_h[cell.second * n_col + cell.first] = 'B';
        h_burned_size++;
        h_burned_ids.push_back(cell);
    }
    h_burned_ids_steps.push_back(h_burned_size);

    // Copy burned_ids to device
    short int* burned_ids_temp = new short int[2 * h_burned_size];
    for (unsigned int i = 0; i < h_burned_size; i++) {
        burned_ids_temp[2*i] = h_burned_ids[i].first;
        burned_ids_temp[2*i+1] = h_burned_ids[i].second;
    }
    CUDA_CHECK(cudaMemcpy(&d_burned_ids, burned_ids_temp, 2 * h_burned_size * sizeof(unsigned int), cudaMemcpyHostToDevice));
    //CUDA_CHECK(cudaMemcpy(&d_burned_ids_steps, burned_ids_steps.back(), sizeof(unsigned int), cudaMemcpyHostToDevice));
    delete[] burned_ids_temp;

    // Copy initial cell states to DEVICE    
    CUDA_CHECK(cudaMemcpy(cell_states_initial_d, cell_states_initial_h, n_col * n_row * sizeof(char), cudaMemcpyHostToDevice));

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

    // Initialize random states
    int block_size = 256;
    int num_blocks = (n_col * n_row + block_size - 1) / block_size;
    init_random_states<<<num_blocks, block_size>>>(d_states, time(NULL));
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    unsigned int iteration = 1;

    //short int d_last_burned_size = d_burned_size;
    
    while (h_burned) {

        // Launch kernel
        calculate_spread_probabilities<<<num_blocks, block_size>>>(
            d_landscape, d_burned_ids, d_burned_size,
            n_col, n_row, d_params, distance, elevation_mean, elevation_sd,
            upper_limit, d_states, burned_d, cell_states_initial_d, cell_states_final_d
        );

        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());


        // Add the number of cells that burned in the last iteration to burned_ids_steps
        int old_burned_size = h_burned_size;
        CUDA_CHECK(cudaMemcpy(&h_burned_size, d_burned_size, sizeof(unsigned int), cudaMemcpyDeviceToHost));
        if (h_burned_size == old_burned_size) {
            // No new cells burned, exit loop
            h_burned = false;
            break;
        }
        h_burned_ids_steps.push_back(h_burned_size);
        

        // Swap initial and final states
        char* temp = cell_states_initial_d;
        cell_states_initial_d = cell_states_final_d;
        cell_states_final_d = temp;

        // Copy burned flag back to host
        CUDA_CHECK(cudaMemcpy(&h_burned, burned_d, sizeof(bool), cudaMemcpyDeviceToHost))

        // Agregar celdas incendiadas a steps
        CUDA_CHECK(cudaMemcpy(d_burned_ids_steps + iteration * sizeof(unsigned short), &d_burned_iteration, sizeof(short int), cudaMemcpyDeviceToDevice));

        // PASAR QUEMANDOSE A QUEMADOS CAPAZ

        // // --- DEBUG ---
        // // Get number of new burned cells
        // CUDA_CHECK(cudaMemcpy(&n_new_burned, d_n_new_burned, sizeof(unsigned int), cudaMemcpyDeviceToHost));
        // CUDA_CHECK(cudaDeviceSynchronize());

        // // Debug print
        // fprintf(stderr, "n_new_burned = %u\n", n_new_burned);
        // if (n_new_burned > 2 * n_col * n_row) {
        //     fprintf(stderr, "ERROR: n_new_burned (%u) exceeds max possible cells (%u)!\n", n_new_burned, 2 * n_col * n_row);
        //     exit(1);
        // }
        // // --- END OF DEBUG ---

        ++iteration;
    }

    CUDA_CHECK(cudaMemcpy(&cell_states_final_h, cell_states_final_d, n_col * n_row * sizeof(char), cudaMemcpyDeviceToHost));
    for (unsigned int i = 0; i < n_col * n_row; i++) {
        unsigned int x = i % n_col;
        unsigned int y = i / n_col;
        burned_bin[x][y] = (cell_states_final_h[i] != 'U');
    }

    short int* h_burned_ids_aux = new short int[2 * n_col * n_row];
    CUDA_CHECK(cudaMemcpy(&h_burned_ids_aux, d_burned_ids, d_burned_ids_size, cudaMemcpyDeviceToHost));
    for (unsigned int i = 0; i < 2 * n_col * n_row; i += 2) {
        h_burned_ids.push_back({h_burned_ids_aux[i], h_burned_ids_aux[i + 1]});
    }

    short int* h_burned_ids_steps_aux = new short int[n_col * n_row];
    CUDA_CHECK(cudaMemcpy(&h_burned_ids_steps_aux, d_burned_ids_steps, d_burned_ids_steps_size, cudaMemcpyDeviceToHost));
    for (unsigned int i = 0; i < n_col * n_row; i++)
        h_burned_ids_steps.push_back(h_burned_ids_steps_aux[i]);

    // Free device memory
    if (d_landscape) CUDA_CHECK(cudaFree(d_landscape));
    if (d_states) CUDA_CHECK(cudaFree(d_states));
    if (d_params) CUDA_CHECK(cudaFree(d_params));
    if (cell_states_initial_d) CUDA_CHECK(cudaFree(cell_states_initial_d));
    if (cell_states_final_d) CUDA_CHECK(cudaFree(cell_states_final_d));
    if (burned_d) CUDA_CHECK(cudaFree(burned_d));
    if (d_burned_ids) CUDA_CHECK(cudaFree(d_burned_ids));
    if (d_burned_ids_steps) CUDA_CHECK(cudaFree(d_burned_ids_steps));
    if (d_burned_size) CUDA_CHECK(cudaFree(d_burned_size));
    
    // Get elapsed time using CUDA events
    CUDA_CHECK(cudaEventRecord(stop_event));
    CUDA_CHECK(cudaEventSynchronize(stop_event));
    float milliseconds = 0;
    CUDA_CHECK(cudaEventElapsedTime(&milliseconds, start_event, stop_event));
    double seconds = milliseconds / 1000.0;

    // fprintf(stderr, "Celdas incendiadas: %ld\n", burned_ids.size());
    // fprintf(stderr, "celdas incendiadas por microsegundo: %lf\n", burned_ids.size() / (1E06 * seconds));

    // Clean up CUDA events
    CUDA_CHECK(cudaEventDestroy(start_event));
    CUDA_CHECK(cudaEventDestroy(stop_event));

    delete[] landscape_data;
    delete[] h_burned_ids_aux;
    delete[] h_burned_ids_steps_aux;

    return { n_col, n_row, burned_bin, h_burned_ids, h_burned_ids_steps };
}
