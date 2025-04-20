#include <immintrin.h>
#include <cmath>
#include <random>
#include <vector>
#include <omp.h>
#include <cstring>

// Helper function for AVX2 exponential approximation
// This is a simple implementation and could be improved for better accuracy
inline __m256 _mm256_exp_ps(__m256 x) {
    // Approximate exp(x) using polynomial approximation
    // exp(x) ≈ 1 + x + x^2/2 + x^3/6 + x^4/24 for small x
    __m256 exp_x = x;
    __m256 exp_x_squared = _mm256_mul_ps(exp_x, exp_x);
    __m256 exp_x_cubed = _mm256_mul_ps(exp_x_squared, exp_x);
    __m256 exp_x_fourth = _mm256_mul_ps(exp_x_cubed, exp_x);
    
    __m256 exp_term1 = exp_x;
    __m256 exp_term2 = _mm256_div_ps(exp_x_squared, _mm256_set1_ps(2.0f));
    __m256 exp_term3 = _mm256_div_ps(exp_x_cubed, _mm256_set1_ps(6.0f));
    __m256 exp_term4 = _mm256_div_ps(exp_x_fourth, _mm256_set1_ps(24.0f));
    
    return _mm256_add_ps(_mm256_set1_ps(1.0f), 
                        _mm256_add_ps(exp_term1, 
                                     _mm256_add_ps(exp_term2, 
                                                  _mm256_add_ps(exp_term3, exp_term4))));
}

// Vectorized version using AVX2
// This function processes 8 elements at once using AVX2
#pragma omp declare simd
void spread_probability_vectorized(
    const float* burning_cell_elevation, 
    const float* burning_cell_wind_direction, 
    const float* neighbour_fwi, 
    const float* neighbour_aspect, 
    const float* neighbour_elevation,
    const float* linpred,
    float fwi_pred, 
    float aspect_pred, 
    float wind_pred,
    float elevation_pred, 
    float slope_pred, 
    const float* angle, 
    float distance, 
    float elevation_mean,
    float elevation_sd, 
    float upper_limit,
    float* result,
    int size
) {
    // Process 8 elements at a time (AVX2 with float)
    #pragma omp simd
    for (int i = 0; i < size; i += 8) {
        // Load 8 elements into AVX2 registers
        __m256 burning_elev = _mm256_loadu_ps(&burning_cell_elevation[i]);
        __m256 burning_wind = _mm256_loadu_ps(&burning_cell_wind_direction[i]);
        __m256 neigh_fwi = _mm256_loadu_ps(&neighbour_fwi[i]);
        __m256 neigh_aspect = _mm256_loadu_ps(&neighbour_aspect[i]);
        __m256 neigh_elev = _mm256_loadu_ps(&neighbour_elevation[i]);
        __m256 lin_pred = _mm256_loadu_ps(&linpred[i]);
        __m256 angles = _mm256_loadu_ps(&angle[i]);
        
        // Constants
        __m256 dist = _mm256_set1_ps(distance);
        __m256 elev_mean = _mm256_set1_ps(elevation_mean);
        __m256 elev_sd = _mm256_set1_ps(elevation_sd);
        __m256 fwi_p = _mm256_set1_ps(fwi_pred);
        __m256 aspect_p = _mm256_set1_ps(aspect_pred);
        __m256 wind_p = _mm256_set1_ps(wind_pred);
        __m256 elev_p = _mm256_set1_ps(elevation_pred);
        __m256 slope_p = _mm256_set1_ps(slope_pred);
        __m256 upper = _mm256_set1_ps(upper_limit);
        __m256 one = _mm256_set1_ps(1.0f);
        
        // Calculate slope_term = sin(atan((neighbour_elevation - burning_cell_elevation) / distance))
        __m256 elev_diff = _mm256_sub_ps(neigh_elev, burning_elev);
        __m256 elev_diff_div_dist = _mm256_div_ps(elev_diff, dist);
        
        // Use AVX2 intrinsics for vectorized math operations
        // For atan and sin, we'll use the AVX2 approximations
        
        // Approximate atan(x) using polynomial approximation
        // atan(x) ≈ x * (0.9969 - 0.3173 * x^2) for |x| < 1
        // For larger values, we use atan(x) ≈ π/2 - 1/x * (0.9969 - 0.3173 / x^2)
        __m256 abs_x = _mm256_and_ps(elev_diff_div_dist, _mm256_castsi256_ps(_mm256_set1_epi32(0x7FFFFFFF)));
        __m256 x_squared = _mm256_mul_ps(elev_diff_div_dist, elev_diff_div_dist);
        __m256 coef1 = _mm256_set1_ps(0.9969f);
        __m256 coef2 = _mm256_set1_ps(0.3173f);
        __m256 term = _mm256_sub_ps(coef1, _mm256_mul_ps(coef2, x_squared));
        __m256 atan_result = _mm256_mul_ps(elev_diff_div_dist, term);
        
        // Approximate sin(x) using polynomial approximation
        // sin(x) ≈ x - x^3/6 + x^5/120 for small x
        __m256 x_cubed = _mm256_mul_ps(atan_result, _mm256_mul_ps(atan_result, atan_result));
        __m256 x_fifth = _mm256_mul_ps(x_cubed, _mm256_mul_ps(atan_result, atan_result));
        __m256 term1 = _mm256_div_ps(x_cubed, _mm256_set1_ps(6.0f));
        __m256 term2 = _mm256_div_ps(x_fifth, _mm256_set1_ps(120.0f));
        __m256 slope_term = _mm256_sub_ps(atan_result, _mm256_sub_ps(term1, term2));
        
        // Calculate wind_term = cos(angle - burning_cell_wind_direction)
        __m256 angle_diff = _mm256_sub_ps(angles, burning_wind);
        
        // Approximate cos(x) using polynomial approximation
        // cos(x) ≈ 1 - x^2/2 + x^4/24 for small x
        // First, normalize angle_diff to [-π, π]
        __m256 pi = _mm256_set1_ps(3.14159265358979323846f);
        __m256 two_pi = _mm256_set1_ps(2.0f * 3.14159265358979323846f);
        __m256 normalized_angle = angle_diff;
        normalized_angle = _mm256_sub_ps(normalized_angle, 
                                        _mm256_mul_ps(_mm256_round_ps(_mm256_div_ps(normalized_angle, two_pi), _MM_FROUND_TO_NEAREST_INT), 
                                                     two_pi));
        
        __m256 angle_squared = _mm256_mul_ps(normalized_angle, normalized_angle);
        __m256 angle_fourth = _mm256_mul_ps(angle_squared, angle_squared);
        __m256 cos_term1 = _mm256_div_ps(angle_squared, _mm256_set1_ps(2.0f));
        __m256 cos_term2 = _mm256_div_ps(angle_fourth, _mm256_set1_ps(24.0f));
        __m256 wind_term = _mm256_sub_ps(_mm256_set1_ps(1.0f), _mm256_sub_ps(cos_term1, cos_term2));
        
        // Calculate elev_term = (neighbour_elevation - elevation_mean) / elevation_sd
        __m256 elev_term = _mm256_div_ps(_mm256_sub_ps(neigh_elev, elev_mean), elev_sd);
        
        // Update linpred
        lin_pred = _mm256_add_ps(lin_pred, _mm256_mul_ps(fwi_p, neigh_fwi));
        lin_pred = _mm256_add_ps(lin_pred, _mm256_mul_ps(aspect_p, neigh_aspect));
        
        __m256 wind_term_contrib = _mm256_mul_ps(wind_term, wind_p);
        __m256 elev_term_contrib = _mm256_mul_ps(elev_term, elev_p);
        __m256 slope_term_contrib = _mm256_mul_ps(slope_term, slope_p);
        
        lin_pred = _mm256_add_ps(lin_pred, wind_term_contrib);
        lin_pred = _mm256_add_ps(lin_pred, elev_term_contrib);
        lin_pred = _mm256_add_ps(lin_pred, slope_term_contrib);
        
        // Calculate prob = upper_limit / (1 + exp(-linpred))
        __m256 neg_lin_pred = _mm256_sub_ps(_mm256_setzero_ps(), lin_pred);
        
        // Approximate exp(x) using polynomial approximation
        __m256 exp_result = _mm256_exp_ps(neg_lin_pred);
        
        // Clamp exp_result to avoid overflow
        __m256 max_val = _mm256_set1_ps(20.0f);  // Limit to a reasonable range
        exp_result = _mm256_min_ps(exp_result, _mm256_exp_ps(max_val));
        
        __m256 denominator = _mm256_add_ps(one, exp_result);
        __m256 prob = _mm256_div_ps(upper, denominator);
        
        // Store the result
        _mm256_storeu_ps(&result[i], prob);
    }

}
