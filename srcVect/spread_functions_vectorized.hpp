#pragma once

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
);