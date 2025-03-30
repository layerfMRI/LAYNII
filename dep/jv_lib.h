#include <stdio.h>
#include <iostream>
#include <limits>
#include <cmath>
#include <tuple>
// #include "./nifti2_io.h"

// ============================================================================
// Command-line log messages
// ============================================================================
void jv_log_welcome(const char* programname);
// void jv_log_nifti_descriptives(nifti_image* nii);

// ============================================================================
// Utility functions
// ============================================================================
float jv_gaussian(float distance, float sigma);

std::tuple<uint32_t, uint32_t, uint32_t> jv_ind2sub_3D(const uint32_t linear_index,
                                                       const uint32_t nx, const uint32_t ny);
std::tuple<uint32_t, uint32_t, uint32_t, uint32_t> jv_ind2sub_4D(const uint32_t linear_index, 
                                                                 const uint32_t nx, const uint32_t ny, const uint32_t nz);

uint32_t jv_sub2ind_3D(const uint32_t x, const uint32_t y, const uint32_t z,
                       const uint32_t nx, const uint32_t ny);
uint32_t jv_sub2ind_4D(const uint32_t x, const uint32_t y, const uint32_t z, const uint32_t t,
                       const uint32_t nx, const uint32_t ny, const uint32_t nz);

// ============================================================================
// Data manipulations
// ============================================================================
void jv_normalize_to_zero_one(float* data, int data_size);

void jv_smooth_gaussian_iterative_3D(float* data_in,
                                     const int   nx, const int   ny, const int   nz, const int nt,
                                     const float dx, const float dy, const float dz,
                                     const float FWHM_val, const int nr_iterations, const bool log = true);

void jv_compute_gradients_3D(const float* data, float* data_grad_x, float* data_grad_y, float* data_grad_z, 
                             const int nx, const int ny, const int nz, const int nt);


void jv_compute_gradients_3D_over_x(const float* data, float* data_out, 
                                    const int nx, const int ny, const int nz, const int nt);
void jv_compute_gradients_3D_over_y(const float* data, float* data_out, 
                                    const int nx, const int ny, const int nz, const int nt);
void jv_compute_gradients_3D_over_z(const float* data, float* data_out, 
                                    const int nx, const int ny, const int nz, const int nt);

void jv_compute_hessian_3D(const float* data, float* data_shorthessian,
                           const int nx, const int ny, const int nz, const int nt);

void jv_compute_eigen_values_3D(const float* data_shorthessian, float* data_eigvals,
                                const int nx, const int ny, const int nz, const int nt);

void jv_compute_vesselness_3D_v2(const float* data_eigvals, float* data_out, const int data_size,
                                 const float A, const float B, const float C,
                                 const bool include_dark, const bool include_bright);
