
#include <stdio.h>
//#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <tuple>
#include "./nifti2_io.h"

using namespace std;

// ============================================================================
// Declarations
// ============================================================================
double ren_average(double arr[], int size);
double ren_stdev(double arr[], int size);
double ren_correl(double arr1[], double arr2[], int size);
double ren_skew(double arr[], int size);
double ren_kurt(double arr[], int size);
double ren_autocor(double arr[], int size);
int ren_most_occurred_number(int nums[], int size);
//int ren_add_if_new(int arr[], int size); 

float dist(float x1, float y1, float z1, float x2, float y2, float z2,
           float dX, float dY, float dZ);
float dist2d(float x1, float y1, float x2, float y2);
float angle(float a, float b, float c);
float gaus(float distance, float sigma);

void log_welcome(const char* programname);
void log_output(const char* filename);
void log_nifti_descriptives(nifti_image* nii);

void save_output_nifti(string filename, string prefix, nifti_image* nii,
                       bool log = true, bool use_outpath = false);

nifti_image* copy_nifti_as_double(nifti_image* nii);
nifti_image* copy_nifti_as_float32(nifti_image* nii);
nifti_image* copy_nifti_as_float16(nifti_image* nii);
nifti_image* copy_nifti_as_int32(nifti_image* nii);
nifti_image* copy_nifti_as_int16(nifti_image* nii);
nifti_image* copy_nifti_as_float32_with_scl_slope_and_scl_inter(nifti_image* nii);

std::tuple<uint32_t, uint32_t, uint32_t> ind2sub_3D(
    const uint32_t linear_index,
    const uint32_t size_x,
    const uint32_t size_y);

std::tuple<uint32_t, uint32_t, uint32_t, uint32_t> ind2sub_4D(
    const uint32_t linear_index,
    const uint32_t size_x,
    const uint32_t size_y,
    const uint32_t size_z);

uint32_t sub2ind_3D(const uint32_t x,
                    const uint32_t y,
                    const uint32_t z,
                    const uint32_t size_x,
                    const uint32_t size_y);

uint32_t sub2ind_4D(const uint32_t x,
                    const uint32_t y,
                    const uint32_t z,
                    const uint32_t t,
                    const uint32_t size_x,
                    const uint32_t size_y,
                    const uint32_t size_z);

std::tuple<float, float> simplex_closure_2D(float x, float y);
std::tuple<float, float> simplex_perturb_2D(float x, float y, float a, float b);

nifti_image* iterative_smoothing(nifti_image* nii_in, int iter_smooth,
                                 nifti_image* nii_mask, int32_t mask_value);

// ============================================================================
// Preprocessor macros.
// ============================================================================
#define PI 3.14159265;
