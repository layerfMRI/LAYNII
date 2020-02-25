#pragma once

#include <stdio.h>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <tuple>
#include "./nifti2_io.h"

using namespace std;

// ============================================================================
// ============================================================================
// ============================================================================
// Declarations
double ren_average(double arr[], int size);
double ren_stdev(double arr[], int size);
double ren_correl(double arr1[], double arr2[], int size);
double ren_skew(double arr[], int size);
double ren_kurt(double arr[], int size);
double ren_autocor(double arr[], int size);

void log_welcome(const char* programname);
void log_output(const char* filename);
void log_nifti_descriptives(nifti_image* nii);

nifti_image* copy_nifti_header_as_float(nifti_image* nii);
nifti_image* copy_nifti_header_as_int(nifti_image* nii);
nifti_image* copy_nifti_header_as_uint(nifti_image* nii);

std::tuple<int, int, int> ind2sub_3D(const int linear_index, const int size_x,
                                     const int size_y);

int sub2ind_3D(const int x, const int y, const int z,
               const int size_x, const int size_y);

// ============================================================================
// ============================================================================
// ============================================================================
// Preprocessor macros.
#define PI 3.14159265;

// NOTE: Do not put any characters after `\` for multiline defines to work.
// Loop voxels 4D (t=time steps, z=slices, y=read steps, x=phase steps)
#define FOR_EACH_VOXEL_TZYX for (int it = 0; it < size_t; ++it) {\
                            for (int iz = 0; iz < size_z; ++iz) {\
                            for (int iy = 0; iy < size_y; ++iy) {\
                            for (int ix = 0; ix < size_x; ++ix) {

#define END_FOR_EACH_VOXEL_TZYX }}}}

// Used inside voxel loops
#define VOXEL_ID_3D (nxy * iz + nx * iy + ix)
#define VOXEL_ID (nxyz * it + nxy * iz + nx * iy + ix)
