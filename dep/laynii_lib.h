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

float dist(float x1, float y1, float z1, float x2, float y2, float z2,
           float dX, float dY, float dZ);
float angle(float a, float b, float c);
float gaus(float distance, float sigma);

void log_welcome(const char* programname);
void log_output(const char* filename);
void log_nifti_descriptives(nifti_image* nii);

void save_output_nifti(string filename, string prefix, nifti_image* nii,
                       bool log = true);

nifti_image* copy_nifti_as_float32(nifti_image* nii);
nifti_image* copy_nifti_as_int32(nifti_image* nii);

std::tuple<uint32_t, uint32_t, uint32_t> ind2sub_3D(
    const uint32_t linear_index, const uint32_t size_x, const uint32_t size_y);

uint32_t sub2ind_3D(const uint32_t x, const uint32_t y, const uint32_t z,
                    const uint32_t size_x, const uint32_t size_y);

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
