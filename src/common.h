#pragma once

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <string>

#include "../deps/nifti2_io.h"
// #include "../deps/nifti2.h"
// #include "../deps/nifti1.h"
// #include "../deps/nifticdf.h"
// #include "../deps/nifti_tool.h"
// #include< gsl/gsl_multifit.h>
// #include< gsl/gsl_statistics_double.h>

using namespace std;

#define PI 3.14159265;

// #include "utils.hpp"

// ============================================================================
// Preprocessor macros.
// NOTE: Do not put any characters after `\` for multiline defines to work.

// Loop voxels 4D (t=time steps, z=slices, y=read steps, x=phase steps)
#define FOR_EACH_VOXEL_TZYX for (int it = 0; it < size_t; ++it) {\
                            for (int iz = 0; iz < size_z; ++iz) {\
                            for (int iy = 0; iy < size_y; ++iy) {\
                            for (int ix = 0; ix < size_x; ++ix) {

#define END_FOR_EACH_VOXEL_TZYX }}}}

// Used inside voxel loops
#define VOXEL_ID (nxyz * it + nxy * iz + nx * ix + iy)
#define VOXEL_ID_3D (nxy * iz + nx * ix + iy)

// ============================================================================
