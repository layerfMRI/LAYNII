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
// Preprocessor macros
// NOTE: Do not put any characters after `\`.
#define FOR_EACH_VOXEL for (int t = 0; t < size_t; t++) {\
                       for (int z = 0; z < size_z; z++) {\
                       for (int y = 0; y < size_y; y++) {\
                       for (int x = 0; x < size_x; x++) {
#define END_FOR_EACH_VOXEL }}}}
