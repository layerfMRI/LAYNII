
#include "../src/common.h"
#include "../src/utils.h"

// Definitions
void log_output(const char* filename) {
    cout << "  Writing output as:" << endl;
    cout << "    " << filename << endl;
}

// float* typecast_voxels_float_to_float(float* nim_in, int size_x, int size_y,
//     int size_z, int size_t) {
//
//     float* nim_new = NULL;
//     int nx = size_x;
//     int nxy = size_x * size_y;
//     int nxyz = size_x * size_y * size_z;
//
//     for (int t = 0; t < size_t; ++t) {  // time
//         for (int z = 0; z < size_z; ++z) {  // slice
//             for (int y = 0; y < size_y; ++y) {  // phase
//                 for (int x = 0; x < size_x; ++x) {  // read
//                     *(nim_new + nxyz * t + nxy * z + nx * x + y) =
//                         static_cast<float>(*(nim_in + nxyz * t
//                             + nxy * z + nx * x + y));
//                 }
//             }
//         }
//     }
//     return nim_new;
// }
//
// float* typecast_voxels_int16_to_float(int16_t* nim_in, int size_x, int size_y,
//     int size_z, int size_t) {
//
//     float* nim_new = NULL;
//     int nx = size_x;
//     int nxy = size_x * size_y;
//     int nxyz = size_x * size_y * size_z;
//
//     for (int t = 0; t < size_t; ++t) {  // time
//         for (int z = 0; z < size_z; ++z) {  // slice
//             for (int y = 0; y < size_y; ++y) {  // phase
//                 for (int x = 0; x < size_x; ++x) {  // read
//                     *(nim_new + nxyz * t + nxy * z + nx * x + y) =
//                         static_cast<float>(*(nim_in + nxyz * t
//                             + nxy * z + nx * x + y));
//                 }
//             }
//         }
//     }
//     return nim_new;
// }
