#include "./idalib.h"


// ====================================================================================================================
// Command-line log messages
// ====================================================================================================================

void ida_log_welcome(const char* program_name) {
    // Print the main program name, version number, and sub-program name
    std::printf("========================\n");
    std::printf("LayNii IDA    v0.0.0    \n");
    std::printf("========================\n");
    std::printf("%s\n", program_name);
    std::printf("\n");
}

// ====================================================================================================================
// Utility functions
// ====================================================================================================================

std::tuple<uint64_t, uint64_t, uint64_t, uint64_t> ida_ind2sub_4D(
    const uint64_t linear_index,
    const uint64_t size_x,
    const uint64_t size_y,
    const uint64_t size_z) {

    uint64_t size_xy = size_x * size_y;
    uint64_t size_xyz = size_x * size_y * size_z;

    uint64_t t = linear_index / size_xyz;
    uint64_t temp = linear_index % size_xyz;
    uint64_t z = temp / size_xy;
    temp = linear_index % size_xy;
    uint64_t y = temp / size_x;
    uint64_t x = temp % size_x;

    return std::make_tuple(x, y, z, t);
}


std::tuple<uint64_t, uint64_t, uint64_t, uint64_t> ida_ind2sub_4D_Tmajor(
    const uint64_t linear_index,
    const uint64_t size_time,
    const uint64_t size_x,
    const uint64_t size_y) {

    uint64_t size_tx = size_time * size_x;
    uint64_t size_txy = size_time * size_x * size_y;

    uint64_t z = linear_index / size_txy;
    uint64_t rem1 = linear_index % size_txy;
    uint64_t y = rem1 / size_tx;
    uint64_t rem2 = rem1 % size_tx;
    uint64_t x = rem2 / size_time;
    uint64_t t = rem2 % size_time;

    return std::make_tuple(t, x, y, z);
}


uint64_t ida_sub2ind_4D(const uint64_t x, const uint64_t y, const uint64_t z, const uint64_t t,
                       const uint64_t size_x, const uint64_t size_y, const uint64_t size_z) {
    return size_x * size_y * size_z * t + size_x * size_y * z + size_x * y + x;
}

uint64_t ida_sub2ind_4D_Tmajor(const uint64_t t, const uint64_t x, const uint64_t y, const uint64_t z,
                              const uint64_t size_time, const uint64_t size_x, const uint64_t size_y) {
    return t + size_time * (x + size_x * (y + size_y * z));
}

