#include <stdio.h>
#include <iostream>
#include <limits>
#include <cmath>
#include <tuple>

// ============================================================================
// Command-line log messages
// ============================================================================
void ida_log_welcome(const char* programname);

// ============================================================================
// Utility functions
// ============================================================================
std::tuple<uint64_t, uint64_t, uint64_t, uint64_t> ida_ind2sub_4D(
    const uint64_t linear_index, const uint64_t nx, const uint64_t ny, const uint64_t nz);

std::tuple<uint64_t, uint64_t, uint64_t, uint64_t> ida_ind2sub_4D_Tmajor(
    const uint64_t linear_index,const uint64_t nt, const uint64_t nx, const uint64_t ny);

uint64_t ida_sub2ind_4D(const uint64_t x, const uint64_t y, const uint64_t z, const uint64_t t,
                        const uint64_t nx, const uint64_t ny, const uint64_t nz);

uint64_t ida_sub2ind_4D_Tmajor(const uint64_t t, const uint64_t x, const uint64_t y, const uint64_t z,
                               const uint64_t nt, const uint64_t nx, const uint64_t ny);

