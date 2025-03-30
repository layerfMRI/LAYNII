
#include "./jv_lib.h"

// ====================================================================================================================
// Command-line log messages
// ====================================================================================================================

void jv_log_welcome(const char* program_name) {
    // Print the main program name, version number, and sub-program name
    std::printf("========================\n");
    std::printf("JV Lib v0.0.0           \n");
    std::printf("========================\n");
    std::printf("%s\n", program_name);
    std::printf("\n");
}


// void jv_log_nifti_descriptives(nifti_image* nii) {
//     // Print nifti descriptives to command line
//     std::printf("  Input information:\n");
//     std::printf("    File name  = %s\n", nii->fname);
//     std::printf("    Image size = [%lli, %lli, %lli, %lli] [nx, ny, nz, nt]\n", nii->nx, nii->ny, nii->nz, nii->nt);
//     std::printf("    Voxel size = [%f, %f, %f] [dx, dy, dz]\n", nii->pixdim[1], nii->pixdim[2], nii->pixdim[3]);
//     std::printf("    Datatype   = %d\n", nii->datatype);
//     std::printf("\n");
// }

// ====================================================================================================================
// Utility functions
// ====================================================================================================================

float jv_gaussian(float distance, float sigma) {
    return 1. / (sigma * std::sqrt(2. * 3.141592))
           * std::exp(-0.5 * distance * distance / (sigma * sigma));
}

std::tuple<uint32_t, uint32_t, uint32_t> jv_ind2sub_3D(
    const uint32_t linear_index,
    const uint32_t size_x,
    const uint32_t size_y) {

    uint32_t z = linear_index / (size_x * size_y);
    uint32_t temp = linear_index % (size_x * size_y);
    uint32_t y = temp / size_x;
    uint32_t x = temp % size_x;

    return std::make_tuple(x, y, z);
}

std::tuple<uint32_t, uint32_t, uint32_t, uint32_t> jv_ind2sub_4D(
    const uint32_t linear_index,
    const uint32_t size_x,
    const uint32_t size_y,
    const uint32_t size_z) {

    uint32_t t = linear_index / (size_x * size_y * size_z);
    uint32_t temp = linear_index % (size_x * size_y * size_z);
    uint32_t z = temp / (size_x * size_y);
    temp = linear_index % (size_x * size_y);
    uint32_t y = temp / size_x;
    uint32_t x = temp % size_x;

    return std::make_tuple(x, y, z, t);
}


uint32_t jv_sub2ind_3D(const uint32_t x, const uint32_t y, const uint32_t z,
                    const uint32_t size_x, const uint32_t size_y) {
    return size_x * size_y * z + size_x * y + x;
}

uint32_t jv_sub2ind_4D(const uint32_t x, const uint32_t y, const uint32_t z, const uint32_t t,
                    const uint32_t size_x, const uint32_t size_y, const uint32_t size_z) {
    return size_x * size_y * size_z * t + size_x * size_y * z + size_x * y + x;
}


// ====================================================================================================================
// Data manipulations
// ====================================================================================================================

void jv_normalize_to_zero_one(float* data, int data_size) {

    // Find minimum and maximum
    float temp_max = std::numeric_limits<float>::min();
    float temp_min = std::numeric_limits<float>::max();
    for (int i = 0; i != data_size; ++i) {
        if (*(data + i) > temp_max) {
            temp_max = *(data + i);
        }
        if (*(data + i) < temp_min) {
            temp_min = *(data + i);
        }
    }

    // Translate minimum to zero
    for (int i = 0; i != data_size; ++i) {
        *(data + i) -= temp_min;
    }

    // Scale maximum to one
    temp_max -= temp_min;
    for (int i = 0; i != data_size; ++i) {
        *(data + i) /= temp_max;
    }    
}


void jv_smooth_gaussian_iterative_3D(float* data_in,
                                     const int   nx, const int   ny, const int   nz, const int nt,
                                     const float dx, const float dy, const float dz,
                                     const float fwhm, const int nr_iterations, const bool log) {
    // NOTE: Overwrites the input with the smoothed image at the end

    int data_size = nx * ny * nz * nt;

    float* data_temp = (float*)malloc(data_size * sizeof(float));

    // Compute Gaussian weights
    float w_0 = jv_gaussian(0, fwhm);
    float w_dX = jv_gaussian(dx, fwhm);
    float w_dY = jv_gaussian(dy, fwhm);
    float w_dZ = jv_gaussian(dz, fwhm);

    // Loop over every data point
    int ix, iy, iz, it, j;
    for (int n = 0; n != nr_iterations; ++n) {

        if (log) std::printf("    Iteration: %i/%i\n", n+1, nr_iterations);

        for (int i = 0; i != data_size; ++i) {

            std::tie(ix, iy, iz, it) = jv_ind2sub_4D(i, nx, ny, nz);
            float new_val = 0, total_weight = 0;

            // Start with the voxel itself
            new_val += *(data_in + i) * w_0;
            total_weight += w_0;

            // ----------------------------------------------------------------
            // 1-jump neighbours
            // ----------------------------------------------------------------
            if (ix > 0) {
                j = jv_sub2ind_4D(ix-1, iy, iz, it, nx, ny, nz);
                new_val += *(data_in + j) * w_dX;
                total_weight += w_dX;
            }
            if (ix < nx-1) {
                j = jv_sub2ind_4D(ix+1, iy, iz, it, nx, ny, nz);
                new_val += *(data_in + j) * w_dX;
                total_weight += w_dX;
            }
            if (iy > 0) {
                j = jv_sub2ind_4D(ix, iy-1, iz, it, nx, ny, nz);
                new_val += *(data_in + j) * w_dY;
                total_weight += w_dY;
            }
            if (iy < ny-1) {
                j = jv_sub2ind_4D(ix, iy+1, iz, it, nx, ny, nz);
                new_val += *(data_in + j) * w_dY;
                total_weight += w_dY;
            }
            if (iz > 0) {
                j = jv_sub2ind_4D(ix, iy, iz-1, it, nx, ny, nz);
                new_val += *(data_in + j) * w_dZ;
                total_weight += w_dZ;
            }
            if (iz < nz-1) {
                j = jv_sub2ind_4D(ix, iy, iz+1, it, nx, ny, nz);
                new_val += *(data_in + j) * w_dZ;
                total_weight += w_dZ;
            }

            // ----------------------------------------------------------------
            *(data_temp + i) = new_val / total_weight;
        }

        // Swap image data for the next iteration
        for (int i = 0; i != data_size; ++i) {
            *(data_in + i) = *(data_temp + i);
        }
    }
    free(data_temp);
}


void jv_compute_gradients_3D(const float* data_in, float* data_grad_x, float* data_grad_y, float* data_grad_z, 
                             const int nx, const int ny, const int nz, const int nt) {

    int data_size = nx * ny * nz * nt;

    // Loop over every data point
    int ix, iy, iz, it, j, k;
    for (int i = 0; i != data_size; ++i) {
        std::tie(ix, iy, iz, it) = jv_ind2sub_4D(i, nx, ny, nz);

        // --------------------------------------------------------------------
        // 1-jump neighbours
        // --------------------------------------------------------------------
        if (ix > 0 && ix < nx-1) {
            j = jv_sub2ind_4D(ix-1, iy, iz, it, nx, ny, nz);
            k = jv_sub2ind_4D(ix+1, iy, iz, it, nx, ny, nz);
            *(data_grad_x + i) = *(data_in + j) - *(data_in + k);
        }
        if (iy > 0 && iy < ny-1) {
            j = jv_sub2ind_4D(ix, iy-1, iz, it, nx, ny, nz);
            k = jv_sub2ind_4D(ix, iy+1, iz, it, nx, ny, nz);
            *(data_grad_y + i)  = *(data_in + j) - *(data_in + k);
        }
        if (iz > 0 && iz < nz-1) {
            j = jv_sub2ind_4D(ix, iy, iz-1, it, nx, ny, nz);
            k = jv_sub2ind_4D(ix, iy, iz+1, it, nx, ny, nz);
            *(data_grad_z + i)  = *(data_in + j) - *(data_in + k);
        }
    }
}


void jv_compute_gradients_3D_over_x(const float* data_in, float* data_out, 
                                    const int nx, const int ny, const int nz, const int nt) {

    int data_size = nx * ny * nz * nt;

    // Loop over every data point
    int ix, iy, iz, it, j, k;
    for (int i = 0; i != data_size; ++i) {
        std::tie(ix, iy, iz, it) = jv_ind2sub_4D(i, nx, ny, nz);
        if (ix > 0 && ix < nx-1) {
            j = jv_sub2ind_4D(ix-1, iy, iz, it, nx, ny, nz);
            k = jv_sub2ind_4D(ix+1, iy, iz, it, nx, ny, nz);
            *(data_out + i) = *(data_in + j) - *(data_in + k);
        }
    }
}


void jv_compute_gradients_3D_over_y(const float* data_in, float* data_out, 
                                    const int nx, const int ny, const int nz, const int nt) {

    int data_size = nx * ny * nz * nt;

    // Loop over every data point
    int ix, iy, iz, it, j, k;
    for (int i = 0; i != data_size; ++i) {
        std::tie(ix, iy, iz, it) = jv_ind2sub_4D(i, nx, ny, nz);
        if (iy > 0 && iy < ny-1) {
            j = jv_sub2ind_4D(ix, iy-1, iz, it, nx, ny, nz);
            k = jv_sub2ind_4D(ix, iy+1, iz, it, nx, ny, nz);
            *(data_out + i)  = *(data_in + j) - *(data_in + k);
        }
    }
}


void jv_compute_gradients_3D_over_z(const float* data_in, float* data_out, 
                                    const int nx, const int ny, const int nz, const int nt) {

    int data_size = nx * ny * nz * nt;

    // Loop over every data point
    int ix, iy, iz, it, j, k;
    for (int i = 0; i != data_size; ++i) {
        std::tie(ix, iy, iz, it) = jv_ind2sub_4D(i, nx, ny, nz);
        if (iz > 0 && iz < nz-1) {
            j = jv_sub2ind_4D(ix, iy, iz-1, it, nx, ny, nz);
            k = jv_sub2ind_4D(ix, iy, iz+1, it, nx, ny, nz);
            *(data_out + i)  = *(data_in + j) - *(data_in + k);
        }
    }
}


void jv_compute_hessian_3D(const float* data_in, float* data_shorthessian,
                           const int nx, const int ny, const int nz, const int nt) {
    // NOTE: Hessian data should be 6 times larger than the input
    // NOTE: Hessian values are saved consecutively for each voxel
    //     *(data_shorthessian + i*6 + 0) = 2nd derivative xx
    //     *(data_shorthessian + i*6 + 1) = 2nd derivative xy yx
    //     *(data_shorthessian + i*6 + 2) = 2nd derivative xz zx
    //     *(data_shorthessian + i*6 + 3) = 2nd derivative yy
    //     *(data_shorthessian + i*6 + 4) = 2nd derivative yz zy
    //     *(data_shorthessian + i*6 + 5) = 2nd derivative zz


    int data_size = nx * ny * nz * nt;

    // Allocate memory (NOTE: I have prioritized RAM optimization)
    float* data_grad_1st = (float*)malloc(data_size * sizeof(float));
    float* data_grad_2nd = (float*)malloc(data_size * sizeof(float));

    // x 
    jv_compute_gradients_3D_over_x(data_in, data_grad_1st, nx, ny, nz, nt);

    // xx
    jv_compute_gradients_3D_over_x(data_grad_1st, data_grad_2nd, nx, ny, nz, nt);
    for (uint32_t i = 0; i != data_size; ++i) {
        *(data_shorthessian + i*6 + 0) =  *(data_grad_2nd + i);
    }

    // xy
    jv_compute_gradients_3D_over_y(data_grad_1st, data_grad_2nd, nx, ny, nz, nt);
    for (uint32_t i = 0; i != data_size; ++i) {
        *(data_shorthessian + i*6 + 1) =  *(data_grad_2nd + i);
    }

    // xz
    jv_compute_gradients_3D_over_z(data_grad_1st, data_grad_2nd, nx, ny, nz, nt);
    for (uint32_t i = 0; i != data_size; ++i) {
        *(data_shorthessian + i*6 + 2) =  *(data_grad_2nd + i);
    }

    // y
    jv_compute_gradients_3D_over_y(data_in, data_grad_1st, nx, ny, nz, nt);

    // yy
    jv_compute_gradients_3D_over_y(data_grad_1st, data_grad_2nd, nx, ny, nz, nt);
    for (uint32_t i = 0; i != data_size; ++i) {
        *(data_shorthessian + i*6 + 3) =  *(data_grad_2nd + i);
    }

    // yz
    jv_compute_gradients_3D_over_z(data_grad_1st, data_grad_2nd, nx, ny, nz, nt);
    for (uint32_t i = 0; i != data_size; ++i) {
        *(data_shorthessian + i*6 + 4) =  *(data_grad_2nd + i);
    }

    // z
    jv_compute_gradients_3D_over_z(data_in, data_grad_1st, nx, ny, nz, nt);

    // zz
    jv_compute_gradients_3D_over_z(data_grad_1st, data_grad_2nd, nx, ny, nz, nt);
    for (uint32_t i = 0; i != data_size; ++i) {
        *(data_shorthessian + i*6 + 5) =  *(data_grad_2nd + i);
    } 

    free(data_grad_1st);
    free(data_grad_2nd);
}

void jv_compute_eigen_values_3D(const float* data_hessian, float* data_eigvals,
                                const int nx, const int ny, const int nz, const int nt) {
    // NOTE: Implementing Delledalle et al. 2017, Hal.
    // NOTE: I simplified complex conjugates as I do not have complex values.

    int data_size = nx * ny * nz * nt;

    for (uint32_t i = 0; i != data_size; ++i) {
        float a = *(data_hessian + i*6 + 0);  // xx
        float b = *(data_hessian + i*6 + 3);  // yy
        float c = *(data_hessian + i*6 + 5);  // zz
        float d = *(data_hessian + i*6 + 1);  // xy, yx
        float e = *(data_hessian + i*6 + 4);  // yz, zy
        float f = *(data_hessian + i*6 + 2);  // xz, zx

        // --------------------------------------------------------------------
        // Compute eigen values 3D
        // --------------------------------------------------------------------
        float x1 = a*a + b*b + c*c - a*b - a*c - b*c + 3 * (d*d + f*f + e*e);

        float t1 = 2*a - b - c;
        float t2 = 2*b - a - c;
        float t3 = 2*c - a - b;
        float x2 = - t1 * t2 * t3 + 9*( t3*(d*d) + t2*(f*f) + t1*(e*e) ) - 54*( d*e*f );

        float phi;
        if (x2 > 0) {
            phi = std::atan2( 0, std::sqrt( 4.*( (x1*x1)*x1 ) - x2*x2 ) / x2 );
        }   else if (x2 < 0) {
            phi = std::atan2( 0, std::sqrt( 4.*( (x1*x1)*x1 ) - x2*x2 ) / x2 ) + M_PI;
        } else {
            phi = M_PI / 2.;
        }

        float lambda1 = ( a + b + c - 2.*std::sqrt( x1 ) * std::cos( phi/3. ) ) / 3.;
        float lambda2 = ( a + b + c + 2.*std::sqrt( x1 ) * std::cos( (phi - M_PI)/3. ) ) / 3.;
        float lambda3 = ( a + b + c + 2.*std::sqrt( x1 ) * std::cos( (phi + M_PI)/3. ) ) / 3.;

        if (std::isnan(lambda1)) {
            lambda1 = 0.;
        }
        if (std::isnan(lambda2)) {
            lambda2 = 0.;
        }
        if (std::isnan(lambda3)) {
            lambda3 = 0.;
        }

        // --------------------------------------------------------------------
        // Sort by magnitude while preserving sign of the eigenvalues
        // --------------------------------------------------------------------
        bool lambda1_sign = std::signbit(lambda1);
        bool lambda2_sign = std::signbit(lambda2);
        bool lambda3_sign = std::signbit(lambda3);

        float lambda1_abs = std::abs(lambda1);
        float lambda2_abs = std::abs(lambda2);
        float lambda3_abs = std::abs(lambda3);

        *(data_eigvals + i*3) = std::min( std::min( lambda1_abs, lambda2_abs ), lambda3_abs );
        *(data_eigvals + i*3+2) = std::max( std::max( lambda1_abs, lambda2_abs ), lambda3_abs );
        *(data_eigvals + i*3+1) = lambda1_abs + lambda2_abs + lambda3_abs - ( *(data_eigvals + i*3) + *(data_eigvals + i*3+2) );
 
        // Put back the signbit
        if (*(data_eigvals + i*3) == lambda1_abs) {
            *(data_eigvals + i*3) *= static_cast<float>(lambda1_sign) * 2 - 1;
        } else if (*(data_eigvals + i*3) == lambda2_abs) {
            *(data_eigvals + i*3) *= static_cast<float>(lambda2_sign) * 2 - 1;
        } else if (*(data_eigvals + i*3) == lambda3_abs) {
            *(data_eigvals + i*3) *= static_cast<float>(lambda3_sign) * 2 - 1;
        }

        if (*(data_eigvals + i*3+1) == lambda1_abs) {
            *(data_eigvals + i*3+1) *= static_cast<float>(lambda1_sign) * 2 - 1;
        } else if (*(data_eigvals + i*3+1) == lambda2_abs) {
            *(data_eigvals + i*3+1) *= static_cast<float>(lambda2_sign) * 2 - 1;
        } else if (*(data_eigvals + i*3+1) == lambda3_abs) {
            *(data_eigvals + i*3+1) *= static_cast<float>(lambda3_sign) * 2 - 1;
        }

        if (*(data_eigvals + i*3+2) == lambda1_abs) {
            *(data_eigvals + i*3+2) *= static_cast<float>(lambda1_sign) * 2 - 1;
        } else if (*(data_eigvals + i*3+2) == lambda2_abs) {
            *(data_eigvals + i*3+2) *= static_cast<float>(lambda2_sign) * 2 - 1;
        } else if (*(data_eigvals + i*3+2) == lambda3_abs) {
            *(data_eigvals + i*3+2) *= static_cast<float>(lambda3_sign) * 2 - 1;
        }
    }
}


void jv_compute_vesselness_3D_v2(const float* data_eigvals, float* data_out, const int data_size,
                                 const float A, const float B, const float C,
                                 const bool include_dark, const bool include_bright) {

    for (uint32_t i = 0; i != data_size; ++i) {

        float lambda1 = *(data_eigvals + i*3  );
        float lambda2 = *(data_eigvals + i*3+1);
        float lambda3 = *(data_eigvals + i*3+2);

        // NOTE: 3D Eigen value sign interpretations from Frangi et al. 1998
        // | e1 | e2 | e3 |      Orientation Pattern      |
        // |----|----|----|-------------------------------|
        // | N  | N  | N  | Noisy, no preferred direction |
        // | L  | L  | H- | Plate-like structure, bright  |
        // | L  | L  | H+ | Plate-like structure, dark    |
        // | L  | H- | H- | Tubular structure, bright     |
        // | L  | H+ | H+ | Tubular structure, dark       |
        // | H- | H- | H- | Blob-like structure, bright   |
        // | H+ | H+ | H+ | Blob-like structure, dark     |

         *(data_out + i) = 0;
        float lambda1_abs = std::abs(lambda1);
        float lambda2_abs = std::abs(lambda2);
        float lambda3_abs = std::abs(lambda3);

        float RA = lambda2_abs / lambda3_abs;
        float RB = lambda1_abs / std::sqrt(lambda2_abs*lambda3_abs);
        float S =  std::sqrt( lambda1*lambda1 + lambda2*lambda2 + lambda3*lambda3 );

        float t1 = 1 - std::exp( -( RA*RA ) / ( 2*(A*A) ) );
        float t2 =     std::exp( -( RB*RB ) / ( 2*(B*B) ) );
        float t3 = 1 - std::exp( -(  S*S  ) / ( 2*(C*C) ) );

        *(data_out + i) = t1 * t2 * t3;  

        // NOTE: I could not figure out a way to not compute vesselness based on inclusion criteria. This part can be
        // improved in the future.
        if (include_dark == false) {
            if (lambda2 < 0 || lambda3 < 0) {
                *(data_out + i) = 0;
            }
        }
        if (include_bright == false) {
            if (lambda2 > 0 || lambda3 > 0) {
                *(data_out + i) = 0;
            }            
        }          
    }
}