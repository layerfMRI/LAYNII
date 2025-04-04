#pragma once

#include <iostream>
#include <cstdint>
#include <string>
#include <zlib.h>
#include <chrono>

#include <SDL.h>
#if defined(IMGUI_IMPL_OPENGL_ES2)
#include <SDL_opengles2.h>
#else
#include <SDL_opengl.h>
#endif


namespace IDA_IO
{
    struct NiftiInfo {
        int   sizeof_hdr;       // MUST be 348
        char  data_type[10];    // Unused
        char  db_name[18];      // Unused
        int   extents;          // Unused
        short session_error;    // Unused
        char  regular;          // Unused
        char  dim_info;         // MRI slice ordering
        short dim[8];           // Data array dimensions
        float intent_p1;        // 1st intent parameter
        float intent_p2;        // 2nd intent parameter
        float intent_p3;        // 3rd intent parameter
        short intent_code;      // NIFTI_INTENT code 
        short datatype;         // Defines data type
        short bitpix;           // Number bits/voxel 
        short slice_start;      // First slice index   
        float pixdim[8];        // Grid spacings    
        float vox_offset;       // Offset into .nii file
        float scl_slope;        // Data scaling: slope
        float scl_inter;        // Data scaling: offset
        short slice_end;        // Last slice index.
        char  slice_code;       // Slice timing order.  
        char  xyzt_units;       // Units of pixdim[1..4]
        float cal_max;          // Max display intensity
        float cal_min;          // Min display intensity
        float slice_duration;   // Time for 1 slice.
        float toffset;          // Time axis shift.
        int   glmax;            // Unused
        int   glmin;            // Unused
        char  descrip[80];      // Any text
        char  aux_file[24];     // Auxiliary filename
        short qform_code;       // NIFTI_XFORM code
        short sform_code;       // NIFTI_XFORM code
        float quatern_b;        // Quaternion b param
        float quatern_c;        // Quaternion c param
        float quatern_d;        // Quaternion d param
        float qoffset_x;        // Quaternion x shift
        float qoffset_y;        // Quaternion y shift
        float qoffset_z;        // Quaternion z shift
        float srow_x[4];        // 1st row affine transform
        float srow_y[4];        // 2nd row affine transform
        float srow_z[4];        // 3rd row affine transform
        char  intent_name[16];  // 'name' or meaning of data
        char  magic[4];         // MUST be "ni1\0" or "n+1\0"
        char  extension[4];     // All bytes are 0 if no extension
    };

    struct FileInfo {
        std::string name;
        std::string path;
        NiftiInfo   header;
        int         dim_k;              // File format independent image dimension
        int         dim_j;              // File format independent image dimension
        int         dim_i;              // File format independent image dimension
        int         dim_t;              // File format independent image dimension
        uint64_t    nr_voxels;          // Useful in general
        float       pixdim_i;           // File format independent image dimension
        float       pixdim_j;           // File format independent image dimension
        float       pixdim_k;           // File format independent image dimension
        float       voxel_volume;       // I often need to compute this
        float*      p_data_float;       // holds image data (voxels)
        int         display_k;          // Display slice number
        int         display_j;          // Display slice number
        int         display_i;          // Display slice number
        int         display_t;          // Display slice number
        float*      p_sliceK_float;     // Data slice that holds high precision data
        float*      p_sliceJ_float;     // Data slice that holds high precision data
        float*      p_sliceI_float;     // Data slice that holds high precision data
        uint8_t*    p_sliceK_uint8;     // Data slice that holds display data
        uint8_t*    p_sliceJ_uint8;     // Data slice that holds display data
        uint8_t*    p_sliceI_uint8;     // Data slice that holds display data
        GLuint      textureIDk;         // OpenGL needs this
        GLuint      textureIDj;         // OpenGL needs this
        GLuint      textureIDi;         // OpenGL needs this
        float*      p_time_course_float; // One voxel's time course
        float       time_course_min;     // Minimum of voxel's time course data
        float       time_course_max;     // Maximum of voxel's time course data
        float       data_min;           // Minimum data value
        float       data_max;           // Maximum data value
        float       display_min;        // Minimum displayed value
        float       display_max;        // Maximum displayed value
        float       display_scale;      // Scaling factor for zooming in or out
        float       display_k_offset_x; // Offset displayed image area
        float       display_k_offset_y; // Offset displayed image area
        float       display_j_offset_x; // Offset displayed image area
        float       display_j_offset_y; // Offset displayed image area
        float       display_i_offset_x; // Offset displayed image area
        float       display_i_offset_y; // Offset displayed image area
        // Overlay related --------------------------------------------------------------------------------------------
        float       overlay_min;            // Minimum data value for masking
        float       overlay_max;            // Maximum data value for masking
        uint8_t*    p_sliceK_RGB_uint8;     // RGB data slice that holds display data
        uint8_t*    p_sliceJ_RGB_uint8;     // RGB data slice that holds display data
        uint8_t*    p_sliceI_RGB_uint8;     // RGB data slice that holds display data
        GLuint      textureIDk_RGB;         // OpenGL needs this
        GLuint      textureIDj_RGB;         // OpenGL needs this
        GLuint      textureIDi_RGB;         // OpenGL needs this
        // Correlation related ----------------------------------------------------------------------------------------
        uint64_t    voxel_k;               // A selected or hovered over voxel index
        uint64_t    voxel_j;               // A selected or hovered over voxel index
        uint64_t    voxel_i;               // A selected or hovered over voxel index
        uint64_t    voxel_t;               // A selected or hovered over voxel index
        uint64_t    voxel_index4D;         // A selected or hovered over voxel index
        int         time_course_onset;     // Omit volumes from start until this number
        int         time_course_offset;    // Omit volumes from end until this number
        int         time_course_shift;     // Shift time course data by this amount
        float*      p_sliceK_float_corr;   // Holds correlation data
        float*      p_sliceJ_float_corr;   // Holds correlation data
        float*      p_sliceI_float_corr;   // Holds correlation data
        int         visualization_mode;    // 0: grayscale, 1: RGB red overlay, 3: RGB correlation overlay
        // Add more file-related information as needed
    };

    struct FileList
    {
        std::vector<FileInfo> files;

        // ============================================================================================================
        // Procedure to add a new file to the list
        // ============================================================================================================
        void addFile(const std::string& name, const std::string& path)
        {
            // Read only the header information
            NiftiInfo nii_header;
            loadNifti(path, nii_header);
            printNiftiInfo(nii_header);  // Optional
            // Save relevant information here
            files.push_back({name, path, nii_header});
        }

        // ============================================================================================================
        // Procedure to remove an element from the list
        // ============================================================================================================
        void removeFile(size_t indexToRemove)
        {
            if (indexToRemove < files.size()) {
                files.erase(files.begin() + indexToRemove);
            }
        }

        // ============================================================================================================
        // Procedure to save format independent image information
        // ============================================================================================================
        void fillFileInfo(FileInfo& fi)
        {
            fi.dim_i = fi.header.dim[1];
            fi.dim_j = fi.header.dim[2];
            fi.dim_k = fi.header.dim[3];
            fi.dim_t = fi.header.dim[4];
            fi.pixdim_i = fi.header.pixdim[1];
            fi.pixdim_j = fi.header.pixdim[2];
            fi.pixdim_k = fi.header.pixdim[3];
            fi.nr_voxels = static_cast<uint64_t>(fi.header.dim[1]) * fi.header.dim[2] * fi.header.dim[3];
            fi.voxel_volume = fi.header.pixdim[1] * fi.header.pixdim[2] * fi.header.pixdim[3];
            fi.display_scale = 1.0;
            fi.visualization_mode = 0;

        }

        // ============================================================================================================
        // Function to access nifti files
        // ============================================================================================================
        int loadNifti(std::string filepath, NiftiInfo& header)
        {
            // Open the gzip-compressed file using zlib
            const char* cString = filepath.c_str();
            const gzFile file = gzopen(cString, "rb");

            // Check if the file is successfully opened
            if (file == nullptr) {
                printf("Error opening file: %s\n", cString);
                return 1;
            }

            // Read into the structure byte by byte
            gzread(file, &header.sizeof_hdr, 4);  // Must be 348

            // Unused fields
            gzread(file, &header.data_type, 10);     // unused
            gzread(file, &header.db_name, 18);       // unused
            gzread(file, &header.extents, 4);        // unused
            gzread(file, &header.session_error, 2);  // unused
            gzread(file, &header.regular, 1);        // unused

            // 1-byte
            gzread(file, &header.dim_info, 1);

            // 2-bytes
            gzread(file, &header.dim, 16);  // 2*8

            // Unused fields
            gzread(file, &header.intent_p1, 4);  // unused
            gzread(file, &header.intent_p2, 4);  // unused
            gzread(file, &header.intent_p3, 4);  // unused

            // 2-bytes
            gzread(file, &header.intent_code, 2);
            gzread(file, &header.datatype, 2);
            gzread(file, &header.bitpix, 2);
            gzread(file, &header.slice_start, 2);

            // 4-bytes
            gzread(file, &header.pixdim, 32);  // 4*8
            gzread(file, &header.vox_offset, 4);
            gzread(file, &header.scl_slope, 4);
            gzread(file, &header.scl_inter, 4);

            // 2-bytes
            gzread(file, &header.slice_end, 2);

            // 1-byte
            gzread(file, &header.slice_code, 1);
            gzread(file, &header.xyzt_units, 1);

            // 4-bytes
            gzread(file, &header.cal_max, 4);
            gzread(file, &header.cal_min, 4);
            gzread(file, &header.slice_duration, 4);
            gzread(file, &header.toffset, 4);
            gzread(file, &header.glmax, 4);
            gzread(file, &header.glmin, 4);

            // 1-byte
            gzread(file, &header.descrip, 80);  // 1*80
            gzread(file, &header.aux_file, 24);  // 1*24

            // 2-bytes
            gzread(file, &header.qform_code, 2);
            gzread(file, &header.sform_code, 2);

            // 4-bytes
            gzread(file, &header.quatern_b, 4);
            gzread(file, &header.quatern_c, 4);
            gzread(file, &header.quatern_d, 4);
            gzread(file, &header.qoffset_x, 4);
            gzread(file, &header.qoffset_y, 4);
            gzread(file, &header.qoffset_z, 4);
            gzread(file, &header.srow_x, 16);  // 4*4
            gzread(file, &header.srow_y, 16);  // 4*4
            gzread(file, &header.srow_z, 16);  // 4*4

            // 1-byte
            gzread(file, &header.intent_name, 16);  // 1*16
            gzread(file, &header.magic, 4);  // 1*4


            // NOTE[taken from nifti-1 documentation]: After the end of the 348 byte header (after the magic field),
            // the next 4 bytes are a char array field named "extension". By default, all 4 bytes of this array should 
            // be set to zero. In a .nii file, these 4 bytes will always be present, since the earliest start point for 
            // the image data is byte #352.
            gzread(file, &header.extension, 4);  // 1*4

            // Close the compressed file
            gzclose(file);

            return 0;
        }

        // ============================================================================================================
        // Procedure to print nifti header in the terminal
        // ============================================================================================================
        void printNiftiInfo(const NiftiInfo& header)
        {
            printf("sizeof_hdr    : %d\n", header.sizeof_hdr);

            // Unused fields
            // NOTE: These are left out

            printf("dim_info      : %u\n", header.dim_info);
            printf("dim           : [ %hd %hd %hd %hd %hd %hd %hd %hd ]\n", 
                header.dim[0], header.dim[1], header.dim[2], header.dim[3], 
                header.dim[4], header.dim[5], header.dim[6], header.dim[7]);

            // Unused fields
            // NOTE: These are left out

            printf("intent_code   : %hd\n", header.intent_code);
            printf("datatype      : %hd\n", header.datatype);
            printf("bitpix        : %hd\n", header.bitpix);
            printf("slice_start   : %hd\n", header.slice_start);
            printf("pixdim        : [ %f %f %f %f %f %f %f %f ]\n", 
                header.pixdim[0], header.pixdim[1], header.pixdim[2], header.pixdim[3], 
                header.pixdim[4], header.pixdim[5], header.pixdim[6], header.pixdim[7]);
            printf("vox_offset    : %f\n", header.vox_offset);
            printf("scl_slope     : %f\n", header.scl_slope);
            printf("scl_inter     : %f\n", header.scl_inter);
            printf("slice_end     : %hd\n", header.slice_end);
            printf("slice_code    : %u\n", header.slice_code);
            printf("xyzt_units    : %u\n", header.xyzt_units);
            printf("cal_max       : %f\n", header.cal_max);
            printf("cal_min       : %f\n", header.cal_min);
            printf("slice_duration: %f\n", header.slice_duration);
            printf("toffset       : %f\n", header.toffset);
            printf("glmax         : %d\n", header.glmax);
            printf("glmin         : %d\n", header.glmin);
            printf("descrip       : \"%s\"\n", header.descrip);
            printf("aux_file      : \"%s\"\n", header.aux_file);
            printf("qform_code    : %hd\n", header.qform_code);
            printf("sform_code    : %hd\n", header.sform_code);
            printf("quatern_b     : %f\n", header.quatern_b);
            printf("quatern_c     : %f\n", header.quatern_c);
            printf("quatern_d     : %f\n", header.quatern_d);
            printf("qoffset_x     : %f\n", header.qoffset_x);
            printf("qoffset_y     : %f\n", header.qoffset_y);
            printf("qoffset_z     : %f\n", header.qoffset_z);
            printf("srow_x        : [ %f %f %f %f ]\n",  // 1st row affine transform
                header.srow_x[0], header.srow_x[1], header.srow_x[2], header.srow_x[3]);
            printf("srow_y        : [ %f %f %f %f ]\n",  // 2nd row affine transform
                header.srow_y[0], header.srow_y[1], header.srow_y[2], header.srow_y[3]);
            printf("srow_z        : [ %f %f %f %f ]\n",  // 3rd row affine transform
                header.srow_z[0], header.srow_z[1], header.srow_z[2], header.srow_z[3]);
            printf("intent_name   : \"%s\"\n", header.intent_name);
            printf("magic         : %s\n", header.magic);  // Must be "ni1\0" or "n+1\0"
            printf("extension     : [ %d %d %d %d ]\n", 
                header.extension[0], header.extension[1], header.extension[2], header.extension[3]);
            printf("\n");
        }

        // ============================================================================================================
        // Procedure to load nifti data into memory
        // ============================================================================================================
        void loadNiftiDataTest(FileInfo& fi)
        {
            printf("\rLoading data...\n"); 

            // Free memory that could be filled if the function is run before
            // NOTE: I am not completely sure that whether this design migth cause crashes. It seems to stop memory
            // accumulation upon multiple "load data" clicks
            free(fi.p_data_float);
            free(fi.p_sliceK_float);
            free(fi.p_sliceJ_float);
            free(fi.p_sliceI_float);
            free(fi.p_sliceK_uint8);
            free(fi.p_sliceJ_uint8);
            free(fi.p_sliceI_uint8);
            free(fi.p_time_course_float);

            // Open the gzip-compressed file using zlib
            const char* cString = fi.path.c_str();
            const gzFile file = gzopen(cString, "rb");

            // Check if the file is successfully opened
            if (file == nullptr) {
                printf("Error opening file: %s\n", cString);
            }

            const uint64_t skipBytes = static_cast<uint64_t>(fi.header.vox_offset);
            gzseek(file, skipBytes, SEEK_SET);  // Skip header
            const uint64_t nr_data_points = fi.nr_voxels * fi.header.dim[4];
            
            // Allocate memory to float data
            const uint64_t array_size = nr_data_points * sizeof(float);
            fi.p_data_float = (float*)malloc(array_size);

            // Read buffer size
            const int READ_LIMIT = 1024 * 1024 * 1024;  // 1 GB in bytes

            // --------------------------------------------------------------------------------------------------------
            // Handle different data types. 
            // --------------------------------------------------------------------------------------------------------
            // | CODE | NAME       | NAME2              |     BIT |    BYTE |
            // |------|------------|--------------------|---------|---------|
            // |    2 | UINT8      | unsigned char      |   8 bit |  1 byte |
            // |  256 | INT8       | signed char        |   8 bit |  1 byte |
            // |    4 | INT16      | signed short       |  16 bit |  2 byte |
            // |  512 | UINT16     | unsigned short     |  16 bit |  2 byte |
            // |    8 | INT32      | signed int         |  32 bit |  4 byte |
            // |   16 | FLOAT32    | float              |  32 bit |  4 byte |
            // |  768 | UINT32     | unsigned int       |  32 bit |  4 byte |
            // |   64 | FLOAT64    | double             |  64 bit |  8 byte |
            // | 1024 | INT64      | long long          |  64 bit |  8 byte |
            // | 1280 | UINT64     | unsigned long long |  64 bit |  8 byte |
            // | 1536 | FLOAT128   | long double        | 128 bit | 16 byte |
            // |   32 | COMPLEX64  | complex            |  64 bit |  8 byte |  << Not implemented
            // | 1792 | COMPLEX128 | double pair        | 128 bit | 16 byte |  << Not implemented
            // | 2048 | COMPLEX256 | long double pair   | 256 bit | 32 byte |  << Not implemented
            // |  128 | RGB24      | RGB                |  24 bit |  3 byte |  << Not implemented
            // | 2304 | RGBA32     | RGBA               |  32 bit |  4 byte |  << Not implemented
            if (fi.header.datatype == 2) {
                uint64_t size_all_voxels = nr_data_points * sizeof(unsigned char);
                uint64_t nr_chunks = size_all_voxels / READ_LIMIT;  // Integer division, quotient
                uint64_t remainder = size_all_voxels % READ_LIMIT;  // Remainder
                uint64_t nr_voxels_per_chunk = READ_LIMIT / sizeof(unsigned char);
                uint64_t nr_voxels_remainder = remainder / sizeof(unsigned char);
                unsigned char* buffer = (unsigned char*)malloc(READ_LIMIT);
                for (uint64_t j = 0; j < nr_chunks; ++j) {
                    printf("\r  Reading chunks (%llu/%llu)...     ", j+1, nr_chunks+1); fflush(stdout);
                    gzread(file, buffer, READ_LIMIT);
                    for (uint64_t i = 0; i < nr_voxels_per_chunk; ++i) {
                        fi.p_data_float[j*nr_voxels_per_chunk + i] = (float)buffer[i];
                    }
                }
                printf("\r  Reading chunks (%llu/%llu)...     \n", nr_chunks+1, nr_chunks+1); 
                gzread(file, buffer, remainder);
                for (uint64_t i = 0; i < nr_voxels_remainder; ++i) {
                    fi.p_data_float[nr_chunks*nr_voxels_per_chunk + i] = (float)buffer[i];
                }
                gzclose(file);
                free(buffer);
            } else if (fi.header.datatype == 256) {
                uint64_t size_all_voxels = nr_data_points * sizeof(signed char);
                uint64_t nr_chunks = size_all_voxels / READ_LIMIT;  // Integer division, quotient
                uint64_t remainder = size_all_voxels % READ_LIMIT;  // Remainder
                uint64_t nr_voxels_per_chunk = READ_LIMIT / sizeof(signed char);
                uint64_t nr_voxels_remainder = remainder / sizeof(signed char);
                signed char* buffer = (signed char*)malloc(READ_LIMIT);
                for (uint64_t j = 0; j < nr_chunks; ++j) {
                    printf("\r  Reading chunks (%llu/%llu)...     ", j+1, nr_chunks+1); fflush(stdout);
                    gzread(file, buffer, READ_LIMIT);
                    for (uint64_t i = 0; i < nr_voxels_per_chunk; ++i) {
                        fi.p_data_float[j*nr_voxels_per_chunk + i] = (float)buffer[i];
                    }
                }
                printf("\r  Reading chunks (%llu/%llu)...     \n", nr_chunks+1, nr_chunks+1); 
                gzread(file, buffer, remainder);
                for (uint64_t i = 0; i < nr_voxels_remainder; ++i) {
                    fi.p_data_float[nr_chunks*nr_voxels_per_chunk + i] = (float)buffer[i];
                }
                gzclose(file);
                free(buffer);
            } else if (fi.header.datatype == 4) {
                uint64_t size_all_voxels = nr_data_points * sizeof(signed short);
                uint64_t nr_chunks = size_all_voxels / READ_LIMIT;  // Integer division, quotient
                uint64_t remainder = size_all_voxels % READ_LIMIT;  // Remainder
                uint64_t nr_voxels_per_chunk = READ_LIMIT / sizeof(signed short);
                uint64_t nr_voxels_remainder = remainder / sizeof(signed short);
                signed short* buffer = (signed short*)malloc(READ_LIMIT);
                for (uint64_t j = 0; j < nr_chunks; ++j) {
                    printf("\r  Reading chunks (%llu/%llu)...     ", j+1, nr_chunks+1); fflush(stdout);
                    gzread(file, buffer, READ_LIMIT);
                    for (uint64_t i = 0; i < nr_voxels_per_chunk; ++i) {
                        fi.p_data_float[j*nr_voxels_per_chunk + i] = (float)buffer[i];
                    }
                }
                printf("\r  Reading chunks (%llu/%llu)...     \n", nr_chunks+1, nr_chunks+1); 
                gzread(file, buffer, remainder);
                for (uint64_t i = 0; i < nr_voxels_remainder; ++i) {
                    fi.p_data_float[nr_chunks*nr_voxels_per_chunk + i] = (float)buffer[i];
                }
                gzclose(file);
                free(buffer);
            } else if (fi.header.datatype == 512) {
                uint64_t size_all_voxels = nr_data_points * sizeof(unsigned short);
                uint64_t nr_chunks = size_all_voxels / READ_LIMIT;  // Integer division, quotient
                uint64_t remainder = size_all_voxels % READ_LIMIT;  // Remainder
                uint64_t nr_voxels_per_chunk = READ_LIMIT / sizeof(unsigned short);
                uint64_t nr_voxels_remainder = remainder / sizeof(unsigned short);
                unsigned short* buffer = (unsigned short*)malloc(READ_LIMIT);
                for (uint64_t j = 0; j < nr_chunks; ++j) {
                    printf("\r  Reading chunks (%llu/%llu)...     ", j+1, nr_chunks+1); fflush(stdout);
                    gzread(file, buffer, READ_LIMIT);
                    for (uint64_t i = 0; i < nr_voxels_per_chunk; ++i) {
                        fi.p_data_float[j*nr_voxels_per_chunk + i] = (float)buffer[i];
                    }
                }
                printf("\r  Reading chunks (%llu/%llu)...     \n", nr_chunks+1, nr_chunks+1); 
                gzread(file, buffer, remainder);
                for (uint64_t i = 0; i < nr_voxels_remainder; ++i) {
                    fi.p_data_float[nr_chunks*nr_voxels_per_chunk + i] = (float)buffer[i];
                }
                gzclose(file);
                free(buffer);
            } else if (fi.header.datatype == 8) {
                uint64_t size_all_voxels = nr_data_points * sizeof(signed int);
                uint64_t nr_chunks = size_all_voxels / READ_LIMIT;  // Integer division, quotient
                uint64_t remainder = size_all_voxels % READ_LIMIT;  // Remainder
                uint64_t nr_voxels_per_chunk = READ_LIMIT / sizeof(signed int);
                uint64_t nr_voxels_remainder = remainder / sizeof(signed int);
                signed int* buffer = (signed int*)malloc(READ_LIMIT);
                for (uint64_t j = 0; j < nr_chunks; ++j) {
                    printf("\r  Reading chunks (%llu/%llu)...     ", j+1, nr_chunks+1); fflush(stdout);
                    gzread(file, buffer, READ_LIMIT);
                    for (uint64_t i = 0; i < nr_voxels_per_chunk; ++i) {
                        fi.p_data_float[j*nr_voxels_per_chunk + i] = (float)buffer[i];
                    }
                }
                printf("\r  Reading chunks (%llu/%llu)...     \n", nr_chunks+1, nr_chunks+1); 
                gzread(file, buffer, remainder);
                for (uint64_t i = 0; i < nr_voxels_remainder; ++i) {
                    fi.p_data_float[nr_chunks*nr_voxels_per_chunk + i] = (float)buffer[i];
                }
                gzclose(file);
                free(buffer);
            } else if (fi.header.datatype == 16) {
                uint64_t size_all_voxels = nr_data_points * sizeof(float);
                uint64_t nr_chunks = size_all_voxels / READ_LIMIT;  // Integer division, quotient
                uint64_t remainder = size_all_voxels % READ_LIMIT;  // Remainder
                uint64_t nr_voxels_per_chunk = READ_LIMIT / sizeof(float);
                uint64_t nr_voxels_remainder = remainder / sizeof(float);
                float* buffer = (float*)malloc(READ_LIMIT);
                for (uint64_t j = 0; j < nr_chunks; ++j) {
                    printf("\r  Reading chunks (%llu/%llu)...     ", j+1, nr_chunks+1); fflush(stdout);
                    gzread(file, buffer, READ_LIMIT);
                    for (uint64_t i = 0; i < nr_voxels_per_chunk; ++i) {
                        fi.p_data_float[j*nr_voxels_per_chunk + i] = (float)buffer[i];
                    }
                }
                printf("\r  Reading chunks (%llu/%llu)...     \n", nr_chunks+1, nr_chunks+1); 
                gzread(file, buffer, remainder);
                for (uint64_t i = 0; i < nr_voxels_remainder; ++i) {
                    fi.p_data_float[nr_chunks*nr_voxels_per_chunk + i] = (float)buffer[i];
                }
                gzclose(file);
                free(buffer);
            } else if (fi.header.datatype == 768) {
                uint64_t size_all_voxels = nr_data_points * sizeof(unsigned int);
                uint64_t nr_chunks = size_all_voxels / READ_LIMIT;  // Integer division, quotient
                uint64_t remainder = size_all_voxels % READ_LIMIT;  // Remainder
                uint64_t nr_voxels_per_chunk = READ_LIMIT / sizeof(unsigned int);
                uint64_t nr_voxels_remainder = remainder / sizeof(unsigned int);
                unsigned int* buffer = (unsigned int*)malloc(READ_LIMIT);
                for (uint64_t j = 0; j < nr_chunks; ++j) {
                    printf("\r  Reading chunks (%llu/%llu)...     ", j+1, nr_chunks+1); fflush(stdout);
                    gzread(file, buffer, READ_LIMIT);
                    for (uint64_t i = 0; i < nr_voxels_per_chunk; ++i) {
                        fi.p_data_float[j*nr_voxels_per_chunk + i] = (float)buffer[i];
                    }
                }
                printf("\r  Reading chunks (%llu/%llu)...     \n", nr_chunks+1, nr_chunks+1); 
                gzread(file, buffer, remainder);
                for (uint64_t i = 0; i < nr_voxels_remainder; ++i) {
                    fi.p_data_float[nr_chunks*nr_voxels_per_chunk + i] = (float)buffer[i];
                }
                gzclose(file);
                free(buffer);
            } else if (fi.header.datatype == 64) {
                uint64_t size_all_voxels = nr_data_points * sizeof(double);
                uint64_t nr_chunks = size_all_voxels / READ_LIMIT;  // Integer division, quotient
                uint64_t remainder = size_all_voxels % READ_LIMIT;  // Remainder
                uint64_t nr_voxels_per_chunk = READ_LIMIT / sizeof(double);
                uint64_t nr_voxels_remainder = remainder / sizeof(double);
                double* buffer = (double*)malloc(READ_LIMIT);
                for (uint64_t j = 0; j < nr_chunks; ++j) {
                    printf("\r  Reading chunks (%llu/%llu)...     ", j+1, nr_chunks+1); fflush(stdout);
                    gzread(file, buffer, READ_LIMIT);
                    for (uint64_t i = 0; i < nr_voxels_per_chunk; ++i) {
                        fi.p_data_float[j*nr_voxels_per_chunk + i] = (float)buffer[i];
                    }
                }
                printf("\r  Reading chunks (%llu/%llu)...     \n", nr_chunks+1, nr_chunks+1); 
                gzread(file, buffer, remainder);
                for (uint64_t i = 0; i < nr_voxels_remainder; ++i) {
                    fi.p_data_float[nr_chunks*nr_voxels_per_chunk + i] = (float)buffer[i];
                }
                gzclose(file);
                free(buffer);
            } else if (fi.header.datatype == 1024) {
                uint64_t size_all_voxels = nr_data_points * sizeof(long long);
                uint64_t nr_chunks = size_all_voxels / READ_LIMIT;  // Integer division, quotient
                uint64_t remainder = size_all_voxels % READ_LIMIT;  // Remainder
                uint64_t nr_voxels_per_chunk = READ_LIMIT / sizeof(long long);
                uint64_t nr_voxels_remainder = remainder / sizeof(long long);
                long long* buffer = (long long*)malloc(READ_LIMIT);
                for (uint64_t j = 0; j < nr_chunks; ++j) {
                    printf("\r  Reading chunks (%llu/%llu)...     ", j+1, nr_chunks+1); fflush(stdout);
                    gzread(file, buffer, READ_LIMIT);
                    for (uint64_t i = 0; i < nr_voxels_per_chunk; ++i) {
                        fi.p_data_float[j*nr_voxels_per_chunk + i] = (float)buffer[i];
                    }
                }
                printf("\r  Reading chunks (%llu/%llu)...     \n", nr_chunks+1, nr_chunks+1); 
                gzread(file, buffer, remainder);
                for (uint64_t i = 0; i < nr_voxels_remainder; ++i) {
                    fi.p_data_float[nr_chunks*nr_voxels_per_chunk + i] = (float)buffer[i];
                }
                gzclose(file);
                free(buffer);
            } else if (fi.header.datatype == 1280) {
                uint64_t size_all_voxels = nr_data_points * sizeof(unsigned long long);
                uint64_t nr_chunks = size_all_voxels / READ_LIMIT;  // Integer division, quotient
                uint64_t remainder = size_all_voxels % READ_LIMIT;  // Remainder
                uint64_t nr_voxels_per_chunk = READ_LIMIT / sizeof(unsigned long long);
                uint64_t nr_voxels_remainder = remainder / sizeof(unsigned long long);
                unsigned long long* buffer = (unsigned long long*)malloc(READ_LIMIT);
                for (uint64_t j = 0; j < nr_chunks; ++j) {
                    printf("\r  Reading chunks (%llu/%llu)...     ", j+1, nr_chunks+1); fflush(stdout);
                    gzread(file, buffer, READ_LIMIT);
                    for (uint64_t i = 0; i < nr_voxels_per_chunk; ++i) {
                        fi.p_data_float[j*nr_voxels_per_chunk + i] = (float)buffer[i];
                    }
                }
                printf("\r  Reading chunks (%llu/%llu)...     \n", nr_chunks+1, nr_chunks+1); 
                gzread(file, buffer, remainder);
                for (uint64_t i = 0; i < nr_voxels_remainder; ++i) {
                    fi.p_data_float[nr_chunks*nr_voxels_per_chunk + i] = (float)buffer[i];
                }
                gzclose(file);
                free(buffer);
            } else if (fi.header.datatype == 1536) {
                uint64_t size_all_voxels = nr_data_points * sizeof(long double);
                uint64_t nr_chunks = size_all_voxels / READ_LIMIT;  // Integer division, quotient
                uint64_t remainder = size_all_voxels % READ_LIMIT;  // Remainder
                uint64_t nr_voxels_per_chunk = READ_LIMIT / sizeof(long double);
                uint64_t nr_voxels_remainder = remainder / sizeof(long double);
                long double* buffer = (long double*)malloc(READ_LIMIT);
                for (uint64_t j = 0; j < nr_chunks; ++j) {
                    printf("\r  Reading chunks (%llu/%llu)...     ", j+1, nr_chunks+1); fflush(stdout);
                    gzread(file, buffer, READ_LIMIT);
                    for (uint64_t i = 0; i < nr_voxels_per_chunk; ++i) {
                        fi.p_data_float[j*nr_voxels_per_chunk + i] = (float)buffer[i];
                    }
                }
                printf("\r  Reading chunks (%llu/%llu)...     \n", nr_chunks+1, nr_chunks+1); 
                gzread(file, buffer, remainder);
                for (uint64_t i = 0; i < nr_voxels_remainder; ++i) {
                    fi.p_data_float[nr_chunks*nr_voxels_per_chunk + i] = (float)buffer[i];
                }
                gzclose(file);
                free(buffer);
            }

            // --------------------------------------------------------------------------------------------------------
            // Save useful format independent data here
            // --------------------------------------------------------------------------------------------------------
            // NOTE: Not sure if this is the best design
            findMinMax(fi);
            fi.display_i = fi.header.dim[1] / 2;
            fi.display_j = fi.header.dim[2] / 2;
            fi.display_k = fi.header.dim[3] / 2;
            fi.display_t = 0;

            // Initialize the visualization data slices
            fi.p_sliceK_float = (float*)malloc(fi.dim_i*fi.dim_j * sizeof(float));
            fi.p_sliceJ_float = (float*)malloc(fi.dim_i*fi.dim_k * sizeof(float));
            fi.p_sliceI_float = (float*)malloc(fi.dim_j*fi.dim_k * sizeof(float));

            // Initialize the visualization display slices
            fi.p_sliceK_uint8 = (uint8_t*)malloc(fi.dim_i*fi.dim_j * sizeof(uint8_t));
            fi.p_sliceJ_uint8 = (uint8_t*)malloc(fi.dim_i*fi.dim_k * sizeof(uint8_t));
            fi.p_sliceI_uint8 = (uint8_t*)malloc(fi.dim_j*fi.dim_k * sizeof(uint8_t));

            // Initialize the voxel data for timecourse visualizations;
            fi.p_time_course_float = (float*)malloc(fi.dim_t * sizeof(float));

            // Initialize hovered over or selected voxel index
            fi.voxel_i = static_cast<uint64_t>(fi.display_i);
            fi.voxel_j = static_cast<uint64_t>(fi.display_j);
            fi.voxel_k = static_cast<uint64_t>(fi.display_k);
            fi.voxel_t = static_cast<uint64_t>(fi.display_t);

            // Initialize time course parameters
            fi.time_course_onset  = 0;
            fi.time_course_offset = fi.header.dim[4];
            fi.time_course_shift = 0;
        }

        // ============================================================================================================
        // Procedures to load slices into memory as float
        // ============================================================================================================
        void loadSliceK_float(FileInfo& fi)
        {
            uint64_t ni = static_cast<uint64_t>(fi.dim_i);
            uint64_t nj = static_cast<uint64_t>(fi.dim_j);
            uint64_t k = static_cast<uint64_t>(fi.display_k);
            uint64_t t = static_cast<uint64_t>(fi.display_t);
            for (uint64_t i = 0; i < ni; ++i) {
                for (uint64_t j = 0; j < nj; ++j) {
                    uint64_t index4D = i + j*ni + k*ni*nj + fi.nr_voxels*t;
                    uint64_t index2D = i + j*ni;
                    fi.p_sliceK_float[index2D] = fi.p_data_float[index4D];
                }
            }
        }

        void loadSliceJ_float(FileInfo& fi)
        {
            uint64_t ni = static_cast<uint64_t>(fi.dim_i);
            uint64_t nj = static_cast<uint64_t>(fi.dim_j);
            uint64_t nk = static_cast<uint64_t>(fi.dim_k);
            uint64_t j = static_cast<uint64_t>(fi.display_j);
            uint64_t t = static_cast<uint64_t>(fi.display_t);
            for (uint64_t i = 0; i < ni; i++) {
                for (uint64_t k = 0; k < nk; k++) {
                    uint64_t index4D = i + j*ni + k*ni*nj + fi.nr_voxels*t;
                    uint64_t index2D = i + k*ni;
                    fi.p_sliceJ_float[index2D] = fi.p_data_float[index4D];
                }
            }
        }

        void loadSliceI_float(FileInfo& fi)
        {
            uint64_t ni = static_cast<uint64_t>(fi.dim_i);
            uint64_t nj = static_cast<uint64_t>(fi.dim_j);
            uint64_t nk = static_cast<uint64_t>(fi.dim_k);
            uint64_t i = static_cast<uint64_t>(fi.display_i);
            uint64_t t = static_cast<uint64_t>(fi.display_t);
            for (uint64_t j = 0; j < nj; j++) {
                for (uint64_t k = 0; k < nk; k++) {
                    uint64_t index4D = i + j*ni + k*ni*nj + fi.nr_voxels*t;
                    uint64_t index2D = j + k*nj;
                    fi.p_sliceI_float[index2D] = fi.p_data_float[index4D];
                }
            }
        }

        // ============================================================================================================
        // Procedure to type cast float into uint8 while clipping the range of uint8 display
        // ============================================================================================================
        // NOTE: I have decided to keep the functions separate because this might be more useful for covering 2D cases.
        // Not completely sure but keeping it simple for now.
        void loadSliceK_uint8(FileInfo& fi)
        {
            int nr_pixels = fi.dim_i * fi.dim_j;
            float thr_min = fi.display_min;
            float thr_max = fi.display_max;
            for (int i = 0; i < nr_pixels; i++) {
                if (fi.p_sliceK_float[i] > thr_max) {
                    fi.p_sliceK_uint8[i] = 255;
                } else if (fi.p_sliceK_float[i] < thr_min) {
                    fi.p_sliceK_uint8[i] = 0;
                } else {
                    fi.p_sliceK_uint8[i] = static_cast<uint8_t>((fi.p_sliceK_float[i] - thr_min)
                                                                            / (thr_max - thr_min) * 255);
                }
            }
        }

        void loadSliceJ_uint8(FileInfo& fi)
        {
            int nr_pixels = fi.dim_i * fi.dim_k;
            float thr_min = fi.display_min;
            float thr_max = fi.display_max;
            for (int i = 0; i < nr_pixels; i++) {
                if (fi.p_sliceJ_float[i] > thr_max) {
                    fi.p_sliceJ_uint8[i] = 255;
                } else if (fi.p_sliceJ_float[i] < thr_min) {
                    fi.p_sliceJ_uint8[i] = 0;
                } else {
                    fi.p_sliceJ_uint8[i] = static_cast<uint8_t>((fi.p_sliceJ_float[i] - thr_min)
                                                                            / (thr_max - thr_min) * 255);
                }
            }
        }

        void loadSliceI_uint8(FileInfo& fi)
        {
            int nr_pixels = fi.dim_j * fi.dim_k;
            float thr_min = fi.display_min;
            float thr_max = fi.display_max;
            for (int i = 0; i < nr_pixels; i++) {
                if (fi.p_sliceI_float[i] > thr_max) {
                    fi.p_sliceI_uint8[i] = 255;
                } else if (fi.p_sliceI_float[i] < thr_min) {
                    fi.p_sliceI_uint8[i] = 0;
                } else {
                    fi.p_sliceI_uint8[i] = static_cast<uint8_t>((fi.p_sliceI_float[i] - thr_min) / (thr_max - thr_min) * 255);
                }
            }
        }

        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // TODO: Procedure to load RGB overlay data into memory
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        void prepareRBGSlices(FileInfo& fi)
        {
            free(fi.p_sliceK_RGB_uint8);
            free(fi.p_sliceJ_RGB_uint8);
            free(fi.p_sliceI_RGB_uint8);

            // Initialize the visualization display slices
            fi.p_sliceK_RGB_uint8 = (uint8_t*)malloc(fi.dim_i*fi.dim_j * sizeof(uint8_t) * 3);
            fi.p_sliceJ_RGB_uint8 = (uint8_t*)malloc(fi.dim_i*fi.dim_k * sizeof(uint8_t) * 3);
            fi.p_sliceI_RGB_uint8 = (uint8_t*)malloc(fi.dim_j*fi.dim_k * sizeof(uint8_t) * 3);
        }

        void loadSliceK_RGB_uint8(FileInfo& fi)
        {
            for (uint64_t i = 0; i < fi.dim_i*fi.dim_j; i++) {
                if ( fi.p_sliceK_float[i] >= fi.overlay_min && fi.p_sliceK_float[i] <= fi.overlay_max) {
                    // Bimodal Highlight Blend
                    if (fi.p_sliceK_uint8[i] <= 122) {
                        fi.p_sliceK_RGB_uint8[i*3]   = 255 - fi.p_sliceK_uint8[i];
                        fi.p_sliceK_RGB_uint8[i*3+1] = 0;
                        fi.p_sliceK_RGB_uint8[i*3+2] = 0;

                    } else {
                        fi.p_sliceK_RGB_uint8[i*3]   = fi.p_sliceK_uint8[i];
                        fi.p_sliceK_RGB_uint8[i*3+1] = 0;
                        fi.p_sliceK_RGB_uint8[i*3+2] = 0;
                    }
                } else {
                    fi.p_sliceK_RGB_uint8[i*3]   = fi.p_sliceK_uint8[i];
                    fi.p_sliceK_RGB_uint8[i*3+1] = fi.p_sliceK_uint8[i];
                    fi.p_sliceK_RGB_uint8[i*3+2] = fi.p_sliceK_uint8[i];
                }
            }
        }

        void loadSliceJ_RGB_uint8(FileInfo& fi)
        {
            for (uint64_t i = 0; i < fi.dim_i*fi.dim_k; i++) {
                if ( fi.p_sliceJ_float[i] >= fi.overlay_min && fi.p_sliceJ_float[i] <= fi.overlay_max) {
                    // Bimodal Highlight Blend
                    if (fi.p_sliceJ_uint8[i] <= 122) {
                        fi.p_sliceJ_RGB_uint8[i*3]   = 255 - fi.p_sliceJ_uint8[i];
                        fi.p_sliceJ_RGB_uint8[i*3+1] = 0;
                        fi.p_sliceJ_RGB_uint8[i*3+2] = 0;

                    } else {
                        fi.p_sliceJ_RGB_uint8[i*3]   = fi.p_sliceJ_uint8[i];
                        fi.p_sliceJ_RGB_uint8[i*3+1] = 0;
                        fi.p_sliceJ_RGB_uint8[i*3+2] = 0;
                    }
                } else {
                    fi.p_sliceJ_RGB_uint8[i*3]   = fi.p_sliceJ_uint8[i];
                    fi.p_sliceJ_RGB_uint8[i*3+1] = fi.p_sliceJ_uint8[i];
                    fi.p_sliceJ_RGB_uint8[i*3+2] = fi.p_sliceJ_uint8[i];
                }
            }
        }

        void loadSliceI_RGB_uint8(FileInfo& fi)
        {
            for (uint64_t i = 0; i < fi.dim_j*fi.dim_k; i++) {
                if ( fi.p_sliceI_float[i] >= fi.overlay_min && fi.p_sliceI_float[i] <= fi.overlay_max) {
                    // Bimodal Highlight Blend
                    if (fi.p_sliceI_uint8[i] <= 122) {
                        fi.p_sliceI_RGB_uint8[i*3]   = 255 - fi.p_sliceI_uint8[i];
                        fi.p_sliceI_RGB_uint8[i*3+1] = 0;
                        fi.p_sliceI_RGB_uint8[i*3+2] = 0;

                    } else {
                        fi.p_sliceI_RGB_uint8[i*3]   = fi.p_sliceI_uint8[i];
                        fi.p_sliceI_RGB_uint8[i*3+1] = 0;
                        fi.p_sliceI_RGB_uint8[i*3+2] = 0;
                    }
                } else {
                    fi.p_sliceI_RGB_uint8[i*3]   = fi.p_sliceI_uint8[i];
                    fi.p_sliceI_RGB_uint8[i*3+1] = fi.p_sliceI_uint8[i];
                    fi.p_sliceI_RGB_uint8[i*3+2] = fi.p_sliceI_uint8[i];
                }
            }
        }

        // ============================================================================================================
        // Procedure upload OPENGL data used to visualize the image slice
        // ============================================================================================================
        void uploadTextureDataToOpenGL(int& dim1, int& dim2, GLuint& textureID, uint8_t* p_slice_uint8)
        {
            // Generate an OpenGL texture and bind it
            glGenTextures(1, &textureID);
            glBindTexture(GL_TEXTURE_2D, textureID);
            
            // Ensure row alignment is byte-precise
            glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

            // Upload the texture data to the OpenGL texture
            glTexImage2D(
                GL_TEXTURE_2D,      // Texture target
                0,                  // Mipmap level
                GL_RED,             // Internal format
                dim1, dim2,         // Width and height of the image
                0,                  // Border (must be 0)
                GL_RED,             // Format of the pixel data (only red channel)
                GL_UNSIGNED_BYTE,   // Data type of pixel data
                p_slice_uint8       // Pointer to the data
            );

            // Set texture filtering and wrapping parameters
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

            // Use single channel data to make grayscale
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_R, GL_RED);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_G, GL_RED);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_B, GL_RED);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_A, GL_ONE);  // Alpha always 1
            
            // Unbind the texture (TODO: Why did I do this? Understand.)
            glBindTexture(GL_TEXTURE_2D, 0);
        }

        void uploadTextureDataToOpenGL_RGB(int& dim1, int& dim2, GLuint& textureID, uint8_t* p_slice_uint8)
        {
            // Generate an OpenGL texture and bind it
            glGenTextures(1, &textureID);
            glBindTexture(GL_TEXTURE_2D, textureID);

            // Ensure row alignment is byte-precise
            glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
            
            // Upload the texture data to the OpenGL texture
            glTexImage2D(
                GL_TEXTURE_2D,      // Texture target
                0,                  // Mipmap level
                GL_RGB,             // Internal format
                dim1, dim2,         // Width and height of the image
                0,                  // Border (must be 0)
                GL_RGB,             // Format of the pixel data (only red channel)
                GL_UNSIGNED_BYTE,   // Data type of pixel data
                p_slice_uint8       // Pointer to the data
            );

            // Set texture filtering and wrapping parameters
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

            // Use single channel data to make grayscale
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_SWIZZLE_A, GL_ONE);  // Alpha always 1
            
            // Unbind the texture (TODO: Why did I do this? Understand.)
            glBindTexture(GL_TEXTURE_2D, 0);
        }

        // ============================================================================================================
        // Procedure update OPENGL data with new slice data
        // ============================================================================================================
        void updateTextureDataInOpenGL(int& dim1, int& dim2, GLuint& textureID, uint8_t* p_slice_uint8)
        {
            glBindTexture(GL_TEXTURE_2D, textureID);
            glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, dim1, dim2, GL_RED, GL_UNSIGNED_BYTE, p_slice_uint8);
            glBindTexture(GL_TEXTURE_2D, 0);
        }

        void updateTextureDataInOpenGL_RGB(int& dim1, int& dim2, GLuint& textureID, uint8_t* p_slice_uint8)
        {
            glBindTexture(GL_TEXTURE_2D, textureID);
            glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, dim1, dim2, GL_RGB, GL_UNSIGNED_BYTE, p_slice_uint8);
            glBindTexture(GL_TEXTURE_2D, 0);
        }

        // ============================================================================================================
        // Procedure to find minimum and maximum of the image data
        // ============================================================================================================
        void findMinMax(FileInfo& fi)
        {
            float max_val = std::numeric_limits<float>::min();
            float min_val = std::numeric_limits<float>::max();
            for (uint64_t i = 0; i < fi.nr_voxels; ++i) {
                if (fi.p_data_float[i] < min_val) {
                    min_val = fi.p_data_float[i];
                }
                if (fi.p_data_float[i] > max_val) {
                    max_val = fi.p_data_float[i];
                }
            }
            fi.data_max = max_val;
            fi.data_min = min_val;
            fi.display_max = max_val;
            fi.display_min = min_val;
            std::cout << "  Min: " << min_val << " | Max: " << max_val << std::endl;
        }

        // ============================================================================================================
        // Procedure to compute one voxel's time course correlations to other visible voxels
        // ============================================================================================================
        void computeCorrelationsForSlices_float(FileInfo& fi)
        {
            uint64_t ni = static_cast<uint64_t>(fi.dim_i);
            uint64_t nj = static_cast<uint64_t>(fi.dim_j);
            uint64_t nk = static_cast<uint64_t>(fi.dim_k);
            uint64_t nt = static_cast<uint64_t>(fi.time_course_offset) - fi.time_course_onset;

            // Prepare 1D arrays buffers that will be used to compute correlations
            float* x_arr = (float*)malloc(fi.dim_t * sizeof(float));
            float* y_arr = (float*)malloc(fi.dim_t * sizeof(float));

            // Prepare x array (selected voxel's data)
            for (uint64_t t = 0; t < nt; ++t) {
                uint64_t index4D = fi.voxel_i + fi.voxel_j*ni + fi.voxel_k*ni*nj + fi.nr_voxels*t;
                uint64_t tt = (t + fi.time_course_shift) % nt;
                *(x_arr + tt) = fi.p_data_float[index4D];
            }

            // --------------------------------------------------------------------------------------------------------
            // Compute correlations for slice k
            // --------------------------------------------------------------------------------------------------------
            uint64_t k = static_cast<uint64_t>(fi.display_k);
            for (uint64_t i = 0; i < ni; ++i) {
                for (uint64_t j = 0; j < nj; ++j) {
                    // Prepare y array
                    for (uint64_t t = 0; t < nt; ++t) {
                        uint64_t index4D = i + j*ni + k*ni*nj + fi.nr_voxels*t;
                        *(y_arr + t)= fi.p_data_float[index4D];
                    }

                    // Compute sums for covariance and variances
                    double sum_X = 0, sum_Y = 0, sum_XY = 0;
                    double sum_X2 = 0, sum_Y2 = 0;
                    for (int t = 0; t < nt; t++) {
                        sum_X += x_arr[t];
                        sum_Y += y_arr[t];
                        sum_XY += x_arr[t] * y_arr[t];
                        sum_X2 += x_arr[t] * x_arr[t];
                        sum_Y2 += y_arr[t] * y_arr[t];
                    }

                    // Compute covariance and variances
                    double numerator = sum_XY - (sum_X * sum_Y / nt);
                    double denominator = sqrt((sum_X2 - (sum_X * sum_X / nt)) * (sum_Y2 - (sum_Y * sum_Y / nt)));

                    // Handle edge cases where denominator might be zero
                    float r;
                    if (denominator == 0) {
                        r = 0.0f; // No correlation if variance is zero
                    }

                    // Compute the Pearson correlation coefficient
                    uint64_t index2D = i + j*ni;
                    r = static_cast<float>(numerator / denominator);
                    fi.p_sliceK_float_corr[index2D] = r;
                }
            }

            // --------------------------------------------------------------------------------------------------------
            // Compute correlations for slice j
            // --------------------------------------------------------------------------------------------------------
            uint64_t j = static_cast<uint64_t>(fi.display_j);
            for (uint64_t i = 0; i < ni; ++i) {
                for (uint64_t k = 0; k < nk; ++k) {
                    // Prepare y array
                    for (uint64_t t = 0; t < nt; ++t) {
                        uint64_t index4D = i + j*ni + k*ni*nj + fi.nr_voxels*t;
                        *(y_arr + t)= fi.p_data_float[index4D];
                    }

                    // Compute sums for covariance and variances
                    double sum_X = 0, sum_Y = 0, sum_XY = 0;
                    double sum_X2 = 0, sum_Y2 = 0;
                    for (int t = 0; t < nt; t++) {
                        sum_X += x_arr[t];
                        sum_Y += y_arr[t];
                        sum_XY += x_arr[t] * y_arr[t];
                        sum_X2 += x_arr[t] * x_arr[t];
                        sum_Y2 += y_arr[t] * y_arr[t];
                    }

                    // Compute covariance and variances
                    double numerator = sum_XY - (sum_X * sum_Y / nt);
                    double denominator = sqrt((sum_X2 - (sum_X * sum_X / nt)) * (sum_Y2 - (sum_Y * sum_Y / nt)));

                    // Handle edge cases where denominator might be zero
                    float r;
                    if (denominator == 0) {
                        r = 0.0f; // No correlation if variance is zero
                    }

                    // Compute the Pearson correlation coefficient
                    uint64_t index2D = i + k*ni;
                    r = static_cast<float>(numerator / denominator);
                    fi.p_sliceJ_float_corr[index2D] = r;
                }
            }

            // --------------------------------------------------------------------------------------------------------
            // Compute correlations for slice i
            // --------------------------------------------------------------------------------------------------------
            uint64_t i = static_cast<uint64_t>(fi.display_i);
            for (uint64_t j = 0; j < nj; ++j) {
                for (uint64_t k = 0; k < nk; ++k) {
                    // Prepare y array
                    for (uint64_t t = 0; t < nt; ++t) {
                        uint64_t index4D = i + j*ni + k*ni*nj + fi.nr_voxels*t;
                        *(y_arr + t)= fi.p_data_float[index4D];
                    }

                    // Compute sums for covariance and variances
                    double sum_X = 0, sum_Y = 0, sum_XY = 0;
                    double sum_X2 = 0, sum_Y2 = 0;
                    for (int t = 0; t < nt; t++) {
                        sum_X += x_arr[t];
                        sum_Y += y_arr[t];
                        sum_XY += x_arr[t] * y_arr[t];
                        sum_X2 += x_arr[t] * x_arr[t];
                        sum_Y2 += y_arr[t] * y_arr[t];
                    }

                    // Compute covariance and variances
                    double numerator = sum_XY - (sum_X * sum_Y / nt);
                    double denominator = sqrt((sum_X2 - (sum_X * sum_X / nt)) * (sum_Y2 - (sum_Y * sum_Y / nt)));

                    // Handle edge cases where denominator might be zero
                    float r;
                    if (denominator == 0) {
                        r = 0.0f; // No correlation if variance is zero
                    }

                    // Compute the Pearson correlation coefficient
                    uint64_t index2D = j + k*nj;
                    r = static_cast<float>(numerator / denominator);
                    fi.p_sliceI_float_corr[index2D] = r;
                }
            }
        }

        // ============================================================================================================
        // Correlation maps specific RGB image creation
        // ============================================================================================================
        void loadSliceK_Correlations_uint8(FileInfo& fi)
        {
            for (int i = 0; i < fi.dim_i*fi.dim_j; ++i) {
                float v = static_cast<float>(fi.p_sliceK_uint8[i]);
                float r = fi.p_sliceK_float_corr[i];

                // Red when r is close to 1, and grayscale underlay otwerwise
                if (r > 0) {
                    fi.p_sliceK_RGB_uint8[i*3]   = static_cast<uint8_t>( v * (1-r) + 255 * r );
                    fi.p_sliceK_RGB_uint8[i*3+1] = static_cast<uint8_t>( v * (1-r) );
                    fi.p_sliceK_RGB_uint8[i*3+2] = static_cast<uint8_t>( v * (1-r) );
                // Blue when r is close to -1, and grayscale underlay otwerwise
                } else {
                    fi.p_sliceK_RGB_uint8[i*3]   = static_cast<uint8_t>( v * (1+r) );
                    fi.p_sliceK_RGB_uint8[i*3+1] = static_cast<uint8_t>( v * (1+r) );
                    fi.p_sliceK_RGB_uint8[i*3+2] = static_cast<uint8_t>( v * (1+r) - 255 * r );
                }
            }
        }

        void loadSliceJ_Correlations_uint8(FileInfo& fi)
        {
            for (int i = 0; i < fi.dim_i*fi.dim_k; ++i) {
                float v = static_cast<float>(fi.p_sliceJ_uint8[i]);
                float r = fi.p_sliceJ_float_corr[i];

                // Red when r is close to 1, and grayscale underlay otwerwise
                if (r > 0) {
                    fi.p_sliceJ_RGB_uint8[i*3]   = static_cast<uint8_t>( v * (1-r) + 255 * r );
                    fi.p_sliceJ_RGB_uint8[i*3+1] = static_cast<uint8_t>( v * (1-r) );
                    fi.p_sliceJ_RGB_uint8[i*3+2] = static_cast<uint8_t>( v * (1-r) );
                // Blue when r is close to -1, and grayscale underlay otwerwise
                } else {
                    fi.p_sliceJ_RGB_uint8[i*3]   = static_cast<uint8_t>( v * (1+r) );
                    fi.p_sliceJ_RGB_uint8[i*3+1] = static_cast<uint8_t>( v * (1+r) );
                    fi.p_sliceJ_RGB_uint8[i*3+2] = static_cast<uint8_t>( v * (1+r) - 255 * r );
                }
            }
        }

        void loadSliceI_Correlations_uint8(FileInfo& fi)
        {
            for (int i = 0; i < fi.dim_j*fi.dim_k; ++i) {
                float v = static_cast<float>(fi.p_sliceI_uint8[i]);
                float r = fi.p_sliceI_float_corr[i];

                // Red when r is close to 1, and grayscale underlay otwerwise
                if (r > 0) {
                    fi.p_sliceI_RGB_uint8[i*3]   = static_cast<uint8_t>( v * (1-r) + 255 * r );
                    fi.p_sliceI_RGB_uint8[i*3+1] = static_cast<uint8_t>( v * (1-r) );
                    fi.p_sliceI_RGB_uint8[i*3+2] = static_cast<uint8_t>( v * (1-r) );
                // Blue when r is close to -1, and grayscale underlay otwerwise
                } else {
                    fi.p_sliceI_RGB_uint8[i*3]   = static_cast<uint8_t>( v * (1+r) );
                    fi.p_sliceI_RGB_uint8[i*3+1] = static_cast<uint8_t>( v * (1+r) );
                    fi.p_sliceI_RGB_uint8[i*3+2] = static_cast<uint8_t>( v * (1+r) - 255 * r );
                }
            }
        }
    };
}