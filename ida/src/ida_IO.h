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
        int         nr_voxels;          // Useful in general
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
        // Frangi related ---------------------------------------------------------------------------------------------
        float*      p_data_frangi;              // Copy original data here
        float*      p_data_scale1_eigvals_3D;   // Eigen values are consecutive per voxel (like in RGB)
        float*      p_data_scale2_eigvals_3D;
        float*      p_data_scale3_eigvals_3D;
        float*      p_sliceK_float_scale1_eigvals;  // Data slice that holds high precision data (e1 < e2 < e3)
        float*      p_sliceJ_float_scale1_eigvals;
        float*      p_sliceI_float_scale1_eigvals;
        float*      p_sliceK_float_scale2_eigvals;
        float*      p_sliceJ_float_scale2_eigvals;
        float*      p_sliceI_float_scale2_eigvals;
        float*      p_sliceK_float_scale3_eigvals;
        float*      p_sliceJ_float_scale3_eigvals;
        float*      p_sliceI_float_scale3_eigvals;
        float*      p_sliceK_float_scale1_vessellness;  // Data slice that holds a scalar for vessellness
        float*      p_sliceJ_float_scale1_vessellness;
        float*      p_sliceI_float_scale1_vessellness;
        float*      p_sliceK_float_scale2_vessellness;
        float*      p_sliceJ_float_scale2_vessellness;
        float*      p_sliceI_float_scale2_vessellness;
        float*      p_sliceK_float_scale3_vessellness;
        float*      p_sliceJ_float_scale3_vessellness;
        float*      p_sliceI_float_scale3_vessellness;
        float*      p_sliceK_float_final_vessellness;
        float*      p_sliceJ_float_final_vessellness;
        float*      p_sliceI_float_final_vessellness;
        uint8_t*    p_sliceK_RGB_uint8_frangi;  // RGB data slice that holds display data
        uint8_t*    p_sliceJ_RGB_uint8_frangi;  // RGB data slice that holds display data
        uint8_t*    p_sliceI_RGB_uint8_frangi;  // RGB data slice that holds display data
        GLuint      textureIDk_frangi;          // OpenGL needs this
        GLuint      textureIDj_frangi;          // OpenGL needs this
        GLuint      textureIDi_frangi;          // OpenGL needs this
        // Correlation related ----------------------------------------------------------------------------------------
        int         voxel_k;               // A selected or hovered over voxel index
        int         voxel_j;               // A selected or hovered over voxel index
        int         voxel_i;               // A selected or hovered over voxel index
        int         time_course_onset;     // Omit volumes from start until this number
        int         time_course_offset;    // Omit volumes from end until this number
        float*      p_sliceK_float_corr; // Holds correlation data
        int         visualization_mode;  // 0: grayscale, 1: red overlay, 2: frangi rgb
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
            fi.nr_voxels = fi.header.dim[1] * fi.header.dim[2] * fi.header.dim[3];
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

            const int skipBytes = static_cast<int>(fi.header.vox_offset);
            gzseek(file, skipBytes, SEEK_SET);  // Skip header
            const uint64_t nr_voxels = fi.header.dim[1] * fi.header.dim[2] * fi.header.dim[3];
            const uint64_t nr_timepoints = fi.header.dim[4];
            const uint64_t nr_data_points = nr_voxels * nr_timepoints;
            
            // Allocate memory to float data
            const size_t array_size = nr_data_points * sizeof(float);
            fi.p_data_float = (float*)malloc(array_size);

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
                const size_t        data_size = nr_data_points * sizeof(unsigned char);
                unsigned char*      data_orig = (unsigned char*)malloc(data_size);
                gzread(file, data_orig, data_size); gzclose(file);
                for (uint64_t i = 0; i < nr_data_points; ++i) fi.p_data_float[i] = (float)data_orig[i];
                free(data_orig);
            } else if (fi.header.datatype == 256) {
                const size_t        data_size = nr_data_points * sizeof(signed char);
                signed char*        data_orig = (signed char*)malloc(data_size);
                gzread(file, data_orig, data_size); gzclose(file);
                for (uint64_t i = 0; i < nr_data_points; ++i) fi.p_data_float[i] = (float)data_orig[i];
                free(data_orig);
            } else if (fi.header.datatype == 4) {
                const size_t        data_size = nr_data_points * sizeof(signed short);
                signed short*       data_orig = (signed short*)malloc(data_size);
                gzread(file, data_orig, data_size); gzclose(file);
                for (uint64_t i = 0; i < nr_data_points; ++i) fi.p_data_float[i] = (float)data_orig[i];
                free(data_orig);
            } else if (fi.header.datatype == 512) {
                const size_t        data_size = nr_data_points * sizeof(unsigned short);
                unsigned short*     data_orig = (unsigned short*)malloc(data_size);
                gzread(file, data_orig, data_size); gzclose(file);
                for (uint64_t i = 0; i < nr_data_points; ++i) fi.p_data_float[i] = (float)data_orig[i];
                free(data_orig);
            } else if (fi.header.datatype == 8) {
                const size_t        data_size = nr_data_points * sizeof(signed int);
                signed int*         data_orig = (signed int*)malloc(data_size);
                gzread(file, data_orig, data_size); gzclose(file);
                for (uint64_t i = 0; i < nr_data_points; ++i) fi.p_data_float[i] = (float)data_orig[i];
                free(data_orig);
            } else if (fi.header.datatype == 16) {  // Special case as the target data type is already float
                gzread(file, fi.p_data_float, array_size); gzclose(file);
            } else if (fi.header.datatype == 768) {
                const size_t        data_size = nr_data_points * sizeof(unsigned int);
                unsigned int*       data_orig = (unsigned int*)malloc(data_size);
                gzread(file, data_orig, data_size); gzclose(file);
                for (uint64_t i = 0; i < nr_data_points; ++i) fi.p_data_float[i] = (float)data_orig[i];
                free(data_orig);
            } else if (fi.header.datatype == 64) {
                const size_t        data_size = nr_data_points * sizeof(double);
                double*             data_orig = (double*)malloc(data_size);
                gzread(file, data_orig, data_size); gzclose(file);
                for (uint64_t i = 0; i < nr_data_points; ++i) fi.p_data_float[i] = (float)data_orig[i];
                free(data_orig);
            } else if (fi.header.datatype == 1024) {
                const size_t        data_size = nr_data_points * sizeof(long long);
                long long*          data_orig = (long long*)malloc(data_size);
                gzread(file, data_orig, data_size); gzclose(file);
                for (uint64_t i = 0; i < nr_data_points; ++i) fi.p_data_float[i] = (float)data_orig[i];
                free(data_orig);
            } else if (fi.header.datatype == 1280) {
                const size_t        data_size = nr_data_points * sizeof(unsigned long long);
                unsigned long long* data_orig = (unsigned long long*)malloc(data_size);
                gzread(file, data_orig, data_size); gzclose(file);
                for (uint64_t i = 0; i < nr_data_points; ++i) fi.p_data_float[i] = (float)data_orig[i];
                free(data_orig);
            } else if (fi.header.datatype == 1536) {
                const size_t        data_size = nr_data_points * sizeof(long double);
                long double*        data_orig = (long double*)malloc(data_size);
                gzread(file, data_orig, data_size); gzclose(file);
                for (uint64_t i = 0; i < nr_data_points; ++i) fi.p_data_float[i] = (float)data_orig[i];
                free(data_orig);
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
            fi.voxel_i = fi.display_i;
            fi.voxel_j = fi.display_j;
            fi.voxel_k = fi.display_k;

            // Initialize time course parameters
            fi.time_course_onset  = 0;
            fi.time_course_offset = fi.header.dim[4];
        }

        // ============================================================================================================
        // Procedures to load slices into memory as float
        // ============================================================================================================
        void loadSliceK_float(FileInfo& fi)
        {
            int ni = fi.dim_i;
            int nj = fi.dim_j;
            int k = fi.display_k;
            int t = fi.display_t;
            for (int i = 0; i < ni; i++) {
                for (int j = 0; j < nj; j++) {
                    int index4D = i + j*ni + k*ni*nj + fi.nr_voxels*t;
                    int index2D = i + j*ni;
                    fi.p_sliceK_float[index2D] = fi.p_data_float[index4D];
                }
            }
        }

        void loadSliceJ_float(FileInfo& fi)
        {
            int ni = fi.dim_i;
            int nj = fi.dim_j;
            int nk = fi.dim_k;
            int j = fi.display_j;
            int t = fi.display_t;
            for (int i = 0; i < ni; i++) {
                for (int k = 0; k < nk; k++) {
                    int index4D = i + j*ni + k*ni*nj + fi.nr_voxels*t;
                    int index2D = i + k*ni;
                    fi.p_sliceJ_float[index2D] = fi.p_data_float[index4D];
                }
            }
        }

        void loadSliceI_float(FileInfo& fi)
        {
            int ni = fi.dim_i;
            int nj = fi.dim_j;
            int nk = fi.dim_k;
            int i = fi.display_i;
            int t = fi.display_t;
            for (int j = 0; j < nj; j++) {
                for (int k = 0; k < nk; k++) {
                    int index4D = i + j*ni + k*ni*nj + fi.nr_voxels*t;
                    int index2D = j + k*nj;
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
            for (int i = 0; i < fi.dim_i*fi.dim_j; i++) {
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
            for (int i = 0; i < fi.dim_i*fi.dim_k; i++) {
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
            for (int i = 0; i < fi.dim_j*fi.dim_k; i++) {
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
            for (int i = 0; i < fi.nr_voxels; ++i) {
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
        // Procedure to extract one voxel's time course
        // ============================================================================================================
        void loadVoxelTimeCourse_float(FileInfo& fi, int i, int j, int k)
        {
            int ni = fi.dim_i;
            int nj = fi.dim_j;
            int nt = fi.time_course_offset - fi.time_course_onset;
            for (int t = 0; t < nt; ++t) {
                int index4D = i + j*ni + k*ni*nj + fi.nr_voxels*t;
                *(fi.p_time_course_float + t) = fi.p_data_float[index4D];
            }
            
            // --------------------------------------------------------------------------------------------------------
            // NOTE: I might make a separate function for finding mix max in arbitrary data
            // NOTE: I can implement percent normalization etc here as well
            float max_val = std::numeric_limits<float>::min();
            float min_val = std::numeric_limits<float>::max();
            for (int t = 0; t < nt; ++t) {
                if (*(fi.p_time_course_float + t) < min_val) {
                    min_val = *(fi.p_time_course_float + t);
                }
                if (*(fi.p_time_course_float + t) > max_val) {
                    max_val = *(fi.p_time_course_float + t);
                }
            }
            fi.time_course_min = min_val;
            fi.time_course_max = max_val;
            // --------------------------------------------------------------------------------------------------------
        }

        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // Procedure to compute one voxel's time course correlations to other voxels
        // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        void computeCorrelationsSliceK_float(FileInfo& fi)
        {
            int ni = fi.dim_i;
            int nj = fi.dim_j;
            int nt = fi.time_course_offset - fi.time_course_onset;

            float* x_arr = (float*)malloc(fi.dim_t * sizeof(float));
            float* y_arr = (float*)malloc(fi.dim_t * sizeof(float));

            // Prepare x array
            for (int t = 0; t < nt; ++t) {
                int index4D = fi.voxel_i + fi.voxel_j*ni + fi.voxel_k*ni*nj + fi.nr_voxels*t;
                *(x_arr + t) = fi.p_data_float[index4D];
            }

            int k = fi.display_k;
            for (int i = 0; i < ni; i++) {
                for (int j = 0; j < nj; j++) {

                    // Prepare y array
                    for (int t = 0; t < nt; ++t) {
                        int index4D = i + j*ni + k*ni*nj + fi.nr_voxels*t;
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
                    int index2D = i + j*ni;
                    r = static_cast<float>(numerator / denominator);
                    fi.p_sliceK_float_corr[index2D] = r;
                }
            }
        }

        void loadSliceK_Correlations_uint8(FileInfo& fi)
        {
            for (int i = 0; i < fi.dim_i*fi.dim_j; i++) {
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
    };
}