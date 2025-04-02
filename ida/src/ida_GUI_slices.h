#pragma once

#include "imgui.h"

// ====================================================================================================================
// Procedure to draw voxel inspector
// ====================================================================================================================
void RenderVoxelInspector(IDA_IO::FileInfo& fi, int slice_window, ImVec2 cursor_screen_pos,
                          bool& show_voxel_value, bool& show_voxel_indices, bool& show_voxel_time_course) {

    // Definitions
    ImGuiIO& io = ImGui::GetIO();
    float scl = fi.display_scale;

    uint64_t ni = static_cast<uint64_t>(fi.dim_i);
    uint64_t nj = static_cast<uint64_t>(fi.dim_j);
    uint64_t nk = static_cast<uint64_t>(fi.dim_k);
    uint64_t nt = static_cast<uint64_t>(fi.dim_t);

    if ( ImGui::IsItemHovered(ImGuiHoveredFlags_DelayNone) && ImGui::BeginTooltip() ) {
        // ------------------------------------------------------------------------------------------------------------
        // Compute hovered over voxel index
        // ------------------------------------------------------------------------------------------------------------
        if ( slice_window == 3 ) {
            fi.voxel_i = -static_cast<int>((io.MousePos.x - cursor_screen_pos.x) / scl) + ni-1;
            fi.voxel_j = -static_cast<int>((io.MousePos.y - cursor_screen_pos.y) / scl) + nj-1;
            fi.voxel_k = fi.display_k;
            fi.voxel_t = fi.display_t;
        } else if ( slice_window == 2 ) {
            fi.voxel_i = -static_cast<int>((io.MousePos.x - cursor_screen_pos.x) / scl) + ni-1;
            fi.voxel_j = fi.display_j;
            fi.voxel_k = -static_cast<int>((io.MousePos.y - cursor_screen_pos.y) / scl) + nk-1;
            fi.voxel_t = fi.display_t;
        } else if ( slice_window == 1 ) {
            fi.voxel_i = fi.display_i;
            fi.voxel_j = -static_cast<int>((io.MousePos.x - cursor_screen_pos.x) / scl) + nj-1;
            fi.voxel_k = -static_cast<int>((io.MousePos.y - cursor_screen_pos.y) / scl) + nk-1;
            fi.voxel_t = fi.display_t;
        }

        // Compute data index
        uint64_t i = static_cast<uint64_t>(fi.voxel_i);
        uint64_t j = static_cast<uint64_t>(fi.voxel_j);
        uint64_t k = static_cast<uint64_t>(fi.voxel_k);
        uint64_t t = static_cast<uint64_t>(fi.voxel_t);
        uint64_t index4D = i + j*ni + k*ni*nj + fi.nr_voxels*t;

        // ------------------------------------------------------------------------------------------------------------
        // Pull voxel data from 4D into 1D memory, only if the hovered over voxel changes
        // ------------------------------------------------------------------------------------------------------------
        if ( fi.voxel_index4D != index4D ) {
            fi.voxel_index4D = index4D;

            // Load voxel data
            for (uint64_t t = 0; t < nt; ++t) {
                uint64_t index4D = i + j*ni + k*ni*nj + fi.nr_voxels*t;
                fi.p_time_course_float[t] = fi.p_data_float[index4D];
            }

            // Adjust min max for better visualizing the timecourse
            // NOTE: I might make a separate function for finding mix max in arbitrary data
            // NOTE: I can implement percent normalization etc here as well
            float max_val = std::numeric_limits<float>::min();
            float min_val = std::numeric_limits<float>::max();
            for (int t = 0; t < nt; ++t) {
                if (fi.p_time_course_float[t] < min_val) {
                    min_val = fi.p_time_course_float[t];
                }
                if (fi.p_time_course_float[t] > max_val) {
                    max_val = fi.p_time_course_float[t];
                }
            }
            fi.time_course_min = min_val;
            fi.time_course_max = max_val;
        }

        // ------------------------------------------------------------------------------------------------------------
        // Render inspector fields
        // ------------------------------------------------------------------------------------------------------------
        if ( show_voxel_value ) {
            ImGui::Text("Value : %.6f", fi.p_time_course_float[fi.voxel_t]);                        
        }
        if ( show_voxel_indices ) {
            ImGui::Text("Index : [%llu i, %llu j, %llu k, %llu t]", fi.voxel_i, fi.voxel_j, fi.voxel_k, fi.voxel_t);                        
        }

        if ( show_voxel_time_course ) {
            ImGui::Text("Time Course:");
            ImGui::PlotLines(
                "",                                            // Label
                fi.p_time_course_float,                        // Values
                fi.time_course_offset - fi.time_course_onset,  // Values count
                0,                                             // Values offset
                NULL,                                          // Overlay Text
                fi.time_course_min,                            // Scale min (FLT_MIN for auto)
                fi.time_course_max,                            // Scale max (FLT_MAX for auto)
                ImVec2(0, 100.0f)                              // Plot Size
                );
            ImGui::Text("[y: min-max scaled]");
        }
    ImGui::EndTooltip();
    }
}


// ====================================================================================================================
// Procedure to draw slice images
// ====================================================================================================================
void RenderSlice(int& dim1_vol, int& dim2_vol, int& dim3_vol, float dim1_sli, float dim2_sli, 
                 float& display_scale, float& display_offset_x, float& display_offset_y,
                 int& disp_idx_1, int& disp_idx_2, int& disp_idx_3, int& disp_idx_4,
                 int& visualization_mode, bool& show_slice_crosshair, bool& show_mouse_crosshair,
                 bool& request_image_data_update, 
                 bool& show_voxel_inspector, bool& show_voxel_value, bool& show_voxel_indices, bool& show_voxel_time_course,
                 GLuint& textureID, GLuint& textureID_RGB, int slice_window, IDA_IO::FileInfo& fi) {

    // Definitions
    ImGuiIO& io = ImGui::GetIO();
    ImVec2 cursor_screen_pos;
    ImDrawList* drawList;
    ImU32 cross_mouse_color = IM_COL32(81, 113, 217, 255);  // RGBA
    ImU32 cross_slice_color = IM_COL32(222, 181, 61, 255);  // RGBA
    float cross_thickness = 1.0f;  // Pixels

    // ----------------------------------------------------------------------------------------------------------------
    // Zoom image
    // TODO: Smoother zooming for the inspector can be achieved centering to the mouse position
    // ----------------------------------------------------------------------------------------------------------------
    if ( ImGui::IsWindowHovered() && io.KeyCtrl && io.MouseWheel < 0 ) {
        display_scale -= 0.1;
    } else if ( ImGui::IsWindowHovered() && io.KeyCtrl && io.MouseWheel > 0 ) {
        display_scale += 0.1;
    }

    // ----------------------------------------------------------------------------------------------------------------
    // Scroll through slices (Hover + Wheel)
    // ----------------------------------------------------------------------------------------------------------------
    if ( ImGui::IsWindowHovered() && io.MouseWheel < 0 && disp_idx_3 > 0 && !io.KeyCtrl) {
        disp_idx_3--;
        request_image_data_update = true;

    } else if ( ImGui::IsWindowHovered() && io.MouseWheel > 0 && disp_idx_3 < dim3_vol-1 && !io.KeyCtrl) {
        disp_idx_3++;
        request_image_data_update = true;
    }

    // ----------------------------------------------------------------------------------------------------------------
    // Prepare image position and scaling
    // ----------------------------------------------------------------------------------------------------------------
    // Compute pixel dimension scales
    float pixscl_w = dim1_sli / dim2_sli;
    float scl = display_scale;
    float img_w = static_cast<float>(dim1_vol) * pixscl_w * scl;
    float img_h = static_cast<float>(dim2_vol) * scl;
    ImVec2 img_size = ImVec2(img_w, img_h);

    // Center the image center to window center
    float center_x = -(img_w - ImGui::GetWindowSize()[0]) / 2;
    float center_y = -(img_h - ImGui::GetWindowSize()[1]) / 2;

    // Move image (important to have it before image render for correct mouse inspector indexing)
    if ( ImGui::IsWindowHovered() || ImGui::IsWindowFocused() ) {
        if ( io.KeyCtrl ) {
            display_offset_x += io.MouseDelta.x / scl;
            display_offset_y += io.MouseDelta.y / scl;
        }
    }

    // ----------------------------------------------------------------------------------------------------------------
    // Render image
    // ----------------------------------------------------------------------------------------------------------------
    ImVec2 image_pos = ImVec2(center_x + display_offset_x * scl,
                              center_y + display_offset_y * scl);
    ImGui::SetCursorPos(image_pos);
    cursor_screen_pos = ImGui::GetCursorScreenPos();  // Relative to top-left corner image
    ImVec2 uv_min = ImVec2(1.0f, 1.0f);  // Default (0.0f, 0.0f) is top-left
    ImVec2 uv_max = ImVec2(0.0f, 0.0f);  // Default (1.0f, 1.0f) is lower-right

    // Render texture
    if (visualization_mode == 0) {
        ImGui::Image((void*)(intptr_t)textureID, img_size, uv_min, uv_max);
    } else if (visualization_mode == 1 || visualization_mode == 3) {
        ImGui::Image((void*)(intptr_t)textureID_RGB, img_size, uv_min, uv_max);
    }

    // ----------------------------------------------------------------------------------------------------------------
    // Draw crosshair at mouse position
    // ----------------------------------------------------------------------------------------------------------------
    if ( show_slice_crosshair ) {
        drawList = ImGui::GetWindowDrawList();

        // Slice crosshair
        // TODO: Cross moves in the wrong direction due to OpenGL flip (lower left corner is 0 0)
        float offset_hori = cursor_screen_pos.y + disp_idx_2;
        // Horizontal
        drawList->AddLine(ImVec2(cursor_screen_pos.x        , offset_hori),
                          ImVec2(cursor_screen_pos.x + img_w, offset_hori),
                          cross_slice_color,
                          cross_thickness);
        // Vertical
        drawList->AddLine(ImVec2(cursor_screen_pos.x + disp_idx_1, cursor_screen_pos.y        ),
                          ImVec2(cursor_screen_pos.x + disp_idx_1, cursor_screen_pos.y + img_h),
                          cross_slice_color,
                          cross_thickness);
    }

    // ----------------------------------------------------------------------------------------------------------------
    // Mouse crosshair
    // ----------------------------------------------------------------------------------------------------------------
    if ( show_mouse_crosshair && ImGui::IsItemHovered() ) {
        drawList = ImGui::GetWindowDrawList();

        // Horizontal
        drawList->AddLine(ImVec2(cursor_screen_pos.x - 0.5f  , io.MousePos.y - 0.5f),
                          ImVec2(cursor_screen_pos.x - 0.5f + img_w, io.MousePos.y - 0.5f),
                          cross_mouse_color,
                          cross_thickness);
        // Vertical
        drawList->AddLine(ImVec2(io.MousePos.x - 0.5f, cursor_screen_pos.y - 0.5f),
                          ImVec2(io.MousePos.x - 0.5f, cursor_screen_pos.y - 0.5f + img_h),
                          cross_mouse_color,
                          cross_thickness);

        // Draw the circle using ImDrawList
        drawList->AddCircle(ImVec2(io.MousePos.x, io.MousePos.y),
                            5.0f, cross_mouse_color, 45, 1.0f);
    }

    if ( show_voxel_inspector ) {
        RenderVoxelInspector(fi, slice_window, cursor_screen_pos, show_voxel_value, show_voxel_indices, show_voxel_time_course);
    };

};
