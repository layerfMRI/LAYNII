#pragma once

#include "imgui.h"

// ====================================================================================================================
// Procedure to draw slice images
// ====================================================================================================================
void RenderSlice(int& dim1_vol, int& dim2_vol, int& dim3_vol, float dim1_sli, float dim2_sli, 
                 float& display_scale, float& display_offset_x, float& display_offset_y,
                 int& disp_idx_1, int& disp_idx_2, int& disp_idx_3, int& disp_idx_4,
                 int& visualization_mode, bool& show_slice_crosshair, bool& show_mouse_crosshair,
                 bool& request_image_data_update,
                 bool& show_voxel_inspector, bool& show_voxel_value, bool& show_voxel_indices,
                 GLuint& textureID, GLuint& textureID_RGB) {

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
};

// ====================================================================================================================
// Procedure to draw voxel inspector
// ====================================================================================================================
void RenderVoxelInspector(IDA_IO::FileInfo& fi, int slice_window,
                          bool& show_voxel_value, bool& show_voxel_indices) {

    // Definitions
    ImGuiIO& io = ImGui::GetIO();
    ImVec2 cursor_screen_pos;
    cursor_screen_pos = ImGui::GetCursorScreenPos();  // Relative to top-left corner image
    float scl = fi.display_scale;

    int ni = fi.dim_i;
    int nj = fi.dim_j;
    int nk = fi.dim_k;

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

        // ------------------------------------------------------------------------------------------------------------
        // Render inspector fields
        // ------------------------------------------------------------------------------------------------------------
        if ( show_voxel_value ) {
            uint64_t ni = static_cast<uint64_t>(fi.dim_i);
            uint64_t nj = static_cast<uint64_t>(fi.dim_j);
            uint64_t i = static_cast<uint64_t>(fi.display_i);
            uint64_t j = static_cast<uint64_t>(fi.display_j);
            uint64_t k = static_cast<uint64_t>(fi.display_k);
            uint64_t t = static_cast<uint64_t>(fi.display_t);
            uint64_t index4D = i + j*ni + k*ni*nj + fi.nr_voxels*t;
            ImGui::Text("Value : %.6f", fi.p_data_float[index4D]);                        
        }
        if ( show_voxel_indices ) {
            ImGui::Text("Index : [%d i, %d j, %d k, %d t]", fi.voxel_i, fi.voxel_j, fi.voxel_k, fi.voxel_t);                        
        }
    ImGui::EndTooltip();
    }

}
