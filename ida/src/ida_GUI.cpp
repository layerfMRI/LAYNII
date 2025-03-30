#include "ida_GUI.h"
#include "../dep/idalib.h"

namespace IDA
{
	void RenderUI(bool& show_demo_window, bool& show_file_window, IDA_IO::FileList& fl)
	{
	    // ====================================================================
		// Enable Docking
	    // ====================================================================
        ImGui::DockSpaceOverViewport(ImGui::GetMainViewport());

        // ------------------------------------------------------------------------------------------------------------
        // Input/Output Menu
        // ------------------------------------------------------------------------------------------------------------
        ImGui::Begin("File Menu");

        // Variables
        // static char str_input[4096] = "Enter nifti path";
        // static char str_input[4096] = "/Users/faruk/Git/LayNii/test_data/lo_BOLD_intemp.nii.gz";
        static char str_input[4096] = "/Users/faruk/Documents/test-LN3_IDA/sub-3003_ses-fine_task-fncloc_run-1_part-mag_bold.nii.gz";

        static bool loaded_file          = false;
        static bool show_header_info     = false;
        static bool show_slice_crosshair = false;
        static bool show_mouse_crosshair = true;
        static bool show_voxel_inspector = false;
        static bool show_voxel_indices   = true;
        static bool show_voxel_value     = true;
        static bool show_voxel_magnifier = false;
        static bool show_voxel_time_course  = false;

        static bool request_image_data_update = false;
        static bool request_image_update      = false;

        static int sf = -1;  // Selected file

        static int magnifier_window_size = 30;
        static int magnifier_zoom = 8;

        ImVec2 cursor_screen_pos;
        ImDrawList* drawList;

        ImU32 cross_mouse_color = IM_COL32(81, 113, 217, 255);  // RGBA
        ImU32 cross_slice_color = IM_COL32(222, 181, 61, 255);  // RGBA

        float cross_thickness = 1.0f;  // Pixels

        int nr_files = fl.files.size();
        float spacing = ImGui::GetStyle().ItemInnerSpacing.x;
        ImGuiIO& io = ImGui::GetIO();

        ImGui::Text("FPS %.0f ", ImGui::GetIO().Framerate);
        ImGui::SameLine();
        ImGui::Checkbox("Demo Window", &show_demo_window);  // NOTE: This checkbox is here for development.
        ImGui::SameLine();
        ImGui::Checkbox("Show Full Header", &show_header_info);

        // Enter file path as text
        ImGui::InputText("Input path", str_input, IM_ARRAYSIZE(str_input));

        // Add new file button
        if (ImGui::Button("Add"))
        {
            // Open a file dialog or get the file path in some way
            std::string filePath = str_input;
            // Extract the file name from the path
            std::string fileName = filePath.substr(filePath.find_last_of("/\\") + 1);
            fl.addFile(fileName, filePath);         // Add the new file to the list
            fl.fillFileInfo(fl.files[nr_files]);    // Fill minimum information
        }

        // Add remove file button
        ImGui::SameLine();
        if (ImGui::Button("Remove"))
        {
            fl.removeFile(sf);
        }

        // Add load data button
        ImGui::SameLine();
        if (ImGui::Button("Load Data"))
        {
            if (sf >= 0)
            {
                fl.loadNiftiDataTest(fl.files[sf]);

                fl.loadSliceK_float(fl.files[sf]);
                fl.loadSliceJ_float(fl.files[sf]);
                fl.loadSliceI_float(fl.files[sf]);

                fl.loadSliceK_uint8(fl.files[sf]);
                fl.loadSliceJ_uint8(fl.files[sf]);
                fl.loadSliceI_uint8(fl.files[sf]);

                fl.uploadTextureDataToOpenGL(fl.files[sf].dim_i, fl.files[sf].dim_j,
                                             fl.files[sf].textureIDk, fl.files[sf].p_sliceK_uint8);
                fl.uploadTextureDataToOpenGL(fl.files[sf].dim_i, fl.files[sf].dim_k,
                                             fl.files[sf].textureIDj, fl.files[sf].p_sliceJ_uint8);
                fl.uploadTextureDataToOpenGL(fl.files[sf].dim_j, fl.files[sf].dim_k,
                                             fl.files[sf].textureIDi, fl.files[sf].p_sliceI_uint8);

                fl.files[sf].overlay_min = fl.files[sf].display_min;
                fl.files[sf].overlay_max = fl.files[sf].display_max;

                loaded_file = true;
            }
        }

        // Display list of selectable file names
        for (int n = 0; n < nr_files; n++)
        {
            char buf[4096];
            snprintf(buf, sizeof(buf), "%d: %s", n, fl.files[n].name.c_str());
            if (ImGui::Selectable(buf, sf == n)) {
                sf = n;
            }
        }

        if (sf >= 0) {
            ImGui::Text("");
            ImGui::Text("Selected File:");
            ImGui::Text("  Number of voxels    : %d"  , fl.files[sf].nr_voxels);
            ImGui::Text("  Voxel volume        : %.3f", fl.files[sf].voxel_volume);
            ImGui::Text("  1st Data Axis       : %d"  , fl.files[sf].dim_i);
            ImGui::Text("  2nd Data Axis       : %d"  , fl.files[sf].dim_j);
            ImGui::Text("  3rd Data Axis       : %d"  , fl.files[sf].dim_k);
            ImGui::Text("  4th Data Axis       : %d"  , fl.files[sf].dim_t);
            ImGui::Text("  1st Voxel Dimension : %.3f", fl.files[sf].pixdim_i);
            ImGui::Text("  2nd Voxel Dimension : %.3f", fl.files[sf].pixdim_j);
            ImGui::Text("  3rd Voxel Dimension : %.3f", fl.files[sf].pixdim_k);
            ImGui::TextWrapped("Complete Path: \n%s" , fl.files[sf].path.c_str());
        }

        ImGui::End();

        // ============================================================================================================
        // Image Views
        // ============================================================================================================

        ImGui::Begin("Slice Z axis", nullptr, ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_NoScrollWithMouse);

        if (loaded_file)  // Slice K window
        {
            float data_w = fl.files[sf].dim_i;
            float data_h = fl.files[sf].dim_j;
            float pix_w = fl.files[sf].pixdim_i;
            float pix_h = fl.files[sf].pixdim_j;
            float pixscl_w = pix_w / pix_h;
            float scl = fl.files[sf].display_scale;
            float img_w = data_w * pixscl_w * scl;
            float img_h = data_h * scl;
            ImVec2 img_size = ImVec2(img_w, img_h);

            // Center the image center to window center
            float center_x =  - (img_w - ImGui::GetWindowSize()[0]) / 2;
            float center_y =  - (img_h - ImGui::GetWindowSize()[1]) / 2;
 
            // Move image (important to have it before image render for correct mouse inspector indexing)
            if ( ImGui::IsWindowHovered() || ImGui::IsWindowFocused() ) {
                if ( io.KeyCtrl ) {
                    fl.files[sf].display_k_offset_x += io.MouseDelta.x / scl;
                    fl.files[sf].display_k_offset_y += io.MouseDelta.y / scl;
                }
            }

            // Zoom image
            // TODO: Smoother zooming for the inspector can be achieved centering to the mouse position
            if ( ImGui::IsWindowHovered() && io.KeyCtrl && io.MouseWheel < 0 ) {
                fl.files[sf].display_scale -= 0.1;
            } else if ( ImGui::IsWindowHovered() && io.KeyCtrl && io.MouseWheel > 0 ) {
                fl.files[sf].display_scale += 0.1;
            }
 
            // Render image
            ImVec2 image_pos = ImVec2(center_x + fl.files[sf].display_k_offset_x * scl,
                                      center_y + fl.files[sf].display_k_offset_y * scl);
            ImGui::SetCursorPos(image_pos);
            cursor_screen_pos = ImGui::GetCursorScreenPos();  // Relative to top-left corner image
            ImVec2 uv_min = ImVec2(1.0f, 1.0f);  // Default (0.0f, 0.0f) is top-left
            ImVec2 uv_max = ImVec2(0.0f, 0.0f);  // Default (1.0f, 1.0f) is lower-right

            // Render texture
            if (fl.files[sf].visualization_mode == 0) {
                ImGui::Image((void*)(intptr_t)fl.files[sf].textureIDk       , img_size, uv_min, uv_max);
            } else if (fl.files[sf].visualization_mode == 1 || fl.files[sf].visualization_mode == 3) {
                ImGui::Image((void*)(intptr_t)fl.files[sf].textureIDk_RGB   , img_size, uv_min, uv_max);
            }

            // Draw crosshair at mouse position
            if ( show_slice_crosshair ) {
                drawList = ImGui::GetWindowDrawList();

                // Slice crosshair
                // TODO: Cross moves in the wrong direction due to OpenGL flip (lower left corner is 0 0)
                float offset_hori = cursor_screen_pos.y + fl.files[sf].display_j;
                // Horizontal
                drawList->AddLine(ImVec2(cursor_screen_pos.x        , offset_hori),
                                  ImVec2(cursor_screen_pos.x + img_w, offset_hori),
                                  cross_slice_color,
                                  cross_thickness);
                // Vertical
                drawList->AddLine(ImVec2(cursor_screen_pos.x + fl.files[sf].display_i, cursor_screen_pos.y        ),
                                  ImVec2(cursor_screen_pos.x + fl.files[sf].display_i, cursor_screen_pos.y + img_h),
                                  cross_slice_color,
                                  cross_thickness);
            }

            // Mouse crosshair
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

            // --------------------------------------------------------------------------------------------------------
            // Voxel inspector
            // --------------------------------------------------------------------------------------------------------
            if ( show_voxel_inspector ) {
                if ( ImGui::IsItemHovered(ImGuiHoveredFlags_DelayNone) && ImGui::BeginTooltip() ) {
                    int voxel_i = static_cast<int>((io.MousePos.x - cursor_screen_pos.x) / scl);
                    int voxel_j = static_cast<int>((io.MousePos.y - cursor_screen_pos.y) / scl);
                    int voxel_k = fl.files[sf].display_k;
                    int voxel_t = fl.files[sf].display_t;

                    // NOTE: Important, adjust for OpenGL flips
                    voxel_i = -voxel_i + fl.files[sf].dim_i-1;
                    voxel_j = -voxel_j + fl.files[sf].dim_j-1;

                    // Render inspector fields
                    if ( show_voxel_time_course ) {
                        fl.loadVoxelTimeCourse_float(fl.files[sf], voxel_i, voxel_j, voxel_k);
                        ImGui::PlotLines(
                            "Time Course",                                                     // Label
                            fl.files[sf].p_time_course_float,                                  // Values
                            fl.files[sf].time_course_offset - fl.files[sf].time_course_onset,  // Values count
                            0,                                                                 // Values offset
                            NULL,                                                              // Overlay Text
                            fl.files[sf].time_course_min,                                      // Scale min (FLT_MIN for auto)
                            fl.files[sf].time_course_max,                                      // Scale max (FLT_MAX for auto)
                            ImVec2(0, 100.0f)                                                  // Plot Size
                            );
                    }

                    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    // TODO: Incorporate this voxel index better
                    if ( fl.files[sf].visualization_mode == 3 ) {
                        if ( fl.files[sf].voxel_i != voxel_i || fl.files[sf].voxel_j != voxel_j || fl.files[sf].voxel_k != voxel_k ) {
                            fl.files[sf].voxel_i = voxel_i;
                            fl.files[sf].voxel_j = voxel_j;
                            fl.files[sf].voxel_k = voxel_k;
                            request_image_data_update = true;
                        }
                    }
                    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    if ( show_voxel_value ) {
                        float voxel_val = fl.files[sf].p_sliceK_float[voxel_j * fl.files[sf].dim_i + voxel_i];
                        ImGui::Text("Value : %.6f", voxel_val);                        
                    }
                    if ( show_voxel_indices ) {
                        ImGui::Text("Index : [%d i, %d j, %d k, %d t]", voxel_i, voxel_j, voxel_k, voxel_t);                        
                    }
                    if ( show_voxel_magnifier ) {  // Voxel window magnifier

                        // Step 1: Mouse pos relative to imgui image cursor
                        float subvoxel_i = io.MousePos.x - cursor_screen_pos.x;
                        float subvoxel_j = io.MousePos.y - cursor_screen_pos.y;

                        // Step 2: Adjust for general image scaling
                        subvoxel_i /= scl;
                        subvoxel_j /= scl;

                        // Step 3: Adjust for anisotropic voxels. NOTE: Generalize this for other windows, hmmm.
                        subvoxel_i /= pixscl_w;

                        // Step 4: Adjust for OpenGL flips. NOTE: This can break for unconventional data, hmmm.
                        subvoxel_i = -subvoxel_i + fl.files[sf].dim_i;
                        subvoxel_j = -subvoxel_j + fl.files[sf].dim_j;

                        // Determine magnifier window size
                        float winmag_sz = static_cast<float>(magnifier_window_size);
                        float winmag_x = subvoxel_i - winmag_sz * 0.5f;
                        float winmag_y = subvoxel_j - winmag_sz * 0.5f;
                        float winmag_zoom = static_cast<float>(magnifier_zoom);

                        // --------------------------------------------------------------------------------------------
                        ImVec2 winmag_pix = ImVec2( winmag_sz * winmag_zoom * pixscl_w,
                                                    winmag_sz * winmag_zoom );
                        ImVec2 uv0 = ImVec2( winmag_x / data_w,
                                             winmag_y / data_h );
                        ImVec2 uv1 = ImVec2( (winmag_x + winmag_sz) / data_w,
                                             (winmag_y + winmag_sz) / data_h );
                        ImGui::Image((void*)(intptr_t)fl.files[sf].textureIDk, winmag_pix, uv1, uv0);

                        ImGui::Text("Magnification");
                        ImGui::Text("  Factor  : %.0f [times]", winmag_zoom);
                        ImGui::Text("  Size    : %.0f X %.0f [pixels]", winmag_sz * winmag_zoom * pixscl_w, winmag_sz * winmag_zoom);
                        ImGui::Text("  Size    : %.0f X %.0f [~voxels]", winmag_sz, winmag_sz);
                        ImGui::Text("  Indices :\n    [%.0f:%.0f i, %.0f:%.0f j, %d k]", 
                            winmag_x, winmag_x + winmag_sz, winmag_y, winmag_y + winmag_sz, voxel_k);
                    }
                    ImGui::EndTooltip();
                }
            }
            // --------------------------------------------------------------------------------------------------------
            ImGui::Text("Pointer = %u", fl.files[sf].textureIDk);
            ImGui::SameLine();
            ImGui::Text("Size = %d x %d", fl.files[sf].dim_i, fl.files[sf].dim_j);
            ImGui::Text("Window size = [%.1f %.1f]", ImGui::GetWindowSize().x, ImGui::GetWindowSize().y);
            ImGui::Text("Window pos  = [%.1f %.1f]", ImGui::GetWindowPos().x, ImGui::GetWindowPos().y);

            // Scroll through slices (Hover + Wheel)
            if ( ImGui::IsWindowHovered() && io.MouseWheel < 0 && fl.files[sf].display_k > 0 && !io.KeyCtrl) {
                fl.files[sf].display_k--;
                request_image_data_update = true;

            } else if ( ImGui::IsWindowHovered() && io.MouseWheel > 0 && fl.files[sf].display_k < fl.files[sf].dim_k-1 && !io.KeyCtrl) {
                fl.files[sf].display_k++;
                request_image_data_update = true;
            }
        }
        ImGui::End();

        ImGui::Begin("Slice Y axis", nullptr, ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
        if (loaded_file)  // Slice J window
        {
            // Render texture
            if (fl.files[sf].visualization_mode == 0) {
                ImGui::Image((void*)(intptr_t)fl.files[sf].textureIDj, 
                             ImVec2(fl.files[sf].pixdim_i / fl.files[sf].pixdim_k * fl.files[sf].dim_i,
                                    fl.files[sf].dim_k),
                             ImVec2(1, 1),  // Horizontal
                             ImVec2(0, 0)   // Vertical
                             );
            } else if (fl.files[sf].visualization_mode == 1) {
                ImGui::Image((void*)(intptr_t)fl.files[sf].textureIDj_RGB, 
                             ImVec2(fl.files[sf].pixdim_i / fl.files[sf].pixdim_k * fl.files[sf].dim_i,
                                    fl.files[sf].dim_k),
                             ImVec2(1, 1),  // Horizontal
                             ImVec2(0, 0)   // Vertical
                             );       
            }

            // --------------------------------------------------------------------------------------------------------
            ImGui::Text("Pointer = %u", fl.files[sf].textureIDj);
            ImGui::SameLine();
            ImGui::Text("Size = %d x %d", fl.files[sf].dim_i, fl.files[sf].dim_k);
            ImGui::Text("Window size = [%.1f %.1f]", ImGui::GetWindowSize().x, ImGui::GetWindowSize().y);
            ImGui::Text("Window pos  = [%.1f %.1f]", ImGui::GetWindowPos().x, ImGui::GetWindowPos().y);

            // Mouse wheel browsing controls
            if ( ImGui::IsWindowHovered() && io.MouseWheel < 0 && fl.files[sf].display_j > 0) {
                fl.files[sf].display_j--;
                request_image_data_update = true;
            } else if ( ImGui::IsWindowHovered() && io.MouseWheel > 0 && fl.files[sf].display_j < fl.files[sf].dim_j-1 ) {
                fl.files[sf].display_j++;
                request_image_data_update = true;
            }
        }
        ImGui::End();

        ImGui::Begin("Slice X axis", nullptr, ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
        if (loaded_file)  // Slice I window
        {
            // Render texture
            if (fl.files[sf].visualization_mode == 0) {
                ImGui::Image((void*)(intptr_t)fl.files[sf].textureIDi, 
                             ImVec2(fl.files[sf].dim_j,
                                    fl.files[sf].dim_k),
                             ImVec2(1, 1),  // Horizontal
                             ImVec2(0, 0)   // Vertical
                             ); 
            } else if (fl.files[sf].visualization_mode == 1) {
                ImGui::Image((void*)(intptr_t)fl.files[sf].textureIDi_RGB, 
                             ImVec2(fl.files[sf].dim_j,
                                    fl.files[sf].dim_k),
                             ImVec2(1, 1),  // Horizontal
                             ImVec2(0, 0)   // Vertical
                             ); 
            }

            // --------------------------------------------------------------------------------------------------------
            ImGui::Text("Pointer = %u", fl.files[sf].textureIDi);
            ImGui::SameLine();
            ImGui::Text("Size = %d x %d", fl.files[sf].dim_j, fl.files[sf].dim_k);
            ImGui::Text("Window size = [%.1f %.1f]", ImGui::GetWindowSize().x, ImGui::GetWindowSize().y);
            ImGui::Text("Window pos  = [%.1f %.1f]", ImGui::GetWindowPos().x, ImGui::GetWindowPos().y);

            // Mouse wheel browsing controls
            if ( ImGui::IsWindowHovered() && io.MouseWheel < 0 && fl.files[sf].display_i > 0) {
                fl.files[sf].display_i--;
                request_image_data_update = true;
            } else if ( ImGui::IsWindowHovered() && io.MouseWheel > 0 && fl.files[sf].display_i < fl.files[sf].dim_i-1 ) {
                fl.files[sf].display_i++;
                request_image_data_update = true;
            }
        }
        ImGui::End();


        // ============================================================================================================
        // Navigation Control Menu
        // ============================================================================================================
        ImGui::Begin("Navigation");

        if (loaded_file)
        {
            // --------------------------------------------------------------------------------------------------------
            ImGui::SeparatorText("SLICE BROWSER CONTROLS");
            // --------------------------------------------------------------------------------------------------------
            ImGui::Text("Visualization Mode: %u", fl.files[sf].visualization_mode);

            // Slice k
            ImGui::PushButtonRepeat(true);
            if (ImGui::ArrowButton("##SliceK_decrease", ImGuiDir_Left)) {
                fl.files[sf].display_k = ((fl.files[sf].display_k - 1) % fl.files[sf].dim_k
                                          + fl.files[sf].dim_k) % fl.files[sf].dim_k;  // Always positive
                request_image_data_update = true;
            }
            ImGui::SameLine(0.0f, spacing);
            if (ImGui::ArrowButton("##SliceK_increase", ImGuiDir_Right)) {
                fl.files[sf].display_k = (fl.files[sf].display_k + 1) % fl.files[sf].dim_k;
                request_image_data_update = true;
            }
            ImGui::PopButtonRepeat();
            ImGui::SameLine();
            if (ImGui::SliderInt("##SliceK", &fl.files[sf].display_k, 0, fl.files[sf].dim_k-1, "%i")) {
                request_image_data_update = true;
            }     
            ImGui::SameLine();
            ImGui::Text("k [%u/%u]", fl.files[sf].display_k, fl.files[sf].dim_k-1);

            // Slice j
            ImGui::PushButtonRepeat(true);
            if (ImGui::ArrowButton("##SliceJ_decrease", ImGuiDir_Left)) {
                fl.files[sf].display_j = ((fl.files[sf].display_j - 1) % fl.files[sf].dim_j
                                          + fl.files[sf].dim_j) % fl.files[sf].dim_j;  // Always positive
                request_image_data_update = true;
            }
            ImGui::SameLine(0.0f, spacing);
            if (ImGui::ArrowButton("##SliceJ_increase", ImGuiDir_Right)) {
                fl.files[sf].display_j = (fl.files[sf].display_j + 1) % fl.files[sf].dim_j;
                request_image_data_update = true;
            }
            ImGui::PopButtonRepeat();
            ImGui::SameLine();
            if (ImGui::SliderInt("##SliceJ", &fl.files[sf].display_j, 0, fl.files[sf].dim_j-1, "%i")) {
                request_image_data_update = true;
            }     
            ImGui::SameLine();
            ImGui::Text("j [%u/%u]", fl.files[sf].display_j, fl.files[sf].dim_j-1);

            // Slice i
            ImGui::PushButtonRepeat(true);
            if (ImGui::ArrowButton("##SliceI_decrease", ImGuiDir_Left)) {
                fl.files[sf].display_i = ((fl.files[sf].display_i - 1) % fl.files[sf].dim_i
                                          + fl.files[sf].dim_i) % fl.files[sf].dim_i;  // Always positive
                request_image_data_update = true;
            }
            ImGui::SameLine(0.0f, spacing);
            if (ImGui::ArrowButton("##SliceI_increase", ImGuiDir_Right)) {
                fl.files[sf].display_i = (fl.files[sf].display_i + 1) % fl.files[sf].dim_i;
                request_image_data_update = true;
            }
            ImGui::PopButtonRepeat();
            ImGui::SameLine();
            if (ImGui::SliderInt("##SliceI", &fl.files[sf].display_i, 0, fl.files[sf].dim_i-1, "%i")) {
                request_image_data_update = true;
            }     
            ImGui::SameLine();
            ImGui::Text("i [%u/%u]", fl.files[sf].display_i, fl.files[sf].dim_i-1);

            // Slice t
            if (fl.files[sf].dim_t > 1) {
                ImGui::PushButtonRepeat(true);
                if (ImGui::ArrowButton("##SliceT_decrease", ImGuiDir_Left)) {
                    fl.files[sf].display_t = ((fl.files[sf].display_t - 1) % fl.files[sf].dim_t
                                              + fl.files[sf].dim_t) % fl.files[sf].dim_t;  // Always positive
                    request_image_data_update = true;
                }
                ImGui::SameLine(0.0f, spacing);
                if (ImGui::ArrowButton("##SliceT_increase", ImGuiDir_Right)) {
                    fl.files[sf].display_t = (fl.files[sf].display_t + 1) % fl.files[sf].dim_t;
                    request_image_data_update = true;
                }
                ImGui::PopButtonRepeat();
                ImGui::SameLine();
                if (ImGui::SliderInt("##SliceT", &fl.files[sf].display_t, 0, fl.files[sf].dim_t-1, "%i")) {
                    request_image_data_update = true;
                }     
                ImGui::SameLine();
                ImGui::Text("t [%u/%u]", fl.files[sf].display_t, fl.files[sf].dim_t-1);
            }


            // --------------------------------------------------------------------------------------------------------
            ImGui::SeparatorText("ZOOM CONTROLS");
            // --------------------------------------------------------------------------------------------------------

            ImGui::PushButtonRepeat(true);
            if (ImGui::ArrowButton("##Zoom_decrease", ImGuiDir_Left)) {
                fl.files[sf].display_scale -= 1;
            }
            ImGui::SameLine(0.0f, spacing);
            if (ImGui::ArrowButton("##Zoom_increase", ImGuiDir_Right)) {
                fl.files[sf].display_scale += 1;
            }
            ImGui::PopButtonRepeat();
            ImGui::SameLine();
            ImGui::SliderFloat("##Zoom", &fl.files[sf].display_scale, 0, 10, "%f");
            ImGui::SameLine();
            ImGui::Text("  %.3fX", fl.files[sf].display_scale);

            if (ImGui::Button("Reset")) {
                fl.files[sf].display_scale = 1.f;
                fl.files[sf].display_k_offset_x = 0.f;
                fl.files[sf].display_k_offset_y = 0.f;
                fl.files[sf].display_j_offset_x = 0.f;
                fl.files[sf].display_j_offset_y = 0.f;
                fl.files[sf].display_i_offset_x = 0.f;
                fl.files[sf].display_i_offset_y = 0.f;
            }

            // --------------------------------------------------------------------------------------------------------
            ImGui::SeparatorText("CONTRAST CONTROLS");
            // --------------------------------------------------------------------------------------------------------
            if ( ImGui::SliderFloat("Display Min.", &fl.files[sf].display_min, fl.files[sf].data_min, fl.files[sf].data_max, "%.2f") ) 
            {
                if (fl.files[sf].display_min > fl.files[sf].display_max) {
                    fl.files[sf].display_min = fl.files[sf].display_max;
                }
                request_image_update = true;
            }

            if ( ImGui::SliderFloat("Display Max.", &fl.files[sf].display_max, fl.files[sf].data_min, fl.files[sf].data_max, "%.2f") ) 
            {
                if (fl.files[sf].display_max < fl.files[sf].display_min) {
                    fl.files[sf].display_max = fl.files[sf].display_min;
                }
                request_image_update = true;
            }

            // ========================================================================================================
            ImGui::SeparatorText("MASK CONTROLS");
            // ========================================================================================================
            if ( fl.files[sf].visualization_mode != 1 ) {
                if (ImGui::Button("Enable Overlay") || ImGui::IsKeyPressed(ImGuiKey_S, false)) {
                    fl.prepareRBGSlices(fl.files[sf]);

                    fl.loadSliceK_uint8(fl.files[sf]);
                    fl.loadSliceJ_uint8(fl.files[sf]);
                    fl.loadSliceI_uint8(fl.files[sf]);

                    fl.loadSliceK_RGB_uint8(fl.files[sf]);
                    fl.loadSliceJ_RGB_uint8(fl.files[sf]);
                    fl.loadSliceI_RGB_uint8(fl.files[sf]);

                    fl.uploadTextureDataToOpenGL_RGB(fl.files[sf].dim_i, fl.files[sf].dim_j, 
                                                     fl.files[sf].textureIDk_RGB, fl.files[sf].p_sliceK_RGB_uint8);
                    fl.uploadTextureDataToOpenGL_RGB(fl.files[sf].dim_i, fl.files[sf].dim_k, 
                                                     fl.files[sf].textureIDj_RGB, fl.files[sf].p_sliceJ_RGB_uint8);
                    fl.uploadTextureDataToOpenGL_RGB(fl.files[sf].dim_j, fl.files[sf].dim_k, 
                                                     fl.files[sf].textureIDi_RGB, fl.files[sf].p_sliceI_RGB_uint8);

                    fl.files[sf].visualization_mode = 1;
                }
            } else if ( fl.files[sf].visualization_mode == 1 ) {
                if (ImGui::Button("Disable Overlay") || ImGui::IsKeyPressed(ImGuiKey_S, false)) {
                    request_image_data_update = true;
                    fl.files[sf].visualization_mode = 0;
                }
            }

            if (fl.files[sf].visualization_mode != 1) {
                ImGui::BeginDisabled(true);
            }

            if ( ImGui::SliderFloat("Mask Min.", &fl.files[sf].overlay_min, fl.files[sf].data_min, fl.files[sf].data_max, "%.2f") ) 
            {
                if (fl.files[sf].overlay_min > fl.files[sf].overlay_max) {
                    fl.files[sf].overlay_min = fl.files[sf].overlay_max;
                }
                request_image_update = true;
            }

            if ( ImGui::SliderFloat("Mask Max.", &fl.files[sf].overlay_max, fl.files[sf].data_min, fl.files[sf].data_max, "%.2f") ) 
            {
                if (fl.files[sf].overlay_max < fl.files[sf].overlay_min) {
                    fl.files[sf].overlay_max = fl.files[sf].overlay_min;
                }
                request_image_update = true;
            }

            if (fl.files[sf].visualization_mode != 1) {
                ImGui::EndDisabled();
            }

            // --------------------------------------------------------------------------------------------------------
            ImGui::SeparatorText("CORRELATIONS CONTROLS");
            // --------------------------------------------------------------------------------------------------------
            if (fl.files[sf].visualization_mode != 3) {
                if (ImGui::Button("Enable Correlations") || ImGui::IsKeyPressed(ImGuiKey_S, false)) {
                    fl.prepareRBGSlices(fl.files[sf]);

                    free(fl.files[sf].p_sliceK_float_corr);
                    fl.files[sf].p_sliceK_float_corr = (float*)malloc(fl.files[sf].dim_i*fl.files[sf].dim_j * sizeof(float));

                    fl.loadSliceK_Correlations_uint8(fl.files[sf]);

                    fl.uploadTextureDataToOpenGL_RGB(fl.files[sf].dim_i, fl.files[sf].dim_j, 
                                                     fl.files[sf].textureIDk_RGB, fl.files[sf].p_sliceK_RGB_uint8);

                    fl.files[sf].visualization_mode = 3;
                }
            } else {
                if (ImGui::Button("Disable Correlations") || ImGui::IsKeyPressed(ImGuiKey_S, false)) {
                    request_image_data_update = true;
                    fl.files[sf].visualization_mode = 0;
                }
            }

            if (fl.files[sf].visualization_mode == 3) {
                if ( ImGui::SliderInt("Time Course Onset ", &fl.files[sf].time_course_onset , 0, fl.files[sf].dim_t, "%i") ) {
                    request_image_data_update = true;
                }
                if ( ImGui::SliderInt("Time Course Offset", &fl.files[sf].time_course_offset, 0, fl.files[sf].dim_t, "%i") ) {
                    request_image_data_update = true;
                }
            }

            // --------------------------------------------------------------------------------------------------------
            ImGui::SeparatorText("VOXEL INSPECTOR CONTROLS");
            // --------------------------------------------------------------------------------------------------------
            ImGui::Checkbox("Toggle voxel inspector", &show_voxel_inspector);

            if (!show_voxel_inspector) {
                ImGui::BeginDisabled(true);
            }
            ImGui::Checkbox("Show value", &show_voxel_value); ImGui::SameLine();
            ImGui::Checkbox("Show indices", &show_voxel_indices); ImGui::SameLine();
            ImGui::Checkbox("Show time course", &show_voxel_time_course);

            ImGui::Checkbox("Show magnifier", &show_voxel_magnifier);
            if ( show_voxel_magnifier ) {
                ImGui::SliderInt("##Magnifier Voxel Size", &magnifier_window_size, 1, 100, "%i");
                ImGui::SameLine();
                ImGui::Text("%i X %i [voxels]", magnifier_window_size, magnifier_window_size);

                ImGui::SliderInt("##Magnifier Zoom", &magnifier_zoom, 1, 20, "%i");
                ImGui::SameLine();
                ImGui::Text("%i [times]", magnifier_zoom);
            }

            if (!show_voxel_inspector) {
                ImGui::EndDisabled();
            }

            // --------------------------------------------------------------------------------------------------------
            ImGui::SeparatorText("CROSSHAIR CONTROLS");
            // --------------------------------------------------------------------------------------------------------
            ImGui::Checkbox("Show mouse crosshair", &show_mouse_crosshair);
            // ImGui::Checkbox("Show slice crosshair", &show_slice_crosshair);
        }
        ImGui::End();


        // ============================================================================================================
        // Keyboard and Mouse Controls Menu
        // ============================================================================================================
        ImGui::Begin("Keyboard & Mouse Controls + Debug"); 
        if (loaded_file)
        {
            ImGui::Text("Image controls:");
            ImGui::Text("  Move         : Focus + CRTL + Move");
            ImGui::Text("  Zoom         : Hover + CRTL + Wheel");
            ImGui::Text("  Slice Scroll : Hover + Wheel");
            ImGui::Text("");
            ImGui::Text("Positioning parameters:");
            ImGui::Text("  MousePos                 : [%.3f x, %.3f y]", io.MousePos.x, io.MousePos.y);
            ImGui::Text("  CursorScreenPos          : [%.3f x, %.3f y]", cursor_screen_pos.x, cursor_screen_pos.y);
            // ImGui::Text("  CursorScreenPos Magnifier: [%.3f x, %.3f y]", cursor_pos_winmag.x, cursor_pos_winmag.y);
            ImGui::Text("");
            ImGui::Text("Move parameters:");
            ImGui::Text("  K offset: [%.3f x, %.3f y]", fl.files[sf].display_k_offset_x, fl.files[sf].display_k_offset_y);
            ImGui::Text("  J offset: [%.3f x, %.3f y]", fl.files[sf].display_j_offset_x, fl.files[sf].display_j_offset_y);
            ImGui::Text("  I offset: [%.3f x, %.3f y]", fl.files[sf].display_i_offset_x, fl.files[sf].display_i_offset_y);
            ImGui::Text("");
            ImGui::Text("Crosshair inspector: (WIP)");
            ImGui::Text("  Value : (WIP)");
            ImGui::Text("  Index : (WIP)");
        }
        ImGui::End();

        // ============================================================================================================
        // Display header information
        // ============================================================================================================
        if (show_header_info && sf >= 0) {
            ImGui::Begin("Header Info");
            auto& r_hdr = fl.files[sf].header;  // reference as shorthand

            ImGui::Text("\nNifti header information:");
            ImGui::Text("  sizeof_hdr    : %d\n", r_hdr.sizeof_hdr);
            ImGui::Text("  dim_info      : %u\n", r_hdr.dim_info);
            ImGui::Text("  dim           : [ %hd %hd %hd %hd %hd %hd %hd %hd ]", 
                r_hdr.dim[0], r_hdr.dim[1], r_hdr.dim[2], r_hdr.dim[3], 
                r_hdr.dim[4], r_hdr.dim[5], r_hdr.dim[6], r_hdr.dim[7]);
            ImGui::Text("  intent_code   : %hd\n", r_hdr.intent_code);
            ImGui::Text("  datatype      : %hd\n", r_hdr.datatype);
            ImGui::Text("  bitpix        : %hd\n", r_hdr.bitpix);
            ImGui::Text("  slice_start   : %hd\n", r_hdr.slice_start);
            ImGui::Text("  pixdim        : [ %f %f %f %f %f %f %f %f ]", 
                r_hdr.pixdim[0], r_hdr.pixdim[1], r_hdr.pixdim[2], r_hdr.pixdim[3], 
                r_hdr.pixdim[4], r_hdr.pixdim[5], r_hdr.pixdim[6], r_hdr.pixdim[7]);
            ImGui::Text("  vox_offset    : %f\n", r_hdr.vox_offset);
            ImGui::Text("  scl_slope     : %f\n", r_hdr.scl_slope);
            ImGui::Text("  scl_inter     : %f\n", r_hdr.scl_inter);
            ImGui::Text("  slice_end     : %hd\n", r_hdr.slice_end);
            ImGui::Text("  slice_code    : %u\n", r_hdr.slice_code);
            ImGui::Text("  xyzt_units    : %u\n", r_hdr.xyzt_units);
            ImGui::Text("  cal_max       : %f\n", r_hdr.cal_max);
            ImGui::Text("  cal_min       : %f\n", r_hdr.cal_min);
            ImGui::Text("  slice_duration: %f\n", r_hdr.slice_duration);
            ImGui::Text("  toffset       : %f\n", r_hdr.toffset);
            ImGui::Text("  glmax         : %d\n", r_hdr.glmax);
            ImGui::Text("  glmin         : %d\n", r_hdr.glmin);
            ImGui::Text("  descrip       : \"%s\"\n", r_hdr.descrip);
            ImGui::Text("  aux_file      : \"%s\"\n", r_hdr.aux_file);
            ImGui::Text("  qform_code    : %hd\n", r_hdr.qform_code);
            ImGui::Text("  sform_code    : %hd\n", r_hdr.sform_code);
            ImGui::Text("  quatern_b     : %f\n", r_hdr.quatern_b);
            ImGui::Text("  quatern_c     : %f\n", r_hdr.quatern_c);
            ImGui::Text("  quatern_d     : %f\n", r_hdr.quatern_d);
            ImGui::Text("  qoffset_x     : %f\n", r_hdr.qoffset_x);
            ImGui::Text("  qoffset_y     : %f\n", r_hdr.qoffset_y);
            ImGui::Text("  qoffset_z     : %f\n", r_hdr.qoffset_z);
            ImGui::Text("  srow_x        : [ %f %f %f %f ]", 
                r_hdr.srow_x[0], r_hdr.srow_x[1], r_hdr.srow_x[2], r_hdr.srow_x[3]);
            ImGui::Text("  srow_y        : [ %f %f %f %f ]", 
                r_hdr.srow_y[0], r_hdr.srow_y[1], r_hdr.srow_y[2], r_hdr.srow_y[3]);
            ImGui::Text("  srow_z        : [ %f %f %f %f ]", 
                r_hdr.srow_z[0], r_hdr.srow_z[1], r_hdr.srow_z[2], r_hdr.srow_z[3]);
            ImGui::Text("  intent_name   : \"%s\"\n", r_hdr.intent_name);
            ImGui::Text("  magic         : %s\n", r_hdr.magic);  // Must be "ni1\0" or "n+1\0"
            ImGui::Text("  extension     : [ %d %d %d %d ]", 
                r_hdr.extension[0], r_hdr.extension[1], r_hdr.extension[2], r_hdr.extension[3]);
            ImGui::End();
        }

        // ============================================================================================================
        // Update opengl textures based on the GUI interaction
        // ============================================================================================================
        if (loaded_file) {

            // --------------------------------------------------------------------------------------------------------
            // Grayscale mode
            // --------------------------------------------------------------------------------------------------------
            if ( fl.files[sf].visualization_mode == 0 ) {
                if ( request_image_data_update ) {
                    fl.loadSliceK_float(fl.files[sf]);
                    fl.loadSliceJ_float(fl.files[sf]);
                    fl.loadSliceI_float(fl.files[sf]);
                    request_image_update = true;
                }
                
                if ( request_image_update ) {
                    fl.loadSliceK_uint8(fl.files[sf]);
                    fl.updateTextureDataInOpenGL(fl.files[sf].dim_i, fl.files[sf].dim_j, 
                                                 fl.files[sf].textureIDk, fl.files[sf].p_sliceK_uint8);
                    fl.loadSliceJ_uint8(fl.files[sf]);
                    fl.updateTextureDataInOpenGL(fl.files[sf].dim_i, fl.files[sf].dim_k, 
                                                 fl.files[sf].textureIDj, fl.files[sf].p_sliceJ_uint8);
                    fl.loadSliceI_uint8(fl.files[sf]);
                    fl.updateTextureDataInOpenGL(fl.files[sf].dim_j, fl.files[sf].dim_k, 
                                                 fl.files[sf].textureIDi, fl.files[sf].p_sliceI_uint8);
                }

            // --------------------------------------------------------------------------------------------------------
            // Threshold mask mode
            // --------------------------------------------------------------------------------------------------------
            } else if ( fl.files[sf].visualization_mode == 1 ) {
                if ( request_image_data_update ) {
                    fl.loadSliceK_float(fl.files[sf]);
                    fl.loadSliceJ_float(fl.files[sf]);
                    fl.loadSliceI_float(fl.files[sf]);
                    request_image_update = true;
                }
                
                if ( request_image_update ) {
                    fl.loadSliceK_uint8(fl.files[sf]);
                    fl.loadSliceK_RGB_uint8(fl.files[sf]);
                    fl.updateTextureDataInOpenGL_RGB(fl.files[sf].dim_i, fl.files[sf].dim_j, 
                                                     fl.files[sf].textureIDk_RGB, fl.files[sf].p_sliceK_RGB_uint8);
                    fl.loadSliceJ_uint8(fl.files[sf]);
                    fl.loadSliceJ_RGB_uint8(fl.files[sf]);
                    fl.updateTextureDataInOpenGL_RGB(fl.files[sf].dim_i, fl.files[sf].dim_k, 
                                                     fl.files[sf].textureIDj_RGB, fl.files[sf].p_sliceJ_RGB_uint8);
                    fl.loadSliceI_uint8(fl.files[sf]);
                    fl.loadSliceI_RGB_uint8(fl.files[sf]);
                    fl.updateTextureDataInOpenGL_RGB(fl.files[sf].dim_j, fl.files[sf].dim_k, 
                                                     fl.files[sf].textureIDi_RGB, fl.files[sf].p_sliceI_RGB_uint8);
                }

            // --------------------------------------------------------------------------------------------------------
            // Voxel correlations mode
            // --------------------------------------------------------------------------------------------------------
            } else if ( fl.files[sf].visualization_mode == 3 ) {
                if ( request_image_data_update ) {
                    fl.loadSliceK_float(fl.files[sf]);
                    request_image_update = true;
                }

                if ( request_image_update) {
                    fl.computeCorrelationsSliceK_float(fl.files[sf]);
                    fl.loadSliceK_uint8(fl.files[sf]);
                    fl.loadSliceK_Correlations_uint8(fl.files[sf]);
                    fl.updateTextureDataInOpenGL_RGB(fl.files[sf].dim_i, fl.files[sf].dim_j, 
                                                     fl.files[sf].textureIDk_RGB, fl.files[sf].p_sliceK_RGB_uint8);                    
                }
            }

            request_image_data_update = false;
            request_image_update = false;
        }

        // ============================================================================================================
        // WIP
        // ============================================================================================================
        // Demo window
        if (show_demo_window)
            ImGui::ShowDemoWindow(&show_demo_window);

        // Show file selection window
        if (show_file_window)
        {
            ImGui::Begin("Another Window", &show_file_window);   // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
            ImGui::Text("Hello from another window!");
            if (ImGui::Button("Close Me"))
                show_file_window = false;
            ImGui::End();
        }
	}
}