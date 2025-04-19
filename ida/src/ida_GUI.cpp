#include "ida_GUI.h"
#include "ida_GUI_slices.h"
#include "../dep/idalib.h"

namespace IDA
{
	void RenderUI(bool& show_demo_window, bool& show_file_window, IDA_IO::FileList& fl)
	{
        // ============================================================================================================
        // Variables
        // ============================================================================================================
        // static char str_input[4096] = "Enter nifti path";
        static char str_input[4096] = "/Users/faruk/Documents/test-LN3_IDA/lo_BOLD_intemp.nii.gz";
        // static char str_input[4096] = "/Users/faruk/Documents/test-LN3_IDA/test.nii.gz";
        // static char str_input[4096] = "/Users/faruk/data/data-alard/75um/sub-99_75um_crop.nii.gz";

        static bool loaded_data          = false;
        static bool show_header_info     = false;
        static bool show_mouse_crosshair = true;
        static bool show_slice_crosshair   = true;
        static bool show_voxel_inspector = true;
        static bool show_voxel_indices   = false;
        static bool show_voxel_value     = true;
        static bool show_voxel_time_course = false;

        static bool request_image_data_update = false;
        static bool request_image_update      = false;

        static int sf = -1;  // Selected file

        int nr_files = fl.files.size();
        float spacing = ImGui::GetStyle().ItemInnerSpacing.x;

	    // ============================================================================================================
		// Enable Docking
	    // ============================================================================================================
        ImGui::DockSpaceOverViewport(ImGui::GetMainViewport());

        // ============================================================================================================
        // Keyboard and Mouse Controls Menu
        // ============================================================================================================
        ImGui::Begin("Keyboard & Mouse Controls + Debug"); 
        if (loaded_data)
        {
            ImGui::Text("Image controls");
            ImGui::Text("  Move             : Focus + CRTL + Move");
            ImGui::Text("  Zoom             : Hover + CRTL + Wheel");
            ImGui::Text("  Slice Scroll     : Hover + Wheel");
            ImGui::Text("  Move focus voxel : W, A, S, D, Q, E, Z, X");

            ImGui::Text("Time course controls");
            ImGui::Text("  Update reference : Right click on slice");

            ImGui::Text("Correlations controls");
            ImGui::Text("  Freeze correlations : Left click on slice");
            ImGui::Text("  Lag +1 : .");
            ImGui::Text("  Lag -1 : ,");
        }
        ImGui::End();

        // ------------------------------------------------------------------------------------------------------------
        // Input/Output Menu
        // ------------------------------------------------------------------------------------------------------------
        ImGui::Begin("File Menu");
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

                fl.files[sf].tc_nr = 1;

                fl.files[sf].loaded_data = true;
                loaded_data = fl.files[sf].loaded_data;
            }
        }

        // Display list of selectable file names
        for (int n = 0; n < nr_files; n++)
        {
            char buf[4096];
            snprintf(buf, sizeof(buf), "%d: %s", n, fl.files[n].name.c_str());
            if (ImGui::Selectable(buf, sf == n)) {
                sf = n;
                loaded_data = fl.files[sf].loaded_data;
            }
        }

        if (sf >= 0) {
            ImGui::Text("");
            ImGui::Text("Selected File:");
            ImGui::Text("  Number of voxels    : %llu", fl.files[sf].nr_voxels);
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
        // --------------------------------------------------------------------------------------------------------
        // Keyboard controls
        // --------------------------------------------------------------------------------------------------------
        if (loaded_data) {
            // k axis nativation
            if ( ImGui::IsKeyPressed(ImGuiKey_E, true) ) {
                fl.files[sf].display_k = (fl.files[sf].display_k + 1) % fl.files[sf].dim_k;

                fl.files[sf].voxel_i = static_cast<uint64_t>(fl.files[sf].display_i);
                fl.files[sf].voxel_j = static_cast<uint64_t>(fl.files[sf].display_j);
                fl.files[sf].voxel_k = static_cast<uint64_t>(fl.files[sf].display_k);
                request_image_data_update = true;
            }
            if ( ImGui::IsKeyPressed(ImGuiKey_Q, true) ) {
                if (fl.files[sf].display_k > 0) {
                    fl.files[sf].display_k--;
                } else {
                    fl.files[sf].display_k = fl.files[sf].dim_k - 1;
                }

                fl.files[sf].voxel_i = static_cast<uint64_t>(fl.files[sf].display_i);
                fl.files[sf].voxel_j = static_cast<uint64_t>(fl.files[sf].display_j);
                fl.files[sf].voxel_k = static_cast<uint64_t>(fl.files[sf].display_k);
                request_image_data_update = true;
            }

            // j axis nativation
            if ( ImGui::IsKeyPressed(ImGuiKey_W, true) ) {
                fl.files[sf].display_j = (fl.files[sf].display_j + 1) % fl.files[sf].dim_j;

                fl.files[sf].voxel_i = static_cast<uint64_t>(fl.files[sf].display_i);
                fl.files[sf].voxel_j = static_cast<uint64_t>(fl.files[sf].display_j);
                fl.files[sf].voxel_k = static_cast<uint64_t>(fl.files[sf].display_k);
                request_image_data_update = true;
            }
            if ( ImGui::IsKeyPressed(ImGuiKey_S, true) ) {
                if (fl.files[sf].display_j > 0) {
                    fl.files[sf].display_j--;
                } else {
                    fl.files[sf].display_j = fl.files[sf].dim_j - 1;
                }

                fl.files[sf].voxel_i = static_cast<uint64_t>(fl.files[sf].display_i);
                fl.files[sf].voxel_j = static_cast<uint64_t>(fl.files[sf].display_j);
                fl.files[sf].voxel_k = static_cast<uint64_t>(fl.files[sf].display_k);
                request_image_data_update = true;
            }

            // i axis nativation
            if ( ImGui::IsKeyPressed(ImGuiKey_A, true) ) {
                fl.files[sf].display_i = (fl.files[sf].display_i + 1) % fl.files[sf].dim_i;

                fl.files[sf].voxel_i = static_cast<uint64_t>(fl.files[sf].display_i);
                fl.files[sf].voxel_j = static_cast<uint64_t>(fl.files[sf].display_j);
                fl.files[sf].voxel_k = static_cast<uint64_t>(fl.files[sf].display_k);
                request_image_data_update = true;
            }
            if ( ImGui::IsKeyPressed(ImGuiKey_D, true) ) {
                if (fl.files[sf].display_i > 0) {
                    fl.files[sf].display_i--;
                } else {
                    fl.files[sf].display_i = fl.files[sf].dim_i - 1;
                }

                fl.files[sf].voxel_i = static_cast<uint64_t>(fl.files[sf].display_i);
                fl.files[sf].voxel_j = static_cast<uint64_t>(fl.files[sf].display_j);
                fl.files[sf].voxel_k = static_cast<uint64_t>(fl.files[sf].display_k);
                request_image_data_update = true;
            }

            // 4th axis nativation
            if ( ImGui::IsKeyPressed(ImGuiKey_X, true) ) {
                fl.files[sf].display_t = (fl.files[sf].display_t + 1) % fl.files[sf].dim_t;
                fl.files[sf].voxel_t = static_cast<uint64_t>(fl.files[sf].display_t);
                request_image_data_update = true;
            }
            if ( ImGui::IsKeyPressed(ImGuiKey_Z, true) ) {
                if (fl.files[sf].display_t > 0) {
                    fl.files[sf].display_t--;
                } else {
                    fl.files[sf].display_t = fl.files[sf].dim_t - 1;
                }
                fl.files[sf].voxel_t = static_cast<uint64_t>(fl.files[sf].display_t);
                request_image_data_update = true;
            }

            // --------------------------------------------------------------------------------------------------------
            // Center on focused voxel
            // TODO: Center on slice windows as well
            if (ImGui::IsKeyPressed(ImGuiKey_C, false)) {
                fl.files[sf].display_i = fl.files[sf].voxel_i;
                fl.files[sf].display_j = fl.files[sf].voxel_j;
                fl.files[sf].display_k = fl.files[sf].voxel_k;
                fl.files[sf].display_t = fl.files[sf].voxel_t;
                request_image_data_update = true;
            }

            // --------------------------------------------------------------------------------------------------------
            // Shift/lag (forward-backwad) the focused voxel's timepoints
            if (ImGui::IsKeyPressed(ImGuiKey_Comma, true)) {
                if (fl.files[sf].tc_shift < -fl.files[sf].dim_t + 1) {
                    fl.files[sf].tc_shift = fl.files[sf].dim_t;
                } else {
                    fl.files[sf].tc_shift = fl.files[sf].tc_shift - 1;
                }
                request_image_data_update = true;
            }
            if (ImGui::IsKeyPressed(ImGuiKey_Period, true)) {
                if (fl.files[sf].tc_shift > fl.files[sf].dim_t - 1) {
                    fl.files[sf].tc_shift = -fl.files[sf].dim_t;
                } else {
                    fl.files[sf].tc_shift = fl.files[sf].tc_shift + 1; 
                }
                request_image_data_update = true;
            }
        }

        // ============================================================================================================
        // Image Views
        // ============================================================================================================
        ImGui::Begin("Slice Z axis", nullptr, ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_NoScrollWithMouse);

        if (loaded_data)  // Slice K window
        {
            RenderSlice(
                fl.files[sf].dim_i, fl.files[sf].dim_j, fl.files[sf].dim_k,
                fl.files[sf].pixdim_j, fl.files[sf].pixdim_j, 
                fl.files[sf].display_scale, fl.files[sf].display_k_offset_x, fl.files[sf].display_k_offset_y,
                fl.files[sf].display_i, fl.files[sf].display_j, fl.files[sf].display_k, fl.files[sf].display_t,
                fl.files[sf].visualization_mode, show_slice_crosshair, show_mouse_crosshair,
                request_image_data_update,
                show_voxel_inspector, show_voxel_value, show_voxel_indices, show_voxel_time_course,
                fl.files[sf].textureIDk, fl.files[sf].textureIDk_RGB, 3, fl.files[sf]
                );
        }
        ImGui::End();

        ImGui::Begin("Slice Y axis", nullptr, ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
        if (loaded_data)  // Slice J window
        {
            RenderSlice(
                fl.files[sf].dim_i, fl.files[sf].dim_k, fl.files[sf].dim_j,
                fl.files[sf].pixdim_i, fl.files[sf].pixdim_k, 
                fl.files[sf].display_scale, fl.files[sf].display_j_offset_x, fl.files[sf].display_j_offset_y,
                fl.files[sf].display_i, fl.files[sf].display_k, fl.files[sf].display_j, fl.files[sf].display_t,
                fl.files[sf].visualization_mode, show_slice_crosshair, show_mouse_crosshair,
                request_image_data_update,
                show_voxel_inspector, show_voxel_value, show_voxel_indices, show_voxel_time_course,
                fl.files[sf].textureIDj, fl.files[sf].textureIDj_RGB, 2, fl.files[sf]
                );

        }
        ImGui::End();

        ImGui::Begin("Slice X axis", nullptr, ImGuiWindowFlags_HorizontalScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
        if (loaded_data)  // Slice I window
        {
            RenderSlice(
                fl.files[sf].dim_j, fl.files[sf].dim_k, fl.files[sf].dim_i,
                fl.files[sf].pixdim_j, fl.files[sf].pixdim_k, 
                fl.files[sf].display_scale, fl.files[sf].display_i_offset_x, fl.files[sf].display_i_offset_y,
                fl.files[sf].display_j, fl.files[sf].display_k, fl.files[sf].display_i,  fl.files[sf].display_t,
                fl.files[sf].visualization_mode, show_slice_crosshair, show_mouse_crosshair,
                request_image_data_update,
                show_voxel_inspector, show_voxel_value, show_voxel_indices, show_voxel_time_course,
                fl.files[sf].textureIDi, fl.files[sf].textureIDi_RGB, 1, fl.files[sf]
                );
        }
        ImGui::End();

        // ============================================================================================================
        // Time course view
        // ============================================================================================================
        if (loaded_data && fl.files[sf].dim_t > 1)
        {
            ImGui::Begin("Time Course View"); 

            ImVec2 plot_size = ImGui::GetContentRegionAvail();
            plot_size.y = 75.0f;
            ImVec2 pos = ImGui::GetCursorScreenPos();
            pos.y += 80.0f*0;
            ImDrawList* draw_list = ImGui::GetWindowDrawList();

            // Draw plot background
            draw_list->AddRectFilled(pos, ImVec2(pos.x + plot_size.x, pos.y + plot_size.y),
                                     IM_COL32(37, 50, 75, 255));
            draw_list->AddRect(pos, ImVec2(pos.x + plot_size.x, pos.y + plot_size.y),
                               IM_COL32(42, 55, 80, 255), 0.0f, 0, 2.0f);

            // Draw reference time course
            float max_y = fl.files[sf].tc_refer_max;
            float min_y = fl.files[sf].tc_refer_min;
            auto get_screen_point1 = [&](int i, float* data) {
                float x = pos.x + (i / float(fl.files[sf].dim_t - 1)) * plot_size.x;
                float y = pos.y + (1.0f - (data[i] - min_y) / (max_y - min_y)) * plot_size.y;
                return ImVec2(x, y);
            };

            for (int i = 0; i < fl.files[sf].dim_t - 1; i++) {
                draw_list->AddLine(get_screen_point1(i  , fl.files[sf].p_tc_refer_float),
                                   get_screen_point1(i+1, fl.files[sf].p_tc_refer_float),
                                   IM_COL32(255, 0, 0, 255), 1.0f);
            }

            // Draw focus time course
            max_y = fl.files[sf].tc_focus_max;
            min_y = fl.files[sf].tc_focus_min;
            auto get_screen_point2 = [&](int i, float* data) {
                float x = pos.x + (i / float(fl.files[sf].dim_t - 1)) * plot_size.x;
                float y = pos.y + (1.0f - (data[i] - min_y) / (max_y - min_y)) * plot_size.y;
                return ImVec2(x, y);
            };

            for (int i = 0; i < fl.files[sf].dim_t - 1; i++) {
                draw_list->AddLine(get_screen_point2(i  , fl.files[sf].p_tc_focus_float),
                                   get_screen_point2(i+1, fl.files[sf].p_tc_focus_float),
                                   IM_COL32(255, 255, 255, 255), 1.0f);
            }

            ImGui::End();
        }

        // ============================================================================================================
        // Navigation Control Menu
        // ============================================================================================================
        ImGui::Begin("Navigation");

        if (loaded_data)
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
            ImGui::SliderFloat("##Zoom", &fl.files[sf].display_scale, 0.1, 20, "%f");
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
                if (ImGui::Button("Enable Overlay")) {
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
                if (ImGui::Button("Disable Overlay")) {
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
            ImGui::SeparatorText("CROSSHAIR CONTROLS");
            // --------------------------------------------------------------------------------------------------------
            ImGui::Checkbox("Show mouse crosshair", &show_mouse_crosshair);
            ImGui::Checkbox("Show slice crosshair", &show_slice_crosshair);

            // --------------------------------------------------------------------------------------------------------
            ImGui::SeparatorText("VOXEL INSPECTOR CONTROLS");
            // --------------------------------------------------------------------------------------------------------
            ImGui::Checkbox("Toggle voxel inspector", &show_voxel_inspector);

            if ( !show_voxel_inspector ) {
                ImGui::BeginDisabled(true);
            }
            ImGui::Checkbox("Show value", &show_voxel_value); ImGui::SameLine();
            ImGui::Checkbox("Show indices", &show_voxel_indices); ImGui::SameLine();

            if ( fl.files[sf].dim_t > 1 ) {
                ImGui::Checkbox("Show time course", &show_voxel_time_course);
            }

            if (!show_voxel_inspector) {
                ImGui::EndDisabled();
            }


            if ( fl.files[sf].dim_t > 1 ) {
                // ----------------------------------------------------------------------------------------------------
                ImGui::SeparatorText("CORRELATIONS CONTROLS");
                // ----------------------------------------------------------------------------------------------------
                if ( fl.files[sf].visualization_mode != 3 ) {
                    if ( ImGui::Button("Enable Correlations") ) {
                        show_voxel_inspector = true;
                        fl.files[sf].voxel_i = static_cast<uint64_t>(fl.files[sf].display_i);
                        fl.files[sf].voxel_j = static_cast<uint64_t>(fl.files[sf].display_j);
                        fl.files[sf].voxel_k = static_cast<uint64_t>(fl.files[sf].display_k);
                        fl.files[sf].voxel_t = static_cast<uint64_t>(fl.files[sf].display_t);

                        fl.files[sf].tc_lock = false; 
                        fl.files[sf].tc_focus_voxel_i = fl.files[sf].voxel_i;
                        fl.files[sf].tc_focus_voxel_j = fl.files[sf].voxel_j;
                        fl.files[sf].tc_focus_voxel_k = fl.files[sf].voxel_k;
                        SampleVoxelTimeCourseFocus(fl.files[sf]);

                        fl.prepareRBGSlices(fl.files[sf]);

                        // Prepare data to hold correlation maps
                        free(fl.files[sf].p_sliceK_float_corr);
                        free(fl.files[sf].p_sliceJ_float_corr);
                        free(fl.files[sf].p_sliceI_float_corr);
                        fl.files[sf].p_sliceK_float_corr = (float*)malloc(fl.files[sf].dim_i*fl.files[sf].dim_j * sizeof(float));
                        fl.files[sf].p_sliceJ_float_corr = (float*)malloc(fl.files[sf].dim_i*fl.files[sf].dim_k * sizeof(float));
                        fl.files[sf].p_sliceI_float_corr = (float*)malloc(fl.files[sf].dim_j*fl.files[sf].dim_k * sizeof(float));

                        // Prepare corelations images
                        fl.loadSliceK_Correlations_uint8(fl.files[sf]);
                        fl.loadSliceJ_Correlations_uint8(fl.files[sf]);
                        fl.loadSliceI_Correlations_uint8(fl.files[sf]);

                        // Prepare corelations textures
                        fl.uploadTextureDataToOpenGL_RGB(fl.files[sf].dim_i, fl.files[sf].dim_j, 
                                                         fl.files[sf].textureIDk_RGB, fl.files[sf].p_sliceK_RGB_uint8);
                        fl.uploadTextureDataToOpenGL_RGB(fl.files[sf].dim_i, fl.files[sf].dim_k, 
                                                         fl.files[sf].textureIDj_RGB, fl.files[sf].p_sliceJ_RGB_uint8);
                        fl.uploadTextureDataToOpenGL_RGB(fl.files[sf].dim_j, fl.files[sf].dim_k, 
                                                         fl.files[sf].textureIDi_RGB, fl.files[sf].p_sliceI_RGB_uint8);

                        fl.files[sf].visualization_mode = 3;
                        request_image_data_update = true;
                    }
                } else {
                    if ( ImGui::Button("Disable Correlations") ) {
                        request_image_data_update = true;
                        fl.files[sf].visualization_mode = 0;
                    }

                    ImGui::SameLine();
                    if ( ImGui::Button("Save Map") ) {
                        fl.computeCorrelationsForVolume_float(fl.files[sf]);
                    }

                    // ------------------------------------------------------------------------------------------------
                    // Time course onset adjustments
                    // ------------------------------------------------------------------------------------------------
                    ImGui::PushButtonRepeat(true);
                    if (ImGui::ArrowButton("##TC_onset_decrease", ImGuiDir_Left)) {
                        fl.files[sf].tc_onset -= 1;
                        if ( fl.files[sf].tc_onset < 0 ) {
                            fl.files[sf].tc_onset = 0;
                        }
                        request_image_data_update = true;
                    }
                    ImGui::SameLine(0.0f, spacing);
                    if (ImGui::ArrowButton("##TC_onset_increase", ImGuiDir_Right)) {
                        fl.files[sf].tc_onset += 1;
                        if ( fl.files[sf].tc_onset == fl.files[sf].tc_offset ) {
                            fl.files[sf].tc_onset = fl.files[sf].tc_offset - 1;
                        }
                        request_image_data_update = true;
                    }
                    ImGui::PopButtonRepeat();
                    ImGui::SameLine();
                    if ( ImGui::SliderInt("Onset ", &fl.files[sf].tc_onset , 0, fl.files[sf].dim_t, "%i") ) {
                        if ( fl.files[sf].tc_onset >= fl.files[sf].tc_offset ) {
                            fl.files[sf].tc_onset = fl.files[sf].tc_offset - 1;
                        }
                        request_image_data_update = true;
                    }

                    // ------------------------------------------------------------------------------------------------
                    // Time course offset adjustments
                    // ------------------------------------------------------------------------------------------------
                    ImGui::PushButtonRepeat(true);
                    if (ImGui::ArrowButton("##TC_offset_decrease", ImGuiDir_Left)) {
                        fl.files[sf].tc_offset -= 1;
                        if ( fl.files[sf].tc_offset == fl.files[sf].tc_onset ) {
                            fl.files[sf].tc_offset = fl.files[sf].tc_onset + 1;
                        }
                        request_image_data_update = true;
                    }
                    ImGui::SameLine(0.0f, spacing);
                    if (ImGui::ArrowButton("##TC_offset_increase", ImGuiDir_Right)) {
                        fl.files[sf].tc_offset += 1;
                        if ( fl.files[sf].tc_offset > fl.files[sf].dim_t ) {
                            fl.files[sf].tc_offset = fl.files[sf].dim_t;
                        }
                        request_image_data_update = true;
                    }
                    ImGui::PopButtonRepeat();
                    ImGui::SameLine();
                    if ( ImGui::SliderInt("Offset", &fl.files[sf].tc_offset, 0, fl.files[sf].dim_t, "%i") ) {
                        if ( fl.files[sf].tc_offset <= fl.files[sf].tc_onset ) {
                            fl.files[sf].tc_offset = fl.files[sf].tc_onset + 1;
                        }
                        request_image_data_update = true;
                    }

                    // ------------------------------------------------------------------------------------------------
                    // Time course shift/lag adjustments
                    // ------------------------------------------------------------------------------------------------
                    ImGui::PushButtonRepeat(true);
                    if ( ImGui::ArrowButton("##TC_shift_decrease", ImGuiDir_Left) ) {
                        fl.files[sf].tc_shift -= 1;
                        if ( fl.files[sf].tc_shift < -fl.files[sf].dim_t ) {
                            fl.files[sf].tc_shift = -fl.files[sf].dim_t;
                        }
                        request_image_data_update = true;
                    }
                    ImGui::SameLine(0.0f, spacing);
                    if ( ImGui::ArrowButton("##TC_shift_increase", ImGuiDir_Right) ) {
                        fl.files[sf].tc_shift += 1;
                        if ( fl.files[sf].tc_shift > fl.files[sf].dim_t ) {
                            fl.files[sf].tc_shift = fl.files[sf].dim_t;
                        }
                        request_image_data_update = true;
                    }
                    ImGui::PopButtonRepeat();
                    ImGui::SameLine();
                    if ( ImGui::SliderInt("Shift", &fl.files[sf].tc_shift, -fl.files[sf].dim_t, fl.files[sf].dim_t, "%i") ) {
                        request_image_data_update = true;
                    }
                }
            }
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
        if (loaded_data) {

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
                    fl.loadSliceJ_float(fl.files[sf]);
                    fl.loadSliceI_float(fl.files[sf]);
                    request_image_update = true;
                }

                if ( request_image_update ) {
                    fl.computeCorrelationsForSlices_float(fl.files[sf]);

                    fl.loadSliceK_uint8(fl.files[sf]);
                    fl.loadSliceJ_uint8(fl.files[sf]);
                    fl.loadSliceI_uint8(fl.files[sf]);

                    fl.loadSliceK_Correlations_uint8(fl.files[sf]);
                    fl.loadSliceJ_Correlations_uint8(fl.files[sf]);
                    fl.loadSliceI_Correlations_uint8(fl.files[sf]);

                    fl.updateTextureDataInOpenGL_RGB(fl.files[sf].dim_i, fl.files[sf].dim_j, 
                                                     fl.files[sf].textureIDk_RGB, fl.files[sf].p_sliceK_RGB_uint8);                    
                    fl.updateTextureDataInOpenGL_RGB(fl.files[sf].dim_i, fl.files[sf].dim_k, 
                                                     fl.files[sf].textureIDj_RGB, fl.files[sf].p_sliceJ_RGB_uint8);                    
                    fl.updateTextureDataInOpenGL_RGB(fl.files[sf].dim_j, fl.files[sf].dim_k, 
                                                     fl.files[sf].textureIDi_RGB, fl.files[sf].p_sliceI_RGB_uint8);                    
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
            ImGui::Begin("Another Window", &show_file_window);  // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
            ImGui::Text("Hello from another window!");
            if (ImGui::Button("Close Me"))
                show_file_window = false;
            ImGui::End();
        }
	}
}