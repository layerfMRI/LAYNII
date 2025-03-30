#pragma once

#include "imgui.h"
#include "ida_IO.h"

namespace IDA
{
    void RenderUI(bool& show_demo_window, bool& show_file_window, IDA_IO::FileList& fileList);
}