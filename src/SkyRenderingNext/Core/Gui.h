#pragma once

#include "imgui.h"
#include "ImGuizmo.h"

namespace ImGui
{

template <int N> void RatioTable(const char *title, const char *(&names)[N], int *value)
{
    if (ImGui::BeginTable(title, N))
    {
        for (int i = 0; i < N; i++)
        {
            ImGui::TableNextColumn();
            ImGui::RadioButton(names[i], value, i);
        }
        ImGui::EndTable();
    }
}

}
