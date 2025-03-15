#include "RenderManager.h"

#include <filesystem>

void SetCurrentDirToExe()
{
    static TCHAR buffer[MAX_PATH];
    memset(buffer, 0, sizeof(buffer));
    GetModuleFileName(0, buffer, MAX_PATH);
    auto path = std::filesystem::path(buffer).parent_path();
    std::filesystem::current_path(path);
}

int main(int argc, const char *argv[])
{
    SetCurrentDirToExe();

    GetRenderManager().MainLoop(argc, argv);
    return 0;
}
