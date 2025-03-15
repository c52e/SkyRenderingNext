#include "Log.h"

#include <Windows.h>

static const char *kLogLevelStrings[] = {
    "debug", "info", "warn", "error", "fatal",
};

static const wchar_t *kLogLevelStringsL[] = {
    L"debug", L"info", L"warn", L"error", L"fatal",
};

void Log(const LogLevel level, const std::string_view message)
{
    const auto levelstr = kLogLevelStrings[static_cast<uint32_t>(level)];
    const auto msg = std::format("[{}] {}\n", levelstr, message);
    std::cout << msg;
    OutputDebugStringA(msg.c_str());
}

void Log(const LogLevel level, const std::wstring_view message)
{
    const auto levelstr = kLogLevelStringsL[static_cast<uint32_t>(level)];
    const auto msg = std::format(L"[{}] {}\n", levelstr, message);
    std::wcout << msg;
    OutputDebugStringW(msg.c_str());
}

void AssertInner(const std::string_view expression, const std::source_location source)
{
    LOG_FATAL("Assertion failed: ({}) at file {} (line {})", expression, source.file_name(), source.line());
}
