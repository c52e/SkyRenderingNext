#pragma once

#include <format>
#include <iostream>
#include <source_location>

enum class LogLevel
{
    Debug = 0,
    Info,
    Warn,
    Error,
    Fatal,
};

void Log(const LogLevel level, const std::string_view message);

void Log(const LogLevel level, const std::wstring_view message);

void AssertInner(const std::string_view expression,
                  const std::source_location source = std::source_location::current());

#define LOG_DEBUG(...) Log(LogLevel::Debug, std::format(__VA_ARGS__));
#define LOG_INFO(...) Log(LogLevel::Info, std::format(__VA_ARGS__));
#define LOG_WARN(...) Log(LogLevel::Warn, std::format(__VA_ARGS__));
#define LOG_ERROR(...) Log(LogLevel::Error, std::format(__VA_ARGS__));
#define LOG_FATAL(...) Log(LogLevel::Fatal, std::format(__VA_ARGS__));

#define ASSERT(expression) (void)((!!(expression)) || (AssertInner(#expression), __debugbreak(), 0))
#define ASSERT_EQ(lhs, rhs)                                                                                            \
    (void)((!!(lhs == rhs)) ||                                                                                         \
           (AssertInner(std::format("except {} == {}, but {} != {}", #lhs, #rhs, lhs, rhs)), __debugbreak(), 0))
#define ASSERT_NE(lhs, rhs)                                                                                            \
    (void)((!!(lhs != rhs)) ||                                                                                         \
           (AssertInner(std::format("except {} != {}, but {} == {}", #lhs, #rhs, lhs, rhs)), __debugbreak(), 0))
