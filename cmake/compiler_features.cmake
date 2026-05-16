include(CheckCXXSourceCompiles)

if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(FATAL_ERROR
        "RNAnue requires GCC/G++ 14 or newer with C++23 standard-library support. "
        "The selected C++ compiler is '${CMAKE_CXX_COMPILER_ID}' (${CMAKE_CXX_COMPILER}). "
        "Reconfigure a fresh build directory with "
        "-DCMAKE_CXX_COMPILER=/path/to/g++-14 -DCMAKE_C_COMPILER=/path/to/gcc-14.")
endif()

if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 14)
    message(FATAL_ERROR
        "RNAnue requires GCC/G++ 14 or newer. The selected compiler is "
        "GNU ${CMAKE_CXX_COMPILER_VERSION} (${CMAKE_CXX_COMPILER}). "
        "GCC 13 lacks required C++23 libstdc++ APIs used by RNAnue, including "
        "<print>, std::forward_like, std::ranges::to, and std::ranges::zip_view. "
        "Reconfigure a fresh build directory with "
        "-DCMAKE_CXX_COMPILER=/path/to/g++-14 -DCMAKE_C_COMPILER=/path/to/gcc-14.")
endif()

set(RNANUE_REQUIRED_CXX_FEATURES_SOURCE [[
#include <coroutine>
#include <format>
#include <print>
#include <ranges>
#include <source_location>
#include <span>
#include <string>
#include <utility>
#include <vector>

enum class rnanue_probe_enum
{
    value
};

template <>
struct std::formatter<rnanue_probe_enum>
{
    constexpr auto parse(auto & ctx)
    {
        return ctx.begin();
    }

    auto format(rnanue_probe_enum, auto & ctx) const
    {
        return std::format_to(ctx.out(), "{}", "value");
    }
};

struct rnanue_probe_generator
{
    struct promise_type
    {
        auto get_return_object() -> rnanue_probe_generator
        {
            return {};
        }

        auto initial_suspend() noexcept -> std::suspend_always
        {
            return {};
        }

        auto final_suspend() noexcept -> std::suspend_always
        {
            return {};
        }

        auto yield_value(int) noexcept -> std::suspend_always
        {
            return {};
        }

        void return_void() noexcept {}
        void unhandled_exception() {}
    };
};

consteval auto rnanue_probe_consteval() -> int
{
    return 14;
}

auto rnanue_probe_coroutine() -> rnanue_probe_generator
{
    co_yield rnanue_probe_consteval();
}

auto main() -> int
{
    std::vector<int> xs{1, 2};
    std::vector<int> ys{3, 4};

    auto formatted = std::format("RNAnue {}", rnanue_probe_enum::value);
    std::print("{}", formatted);

    int value = 0;
    auto && forwarded = std::forward_like<int &&>(value);
    (void)forwarded;

    auto strings = xs
                 | std::views::transform([](int elem) { return std::to_string(elem); })
                 | std::ranges::to<std::vector>();

    for (auto pair : std::ranges::zip_view(xs, ys))
    {
        (void)pair;
    }

    std::span<int> span{xs};
    auto location = std::source_location::current();
    auto generator = rnanue_probe_coroutine();

    (void)strings;
    (void)span;
    (void)location;
    (void)generator;

    if (xs.empty())
    {
        std::unreachable();
    }

    return 0;
}
]])

set(RNANUE_PREVIOUS_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -std=c++23")
check_cxx_source_compiles("${RNANUE_REQUIRED_CXX_FEATURES_SOURCE}"
                          RNANUE_HAS_REQUIRED_CXX_FEATURES)
set(CMAKE_REQUIRED_FLAGS "${RNANUE_PREVIOUS_REQUIRED_FLAGS}")
unset(RNANUE_PREVIOUS_REQUIRED_FLAGS)

if(NOT RNANUE_HAS_REQUIRED_CXX_FEATURES)
    message(FATAL_ERROR
        "The selected GCC/G++ compiler does not support the C++20/23 language and "
        "standard-library features required by RNAnue. Required features include "
        "std::format, std::formatter, std::print, std::forward_like, std::ranges::to, "
        "std::ranges::zip_view, std::unreachable, std::source_location, std::span, "
        "consteval, and coroutines. Use GCC/G++ 14 or newer and reconfigure a fresh "
        "build directory with -DCMAKE_CXX_COMPILER=/path/to/g++-14 "
        "-DCMAKE_C_COMPILER=/path/to/gcc-14.")
endif()

message(STATUS "RNAnue C++ compiler feature check passed")
