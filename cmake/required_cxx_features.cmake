include_guard(GLOBAL)

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
