#pragma once
#include <memory>
#include <optional>
#include <ranges>
#include <seqan3/core/range/detail/adaptor_from_functor.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <vector>

#include "SamRecord.hpp"
#include "utility/ConcurrentInput.hpp"

using record_input_t =
    seqan3::sam_file_input<seqan3::sam_file_input_default_traits<>, dataTypes::SamFieldIDs>;

template <std::ranges::range urng_t>
class AsyncSplitReadGroupBufferView
    : public std::ranges::view_interface<AsyncSplitReadGroupBufferView<urng_t>> {
    using Group = std::vector<record_input_t::value_type>;
    std::shared_ptr<utility::ConcurrentInput<Group>> input;

   public:
    AsyncSplitReadGroupBufferView(urng_t range, size_t capacity) {
        input = std::make_shared<utility::ConcurrentInput<Group>>(
            capacity,
            [range = std::move(range),
             it = std::optional<std::ranges::iterator_t<urng_t>>{}]() mutable
                -> std::optional<Group> {
                if (!it) it.emplace(range.begin());
                Group group;
                while (*it != range.end()) {
                    auto& record = **it;
                    if (!group.empty() && group.front().id() != record.id()) break;
                    group.push_back(std::move(record));
                    ++*it;
                }
                if (group.empty()) return std::nullopt;
                return group;
            });
    }
    template <typename R>
        requires(!std::same_as<std::remove_cvref_t<R>, urng_t>)
    AsyncSplitReadGroupBufferView(R&& range, size_t capacity)
        : AsyncSplitReadGroupBufferView(std::views::all(range), capacity) {}
    auto begin() { return input->begin(); }
    auto end() const { return std::default_sentinel; }
    void cancel() { input->cancel(); }
};

template <std::ranges::viewable_range urng_t>
AsyncSplitReadGroupBufferView(urng_t&&, size_t buffer_size)
    -> AsyncSplitReadGroupBufferView<std::views::all_t<urng_t>>;

struct AsyncSplitReadGroupBufferViewFn {
    constexpr auto operator()(size_t const bufferSize) const {
        return seqan3::detail::adaptor_from_functor{*this, bufferSize};
    }

    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t&& urange, size_t const buffer_size) const {
        static_assert(std::ranges::input_range<urng_t>,
                      "The range parameter to views::async_input_buffer must be at least a "
                      "std::ranges::input_range.");
        static_assert(std::ranges::viewable_range<urng_t>,
                      "The range parameter to views::async_input_buffer cannot be a temporary of a "
                      "non-view range.");
        static_assert(std::movable<std::ranges::range_value_t<urng_t>>,
                      "The range parameter to views::async_input_buffer must have a value_type "
                      "that is std::movable.");
        static_assert(std::constructible_from<
                          std::ranges::range_value_t<urng_t>,
                          std::remove_reference_t<std::ranges::range_reference_t<urng_t>>&&>,
                      "The range parameter to views::async_input_buffer must have a value_type "
                      "that is constructible by a moved "
                      "value of its reference type.");

        if (buffer_size == 0) {
            throw std::invalid_argument{
                "The buffer_size parameter to views::async_input_buffer must be > 0."};
        }

        return AsyncSplitReadGroupBufferView{std::forward<urng_t>(urange), buffer_size};
    }
};

inline constexpr auto AsyncSplitReadGroupBuffer = AsyncSplitReadGroupBufferViewFn{};
