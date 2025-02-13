#pragma once

// Standard
#include <cassert>
#include <concepts>
#include <cstddef>
#include <iterator>
#include <memory>
#include <ranges>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

// seqan3
#include <seqan3/contrib/parallel/buffer_queue.hpp>
#include <seqan3/core/range/detail/adaptor_from_functor.hpp>
#include <seqan3/io/sam_file/input.hpp>

// Internal
#include "SamRecord.hpp"

using record_input_t =
    seqan3::sam_file_input<seqan3::sam_file_input_default_traits<>, dataTypes::SamFieldIDs>;

template <std::ranges::range urng_t>
class AsyncSplitRecordGroupBufferView
    : public std::ranges::view_interface<AsyncSplitRecordGroupBufferView<urng_t>> {
   private:
    static_assert(std::ranges::input_range<urng_t>,
                  "The range parameter to async_input_buffer_view must be at least a "
                  "std::ranges::input_range.");
    static_assert(std::ranges::view<urng_t>,
                  "The range parameter to async_input_buffer_view must model std::ranges::view.");
    static_assert(std::movable<std::ranges::range_value_t<urng_t>>,
                  "The range parameter to async_input_buffer_view must have a value_type that is "
                  "std::movable.");
    static_assert(
        std::constructible_from<std::ranges::range_value_t<urng_t>,
                                std::remove_reference_t<std::ranges::range_reference_t<urng_t>>&&>,
        "The range parameter to async_input_buffer_view must have a value_type that is "
        "constructible by a moved "
        "value of its reference type.");
    static_assert(std::same_as<urng_t, std::ranges::ref_view<record_input_t>>,
                  "Range type must be sam record type");

    using urng_iterator_type = std::ranges::iterator_t<urng_t>;

    struct state {
        urng_t urange;

        seqan3::contrib::fixed_buffer_queue<std::vector<record_input_t::value_type>> buffer;

        std::thread producer;
    };

    std::shared_ptr<state> statePtr = nullptr;

    class iterator {
        using sentinel_type = std::default_sentinel_t;

        seqan3::contrib::fixed_buffer_queue<std::vector<record_input_t::value_type>>* bufferPtr =
            nullptr;

        mutable std::vector<record_input_t::value_type> cached_value;

        bool at_end = false;

       public:
        using difference_type = std::iter_difference_t<urng_iterator_type>;
        using value_type = std::vector<record_input_t::value_type>;
        using pointer = std::vector<record_input_t::value_type>*;
        using reference = std::vector<record_input_t::value_type>&;
        using iterator_category = std::input_iterator_tag;
        using iterator_concept = iterator_category;

        iterator() = default;  //!< Defaulted.
        // TODO: delete:
        iterator(iterator const& rhs) = default;  //!< Defaulted.
        iterator(iterator&& rhs) = default;       //!< Defaulted.
        // TODO: delete:
        auto operator=(iterator const& rhs) -> iterator& = default;  //!< Defaulted.
        auto operator=(iterator&& rhs) -> iterator& = default;       //!< Defaulted.
        ~iterator() noexcept = default;                              //!< Defaulted.

        iterator(seqan3::contrib::fixed_buffer_queue<std::vector<record_input_t::value_type>>&
                     buffer) noexcept
            : bufferPtr{&buffer} {
            ++(*this);
        }

        auto operator*() const noexcept -> reference { return cached_value; }

        auto operator->() const noexcept -> pointer { return std::addressof(cached_value); }

        auto operator++() noexcept -> iterator& {
            if (at_end) {
                return *this;
            }

            assert(bufferPtr != nullptr);

            if (bufferPtr->wait_pop(cached_value) == seqan3::contrib::queue_op_status::closed) {
                at_end = true;
            }

            return *this;
        }

        void operator++(int) noexcept { ++(*this); }

        friend constexpr auto operator==(iterator const& lhs,
                                         std::default_sentinel_t const&) noexcept -> bool {
            return lhs.at_end;
        }

        //!\copydoc operator==
        friend constexpr auto operator==(std::default_sentinel_t const&,
                                         iterator const& rhs) noexcept -> bool {
            return rhs == std::default_sentinel_t{};
        }

        //!\brief Compares for inequality with sentinel.
        friend constexpr auto operator!=(iterator const& lhs,
                                         std::default_sentinel_t const&) noexcept -> bool {
            return !(lhs == std::default_sentinel_t{});
        }

        //!\copydoc operator!=
        friend constexpr auto operator!=(std::default_sentinel_t const&,
                                         iterator const& rhs) noexcept -> bool {
            return rhs != std::default_sentinel_t{};
        }
    };

   public:
    AsyncSplitRecordGroupBufferView(urng_t _urng, size_t const bufferSize) {
        auto deleter = [](state* pointer) {
            if (pointer != nullptr) {
                pointer->buffer.close();
                pointer->producer.join();
                delete pointer;
            }
        };

        statePtr = std::shared_ptr<state>(
            new state{std::move(_urng),
                      seqan3::contrib::fixed_buffer_queue<std::vector<record_input_t::value_type>>{
                          bufferSize},
                      std::thread{}},  // thread is set/started below, needs rest of state
            deleter);

        auto runner = [&state = *statePtr]() {
            std::vector<record_input_t::value_type> recordGroup;
            recordGroup.reserve(10);

            std::string currentRecordGroupID;

            bool firstRecord = true;

            for (auto&& record : state.urange) {
                if (firstRecord) {
                    firstRecord = false;
                    currentRecordGroupID = record.id();
                }

                if (currentRecordGroupID == record.id()) {
                    recordGroup.emplace_back(std::move(record));
                } else {
                    const auto status = state.buffer.wait_push(recordGroup);
                    recordGroup.clear();
                    currentRecordGroupID = record.id();
                    recordGroup.emplace_back(std::move(record));

                    if (status == seqan3::contrib::queue_op_status::closed) {
                        break;
                    }
                }
            }

            if (!recordGroup.empty() && !state.buffer.is_closed()) {
                state.buffer.wait_push(recordGroup);
            }

            state.buffer.close();
        };

        statePtr->producer = std::thread{runner};
    }

    template <typename other_urange_t>
        requires(!std::same_as<std::remove_cvref_t<other_urange_t>,
                               AsyncSplitRecordGroupBufferView>) &&
                std::ranges::viewable_range<other_urange_t> &&
                std::constructible_from<
                    urng_t, std::ranges::ref_view<std::remove_reference_t<other_urange_t>>>
    AsyncSplitRecordGroupBufferView(other_urange_t&& _urng, size_t const bufferSize)
        : AsyncSplitRecordGroupBufferView{std::views::all(_urng), bufferSize} {}

    auto begin() -> iterator {
        assert(statePtr != nullptr);
        return {statePtr->buffer};
    }

    auto begin() const -> iterator = delete;

    auto end() -> std::default_sentinel_t { return std::default_sentinel; }

    [[nodiscard]] auto end() const -> std::default_sentinel_t = delete;
};

template <std::ranges::viewable_range urng_t>
AsyncSplitRecordGroupBufferView(urng_t&&, size_t buffer_size)
    -> AsyncSplitRecordGroupBufferView<std::views::all_t<urng_t>>;

struct AsyncSplitRecordGroupBufferViewFn {
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

        return AsyncSplitRecordGroupBufferView{std::forward<urng_t>(urange), buffer_size};
    }
};

inline constexpr auto AsyncSplitRecordGroupBuffer = AsyncSplitRecordGroupBufferViewFn{};
