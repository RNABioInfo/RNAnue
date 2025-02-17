#pragma once

// Standard
#include <cstddef>
#include <iterator>
#include <ranges>
#include <vector>

// seqan3
#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/range/detail/adaptor_from_functor.hpp>
#include <seqan3/io/sam_file/sam_flag.hpp>

/*!
 * @brief An iterator that wraps an underlying random access iterator for nucleotide sequences.
 *
 * This iterator can traverse the given nucleotide sequence in either forward or reverse order,
 * depending on the SAM flag provided at construction time.
 * In reverse mode the iterator acts much like std::reverse_iterator.
 *
 * @tparam urng_t Type of the underlying range, which must be a random_access_range.
 */
template <std::ranges::random_access_range urng_t>
class UnderlyingSequenceIterator {
   public:
    using base_t = std::ranges::iterator_t<urng_t>;
    using value_type = typename std::iterator_traits<base_t>::value_type;
    using difference_type = typename std::iterator_traits<base_t>::difference_type;
    using pointer = typename std::iterator_traits<base_t>::pointer;
    using reference = value_type;
    // Change from forward to bidirectional iterator tag (since we use -- in reverse mode).
    using iterator_category = std::forward_iterator_tag;

    // Default constructor.
    UnderlyingSequenceIterator() = default;

    /*!
     * @brief Constructs an iterator wrapping the given range.
     *
     * In forward mode, initializes to the container’s beginning.
     * In reverse mode, initializes to the container’s end.
     */
    explicit UnderlyingSequenceIterator(urng_t &seq, const seqan3::sam_flag &recordFlags)
        : is_on_reverse{static_cast<bool>(seqan3::sam_flag::on_reverse_strand & recordFlags)} {
        if (!is_on_reverse) {
            iter = std::ranges::begin(seq);
        } else {
            iter = std::ranges::end(seq);
        }
    }

    /*!
     * @brief Constructs an iterator from an already available base iterator.
     *
     * Used by the view to create the end iterator.
     */
    UnderlyingSequenceIterator(base_t iter, bool is_reverse)
        : iter(iter), is_on_reverse(is_reverse) {}

    /*!
     * @brief Pre-increment operator.
     *
     * In forward mode, advances to the next element.
     * In reverse mode, moves to the previous element.
     */
    inline auto operator++() -> UnderlyingSequenceIterator & {
        if (!is_on_reverse) {
            ++iter;
        } else {
            --iter;
        }
        return *this;
    }

    /*!
     * @brief Post-increment operator.
     */
    inline auto operator++(int) -> UnderlyingSequenceIterator {
        UnderlyingSequenceIterator tmp = *this;
        ++(*this);
        return tmp;
    }

    /*!
     * @brief Dereference operator.
     *
     * In forward mode, returns the current element.
     * In reverse mode, returns the complement of the element preceding the current iterator.
     */
    inline auto operator*() const -> value_type {
        if (!is_on_reverse) {
            return *iter;
        }

        auto tmp = iter;
        --tmp;
        return seqan3::complement(*tmp);
    }

    /*!
     * @brief Equality comparison operator.
     *
     * Two iterators compare equal if they have the same reverse-mode flag and their underlying
     * iterators are equal.
     */
    inline auto operator==(UnderlyingSequenceIterator const &rhs) const -> bool {
        return (is_on_reverse == rhs.is_on_reverse) && (iter == rhs.iter);
    }
    inline auto operator!=(UnderlyingSequenceIterator const &rhs) const -> bool {
        return !(*this == rhs);
    }

    /// Returns the underlying base iterator.
    [[nodiscard]] inline auto base() const -> base_t { return iter; }

    /// Queries whether the iterator is in reverse mode.
    [[nodiscard]] constexpr auto is_reverse() const -> bool { return is_on_reverse; }

   private:
    base_t iter{};              ///< The underlying iterator.
    bool is_on_reverse{false};  ///< Indicates whether reverse mode is active.
};

static_assert(std::forward_iterator<UnderlyingSequenceIterator<std::vector<seqan3::dna5>>>);

/*!
 * @brief A view that wraps a nucleotide sequence and exposes an iterator that can traverse it
 * either in forward or reverse order.
 *
 * The view is constructed from a viewable range; it always wraps the underlying range via
 * std::views::all so that if an lvalue is provided the view reflects modifications to the original
 * container.
 *
 * @tparam urng_t Type of the underlying range (after applying std::views::all).
 */
template <std::ranges::random_access_range urng_t>
class UnderlyingSequenceView {
   public:
    using iterator = UnderlyingSequenceIterator<urng_t>;
    using const_iterator = UnderlyingSequenceIterator<urng_t const>;

    // Remove the by-value constructor – use only the viewable_range overload.
    template <std::ranges::viewable_range orng_t>
    explicit UnderlyingSequenceView(orng_t &&seq, const seqan3::sam_flag &recordFlags)
        : seq(std::views::all(std::forward<orng_t>(seq))),
          is_on_reverse{static_cast<bool>(seqan3::sam_flag::on_reverse_strand & recordFlags)} {}

    /*!
     * @brief Returns an iterator to the beginning of the sequence.
     *
     * In forward mode, this is the first element.
     * In reverse mode, the iterator is initialized to the container's end.
     */
    inline auto begin() -> iterator {
        return iterator{seq,
                        is_on_reverse ? seqan3::sam_flag::on_reverse_strand : seqan3::sam_flag{}};
    }

    inline auto begin() const -> const_iterator {
        return const_iterator{
            seq, is_on_reverse ? seqan3::sam_flag::on_reverse_strand : seqan3::sam_flag{}};
    }

    inline auto cbegin() const -> const_iterator { return begin(); }

    /*!
     * @brief Returns an iterator representing the end of the sequence.
     *
     * For forward iteration this wraps std::ranges::end(seq);
     * For reverse iteration this wraps std::ranges::begin(seq).
     */
    inline auto end() -> iterator {
        if (!is_on_reverse) {
            return iterator{std::ranges::end(seq), false};
        }
        return iterator{std::ranges::begin(seq), true};
    }

    inline auto end() const -> const_iterator {
        if (!is_on_reverse) {
            return const_iterator{std::ranges::end(seq), false};
        }
        return const_iterator{std::ranges::begin(seq), true};
    }

    inline auto cend() const -> const_iterator { return end(); }

    [[nodiscard]] inline constexpr auto size() const -> size_t { return std::ranges::size(seq); }

   private:
    decltype(std::views::all(std::declval<urng_t>())) seq;  ///< The underlying nucleotide sequence.
    bool is_on_reverse{false};                              ///< Indicates reverse-mode iteration.
};

template <std::ranges::viewable_range orng_t>
UnderlyingSequenceView(orng_t &&, const seqan3::sam_flag &recordFlags)
    -> UnderlyingSequenceView<std::views::all_t<orng_t>>;

struct UnderlyingSequenceViewFn {
    constexpr auto operator()(seqan3::sam_flag recordFlags) const {
        return seqan3::detail::adaptor_from_functor{*this, recordFlags};
    }

    template <std::ranges::input_range urng_t>
    auto operator()(urng_t &&seq, seqan3::sam_flag recordFlags) const {
        static_assert(std::ranges::input_range<urng_t>,
                      "The range parameter to views::underlying_sequence must be at least a "
                      "std::ranges::input_range.");
        static_assert(std::ranges::viewable_range<urng_t>,
                      "The range parameter to views::underlying_sequence cannot be a temporary of "
                      "a non-view range.");

        return UnderlyingSequenceView{std::forward<urng_t>(seq), recordFlags};
    }
};

namespace views {
inline constexpr UnderlyingSequenceViewFn underlying_sequence{};
}
