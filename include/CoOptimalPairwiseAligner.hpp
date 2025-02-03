#pragma once

// Standard
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <optional>
#include <ranges>
#include <type_traits>
#include <utility>
#include <vector>

// seqan3
#include <seqan3/alignment/configuration/align_config_debug.hpp>
#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/configuration/align_config_output.hpp>
#include <seqan3/alignment/configuration/align_config_scoring_scheme.hpp>
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/alignment/matrix/detail/aligned_sequence_builder.hpp>
#include <seqan3/alignment/matrix/detail/debug_matrix.hpp>
#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/alignment/matrix/detail/row_wise_matrix.hpp>
#include <seqan3/alignment/matrix/detail/trace_directions.hpp>
#include <seqan3/alignment/matrix/detail/trace_iterator.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/core/configuration/detail/configuration_element_debug_mode.hpp>

class CoOptimalPairwiseAligner {
   public:
    explicit CoOptimalPairwiseAligner(
        seqan3::nucleotide_scoring_scheme<int8_t> scoringScheme) noexcept
        : cfg(getAlignmentConfig(scoringScheme)) {}

    struct Result {
        int score;
        double complementarity;
        double fraction;
        std::pair<std::size_t, std::size_t> beginPositions;
        std::pair<std::size_t, std::size_t> endPositions;
    };

    template <typename sequence_pair_t>
    [[nodiscard]] auto getLocalAlignments(const sequence_pair_t &sequencePair) const noexcept
        -> std::vector<Result> {
        auto alignResults = seqan3::align_pairwise(sequencePair, cfg);

        if (alignResults.begin() == alignResults.end()) {
            return {};
        }

        const auto &alignResult = *alignResults.begin();

        std::vector<Result> results;
        constexpr size_t MAX_EXPECTED_RESULTS_COUNT = 5;
        results.reserve(MAX_EXPECTED_RESULTS_COUNT);

        using TraceMatrix =
            seqan3::detail::two_dimensional_matrix<std::optional<seqan3::detail::trace_directions>>;
        using ScoreMatrix = seqan3::detail::row_wise_matrix<std::optional<int>>;
        const TraceMatrix traceMatrix{alignResult.trace_matrix()};
        const ScoreMatrix scoreMatrix{alignResult.score_matrix()};

        const std::optional<int> score = alignResult.score();

        if (!score.has_value()) {
            return {};
        }

        auto iterator = std::ranges::find(scoreMatrix, score);
        while (iterator != scoreMatrix.end()) {
            size_t row = std::distance(scoreMatrix.begin(), iterator) / scoreMatrix.cols();
            size_t col = std::distance(scoreMatrix.begin(), iterator) % scoreMatrix.cols();

            seqan3::detail::matrix_coordinate traceBegin{seqan3::detail::row_index_type{row},
                                                         seqan3::detail::column_index_type{col}};
            results.push_back(makeResult(sequencePair, score.value(), traceBegin, traceMatrix));

            iterator = std::ranges::find(iterator + 1, scoreMatrix.end(), score);
        }

        return results;
    };

   private:
    using Configuration = seqan3::configuration<
        seqan3::align_cfg::method_local,
        seqan3::align_cfg::scoring_scheme<seqan3::nucleotide_scoring_scheme<signed char>>,
        seqan3::align_cfg::gap_cost_affine, seqan3::align_cfg::output_score,
        seqan3::align_cfg::output_begin_position, seqan3::align_cfg::output_end_position,
        seqan3::align_cfg::output_alignment,
        seqan3::detail::debug_mode<std::integral_constant<seqan3::detail::align_config_id,
                                                          seqan3::detail::align_config_id::debug>>>;

    Configuration cfg;

    [[nodiscard]] static constexpr auto getAlignmentConfig(
        const seqan3::nucleotide_scoring_scheme<int8_t> &scoringScheme) noexcept -> Configuration {
        return seqan3::align_cfg::method_local{} |
               seqan3::align_cfg::scoring_scheme{scoringScheme} |
               seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-2},
                                                  seqan3::align_cfg::extension_score{-2}} |
               seqan3::align_cfg::output_score{} | seqan3::align_cfg::output_begin_position{} |
               seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::output_alignment{} |
               seqan3::align_cfg::detail::debug{};
    };

    static auto tracePath(
        const seqan3::detail::matrix_coordinate &trace_begin,
        const seqan3::detail::two_dimensional_matrix<seqan3::detail::trace_directions>
            &complete_matrix) noexcept {
        using matrix_t = seqan3::detail::two_dimensional_matrix<seqan3::detail::trace_directions>;
        using matrix_iter_t = std::ranges::iterator_t<matrix_t const>;
        using trace_iterator_t = seqan3::detail::trace_iterator<matrix_iter_t>;
        using path_t = std::ranges::subrange<trace_iterator_t, std::default_sentinel_t>;

        seqan3::detail::matrix_offset const offset{trace_begin};

        return path_t{trace_iterator_t{complete_matrix.begin() + offset}, std::default_sentinel};
    };

    template <typename sequence_pair_t, typename score_t, typename matrix_coordinate_t>
    auto makeResult(const sequence_pair_t &sequencePair, score_t score,
                    matrix_coordinate_t endPositions, auto const &alignmentMatrix) const noexcept
        -> Result {
        const size_t elementsN = alignmentMatrix.rows() * alignmentMatrix.cols();
        std::vector<seqan3::detail::trace_directions> traceDirections;
        traceDirections.reserve(elementsN);

        std::transform(alignmentMatrix.begin(), alignmentMatrix.end(),
                       std::back_inserter(traceDirections), [](const auto &direction) {
                           return direction.value_or(seqan3::detail::trace_directions::none);
                       });

        seqan3::detail::number_rows rowsN{alignmentMatrix.rows()};
        seqan3::detail::number_cols colsN{alignmentMatrix.cols()};

        seqan3::detail::two_dimensional_matrix<seqan3::detail::trace_directions> traceMatrix(
            rowsN, colsN, traceDirections);

        using std::get;
        seqan3::detail::aligned_sequence_builder builder{get<0>(sequencePair),
                                                         get<1>(sequencePair)};

        const auto traceResult =
            builder(CoOptimalPairwiseAligner::tracePath(endPositions, traceMatrix));

        const auto &alignmentResult = traceResult.alignment;

        const int gapCount = std::ranges::count_if(get<0>(alignmentResult),
                                                   [](auto n) { return n == seqan3::gap{}; }) +
                             std::ranges::count_if(get<1>(alignmentResult),
                                                   [](auto n) { return n == seqan3::gap{}; });

        const size_t alignmentSize = get<0>(alignmentResult).size();

        assert(alignmentSize == get<1>(alignmentResult).size());

        // TODO The number to add depends on the scoring scheme and gap penalties
        const size_t matchingCount = gapCount > 0 ? ((score + alignmentSize + gapCount + 2) / 2)
                                                  : ((score + alignmentSize) / 2);
        double complementarity = 0.0;
        if (alignmentSize != 0) {
            complementarity = (static_cast<double>(matchingCount) / alignmentSize);
        }

        const int matchCount = alignmentSize - gapCount;
        const double fraction = static_cast<double>(matchCount) /
                                std::min(get<0>(sequencePair).size(), get<1>(sequencePair).size());

        return {score, complementarity, fraction,
                std::make_pair(traceResult.first_sequence_slice_positions.first,
                               traceResult.second_sequence_slice_positions.first),
                std::make_pair(endPositions.col, endPositions.row)};
    };
};
