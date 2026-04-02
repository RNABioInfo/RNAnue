#include "GeneCopyMasker.hpp"

// Standard
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <ranges>
#include <span>
#include <utility>

// Seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>

// Internal
#include "GenomicFeatureGroup.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStrand.hpp"
#include "MaskedFeatureCluster.hpp"
#include "seqan3/alphabet/views/complement.hpp"

namespace pipelines::align {

auto GeneCopyMasker::addNewCluster(seqan3::dna5_vector&& representativeSequence,
                                   GenomicFeatureGroup const& featureGroup) -> void {
    const auto newClusterIndex = clusters.size();

    clusters.emplace_back(MaskedFeatureCluster{.baseFeatureGroup = featureGroup,
                                               .baseSequence = std::move(representativeSequence),
                                               .subFeatureGroups = {}});

    clusterIndicesByLength[clusters[newClusterIndex].baseSequence.size()].push_back(
        newClusterIndex);
}

auto GeneCopyMasker::extractRegionSequenceSpan(const GenomicRegion& region) const
    -> seqan3::dna5_vector {
    const auto referenceSpan =
        referenceGenome.getReferenceSequenceSpan(region.getReferenceIDIndex());

    const auto startPos = static_cast<std::size_t>(region.getStart());
    const auto regionLength = static_cast<std::size_t>(region.length());

    if (startPos >= referenceSpan.size()) {
        return {};
    }

    const auto maxLength = referenceSpan.size() - startPos;
    const auto safeLength = std::min(regionLength, maxLength);

    seqan3::dna5_vector result;
    result.reserve(safeLength);

    auto const subseq = referenceSpan.subspan(startPos, safeLength);

    if (region.getStrand() == GenomicStrand::FORWARD) {
        std::ranges::copy(subseq, std::back_inserter(result));
    } else {
        auto view = subseq | std::views::reverse | seqan3::views::complement;
        std::ranges::copy(view, std::back_inserter(result));
    }

    return result;
}
auto GeneCopyMasker::maskRegion(const GenomicRegion& region) -> void {
    auto maskSequence = seqan3::dna5_vector(region.length(), kMaskNucleotide);
    maskedGenome.modify(region.getReferenceIDIndex(), region.getRegion(), maskSequence);
}
auto GeneCopyMasker::allowedLengthRange(std::size_t length, double minIdentity)
    -> std::pair<std::size_t, std::size_t> {
    const auto identityThreshold = std::clamp(minIdentity, kMinIdentityBound, kMaxIdentityBound);
    if (identityThreshold <= kMinIdentityBound) {
        return {kZeroEdits, std::numeric_limits<std::size_t>::max()};
    }

    const auto lengthAsDouble = static_cast<double>(length);
    const auto minLength = static_cast<std::size_t>(std::floor(identityThreshold * lengthAsDouble));
    const auto maxLength = static_cast<std::size_t>(std::ceil(lengthAsDouble / identityThreshold));

    return {minLength, std::max(minLength, maxLength)};
}
auto GeneCopyMasker::maxAllowedEdits(std::size_t lengthFirst, std::size_t lengthSecond,
                                     double minIdentity) -> std::size_t {
    const auto identityThreshold = std::clamp(minIdentity, kMinIdentityBound, kMaxIdentityBound);
    const auto longestLength = static_cast<double>(std::max(lengthFirst, lengthSecond));
    const auto editsAsDouble = (kMaxIdentityBound - identityThreshold) * longestLength;

    return editsAsDouble <= kMinIdentityBound ? kZeroEdits
                                              : static_cast<std::size_t>(std::floor(editsAsDouble));
}
auto GeneCopyMasker::editDistanceToIdentity(std::size_t editCount, std::size_t lengthFirst,
                                            std::size_t lengthSecond) -> double {
    const auto longestLength = static_cast<double>(std::max(lengthFirst, lengthSecond));
    if (longestLength <= kMinIdentityBound) {
        return kFullIdentity;
    }
    return kMaxIdentityBound - (static_cast<double>(editCount) / longestLength);
}
}  // namespace pipelines::align
