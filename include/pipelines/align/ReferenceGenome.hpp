#pragma once

// std
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <format>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

// internal
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "ReferenceIndexMapping.hpp"
#include "Region.hpp"
#include "seqan3/alphabet/nucleotide/dna5.hpp"

// Seqan3
#include <seqan3/io/record.hpp>
#include <seqan3/io/sequence_file/input.hpp>

namespace pipelines::align {

class ReferenceGenome {
   public:
    using selected_fields_t = seqan3::fields<seqan3::field::id, seqan3::field::seq>;
    using reference_file_t =
        seqan3::sequence_file_input<seqan3::sequence_file_input_default_traits_dna,
                                    selected_fields_t>;

    using record_t = typename reference_file_t::record_type;
    using sequence_t = std::remove_cvref_t<decltype(std::declval<record_t&>().sequence())>;

    explicit ReferenceGenome(const std::filesystem::path& referencePath) {
        auto parsed = parseReferenceFile(referencePath);
        referenceIndexMapping = std::move(parsed.mapping);
        referenceSequences = std::move(parsed.sequences);
    }

    [[nodiscard]] auto getReferenceIndexMapping() const noexcept
        -> const annotation::ReferenceIndexMapping& {
        return referenceIndexMapping;
    }

    [[nodiscard]] auto getReferenceSequences() const noexcept -> const std::vector<sequence_t>& {
        return referenceSequences;
    }

    [[nodiscard]] auto getSequenceCount() const noexcept -> std::size_t {
        return referenceSequences.size();
    }

    [[nodiscard]] auto getSequenceID(int index) const noexcept -> std::optional<std::string_view> {
        return referenceIndexMapping.findID(index);
    }

    [[nodiscard]] auto getSequence(int index) const noexcept -> const sequence_t& {
        return referenceSequences[index];
    }

    [[nodiscard]] auto findIndex(std::string_view identifier) const noexcept -> std::optional<int> {
        return referenceIndexMapping.findIndex(identifier);
    }

    [[nodiscard]] auto containsIdentifier(std::string_view identifier) const noexcept -> bool {
        return referenceIndexMapping.findIndex(identifier).has_value();
    }

    [[nodiscard]] auto getReferenceSequence(std::string_view identifier) const
        -> const sequence_t& {
        const auto idxOpt = referenceIndexMapping.findIndex(identifier);
        if (!idxOpt) {
            Logger::log<SourceLocation{}, LogLevel::ERROR>(std::format(
                "Reference sequence not found for identifier {}", std::string(identifier)));
        }
        return getReferenceSequence(*idxOpt);
    }

    [[nodiscard]] auto getReferenceSequence(int referenceIndex) const -> const sequence_t& {
        if (referenceIndex < 0 ||
            static_cast<std::size_t>(referenceIndex) >= referenceSequences.size()) {
            Logger::log<SourceLocation{}, LogLevel::ERROR>(
                std::format("Reference sequence not found for index {}", referenceIndex));
        }
        return referenceSequences[static_cast<std::size_t>(referenceIndex)];
    }

    [[nodiscard]] auto getReferenceSequenceSpan(int referenceIndex) const noexcept
        -> std::span<const seqan3::dna5> {
        const auto& sequence = getReferenceSequence(referenceIndex);
        return {sequence.data(), sequence.size()};
    }

    auto modify(std::string_view identifier, const dataTypes::Region& region,
                const sequence_t& sequence) -> void {
        const auto idxOpt = referenceIndexMapping.findIndex(identifier);
        if (!idxOpt) {
            Logger::log<SourceLocation{}, LogLevel::ERROR>(std::format(
                "Reference sequence not found for identifier {}", std::string(identifier)));
        }
        modify(*idxOpt, region, sequence);
    }

    auto modify(int referenceIndex, const dataTypes::Region& region, const sequence_t& sequence)
        -> void {
        if (region.length() != static_cast<std::int32_t>(sequence.size())) {
            Logger::log<SourceLocation{}, LogLevel::ERROR>(
                "Sequence length does not match region length.");
        }

        if (referenceIndex < 0 ||
            static_cast<std::size_t>(referenceIndex) >= referenceSequences.size()) {
            Logger::log<SourceLocation{}, LogLevel::ERROR>(
                std::format("Reference sequence not found for index {}", referenceIndex));
        }

        auto& refSeq = referenceSequences[static_cast<std::size_t>(referenceIndex)];

        const auto start = static_cast<std::size_t>(region.startPosition);
        if (start + sequence.size() > refSeq.size()) {
            Logger::log<SourceLocation{}, LogLevel::ERROR>(
                "Modify region exceeds reference sequence bounds.");
        }

        std::ranges::copy(sequence, refSeq.begin() + static_cast<std::ptrdiff_t>(start));
    }

   private:
    struct ParsedReferenceData {
        annotation::ReferenceIndexMapping mapping{
            annotation::ReferenceIndexMapping::defaultCreate()};
        std::vector<sequence_t> sequences;
    };

    annotation::ReferenceIndexMapping referenceIndexMapping{
        annotation::ReferenceIndexMapping::defaultCreate()};
    std::vector<sequence_t> referenceSequences;

    static auto parseReferenceFile(const std::filesystem::path& referencePath)
        -> ParsedReferenceData {
        Logger::log("Parsing reference file");

        reference_file_t referenceFile{referencePath};

        ParsedReferenceData out{};
        std::size_t refSeqCount = 0;

        for (auto&& record : referenceFile) {
            ++refSeqCount;

            const auto identifier = extractID(record.id());
            const int idx = out.mapping.getIndex(identifier);

            const auto uidx = static_cast<std::size_t>(idx);
            if (uidx >= out.sequences.size()) {
                out.sequences.resize(uidx + 1);
            }

            out.sequences[uidx] = std::move(record.sequence());
        }

        Logger::log(std::format("Parsed {} reference sequences.", refSeqCount));

        return out;
    }

    static auto extractID(const std::string& definitionLine) -> std::string {
        const auto pos = definitionLine.find(' ');
        if (pos == std::string::npos) {
            return definitionLine;
        }
        return definitionLine.substr(0, pos);
    }
};

}  // namespace pipelines::align
