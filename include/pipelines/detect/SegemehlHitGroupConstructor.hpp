#pragma once

// Standard
#include <cassert>
#include <cstddef>
#include <iterator>
#include <memory>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

// seqan3
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
#include <seqan3/utility/views/slice.hpp>

// Internal
#include "CustomSamTags.hpp"  // IWYU pragma: keep
#include "EvaluationContext.hpp"
#include "Generator.hpp"
#include "HitGroup.hpp"
#include "HitGroupFailureReason.hpp"
#include "ReadGroupPreprocessor.hpp"
#include "SamRecord.hpp"
#include "SplitRecords.hpp"
#include "VariantUnion.hpp"
#include "seqan3/io/sam_file/sam_flag.hpp"

namespace pipelines::detect {

using namespace dataTypes;

struct SegemehlHitGroupConstructionParameters {
    size_t minimumMapQuality;
    size_t minimumFragmentLength;
    bool excludeSoftClipping;
};

class SegemehlHitGroupConstructor {
   public:
    /*
     * Constructor moves the record and extracts:
     * hit group type, expected record count and hit group tag.
     */
    explicit SegemehlHitGroupConstructor(SamRecord&& record,
                                         SegemehlHitGroupConstructionParameters parameters)
        : params(parameters),
          hitGroupType(extractHitGroupType(record)),
          expectedRecordsCount(extractExpectedSplitRecords(record)),
          hitGroupTag(extractHitGroupTag(record)) {
        records.reserve(expectedRecordsCount);
        insert(std::move(record));
    }

    void insert(SamRecord&& record) {
        // If already failed, simply add record.
        if (!isValidHitGroup()) {
            records.emplace_back(std::move(record));
            return;
        }

        if (not record.reference_position()) {
            records.emplace_back(std::move(record));
            failureReason = HitGroupFailureReason::UNMAPPED;
            return;
        }

        // Check mapping quality and fragment length first.
        if (record.mapping_quality() < params.minimumMapQuality) {
            records.emplace_back(std::move(record));
            failureReason = HitGroupFailureReason::MAPPING_QUALITY;
            return;
        }

        if (record.sequence().size() < params.minimumFragmentLength) {
            records.emplace_back(std::move(record));
            failureReason = HitGroupFailureReason::FRAGMENT_LENGTH;
            return;
        }

        const size_t expectedRecordCount = record.tags().contains("XH"_tag)
                                               ? static_cast<size_t>(record.tags().get<"XH"_tag>())
                                               : 1;

        std::vector<SamRecord> subRecords;
        subRecords.reserve(expectedRecordCount);

        auto decomposer = decomposeRecord(std::forward<SamRecord>(record));

        for (auto&& subRecord : decomposer) {
            if (subRecord.sequence().size() < params.minimumFragmentLength) {
                failureReason = HitGroupFailureReason::FRAGMENT_LENGTH;
            }

            subRecords.emplace_back(std::move(subRecord));
        }

        if (subRecords.size() != expectedRecordCount) {
            failureReason = HitGroupFailureReason::MALFORMED_RECORD;
        }

        records.insert(records.end(), std::make_move_iterator(subRecords.begin()),
                       std::make_move_iterator(subRecords.end()));
    }

    [[nodiscard]] auto getConstructedEvalContext() noexcept -> ConstructedEvaluationContextVariant {
        if (!failureReason && !hasValidSplitRecordsCount()) {
            // TODO: If read should be singleton but has more records due to being paired alignment
            // -> merge or filter

            failureReason = HitGroupFailureReason::MALFORMED_HIT_GROUP;
        }

        // TODO: Remove this as soon as paired end is handled
        if (static_cast<bool>(records.front().flag() & seqan3::sam_flag::paired)) {
            failureReason = HitGroupFailureReason::NOT_IMPLEMENTED;
        }

        switch (hitGroupType) {
            case SplitRecordType::SINGLETON:
                if (failureReason) {
                    return EvaluationContext{
                        std::make_unique<SingletonHitGroup>(
                            SingletonRecord{std::move(records.front())}, hitGroupTag),
                        std::make_tuple(*failureReason)};
                }

                return EvaluationContext{std::make_unique<SingletonHitGroup>(
                    SingletonRecord{std::move(records.front())}, hitGroupTag)};
            case SplitRecordType::CHIMERIC:
                if (failureReason) {
                    return EvaluationContext{
                        std::make_unique<ChimericHitGroup>(
                            ChimericRecords{std::move(records.front()), std::move(records.back())},
                            hitGroupTag),
                        std::make_tuple(*failureReason)};
                }

                return EvaluationContext{std::make_unique<ChimericHitGroup>(
                    ChimericRecords{std::move(records.front()), std::move(records.back())},
                    hitGroupTag)};

            case SplitRecordType::MULTIMERIC:
                // TODO: Implement multimeric hit group handling
                return EvaluationContext{
                    std::make_unique<MultimericHitGroup>(MultimericRecords{records}, hitGroupTag),
                    std::make_tuple(HitGroupFailureReason::NOT_IMPLEMENTED)};

                // if (failureReason) {
                //     return EvaluationContext{std::make_unique<MultimericHitGroup>(
                //                                  MultimericRecords{records}, hitGroupTag),
                //                              std::make_tuple(*failureReason)};
                // }

                // return EvaluationContext{
                //     std::make_unique<MultimericHitGroup>(MultimericRecords{records},
                //     hitGroupTag)};
        }

        std::unreachable();
    }

   private:
    SegemehlHitGroupConstructionParameters params;

    SplitRecordType hitGroupType;
    size_t expectedRecordsCount;
    int hitGroupTag;
    std::vector<SamRecord> records;

    std::optional<HitGroupFailureReason> failureReason = std::nullopt;

    [[nodiscard]] auto splitRecordCount() const noexcept -> size_t { return records.size(); }

    [[nodiscard]] auto isValidHitGroup() const noexcept -> bool { return !failureReason; }

    [[nodiscard]] static auto extractHitGroupType(const SamRecord& record) -> SplitRecordType {
        if (!record.tags().contains("XJ"_tag)) {
            return SplitRecordType::SINGLETON;
        }

        const int totalSplitAlignmentCount = record.tags().get<"XJ"_tag>();

        switch (totalSplitAlignmentCount) {
            case 1:
                return SplitRecordType::SINGLETON;
            case 2:
                return SplitRecordType::CHIMERIC;
            default:
                return SplitRecordType::MULTIMERIC;
        }
    }

    [[nodiscard]] static auto extractExpectedSplitRecords(const SamRecord& record) -> size_t {
        if (!record.tags().contains("XJ"_tag)) {
            return 1;
        }

        return static_cast<size_t>(record.tags().get<"XJ"_tag>());
    }

    [[nodiscard]] static auto extractHitGroupTag(const SamRecord& record) -> int {
        assert(record.tags().contains("HI"_tag) && "Record does not contain HI tag.");

        return record.tags().get<"HI"_tag>();
    }

    // Checks if the hit group type is compatible with the split record count
    [[nodiscard]] auto hasValidSplitRecordsCount() const -> bool {
        switch (hitGroupType) {
            case SplitRecordType::SINGLETON:
                return splitRecordCount() == 1;
            case SplitRecordType::CHIMERIC:
                return splitRecordCount() == 2;
            case SplitRecordType::MULTIMERIC:
                return splitRecordCount() > 2;
        }

        // Should never be reached
        std::unreachable();
    }

    auto decomposeRecord(SamRecord&& record) -> Generator<SamRecord> {
        std::vector<seqan3::cigar> currentCigar{};
        constexpr size_t initialCigarCapacity = 10;
        currentCigar.reserve(initialCigarCapacity);

        size_t referencePosition = record.reference_position().value_or(0);
        size_t startPosRead{};
        size_t endPosRead{};  // absolute position in read/alignment (e.g., 1 to end)
        size_t startPosSplit{};
        size_t endPosSplit{};  // position in split read (e.g., XX:i to XY:i / 14 to 20)

        int nextSplitReferenceShift = 0;

        auto const getCurrentSplit = [&] [[nodiscard]] () -> SamRecord {
            const auto splitSeq =
                record.sequence() | seqan3::views::slice(static_cast<ptrdiff_t>(startPosRead),
                                                         static_cast<ptrdiff_t>(endPosRead));
            const auto splitQual =
                record.base_qualities() | seqan3::views::slice(static_cast<ptrdiff_t>(startPosRead),
                                                               static_cast<ptrdiff_t>(endPosRead));

            seqan3::sam_tag_dictionary tags = record.tags();
            tags.get<"XX"_tag>() = static_cast<int>(startPosSplit);
            tags.get<"XY"_tag>() = static_cast<int>(endPosSplit);
            tags.get<"XN"_tag>() = static_cast<int>(splitRecordCount());

            return SamRecord{record.id(),
                             record.flag(),
                             record.reference_id(),
                             referencePosition,
                             record.mapping_quality(),
                             currentCigar,
                             seqan3::dna5_vector(splitSeq.begin(), splitSeq.end()),
                             std::vector<seqan3::phred42>(splitQual.begin(), splitQual.end()),
                             std::move(tags)};
        };

        auto const addOtherCigar = [&](const auto& cigar) {
            const auto cigarValue = get<0>(cigar);
            endPosRead += cigarValue;
            endPosSplit += cigarValue;
            currentCigar.push_back(cigar);
        };

        auto const addInsertionCigar = [&](const auto& cigar) {
            const auto cigarValue = get<0>(cigar);
            endPosRead += cigarValue;
            endPosSplit += cigarValue;
            nextSplitReferenceShift -= cigarValue;
            currentCigar.push_back(cigar);
        };

        auto const addDeletionCigar = [&](const auto& cigar) {
            currentCigar.push_back(cigar);
            nextSplitReferenceShift += get<0>(cigar);
        };

        auto const addSoftClipCigar = [&](const auto& cigar) {
            const auto cigarValue = get<0>(cigar);
            if (!params.excludeSoftClipping) {
                nextSplitReferenceShift -= cigarValue;
                addOtherCigar(cigar);
                return;
            }

            /* If current cigar is empty, we are at the beginning of the read in case
            of soft clipping at the end of the read it is just ignored */
            if (currentCigar.empty()) {
                nextSplitReferenceShift -= cigarValue;
                startPosRead += cigarValue;
                endPosRead += cigarValue;
                startPosSplit += cigarValue;
                endPosSplit += cigarValue;
            }
        };

        auto const addSkipCigar = [&] [[nodiscard]] (const auto& cigar) -> SamRecord {
            if (currentCigar.empty()) {
                return {};
            }

            SamRecord currentSplit = getCurrentSplit();

            // Set up positions for the next split
            const auto cigarValue = get<0>(cigar);
            assert((cigarValue + endPosRead + nextSplitReferenceShift) >= 0);
            referencePosition += cigarValue + endPosRead + nextSplitReferenceShift;
            startPosSplit = endPosSplit;
            startPosRead = endPosRead;
            nextSplitReferenceShift = 0;
            currentCigar.clear();

            return currentSplit;
        };

        for (const auto& cigar : record.cigar_sequence()) {
            if (cigar == 'M'_cigar_operation || cigar == '='_cigar_operation ||
                cigar == 'X'_cigar_operation) {
                addOtherCigar(cigar);
            } else if (cigar == 'I'_cigar_operation) {
                addInsertionCigar(cigar);
            } else if (cigar == 'D'_cigar_operation) {
                addDeletionCigar(cigar);
            } else if (cigar == 'S'_cigar_operation) {
                addSoftClipCigar(cigar);
            } else if (cigar == 'N'_cigar_operation) {
                co_yield addSkipCigar(cigar);
            }
        }

        co_yield getCurrentSplit();
    };
};

}  // namespace pipelines::detect
