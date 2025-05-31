#pragma once

// Standard
#include <cstddef>

// Internal
#include <vector>

#include "SamRecord.hpp"
#include "SplitRecords.hpp"
#include "SplitRecordsVariantGroups.hpp"

using namespace dataTypes;

class SplitRecordsConstructor {
   public:
    SplitRecordsConstructor() = delete;
    SplitRecordsConstructor(const SplitRecordsConstructor &) = delete;
    SplitRecordsConstructor(SplitRecordsConstructor &&) = delete;
    auto operator=(const SplitRecordsConstructor &) -> SplitRecordsConstructor & = delete;
    auto operator=(SplitRecordsConstructor &&) -> SplitRecordsConstructor & = delete;
    ~SplitRecordsConstructor() = delete;

    static auto getSplitRecords(const std::vector<SamRecord> &recordGroup)
        -> SplitRecordsVariantGroups;

   private:
    static auto constructSplitRecordsForHitGroup(const std::vector<SamRecord> &hitGroup)
        -> MultimericRecords;

    static auto constructSplitRecords(const SamRecord &record) -> MultimericRecords;
};
