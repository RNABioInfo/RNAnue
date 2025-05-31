#pragma once

// Standard
#include <variant>

// Internal
#include "HitGroup.hpp"
#include "HitGroupFailureReason.hpp"
#include "SplitRecords.hpp"

namespace dataTypes {

template <RecordsContainer T>
class FailedHitGroup : public HitGroup<T> {
   public:
    FailedHitGroup(T&& records, const int hitGroupTag, const HitGroupFailureReason failureReason)
        : HitGroup<T>(std::move(records), hitGroupTag), failureReason(failureReason) {}

   protected:
    HitGroupFailureReason failureReason;
};

using FailedHitGroupVariant =
    std::variant<FailedHitGroup<SingletonRecord>, FailedHitGroup<ChimericRecords>,
                 FailedHitGroup<MultimericRecords>>;

}  // namespace dataTypes
