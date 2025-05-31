#pragma once

#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdlib>
#include <ostream>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

// Internal
#include "Generator.hpp"
#include "SplitRecords.hpp"

namespace dataTypes {

template <RecordsContainer T>
class HitGroup {
   public:
    HitGroup(T&& records, const int hitGroupTag)
        : records(std::move(records)), hitGroupTag(hitGroupTag) {}

    HitGroup(const T& records, const int hitGroupTag)
        : records(records), hitGroupTag(hitGroupTag) {}

    using RecordsContainerT = T;

    static constexpr SplitRecordType hitGroupSplitRecordType{T::splitRecordType};
    [[nodiscard]] constexpr auto getHitGroupTag() const noexcept -> int { return hitGroupTag; };
    [[nodiscard]] constexpr auto getRecordCount() const noexcept -> size_t {
        return records.size();
    }
    [[nodiscard]] constexpr auto getRecordContainer() const -> const T& { return records; }
    [[nodiscard]] constexpr auto getRecordContainer() -> T& { return records; }
    [[nodiscard]] constexpr auto readID() const noexcept -> const std::string& {
        return records.recordID();
    }

    [[nodiscard]] auto getDeconstructedHitGroups() -> std::vector<HitGroup<SingletonRecord>>
        requires(!std::same_as<T, SingletonRecord>)
    {
        std::vector<HitGroup<SingletonRecord>> deconstructedHitGroups;
        deconstructedHitGroups.reserve(getRecordCount());

        for (const auto& record : getRecordContainer()) {
            deconstructedHitGroups.emplace_back(SingletonRecord{record}, hitGroupTag);
        }

        return deconstructedHitGroups;
    }

   protected:
    T records;
    int hitGroupTag;
};

template <typename RecordT>
auto operator<<(std::ostream& out, const HitGroup<RecordT>& hitGroup) -> std::ostream& {
    out << "RecordID: " << hitGroup.readID() << ", Type: "
        << std::remove_cvref_t<decltype(hitGroup.getRecordContainer())>::splitRecordType
        << ", Tag: " << hitGroup.getHitGroupTag() << ", Records: " << hitGroup.getRecordContainer();

    return out;
}

using SingletonHitGroup = HitGroup<SingletonRecord>;
using ChimericHitGroup = HitGroup<ChimericRecords>;
using MultimericHitGroup = HitGroup<MultimericRecords>;

using HitGroupVariant =
    std::variant<HitGroup<SingletonRecord>, HitGroup<ChimericRecords>, HitGroup<MultimericRecords>>;

}  // namespace dataTypes
