#pragma once

// Standard
#include <algorithm>
#include <cstddef>
#include <deque>
#include <iterator>
#include <seqan3/io/sam_file/input.hpp>
#include <string>
#include <vector>

namespace pipelines::detect {

struct SamReference {
    explicit SamReference(auto samFileHeader)
        : referenceLengths(getReferenceLengths(samFileHeader)),
          referenceIDs(samFileHeader.ref_ids()) {}

    std::vector<size_t> referenceLengths;
    std::deque<std::string> referenceIDs;

   private:
    [[nodiscard]] static auto getReferenceLengths(const auto& samFileHeader)
        -> std::vector<size_t> {
        std::vector<size_t> referenceLengths{};
        std::ranges::transform(samFileHeader.ref_id_info, std::back_inserter(referenceLengths),
                               [](auto const& info) { return std::get<0>(info); });

        return referenceLengths;
    }
};

}  // namespace pipelines::detect
