#include "SequenceFileUtility.hpp"

// Standard
#include <cstddef>
#include <filesystem>
#include <iterator>

// seqan3
#include <seqan3/io/record.hpp>
#include <seqan3/io/sequence_file/input.hpp>

namespace SequenceFileUtility {

auto countEntries(const fs::path& path) -> std::size_t {
    seqan3::sequence_file_input fin{path.string(), seqan3::fields<>{}};
    return std::ranges::distance(fin.begin(), fin.end());
};

auto hasAtLeastEntries(const fs::path& path, size_t requiredCount) -> bool {
    seqan3::sequence_file_input fin{path, seqan3::fields<>{}};

    size_t count{};

    for ([[maybe_unused]] auto const& entry : fin) {
        count++;

        if (count >= requiredCount) {
            return true;
        }
    }

    return false;
};

}  // namespace SequenceFileUtility
