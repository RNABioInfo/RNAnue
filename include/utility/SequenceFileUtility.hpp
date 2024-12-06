#pragma once

// Standard
#include <cstddef>
#include <filesystem>

// seqan3
#include <seqan3/io/sequence_file/input.hpp>

namespace SequenceFileUtility {

namespace fs = std::filesystem;

auto countEntries(const fs::path& path) -> std::size_t;

auto hasAtLeastEntries(const fs::path& path, size_t requiredCount) -> bool;

}  // namespace SequenceFileUtility
