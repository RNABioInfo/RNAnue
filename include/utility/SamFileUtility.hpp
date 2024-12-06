#pragma once

// Standard
#include <cstddef>
#include <filesystem>

namespace SamFileUtility {

namespace fs = std::filesystem;

auto countEntries(fs::path& path) -> std::size_t;

auto hasAtLeastEntries(fs::path& path, size_t requiredCount) -> bool;

}  // namespace SamFileUtility
