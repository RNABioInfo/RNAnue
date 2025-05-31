#pragma once

// Standard
#include <filesystem>

namespace pipelines::detect {

namespace fs = std::filesystem;

struct TempOutputDirs {
    fs::path outputTmpSplitsDir;
    fs::path outputTmpMultisplitsDir;
    fs::path outputTmpUnassignedSingletonDir;
};

}  // namespace pipelines::detect
