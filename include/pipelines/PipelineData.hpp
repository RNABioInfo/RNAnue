#pragma once

// Standard
#include <concepts>
#include <cstddef>
#include <filesystem>
#include <ranges>
#include <string>
#include <vector>

// Internal
#include "Logger.hpp"

namespace pipelines {

namespace fs = std::filesystem;

template <typename Container>
concept StringContainer = std::ranges::range<Container> &&
                          std::same_as<std::ranges::range_value_t<Container>, std::string>;

static const std::string treatmentSampleGroup = "treatment";
static const std::string controlSampleGroup = "control";

struct PipelineData {
   protected:
    static auto isHidden(const std::filesystem::path &path) -> bool {
        std::string filename = path.filename().string();
        return !filename.empty() && filename[0] == '.';
    }

    static auto getSampleName(const fs::path &filePath) -> std::string {
        return filePath.parent_path().stem().string();
    }

    static auto contains(const std::string &fullString, const std::string &substring) -> bool {
        return fullString.find(substring) != std::string::npos;
    };

    template <StringContainer Container>
    static auto containsAny(const std::string &fullString, const Container &substrings) -> bool {
        for (const auto &substring : substrings) {
            if (contains(fullString, substring)) {
                return true;
            }
        }
        return false;
    };

    static auto getSubDirectories(const fs::path &parentDir) -> std::vector<fs::path> {
        std::vector<fs::path> directories;
        for (const auto &entry : fs::directory_iterator(parentDir)) {
            if (entry.is_directory()) {
                directories.push_back(entry.path());
            } else if (entry.is_regular_file() && !isHidden(entry)) {
                Logger::log(LogLevel::WARNING,
                            "Found file in directory, which is not allowed: ", entry.path());
            }
        }
        return directories;
    };
};
}  // namespace pipelines
