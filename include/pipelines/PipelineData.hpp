#pragma once

// Standard
#include <algorithm>
#include <concepts>
#include <cstddef>
#include <filesystem>
#include <ranges>
#include <string>
#include <vector>

// Internal
#include "LogLevel.hpp"
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
    static auto validateDirEmpty(const std::filesystem::path &path, const bool forceOverwrite) {
        if (fs::exists(path)) {
            Logger::log<LogLevel::WARNING>(
                "Output directory already exists. Results will be overwritten.");

            if (forceOverwrite) {
                Logger::log<LogLevel::WARNING>("Force overwriting existing output");
                fs::remove_all(path);
            } else {
                Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
                    "Cancelling. If you want to overwrite existing results set the -f or --force "
                    "flag.");
            }
        }
    }

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
        return std::ranges::any_of(substrings, [&fullString](const std::string &substring) {
            return contains(fullString, substring);
        });
    };

    static auto getSubDirectories(const fs::path &parentDir) -> std::vector<fs::path> {
        std::vector<fs::path> directories;
        for (const auto &entry : fs::directory_iterator(parentDir)) {
            if (entry.is_directory()) {
                directories.push_back(entry.path());
            } else if (entry.is_regular_file() && !isHidden(entry)) {
                Logger::log<LogLevel::WARNING>(
                    "Found file in directory, expected only sub-directories: ", entry.path());
            }
        }
        return directories;
    };
};
}  // namespace pipelines
