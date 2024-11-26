#pragma once

#include <filesystem>
#include <string>

inline auto getTestFilePath(const std::string& fileName) -> std::string {
    return (std::filesystem::path{__FILE__}.parent_path() / "test_data" / fileName).string();
}
