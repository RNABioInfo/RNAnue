#pragma once

#include <filesystem>
#include <string_view>

namespace utility::star {

[[nodiscard]] auto enabled() noexcept -> bool;
[[nodiscard]] auto version() noexcept -> std::string_view;
[[nodiscard]] auto currentExecutablePath() -> std::filesystem::path;
// The explicit anchor also allows launchers/tests to resolve a relocated installation.
// Never searches PATH or the current working directory.
[[nodiscard]] auto executablePath(const std::filesystem::path& hostExecutable) -> std::filesystem::path;
[[nodiscard]] auto executablePath() -> std::filesystem::path;

}  // namespace utility::star
