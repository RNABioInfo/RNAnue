#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace utility {

class TemporaryDirectory {
   public:
    explicit TemporaryDirectory(const std::filesystem::path& parent);
    ~TemporaryDirectory();
    TemporaryDirectory(const TemporaryDirectory&) = delete;
    auto operator=(const TemporaryDirectory&) -> TemporaryDirectory& = delete;
    [[nodiscard]] auto path() const -> const std::filesystem::path& { return directory; }

   private:
    std::filesystem::path directory;
};

// Executes an absolute executable without a shell. Captured diagnostics are bounded.
[[nodiscard]] auto runExternalProcess(const std::filesystem::path& executable,
                                      const std::vector<std::string>& arguments,
                                      const std::filesystem::path& workspace) -> std::string;
void validateStarExecutable();

}  // namespace utility
