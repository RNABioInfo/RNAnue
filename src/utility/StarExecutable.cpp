#include "StarExecutable.hpp"
#include "StarBuildConfig.hpp"

#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>
#include <unistd.h>
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

namespace utility::star {
auto enabled() noexcept -> bool { return build::enabled; }
auto version() noexcept -> std::string_view { return build::version; }

auto currentExecutablePath() -> std::filesystem::path {
#ifdef __APPLE__
    uint32_t size = 0;
    _NSGetExecutablePath(nullptr, &size);
    std::vector<char> buffer(size);
    if (_NSGetExecutablePath(buffer.data(), &size) != 0)
        throw std::runtime_error("Could not determine the RNAnue executable path");
    return std::filesystem::canonical(buffer.data());
#elif defined(__linux__)
    return std::filesystem::canonical("/proc/self/exe");
#else
    throw std::runtime_error("Bundled STAR executable discovery supports Linux and macOS only");
#endif
}

auto executablePath(const std::filesystem::path& hostExecutable) -> std::filesystem::path {
    if (!enabled()) throw std::runtime_error("This RNAnue build does not include STAR (RNANUE_BUILD_STAR=OFF)");
    if (!hostExecutable.is_absolute())
        throw std::invalid_argument("STAR discovery requires an absolute host executable path");
    const auto directory = std::filesystem::canonical(hostExecutable).parent_path();
    auto candidate = directory / build::installRelativePath;
    if (directory == std::filesystem::weakly_canonical(build::hostBuildDirectory))
        candidate = build::buildExecutable;
    std::error_code error;
    if (!std::filesystem::is_regular_file(candidate, error) || access(candidate.c_str(), X_OK) != 0)
        throw std::runtime_error("Bundled STAR " + std::string(version()) +
                                 " is missing or not executable: " + candidate.string() +
                                 ". Rebuild/install RNAnue with RNANUE_BUILD_STAR=ON.");
    return std::filesystem::canonical(candidate);
}

auto executablePath() -> std::filesystem::path { return executablePath(currentExecutablePath()); }
}  // namespace utility::star
