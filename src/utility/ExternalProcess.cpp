#include "ExternalProcess.hpp"

#include <fcntl.h>
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>

#include <cerrno>
#include <csignal>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <vector>

#include "StarExecutable.hpp"

extern char** environ;

namespace utility {
namespace {
volatile sig_atomic_t interrupted = 0;
void interruptProcess(int signal) { interrupted = signal; }
struct SignalScope {
    struct sigaction oldInt{}, oldTerm{};
    SignalScope() {
        interrupted = 0;
        struct sigaction action{};
        action.sa_handler = interruptProcess;
        sigemptyset(&action.sa_mask);
        sigaction(SIGINT, &action, &oldInt);
        sigaction(SIGTERM, &action, &oldTerm);
    }
    ~SignalScope() {
        sigaction(SIGINT, &oldInt, nullptr);
        sigaction(SIGTERM, &oldTerm, nullptr);
    }
};
}  // namespace

TemporaryDirectory::TemporaryDirectory(const std::filesystem::path& parent) {
    auto pattern = (std::filesystem::absolute(parent) / ".rnanue-star-XXXXXX").string();
    std::vector<char> buffer(pattern.begin(), pattern.end());
    buffer.push_back('\0');

    if (mkdtemp(buffer.data()) == nullptr) {
        throw std::runtime_error("Cannot create STAR workspace in " + parent.string() + ": " +
                                 std::strerror(errno));
    }
    directory = buffer.data();
}
TemporaryDirectory::~TemporaryDirectory() {
    std::error_code error;
    std::filesystem::remove_all(directory, error);
}

auto runExternalProcess(const std::filesystem::path& executable,
                        const std::vector<std::string>& arguments,
                        const std::filesystem::path& workspace) -> std::string {
    if (!executable.is_absolute())
        throw std::invalid_argument("External executable must be absolute");
    const auto logPath = workspace / "process.log";
    posix_spawn_file_actions_t actions;
    posix_spawnattr_t attributes;
    int error = posix_spawn_file_actions_init(&actions);

    if (error != 0) {
        throw std::runtime_error("Cannot initialize process actions: " +
                                 std::string(std::strerror(error)));
    }

    error = posix_spawnattr_init(&attributes);
    if (error != 0) {
        posix_spawn_file_actions_destroy(&actions);
        throw std::runtime_error("Cannot initialize process attributes: " +
                                 std::string(std::strerror(error)));
    }

    struct SpawnResources {
        posix_spawn_file_actions_t* actions;
        posix_spawnattr_t* attributes;
        ~SpawnResources() {
            posix_spawn_file_actions_destroy(actions);
            posix_spawnattr_destroy(attributes);
        }
    } resources{&actions, &attributes};

    error = posix_spawn_file_actions_addopen(&actions, STDOUT_FILENO, logPath.c_str(),
                                             O_WRONLY | O_CREAT | O_TRUNC, 0600);
    if (error == 0) {
        error = posix_spawn_file_actions_adddup2(&actions, STDOUT_FILENO, STDERR_FILENO);
    }
    if (error == 0) {
        error = posix_spawnattr_setflags(&attributes, POSIX_SPAWN_SETPGROUP);
    }
    if (error == 0) {
        error = posix_spawnattr_setpgroup(&attributes, 0);
    }
    if (error != 0) {
        throw std::runtime_error("Cannot configure STAR process: " +
                                 std::string(std::strerror(error)));
    }
    std::vector<std::string> owned{executable.string()};
    owned.insert(owned.end(), arguments.begin(), arguments.end());
    std::vector<char*> argv;
    argv.reserve(owned.size());

    for (auto& argument : owned) {
        argv.push_back(argument.data());
    }

    argv.push_back(nullptr);
    SignalScope signals;
    pid_t child{};

    error = posix_spawn(&child, executable.c_str(), &actions, &attributes, argv.data(), environ);

    if (error != 0) {
        throw std::runtime_error("Cannot execute STAR: " + std::string(std::strerror(error)));
    }
    int status{};

    while (true) {
        if (interrupted != 0) {
            kill(-child, SIGKILL);
        }

        const auto result = waitpid(child, &status, 0);
        if (result == child) {
            break;
        }

        if (result < 0 && errno == EINTR) {
            continue;
        }

        if (result < 0) {
            const auto reason = std::string(std::strerror(errno));
            kill(-child, SIGKILL);
            while (waitpid(child, &status, 0) < 0 && errno == EINTR) {
            }
            throw std::runtime_error("Cannot wait for STAR: " + reason);
        }
    }
    std::ifstream log(logPath, std::ios::binary);
    const auto size = std::filesystem::file_size(logPath);
    constexpr size_t diagnosticLimit = 65536;

    if (size > diagnosticLimit) {
        log.seekg(static_cast<std::streamoff>(size - diagnosticLimit));
    }

    std::string diagnostics{std::istreambuf_iterator<char>{log}, {}};
    if (interrupted != 0) {
        throw std::runtime_error("STAR interrupted by signal " + std::to_string(interrupted) +
                                 "\n" + diagnostics);
    }

    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        const auto reason = WIFSIGNALED(status)
                                ? "signal " + std::to_string(WTERMSIG(status))
                                : "exit status " + std::to_string(WEXITSTATUS(status));
        throw std::runtime_error("STAR failed (" + reason + ")\n" + diagnostics);
    }

    return diagnostics;
}

void validateStarExecutable() {
    const auto executable = star::executablePath();
    TemporaryDirectory workspace{std::filesystem::temp_directory_path()};
    auto result = runExternalProcess(executable, {"--version"}, workspace.path());
    const auto end = result.find_last_not_of(" \r\n\t");
    result.erase(end == std::string::npos ? 0 : end + 1);
    if (result != star::version()) {
        throw std::runtime_error("Bundled STAR version mismatch: expected " +
                                 std::string(star::version()) + ", got " + result);
    }
}
}  // namespace utility
