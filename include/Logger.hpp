#pragma once

// Standard
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <map>
#include <mutex>
#include <source_location>
#include <string>
#include <utility>

// seqan3
#include <seqan3/core/debug_stream.hpp>

enum class LogLevel : std::uint8_t { DEBUG, INFO, WARNING, ERROR };

template <typename T>
concept IsOutputStreamable = requires() { T::operator<<(); };

constexpr auto IsErrorLogLevel(const LogLevel &logLevel) -> bool {
    return logLevel == LogLevel::ERROR;
};

struct X {};
inline constexpr auto IncludeSourceLocation = X{};

struct SourceLocation {
    SourceLocation(const SourceLocation &) = default;
    constexpr SourceLocation(std::source_location loc = std::source_location::current()) {
        auto input = loc.file_name();
        for (auto out = fileName; *input++; *out++ = *input);
        line = loc.line();
    }
    SourceLocation(SourceLocation &&) = delete;
    auto operator=(const SourceLocation &) -> SourceLocation & = default;
    auto operator=(SourceLocation &&) -> SourceLocation & = delete;
    constexpr SourceLocation(X /*unused*/,
                             std::source_location loc = std::source_location::current())
        : SourceLocation(loc) {}
    ~SourceLocation() = default;

    char fileName[256] = {};
    uint_least32_t line{};
};

class Logger {
   public:
    Logger(Logger &&) = delete;
    auto operator=(Logger &&) -> Logger & = delete;
    Logger(const Logger &) = delete;
    auto operator=(const Logger &) -> Logger & = delete;
    ~Logger() = default;

    static auto getInstance() -> Logger & {
        static Logger instance;  // Singleton instance
        return instance;
    }

    static void setLogLevel(const std::string &logLevelString) {
        const std::map<std::string, LogLevel> stringToLogLevelMap{{"debug", LogLevel::DEBUG},
                                                                  {"info", LogLevel::INFO},
                                                                  {"warning", LogLevel::WARNING},
                                                                  {"error", LogLevel::ERROR}};

        auto iterator = stringToLogLevelMap.find(logLevelString);
        if (iterator != stringToLogLevelMap.end()) {
            getInstance().logLevel = iterator->second;
        } else {
            log<IncludeSourceLocation, LogLevel::ERROR>("Invalid log level: ", logLevelString);
        }
    }

    static void setLogLevel(LogLevel level) { getInstance().logLevel = level; }

    template <LogLevel level = LogLevel::INFO, typename... Args>
        requires(not IsErrorLogLevel(level))
    static void log(Args &&...args) {
        std::lock_guard<std::mutex> lock(getInstance().logMutex);
        if (level >= getInstance().logLevel) {
            std::string levelStr;
            switch (level) {
                case LogLevel::DEBUG:
                    levelStr = "[DEBUG]";
                    break;
                case LogLevel::INFO:
                    levelStr = "[INFO]";
                    break;
                case LogLevel::WARNING:
                    levelStr = "[WARNING]";
                    break;
                case LogLevel::ERROR:
                    levelStr = "[ERROR]";
                    break;
            }
            seqan3::debug_stream << levelStr << " " << getTime() << " ";

            (seqan3::debug_stream << ... << std::forward<Args>(args)) << "\n";

            if (level == LogLevel::ERROR) {
                exit(EXIT_FAILURE);
            }
        }
    }

    template <SourceLocation Source, LogLevel level = LogLevel::INFO, typename... Args>
    static void log(Args &&...args) {
        std::lock_guard<std::mutex> lock(getInstance().logMutex);
        if (level >= getInstance().logLevel) {
            std::string levelStr;
            switch (level) {
                case LogLevel::DEBUG:
                    levelStr = "[DEBUG]";
                    break;
                case LogLevel::INFO:
                    levelStr = "[INFO]";
                    break;
                case LogLevel::WARNING:
                    levelStr = "[WARNING]";
                    break;
                case LogLevel::ERROR:
                    levelStr = "[ERROR]";
                    break;
            }
            seqan3::debug_stream << levelStr << " " << getTime() << " ";

            (seqan3::debug_stream << ... << std::forward<Args>(args)) << "; ";

            if (level == LogLevel::ERROR) {
                seqan3::debug_stream << "File: " << Source.fileName << "; Line: " << Source.line
                                     << "\n";

                exit(EXIT_FAILURE);
            }
        }
    }

   private:
    Logger() = default;

    LogLevel logLevel{LogLevel::INFO};
    std::mutex logMutex;

    static auto getTime() -> std::string {
        const auto now = std::chrono::system_clock::now();
        const std::time_t current_time = std::chrono::system_clock::to_time_t(now);

        std::ostringstream time_stream;
        time_stream << std::put_time(std::localtime(&current_time), "[%Y-%m-%d %H:%M:%S]") << " ";

        return time_stream.str();
    };
};
