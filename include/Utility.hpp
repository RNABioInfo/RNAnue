#pragma once

// Standard
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

// Boost
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// seqan3
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/sam_file/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>

// Internal
#include "Logger.hpp"

namespace helper {

constexpr auto RELATIVE_DIFFERENCE_FACTOR = 0.0001;

inline auto isApproxEqual(double lhs, double rhs,
                          double relativeDifferenceFactor = RELATIVE_DIFFERENCE_FACTOR) -> bool {
    const auto greaterMagnitude = std::max(std::fabs(lhs), std::fabs(rhs));
    return fabs(lhs - rhs) < relativeDifferenceFactor * greaterMagnitude;
}

inline auto vectorsApproxEqual(const std::vector<double> &first, const std::vector<double> &second,
                               double relativeDifferenceFactor = RELATIVE_DIFFERENCE_FACTOR)
    -> bool {
    if (first.size() != second.size()) {
        return false;
    }

    return std::ranges::equal(first, second, [&](double lhs, double rhs) {
        return isApproxEqual(lhs, rhs, relativeDifferenceFactor);
    });
}

void crashHandler(int signal);

template <typename Container>
concept StringContainer = std::ranges::range<Container> &&
                          std::same_as<std::ranges::range_value_t<Container>, std::string>;

namespace fs = std::filesystem;

void createTmpDir(const fs::path &subpath);
void deleteDir(const fs::path &path);

inline auto hasSuffix(const std::string &fullString, const std::string &ending) -> bool {
    if (fullString.length() >= ending.length()) {
        return (0 ==
                fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
    }
    return false;
};
template <StringContainer Container>
auto hasAnySuffix(const std::string &fullString, const Container &endings) -> bool {
    if (endings.empty()) {
        return true;
    }

    for (const auto &ending : endings) {
        if (hasSuffix(fullString, ending)) {
            return true;
        }
    }
    return false;
};

inline auto hasPrefix(const std::string &fullString, const std::string &prefix) -> bool {
    if (fullString.length() >= prefix.length()) {
        return (0 == fullString.compare(0, prefix.length(), prefix));
    }
    return false;
};

template <StringContainer Container>
auto hasAnyPrefix(const std::string &fullString, const Container &prefixes) -> bool {
    if (prefixes.empty()) {
        return true;
    }

    for (const auto &prefix : prefixes) {
        if (hasPrefix(fullString, prefix)) {
            return true;
        }
    }
    return false;
};

template <StringContainer Container = std::vector<std::string>>
auto getValidFilePaths(const fs::path &directory,
                       const Container &allowedSuffixes = Container{},  // NOLINT
                       const Container &allowedPrefixes = Container{}) -> std::vector<fs::path> {
    std::vector<fs::path> filePaths;

    for (const auto &entry : fs::directory_iterator(directory)) {
        if (!entry.is_regular_file()) {
            Logger::log<LogLevel::WARNING>("Found not supported type in directory: ", entry.path());
            continue;
        }

        if (entry.path().filename().string().front() == '.') {
            Logger::log("Ignoring hidden file: ", entry);
            continue;
        }

        Logger::log<LogLevel::DEBUG>("Found file: ", entry);

        const auto &filePathStr = entry.path().string();

        bool hasValidSuffix = hasAnySuffix(filePathStr, allowedSuffixes);
        bool hasValidPrefix = hasAnyPrefix(filePathStr, allowedPrefixes);

        if (hasValidSuffix && hasValidPrefix) {
            filePaths.push_back(entry.path());
        }
    }

    return filePaths;
};

auto getDirIfExists(const fs::path &path) -> std::optional<fs::path>;

auto getUUID() -> std::string;

void mergeSamFiles(const std::vector<fs::path> &inputPaths, const fs::path &outputPath);
void mergeFastqFiles(const std::vector<fs::path> &inputPaths, const fs::path &outputPath);

auto getFilePathsInDir(const fs::path &dir) -> std::vector<fs::path>;

void concatAndDeleteFilesInTmpDir(const fs::path &tmpDir, const fs::path &outPath);

auto countUniqueSamEntries(fs::path &path) -> std::size_t;
auto countSamEntries(fs::path &path) -> std::size_t;
auto countSamEntriesSeqAn(fs::path &path) -> std::size_t;

/** Checks if a value is contained within a specified range, with a given tolerance.
 *
 * @param value The value to check.
 * @param comparisonValue The value to compare against.
 * @param tolerance The tolerance within which the values are considered equal.
 * @return true if the value is contained within the range, false otherwise.
 **/
auto isContained(int32_t value, int32_t comparisonValue, int32_t tolerance) -> bool;

template <typename T>
auto calculateMedian(std::vector<T> values) -> T {
    std::sort(values.begin(), values.end());
    const auto size = values.size();
    if (size % 2 == 0) {
        return (values[size / 2 - 1] + values[size / 2]) / 2;
    } else {
        return values[size / 2];
    }
}

auto generateRandomHexColor() -> std::string;

auto getTime() -> std::string;  // reports the current time

class Timer {
   public:
    Timer() : start(std::chrono::high_resolution_clock::now()) {}
    Timer(const Timer &) = default;
    Timer(Timer &&) = delete;
    auto operator=(const Timer &) -> Timer & = default;
    auto operator=(Timer &&) -> Timer & = delete;
    ~Timer() { stop(); }

    void stop();

   private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start;
};
}  // namespace helper

template <class F, class T, class I, class U>
concept binary_left_foldable_impl =
    std::movable<T> && std::movable<U> && std::convertible_to<T, U> &&
    std::invocable<F &, U, std::iter_reference_t<I>> &&
    std::assignable_from<U &, std::invoke_result_t<F &, U, std::iter_reference_t<I>>>;

template <class F, class T, class I>
concept binary_left_foldable =
    std::copy_constructible<F> && std::indirectly_readable<I> &&
    std::invocable<F &, T, std::iter_reference_t<I>> &&
    std::convertible_to<std::invoke_result_t<F &, T, std::iter_reference_t<I>>,
                        std::decay_t<std::invoke_result_t<F &, T, std::iter_reference_t<I>>>> &&
    binary_left_foldable_impl<F, T, I,
                              std::decay_t<std::invoke_result_t<F &, T, std::iter_reference_t<I>>>>;

struct fold_left_fn {
    template <std::input_iterator I, std::sentinel_for<I> S, class T = std::iter_value_t<I>,
              binary_left_foldable<T, I> F>
    constexpr auto operator()(I first, S last, T init, F f) const {
        using U = std::decay_t<std::invoke_result_t<F &, T, std::iter_reference_t<I>>>;
        if (first == last) return U(std::move(init));
        U accum = std::invoke(f, std::move(init), *first);
        for (++first; first != last; ++first) {
            accum = std::invoke(f, std::move(accum), *first);
        }
        return std::move(accum);
    }

    template <std::ranges::input_range R, class T = std::ranges::range_value_t<R>,
              binary_left_foldable<T, std::ranges::iterator_t<R>> F>
    constexpr auto operator()(R &&r, T init, F f) const {
        return (*this)(std::ranges::begin(r), std::ranges::end(r), std::move(init), std::ref(f));
    }
};

inline constexpr fold_left_fn fold_left;

// inline auto meanQualityScore(std::vector<seqan3::phred42> &qualities) {
//     constexpr double PHRED_SCALE_BASE = 10;

//     const double sumQualities =
//         fold_left(qualities | std::views::transform([](seqan3::phred42 &qual) {
//                       return std::pow(PHRED_SCALE_BASE,
//                                       -static_cast<double>(qual.to_phred()) / PHRED_SCALE_BASE);
//                   }),
//                   0.0, std::plus<>());
//     const double meanErrorProbability =
//         static_cast<double>(sumQualities) / static_cast<double>(qualities.size());

//     return (-PHRED_SCALE_BASE * std::log10(meanErrorProbability));
// }

// inline auto meanQualityScore(const std::span<seqan3::phred42> &qualities) {
//     constexpr double PHRED_SCALE_BASE = 10;
//     const double sumQualities =
//         fold_left(qualities | std::views::transform([](seqan3::phred42 &qual) {
//                       return std::pow(PHRED_SCALE_BASE,
//                                       -static_cast<double>(qual.to_phred()) / PHRED_SCALE_BASE);
//                   }),
//                   0.0, std::plus<>());
//     const double meanErrorProbability =
//         static_cast<double>(sumQualities) / static_cast<double>(qualities.size());

//     return (-PHRED_SCALE_BASE * std::log10(meanErrorProbability));
// }
