#pragma once

#include <filesystem>
#include <optional>
#include <seqan3/io/sequence_file/input.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "FastqRecord.hpp"

namespace pipelines::preprocess {

struct MateName {
    std::string_view key;
    int mate{};
};

inline auto mateName(std::string_view header) -> MateName {
    const auto separator = header.find_first_of(" \t");
    auto key = header.substr(0, separator);
    int mate = 0;
    if (key.ends_with("/1") || key.ends_with("/2")) {
        mate = key.back() - '0';
        key.remove_suffix(2);
    }
    if (separator != std::string_view::npos) {
        const auto start = header.find_first_not_of(" \t", separator);
        if (start != std::string_view::npos) {
            const auto comment = header.substr(start);
            if (comment.starts_with("1:") || comment.starts_with("2:")) {
                const int label = comment.front() - '0';
                if (mate && mate != label) throw std::runtime_error("Conflicting mate labels");
                mate = label;
            }
        }
    }
    return {key, mate};
}

inline void validateMateNames(std::string_view forward, std::string_view reverse) {
    const auto first = mateName(forward);
    const auto second = mateName(reverse);
    if (first.key.empty() || first.key != second.key || first.mate == 2 || second.mate == 1)
        throw std::runtime_error("Mate identifiers do not match their pair/orientation");
}

class CheckedFastqReader {
    using Input = seqan3::sequence_file_input<>;
    using Iterator = decltype(std::declval<Input&>().begin());
    std::filesystem::path path;
    Input input;
    std::optional<Iterator> it;
    size_t ordinal{};

   public:
    explicit CheckedFastqReader(const std::filesystem::path& path) : path(path), input(path) {}
    auto next() -> std::optional<dataTypes::FastqRecord> {
        try {
            if (!it)
                it.emplace(input.begin());
            else
                ++*it;
            if (*it == input.end()) return std::nullopt;
            ++ordinal;
            return std::move(**it);
        } catch (const std::exception& error) {
            throw std::runtime_error(path.string() + ": FASTQ record " +
                                     std::to_string(ordinal + 1) + ": " + error.what());
        }
    }
};

class CheckedFastqPairReader {
    CheckedFastqReader forward;
    CheckedFastqReader reverse;
    std::filesystem::path forwardPath, reversePath;
    size_t ordinal{};

   public:
    CheckedFastqPairReader(const std::filesystem::path& first, const std::filesystem::path& second)
        : forward(first), reverse(second), forwardPath(first), reversePath(second) {}
    auto next() -> std::optional<dataTypes::PairedFastqRecords> {
        ++ordinal;
        std::optional<dataTypes::FastqRecord> first, second;
        try {
            first = forward.next();
            second = reverse.next();
            if (!first && !second) return std::nullopt;
            if (!first || !second) throw std::runtime_error("Unequal mate counts (unexpected EOF)");
            validateMateNames(first->id(), second->id());
        } catch (const std::exception& error) {
            throw std::runtime_error("FASTQ pair " + std::to_string(ordinal) + " in " +
                                     forwardPath.string() + " and " + reversePath.string() + ": " +
                                     error.what() + "; IDs: " + (first ? first->id() : "<EOF>") +
                                     " / " + (second ? second->id() : "<EOF>"));
        }
        return dataTypes::PairedFastqRecords{std::move(*first), std::move(*second)};
    }
};
}  // namespace pipelines::preprocess
