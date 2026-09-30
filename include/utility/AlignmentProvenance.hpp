#pragma once

#include <algorithm>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "AlignmentBackend.hpp"

namespace utility {
using pipelines::align::AlignmentBackend;
inline constexpr auto provenancePrefix = "RNAnue alignment backend=";

template <typename Header>
auto alignmentBackend(const Header& header) -> AlignmentBackend {
    std::optional<AlignmentBackend> result;
    for (const auto& comment : header.comments) {
        if (!comment.starts_with(provenancePrefix)) {
            continue;
        }
        AlignmentBackend current;

        if (comment == "RNAnue alignment backend=star;adapter=1;version=2.7.11b") {
            current = AlignmentBackend::Star;
        } else if (comment == "RNAnue alignment backend=segemehl;adapter=1") {
            current = AlignmentBackend::Segemehl;
        } else {
            throw std::runtime_error("Unsupported RNAnue alignment provenance: " + comment);
        }

        if (result && *result != current) {
            throw std::runtime_error("Conflicting alignment backend provenance");
        }

        result = current;
    }
    for (const auto& program : header.program_infos) {
        if (program.name == "STAR" || program.id == "STAR") {
            if (!result) {
                throw std::runtime_error(
                    "STAR BAM lacks supported RNAnue alignment provenance; run RNAnue align first");
            }

            if (*result != AlignmentBackend::Star ||
                (!program.version.empty() && program.version != "2.7.11b")) {
                throw std::runtime_error(
                    "Conflicting STAR program and RNAnue alignment provenance");
            }
        }
    }
    return result.value_or(AlignmentBackend::Segemehl);
}

template <typename Header>
void stampAlignment(Header& header, AlignmentBackend backend, const std::string& command = {}) {
    const std::string marker = backend == AlignmentBackend::Star
                                   ? "RNAnue alignment backend=star;adapter=1;version=2.7.11b"
                                   : "RNAnue alignment backend=segemehl;adapter=1";
    std::erase_if(header.comments,
                  [](const auto& value) { return value.starts_with(provenancePrefix); });
    header.comments.push_back(marker);
    if (!command.empty()) {
        std::string id = "RNAnue.align";
        while (std::ranges::any_of(header.program_infos, [&](const auto& p) { return p.id == id; }))
            id += ".1";
        header.program_infos.push_back({.id = id,
                                        .name = "RNAnue",
                                        .command_line_call = command,
                                        .previous = {},
                                        .description = marker,
                                        .version = {},
                                        .user_tags = {}});
    }
}

void stampAlignmentFile(const std::filesystem::path& path, AlignmentBackend backend,
                        const std::string& command = {});
}  // namespace utility
