#pragma once

// Standard
#include <istream>
#include <sstream>
#include <string>
#include <unordered_map>

namespace pipelines::analyze {

using TranscriptContributionsByID = std::unordered_map<std::string, float>;

inline void parseTranscriptContributionStream(std::istream &transcriptCountsIn,
                                              TranscriptContributionsByID &transcriptCounts) {
    std::string line;
    while (std::getline(transcriptCountsIn, line)) {
        std::string transcriptID;
        std::string contributionToken;
        std::istringstream lineStream(line);

        if (!std::getline(lineStream, transcriptID, '\t') ||
            !std::getline(lineStream, contributionToken, '\t') || transcriptID.empty()) {
            continue;
        }

        transcriptCounts[transcriptID] += std::stof(contributionToken);
    }
}

}  // namespace pipelines::analyze
