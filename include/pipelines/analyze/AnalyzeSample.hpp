#pragma once

// Standard
#include <filesystem>
#include <string>

namespace pipelines::analyze {
namespace fs = std::filesystem;

struct AnalyzeInput {
    std::string sampleName;

    fs::path splitAlignmentsPath;
    fs::path multisplitAlignmentsPath;
    fs::path unassignedContiguousAlignmentsPath;

    fs::path contiguousAlignmentsTranscriptCountsPath;

    fs::path sampleFragmentCountsPath;
};

struct AnalyzeOutput {
    fs::path interactionsPath;
    fs::path interactionsReadIDsPath;
    fs::path interactionsTranscriptCountsPath;

    fs::path interactionsBEDPath;
    fs::path interactionsBEDARCPath;
    fs::path interactionsArmCoverageBedGraphPath;

    fs::path supplementaryFeaturesPath;
};

struct AnalyzeSample {
    AnalyzeInput input;
    AnalyzeOutput output;
};

}  // namespace pipelines::analyze
