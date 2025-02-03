#pragma once

// Standard
#include <filesystem>
#include <variant>

// Internal
#include "DeduplicationOutput.hpp"

namespace fs = std::filesystem;

struct DeduplicationBySequenceSingleConfig {
    using OutputType = DeduplicationOutputSingle;

    fs::path recordsPath;
};

struct DeduplicationBySequencePairedConfig {
    using OutputType = DeduplicationOutputPaired;

    fs::path recordsPathFwd;
    fs::path recordsPathRev;
};

using DeduplicationBySequenceInput =
    std::variant<DeduplicationBySequenceSingleConfig, DeduplicationBySequencePairedConfig>;
