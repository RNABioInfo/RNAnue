#pragma once

// Standard
#include <vector>

// Internal
#include "FastqRecord.hpp"

using namespace dataTypes;

struct DeduplicationOutputSingle {
    std::vector<FastqRecord> records;
};

struct DeduplicationOutputPaired {
    std::vector<std::pair<FastqRecord, FastqRecord>> recordPairs;
};

using DeduplicationOutput = std::variant<DeduplicationOutputSingle, DeduplicationOutputPaired>;
