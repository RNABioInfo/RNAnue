#pragma once

// Standard
#include <set>
#include <string>
#include <variant>
#include <vector>

struct DeduplicationOutputSingle {
    std::set<std::string> validRecordIDs;
};

struct DeduplicationOutputPaired {
    std::set<std::string> validRecordIDs;
};

using DeduplicationOutput = std::variant<DeduplicationOutputSingle, DeduplicationOutputPaired>;
