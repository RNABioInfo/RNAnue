#pragma once

// Standard
#include <set>
#include <cstddef>
#include <variant>

struct DeduplicationOutputSingle {
    std::set<size_t> validRecordOrdinals;
};

struct DeduplicationOutputPaired {
    std::set<size_t> validRecordOrdinals;
};

using DeduplicationOutput = std::variant<DeduplicationOutputSingle, DeduplicationOutputPaired>;
