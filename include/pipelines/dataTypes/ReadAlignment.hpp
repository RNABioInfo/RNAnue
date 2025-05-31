#pragma once

// Standard
#include <variant>

// Internal

namespace dataTypes {
struct SingletonAlignment {};

struct SplitAlignment {};

using ReadAlignmentVariant = std::variant<SingletonAlignment, SplitAlignment>;

}  // namespace dataTypes
