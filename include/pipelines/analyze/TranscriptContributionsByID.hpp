#pragma once

// Standard
#include <string>
#include <unordered_map>

namespace pipelines::analyze {

using TranscriptContributionsByID = std::unordered_map<std::string, float>;

}
