#pragma once

#include <cstddef>
#include <utility>

using NucleotidePairPositions = std::pair<size_t, size_t>;
using NucleotidePositionsWindow = std::pair<NucleotidePairPositions, NucleotidePairPositions>;
