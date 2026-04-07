#pragma once

// Standard
#include <cstddef>
#include <limits>
#include <string_view>

// Include
#include "FileFormat.hpp"

namespace plotting::defaults {

static constexpr size_t widthPixel = 1920;
static constexpr size_t heightPixel = 1080;
static constexpr size_t maxDataPointsPerScatterPlot = 50000;
static constexpr size_t maxDataPointsPerHistogram = std::numeric_limits<size_t>::max();
static constexpr std::string_view fileFormat = FileFormat::svg;

}  // namespace plotting::defaults
