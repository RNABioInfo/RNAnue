#pragma once

// Standard
#include <cstddef>
#include <string_view>

// Include
#include "FileFormat.hpp"

namespace plotting::defaults {

static constexpr size_t widthPixel = 1920;
static constexpr size_t heightPixel = 1080;
static constexpr size_t maxDatapointsPerPlot = 50000;
static constexpr std::string_view fileFormat = FileFormat::svg;

}  // namespace plotting::defaults
