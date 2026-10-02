#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace plotting {

// Keep Freedman-Diaconis bin widths, but supply explicit edges to Matplot.
// Its automatic bin picker reserves a negative range for a single bin and
// its quartile indexing is out of bounds for two or three observations.
inline auto histogramBinEdges(std::vector<double> values) -> std::vector<double> {
    if (values.empty()) return {};
    if (std::ranges::any_of(values, [](double value) { return !std::isfinite(value); })) {
        throw std::invalid_argument("Histogram observations must be finite");
    }
    std::ranges::sort(values);
    const double minimum = values.front();
    const double maximum = values.back();
    if (minimum == maximum) {
        const double padding = std::max(0.5, std::abs(minimum) * 0.01);
        const double limit = std::numeric_limits<double>::max();
        return {std::max(-limit, minimum - padding), std::min(limit, maximum + padding)};
    }

    const double range = maximum - minimum;
    if (!std::isfinite(range)) return {minimum, maximum};
    const size_t lowerQuartile = values.size() / 4;
    const size_t upperQuartile = std::min(values.size() - 1, values.size() - lowerQuartile);
    const double interquartileRange = values[upperQuartile] - values[lowerQuartile];
    const double width = 2.0 * std::max(interquartileRange, range / 10.0) /
                         std::cbrt(static_cast<double>(values.size()));
    const size_t binCount = static_cast<size_t>(std::clamp(std::ceil(range / width), 1.0, 65536.0));

    std::vector<double> edges;
    edges.reserve(binCount + 1);
    edges.push_back(minimum);
    for (size_t bin = 1; bin < binCount; ++bin) {
        const double edge = minimum + range * (static_cast<double>(bin) / binCount);
        // Floating-point rounding can collapse adjacent edges in narrow ranges.
        if (edge > edges.back() && edge < maximum) edges.push_back(edge);
    }
    edges.push_back(maximum);
    return edges;
}

}  // namespace plotting
