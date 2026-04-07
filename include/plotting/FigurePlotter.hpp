#pragma once

// Standard
#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <concepts>
#include <cstddef>
#include <filesystem>
#include <format>
#include <optional>
#include <random>
#include <string>
#include <string_view>
#include <vector>

// matplot
#include <matplot/matplot.h>
#include <matplot/util/handle_types.h>

// Internal
#include "Logger.hpp"
#include "NumericConcept.hpp"
#include "PlottingDefaults.hpp"
#include "Subsampling.hpp"
#include "matplot/core/figure_registry.h"
#include "matplot/freestanding/axes_functions.h"
#include "matplot/util/keywords.h"

namespace plotting {

namespace fs = std::filesystem;

struct FigureConfig {
    std::string title;
    size_t widthPixel;
    size_t heightPixel;

    size_t maxDataPointsPerScatterPlot;
    size_t maxDataPointsPerHistogram;

    std::string_view fileFormat;
    fs::path outDirPath;

    static auto makeDefault(const std::string& title, const fs::path& outDirPath) -> FigureConfig {
        return FigureConfig{.title = title,
                            .widthPixel = defaults::widthPixel,
                            .heightPixel = defaults::heightPixel,
                            .maxDataPointsPerScatterPlot = defaults::maxDataPointsPerScatterPlot,
                            .maxDataPointsPerHistogram = defaults::maxDataPointsPerHistogram,
                            .fileFormat = defaults::fileFormat,
                            .outDirPath = outDirPath};
    }
};

struct PlotConfig {
    std::string title;

    std::string xlabel;
    std::string ylabel;
};

template <NumericType T>
struct StackedBarData {
    std::vector<std::vector<T>> data;
    std::string_view legendTitle;
    std::vector<std::string> legendLabels;
    std::optional<std::vector<std::string>> groupLabels;
};

template <NumericType T>
struct HistogramData {
    std::vector<T> data;

    std::string datapointLabel;

    std::optional<std::string> cutoffLabel;
    std::optional<double> cutoffValue;
};

template <NumericType T>
struct ScatterData {
    std::vector<T> x_vals;
    std::vector<T> y_vals;

    std::string datapointLabel;
};

class FigurePlotter {
   public:
    FigurePlotter(const FigureConfig& config)
        : config(config), outFilePath(makeFilePath(config)), figure(setupFigure(config)) {}

    void save() {
        Logger::log("Saving figure to: ", outFilePath);
        figure->save(outFilePath);
    }

    template <NumericType T>
    void addStackedBar(const PlotConfig& config, const StackedBarData<T>& data) {
        auto axes = addAxes();
        auto groupLabels =
            data.groupLabels.value_or(std::vector<std::string>(data.data.size(), ""));
        auto plot = axes->barstacked(data.data);

        addLabels(axes, config);
        axes->x_axis().ticklabels(groupLabels);
        addLegend(axes, data.legendTitle, data.legendLabels);
    }

    template <NumericType T>
    void addHistogram(const PlotConfig& config,
                      const std::same_as<HistogramData<T>> auto&... datas) {
        static constexpr double face_alpha = 0.7;
        static_assert(sizeof...(datas) > 0,
                      "addHistogram() would like to have at least one data object...");

        auto axes = addAxes();

        std::vector<std::string> datapointLabels;
        datapointLabels.reserve(sizeof...(datas));

        for (const auto& data : {datas...}) {
            const auto values =
                subsample(data.data, this->config.maxDataPointsPerHistogram, randomDevice);
            auto plot = axes->hist(values, matplot::histogram::binning_algorithm::fd);
            plot->face_alpha(face_alpha);
            datapointLabels.push_back(data.datapointLabel);
            matplot::hold(matplot::on);
            if (data.cutoffValue) {
                addVerticalLine(axes, plot, data.cutoffLabel, *data.cutoffValue);
            }
        }

        matplot::hold(matplot::off);
        addLabels(axes, config);

        addLegend(axes, defaultLegendTitle, datapointLabels);
    }

    template <NumericType T>
    void addScatter(const PlotConfig& config, const ScatterData<T>& data) {
        const auto [x_vals, y_vals] =
            subsample(data.x_vals, data.y_vals, this->config.maxDataPointsPerScatterPlot);
        auto axes = addAxes();

        constexpr double dotSize = 3;
        auto plot =
            axes->scatter(matplot::to_vector_1d(x_vals), matplot::to_vector_1d(y_vals), dotSize);

        addLabels(axes, config);
        addLegend(axes, defaultLegendTitle, {data.datapointLabel});
    }

   private:
    static constexpr std::string_view defaultLegendTitle = "Datapoints";
    std::mt19937 randomDevice = std::mt19937{std::random_device{}()};

    FigureConfig config;
    fs::path outFilePath;

    matplot::figure_handle figure;
    std::vector<matplot::axes_handle> axis;

    [[nodiscard]] static auto setupFigure(const FigureConfig& config) -> matplot::figure_handle {
        auto figure = matplot::figure(true);

        figure->backend()->run_command("unset warnings");
        figure->size(config.widthPixel, config.heightPixel);
        figure->title(config.title);

        return figure;
    }

    [[nodiscard]] static auto makeFilePath(const FigureConfig& config) -> fs::path {
        std::string filename = config.title;

        // Make all lowercase
        std::ranges::transform(filename, filename.begin(),
                               [](unsigned char character) { return std::tolower(character); });

        static constexpr std::array<char, 10> toRemoveCharacters{'/', '.', '-', '?', '%',
                                                                 '!', '<', '>', '|', ' '};

        for (const char character : toRemoveCharacters) {
            std::ranges::replace(filename, character, '_');
        }

        filename += ".";
        filename += config.fileFormat;

        return config.outDirPath / filename;
    }

    [[nodiscard]] auto addAxes() -> matplot::axes_handle {
        axis.emplace_back(figure->nexttile());

        return axis.back();
    }

    static void addLabels(const matplot::axes_handle& axes, const PlotConfig& config) {
        axes->title(config.title);
        axes->xlabel(config.xlabel);
        axes->ylabel(config.ylabel);
    }

    static void addLegend(const matplot::axes_handle& axes, const std::string_view legendTitle,
                          const std::vector<std::string>& legendValues) {
        axes->legend(legendValues);
        axes->legend()->title(legendTitle);
        axes->legend()->box(false);
    }

    static void addVerticalLine(const matplot::axes_handle& axes,
                                const matplot::axes_object_handle& plot,
                                const std::optional<std::string>& cutoffLabel, double cutoffValue) {
        const double min = plot->ymin();
        const double max = plot->ymax();

        auto line = axes->line(cutoffValue, min, cutoffValue, max);
        line->color("red");
        line->line_style("--");
        line->impulse(true);

        if (cutoffLabel) {
            const double offset = (plot->xmax() - plot->xmin()) * 0.01;
            const auto label = *cutoffLabel + std::format("{:.2f}", cutoffValue);

            axes->text(cutoffValue + offset, max - offset, label)->font("Times New Roman");
        }
    }
};

}  // namespace plotting
