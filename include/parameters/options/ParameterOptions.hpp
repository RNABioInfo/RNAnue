#pragma once

// Boost
#include <boost/program_options/options_description.hpp>
#include <string>

namespace po = boost::program_options;

class ParameterOptions {
   public:
    ParameterOptions() = delete;
    ParameterOptions(const ParameterOptions &) = default;
    ParameterOptions(ParameterOptions &&) = delete;
    auto operator=(const ParameterOptions &) -> ParameterOptions & = default;
    auto operator=(ParameterOptions &&) -> ParameterOptions & = delete;
    ~ParameterOptions() = delete;

    template <typename T>
    static auto getOptions() -> po::options_description {
        po::options_description options{std::string(T::optionsDescription)};

        std::apply([&options](const auto &...option) { (option.addOptionTo(options), ...); },
                   T::allOptions);

        return options;
    };
};
