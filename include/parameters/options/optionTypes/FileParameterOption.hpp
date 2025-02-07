#pragma once

// Standard
#include <filesystem>
#include <string_view>
#include <utility>

// Boost
#include <boost/program_options/variables_map.hpp>

// Internal
#include "ParameterNames.hpp"
#include "ParameterOption.hpp"
#include "ParameterValidator.hpp"

/**
 * @brief Specialization for parameters representing file paths.
 */
class FileParameterOption : public ParameterOption<std::filesystem::path> {
   public:
    /**
     * @brief Constructs a new FileParameterOption object.
     *
     * @param names The names for the file parameter.
     * @param description A description of the file parameter.
     */
    consteval explicit FileParameterOption(ParameterOptionNames&& names,
                                           std::string_view description)
        : ParameterOption<std::filesystem::path>(std::move(names), description) {}

    /**
     * @brief Extracts and validates the file path from the variables_map.
     *
     * @param variables The Boost variables_map.
     * @return std::filesystem::path The validated file path.
     */
    [[nodiscard]] auto extractValue(const po::variables_map& variables) const
        -> std::filesystem::path override {
        return ParameterValidator::validateFilePath(variables, names.longName);
    }
};
