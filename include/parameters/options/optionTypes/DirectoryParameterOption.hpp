#pragma once

// Standard
#include <filesystem>
#include <optional>
#include <string_view>
#include <type_traits>
#include <utility>

// Boost
#include <boost/program_options/variables_map.hpp>

// Internal
#include "ParameterNames.hpp"
#include "ParameterOption.hpp"
#include "ParameterValidator.hpp"

/**
 * @brief Specialization for parameters representing directories.
 */
template <bool isOptional = false>
class DirectoryParameterOption : public ParameterOption<std::filesystem::path, isOptional> {
   public:
    using valueType =
        std::conditional_t<isOptional, std::optional<std::filesystem::path>, std::filesystem::path>;

    /**
     * @brief Constructs a new DirectoryParameterOption object.
     *
     * @param names The names for the directory parameter.
     * @param description A description of the directory parameter.
     */
    consteval explicit DirectoryParameterOption(ParameterOptionNames&& names,
                                                std::string_view description)
        : ParameterOption<std::filesystem::path, isOptional>(std::move(names), description) {}

    /**
     * @brief Extracts and validates the directory path from the variables_map.
     *
     * @param variables The Boost variables_map.
     * @return std::filesystem::path The validated directory path.
     */
    [[nodiscard]] auto extractValue(const po::variables_map& variables) const
        -> valueType override {
        if constexpr (isOptional) {
            return ParameterValidator::validateOptionalDirectory(variables, this->getLongName());
        } else {
            return ParameterValidator::validateDirectory(variables, this->getLongName());
        }
    }
};
