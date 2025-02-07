#pragma once

// Standard
#include <cstddef>
#include <cstdlib>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

// Boost
#include <boost/program_options/errors.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>

// Internal
#include "ConstexprStringFrom.hpp"
#include "ParameterNames.hpp"
#include "ParameterValidator.hpp"

namespace po = boost::program_options;

using namespace std::string_view_literals;

/**
 * @brief Base class representing a command-line parameter option.
 *
 * @tparam T The type of the parameter's value.
 */
template <typename T, bool isOptional = false>
class ParameterOption {
   public:
    /**
     * @brief Constructs a new ParameterOption object.
     *
     * @param names The names (short, long, and options name) for the parameter.
     * @param description A brief description of the parameter.
     */
    consteval explicit ParameterOption(ParameterOptionNames&& names, std::string_view description)
        : names(std::move(names)), description(description) {}

    using valueType = std::conditional_t<isOptional, std::optional<T>, T>;
    using underlyingType = T;

    /**
     * @brief Retrieves the short name of the parameter.
     *
     * @return char The parameter's short name.
     */
    [[nodiscard]] constexpr auto getShortName() const noexcept -> std::optional<char> {
        return names.shortName;
    }

    /**
     * @brief Retrieves the long name of the parameter.
     *
     * @return const std::string& The parameter's long name.
     */
    [[nodiscard]] constexpr auto getLongName() const noexcept -> const std::string& {
        return names.longName;
    }

    /**
     * @brief Retrieves the description of the parameter.
     *
     * @return std::string A string containing the parameter's description.
     */
    [[nodiscard]] virtual auto getDescription() const -> std::string {
        return std::string{description};
    }

    /**
     * @brief Adds the parameter option to a Boost options_description.
     *
     * For boolean parameters, a bool_switch is used.
     *
     * @param optionsDescription The Boost options_description instance.
     */
    virtual auto addOptionTo(po::options_description& optionsDescription) const noexcept -> void {
        if constexpr (std::is_same_v<T, bool>) {
            optionsDescription.add_options()(names.optionsName().data(), po::bool_switch(),
                                             getDescription().data());
        } else {
            optionsDescription.add_options()(names.optionsName().data(),
                                             po::value<underlyingType>(), getDescription().data());
        }
    }

    /**
     * @brief Extracts and validates the parameter's value from a Boost variables_map.
     *
     * @param variables The Boost variables_map containing parsed options.
     * @return T The validated parameter value.
     */
    [[nodiscard]] virtual auto extractValue(const po::variables_map& variables) const -> valueType {
        if constexpr (isOptional) {
            return ParameterValidator::validateOptional<T>(variables, names.longName);
        } else {
            return ParameterValidator::validate<T>(variables, names.longName);
        }
    }

   protected:
    ParameterOptionNames names;
    std::string_view description;
};
