#pragma once

// Standard
#include <string>
#include <string_view>
#include <utility>

// Boost
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>

// Internal
#include "ParameterNames.hpp"
#include "ParameterOption.hpp"

/**
 * @brief Parameter option that includes a compile-time default value.
 *
 * @tparam T The type of the parameter.
 * @tparam defaultValue The compile-time default value.
 */
template <typename T, T defaultValue>
class DefaultedParameterOption : public ParameterOption<T> {
   public:
    /**
     * @brief Constructs a new DefaultedParameterOption object.
     *
     * @param names The names for the parameter.
     * @param description A description of the parameter.
     */
    consteval explicit DefaultedParameterOption(ParameterOptionNames&& names,
                                                std::string_view description)
        : ParameterOption<T>(std::move(names), description) {}

    /**
     * @brief Retrieves the compile-time default value.
     *
     * @return constexpr T The default value.
     */
    [[nodiscard]] constexpr auto getDefaultValue() const noexcept -> T { return defaultValue; }

    /**
     * @brief Retrieves the description with the default value appended.
     *
     * @return std::string The description including the default value.
     */
    [[nodiscard]] constexpr auto getDescription() const -> std::string override {
        // TODO: Implement default description  string_from<T, defaultValue>::value
        return std::string{this->description} + " (default: " + ")";
    }

    /**
     * @brief Adds the parameter option with its default value to a Boost options_description.
     *
     * @param optionsDescription The Boost options_description instance.
     */
    auto addOptionTo(po::options_description& optionsDescription) const noexcept -> void override {
        if constexpr (std::is_same_v<T, bool>) {
            optionsDescription.add_options()(this->names.optionsName().data(),
                                             po::bool_switch()->default_value(defaultValue),
                                             getDescription().data());
        } else {
            optionsDescription.add_options()(this->names.optionsName().data(),
                                             po::value<T>()->default_value(defaultValue),
                                             getDescription().data());
        }
    }
};
