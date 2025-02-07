#pragma once

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <string>
#include <string_view>
#include <utility>

#include "ParameterNames.hpp"
#include "ParameterOption.hpp"

/**
 * @brief Represents a command-line string parameter option that has a default value.
 *
 * This class extends ParameterOption by providing a default value for a string parameter.
 * It automatically appends the default value to the option's description.
 */
class DefaultedStringParameterOption : public ParameterOption<std::string> {
   public:
    /**
     * @brief Constructs a DefaultedStringParameterOption.
     *
     * @param names The parameter names associated with this option.
     * @param description The description of the option, wrapped in a tag type to avoid mix-ups.
     * @param defaultValue The default value for the option, wrapped in a tag type to avoid mix-ups.
     */
    consteval DefaultedStringParameterOption(ParameterOptionNames&& names,
                                             std::string_view description,  // NOLINT
                                             std::string_view defaultValue)
        : ParameterOption<std::string>(std::move(names), description), defaultValue{defaultValue} {}

    /**
     * @brief Retrieves the default value for this parameter option.
     *
     * @return A std::string_view representing the default value.
     */
    [[nodiscard]] constexpr auto getDefaultValue() const noexcept -> std::string_view {
        return defaultValue;
    }

    /**
     * @brief Gets the complete description of the option, including its default value.
     *
     * @return A std::string containing the description with the default value appended.
     */
    [[nodiscard]] auto getDescription() const -> std::string override {
        return std::string{this->description} + " (default: " + std::string{defaultValue} + ")";
    }

    /**
     * @brief Adds this option to a Boost Program Options description.
     *
     * This function registers the option with Boost's options_description,
     * specifying its name, value semantic, default value, and description.
     *
     * @param optionsDescription A reference to the Boost options_description object.
     */
    auto addOptionTo(po::options_description& optionsDescription) const noexcept -> void override {
        optionsDescription.add_options()(
            this->names.optionsName().data(),
            po::value<std::string>()->default_value(std::string{defaultValue}),
            getDescription().data());
    }

   private:
    std::string_view defaultValue;  ///< The default value for the option.
};
