#pragma once

// Standard
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

// Boost
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>

// Internal
#include "ConstexprStringFrom.hpp"
#include "ParameterNames.hpp"
#include "ParameterOption.hpp"

namespace po = boost::program_options;

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
    [[nodiscard]] auto getDescription() const -> std::string override {
        // TODO: Implement default description  string_from<T, defaultValue>::value

        return std::string{this->description} + " (default: " + std::to_string(defaultValue) + ")";
    }

    [[nodiscard]] auto getInverseDescription() const -> std::string {
        if (!this->names.inverseDescription.empty()) {
            return std::string{this->names.inverseDescription};
        }

        return "set " + this->names.longName + " to false";
    }

    /**
     * @brief Adds the parameter option with its default value to a Boost options_description.
     *
     * @param optionsDescription The Boost options_description instance.
     */
    auto addOptionTo(po::options_description& optionsDescription) const noexcept -> void override {
        if constexpr (std::is_same_v<T, bool>) {
            optionsDescription.add_options()(this->names.optionsName().data(),
                                             po::value<bool>()
                                                 ->default_value(defaultValue)
                                                 ->implicit_value(true),
                                             getDescription().data());

            if (this->names.hasInverseName()) {
                optionsDescription.add_options()(
                    this->names.inverseLongName.data(),
                    po::bool_switch()->default_value(false),
                    getInverseDescription().data());
            }
        } else if constexpr (std::is_floating_point_v<T>) {
            optionsDescription.add_options()(
                this->names.optionsName().data(),
                po::value<T>()->default_value(defaultValue, std::format("{:.{}}", defaultValue, 3)),
                getDescription().data());
        } else {
            optionsDescription.add_options()(this->names.optionsName().data(),
                                             po::value<T>()->default_value(defaultValue),
                                             getDescription().data());
        }
    }

    [[nodiscard]] auto extractValue(const po::variables_map& variables) const -> T override {
        auto value = ParameterOption<T>::extractValue(variables);

        if constexpr (std::is_same_v<T, bool>) {
            if (this->names.hasInverseName() &&
                variables.contains(std::string{this->names.inverseLongName}) &&
                variables.at(std::string{this->names.inverseLongName}).template as<bool>()) {
                return false;
            }
        }

        return value;
    }
};
