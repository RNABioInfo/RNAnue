#pragma once

// Standard
#include <string_view>
#include <utility>

// Boost
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>

// Internal
#include "DefaultedParameterOption.hpp"
#include "ParameterNames.hpp"
#include "ParameterValidator.hpp"
#include "SatisfiesArithmeticParameterConstraints.hpp"
#include "ValueBounds.hpp"

/**
 * @brief Parameter option for arithmetic types with default values and bounds.
 *
 * @tparam T The arithmetic type of the parameter.
 * @tparam defaultValue The compile-time default value.
 * @tparam valueBounds The bounds for the parameter's value.
 */
template <typename T, T defaultValue, ValueBounds<T> valueBounds>
    requires SatisfiesArithmeticParameterConstraints<T, defaultValue, valueBounds>
class ArithmeticParameterOption : public DefaultedParameterOption<T, defaultValue> {
   public:
    /**
     * @brief Constructs a new ArithmeticParameterOption object.
     *
     * @param names The names for the arithmetic parameter.
     * @param description A description of the arithmetic parameter.
     */
    consteval explicit ArithmeticParameterOption(ParameterOptionNames&& names,
                                                 std::string_view description)
        : DefaultedParameterOption<T, defaultValue>(std::move(names), description) {}

    /**
     * @brief Extracts and validates the arithmetic value from the variables_map.
     *
     * @param variables The Boost variables_map.
     * @return T The validated arithmetic value.
     */
    [[nodiscard]] auto extractValue(const po::variables_map& variables) const -> T override {
        return ParameterValidator::validateArithmetic(variables, this->getLongName(), valueBounds);
    }
};

// Example instantiation and static check.
static constexpr ArithmeticParameterOption<int, 10, {.lowerBound = 0, .upperBound = 10}> testParam{
    {.shortName = 't', .longName = "test"}, "This is a very long test description..."sv};

static_assert(testParam.getShortName() == 't', "Test parameter short name should be 't'");
