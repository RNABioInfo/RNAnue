#pragma once

// Standard
#include <cstddef>

// Boost
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <cstdlib>
#include <filesystem>
#include <optional>
#include <string>

// Boost
#include "LogLevel.hpp"
#include "boost/program_options/errors.hpp"

// Internal
#include "Logger.hpp"
#include "ValueBounds.hpp"
#include "concepts/ArithmeticConcept.hpp"

namespace po = boost::program_options;

struct ParameterValidator {
    template <typename T>
    static auto validateOptional(const po::variables_map& variables, const std::string& optionName)
        -> std::optional<T> {
        if (!variables.contains(optionName)) {
            return std::nullopt;
        }

        return validate<T>(variables, optionName);
    }

    template <typename T>
    static auto validate(const po::variables_map& variables, const std::string& optionName) -> T {
        if (!variables.contains(optionName)) {
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
                "Parsed options do not contain option named: ", optionName);
        }

        T value;

        try {
            value = variables.at(optionName).as<T>();
        } catch (const po::required_option& e) {
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(optionName,
                                                                " is a required parameter.");
            exit(EXIT_FAILURE);
        } catch (const po::invalid_option_value& e) {
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>("Invalid value for ", optionName,
                                                                ". Must be an ", typeid(T).name());
            exit(EXIT_FAILURE);
        } catch (const po::error& e) {
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
                "Unknown error occurred while parsing ", optionName, ". ", std::string(e.what()));
            exit(EXIT_FAILURE);
        }

        return value;
    }

    template <typename T>
        requires arithmetic<T>
    static auto validateArithmetic(const po::variables_map& params, const std::string& paramName,
                                   const ValueBounds<T> valueBounds) -> T {
        const auto value = validate<T>(params, paramName);

        if (!valueBounds.isValid(value)) {
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
                paramName + " must be an integer between " +
                std::to_string(valueBounds.lowerBound) + " and " +
                std::to_string(valueBounds.upperBound));
        }

        return value;
    }

    static auto validateFilePath(const po::variables_map& params, const std::string& paramName)
        -> std::filesystem::path {
        auto filePath = validate<std::filesystem::path>(params, paramName);

        if (!std::filesystem::exists(filePath) || std::filesystem::is_directory(filePath)) {
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
                "Check parameter '", paramName, "': ", filePath, " is not a valid file path.");
        }

        return filePath;
    }

    static auto validateDirectory(const po::variables_map& params, const std::string& paramName)
        -> std::filesystem::path {
        auto dirPath = validate<std::filesystem::path>(params, paramName);

        if (!std::filesystem::is_directory(dirPath)) {
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
                "Check parameter '", paramName, "': ", dirPath, " is not a valid directory.");
        }

        return dirPath;
    }

    static auto validateOptionalDirectory(const po::variables_map& params,
                                          const std::string& paramName)
        -> std::optional<std::filesystem::path> {
        if (!params.contains(paramName)) {
            return std::nullopt;
        }

        return validateDirectory(params, paramName);
    }
};
