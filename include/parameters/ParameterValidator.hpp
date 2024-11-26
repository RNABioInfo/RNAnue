#pragma once

// Standard
#include <concepts>
#include <cstddef>

// Boost
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <cstdlib>
#include <filesystem>
#include <string>

// Internal
#include "Logger.hpp"
#include "boost/program_options/errors.hpp"

namespace po = boost::program_options;

template <typename T>
concept arithmetic = std::integral<T> or std::floating_point<T>;

struct ParameterValidator {
    template <typename T>
        requires arithmetic<T>
    static auto validateArithmetic(const po::variables_map& params, const std::string& paramName,
                                   const T lowerLimit, const T upperLimit) -> T {
        Logger::log<LogLevel::DEBUG>(
            "Validating ", paramName, " parameter. Type: ", typeid(T).name(),
            ". Lower limit: ", lowerLimit, ". Upper limit: ", upperLimit, ".");

        T value;

        try {
            value = params[paramName].as<T>();
        } catch (const po::required_option& e) {
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(paramName,
                                                                " is a required parameter.");
            exit(EXIT_FAILURE);
        } catch (const po::invalid_option_value& e) {
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
                paramName, " must be an integer between ", std::to_string(lowerLimit), " and ",
                std::to_string(upperLimit));
            exit(EXIT_FAILURE);
        } catch (const po::error& e) {
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
                "Unknown error occurred while parsing ", paramName, ". ", std::string(e.what()));
            exit(EXIT_FAILURE);
        }

        if (value <= lowerLimit && value >= upperLimit) {
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
                paramName + " must be an integer between " + std::to_string(lowerLimit) + " and " +
                std::to_string(upperLimit));
        }

        return value;
    }

    static auto validateFilePath(const po::variables_map& params, const std::string& paramName)
        -> std::filesystem::path {
        const std::string filePathStr = params[paramName].as<std::string>();
        std::filesystem::path filePath = std::filesystem::path(filePathStr);

        if (!std::filesystem::exists(filePath) || std::filesystem::is_directory(filePath)) {
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
                "Check parameter '", paramName, "': ", filePath, " is not a valid file path.");
        }

        return filePath;
    }

    static auto validateDirectory(const po::variables_map& params, const std::string& paramName)
        -> std::filesystem::path {
        const std::string dirPathStr = params[paramName].as<std::string>();
        std::filesystem::path dirPath = std::filesystem::path(dirPathStr);

        if (!std::filesystem::is_directory(dirPath)) {
            Logger::log<IncludeSourceLocation, LogLevel::ERROR>(
                "Check parameter '", paramName, "': ", dirPath, " is not a valid directory.");
        }

        return dirPath;
    }
};
