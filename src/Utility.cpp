#include "Utility.hpp"

// Standard
#include <execinfo.h>
#include <unistd.h>

#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

// seqan3
#include <seqan3/io/record.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

// boost
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>

// Internal
#include "Config.hpp"
#include "Logger.hpp"

namespace helper {

void crashHandler(int sig) {
    constexpr size_t MAX_FRAMES = 10;
    std::array<void*, MAX_FRAMES> array{};
    int size = backtrace(array.data(), MAX_FRAMES);

    // print out all the frames to stderr
    std::cerr << "Error: signal " << sig << ":" << '\n';
    backtrace_symbols_fd(array.data(), size, STDERR_FILENO);
    exit(1);
}

void createTmpDir(const fs::path& path) {
    Logger::log("Create temporary directory " + path.string());
    deleteDir(path);
    fs::create_directory(path);
}

// delete folder
void deleteDir(const fs::path& path) {
    if (fs::exists(path)) {
        fs::remove_all(path);
    }
}

auto getDirIfExists(const fs::path& path) -> std::optional<fs::path> {
    if (fs::is_directory(path)) {
        return path;
    }

    return std::nullopt;
}

auto getUUID() -> std::string {
    boost::uuids::random_generator uuidGenerator;
    return boost::uuids::to_string(uuidGenerator());
}

void mergeSamFiles(const std::vector<fs::path>& inputPaths, const fs::path& outputPath) {
    if (inputPaths.empty()) {
        Logger::log<LogLevel::WARNING>("No input files to merge");
        return;
    }

    seqan3::sam_file_output outputFile{outputPath};

    for (const auto& inputPath : inputPaths) {
        if (!fs::exists(inputPath) || fs::file_size(inputPath) == 0) {
            continue;
        }

        seqan3::sam_file_input inputFile{inputPath};

        if (inputFile.begin() == inputFile.end()) {
            continue;
        }

        inputFile | outputFile;
    }
}

void mergeFastqFiles(const std::vector<fs::path>& inputPaths, const fs::path& outputPath) {
    if (inputPaths.empty()) {
        Logger::log<LogLevel::WARNING>("No input files to merge");
        return;
    }

    seqan3::sequence_file_output outputFile{outputPath};

    constexpr size_t MIN_ZIPPED_FILE_SIZE = 24;

    for (const auto& inputPath : inputPaths) {
        if (!fs::exists(inputPath) || inputPath.extension() == ".gz"
                ? fs::file_size(inputPath) <= MIN_ZIPPED_FILE_SIZE
                : fs::file_size(inputPath) == 0) {
            Logger::log<LogLevel::DEBUG>(inputPath);
            continue;
        }
        seqan3::sequence_file_input inputFile{inputPath};
        inputFile | outputFile;
    }
}

auto getFilePathsInDir(const fs::path& dir) -> std::vector<fs::path> {
    std::vector<fs::path> filePaths;
    for (const auto& entry : fs::directory_iterator(dir)) {
        if (fs::is_regular_file(entry)) {
            filePaths.push_back(entry.path());
        }
    }

    return filePaths;
}

void concatAndDeleteFilesInTmpDir(const fs::path& tmpDir, const fs::path& outPath) {
    std::ofstream outStream(outPath, std::ios::binary);

    for (const auto& entry : fs::directory_iterator(tmpDir)) {
        std::ifstream inStream(entry.path(), std::ios::binary);
        outStream << inStream.rdbuf();
        inStream.close();
    }

    outStream.close();

    deleteDir(tmpDir);
}

auto isContained(const int32_t value, const int32_t comparisonValue, const int32_t tolerance)
    -> bool {
    return value >= comparisonValue - tolerance && value <= comparisonValue + tolerance;
}

auto countUniqueSamEntries(fs::path& path) -> std::size_t {
    std::ifstream infile(path);
    if (!infile.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    std::unordered_set<std::string> uniqueIDs;
    std::string line;
    while (std::getline(infile, line)) {
        if (!line.empty() && line[0] != '@') {
            std::string identifier = line.substr(0, line.find('\t'));
            uniqueIDs.insert(identifier);
        }
    }

    return uniqueIDs.size();
}

auto countSamEntriesSeqAn(fs::path& path) -> std::size_t {
    seqan3::sam_file_input fin{path.string(), seqan3::fields<>{}};
    return std::ranges::distance(fin.begin(), fin.end());
}

auto countSamEntries(fs::path& path) -> size_t {
    std::ifstream infile(path);
    if (!infile.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    int count = 0;
    std::string line;
    while (std::getline(infile, line)) {
        if (!line.empty() && line[0] != '@') {
            ++count;
        }
    }

    return count;
}

auto generateRandomHexColor() -> std::string {
    std::random_device randomDevice;
    std::mt19937 gen(randomDevice());
    std::uniform_int_distribution<> distr(0, 255);  // NOLINT

    std::stringstream stringStream;
    stringStream << "#";
    for (int i = 0; i < 3; ++i) {
        stringStream << std::setfill('0') << std::setw(2) << std::hex << distr(gen);
    }

    return stringStream.str();  // Return the hex color code as a string
}

auto getTime() -> std::string {
    const auto now = std::chrono::system_clock::now();
    const std::time_t current_time = std::chrono::system_clock::to_time_t(now);

    std::ostringstream time_stream;
    time_stream << std::put_time(std::localtime(&current_time), "[%Y-%m-%d %H:%M:%S]") << " ";

    return time_stream.str();
}

void Timer::stop() {
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";
}

}  // namespace helper
