#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "AlignParameters.hpp"
#include "SamReference.hpp"

namespace pipelines::align {
namespace fs = std::filesystem;
class SegemehlAligner {
   public:
    explicit SegemehlAligner(AlignParameters params) : parameters(std::move(params)) {}
    void buildIndex();
    void alignSingleReads(const fs::path& input, const fs::path& output) const;
    void alignPairedReads(const fs::path& forward, const fs::path& reverse,
                          const fs::path& output) const;

   private:
    AlignParameters parameters;
    fs::path indexPath;
    [[nodiscard]] auto findIndex(const fs::path& reference) const -> std::optional<fs::path>;
    [[nodiscard]] auto threadsAdaptedToEntries(const fs::path& input) const -> size_t;
    [[nodiscard]] auto getGeneralAlignmentArgs(size_t threads) const -> std::vector<std::string>;
    [[nodiscard]] auto referenceFromGenome() const -> dataTypes::SamReference;
    void writeEmptyAlignments(const fs::path& output, const fs::path& input) const;
    static void runSegemehlAlignment(std::vector<std::string> args, const fs::path& output,
                                     const std::string& error);
    static auto convertToCStrings(std::vector<std::string>& args) -> std::vector<char*>;
};
}  // namespace pipelines::align
