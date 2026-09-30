#pragma once

#include <filesystem>
#include <optional>
#include <vector>
#include "AlignParameters.hpp"
#include "SamReference.hpp"

namespace pipelines::align {
class StarAligner {
   public:
    explicit StarAligner(AlignParameters params)
        : parameters(std::move(params)), reference({}, {}) {}
    void buildIndex();
    void alignSingleReads(const std::filesystem::path& input,
                          const std::filesystem::path& output) const;
    void alignPairedReads(const std::filesystem::path& forward,
                          const std::filesystem::path& reverse,
                          const std::filesystem::path& output) const;
    [[nodiscard]] auto alignmentArguments() const -> std::vector<std::string>;

   private:
    AlignParameters parameters;
    std::filesystem::path indexPath;
    dataTypes::SamReference reference;
    void align(const std::filesystem::path& forward,
               const std::optional<std::filesystem::path>& reverse,
               const std::filesystem::path& output) const;
};
}  // namespace pipelines::align
