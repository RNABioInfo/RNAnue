#pragma once

#include <filesystem>
#include <utility>
#include <variant>

#include "AlignData.hpp"
#include "AlignParameters.hpp"
#include "SegemehlAligner.hpp"
#include "StarAligner.hpp"

namespace pipelines::align {
class Align {
   public:
    explicit Align(AlignParameters params)
        : parameters(params), backend(SegemehlAligner{std::move(params)}) {}
    void process(const AlignData& data);

   private:
    AlignParameters parameters;
    std::variant<SegemehlAligner, StarAligner> backend;
    void preprocessReferences();
    void processSample(const AlignSampleType& sample);
    void processSingleEnd(const AlignSampleSingle& sample);
    void processMergedPairedEnd(const AlignSampleMergedPaired& sample);
    void alignSingleReads(const fs::path& input, const fs::path& output);
    void alignPairedReads(const fs::path& forward, const fs::path& reverse, const fs::path& output);
    [[nodiscard]] auto referenceFromGenome() const -> dataTypes::SamReference;
    void sortAlignmentsByQueryName(const fs::path& alignmentsPath,
                                   const fs::path& sortedAlignmentsPath) const;
};
}  // namespace pipelines::align
