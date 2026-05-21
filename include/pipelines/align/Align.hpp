#pragma once

// Standard
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <optional>
#include <string>
#include <utility>
#include <vector>

// Boost
#include <boost/program_options.hpp>

// segemehl
extern "C" {
#include <segemehl.h>
}  // segemehl

// Class
#include "AlignData.hpp"
#include "AlignParameters.hpp"
#include "AlignSample.hpp"
#include "SamReference.hpp"

namespace pipelines::align {

class Align {
   public:
    explicit Align(AlignParameters params) : parameters(std::move(params)) {};

    void process(const AlignData &data);

   private:
    AlignParameters parameters;
    fs::path indexPath;

    [[nodiscard]] auto threadsAdaptedToEntries(const fs::path &inputPath) const -> size_t;

    void preprocessReferences();

    void processSample(const AlignSampleType &sample);

    void processSingleEnd(const AlignSampleSingle &sample);
    void processMergedPairedEnd(const AlignSampleMergedPaired &sample);

    [[nodiscard]] auto findIndex(const fs::path &referenceGenomePath) const
        -> std::optional<fs::path>;

    void buildIndex();

    [[nodiscard]] auto referenceFromGenome() const -> dataTypes::SamReference;
    void writeEmptyAlignments(const fs::path &alignmentsOutPath,
                              const fs::path &emptyInputPath) const;

    [[nodiscard]] auto getGeneralAlignmentArgs(size_t threadCount) const
        -> std::vector<std::string>;
    void runSegemehlAlignment(std::vector<std::string> args, const fs::path &outputPath,
                              const std::string &errorMessage) const;
    void alignReads(const std::string &query, const std::string &mate,
                    const std::string &matched) const;
    void alignSingleReads(const fs::path &queryFastqInPath,
                          const fs::path &alignmentsFastqOutPath) const;
    void alignPairedReads(const fs::path &queryForwardFastqInPath,
                          const fs::path &queryReverseFastqInPath,
                          const fs::path &alignmentsFastqOutPath) const;

    void sortAlignmentsByQueryName(const fs::path &alignmentsPath,
                                   const fs::path &sortedAlignmentsPath) const;

    static auto convertToCStrings(std::vector<std::string> &args) -> std::vector<char *>;
};

}  // namespace pipelines::align
