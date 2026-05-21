#pragma once

// Standard
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>

// Internal
#include "AlignParameters.hpp"
#include "FeatureParser.hpp"
#include "FeatureWriter.hpp"
#include "FileType.hpp"
#include "GeneCopyMasker.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "ReferenceGenome.hpp"
#include "seqan3/io/sequence_file/output.hpp"

namespace pipelines::align {

namespace fs = std::filesystem;

class GenomePreprocessor {
   public:
    GenomePreprocessor(const AlignParameters& params) : params(params) {}

    struct Input {
        fs::path refGenomePath;
        fs::path annotationPath;
    };

    struct Output {
        fs::path preprocessedGenomePath;
        fs::path preprocessedAnnotationPath;
    };

    void process(const Input& input, const Output& output) const {
        ReferenceGenome referenceGenome = ReferenceGenome{input.refGenomePath};
        const auto parser = annotation::FeatureParser::defaultAllTranscript();

        GeneCopyMasker geneCopyMasker{
            params.geneCopyMaskerParams,
            parser.parseGroupedByHierarchy(input.annotationPath,
                                           referenceGenome.getReferenceIndexMapping()),
            std::move(referenceGenome)};
        auto result = geneCopyMasker.process();
        writeResults(output, result);
    };

   private:
    struct Parameters {
        explicit Parameters(const AlignParameters& params) : geneCopyMaskerParams(params) {}

        GeneCopyMasker::Parameters geneCopyMaskerParams;
    };

    const Parameters params;

    static auto writeResults(const Output& output, const GeneCopyMasker::Result& result) -> void {
        annotation::FeatureWriter::write(
            result.maskedClusters,
            result.maskedGenome.getReferenceIndexMapping().sortedReferenceIDs(),
            output.preprocessedAnnotationPath, annotation::FileType::GFF);

        std::ofstream outputGenome(output.preprocessedGenomePath);

        if (!outputGenome.is_open()) {
            Logger::log<SourceLocation{}, LogLevel::ERROR>("Failed to open output genome file");
        }

        seqan3::sequence_file_output genomeOutput{output.preprocessedGenomePath};

        for (int index = 0; index < static_cast<int>(result.maskedGenome.getSequenceCount());
             ++index) {
            const auto& sequence = result.maskedGenome.getSequence(index);
            const auto& referenceID = result.maskedGenome.getSequenceID(index);

            if (referenceID) {
                genomeOutput.emplace_back(sequence, std::string(referenceID.value()));
            }
        }
    };
};

}  // namespace pipelines::align
