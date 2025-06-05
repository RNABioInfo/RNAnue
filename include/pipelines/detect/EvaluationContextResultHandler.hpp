#pragma once

// Standard
#include <cstddef>
#include <cstdlib>
#include <deque>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

// seqan3
#include <seqan3/io/sam_file/format_bam.hpp>
#include <seqan3/io/sam_file/format_sam.hpp>
#include <seqan3/io/sam_file/output.hpp>

// Internal
#include "AnnotationStep.hpp"
#include "EvaluationContext.hpp"
#include "HitGroup.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "OneOf.hpp"
#include "ReadGroupPostScoringStep.hpp"
#include "SamRecord.hpp"
#include "SamReference.hpp"
#include "TempOutputDirs.hpp"
#include "Utility.hpp"

namespace pipelines::detect {

namespace fs = std::filesystem;

struct EvaluationContextResultHandler {
   public:
    using TranscriptCounts = std::unordered_map<std::string, double>;
    using OutputT =
        seqan3::sam_file_output<SamFieldIDs,
                                seqan3::type_list<seqan3::format_sam, seqan3::format_bam>,
                                std::deque<std::string>>;

    EvaluationContextResultHandler(const TempOutputDirs& outDirs, SamReference& reference,
                                   size_t bufferSize) noexcept
        : bufferSize(bufferSize),
          singletonUnassignedOutput{makePath(outDirs.outputTmpUnassignedSingletonDir),
                                    reference.referenceIDs, reference.referenceLengths,
                                    SamFieldIDs{}},
          chimericOutput{makePath(outDirs.outputTmpSplitsDir), reference.referenceIDs,
                         reference.referenceLengths, SamFieldIDs{}},
          multimericOutput{makePath(outDirs.outputTmpMultisplitsDir), reference.referenceIDs,
                           reference.referenceLengths, SamFieldIDs{}} {
        singletonRecordBuffer.reserve(bufferSize);
        chimericRecordBuffer.reserve(bufferSize);
    }

    [[nodiscard]] auto getSingletonTranscriptCounts() const -> const TranscriptCounts& {
        return singletonTranscriptCounts;
    }

    template <typename... Results>
        requires(one_of<AnnotationStepResult, Results...> &&
                 one_of<ReadGroupPostScoringResult, Results...>)
    void operator()(const EvaluationContext<SingletonHitGroup, Results...>& context) noexcept {
        const auto& annotation = context.template get<AnnotationStepResult>();
        if (annotation.feature) {
            const auto& score = context.template get<ReadGroupPostScoringResult>();
            singletonTranscriptCounts[annotation.feature->getAnnotationID()] +=
                score.contributionScore;
        } else {
            const auto& records = context.getRecords();
            singletonRecordBuffer.insert(singletonRecordBuffer.end(), records.begin(),
                                         records.end());
        }

        checkWriteBuffers();
    }

    template <typename... Results>
    void operator()(const EvaluationContext<ChimericHitGroup, Results...>& context) noexcept {
        const auto& records = context.getRecords();
        chimericRecordBuffer.insert(chimericRecordBuffer.end(), records.begin(), records.end());

        checkWriteBuffers();
    }

    void operator()(auto const& /* unused */) noexcept {
        Logger::log<SourceLocation{}, LogLevel::ERROR>("Not all context cases have been handled!");
    }

    void save() noexcept {
        for (SamRecord& record : singletonRecordBuffer) {
            singletonUnassignedOutput.push_back(record);
        }

        singletonRecordBuffer.clear();

        for (SamRecord& record : chimericRecordBuffer) {
            chimericOutput.push_back(record);
        }

        chimericRecordBuffer.clear();
    }

   private:
    // Generates output file path for each category
    static auto makePath(const fs::path& dir) noexcept -> fs::path {
        return dir / (helper::getUUID() + ".bam");
    }

    void checkWriteBuffers() noexcept {
        if (singletonRecordBuffer.size() + chimericRecordBuffer.size() > bufferSize) {
            save();
        }
    }

    size_t bufferSize;
    OutputT singletonUnassignedOutput;
    OutputT chimericOutput;
    OutputT multimericOutput;
    TranscriptCounts singletonTranscriptCounts;
    std::vector<SamRecord> singletonRecordBuffer;
    std::vector<SamRecord> chimericRecordBuffer;
};

}  // namespace pipelines::detect
