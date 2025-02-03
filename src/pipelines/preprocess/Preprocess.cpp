#include "Preprocess.hpp"

// Standard
#include <cassert>
#include <utility>
#include <variant>

// seqan3
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

// Internal
#include "Constants.hpp"
#include "Logger.hpp"
#include "PairedEndPreprocessor.hpp"
#include "PreprocessData.hpp"
#include "PreprocessParameters.hpp"
#include "PreprocessSample.hpp"
#include "SingleEndPreprocessor.hpp"
#include "VariantOverload.hpp"

namespace pipelines::preprocess {
Preprocess::Preprocess(PreprocessParameters params) : parameters(std::move(params)) {}

void Preprocess::process(const PreprocessData &data) const {
    Logger::log(constants::pipelines::PROCESSING_TREATMENT_MESSAGE);

    for (const auto &sample : data.treatmentSamples) {
        processSample(sample);
    }

    if (!data.controlSamples.has_value()) {
        Logger::log("No control samples provided");
        return;
    }

    Logger::log(constants::pipelines::PROCESSING_CONTROL_MESSAGE);

    for (const auto &sample : data.controlSamples.value()) {
        processSample(sample);
    }
};

void Preprocess::processSample(const PreprocessSampleType &sample) const {
    std::visit(
        overloaded{[this](const PreprocessSampleSingle &sample) { processSingleEnd(sample); },
                   [this](const PreprocessSamplePaired &sample) { processPairedEnd(sample); }},
        sample);
};

/**
 * Processes a single-end sample.
 *
 * @param sample The sample to be processed.
 */
void Preprocess::processSingleEnd(const PreprocessSampleSingle &sample) const {
    const SingleEndPreprocessor preprocessor{parameters};
    preprocessor.process(sample);
};

/**
 * Process a paired-end sample.
 *
 * @param sample The sample to be processed.
 */
void Preprocess::processPairedEnd(const PreprocessSamplePaired &sample) const {
    const PairedEndPreprocessor preprocessor{parameters};
    preprocessor.process(sample);
};

}  // namespace pipelines::preprocess
