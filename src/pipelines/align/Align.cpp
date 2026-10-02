#include "Align.hpp"

#include "AlignmentProvenance.hpp"

// Standard
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <optional>
#include <variant>
#include <vector>

// segmehl
#include "AnnotationFilePicker.hpp"
#include "GenomePreprocessor.hpp"
#include "ReferenceGenomeFilePicker.hpp"

// Internal
#include "AlignData.hpp"
#include "AlignSample.hpp"
#include "Constants.hpp"
#include "Logger.hpp"
#include "ReferenceGenome.hpp"
#include "SamFileUtility.hpp"
#include "SequenceFileUtility.hpp"
#include "Utility.hpp"
#include "VariantOverload.hpp"

namespace pipelines::align {

void Align::process(const AlignData &data) {
    Logger::log(constants::pipelines::PROCESSING_TREATMENT_MESSAGE);

    preprocessReferences();
    if (parameters.aligner == AlignmentBackend::Star) {
        backend.emplace<StarAligner>(parameters);
    } else {
        backend.emplace<SegemehlAligner>(parameters);
    }
    std::visit([](auto &mapper) { mapper.buildIndex(); }, backend);

    for (const auto &sample : data.treatmentSamples) {
        processSample(sample);
    }

    if (!data.controlSamples) {
        return;
    }

    Logger::log(constants::pipelines::PROCESSING_CONTROL_MESSAGE);

    for (const auto &sample : *data.controlSamples) {
        processSample(sample);
    }
}

void Align::preprocessReferences() {
    if (!parameters.maskMultiCopyGenes) {
        Logger::log("Gene copy masking disabled.");
        return;
    }

    const auto genomePreprocessor = GenomePreprocessor{parameters};
    const GenomePreprocessor::Output genomePreprocessingOutput = {
        .preprocessedGenomePath = utility::ReferenceGenomeFilePicker::maskedFilePath(parameters),
        .preprocessedAnnotationPath = utility::AnnotationFilePicker::maskedFilePath(parameters)};

    genomePreprocessor.process(
        {.refGenomePath = parameters.referenceGenome, .annotationPath = parameters.featuresInPath},
        genomePreprocessingOutput);

    Logger::log("Replacing provided reference genome and annotation with masked.");
    parameters.referenceGenome = genomePreprocessingOutput.preprocessedGenomePath;
    parameters.featuresInPath = genomePreprocessingOutput.preprocessedAnnotationPath;
}

void Align::processSample(const AlignSampleType &sample) {
    std::visit(overloaded{[this](const AlignSampleSingle &sample) { processSingleEnd(sample); },
                          [this](const AlignSampleMergedPaired &sample) {
                              processMergedPairedEnd(sample);
                          }},
               sample);
}

void Align::processSingleEnd(const AlignSampleSingle &sample) {
    Logger::log("Processing single end reads");

    alignSingleReads(sample.input.inputFastqPath, sample.output.outputAlignmentsPath);
    sortAlignmentsByQueryName(sample.output.outputAlignmentsPath,
                              sample.output.outputAlignmentsPath);
}

void Align::processMergedPairedEnd(const AlignSampleMergedPaired &sample) {
    Logger::log("Processing merged paired end reads");

    alignSingleReads(sample.input.inputMergedFastqPath,
                     sample.output.outputAlignmentsMergedReadsPath);

    alignSingleReads(sample.input.inputSingletonForwardFastqPath,
                     sample.output.outputAlignmentsSingletonForwardReadsPath);

    alignSingleReads(sample.input.inputSingletonReverseFastqPath,
                     sample.output.outputAlignmentsSingletonReverseReadsPath);

    alignPairedReads(sample.input.inputPairedForwardFastqPath,
                     sample.input.inputPairedReverseFastqPath,
                     sample.output.outputAlignmentsPairedReadsPath);

    std::vector<fs::path> samFiles{sample.output.outputAlignmentsMergedReadsPath,
                                   sample.output.outputAlignmentsSingletonForwardReadsPath,
                                   sample.output.outputAlignmentsSingletonReverseReadsPath,
                                   sample.output.outputAlignmentsPairedReadsPath};

    Logger::log("Merging alignment files");

    helper::mergeSamFiles(samFiles, sample.output.outputAlignmentsPath, referenceFromGenome());

    sortAlignmentsByQueryName(sample.output.outputAlignmentsPath,
                              sample.output.outputAlignmentsPath);
}

void Align::alignSingleReads(const fs::path &input, const fs::path &output) {
    std::visit([&](auto &mapper) { mapper.alignSingleReads(input, output); }, backend);
    if (parameters.aligner == AlignmentBackend::Segemehl) {
        utility::stampAlignmentFile(output, parameters.aligner);
    }
}

void Align::alignPairedReads(const fs::path &forward, const fs::path &reverse,
                             const fs::path &output) {
    std::visit([&](auto &mapper) { mapper.alignPairedReads(forward, reverse, output); }, backend);
    if (parameters.aligner == AlignmentBackend::Segemehl) {
        utility::stampAlignmentFile(output, parameters.aligner);
    }
}

auto Align::referenceFromGenome() const -> dataTypes::SamReference {
    const ReferenceGenome genome{parameters.referenceGenome};
    auto ids = genome.getReferenceIndexMapping().sortedReferenceIDs();
    std::vector<size_t> lengths;
    lengths.reserve(ids.size());
    for (size_t i = 0; i < ids.size(); ++i) {
        lengths.push_back(genome.getSequence(static_cast<int>(i)).size());
    }
    return dataTypes::SamReference{std::move(ids), std::move(lengths)};
}

void Align::sortAlignmentsByQueryName(const fs::path &alignmentsPath,
                                      const fs::path &sortedAlignmentsPath) const {
    Logger::log("Sorting alignments");
    SamFileUtility::sortByQueryName(alignmentsPath, sortedAlignmentsPath, parameters.threadCount);
    Logger::log("Sorting alignments done");
}

}  // namespace pipelines::align
