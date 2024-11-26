#include "Runner.hpp"

#include "Align.hpp"
#include "Analyze.hpp"
#include "AnalyzeData.hpp"
#include "Closing.hpp"
#include "Detect.hpp"
#include "ParameterParser.hpp"
#include "Preprocess.hpp"

void Runner::runPipeline(int argc, const char *const argv[]) {  // NOLINT
    const auto parameters = ParameterParser::getParameters(argc, argv);

    std::visit(Pipeline(), parameters);

    Closing::printQuote();
}

void Runner::runPreprocessPipeline(const preprocess::PreprocessParameters &parameters) {
    Logger::log("Running preprocess pipeline");

    if (!parameters.preprocessEnabled) {
        Logger::log("Preprocess pipeline is disabled in the parameters");
        return;
    }

    const auto data = preprocess::PreprocessData(parameters.outputDir, parameters.treatmentsDir,
                                                 parameters.controlDir);

    const auto pipeline = preprocess::Preprocess(parameters);
    pipeline.process(data);
}

void Runner::runAlignPipeline(const align::AlignParameters &parameters) {
    Logger::log("Running align pipeline");

    const auto inputDirs = InputDirectories(parameters.outputDir, preprocess::pipelinePrefix);

    const auto data = align::AlignData(parameters.outputDir, inputDirs.treatmentInputDir,
                                       inputDirs.controlInputDir);

    auto pipeline = align::Align(parameters);
    pipeline.process(data);
}

void Runner::runDetectPipeline(const detect::DetectParameters &parameters) {
    Logger::log("Running detect pipeline");

    const auto inputDirs = InputDirectories(parameters.outputDir, align::pipelinePrefix);

    const auto data = detect::DetectData(parameters.outputDir, inputDirs.treatmentInputDir,
                                         inputDirs.controlInputDir);

    auto pipeline = detect::Detect(parameters);
    pipeline.process(data);
}

void Runner::runAnalyzePipeline(const analyze::AnalyzeParameters &parameters) {
    Logger::log("Running analyze pipeline");

    const auto inputDirs = InputDirectories(parameters.outputDir, detect::pipelinePrefix);

    const auto data = analyze::AnalyzeData(parameters.outputDir, inputDirs.treatmentInputDir,
                                           inputDirs.controlInputDir);

    auto pipeline = analyze::Analyze(parameters);
    pipeline.process(data);
}

void Runner::runCompletePipeline(const CompleteParameters &parameters) {
    Logger::log("Running complete pipeline");

    runPreprocessPipeline(parameters.preprocessParameters);
    runAlignPipeline(parameters.alignParameters);
    runDetectPipeline(parameters.detectParameters);
    runAnalyzePipeline(parameters.analyzeParameters);
}

void Runner::Pipeline::operator()(const preprocess::PreprocessParameters &params) {
    Runner::runPreprocessPipeline(params);
};
void Runner::Pipeline::operator()(const align::AlignParameters &params) {
    Runner::runAlignPipeline(params);
};
void Runner::Pipeline::operator()(const detect::DetectParameters &params) {
    Runner::runDetectPipeline(params);
};
void Runner::Pipeline::operator()(const analyze::AnalyzeParameters &params) {
    Runner::runAnalyzePipeline(params);
};
void Runner::Pipeline::operator()(const CompleteParameters &params) {
    Runner::runCompletePipeline(params);
};
