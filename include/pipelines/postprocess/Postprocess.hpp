#pragma once

// Standard
#include <utility>

// Internal
#include "GenomicStrandSpecificity.hpp"
#include "Logger.hpp"
#include "PostprocessData.hpp"
#include "PostprocessParameters.hpp"
#include "postprocess/InteractionGenerator.hpp"
#include "postprocess/InteractionParser.hpp"

namespace pipelines::postprocess {

class Postprocess {
   public:
    explicit Postprocess(PostprocessParameters params) : parameters(std::move(params)) {};

    void process(const PostprocessData& data) {
        auto parsingResult = parseInteractions(data.samples);

        Logger::log("Parsed Interactions: ", parsingResult.interactions.size());

        Logger::log("Merging Interactions");
        InteractionGenerator generator{
            {.minSegmentOverlapFraction = parameters.minSegmentOverlapFraction,
             .mergingStrandSpecificity = GenomicStrandSpecificity::SPECIFIC}};

        auto mergingResult = generator.merge(std::exchange(parsingResult.interactions, {}));

        Logger::log("Merged Interactions: ", mergingResult.superInteractions.size());
        for (const auto& interaction : mergingResult.superInteractions) {
            Logger::log(interaction);
        }
    };

   private:
    PostprocessParameters parameters;

    [[nodiscard]] static auto parseInteractions(
        const PostprocessData::SampleByConditionMap& samplesMap) -> InteractionParser::Result {
        InteractionParser parser{};

        for (const auto& [condition, samples] : samplesMap) {
            for (const auto& sample : samples) {
                parser.parse(sample.input.interactionsPath, sample.input.sampleName);
            }
        }

        return parser.getResultAndReset();
    }
};

}  // namespace pipelines::postprocess
