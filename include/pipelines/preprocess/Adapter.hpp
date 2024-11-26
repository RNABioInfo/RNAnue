#pragma once

// Standard
#include <ostream>
#include <string>
#include <vector>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/utility/range/to.hpp>

// Internal
#include "PreprocessParameters.hpp"
#include "TrimConfig.hpp"
#include "VariantOverload.hpp"

namespace pipelines::preprocess {

struct Adapter {
    const seqan3::dna5_vector sequence;
    const double maxMissMatchFraction;
    const TrimConfig::Mode trimmingMode;

    static auto loadAdapters(const AdapterInput &adapterInput, const double maxMissMatchFraction,
                             const TrimConfig::Mode trimmingMode) -> std::vector<Adapter> {
        std::vector<Adapter> adapters;

        auto addAdapterFromSequence = [&adapters, &trimmingMode,
                                       maxMissMatchFraction](const seqan3::dna5_vector &sequence) {
            if (sequence.size() == 0) {
                return;
            }
            adapters.push_back(Adapter{sequence, maxMissMatchFraction, trimmingMode});
        };

        auto addAdapterFromFile = [&adapters, &trimmingMode,
                                   maxMissMatchFraction](const std::string &filePath) {
            seqan3::sequence_file_input adapterFile{filePath};
            for (const auto &record : adapterFile) {
                adapters.push_back(Adapter{record.sequence(), maxMissMatchFraction, trimmingMode});
            }
        };

        std::visit(
            overloaded{[](const std::monostate &) {}, addAdapterFromSequence, addAdapterFromFile},
            adapterInput);

        for (const auto &adapter : adapters) {
            Logger::log("Loaded adapter: ", adapter, " with ",
                        std::to_string(adapter.maxMissMatchFraction),
                        " allowed miss-match fraction.");
        }

        return adapters;
    }
};

inline std::ostream &operator<<(std::ostream &outputStream, const Adapter &adapter) {
    outputStream << (adapter.sequence | seqan3::views::to_char | seqan3::ranges::to<std::string>())
                 << " ";
    return outputStream;
};

}  // namespace pipelines::preprocess
