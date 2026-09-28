#pragma once

// Standard
#include <cstddef>
#include <filesystem>
#include <ranges>
#include <set>
#include <set>
#include <utility>
#include <vector>

// seqan3
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/views/async_input_buffer.hpp>

// Internal
#include "Adapter.hpp"
#include "PreprocessParameters.hpp"
#include "PreprocessSample.hpp"
#include "TrimConfig.hpp"

namespace pipelines::preprocess {

class SingleEndPreprocessor {
   public:
    SingleEndPreprocessor(PreprocessParameters parameters)
        : parameters(std::move(parameters)),
          adapters5(Adapter::loadAdapters(this->parameters.adapter5Forward,
                                          this->parameters.maxMissMatchFractionTrimming,
                                          TrimConfig::Mode::FIVE_PRIME)),
          adapters3(Adapter::loadAdapters(this->parameters.adapter3Forward,
                                          this->parameters.maxMissMatchFractionTrimming,
                                          TrimConfig::Mode::THREE_PRIME)) {};

    void process(const PreprocessSampleSingle& sample) const;

   private:
    struct ChunkResult {
        size_t polyGChangedRecords{};
        size_t polyGRemovedBases{};
        void operator+=(const ChunkResult& other);

        [[nodiscard]] auto getPassedRecords() const -> size_t { return passedRecords; };
        [[nodiscard]] auto getFailedRecords() const -> size_t { return failedRecords; };

        void incrementPassedRecords() { ++passedRecords; };
        void incrementFailedRecords() { ++failedRecords; };

       private:
        size_t passedRecords{0};
        size_t failedRecords{0};
    };

    PreprocessParameters parameters;

    std::vector<Adapter> adapters5;
    std::vector<Adapter> adapters3;

    [[nodiscard]] auto processWithDeduplication(const PreprocessSampleSingle& sample) const
        -> ChunkResult;
    [[nodiscard]] auto processWithoutDeduplication(const PreprocessSampleSingle& sample) const
        -> ChunkResult;

    auto processInput(const PreprocessSampleSingle& sample,
                      const std::set<size_t>* retained) const -> ChunkResult;

    template <typename T>
    auto processChunk(T& recordIterator, const fs::path& tmpOutDir) const -> ChunkResult;
};

}  // namespace pipelines::preprocess
