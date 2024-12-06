#pragma once

// Standard
#include <cstddef>

// boost
#include <boost/program_options.hpp>

// seqan3
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/views/async_input_buffer.hpp>

// Class
#include "PreprocessData.hpp"
#include "PreprocessParameters.hpp"
#include "PreprocessSample.hpp"

using seqan3::operator""_dna5;

namespace pipelines::preprocess {
class Preprocess {
   public:
    explicit Preprocess(PreprocessParameters params);

    void process(const PreprocessData &data) const;

   private:
    PreprocessParameters parameters;

    void processSample(const PreprocessSampleType &sample) const;

    void processSingleEnd(const PreprocessSampleSingle &sample) const;
    void processPairedEnd(const PreprocessSamplePaired &sample) const;
};

}  // namespace pipelines::preprocess
