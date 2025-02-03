#pragma once

// seqan3
#include <cstdint>
#include <seqan3/alignment/configuration/align_config_method.hpp>

namespace pipelines::preprocess {

struct TrimConfig {
   public:
    enum class Mode : std::uint8_t { FIVE_PRIME, THREE_PRIME };

    /**
     * @brief Returns the semi-global alignment configuration for the given mode, assumes adapter as
     * sequence1.
     *
     * @param mode The mode to use for trimming.
     * @return seqan3::align_cfg::method_global The semi-global alignment configuration.
     */
    static auto alignmentConfigFor(TrimConfig::Mode mode) -> seqan3::align_cfg::method_global {
        seqan3::align_cfg::method_global config;

        switch (mode) {
            case TrimConfig::Mode::FIVE_PRIME:
                config.free_end_gaps_sequence1_leading = true;
                config.free_end_gaps_sequence2_leading = true;
                config.free_end_gaps_sequence1_trailing = false;
                config.free_end_gaps_sequence2_trailing = true;
                break;
            case TrimConfig::Mode::THREE_PRIME:
                config.free_end_gaps_sequence1_leading = false;
                config.free_end_gaps_sequence2_leading = true;
                config.free_end_gaps_sequence1_trailing = true;
                config.free_end_gaps_sequence2_trailing = true;
                break;
        }

        return config;
    }
};
}  // namespace pipelines::preprocess
