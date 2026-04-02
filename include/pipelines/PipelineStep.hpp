#pragma once

#include <cstdint>
namespace pipelines {
enum class PipelineStep : std::uint8_t { PREPROCESS, ALIGN, DETECT, ASSEMBLE, EVALUATE, COMPLETE };
}
