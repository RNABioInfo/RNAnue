#pragma once

#include <concepts>
#include <vector>

#include "VariantMeta.hpp"

namespace pipelines::detect {

/// A Post-processor is any callable T which,
/// given a vector of some variant type InVariant,
/// returns a vector whose element-type is exactly
/// detail::meta::transform_variant_t<T,InVariant>.
template <typename T, typename InVariant>
concept ReadGroupPostprocessor = requires(T postprocessor, std::vector<InVariant>&& contexts) {
    {
        postprocessor(std::move(contexts))
    } -> std::same_as<std::vector<detail::meta::transform_variant_t<T, InVariant> > >;
};

}  // namespace pipelines::detect
