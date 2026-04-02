#pragma once

// Standard
#include <cstddef>
#include <cstdint>
#include <deque>
#include <iosfwd>
#include <span>
#include <string>
#include <string_view>
#include <vector>

// Internal
#include "FeatureAnnotator.hpp"
#include "FileType.hpp"
#include "GenomicFeature.hpp"
#include "MaskedFeatureCluster.hpp"

namespace annotation {

class FeatureWriter {
   public:
    FeatureWriter() = delete;
    FeatureWriter(const FeatureWriter&) = delete;
    FeatureWriter(FeatureWriter&&) = delete;
    auto operator=(const FeatureWriter&) -> FeatureWriter& = delete;
    auto operator=(FeatureWriter&&) -> FeatureWriter& = delete;
    ~FeatureWriter() = delete;

    struct Attribute {
        std::string_view key;
        std::string_view value;
    };

    static auto write(const FeatureTreeMap& featureTreeMap,
                      const std::deque<std::string>& sortedReferenceIDs,
                      const std::string& outputPath, FileType::Value fileType) -> void;

    static auto write(const std::vector<MaskedFeatureCluster>& featureClusters,
                      const std::deque<std::string>& sortedReferenceIDs,
                      const std::string& outputPath, FileType::Value fileType) -> void;

   private:
    static auto writeHeader(std::ostream& outputStream, FileType::Value fileType) -> void;

    // Reusable building block: formats one feature line into `lineBuffer` and returns a view.
    static auto formatFeatureLine(std::string& lineBuffer, std::string_view referenceId,
                                  const GenomicFeature& feature, FileType::Value fileType)
        -> std::string_view;

    static auto formatFeatureLine(std::string& lineBuffer, std::string_view referenceId,
                                  const GenomicFeature& feature, FileType::Value fileType,
                                  std::span<const Attribute> extraAttributes) -> std::string_view;

    static auto appendUnsigned(std::string& lineBuffer, std::uint64_t value) -> void;

    // New: reusable method for arbitrary attributes.
    static auto appendAttribute(std::string& lineBuffer, FileType::Value fileType,
                                std::string_view key, std::string_view value, bool isFirstAttribute)
        -> void;

    [[nodiscard]] static auto getReferenceID(int referenceIndex,
                                             const std::deque<std::string>& sortedReferenceIDs)
        -> std::string_view;

    [[noreturn]] static auto throwLogged(std::string_view message) -> void;
};

}  // namespace annotation
