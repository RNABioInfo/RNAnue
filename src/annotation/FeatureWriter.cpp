// FeatureWriter.cpp
#include "FeatureWriter.hpp"

#include "FeatureAnnotator.hpp"
#include "FileType.hpp"
#include "GenomicFeature.hpp"
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "MaskedFeatureCluster.hpp"

// Standard
#include <array>
#include <cassert>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <fstream>
#include <ios>
#include <ostream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

namespace annotation {
namespace {
constexpr auto kFieldSeparator = '\t';
constexpr auto kLineSeparator = '\n';
constexpr auto kDotField = std::string_view{"."};
constexpr auto kOneBasedOffset = std::uint64_t{1};
constexpr auto kUnsignedBuffer = std::size_t{32};
constexpr auto kLineBufferSize = std::size_t{256};
constexpr auto kStreamBufferSizeBytes = std::size_t{1U << 20};  // 1 MiB
}  // namespace

auto FeatureWriter::throwLogged(const std::string_view message) -> void {
    Logger::log<SourceLocation{}, LogLevel::ERROR>(std::string{message});
    throw std::runtime_error(std::string{message});
}

auto FeatureWriter::writeHeader(std::ostream& outputStream, const FileType::Value fileType)
    -> void {
    if (fileType == FileType::GFF) {
        outputStream << "##gff-version 3\n";
        return;
    }
    if (fileType == FileType::GTF) {
        outputStream << "##gtf-version 2.2\n";
        return;
    }
    throwLogged("Unsupported file type.");
}

auto FeatureWriter::appendUnsigned(std::string& lineBuffer, const std::uint64_t value) -> void {
    std::array<char, kUnsignedBuffer> buffer{};
    auto* const first = buffer.data();
    auto* const last = buffer.data() + buffer.size();

    const auto [ptr, errorCode] = std::to_chars(first, last, value);
    if (errorCode != std::errc{}) {
        throwLogged("Failed to format integer.");
    }

    lineBuffer.append(first, static_cast<std::size_t>(ptr - first));
}

auto FeatureWriter::appendAttribute(std::string& lineBuffer, const FileType::Value fileType,
                                    const std::string_view key, const std::string_view value,
                                    const bool isFirstAttribute) -> void {
    if (fileType == FileType::GFF) {
        if (!isFirstAttribute) {
            lineBuffer.push_back(';');
        }
        lineBuffer.append(key);
        lineBuffer.push_back('=');
        lineBuffer.append(value);
        return;
    }

    if (fileType == FileType::GTF) {
        if (!isFirstAttribute) {
            lineBuffer.push_back(' ');
        }
        lineBuffer.append(key);
        lineBuffer.append(" \"");
        lineBuffer.append(value);
        lineBuffer.append("\";");
        return;
    }

    throwLogged("Unsupported file type.");
}

auto FeatureWriter::formatFeatureLine(std::string& lineBuffer, const std::string_view referenceId,
                                      const GenomicFeature& feature, const FileType::Value fileType)
    -> std::string_view {
    return formatFeatureLine(lineBuffer, referenceId, feature, fileType, {});
}

auto FeatureWriter::formatFeatureLine(std::string& lineBuffer, const std::string_view referenceId,
                                      const GenomicFeature& feature, const FileType::Value fileType,
                                      const std::span<const Attribute> extraAttributes)
    -> std::string_view {
    lineBuffer.clear();

    const auto& region = feature.getGenomicRegion();

    // seqname, source, feature, start, end, score, strand, frame, attributes
    lineBuffer.append(referenceId);
    lineBuffer.push_back(kFieldSeparator);

    lineBuffer.append(kDotField);
    lineBuffer.push_back(kFieldSeparator);

    lineBuffer.append(feature.getType());
    lineBuffer.push_back(kFieldSeparator);

    appendUnsigned(lineBuffer, static_cast<std::uint64_t>(region.getStart()) + kOneBasedOffset);
    lineBuffer.push_back(kFieldSeparator);

    appendUnsigned(lineBuffer, static_cast<std::uint64_t>(region.getEnd()));
    lineBuffer.push_back(kFieldSeparator);

    lineBuffer.append(kDotField);
    lineBuffer.push_back(kFieldSeparator);

    lineBuffer.push_back(region.getStrand());
    lineBuffer.push_back(kFieldSeparator);

    lineBuffer.append(kDotField);
    lineBuffer.push_back(kFieldSeparator);

    auto isFirstAttribute = true;

    if (fileType == FileType::GFF) {
        appendAttribute(lineBuffer, fileType, "ID", feature.getID(), isFirstAttribute);
    } else if (fileType == FileType::GTF) {
        appendAttribute(lineBuffer, fileType, "gene_id", feature.getID(), isFirstAttribute);
    } else {
        throwLogged("Unsupported file type.");
    }
    isFirstAttribute = false;

    for (const auto& [key, value] : feature.getAttributes()) {
        appendAttribute(lineBuffer, fileType, key, value, isFirstAttribute);
    }

    for (const auto& attribute : extraAttributes) {
        appendAttribute(lineBuffer, fileType, attribute.key, attribute.value, isFirstAttribute);
    }

    lineBuffer.push_back(kLineSeparator);
    return std::string_view{lineBuffer};
}

auto FeatureWriter::write(const FeatureTreeMap& featureTreeMap,
                          const std::deque<std::string>& sortedReferenceIDs,
                          const std::string& outputPath, const FileType::Value fileType) -> void {
    std::ofstream outputFile(outputPath, std::ios::binary);
    if (!outputFile.is_open()) {
        throwLogged(std::string{"Could not open file for writing: "} + outputPath);
    }

    std::string streamBuffer;
    streamBuffer.resize(kStreamBufferSizeBytes);
    outputFile.rdbuf()->pubsetbuf(streamBuffer.data(),
                                  static_cast<std::streamsize>(streamBuffer.size()));

    writeHeader(outputFile, fileType);

    std::string lineBuffer;
    lineBuffer.reserve(kLineBufferSize);

    for (const auto& [referenceIndex, tree] : featureTreeMap) {
        const auto& referenceId = getReferenceID(referenceIndex, sortedReferenceIDs);

        for (const auto& interval : tree.intervals()) {
            const auto& feature = interval.data;

            const auto lineView = formatFeatureLine(lineBuffer, referenceId, feature, fileType);
            outputFile.write(lineView.data(), static_cast<std::streamsize>(lineView.size()));

            if (!outputFile.good()) {
                throwLogged(std::string{"I/O error while writing: "} + outputPath);
            }
        }
    }
}

auto FeatureWriter::write(const std::vector<MaskedFeatureCluster>& featureClusters,
                          const std::deque<std::string>& sortedReferenceIDs,
                          const std::string& outputPath, FileType::Value fileType) -> void {
    std::ofstream outputFile(outputPath, std::ios::binary);
    if (!outputFile.is_open()) {
        throwLogged(std::string{"Could not open file for writing: "} + outputPath);
    }

    std::string streamBuffer;
    streamBuffer.resize(kStreamBufferSizeBytes);
    outputFile.rdbuf()->pubsetbuf(streamBuffer.data(),
                                  static_cast<std::streamsize>(streamBuffer.size()));

    writeHeader(outputFile, fileType);

    std::string lineBuffer;
    lineBuffer.reserve(kLineBufferSize);

    for (const auto& cluster : featureClusters) {
        const auto& baseFeatureGroup = cluster.baseFeatureGroup;
        const auto& baseFeature = baseFeatureGroup.getRoot().feature;

        const auto& rootReferenceID = getReferenceID(
            baseFeature.getGenomicRegion().getReferenceIDIndex(), sortedReferenceIDs);

        if (cluster.isMultiCopy()) {
            const auto maskedFeatureString = cluster.csvSubFeatureRootIDs();
            std::array attributes = {
                Attribute{.key = "MaskedFeatures", .value = maskedFeatureString}};

            const auto lineView =
                formatFeatureLine(lineBuffer, rootReferenceID, baseFeature, fileType, attributes);
            outputFile.write(lineView.data(), static_cast<std::streamsize>(lineView.size()));
        } else {
            const auto lineView =
                formatFeatureLine(lineBuffer, rootReferenceID, baseFeature, fileType);
            outputFile.write(lineView.data(), static_cast<std::streamsize>(lineView.size()));
        }

        if (!outputFile.good()) {
            throwLogged(std::string{"I/O error while writing: "} + outputPath);
        }

        // Write child features
        for (const auto& child : baseFeatureGroup.getNodes()) {
            if (&child == &baseFeatureGroup.getRoot()) {  // if getNodes() includes root
                continue;
            }

            const auto& parent = baseFeatureGroup.getNode(child.parentIndex);

            const auto& feature = child.feature;

            const auto& referenceID = getReferenceID(
                feature.getGenomicRegion().getReferenceIDIndex(), sortedReferenceIDs);

            std::array attributes = {Attribute{.key = "Parent", .value = parent.feature.getID()}};

            const auto lineView =
                formatFeatureLine(lineBuffer, referenceID, feature, fileType, attributes);
            outputFile.write(lineView.data(), static_cast<std::streamsize>(lineView.size()));

            if (!outputFile.good()) {
                throwLogged(std::string{"I/O error while writing: "} + outputPath);
            }
        }
    }
}

auto FeatureWriter::getReferenceID(int referenceIndex,
                                   const std::deque<std::string>& sortedReferenceIDs)
    -> std::string_view {
    if (referenceIndex < 0 ||
        static_cast<std::size_t>(referenceIndex) >= sortedReferenceIDs.size()) {
        throwLogged("Reference index out of bounds for sortedReferenceIDs.");
    }

    return sortedReferenceIDs[referenceIndex];
}
}  // namespace annotation
