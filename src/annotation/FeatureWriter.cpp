#include "FeatureWriter.hpp"

// Standard
#include <fstream>

namespace annotation {

void FeatureWriter::write(const FeatureTreeMap &featureTreeMap, const std::string &outputPath,
                          const FileType::Value fileType) {
    std::ofstream outputFile(outputPath);
    if (!outputFile.is_open()) {
        throw std::runtime_error("Could not open file for writing: " + outputPath);
    }

    // Write the file header based on the fileType
    if (fileType == FileType::GFF) {
        outputFile << "##gff-version 3\n";
    } else if (fileType == FileType::GTF) {
        outputFile << "##gtf-version 2.2\n";
    }

    for (const auto &[referenceID, tree] : featureTreeMap) {
        for (const auto &interval : tree.intervals()) {
            const auto &feature = interval.data;
            outputFile << referenceID << '\t' << "." << '\t' << feature.type << '\t'
                       << feature.startPosition + 1 << '\t' << feature.endPosition + 1 << '\t'
                       << "." << '\t'
                       << (feature.strand == dataTypes::GenomicStrand::FORWARD ? '+' : '-') << '\t'
                       << "." << '\t';

            // Attributes field
            if (fileType == FileType::GFF) {
                outputFile << "ID=" << feature.id;
            } else if (fileType == FileType::GTF) {
                outputFile << "gene_id \"" << feature.id << "\"; ";
            }
            outputFile << '\n';
        }
    }

    outputFile.close();
}

}  // namespace annotation
