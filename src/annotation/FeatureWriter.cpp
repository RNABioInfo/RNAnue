#include "FeatureWriter.hpp"

#include "FeatureAnnotator.hpp"
#include "FileType.hpp"

// Standard
#include <fstream>
#include <stdexcept>
#include <string>

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
                       << feature.genomicRegion.getStart() + 1 << '\t'
                       << feature.genomicRegion.getEnd() << '\t' << "." << '\t'
                       << feature.genomicRegion.getStrand() << '\t' << "." << '\t';

            // Attributes field
            if (fileType == FileType::GFF) {
                outputFile << "ID=" << feature.featureID;
            } else if (fileType == FileType::GTF) {
                outputFile << "gene_id \"" << feature.featureID << "\"; ";
            }
            outputFile << '\n';
        }
    }

    outputFile.close();
}

}  // namespace annotation
