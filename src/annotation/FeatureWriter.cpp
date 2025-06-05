#include "FeatureWriter.hpp"

#include "FeatureAnnotator.hpp"
#include "FileType.hpp"

// Standard
#include <cassert>
#include <cstddef>
#include <deque>
#include <fstream>
#include <stdexcept>
#include <string>

namespace annotation {

void FeatureWriter::write(const FeatureTreeMap &featureTreeMap,
                          const std::deque<std::string> &sortedReferenceIDs,
                          const std::string &outputPath, const FileType::Value fileType) {
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

    for (const auto &[referenceIndex, tree] : featureTreeMap) {
        assert(sortedReferenceIDs.size() > size_t(referenceIndex));

        for (const auto &interval : tree.intervals()) {
            const auto &feature = interval.data;
            const auto &referenceID = sortedReferenceIDs[referenceIndex];
            outputFile << referenceID << '\t' << "." << '\t' << feature.getType() << '\t'
                       << feature.getGenomicRegion().getStart() + 1 << '\t'
                       << feature.getGenomicRegion().getEnd() << '\t' << "." << '\t'
                       << feature.getGenomicRegion().getStrand() << '\t' << "." << '\t';

            // Attributes field
            if (fileType == FileType::GFF) {
                outputFile << "ID=" << feature.getID();
            } else if (fileType == FileType::GTF) {
                outputFile << "gene_id \"" << feature.getID() << "\"; ";
            }
            outputFile << '\n';
        }
    }

    outputFile.close();
}

}  // namespace annotation
