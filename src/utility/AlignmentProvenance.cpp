#include "AlignmentProvenance.hpp"

#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sam_file/output.hpp>

#include "AlignmentFileInput.hpp"
#include "ExternalProcess.hpp"
#include "SamFileUtility.hpp"

namespace utility {
void stampAlignmentFile(const std::filesystem::path& path, AlignmentBackend backend,
                        const std::string& command) {
    TemporaryDirectory workspace{path.parent_path()};
    const auto temporary = workspace.path() / "tagged.bam";
    const auto inspection = SamFileUtility::inspect(path);
    if (!inspection.isReadable()) {
        throw std::runtime_error("Cannot stamp unreadable alignment BAM");
    }

    {
        AlignmentFileInput input{path};
        const auto& header = input.header();
        dataTypes::SamReference reference{header};
        seqan3::sam_file_output output{temporary, reference.referenceIDs,
                                       reference.referenceLengths};
        output.header().program_infos = header.program_infos;
        output.header().comments = header.comments;
        output.header().sorting = header.sorting;
        stampAlignment(output.header(), backend, command);

        if (inspection.hasRecords()) {
            for (auto& record : input) {
                output.push_back(record);
            }
        }
    }

    if (!SamFileUtility::inspect(temporary).isReadable()) {
        throw std::runtime_error("Cannot stamp alignment BAM");
    }
    std::filesystem::rename(temporary, path);
}
}  // namespace utility
