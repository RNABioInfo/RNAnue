#include "StarAlignmentDecoder.hpp"

#include <algorithm>
#include <charconv>
#include <climits>
#include <map>
#include <stdexcept>

namespace pipelines::align {
namespace {
using namespace dataTypes;
auto flag(const SamRecord& r, unsigned bit) -> bool {
    return (static_cast<unsigned>(r.flag()) & bit) != 0;
}
auto mate(const SamRecord& r) -> size_t { return flag(r, 0x80) ? 1 : 0; }

[[noreturn]] void invalid(const SamRecord& r, const std::string& reason) {
    throw std::runtime_error("Malformed STAR alignment for '" + r.id() + "': " + reason);
}

auto cigarString(const SamRecord& r) -> std::string {
    std::string result;

    for (const auto& c : r.cigar_sequence()) {
        result += std::to_string(get<0>(c)) + seqan3::to_char(get<1>(c));
    }

    return result;
}
struct MdRun {
    size_t count;
    bool deletion;
    bool mismatch;
};
class MdCursor {
    std::vector<MdRun> runs;
    size_t next{};

   public:
    explicit MdCursor(const SamRecord& r) {
        const auto& md = r.tags().get<"MD"_tag>();
        size_t pos = 0;

        while (pos < md.size()) {
            if (md[pos] >= '0' && md[pos] <= '9') {
                size_t count{};
                auto parsed = std::from_chars(md.data() + pos, md.data() + md.size(), count);
                if (parsed.ec != std::errc{}) invalid(r, "invalid MD length");

                pos = static_cast<size_t>(parsed.ptr - md.data());
                if (count != 0U) {
                    runs.push_back({count, false, false});
                }

            } else if (md[pos] == '^') {
                const auto start = ++pos;
                while (pos < md.size() && md[pos] >= 'A' && md[pos] <= 'Z') {
                    ++pos;
                }

                if (pos == start) {
                    invalid(r, "empty MD deletion");
                }

                runs.push_back({pos - start, true, false});
            } else if (md[pos] >= 'A' && md[pos] <= 'Z') {
                ++pos;
                runs.push_back({1, false, true});
            } else {
                invalid(r, "invalid MD token");
            }
        }
    }
    auto consume(size_t count, bool deletion, const SamRecord& r) -> size_t {
        size_t mismatches = 0;

        while (count != 0U) {
            if (next == runs.size() || runs[next].deletion != deletion) {
                invalid(r, "MD/CIGAR disagree");
            }

            const auto used = std::min(count, runs[next].count);
            if (runs[next].mismatch) {
                mismatches += used;
            }

            count -= used;
            runs[next].count -= used;
            if (runs[next].count == 0U) {
                ++next;
            }
        }
        return mismatches;
    }
    auto complete() const -> bool { return next == runs.size(); }
};
struct Piece {
    SamRecord record;
    size_t mateIndex{}, begin{}, end{}, originalBegin{}, originalEnd{}, rawIndex{};
};
auto pieces(const SamRecord& raw) -> std::vector<Piece> {
    MdCursor md{raw};
    std::vector<Piece> result;

    size_t q = 0;
    size_t begin = 0;
    size_t editDistance = 0;
    size_t totalEdits = 0;
    int64_t ref = *raw.reference_position();
    int64_t start = ref;
    std::vector<seqan3::cigar> cigar;

    const auto finish = [&] {
        if (cigar.empty() || begin == q) {
            invalid(raw, "empty alignment block");
        }

        auto r = raw;
        r.reference_position() = static_cast<int32_t>(start);
        r.cigar_sequence() = cigar;
        r.sequence() =
            seqan3::dna5_vector(raw.sequence().begin() + begin, raw.sequence().begin() + q);

        if (!raw.base_qualities().empty()) {
            r.base_qualities() = std::vector<seqan3::phred42>(raw.base_qualities().begin() + begin,
                                                              raw.base_qualities().begin() + q);
        }

        for (const auto tag : {"SA"_tag, "MD"_tag, "AS"_tag, "nM"_tag, "ch"_tag}) {
            r.tags().erase(tag);
        }

        r.tags()["NM"_tag] = static_cast<int32_t>(editDistance);
        const auto originalBegin = flag(raw, 0x10) ? raw.sequence().size() - q : begin;
        const auto originalEnd = flag(raw, 0x10) ? raw.sequence().size() - begin : q;
        r.tags()["XX"_tag] = static_cast<int32_t>(originalBegin);
        r.tags()["XY"_tag] = static_cast<int32_t>(originalEnd);
        result.push_back({std::move(r), mate(raw), begin, q, originalBegin, originalEnd});
        totalEdits += editDistance;
        cigar.clear();
        editDistance = 0;
    };

    for (size_t i = 0; i < raw.cigar_sequence().size(); ++i) {
        const auto& c = raw.cigar_sequence()[i];
        const size_t n = get<0>(c);
        const char op = seqan3::to_char(get<1>(c));

        if (n == 0U) {
            invalid(raw, "zero-length CIGAR operation");
        }

        if (op == 'S') {
            if (i != 0 && i + 1 != raw.cigar_sequence().size()) {
                invalid(raw, "internal soft clipping");
            }

            if (i == 0) {
                q += n;
                begin = q;
            }
            // Trailing clipping is attached later, only at the exterior of the entire hypothesis.
        } else if (op == 'N') {
            finish();
            ref += n;
            start = ref;
            begin = q;
        } else if (op == 'M' || op == '=' || op == 'X') {
            const auto mismatches = md.consume(n, false, raw);

            if ((op == '=' && mismatches) || (op == 'X' && mismatches != n)) {
                invalid(raw, "MD/CIGAR match type disagrees");
            }

            editDistance += mismatches;
            q += n;
            ref += n;
            cigar.push_back(c);
        } else if (op == 'I') {
            editDistance += n;
            q += n;
            cigar.push_back(c);
        } else if (op == 'D') {
            md.consume(n, true, raw);
            editDistance += n;
            ref += n;
            cigar.push_back(c);
        } else if (op == 'P') {
            cigar.push_back(c);
        } else {
            invalid(
                raw,
                "unsupported CIGAR operation (STAR must emit full sequences with soft clipping)");
        }

        if (q > raw.sequence().size() || ref > INT_MAX) {
            invalid(raw, "CIGAR coordinate overflow");
        }
    }
    finish();

    if (!md.complete() || totalEdits != static_cast<size_t>(raw.tags().get<"NM"_tag>())) {
        invalid(raw, "NM/MD/CIGAR disagree");
    }
    return result;
}
}  // namespace

auto StarAlignmentDecoder::decode(const std::vector<SamRecord>& records,
                                  bool excludeSoftClipping) const
    -> std::vector<StarAlignmentHypothesis> {
    std::map<int, std::vector<size_t>> groups;
    for (size_t i = 0; i < records.size(); ++i) {
        const auto& r = records[i];

        if (r.id() != records.front().id()) {
            invalid(r, "different query names in one read group");
        }

        if (flag(r, 4)) {
            continue;
        }

        if (!r.reference_id() || !r.reference_position() || *r.reference_position() < 0 ||
            *r.reference_id() < 0 ||
            static_cast<size_t>(*r.reference_id()) >= referenceIDs.size()) {
            invalid(r, "invalid reference coordinates");
        }

        for (const auto tag : {"HI"_tag, "NH"_tag, "AS"_tag, "NM"_tag, "MD"_tag}) {
            if (!r.tags().contains(tag)) {
                invalid(r, "missing required STAR attribute");
            }
        }

        for (const auto tag : {"HI"_tag, "NH"_tag, "AS"_tag, "NM"_tag}) {
            if (!std::holds_alternative<int32_t>(r.tags().at(tag))) {
                invalid(r, "noninteger STAR attribute");
            }
        }

        if (!std::holds_alternative<std::string>(r.tags().at("MD"_tag)) ||
            (r.tags().contains("SA"_tag) &&
             !std::holds_alternative<std::string>(r.tags().at("SA"_tag)))) {
            invalid(r, "invalid MD/SA attribute type");
        }

        if (r.tags().get<"HI"_tag>() < 1 || r.tags().get<"NH"_tag>() < r.tags().get<"HI"_tag>() ||
            r.tags().get<"NM"_tag>() < 0) {
            invalid(r, "invalid HI/NH/NM");
        }

        if (r.tags().contains("ch"_tag) && (!std::holds_alternative<char>(r.tags().at("ch"_tag)) ||
                                            std::get<char>(r.tags().at("ch"_tag)) != '1')) {
            invalid(r, "invalid ch marker");
        }

        if (!flag(r, 1) && (flag(r, 0x40) || flag(r, 0x80))) {
            invalid(r, "mate flags on unpaired record");
        }

        if (!referenceLengths.empty()) {
            uint64_t end = *r.reference_position();
            for (const auto& c : r.cigar_sequence()) {
                const auto op = seqan3::to_char(get<1>(c));
                if (op == 'M' || op == '=' || op == 'X' || op == 'D' || op == 'N') end += get<0>(c);
            }

            if (end > referenceLengths.at(*r.reference_id())) {
                invalid(r, "alignment exceeds reference length");
            }
        }

        if (flag(r, 1) && flag(r, 0x40) == flag(r, 0x80)) {
            invalid(r, "invalid mate flags");
        }

        if (r.sequence().empty() ||
            (!r.base_qualities().empty() && r.base_qualities().size() != r.sequence().size())) {
            invalid(r, "invalid sequence/qualities");
        }

        try {
            if (dataTypes::alignmentLength(r) == 0U) {
                invalid(r, "CIGAR has no aligned query bases");
            }
        } catch (const std::exception& e) {
            invalid(r, e.what());
        }

        groups[r.tags().get<"HI"_tag>()].push_back(i);
    }

    const auto descriptor = [&](const SamRecord& r) {
        return referenceIDs[*r.reference_id()] + "," + std::to_string(*r.reference_position() + 1) +
               "," + (flag(r, 0x10) ? "-" : "+") + "," + cigarString(r) + "," +
               std::to_string(r.mapping_quality()) + "," +
               std::to_string(r.tags().get<"NM"_tag>()) + ";";
    };

    std::vector<StarAlignmentHypothesis> result;
    for (const auto& [hi, indices] : groups) {
        StarAlignmentHypothesis hypothesis{.hitIndex = hi, .rawIndices = indices, .fragments = {}};
        std::vector<Piece> decoded;
        std::array<size_t, 2> first{SIZE_MAX, SIZE_MAX};
        std::array<size_t, 2> last{};

        const auto& representative = records[indices.front()];

        for (const auto index : indices) {
            const auto& r = records[index];

            if (flag(r, 1) != flag(representative, 1) ||
                r.tags().get<"NH"_tag>() != representative.tags().get<"NH"_tag>()) {
                invalid(r, "inconsistent hypothesis metadata");
            }

            hypothesis.paired = flag(r, 1);
            const auto m = mate(r);

            if ((hypothesis.queryLengths[m] != 0U) &&
                hypothesis.queryLengths[m] != r.sequence().size()) {
                invalid(r, "inconsistent query lengths");
            }

            hypothesis.queryLengths[m] = r.sequence().size();
            if (flag(r, 0x800) && (!r.tags().contains("SA"_tag) || !r.tags().contains("ch"_tag))) {
                invalid(r, "supplementary record lacks SA/ch");
            }

            size_t sameMate = 0;
            size_t representatives = 0;

            for (auto other : indices) {
                if (mate(records[other]) == m) {
                    ++sameMate;
                    if (!flag(records[other], 0x800)) {
                        ++representatives;
                    }
                }
            }
            if (representatives != 1) {
                invalid(r, "hypothesis requires one representative per mate");
            }

            if (sameMate > 1 && (!r.tags().contains("ch"_tag) || !r.tags().contains("SA"_tag))) {
                invalid(r, "ambiguous records sharing HI");
            }

            if (!hypothesis.paired && r.tags().contains("ch"_tag) && sameMate < 2) {
                invalid(r, "incomplete chimera");
            }

            if (r.tags().contains("SA"_tag)) {
                const auto& sa = r.tags().get<"SA"_tag>();
                size_t pos = 0;
                size_t links = 0;

                while (pos < sa.size()) {
                    const auto end = sa.find(';', pos);
                    if (end == std::string::npos) {
                        invalid(r, "unterminated SA");
                    }

                    const auto entry = sa.substr(pos, end - pos + 1);
                    size_t matches = 0;

                    for (const auto other : indices) {
                        if (other == index || mate(records[other]) != m ||
                            descriptor(records[other]) != entry) {
                            continue;
                        }
                        const auto& partner = records[other];

                        if (!partner.tags().contains("SA"_tag) ||
                            (";" + partner.tags().get<"SA"_tag>()).find(";" + descriptor(r)) ==
                                std::string::npos) {
                            invalid(r, "nonreciprocal SA");
                        }

                        ++matches;
                    }

                    if (matches != 1) {
                        invalid(r, "SA partner missing or ambiguous");
                    }

                    ++links;
                    pos = end + 1;
                }

                if (links != sameMate - 1) {
                    invalid(r, "incomplete SA relationships");
                }
            }

            for (auto& piece : pieces(r)) {
                piece.rawIndex = index;
                first[m] = std::min(first[m], piece.originalBegin);
                last[m] = std::max(last[m], piece.originalEnd);
                decoded.push_back(std::move(piece));
            }
        }
        for (size_t m = 0; m < 2; ++m) {
            std::vector<std::pair<size_t, size_t>> intervals;
            for (const auto& piece : decoded) {
                if (piece.mateIndex == m) {
                    intervals.emplace_back(piece.originalBegin, piece.originalEnd);
                }
            }
            std::ranges::sort(intervals);
            size_t end = 0;
            for (const auto& [begin, stop] : intervals) {
                if (stop > end) {
                    hypothesis.coveredBases[m] += stop - std::max(begin, end);
                }

                end = std::max(end, stop);
            }
        }
        for (auto& piece : decoded) {
            auto& fragment = piece.record;

            if (!excludeSoftClipping) {
                // Only the two outside ends can retain clipping; another arm's bases never can.
                const auto& raw = records[piece.rawIndex];
                const bool reverse = flag(raw, 0x10);
                const bool extendLeft = reverse ? piece.originalEnd == last[piece.mateIndex]
                                                : piece.originalBegin == first[piece.mateIndex];
                const bool extendRight = reverse ? piece.originalBegin == first[piece.mateIndex]
                                                 : piece.originalEnd == last[piece.mateIndex];
                const auto begin = extendLeft ? 0 : piece.begin;
                const auto end = extendRight ? raw.sequence().size() : piece.end;

                if (begin < piece.begin) {
                    fragment.cigar_sequence().insert(
                        fragment.cigar_sequence().begin(),
                        seqan3::cigar{static_cast<uint32_t>(piece.begin - begin),
                                      'S'_cigar_operation});
                }

                if (end > piece.end) {
                    fragment.cigar_sequence().push_back(
                        seqan3::cigar{static_cast<uint32_t>(end - piece.end), 'S'_cigar_operation});
                }

                fragment.sequence() =
                    seqan3::dna5_vector(raw.sequence().begin() + static_cast<long>(begin),
                                        raw.sequence().begin() + static_cast<long>(end));

                if (!raw.base_qualities().empty()) {
                    fragment.base_qualities() = std::vector<seqan3::phred42>(
                        raw.base_qualities().begin() + static_cast<long>(begin),
                        raw.base_qualities().begin() + static_cast<long>(end));
                }
            }
            fragment.tags()["XN"_tag] = static_cast<int32_t>(hypothesis.fragments.size());
            hypothesis.fragments.push_back(std::move(fragment));
        }
        result.push_back(std::move(hypothesis));
    }
    return result;
}
}  // namespace pipelines::align
