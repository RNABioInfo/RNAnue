#include <gtest/gtest.h>

#include <algorithm>
#include <optional>
#include <random>
#include <span>
#include <string>
#include <tuple>
#include <vector>

#include "pipelines/preprocess/PreprocessFilter.hpp"
#include "pipelines/preprocess/RecordTrimmer.hpp"

using namespace pipelines::preprocess;
using dataTypes::FastqRecord;

namespace {
constexpr auto adapterText = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";

auto dna(const std::string& text) -> seqan3::dna5_vector {
    seqan3::dna5_vector result;
    for (char c : text) result.push_back(seqan3::assign_char_to(c, seqan3::dna5{}));
    return result;
}

auto record(const std::string& text) -> FastqRecord {
    FastqRecord result;
    result.sequence() = dna(text);
    result.id() = "probe";
    for (std::size_t i = 0; i < text.size(); ++i)
        result.base_qualities().push_back(seqan3::assign_phred_to(i % 42, seqan3::phred42{}));
    return result;
}

auto adapter(TrimConfig::Mode mode, double fraction = 0.05,
             const std::string& sequence = adapterText) -> Adapter {
    return {dna(sequence), fraction, mode};
}

void expectTrim(const std::string& input, const std::string& expected, TrimConfig::Mode mode,
                double fraction = 0.05, std::size_t minOverlap = 5) {
    auto r = record(input);
    const auto originalQualities = r.base_qualities();
    const auto a = adapter(mode, fraction);
    RecordTrimmer::trimAdapter(a, r, minOverlap);
    EXPECT_EQ(r.sequence(), dna(expected));
    ASSERT_EQ(r.base_qualities().size(), expected.size());
    const auto offset = mode == TrimConfig::Mode::FIVE_PRIME ? input.size() - expected.size() : 0;
    EXPECT_TRUE(std::equal(r.base_qualities().begin(), r.base_qualities().end(),
                           originalQualities.begin() + offset));
    RecordTrimmer::trimAdapter(a, r, minOverlap);
    EXPECT_EQ(r.sequence(), dna(expected));
}

// Intentionally exhaustive: no optimized score bounds, bands, batching, or
// ungapped assumptions. Every permitted fixed interval is scored by SeqAn3.
auto oracle(const Adapter& a, const seqan3::dna5_vector& read, std::size_t minimum)
    -> std::optional<std::pair<std::size_t, std::size_t>> {
    using View = std::span<const seqan3::dna5>;
    struct Candidate { View adapter; View read; std::size_t begin; };
    std::vector<Candidate> candidates;
    const View sequence{a.sequence};
    const View input{read};
    const auto minOverlap = std::max(std::size_t{1}, minimum);
    if (sequence.empty() || read.size() < minOverlap) return std::nullopt;
    for (std::size_t begin = 0; begin + minOverlap <= read.size(); ++begin) {
        for (std::size_t end = begin + minOverlap; end <= read.size(); ++end)
            candidates.push_back({sequence, input.subspan(begin, end - begin), begin});
        for (std::size_t length = 1; length < sequence.size(); ++length) {
            if (a.trimmingMode == TrimConfig::Mode::THREE_PRIME)
                candidates.push_back({sequence.first(length), input.subspan(begin), begin});
            else if (begin == 0)
                for (std::size_t end = minOverlap; end <= read.size(); ++end)
                    candidates.push_back({sequence.last(length), input.first(end), 0});
        }
    }
    const auto config = seqan3::align_cfg::method_global{} |
                        seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{
                            seqan3::match_score{1}, seqan3::mismatch_score{-1}}} |
                        seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-2},
                                                        seqan3::align_cfg::extension_score{-4}} |
                        seqan3::align_cfg::output_score{};
    auto pairs = candidates | std::views::transform([](const Candidate& c) {
                     return std::tuple{c.adapter, c.read};
                 });
    std::optional<std::pair<std::size_t, std::size_t>> best;
    int bestScore = 0;
    std::size_t index = 0;
    for (const auto& result : seqan3::align_pairwise(pairs, config)) {
        const auto& c = candidates[index++];
        const auto length = c.read.size();
        const int threshold = double(length) - ((length * a.maxMissMatchFraction) * 2);
        if (result.score() < threshold) continue;
        const auto end = c.begin + length;
        const bool outermost = best && (a.trimmingMode == TrimConfig::Mode::FIVE_PRIME
                                           ? c.begin > best->first : c.begin < best->first);
        if (!best || outermost || (c.begin == best->first &&
                                   (result.score() > bestScore ||
                                    (result.score() == bestScore && end > best->second)))) {
            best = {c.begin, end};
            bestScore = result.score();
        }
    }
    return best;
}
}  // namespace

TEST(AdapterTrimming, RemovesRepeatedAndSeparatedAdaptersAtBothEnds) {
    const std::string a{adapterText};
    const std::string insert{"TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGCTCTCACCG"};
    expectTrim(insert + a + a + a.substr(0, 9), insert, TrimConfig::Mode::THREE_PRIME);
    expectTrim(insert + a + "CCCCCCCC" + a, insert, TrimConfig::Mode::THREE_PRIME);
    expectTrim(a + a + insert, insert, TrimConfig::Mode::FIVE_PRIME);
    expectTrim(a + "CCCCCCCC" + a + insert, insert, TrimConfig::Mode::FIVE_PRIME);
}

TEST(AdapterTrimming, UsesPositionBeforeScore) {
    const std::string a{adapterText};
    auto imperfect = a;
    imperfect[15] = 'T';
    const std::string insert{"CCCCCCCCCCCCCCCCCCCC"};
    expectTrim(insert + imperfect + a, insert, TrimConfig::Mode::THREE_PRIME);
    expectTrim(a + imperfect + insert, insert, TrimConfig::Mode::FIVE_PRIME);
}

TEST(AdapterTrimming, AcceptsTerminalPartialsButProtectsInternalMotifs) {
    const std::string a{adapterText};
    const std::string insert{"CCCCCCCCAGATCCCCCCCCCCCCCCCCCC"};
    expectTrim(insert + a.substr(0, 12), insert, TrimConfig::Mode::THREE_PRIME);
    expectTrim(a.substr(a.size() - 12) + insert, insert, TrimConfig::Mode::FIVE_PRIME);
    expectTrim(insert, insert, TrimConfig::Mode::THREE_PRIME);
    expectTrim(insert, insert, TrimConfig::Mode::FIVE_PRIME);
    expectTrim(insert + a.substr(0, 4), insert + a.substr(0, 4), TrimConfig::Mode::THREE_PRIME);
}

TEST(AdapterTrimming, RejectsExcessiveErrorsAndHonorsOverlapAndZero) {
    const std::string a{adapterText};
    auto bad = a;
    bad[12] = 'T';
    bad[15] = 'T';
    bad[18] = bad[18] == 'T' ? 'A' : 'T';
    expectTrim("CCCCCCCC" + bad, "CCCCCCCC" + bad, TrimConfig::Mode::THREE_PRIME);
    expectTrim("CCCCCCCC" + a, "CCCCCCCC" + a, TrimConfig::Mode::THREE_PRIME, .05, 40);
    auto r = record("");
    RecordTrimmer::trimAdapter(adapter(TrimConfig::Mode::THREE_PRIME), r, 0);
    EXPECT_TRUE(r.sequence().empty());
    r = record("CCCCCCCC");
    RecordTrimmer::trimAdapter(adapter(TrimConfig::Mode::THREE_PRIME, .05, ""), r, 0);
    EXPECT_EQ(r.sequence(), dna("CCCCCCCC"));
}

TEST(AdapterTrimming, RemovesDimersAndRejectsEmptyRecordsWithZeroLengthThreshold) {
    const std::string a{adapterText};
    expectTrim(a + a, "", TrimConfig::Mode::THREE_PRIME);
    expectTrim(a + a, "", TrimConfig::Mode::FIVE_PRIME);
    auto r = record(a + a);
    RecordTrimmer::trimAdapter(adapter(TrimConfig::Mode::THREE_PRIME), r, 5);
    EXPECT_FALSE(PreprocessFilter::passes({0, 0}, r));
}

TEST(AdapterTrimming, TrimsBothEndsAndSupportsMultipleAdapters) {
    const std::string a{adapterText};
    const std::string b{"TGTAGATCTCGGTGGTCGCCGTATCATT"};
    auto r = record(a + a + "CCCCCCCCCCCCCCCCCCCC" + b + b);
    RecordTrimmer::trimAdapter(adapter(TrimConfig::Mode::FIVE_PRIME), r, 5);
    RecordTrimmer::trimAdapter(adapter(TrimConfig::Mode::THREE_PRIME, .05, b), r, 5);
    EXPECT_EQ(r.sequence(), dna("CCCCCCCCCCCCCCCCCCCC"));
    EXPECT_EQ(r.base_qualities().size(), r.sequence().size());
}

TEST(AdapterTrimming, RecoversAllSimulationProbeInsertsExactly) {
    const std::string source =
        "TCCCTGGTGGTCTAGTGGTTAGGATTCGGCGCTCTCACCGCCGCGGCCCGGGTTCGATTCCCGGTCAGGGAAAGA";
    const std::string a{adapterText};
    for (std::size_t length : {40, 60, 73, 75, 84, 90, 100, 120}) {
        SCOPED_TRACE(length);
        const auto insert = (source + source).substr(0, length);
        std::string input = insert;
        while (input.size() < 150) input += a;
        input.resize(150);
        expectTrim(input, insert, TrimConfig::Mode::THREE_PRIME);
    }
}

TEST(AdapterTrimming, MatchesExhaustiveSeqAnOracleAcrossConfigurations) {
    std::mt19937 random{271828};
    const std::string alphabet{"ACGTN"};
    auto randomText = [&](std::size_t length) {
        std::string text;
        while (text.size() < length) text += alphabet[random() % alphabet.size()];
        return text;
    };
    const std::array fractions{0.0, .05, .1, .2, .4, .5, 1.0};
    for (int trial = 0; trial < 210; ++trial) {
        const auto atext = randomText(4 + random() % 13);
        auto input = randomText(1 + random() % 18);
        if (trial % 3 == 0) input.insert(random() % (input.size() + 1), atext);
        if (trial % 3 == 1) {
            auto edited = atext;
            if (trial % 2) edited.erase(edited.begin() + random() % edited.size());
            else edited.insert(random() % edited.size(), "N");
            input += edited;
        }
        const auto overlap = std::size_t(trial % 7);
        const auto fraction = fractions[trial % fractions.size()];
        for (auto mode : {TrimConfig::Mode::FIVE_PRIME, TrimConfig::Mode::THREE_PRIME}) {
            SCOPED_TRACE(std::to_string(trial) + ": " + atext + " / " + input +
                         " / " + std::to_string(fraction) + " / " + std::to_string(int(mode)));
            const auto a = adapter(mode, fraction, atext);
            auto r = record(input);
            const auto q = r.base_qualities();
            const auto expected = oracle(a, r.sequence(), overlap);
            std::size_t begin = 0, end = input.size();
            if (expected) {
                if (mode == TrimConfig::Mode::FIVE_PRIME) begin = expected->second;
                else end = expected->first;
            }
            RecordTrimmer::trimAdapter(a, r, overlap);
            EXPECT_EQ(r.sequence(), dna(input.substr(begin, end - begin)));
            EXPECT_EQ(r.base_qualities(), (std::vector<seqan3::phred42>{q.begin() + begin,
                                                                      q.begin() + end}));
        }
    }
}

TEST(AdapterTrimming, GappedBatchBoundariesAndLongAdapterAgreeWithOracle) {
    const std::string a = std::string(adapterText) + adapterText;
    for (auto mode : {TrimConfig::Mode::FIVE_PRIME, TrimConfig::Mode::THREE_PRIME}) {
        for (const auto& input : {a.substr(0, 30) + "T" + a.substr(30) + "CCCCCC",
                                 std::string("CCCCCC") + a.substr(0, 30) + a.substr(31),
                                 std::string("CCCCCC") + a + a}) {
            auto r = record(input);
            const auto ad = adapter(mode, .15, a);
            const auto expected = oracle(ad, r.sequence(), 20);
            ASSERT_TRUE(expected);
            RecordTrimmer::trimAdapter(ad, r, 20);
            const auto retained = mode == TrimConfig::Mode::FIVE_PRIME
                                      ? input.substr(expected->second) : input.substr(0, expected->first);
            EXPECT_EQ(r.sequence(), dna(retained));
        }
    }
}
