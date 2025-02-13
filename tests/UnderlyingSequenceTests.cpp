// NOLINTBEGIN
#include <gtest/gtest.h>

// Standard
#include <cstddef>
#include <iterator>
#include <ranges>
#include <string>
#include <vector>

// seqan3
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sam_file/sam_flag.hpp>

// UnderlyingSequence (iterator and view) under test.
#include "UnderlyingSequence.hpp"

// Convenience: using directive for seqan3 nucleotides.
using seqan3::dna5;
using seqan3::sam_flag;
using seqan3::operator""_dna5;

// Helper function to convert a range of dna5 elements to a std::string (using to_char).
auto to_string(std::ranges::range auto const& range) -> std::string {
    std::string ret;
    for (auto nuc : range) {
        ret.push_back(seqan3::to_char(nuc));
    }
    return ret;
}

// Given a vector of nucleotides and a flag, returns a std::vector<char> with the iterated letters.
auto collect_letters(auto&& view) -> std::vector<char> {
    std::vector<char> letters;
    for (auto nuc : view) {
        letters.push_back(seqan3::to_char(nuc));
    }
    return letters;
}

TEST(UnderlyingSequenceTest, ForwardIteration) {
    // Create an example nucleotide sequence.
    std::vector<dna5> seq = {"ATGC"_dna5};

    // Create UnderlyingSequenceView in forward mode (i.e. SAM flag not set).
    UnderlyingSequenceView seq_view(seq, sam_flag{});

    // Test that begin() dereferences to the first element.
    auto iter = seq_view.begin();
    EXPECT_EQ(seqan3::to_char(*iter), 'A');

    // Collect the letters via range-based for loop.
    std::vector<char> letters = collect_letters(seq_view);
    // In forward mode, the iterator returns the underlying value.
    std::vector<char> expected = {'A', 'T', 'G', 'C'};
    EXPECT_EQ(letters, expected);

    // Test size() and that cbegin() equals begin().
    EXPECT_EQ(seq_view.size(), seq.size());
    auto cit = seq_view.cbegin();
    EXPECT_EQ(seqan3::to_char(*cit), 'A');
}

TEST(UnderlyingSequenceTest, ReverseIteration) {
    // Create a nucleotide sequence that is not self-complementary.
    // For instance "A", "C", "G", "T" will be reversed to "T", "G", "C", "A"
    // but then complemented becomes: complement(T)=A, complement(G)=C, complement(C)=G,
    // complement(A)=T.
    std::vector<dna5> seq = {'A'_dna5, 'C'_dna5, 'G'_dna5, 'T'_dna5};

    // Create UnderlyingSequenceView in reverse mode.
    UnderlyingSequenceView rev_view(seq, sam_flag::on_reverse_strand);

    // Test that begin() dereferences correctly:
    auto iter = rev_view.begin();
    // In reverse mode, first element is complement(last): last in seq is 'T' -> complement('T') =
    // 'A'.
    EXPECT_EQ(seqan3::to_char(*iter), 'A');

    // Advance the iterator and check subsequent elements:
    ++iter;
    // Second element: should be complement(seq[2]): 'G' -> complement('G') = 'C'
    EXPECT_EQ(seqan3::to_char(*iter), 'C');
    ++iter;
    EXPECT_EQ(seqan3::to_char(*iter), 'G');  // complement(seq[1]): 'C' -> 'G'
    ++iter;
    EXPECT_EQ(seqan3::to_char(*iter), 'T');  // complement(seq[0]): 'A' -> 'T'
    ++iter;
    // Ensure that iteration stops at the reverse end.
    EXPECT_EQ(iter, rev_view.end());

    // Collect the letters via range-based for loop.
    std::vector<char> letters = collect_letters(rev_view);
    // Expected reverse view sequence: "A", "C", "G", "T"
    std::vector<char> expected = {'A', 'C', 'G', 'T'};
    EXPECT_EQ(letters, expected);

    // Verify size is unchanged.
    EXPECT_EQ(rev_view.size(), seq.size());
}

TEST(UnderlyingSequenceTest, IteratorEqualityAndRandomAccess) {
    // This test exercises the iterator’s equality and arithmetic operators.
    std::vector<dna5> seq = {'A'_dna5, 'T'_dna5, 'G'_dna5, 'C'_dna5, 'G'_dna5};

    // Forward view.
    UnderlyingSequenceView seq_view(seq, sam_flag{});
    auto begin_it = seq_view.begin();
    auto end_it = seq_view.end();

    // Test operator++ and equality.
    auto iter = begin_it;
    size_t count = 0;
    while (iter != end_it) {
        ++iter;
        ++count;
    }
    EXPECT_EQ(count, seq.size());

    // Test random access operations on the iterator. (Note: our original iterator is only forward.
    // If you later adopt the random access version, then these tests can check operator+ etc.)
    // For now, use std::next and operator==.
    auto it2 = std::next(begin_it, 3);
    // In forward mode, character at position 3 should be 'C'
    EXPECT_EQ(seqan3::to_char(*it2), 'C');
}

TEST(UnderlyingSequenceTest, ConstViewAndIterator) {
    // Create a const sequence and test const view iteration.
    const std::vector<dna5> seq = {'T'_dna5, 'G'_dna5, 'C'_dna5, 'A'_dna5};

    // Create a const UnderlyingSequenceView in forward mode.
    const UnderlyingSequenceView const_view(seq, sam_flag{});
    static_cast<void>(const_view.begin());
    std::vector<char> letters;
    for (auto nuc : const_view) {
        letters.push_back(seqan3::to_char(nuc));
    }
    std::vector<char> expected = {'T', 'G', 'C', 'A'};
    EXPECT_EQ(letters, expected);
}

TEST(UnderlyingSequenceTest, AdaptorFromFunctor) {
    // Test that the view adaptor in the views namespace works as expected.
    std::vector<dna5> seq = {'A'_dna5, 'T'_dna5, 'G'_dna5, 'C'_dna5};

    // Use the underlying_sequence adaptor (in forward mode).
    auto view0 = views::underlying_sequence(seq, sam_flag{});
    std::vector<char> letters_forward = collect_letters(view0);
    std::vector<char> expected_forward = {'A', 'T', 'G', 'C'};
    EXPECT_EQ(letters_forward, expected_forward);

    // And test with reverse flag.
    auto view1 = views::underlying_sequence(seq, sam_flag::on_reverse_strand);
    // For our sequence, reverse view yields: first element = complement(last) = complement('C') =
    // 'G', then complement('G') = 'C', then complement('T') = 'A', then complement('A') = 'T'
    std::vector<char> letters_reverse = collect_letters(view1);
    std::vector<char> expected_reverse = {'G', 'C', 'A', 'T'};
    EXPECT_EQ(letters_reverse, expected_reverse);
}

TEST(UnderlyingSequenceTest, MultiplePassRange) {
    // Ensure that the UnderlyingSequenceView can be iterated multiple times.
    std::vector<dna5> seq = {'A'_dna5, 'C'_dna5, 'G'_dna5, 'T'_dna5};
    UnderlyingSequenceView view(seq, sam_flag{});

    std::string first_pass = to_string(view);
    std::string second_pass = to_string(view);

    EXPECT_EQ(first_pass, second_pass);
    EXPECT_EQ(first_pass, "ACGT");
}

TEST(UnderlyingSequenceTest, ModifyUnderlyingSequence) {
    // Verify that when a sequence is passed by value into the view constructor,
    // modifications to the original do not affect the view.
    std::vector<dna5> seq = {'A'_dna5, 'T'_dna5, 'G'_dna5, 'C'_dna5};
    UnderlyingSequenceView view(seq, sam_flag{});
    seq[0] = 'C'_dna5;  // modify original
    // The view holds its own copy (if using std::views::all on an lvalue,
    // it usually stores a reference, so be careful; adjust this test if your design differs).
    std::vector<char> letters = collect_letters(view);
    // In this design, we store the underlying range via std::views::all, so modifications
    // to seq might be seen. For this test assume that we want the view to reflect changes.
    EXPECT_EQ(letters[0], 'C');
}
// NOLINTEND
