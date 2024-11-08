#pragma once

// seqan3
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

using seqan3::operator""_tag;

template <>
struct seqan3::sam_tag_type<"XX"_tag> {
    using type = int32_t;
};
template <>
struct seqan3::sam_tag_type<"XY"_tag> {
    using type = int32_t;
};
template <>
struct seqan3::sam_tag_type<"XJ"_tag> {
    using type = int32_t;
};
template <>
struct seqan3::sam_tag_type<"XH"_tag> {
    using type = int32_t;
};

template <>
struct seqan3::sam_tag_type<"XM"_tag> {
    using type = int32_t;
};  // matches in alignment
template <>
struct seqan3::sam_tag_type<"XL"_tag> {
    using type = int32_t;
};  // length of alignment

template <>
struct seqan3::sam_tag_type<"XN"_tag> {
    using type = float;
};

template <>
struct seqan3::sam_tag_type<"XS"_tag> {
    using type = std::string;
};

template <>
struct seqan3::sam_tag_type<"XC"_tag> {
    using type = float;
};

template <>
struct seqan3::sam_tag_type<"XR"_tag> {
    using type = float;
};

template <>
struct seqan3::sam_tag_type<"XE"_tag> {
    using type = float;
};

template <>
struct seqan3::sam_tag_type<"XO"_tag> {
    using type = int32_t;
};
