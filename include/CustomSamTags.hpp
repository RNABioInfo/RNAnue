#pragma once

// seqan3
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
#include <string>

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
};  // complementarity score

template <>
struct seqan3::sam_tag_type<"XC"_tag> {
    using type = float;
};  // complementarity

template <>
struct seqan3::sam_tag_type<"XR"_tag> {
    using type = float;
};  // complementarity fragment fraction

template <>
struct seqan3::sam_tag_type<"XE"_tag> {
    using type = float;
};  // hybridization energy

template <>
struct seqan3::sam_tag_type<"XO"_tag> {
    using type = int32_t;
};  // number of crosslinking sites

template <>
struct seqan3::sam_tag_type<"XA"_tag> {
    using type = std::string;
};  // Intra-molecular crosslinking sites

template <>
struct seqan3::sam_tag_type<"XI"_tag> {
    using type = std::string;
};  // Inter-molecular crosslinking sites

template <>
struct seqan3::sam_tag_type<"XD"_tag> {
    using type = std::string;
};  // Secondary structure
