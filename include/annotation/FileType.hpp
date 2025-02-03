#pragma once

// Standard
#include <cstdint>
#include <string>

namespace annotation {

class FileType {
   public:
    enum Value : uint8_t { GFF, GTF };
    explicit constexpr FileType(Value p_value) : value(p_value) {};

    [[nodiscard]] static constexpr auto attributeDelimiter() -> char { return ';'; }

    [[nodiscard]] constexpr auto attributeAssignment() const -> char {
        switch (value) {
            case GFF:
                return '=';
            default:
                return ' ';
        }
    }

    [[nodiscard]] constexpr auto defaultIDKey() const -> std::string {
        switch (value) {
            case GFF:
                return "ID";
            default:
                return "gene_id";
        }
    }

    [[nodiscard]] constexpr auto defaultGroupKey() const -> std::string {
        switch (value) {
            case GFF:
                return "Parent";
            default:
                return "transcript_id";
        }
    }

    [[nodiscard]] static constexpr auto defaultGeneNameKey() -> std::string { return "gene"; }

    constexpr operator Value() const { return value; }
    explicit operator bool() const = delete;

   private:
    Value value;
};

}  // namespace annotation
