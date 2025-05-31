#pragma once

// Standard
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <ostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

// Internal
#include "LogLevel.hpp"
#include "Logger.hpp"
#include "SamRecord.hpp"
#include "seqan3/alphabet/views/to_char.hpp"
#include "seqan3/core/debug_stream/default_printer.hpp"
#include "seqan3/utility/range/to.hpp"

namespace dataTypes {

enum class SplitRecordType : std::uint8_t { SINGLETON, CHIMERIC, MULTIMERIC };

inline auto to_string(SplitRecordType direction) -> std::string {
    switch (direction) {
        case SplitRecordType::SINGLETON:
            return "Singleton";
        case SplitRecordType::CHIMERIC:
            return "Chimeric";
        case SplitRecordType::MULTIMERIC:
            return "Multimeric";
        default:
            return "Other";
    }
}

inline auto operator<<(std::ostream& out, SplitRecordType type) -> std::ostream& {
    out << to_string(type);

    return out;
}

template <typename T>
concept RecordsContainer = requires(T value) {
    std::is_same_v<SamRecord, std::iter_value_t<typename T::iterator>>;
    value.begin();
    value.end();
    value.size();
    T::splitRecordType;
};

/**
 * @brief Represents a simple singleton read (not split).
 *
 * Provides a container–like API (begin, end, size, operator[]) so that it
 * behaves similarly to MultisplitRecords and SplitRecords.
 */
class SingletonRecord {
   public:
    using container = std::array<SamRecord, 1>;
    using iterator = container::iterator;
    using const_iterator = container::const_iterator;

    static constexpr auto splitRecordType{SplitRecordType::SINGLETON};

    /**
     * @brief Constructs a SingletonRecord from a single SamRecord.
     *
     * @param record The SamRecord to wrap.
     */
    explicit SingletonRecord(const SamRecord& record) : record{record} {}

    /**
     * @brief Constructs a SingletonRecord from a single SamRecord rvalue reference.
     *
     * @param record The SamRecord to wrap.
     */
    explicit SingletonRecord(SamRecord&& record) : record{std::move(record)} {}

    [[nodiscard]] auto getRecord() const noexcept -> const SamRecord& { return record.front(); }
    [[nodiscard]] auto getRecord() noexcept -> SamRecord { return record.front(); }

    [[nodiscard]] constexpr auto recordID() const noexcept -> const std::string& {
        return getRecord().id();
    }

    /**
     * @brief Returns an iterator to the beginning (the sole element).
     */
    [[nodiscard]] auto begin() noexcept -> iterator { return record.begin(); }
    [[nodiscard]] auto begin() const noexcept -> const_iterator { return record.begin(); }

    /**
     * @brief Returns an iterator to the end (one past the sole element).
     */
    [[nodiscard]] auto end() noexcept -> iterator { return record.end(); }
    [[nodiscard]] auto end() const noexcept -> const_iterator { return record.end(); }

    /**
     * @brief Returns the number of contained records (always 1).
     */
    [[nodiscard]] constexpr auto size() const noexcept -> std::size_t { return record.size(); }

    /**
     * @brief Returns the first (and only) element.
     */
    [[nodiscard]] auto front() noexcept -> SamRecord& { return record.front(); }
    [[nodiscard]] auto front() const noexcept -> const SamRecord& { return record.front(); }

    /**
     * @brief Returns the last (and only) element.
     */
    [[nodiscard]] auto back() noexcept -> SamRecord& { return record.back(); }
    [[nodiscard]] auto back() const noexcept -> const SamRecord& { return record.back(); }

    /**
     * @brief Less-than operator for comparing two SingletonRecord objects.
     *
     * The comparison is done on the underlying SamRecord taken from back().
     * First the reference ids are compared; if they differ then that result is returned.
     * Otherwise, the comparison of the end positions (using recordEndPosition) is returned.
     *
     * @param rhs The SingletonRecord to compare with.
     * @return true if this object's record is less than rhs's record; false otherwise.
     */
    [[nodiscard]] auto operator<(const SingletonRecord& rhs) const -> bool {
        const SamRecord& lhsRecord = record.back();
        const SamRecord& rhsRecord = rhs.record.back();
        if (lhsRecord.reference_id() != rhsRecord.reference_id()) {
            return lhsRecord.reference_id() < rhsRecord.reference_id();
        }
        return recordEndPosition(lhsRecord) < recordEndPosition(rhsRecord);
    }

    /**
     * @brief Greater-than operator.
     *
     * @param rhs The SingletonRecord to compare with.
     * @return true if this object's record is greater than rhs's record.
     */
    [[nodiscard]] auto operator>(const SingletonRecord& rhs) const -> bool { return rhs < *this; }

   private:
    container record;
};

inline auto operator<<(std::ostream& out, const SingletonRecord& record) -> std::ostream& {
    const auto printer = seqan3::default_printer();

    out << "<(";
    printer(out, (record.getRecord().sequence() | seqan3::views::to_char));
    out << ", ";
    printer(out, record.getRecord().cigar_sequence());
    out << ")>";

    return out;
}

static_assert(RecordsContainer<SingletonRecord>);

/**
 * @brief Represents a chimeric read that is split into two records and stored in a sorted fashion.
 *
 * The two records are rearranged on construction so that the order is
 * determined by dataTypes::operator<. This ordering guarantees that the "second"
 * record (accessed via back() or operator[](1)) is always greater or equal to
 * the first, which can later be used for comparisons.
 */
class ChimericRecords {
   public:
    /// The underlying container for two SamRecord objects.
    using container = std::array<SamRecord, 2>;
    using iterator = container::iterator;
    using const_iterator = container::const_iterator;

    static constexpr auto splitRecordType{SplitRecordType::CHIMERIC};

    /**
     * @brief Constructs a SplitRecords object from two SamRecord objects.
     *
     * The records are sorted on construction so that the lesser record (using dataTypes::operator<)
     * is stored in position 0 and the larger in position 1.
     *
     * @param rec1 The first record.
     * @param rec2 The second record.
     */
    ChimericRecords(const SamRecord& rec1, const SamRecord& rec2) : records{rec1, rec2} {
        sortRecords();
    }

    /**
     * @brief Constructs a SplitRecords object from two moved SamRecord objects.
     *
     * The records are sorted on construction so that the lesser record (using dataTypes::operator<)
     * is stored in position 0 and the larger in position 1.
     *
     * @param rec1 The first record.
     * @param rec2 The second record.
     */
    ChimericRecords(SamRecord&& rec1, SamRecord&& rec2)
        : records{std::move(rec1), std::move(rec2)} {
        sortRecords();
    }

    /**
     * @brief Constructs a SplitRecords object from a std::pair of SamRecord objects.
     *
     * The records are sorted on construction.
     *
     * @param recPair A std::pair of SamRecord objects.
     */
    explicit ChimericRecords(const std::pair<SamRecord, SamRecord>& recPair)
        : records{recPair.first, recPair.second} {
        sortRecords();
    }

    /**
     * @brief Returns an iterator to the first element.
     */
    [[nodiscard]] auto begin() noexcept -> iterator { return records.begin(); }
    [[nodiscard]] auto begin() const noexcept -> const_iterator { return records.begin(); }

    /**
     * @brief Returns an iterator pointing past the last element.
     */
    [[nodiscard]] auto end() noexcept -> iterator { return records.end(); }
    [[nodiscard]] auto end() const noexcept -> const_iterator { return records.end(); }

    /**
     * @brief Returns the first element.
     */
    [[nodiscard]] auto first() noexcept -> SamRecord& { return records.front(); }
    [[nodiscard]] auto first() const noexcept -> const SamRecord& { return records.front(); }

    /**
     * @brief Returns the second element.
     */
    [[nodiscard]] auto second() noexcept -> SamRecord& { return records.back(); }
    [[nodiscard]] auto second() const noexcept -> const SamRecord& { return records.back(); }

    /**
     * @brief Returns the number of records (always 2).
     */
    [[nodiscard]] constexpr auto size() const noexcept -> std::size_t { return records.size(); }

    /**
     * @brief Returns the last record (i.e. the one that compares greater).
     *
     * This is useful when comparisons and sorting use the "end" of the record pair.
     */
    [[nodiscard]] constexpr auto front() noexcept -> SamRecord& { return records.front(); }
    [[nodiscard]] constexpr auto front() const noexcept -> const SamRecord& {
        return records.front();
    }

    /**
     * @brief Returns the last record (i.e. the one that compares greater).
     *
     * This is useful when comparisons and sorting use the "end" of the record pair.
     */
    [[nodiscard]] constexpr auto back() noexcept -> SamRecord& { return records.back(); }
    [[nodiscard]] constexpr auto back() const noexcept -> const SamRecord& {
        return records.back();
    }

    [[nodiscard]] constexpr auto recordID() const noexcept -> const std::string& {
        return front().id();
    }

    /**
     * @brief Less-than operator for comparing two SplitRecords objects.
     *
     * The comparison is done by comparing the "last" record of each pair. First
     * the reference ids are compared; if they differ then that result is returned.
     * Otherwise, the comparison of the end positions (using recordEndPosition) is returned.
     * This means the last base (3' direction) of all records determines the larger SplitRecords.
     *
     * Precondition: Both SplitRecords objects have been properly constructed.
     *
     * @param rhs The SplitRecords to compare with.
     * @return true if this object is less than rhs; false otherwise.
     */
    [[nodiscard]] auto operator<(const ChimericRecords& rhs) const -> bool {
        const SamRecord& lhs_last = this->back();
        const SamRecord& rhs_last = rhs.back();

        if (lhs_last.reference_id() != rhs_last.reference_id()) {
            return lhs_last.reference_id() < rhs_last.reference_id();
        }
        return recordEndPosition(lhs_last) < recordEndPosition(rhs_last);
    }

    /**
     * @brief Greater-than operator.
     *
     * @param rhs The SplitRecords to compare with.
     * @return true if this object is greater than rhs; false otherwise.
     */
    [[nodiscard]] auto operator>(const ChimericRecords& rhs) const -> bool { return rhs < *this; }

   private:
    container records;

    /**
     * @brief Sorts the two SamRecord objects.
     *
     * Uses dataTypes::operator< to compare the two records.
     */
    void sortRecords() noexcept {
        if (dataTypes::operator>(records[0], records[1])) {
            std::swap(records[0], records[1]);
        }
    }
};

inline auto operator<<(std::ostream& out, const ChimericRecords& record) -> std::ostream& {
    const auto printer = seqan3::default_printer();

    out << "<";
    for (const auto& record : record) {
        out << "(";
        printer(out, (record.sequence() | seqan3::views::to_char));
        out << ", ";
        printer(out, record.cigar_sequence());
        out << ")";
    }
    out << ">";

    return out;
}

static_assert(RecordsContainer<ChimericRecords>);

/**
 * @brief Represents a multimeric read which has been split into multiple SAM records.
 *
 * This type manages a collection of SamRecord instances. The records are automatically
 * sorted (using std::ranges::sort) upon construction, and the comparison operator (<)
 * uses the last element’s reference id and end position.
 */
class MultimericRecords {
   public:
    using container = std::vector<SamRecord>;
    using iterator = container::iterator;
    using const_iterator = container::const_iterator;

    static constexpr auto splitRecordType{SplitRecordType::MULTIMERIC};

    /**
     * @brief Construct from an existing collection of SamRecord objects.
     *
     * @param records The records to store. They will be copied and then sorted.
     */
    explicit MultimericRecords(const container& records) : records{records} {
        // TODO: Change this to warning if multisplits is implemented
        if (records.size() < 3) {
            Logger::log<LogLevel::DEBUG>(
                "Multisplit records are required to contain more then two records. Record ID: ",
                records.empty() ? "No records" : records.front().id());
        }
        sortRecords();
    }

    /**
     * @brief Construct from an rvalue reference to a collection of SamRecord objects.
     *
     * @param records The records to store. They will be moved and then sorted.
     */
    explicit MultimericRecords(container&& records) : records{std::move(records)} {
        if (records.size() < 3) {
            Logger::log<LogLevel::WARNING>(
                "Multimeric records are required to contain more then two records. Record ID: ",
                records.empty() ? "No records" : records.front().id());
        }
        sortRecords();
    }

    // --- Element access and container interface ---

    /**
     * @brief Returns an iterator to the first element.
     */
    [[nodiscard]] auto begin() noexcept -> iterator { return records.begin(); }
    [[nodiscard]] auto begin() const noexcept -> const_iterator { return records.begin(); }

    /**
     * @brief Returns an iterator to the element following the last element.
     */
    [[nodiscard]] auto end() noexcept -> iterator { return records.end(); }
    [[nodiscard]] auto end() const noexcept -> const_iterator { return records.end(); }

    /**
     * @brief Returns the last record (i.e. the one that compares greater).
     *
     * This is useful when comparisons and sorting use the "end" of the record pair.
     */
    [[nodiscard]] auto front() noexcept -> SamRecord& { return records.front(); }
    [[nodiscard]] auto front() const noexcept -> const SamRecord& { return records.front(); }

    /**
     * @brief Returns the last record (i.e. the one that compares greater).
     *
     * This is useful when comparisons and sorting use the "end" of the record pair.
     */
    [[nodiscard]] auto back() noexcept -> SamRecord& { return records.back(); }
    [[nodiscard]] auto back() const noexcept -> const SamRecord& { return records.back(); }

    /**
     * @brief Returns the number of contained SamRecord objects.
     */
    [[nodiscard]] auto size() const noexcept -> std::size_t { return records.size(); }

    /**
     * @brief Checks whether the container is empty.
     */
    [[nodiscard]] auto empty() const noexcept -> bool { return records.empty(); }

    /**
     * @brief Element access: returns a reference to the element at the specified location.
     *
     * @param pos Position of the element to return.
     */
    [[nodiscard]] auto operator[](std::size_t pos) noexcept -> SamRecord& { return records[pos]; }
    [[nodiscard]] auto operator[](std::size_t pos) const noexcept -> const SamRecord& {
        return records[pos];
    }

    /**
     * @brief Less-than operator for comparing two MultisplitRecords objects.
     *
     * The comparison is performed using the last record in the container.
     * First the record’s reference id is compared; if those are equal, then
     * the (external) function recordEndPosition is used.
     * This means the last base (3' direction) of all records determines the larger SplitRecords.
     *
     * Precondition: There is at least one record in the container.
     *
     * @param rhs The right-hand side of the comparison.
     * @return true if *this is less than rhs, false otherwise.
     */
    [[nodiscard]] auto operator<(const MultimericRecords& rhs) const -> bool {
        const SamRecord& lhs_last = records.back();
        const SamRecord& rhs_last = rhs.records.back();

        if (lhs_last.reference_id() != rhs_last.reference_id()) {
            return lhs_last.reference_id() < rhs_last.reference_id();
        }
        return recordEndPosition(lhs_last) < recordEndPosition(rhs_last);
    }

    /**
     * @brief Greater-than operator.
     *
     * @param rhs The right-hand side of the comparison.
     * @return true if *this is greater than rhs, false otherwise.
     */
    [[nodiscard]] auto operator>(const MultimericRecords& rhs) const -> bool { return rhs < *this; }

    [[nodiscard]] constexpr auto recordID() const noexcept -> const std::string& {
        return front().id();
    }

   private:
    container records;

    /**
     * @brief Sorts the contained SamRecord objects.
     *
     * Skips sorting if there are 0 or 1 records.
     */
    void sortRecords() {
        if (records.size() > 1) {
            std::ranges::sort(records, dataTypes::operator<);
        }
    }
};

inline auto operator<<(std::ostream& out, const MultimericRecords& record) -> std::ostream& {
    const auto printer = seqan3::default_printer();

    out << "<";
    for (const auto& record : record) {
        out << "(";
        printer(out, (record.sequence() | seqan3::views::to_char));
        out << ", ";
        printer(out, record.cigar_sequence());
        out << ")";
    }
    out << ">";

    return out;
}

static_assert(RecordsContainer<MultimericRecords>);

template <typename T>
concept SplitRecordsContainer = RecordsContainer<T> && (std::is_same_v<T, ChimericRecords> ||
                                                        std::is_same_v<T, MultimericRecords>);

}  // namespace dataTypes
