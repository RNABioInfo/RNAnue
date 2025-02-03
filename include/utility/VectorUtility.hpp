#pragma once

// Standard
#include <algorithm>
#include <cstddef>
#include <vector>

namespace helper {

/**
 * @brief Erases selected elements from the vector.
 *
 * Removes elements from the provided vector `vec` at the indices specified in the sorted
 * `selection` vector. After removal, the vector is resized, and the new elements beyond
 * the resized length are initialized with `replaceVal`. It is expected that the `selection`
 * vector is sorted in ascending order to ensure correct behavior.
 *
 * @tparam T The type of elements stored in the vector.
 * @param vec The vector from which elements will be erased.
 * @param selection A sorted vector of indices indicating which elements to erase from `vec`.
 * @param replaceVal The value to assign to the elements beyond the new size of the vector after
 * resizing.
 */
template <class T>
inline void erase_selected(std::vector<T>& vec, const std::vector<size_t>& selection,
                           const T& replaceVal) {
    vec.resize(std::distance(vec.begin(),
                             std::stable_partition(
                                 vec.begin(), vec.end(),
                                 [&selection, &vec](const T& item) {
                                     return !std::ranges::binary_search(
                                         selection, static_cast<int>(static_cast<const T*>(&item) -
                                                                     vec.data()));
                                 })),
               replaceVal);
}

}  // namespace helper
