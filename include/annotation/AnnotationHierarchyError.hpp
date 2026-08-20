#pragma once

#include <stdexcept>
#include <string>
#include <utility>

namespace annotation {

class AnnotationHierarchyError final : public std::runtime_error {
   public:
    explicit AnnotationHierarchyError(std::string message)
        : std::runtime_error(std::move(message)) {}
};

}  // namespace annotation
