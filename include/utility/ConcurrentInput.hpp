#pragma once

#include <condition_variable>
#include <cstddef>
#include <deque>
#include <exception>
#include <future>
#include <iterator>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

namespace utility {

// One reader owns the input iterator. Each consumer owns its popped record.
// Cancellation wakes both sides; producer exceptions are rethrown by consumers.
template <typename T>
class ConcurrentInput {
    std::mutex mutex;
    std::condition_variable ready;
    std::deque<T> queue;
    size_t capacity;
    bool finished{false};
    bool cancelled{false};
    std::exception_ptr error;
    std::thread producer;

    auto pop() -> std::optional<T> {
        std::unique_lock lock(mutex);
        ready.wait(lock, [&] { return cancelled || error || finished || !queue.empty(); });
        if (error) std::rethrow_exception(error);
        if (cancelled || queue.empty()) return std::nullopt;
        T value = std::move(queue.front());
        queue.pop_front();
        ready.notify_all();
        return value;
    }

   public:
    template <typename Reader>
    ConcurrentInput(size_t capacity, Reader reader) : capacity(capacity) {
        if (capacity == 0) throw std::invalid_argument("Input buffer size must be positive");
        producer = std::thread([this, reader = std::move(reader)]() mutable {
            try {
                while (true) {
                    {
                        std::unique_lock lock(mutex);
                        ready.wait(lock,
                                   [&] { return cancelled || queue.size() < this->capacity; });
                        if (cancelled) break;
                    }
                    auto value = reader();
                    if (!value) break;
                    std::lock_guard lock(mutex);
                    if (cancelled) break;
                    queue.push_back(std::move(*value));
                    ready.notify_all();
                }
            } catch (...) {
                std::lock_guard lock(mutex);
                error = std::current_exception();
            }
            {
                std::lock_guard lock(mutex);
                finished = true;
            }
            ready.notify_all();
        });
    }

    ConcurrentInput(const ConcurrentInput&) = delete;
    auto operator=(const ConcurrentInput&) -> ConcurrentInput& = delete;
    ~ConcurrentInput() {
        cancel();
        if (producer.joinable()) producer.join();
    }

    void cancel() {
        std::lock_guard lock(mutex);
        cancelled = true;
        ready.notify_all();
    }

    class iterator {
        ConcurrentInput* source{};
        mutable std::optional<T> value;

       public:
        using value_type = T;
        using difference_type = std::ptrdiff_t;
        using iterator_concept = std::input_iterator_tag;
        using iterator_category = std::input_iterator_tag;
        iterator() = default;
        explicit iterator(ConcurrentInput& source) : source(&source), value(source.pop()) {}
        auto operator*() const -> T& { return *value; }
        auto operator->() const -> T* { return &*value; }
        auto operator++() -> iterator& {
            value = source->pop();
            return *this;
        }
        void operator++(int) { ++*this; }
        auto operator==(std::default_sentinel_t) const -> bool { return !value; }
    };
    auto begin() -> iterator { return iterator{*this}; }
    auto end() const -> std::default_sentinel_t { return {}; }
};

// Always collect every worker before propagating failure. A failed consumer
// cancels the source so a producer blocked on a full queue cannot deadlock.
template <typename Input, typename Consumer>
auto consumeConcurrently(Input& input, size_t workers, Consumer consumer) {
    using Result = decltype(consumer());
    std::vector<std::future<Result>> futures;
    std::exception_ptr error;
    std::optional<Result> total;
    if (workers == 0) throw std::invalid_argument("At least one input consumer is required");
    try {
        for (size_t i = 0; i < workers; ++i) {
            futures.emplace_back(std::async(std::launch::async, [&] {
                try {
                    return consumer();
                } catch (...) {
                    input.cancel();
                    throw;
                }
            }));
        }
    } catch (...) {
        error = std::current_exception();
        input.cancel();
    }
    for (auto& future : futures) {
        try {
            auto result = future.get();
            if (total)
                *total += result;
            else
                total.emplace(std::move(result));
        } catch (...) {
            if (!error) error = std::current_exception();
            input.cancel();
        }
    }
    if (error) std::rethrow_exception(error);
    return std::move(*total);
}
}  // namespace utility
