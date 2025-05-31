#pragma once
#include <cassert>
#include <coroutine>
#include <cstddef>
#include <exception>
#include <iostream>
#include <memory_resource>
#include <new>
#include <optional>

// the default pmr memory resource is a thread-safe allocator that uses the global new and delete
inline static std::pmr::synchronized_pool_resource mem_pool{
    std::pmr::new_delete_resource()};  // default allocator

inline static thread_local std::pmr::memory_resource* pmem_pool =
    &mem_pool;  // never owns any supplied mem resource
inline static void set_pmr_mem_pool(std::pmr::memory_resource* mem_pool_cust) {
    pmem_pool = mem_pool_cust;  // set a custom pmr allocator (but does not take ownership)
}
inline static void reset_default_pmr_mem_pool() {
    pmem_pool = &mem_pool;  // reset using the default allocator (does not take ownership)
}

/**
 * General purpose C++20 coroutine generator template class.
 *
 * @tparam T the type of value that the generator returns to the caller
 */
template <typename T>
class [[nodiscard]] Generator {
   public:
    struct promise_type;
    using coro_handle_type = std::coroutine_handle<promise_type>;

   private:
    coro_handle_type coro;

   public:
    explicit Generator(coro_handle_type handle) : coro{handle} {}
    Generator(const Generator&) = delete;                     // do not allow copy construction
    auto operator=(const Generator&) -> Generator& = delete;  // do not allow copy assignment
    Generator(Generator&& oth) noexcept : coro{std::move(oth.coro)} {
        oth.coro = nullptr;  // insure the other moved handle is null
    }
    auto operator=(Generator&& other) noexcept -> Generator& {
        if (this != &other) {  // ignore assignment to self
            if (coro) {        // destroy self current handle
                coro.destroy();
            }
            coro = std::move(other.coro);  // move other coro handle into self
            other.coro = nullptr;          // insure other moved handle is null
        }
        return *this;
    }
    ~Generator() {
        if (coro) {
            coro.destroy();
            coro = nullptr;
        }
    }

    // API
    [[nodiscard]] auto next() const -> bool {
        if (!coro || coro.done()) {
            return false;  // nothing more to process
        }
        coro.resume();
        return !coro.done();
    }

    auto getValue() noexcept -> std::optional<T> {
        return coro ? std::make_optional(coro.promise().current_value) : std::nullopt;
    }

    // implementation of above opaque declaration promise_type
    struct promise_type {
       public:
        auto operator new(std::size_t size) -> void* {
            assert(pmem_pool != nullptr);
            return pmem_pool->allocate(size);
        }
        void operator delete(void* ptr, std::size_t size) {
            assert(pmem_pool != nullptr);
            pmem_pool->deallocate(ptr, size);
        }

       private:
        T current_value;
        friend class Generator;

       public:
        promise_type() = default;
        ~promise_type() = default;
        promise_type(const promise_type&) = delete;
        promise_type(promise_type&&) = delete;
        auto operator=(const promise_type&) -> promise_type& = delete;
        auto operator=(promise_type&&) -> promise_type& = delete;

        auto get_return_object() { return Generator{coro_handle_type::from_promise(*this)}; }

        auto initial_suspend() { return std::suspend_always{}; }

        auto final_suspend() noexcept { return std::suspend_always{}; }

        void return_void() {}

        auto yield_value(T some_value) {
            current_value = some_value;
            return std::suspend_always{};
        }

        void unhandled_exception() { std::terminate(); }
    };

    struct iterator {
        using difference_type [[maybe_unused]] = std::ptrdiff_t;
        using value_type [[maybe_unused]] = T;
        coro_handle_type hdl = nullptr;
        iterator() = default;
        iterator(coro_handle_type handle) : hdl{handle} {}
        void getNext() {
            if (hdl) {
                hdl.resume();
                if (hdl.done()) {
                    hdl = nullptr;
                }
            }
        }
        auto operator*() const -> T& {
            assert(hdl);
            return hdl.promise().current_value;
        }
        auto operator++() -> iterator& {  // pre-incrementable
            getNext();
            return *this;
        }
        void operator++(int) {  // post-incrementable
            ++*this;
        }
        auto operator==(const iterator& iter) const -> bool = default;
    };

    auto begin() const -> iterator {
        if (!coro || coro.done()) {
            return iterator{nullptr};
        }
        iterator itr{coro};
        itr.getNext();
        return itr;
    }

    auto end() const -> iterator { return iterator{nullptr}; }
};

// /**
//  * Helper class for establishing a pmr memory_resource compliant
//  * allocator that allocates from a supplied fixed size buffer, e.g.,
//  * such as a buffer allocated on the stack. Allocations are not
//  * freed so depends on buffer context being reclaimed when going
//  * out of scope (this class does not take ownership of the buffer).
//  */
// class fixed_buffer_pmr_allocator : public std::pmr::memory_resource {
//    private:
//     void* const buf;
//     size_t buf_size;

//    public:
//     const size_t max_buf_size;
//     fixed_buffer_pmr_allocator(void* buf, size_t buf_size)
//         : buf(buf), buf_size(buf_size), max_buf_size(buf_size) {}
//     fixed_buffer_pmr_allocator() = delete;
//     fixed_buffer_pmr_allocator(const fixed_buffer_pmr_allocator&) = delete;
//     fixed_buffer_pmr_allocator(fixed_buffer_pmr_allocator&&) = delete;
//     auto operator=(const fixed_buffer_pmr_allocator&) -> fixed_buffer_pmr_allocator& = delete;
//     auto operator=(fixed_buffer_pmr_allocator&&) -> fixed_buffer_pmr_allocator& = delete;

//    private:
//     auto do_allocate(size_t bytes, size_t alignment) -> void* override {
//         if (bytes > buf_size) {
//             std::cerr << "requested bytes: " << bytes << ", remaining bytes capacity: " <<
//             buf_size
//                       << '\n';
//             throw std::bad_alloc();
//         }
//         buf_size -= bytes;
//         return buf;
//     }
//     void do_deallocate(void* p, size_t bytes, size_t alignment) override {}
//     auto do_is_equal(const memory_resource& _other) const noexcept -> bool override {
//         return false;
//     }
// };
