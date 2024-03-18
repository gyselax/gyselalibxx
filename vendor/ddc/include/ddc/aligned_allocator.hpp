// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <new>
#include <type_traits>

namespace ddc {

template <class T, std::size_t N>
class AlignedAllocator
{
public:
    using value_type = T;

    template <class U>
    struct rebind
    {
        using other = AlignedAllocator<U, N>;
    };

    constexpr AlignedAllocator() = default;

    constexpr AlignedAllocator(AlignedAllocator const& x) = default;

    constexpr AlignedAllocator(AlignedAllocator&& x) noexcept = default;

    template <class U>
    constexpr explicit AlignedAllocator(AlignedAllocator<U, N> const&) noexcept
    {
    }

    ~AlignedAllocator() = default;

    constexpr AlignedAllocator& operator=(AlignedAllocator const& x) = default;

    constexpr AlignedAllocator& operator=(AlignedAllocator&& x) noexcept = default;

    template <class U>
    constexpr AlignedAllocator& operator=(AlignedAllocator<U, N> const&) noexcept
    {
    }

    [[nodiscard]] T* allocate(std::size_t n) const
    {
        return new (std::align_val_t(N)) value_type[n];
    }

    void deallocate(T* p, std::size_t) const
    {
        operator delete[](p, std::align_val_t(N));
    }
};

template <class T, std::size_t NT, class U, std::size_t NU>
constexpr bool operator==(AlignedAllocator<T, NT> const&, AlignedAllocator<U, NU> const&) noexcept
{
    return std::is_same_v<AlignedAllocator<T, NT>, AlignedAllocator<U, NU>>;
}

template <class T, std::size_t NT, class U, std::size_t NU>
constexpr bool operator!=(AlignedAllocator<T, NT> const&, AlignedAllocator<U, NU> const&) noexcept
{
    return !std::is_same_v<AlignedAllocator<T, NT>, AlignedAllocator<U, NU>>;
}

} // namespace ddc
