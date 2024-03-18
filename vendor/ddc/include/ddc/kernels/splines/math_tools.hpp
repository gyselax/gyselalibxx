// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cmath>

#include <experimental/mdspan>

namespace ddc::detail {
template <typename T>
KOKKOS_INLINE_FUNCTION T sum(T* array, int size)
{
    T val(0.0);
    for (int i(0); i < size; ++i) {
        val += array[i];
    }
    return val;
}

template <class ElementType, class LayoutPolicy, class AccessorPolicy, std::size_t Ext>
KOKKOS_INLINE_FUNCTION ElementType sum(std::experimental::mdspan<
                                       ElementType,
                                       std::experimental::extents<std::size_t, Ext>,
                                       LayoutPolicy,
                                       AccessorPolicy> const& array)
{
    ElementType val(0.0);
    for (std::size_t i(0); i < array.extent(0); ++i) {
        val += array[i];
    }
    return val;
}

template <class ElementType, class LayoutPolicy, class AccessorPolicy, std::size_t Ext>
KOKKOS_INLINE_FUNCTION ElementType
sum(std::experimental::mdspan<
            ElementType,
            std::experimental::extents<std::size_t, Ext>,
            LayoutPolicy,
            AccessorPolicy> const& array,
    int start,
    int end)
{
    ElementType val(0.0);
    for (int i(start); i < end; ++i) {
        val += array[i];
    }
    return val;
}

template <typename T>
KOKKOS_INLINE_FUNCTION T modulo(T x, T y)
{
    return x - y * Kokkos::floor(double(x) / y);
}

KOKKOS_INLINE_FUNCTION double ipow(double a, std::size_t i)
{
    double r(1.0);
    for (std::size_t j(0); j < i; ++j) {
        r *= a;
    }
    return r;
}

KOKKOS_INLINE_FUNCTION double ipow(double a, int i)
{
    double r(1.0);
    if (i > 0) {
        for (int j(0); j < i; ++j) {
            r *= a;
        }
    } else if (i < 0) {
        for (int j(0); j < -i; ++j) {
            r *= a;
        }
        r = 1.0 / r;
    }
    return r;
}

KOKKOS_INLINE_FUNCTION std::size_t factorial(std::size_t f)
{
    std::size_t r = 1;
    for (std::size_t i(2); i < f + 1; ++i) {
        r *= i;
    }
    return r;
}

template <class T, std::size_t D>
KOKKOS_INLINE_FUNCTION T dot_product(std::array<T, D> const& a, std::array<T, D> const& b)
{
    T result = 0;
    for (std::size_t i(0); i < D; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

template <typename T>
KOKKOS_INLINE_FUNCTION T min(T x, T y)
{
    return x < y ? x : y;
}

template <typename T>
KOKKOS_INLINE_FUNCTION T max(T x, T y)
{
    return x > y ? x : y;
}
} // namespace ddc::detail
