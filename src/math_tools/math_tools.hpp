// SPDX-License-Identifier: MIT
#pragma once

#include <algorithm>
#include <cmath>

#include <ddc/ddc.hpp>

#include <Kokkos_Core.hpp>

#include "tensor.hpp"

template <typename T>
KOKKOS_INLINE_FUNCTION T sum(const T* array, int size)
{
    T val(0.0);
    for (int i(0); i < size; ++i) {
        val += array[i];
    }
    return val;
}

template <class ElementType, class LayoutPolicy, class AccessorPolicy, std::size_t Ext>
KOKKOS_INLINE_FUNCTION ElementType sum(Kokkos::mdspan<
                                       ElementType,
                                       Kokkos::extents<std::size_t, Ext>,
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
sum(Kokkos::mdspan<
            ElementType,
            Kokkos::extents<std::size_t, Ext>,
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
inline T modulo(T x, T y)
{
    return x - y * std::floor(double(x) / y);
}

KOKKOS_INLINE_FUNCTION constexpr double ipow(double a, std::size_t i)
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

inline std::size_t factorial(std::size_t f)
{
    std::size_t r = 1;
    for (std::size_t i(2); i < f + 1; ++i) {
        r *= i;
    }
    return r;
}

template <class T, class... Dims>
KOKKOS_INLINE_FUNCTION T
dot_product(Vector<T, Dims...> const& a, Vector<T, typename Dims::Dual...> const& b)
{
    return ((ddcHelper::get<Dims>(a) * ddcHelper::get<typename Dims::Dual>(b)) + ...);
}


template <typename T>
inline T min(T x, T y)
{
    return x < y ? x : y;
}

template <typename T>
inline T max(T x, T y)
{
    return x > y ? x : y;
}

template <std::size_t N, std::size_t M, std::size_t P>
KOKKOS_INLINE_FUNCTION std::array<std::array<double, N>, P> mat_mul(
        std::array<std::array<double, N>, M> const& a,
        std::array<std::array<double, M>, P> const& b)
{
    std::array<std::array<double, N>, P> result;
    for (std::size_t i(0); i < N; ++i) {
        for (std::size_t j(0); j < P; ++j) {
            result[i][j] = 0.0;
            for (std::size_t k(0); k < M; ++k) {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return result;
}

template <std::size_t N, std::size_t M>
KOKKOS_INLINE_FUNCTION std::array<double, N> mat_vec_mul(
        std::array<std::array<double, N>, M> const& a,
        std::array<double, M> const& b)
{
    std::array<double, N> result;
    for (std::size_t i(0); i < N; ++i) {
        result[i] = 0.0;
        for (std::size_t k(0); k < M; ++k) {
            result[i] += a[i][k] * b[k];
        }
    }
    return result;
}

KOKKOS_INLINE_FUNCTION double determinant(std::array<std::array<double, 2>, 2> arr)
{
    return arr[0][0] * arr[1][1] - arr[0][1] * arr[1][0];
}

KOKKOS_INLINE_FUNCTION std::array<std::array<double, 2>, 2> inverse(
        std::array<std::array<double, 2>, 2> arr)
{
    std::array<std::array<double, 2>, 2> inv;
    double det = determinant(arr);
    inv[0][0] = arr[1][1] / det;
    inv[1][0] = -arr[1][0] / det;
    inv[0][1] = -arr[0][1] / det;
    inv[1][1] = arr[0][0] / det;
    return inv;
}
