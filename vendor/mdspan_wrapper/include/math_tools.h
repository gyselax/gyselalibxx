#pragma once

#include <algorithm>
#include <cmath>

#include <experimental/mdspan>

template <typename T>
inline T sum(T* array, int size)
{
    T val(0.0);
    for (int i(0); i < size; ++i) {
        val += array[i];
    }
    return val;
}

template <class ElementType, class LayoutPolicy, class AccessorPolicy, std::size_t Ext>
inline ElementType sum(std::experimental::mdspan<
                       ElementType,
                       std::experimental::extents<Ext>,
                       LayoutPolicy,
                       AccessorPolicy> const& array)
{
    ElementType val(0.0);
    for (int i(0); i < array.extent(0); ++i) {
        val += array[i];
    }
    return val;
}

template <class ElementType, class LayoutPolicy, class AccessorPolicy, std::size_t Ext>
inline ElementType sum(
        std::experimental::mdspan<
                ElementType,
                std::experimental::extents<Ext>,
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

inline double ipow(double a, int i)
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

using std::max;

using std::min;
