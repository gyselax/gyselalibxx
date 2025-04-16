

# File math\_tools.hpp

[**File List**](files.md) **>** [**math\_tools**](dir_3ced5d1c6eac490d7704c2e023d148d8.md) **>** [**math\_tools.hpp**](math__tools_8hpp.md)

[Go to the documentation of this file](math__tools_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <algorithm>
#include <cmath>

#include <ddc/ddc.hpp>

#include <Kokkos_Core.hpp>

#include "indexed_tensor.hpp"
#include "tensor.hpp"
#include "vector_field.hpp"

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

template <class ElementType, class VectorIndexSetType>
KOKKOS_INLINE_FUNCTION ElementType
norm(Tensor<ElementType, VectorIndexSetType, VectorIndexSetType> const& metric,
     Tensor<ElementType, vector_index_set_dual_t<VectorIndexSetType>> const& vec)
{
    return Kokkos::sqrt(tensor_mul(index<'i'>(vec), index<'i', 'j'>(metric), index<'j'>(vec)));
}

template <
        class ExecSpace,
        class IdxRangeType,
        class MetricTensorEvaluator,
        class VectorIndexSetType>
void norm(
        ExecSpace exec_space,
        DField<IdxRangeType, typename ExecSpace::memory_space> norm_vals,
        MetricTensorEvaluator const& get_metric,
        DVectorConstField<IdxRangeType, VectorIndexSetType, typename ExecSpace::memory_space> vals)
{
    using IdxType = typename IdxRangeType::discrete_element_type;

    if constexpr (is_contravariant_vector_index_set_v<VectorIndexSetType>) {
        ddc::parallel_for_each(
                exec_space,
                get_idx_range(vals),
                KOKKOS_LAMBDA(IdxType idx) {
                    Tensor metric = get_metric(ddc::coordinate(idx));
                    norm_vals(idx) = norm(metric, vals(idx));
                });
    } else {
        ddc::parallel_for_each(
                exec_space,
                get_idx_range(vals),
                KOKKOS_LAMBDA(IdxType idx) {
                    Tensor metric = get_metric.inverse(ddc::coordinate(idx));
                    norm_vals(idx) = norm(metric, vals(idx));
                });
    }
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

template <class RowDim1, class RowDim2, class ColDim1, class ColDim2>
KOKKOS_INLINE_FUNCTION double determinant(
        DTensor<VectorIndexSet<RowDim1, RowDim2>, VectorIndexSet<ColDim1, ColDim2>> arr)
{
    return ddcHelper::get<RowDim1, ColDim1>(arr) * ddcHelper::get<RowDim2, ColDim2>(arr)
           - ddcHelper::get<RowDim2, ColDim1>(arr) * ddcHelper::get<RowDim1, ColDim2>(arr);
}

template <class RowDim1, class RowDim2, class RowDim3, class ColDim1, class ColDim2, class ColDim3>
KOKKOS_INLINE_FUNCTION double determinant(
        DTensor<VectorIndexSet<RowDim1, RowDim2, RowDim3>,
                VectorIndexSet<ColDim1, ColDim2, ColDim3>> arr)
{
    return ddcHelper::get<RowDim1, ColDim1>(arr) * ddcHelper::get<RowDim2, ColDim2>(arr)
                   * ddcHelper::get<RowDim3, ColDim3>(arr)
           + ddcHelper::get<RowDim1, ColDim2>(arr) * ddcHelper::get<RowDim2, ColDim3>(arr)
                     * ddcHelper::get<RowDim3, ColDim1>(arr)
           + ddcHelper::get<RowDim1, ColDim3>(arr) * ddcHelper::get<RowDim2, ColDim1>(arr)
                     * ddcHelper::get<RowDim3, ColDim2>(arr)
           - ddcHelper::get<RowDim3, ColDim1>(arr) * ddcHelper::get<RowDim2, ColDim2>(arr)
                     * ddcHelper::get<RowDim1, ColDim3>(arr)
           - ddcHelper::get<RowDim3, ColDim2>(arr) * ddcHelper::get<RowDim2, ColDim3>(arr)
                     * ddcHelper::get<RowDim1, ColDim1>(arr)
           - ddcHelper::get<RowDim3, ColDim3>(arr) * ddcHelper::get<RowDim2, ColDim1>(arr)
                     * ddcHelper::get<RowDim1, ColDim2>(arr);
}

template <class RowDim1, class RowDim2, class ColDim1, class ColDim2>
KOKKOS_INLINE_FUNCTION DTensor<
        VectorIndexSet<typename ColDim1::Dual, typename ColDim2::Dual>,
        VectorIndexSet<typename RowDim1::Dual, typename RowDim2::Dual>>
inverse(DTensor<VectorIndexSet<RowDim1, RowDim2>, VectorIndexSet<ColDim1, ColDim2>> arr)
{
    using OutRowDim1 = typename ColDim1::Dual;
    using OutRowDim2 = typename ColDim2::Dual;
    using OutColDim1 = typename RowDim1::Dual;
    using OutColDim2 = typename RowDim2::Dual;
    DTensor<VectorIndexSet<OutRowDim1, OutRowDim2>, VectorIndexSet<OutColDim1, OutColDim2>> inv;
    double det = determinant(arr);
    ddcHelper::get<OutRowDim1, OutColDim1>(inv) = ddcHelper::get<RowDim2, ColDim2>(arr) / det;
    ddcHelper::get<OutRowDim2, OutColDim1>(inv) = -ddcHelper::get<RowDim2, ColDim1>(arr) / det;
    ddcHelper::get<OutRowDim1, OutColDim2>(inv) = -ddcHelper::get<RowDim1, ColDim2>(arr) / det;
    ddcHelper::get<OutRowDim2, OutColDim2>(inv) = ddcHelper::get<RowDim1, ColDim1>(arr) / det;
    return inv;
}
```


