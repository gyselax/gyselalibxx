

# File gradient.hpp

[**File List**](files.md) **>** [**math\_tools**](dir_3ced5d1c6eac490d7704c2e023d148d8.md) **>** [**gradient.hpp**](gradient_8hpp.md)

[Go to the documentation of this file](gradient_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "metric_tensor_evaluator.hpp"
#include "vector_field.hpp"


template <class MetricTensorType>
class Gradient
{
    using CoordArg = typename MetricTensorType::CoordArg;

    using Dims = ddc::to_type_seq_t<typename MetricTensorType::CoordArg>;
    using Dims_cov = vector_index_set_dual_t<Dims>;

    using DVectorType = typename MetricTensorType::ContravariantVectorType;
    using DVectorCov = typename MetricTensorType::CovariantVectorType;

    template <class IdxRange>
    using DVectorFieldType = DVectorField<IdxRange, Dims>;

    template <class IdxRange>
    using DVectorFieldCovType = DVectorField<IdxRange, Dims_cov>;

    template <class IdxRange>
    using DVectorConstFieldCovType = DVectorConstField<IdxRange, Dims_cov>;

    MetricTensorType const& m_metric_tensor;

public:
    explicit KOKKOS_FUNCTION Gradient(MetricTensorType const& metric_tensor)
        : m_metric_tensor(metric_tensor)
    {
    }

    KOKKOS_INLINE_FUNCTION DVectorCov operator()(DVectorCov const& partial_derivatives) const
    {
        return partial_derivatives;
    }

    template <class IdxRange>
    KOKKOS_INLINE_FUNCTION void operator()(
            DVectorFieldCovType<IdxRange> gradient,
            DVectorConstFieldCovType<IdxRange> const partial_derivatives) const
    {
        ddcHelper::deepcopy(gradient, partial_derivatives);
    }

    KOKKOS_INLINE_FUNCTION DVectorType
    operator()(DVectorCov const& partial_derivatives, CoordArg const& coord) const
    {
        return tensor_mul(
                index<'i', 'j'>(m_metric_tensor.inverse(coord)),
                index<'j'>(partial_derivatives));
    }

    template <class IdxRange>
    KOKKOS_FUNCTION void operator()(
            DVectorFieldType<IdxRange> const gradient,
            DVectorConstFieldCovType<IdxRange> const partial_derivatives) const
    {
        using IdxType = typename IdxRange::discrete_element_type;
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(gradient),
                KOKKOS_CLASS_LAMBDA(IdxType const idx) {
                    DVectorType gradient_proxy
                            = (*this)(partial_derivatives(idx), ddc::coordinate(idx));
                    ddcHelper::assign_vector_field_element(gradient, idx, gradient_proxy);
                });
    }
};
```


