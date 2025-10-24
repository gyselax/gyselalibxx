

# File lie\_poisson\_bracket.hpp

[**File List**](files.md) **>** [**math\_tools**](dir_3ced5d1c6eac490d7704c2e023d148d8.md) **>** [**lie\_poisson\_bracket.hpp**](lie__poisson__bracket_8hpp.md)

[Go to the documentation of this file](lie__poisson__bracket_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "gradient.hpp"
#include "metric_tensor_evaluator.hpp"
#include "static_tensors.hpp"
#include "vector_field.hpp"


template <class Mapping3D, class MappingCoord = typename Mapping3D::CoordArg>
class LiePoissonBracket
{
    static_assert(is_mapping_v<Mapping3D>);

    using BasisSpatial = ddc::to_type_seq_t<typename Mapping3D::CoordArg>;
    using CovBasisSpatial = get_covariant_dims_t<BasisSpatial>;

private:
    Mapping3D m_mapping;
    MetricTensorEvaluator<Mapping3D, MappingCoord> m_metric_tensor;
    Gradient<MetricTensorEvaluator<Mapping3D, MappingCoord>> m_grad;

public:
    explicit LiePoissonBracket(Mapping3D const& mapping)
        : m_mapping(mapping)
        , m_metric_tensor(mapping)
        , m_grad(m_metric_tensor)
    {
    }

    template <
            class TensorType,
            class = std::enable_if_t<is_tensor_type_v<TensorType> && TensorType::rank() == 2>>
    KOKKOS_INLINE_FUNCTION DTensor<ddc::type_seq_element_t<0, typename TensorType::index_set>>
    operator()(
            DTensor<CovBasisSpatial> const& partial_derivatives_f,
            TensorType const& partial_derivatives_g,
            DTensor<BasisSpatial> const& B,
            MappingCoord const& coord) const
    {
        double J = m_mapping.jacobian(coord);
        LeviCivitaTensor<double, BasisSpatial> eps(J);
        DTensor<CovBasisSpatial, CovBasisSpatial> metric_tensor = m_metric_tensor(coord);
        double B_norm = norm(metric_tensor, B);
        return tensor_mul(
                index<'i', 'j', 'k'>(eps),
                index<'i', 'l'>(metric_tensor),
                index<'l'>(B / B_norm),
                index<'j'>(partial_derivatives_f),
                index<'m', 'k'>(partial_derivatives_g));
    }

    KOKKOS_INLINE_FUNCTION double operator()(
            DTensor<CovBasisSpatial> const& partial_derivatives_f,
            DTensor<CovBasisSpatial> const& partial_derivatives_g,
            DTensor<BasisSpatial> const& B,
            MappingCoord const& coord) const
    {
        double J = m_mapping.jacobian(coord);
        LeviCivitaTensor<double, BasisSpatial> eps(J);
        DTensor<CovBasisSpatial, CovBasisSpatial> metric_tensor = m_metric_tensor(coord);
        double B_norm = norm(metric_tensor, B);
        return tensor_mul(
                index<'i', 'j', 'k'>(eps),
                index<'i', 'l'>(metric_tensor),
                index<'l'>(B / B_norm),
                index<'j'>(partial_derivatives_f),
                index<'k'>(partial_derivatives_g));
    }

    template <class ExecSpace, class IdxRange, class MemorySpace>
    void operator()(
            ExecSpace exec_space,
            DField<IdxRange, MemorySpace> poisson_bracket,
            DVectorConstField<IdxRange, CovBasisSpatial, MemorySpace> const partial_derivatives_f,
            DVectorConstField<IdxRange, CovBasisSpatial, MemorySpace> const partial_derivatives_g,
            DVectorConstField<IdxRange, BasisSpatial, MemorySpace> const B)
    {
        static_assert(is_accessible_v<ExecSpace, Mapping3D>);
        static_assert(Kokkos::SpaceAccessibility<ExecSpace, MemorySpace>::accessible);
        using IdxType = typename IdxRange::discrete_element_type;
        ddc::parallel_for_each(
                exec_space,
                get_idx_range(poisson_bracket),
                KOKKOS_CLASS_LAMBDA(IdxType const idx) {
                    poisson_bracket(idx) = (*this)(
                            partial_derivatives_f(idx),
                            partial_derivatives_g(idx),
                            B(idx),
                            ddc::coordinate(idx));
                });
    }
};
```


