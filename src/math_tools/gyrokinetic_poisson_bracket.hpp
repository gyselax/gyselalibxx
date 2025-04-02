// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "metric_tensor_evaluator.hpp"
#include "vector_field.hpp"


/**
 * @brief A class which implements a gradient operator.
 * @tparam MetricTensorType A type representing a metric tensor.
 */
template <class Mapping3D>
class GyrokineticPoissonBracket
{
    static_assert(is_mapping_v<Mapping3D>);
    static_assert(Mapping3D::CoordArg::size() == 3);

    using BasisSpatial = ddc::to_type_seq_t<typename Mapping3D::CoordArg>;
    using CovBasisSpatial = get_covariant_dims_t<BasisSpatial>;

private:
    Mapping3D const& m_mapping;
    Gradient<MetricTensorEvaluator<Mapping3D>> m_grad;

public:
    GyrokineticPoissonBracket(Mapping3D const& mapping)
        : m_mapping(mapping)
        , m_grad(MetricTensorEvaluator<Mapping3D>(mapping))
    {
    }

    KOKKOS_INLINE_FUNCTION double operator()(
            DTensor<CovBasisSpatial> const& partial_derivatives_f,
            DTensor<CovBasisSpatial> const& partial_derivatives_g,
            DTensor<BasisSpatial> const& B,
            CoordArg const& coord) const
    {
        LeviCivitaTensor<double, ValidIndexSet> eps;
        double B_norm = norm(B);
        return tensor_mul(
                       index<'i', 'j', 'k'>(eps),
                       index<'i'>(m_grad(partial_derivatives_f, coord)),
                       index<'j'>(m_grad(partial_derivatives_g, coord)),
                       index<'k'>(B / B_norm))
               / m_mapping.jacobian(coord);
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
        static_assert(Kokkos::SpaceAccessibility<ExecSpace, MemorySpace>);
        using IdxType = typename IdxRange::discrete_element_type;
        ddc::parallel_for_each(
                exec_space,
                get_idx_range(poisson_bracket),
                KOKKOS_CLASS_LAMBDA(IdxType const idx) {
                    DVectorType poisson_bracket_elem = (*this)(
                            partial_derivatives_f(idx),
                            partial_derivatives_g(idx),
                            B(idx),
                            ddc::coordinate(idx));
                    ddcHelper::
                            assign_vector_field_element(poisson_bracket, idx, poisson_bracket_elem);
                });
    }
};
