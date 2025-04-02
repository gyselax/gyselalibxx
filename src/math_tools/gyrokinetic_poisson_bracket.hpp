// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "gradient.hpp"
#include "metric_tensor_evaluator.hpp"
#include "vector_field.hpp"


/**
 * @brief A class which implements a gyrokinetic Poisson bracket operator.
 * The implemented equation is:
 * @f$ \{F, G\} = b\dot(\nabla F \cross \nabla G) @f$
 * @tparam Mapping3D A type representing a mapping in 3 dimensions.
 */
template <class Mapping3D>
class GyrokineticPoissonBracket
{
    static_assert(is_mapping_v<Mapping3D>);
    static_assert(Mapping3D::CoordArg::size() == 3);

    using BasisSpatial = ddc::to_type_seq_t<typename Mapping3D::CoordArg>;
    using CovBasisSpatial = get_covariant_dims_t<BasisSpatial>;

private:
    Mapping3D m_mapping;
    Gradient<MetricTensorEvaluator<Mapping3D>> m_grad;

public:
    /**
     * @brief Build a GyrokineticPoissonBracket operator.
     * @param mapping The mapping describing the system of coordinates on which the
     *              expression is calculated.
     */
    explicit GyrokineticPoissonBracket(Mapping3D const& mapping)
        : m_mapping(mapping)
        , m_grad(MetricTensorEvaluator<Mapping3D>(mapping))
    {
    }

    /**
     * @brief Compute the gyrokinetic Poisson bracket at a given coordinate, from
     * the partial derivatives of the two fields, and the magnetic field.
     *
     * @param[in] partial_derivatives_f A vector containing the partial derivatives
     * of the scalar field f expressed at the given coordinate.
     * @param[in] partial_derivatives_g A vector containing the partial derivatives
     * of the scalar field g expressed at the given coordinate.
     * @param[in] B A vector containing the magnetic field at the given coordinate.
     * @param[in] coord The coordinate where the calculation is carried out.
     */
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

    /**
     * @brief Compute the gyrokinetic Poisson bracket at every point on a grid from the
     * partial derivatives of the fields f and g, and the magnetic field.
     *
     * @param[in] exec_space The space (CPU/GPU) where the calculation should be executed.
     * @param[out] poisson_bracket The result of the calculation of the gyrokinetic Poisson
     * bracket at every point on a grid.
     * @param[in] partial_derivatives_f A vector field containing the partial derivatives
     * of f at each point on the grid.
     * @param[in] partial_derivatives_g A vector field containing the partial derivatives
     * of f at each point on the grid.
     * @param[in] B A vector field describing the magnetic field.
     */
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
