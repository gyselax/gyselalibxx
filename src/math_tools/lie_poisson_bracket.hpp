// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "gradient.hpp"
#include "metric_tensor_evaluator.hpp"
#include "static_tensors.hpp"
#include "vector_field.hpp"


/**
 * @brief A class which implements a gyrokinetic Poisson bracket operator.
 * The implemented equation is:
 * @f$ \{F, G\} = b\dot(\nabla F \cross \nabla G) @f$
 * with @f$ b= \mathbf{B} / B @f$ the unitary magnetic field, i.e:
 * @f$ \{F, G\} = {\cal J}_{\rm x}^{-1}\epsilon^{ijk}\partial_{x^i} F \partial_{x^j} G b_k @f$
 * with @f$ {\cal J}_{\rm x} @f$ the jacobian of the system,
 * @f$ b_k @f$ the covariant components of b and @f$\epsilon^{ijk} @f$ the Levi-Civita symbol.
 * @tparam Mapping3D A type representing a mapping in 3 dimensions.
 * @tparam MappingCoord The type of the coordinate that will be used to evaluate the
 *                  mapping. This coordinate is used to calculate the determinant of
 *                  the Jacobian and the metric tensor. It is almost always equal to
 *                  the argument type of the 3D mapping but axi-symmetry may make it
 *                  useful to evaluate the mapping on a coordinate with a reduced
 *                  dimensionality.
 */
template <concepts::Mapping Mapping3D, class MappingCoord = typename Mapping3D::CoordArg>
class LiePoissonBracket
{
    using BasisSpatial = ddc::to_type_seq_t<typename Mapping3D::CoordArg>;
    using CovBasisSpatial = get_covariant_dims_t<BasisSpatial>;

private:
    Mapping3D m_mapping;
    MetricTensorEvaluator<Mapping3D, MappingCoord> m_metric_tensor;
    Gradient<MetricTensorEvaluator<Mapping3D, MappingCoord>> m_grad;

public:
    /**
     * @brief Build a LiePoissonBracket operator.
     * @param mapping The mapping describing the system of coordinates on which the
     *              expression is calculated.
     */
    explicit LiePoissonBracket(Mapping3D const& mapping)
        : m_mapping(mapping)
        , m_metric_tensor(mapping)
        , m_grad(m_metric_tensor)
    {
    }

    /**
     * @brief Compute the gyrokinetic Poisson bracket at a given coordinate, from
     * the partial derivatives of the two fields, and the magnetic field.
     *
     * @param[in] partial_derivatives_f A vector containing the (covariant) gradient
     * of the scalar field f expressed at the given coordinate.
     * @param[in] partial_derivatives_g A tensor containing the partial derivatives
     * of the vector field G expressed at the given coordinate : @f$ DG_{ij} = \partial_j G_i @f$.
     * @param[in] B A contravariant vector containing the magnetic field at the given coordinate.
     * @param[in] coord The coordinate where the calculation is carried out.
     *
     * @return A vector describing the gyrokinetic Poisson bracket at a given coordinate.
     * The vector is indexed in the same way and has the same covariance as the elements of the vector G.
     */
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

    /**
     * @brief Compute the gyrokinetic Poisson bracket at a given coordinate, from
     * the partial derivatives of the two fields, and the magnetic field.
     *
     * @param[in] partial_derivatives_f A vector containing the (covariant) gradient
     * of the scalar field f expressed at the given coordinate.
     * @param[in] partial_derivatives_g A vector containing the (covariant) gradient
     * of the scalar field g expressed at the given coordinate.
     * @param[in] B A contravariant vector containing the magnetic field at the given coordinate.
     * @param[in] coord The coordinate where the calculation is carried out.
     *
     * @return The gyrokinetic Poisson bracket at a given coordinate.
     */
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
