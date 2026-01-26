// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "gradient.hpp"
#include "metric_tensor_evaluator.hpp"
#include "static_tensors.hpp"
#include "vector_field.hpp"

/**
 * @brief A class which implements a curl operator
 * The implemented equation is:
 * @f$ \nabla \cross \mathbf{F} @f$
 * @f$ \nabla \cross \mathbf{F} = {\cal J}_{\rm x}^{-1}\epsilon^{klm}\partial_{x^l} F_m @f$
 * with @f$ {\cal J}_{\rm x} @f$ the jacobian of the system,
 * @f$ F_m @f$ the covariant components of F and @f$\epsilon^{klm} @f$ the Levi-Civita symbol.
 * @tparam Mapping3D A type representing a mapping in 3 dimensions.
 * @tparam MappingCoord The type of the coordinate that will be used to evaluate the
 *                  mapping. This coordinate is used to calculate the determinant of
 *                  the Jacobian and the metric tensor. It is almost always equal to
 *                  the argument type of the 3D mapping but axi-symmetry may make it
 *                  useful to evaluate the mapping on a coordinate with a reduced
 *                  dimensionality.
 */
template <concepts::Mapping Mapping3D, class MappingCoord = typename Mapping3D::CoordArg>
class Curl
{
    using BasisSpatial = ddc::to_type_seq_t<typename Mapping3D::CoordArg>;
    using CovBasisSpatial = get_covariant_dims_t<BasisSpatial>;

private:
    Mapping3D m_mapping;
    MetricTensorEvaluator<Mapping3D, MappingCoord> m_metric_tensor;

public:
    /**
     * @brief Build a curl operator.
     * @param mapping The mapping describing the system of coordinates on which the
     *              expression is calculated.
     */
    explicit Curl(Mapping3D const& mapping) : m_mapping(mapping), m_metric_tensor(mapping) {}

    /**
     * @brief Compute a curl at a given coordinate, from
     * the partial derivatives of the vector field f
     *
     * @param[in] partial_derivatives_f A tensor containing the partial derivatives
     * of the vector field f expressed at the given coordinate.
     * @param[in] coord The coordinate where the calculation is carried out.
     */
    template <class T, class CovBasisF>
    KOKKOS_INLINE_FUNCTION Tensor<T, get_contravariant_dims_t<CovBasisF>> operator()(
            Tensor<T, CovBasisF, CovBasisSpatial> const& partial_derivatives_f,
            MappingCoord const& coord) const
    {
        static_assert(
                is_covariant_vector_index_set_v<CovBasisF>,
                "Expected the derivative of f expressed in the covariant basis");
        double J = m_mapping.jacobian(coord);
        LeviCivitaTensor<double, BasisSpatial> eps(J);

        // The following code computes: 1/sqrt(g) * eps^{klm} \nabla_{l} F_{m}
        return tensor_mul(index<'k', 'l', 'm'>(eps), index<'m', 'l'>(partial_derivatives_f));
    }
};
