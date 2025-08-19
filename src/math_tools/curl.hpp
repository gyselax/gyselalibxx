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
template <class Mapping3D, class MappingCoord = typename Mapping3D::CoordArg>
class Curl
{
    static_assert(is_mapping_v<Mapping3D>);

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
    explicit Curl(Mapping3D const& mapping)
        : m_mapping(mapping)
        , m_metric_tensor(mapping)
        , m_grad(m_metric_tensor)
    {
    }

    /**
     * @brief Compute the gyrokinetic Poisson bracket at a given coordinate, from
     * the partial derivatives of the two fields, and the magnetic field.
     *
     * @param[in] partial_derivatives_f A tensor containing the partial derivatives
     * of the vector field f expressed at the given coordinate.
     * @param[in] coord The coordinate where the calculation is carried out.
     */
    template <
            class TensorType,
            class = std::enable_if_t<is_tensor_type_v<TensorType> && TensorType::rank() == 2>>
    KOKKOS_INLINE_FUNCTION auto operator()(
            TensorType const& partial_derivatives_f,
            MappingCoord const& coord) const
    {
        double J = m_mapping.jacobian(coord);
        LeviCivitaTensor<double, BasisSpatial> eps(J);
        return tensor_mul(index<'k', 'l', 'm'>(eps), index<'l', 'm'>(partial_derivatives_f));
    }
};
