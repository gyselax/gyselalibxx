// SPDX-License-Identifier: MIT

#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "metric_tensor_evaluator.hpp"


/**
 * @brief A class which implements a gradient operator 
 * @tparam Mapping A mapping.
 * @tparam PositionCoordinate The coordinate type where the gradient can be evaluated.
 */
template <class Mapping, class PositionCoordinate>
class Gradient
{
    using MetricTensorType = MetricTensorEvaluator<Mapping, PositionCoordinate>;

    using ContravariantVectorType = typename MetricTensorType::ContravariantVectorType;
    using CovariantVectorType = typename MetricTensorType::CovariantVectorType;

    using Dims = ddc::to_type_seq_t<typename Mapping::CoordArg>;
    using Dims_cov = vector_index_set_dual_t<Dims>;
    using Dim0_cov = ddc::type_seq_element_t<0, Dims_cov>;
    using Dim1_cov = ddc::type_seq_element_t<1, Dims_cov>;

    MetricTensorType const m_metric_tensor;

public:
    /**
     * @brief Construct an instance of the class Gradient.
     *
     * @param metric_tensor A MetricTensorEvaluator.
     */
    explicit KOKKOS_FUNCTION Gradient(MetricTensorType const metric_tensor)
        : m_metric_tensor(metric_tensor)
    {
    }

    /**
     * @brief Compute the gradient of a scalar field from the partial 
     * derivatives of the field. The gradient is expressed on the covariant 
     * basis. The components of the gradient in the covariant basis are simply 
     * equal to the value of the partial derivatives of the scalar field.
     * See [Differential operators](#docs_mathematical_and_physical_conventions__Differential_operators)
     *
     * @param[in] partial_derivatives Contains the partial derivatives of the scalar field.
     *
     * @return The components of the gradient expressed on the covariant basis.
     */
    KOKKOS_INLINE_FUNCTION CovariantVectorType
    operator()(CovariantVectorType const& partial_derivatives) const
    {
        return partial_derivatives;
    }

    /**
     * @brief Compute the gradient of a scalar field from the partial 
     * derivatives of the field. The gradient is expressed on the contravariant 
     * basis.  
     * See [Differential operators](#docs_mathematical_and_physical_conventions__Differential_operators)
     *
     * @param[in] partial_derivatives Contains the partial derivatives of the scalar field.
     * @param[in] coord The coordinate at which the gradient should be evaluated.
     *
     * @return The components of the gradient expressed on the contravariant basis.
     */
    KOKKOS_INLINE_FUNCTION ContravariantVectorType operator()(
            CovariantVectorType const& partial_derivatives,
            PositionCoordinate const& coord) const
    {
        return tensor_mul(
                index<'i', 'j'>(m_metric_tensor.inverse(coord)),
                index<'j'>(partial_derivatives));
    }
};
