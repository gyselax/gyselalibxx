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
    explicit KOKKOS_FUNCTION Gradient(
            MetricTensorType const metric_tensor)
        : m_metric_tensor(metric_tensor)
    {
    }

    /**
     * @brief Compute the partial derivative of a field in the direction 
     * where the field is represented using 1d splines.
     *
     * @param[out] differentiated_field Contains on output the value of the differentiated field.
     */
    KOKKOS_FUNCTION CovariantVectorType operator()(CovariantVectorType const& partial_derivatives) const
    {
        CovariantVectorType gradient;
        ddcHelper::get<Dim0_cov>(gradient) = ddcHelper::get<Dim0_cov>(partial_derivatives);
        ddcHelper::get<Dim1_cov>(gradient) = ddcHelper::get<Dim1_cov>(partial_derivatives);

        return gradient;
    }

    /**
     * @brief Compute the partial derivative of a field in the direction 
     * where the field is represented using 1d splines.
     *
     * @param[out] differentiated_field Contains on output the value of the differentiated field.
     */
    KOKKOS_FUNCTION ContravariantVectorType operator()(ContravariantVectorType const& partial_derivatives, PositionCoordinate const& coord) const
    {
        ContravariantVectorType gradient;

        return tensor_mul(
                index<'i', 'j'>(m_metric_tensor.inverse(coord)),
                index<'j'>(partial_derivatives));
    }
};
