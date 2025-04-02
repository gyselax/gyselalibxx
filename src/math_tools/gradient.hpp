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
    /**
     * @brief Construct an instance of the class Gradient.
     *
     * @param metric_tensor A MetricTensorEvaluator.
     */
    explicit KOKKOS_FUNCTION Gradient(MetricTensorType const& metric_tensor)
        : m_metric_tensor(metric_tensor)
    {
    }

    /**
     * @brief Compute the gradient of a scalar field at a given coordinate, 
     * from the partial derivatives of the field. The gradient is expressed 
     * on the covariant basis. The components of the gradient in the covariant 
     * basis are simply equal to the value of the partial derivatives of the 
     * scalar field.
     *
     * @param[in] partial_derivatives A vector containing the partial derivatives 
     * of the scalar field expressed at a given coordinate.
     * 
     * @return The components of the gradient at a given coordinate,
     * expressed on the covariant basis.
     */
    KOKKOS_INLINE_FUNCTION DVectorCov operator()(DVectorCov const& partial_derivatives) const
    {
        return partial_derivatives;
    }

    /**
     * @brief Compute the gradient of a scalar field at a given coordinate, 
     * using partial derivatives of the field. The gradient is expressed 
     * on the covariant basis. The components of the gradient in the covariant 
     * basis are simply equal to the value of the partial derivatives of the 
     * scalar field.
     *
     * @param[out] gradient A vector field that contains on output the 
     * value of the gradient components expressed on the contravariant
     * @param[in] partial_derivatives A vector field containing the 
     * partial derivatives of the scalar field.
     * basis. 
     */
    template <class IdxRange>
    KOKKOS_INLINE_FUNCTION void operator()(
            DVectorFieldCovType<IdxRange> gradient,
            DVectorConstFieldCovType<IdxRange> const partial_derivatives) const
    {
        ddcHelper::deepcopy(gradient, partial_derivatives);
    }

    /**
     * @brief Compute the gradient of a scalar field at a given coordinate, 
     * using partial derivatives of the field. The gradient is expressed 
     * on the contravariant basis.  
     *
     * @param[in] partial_derivatives A vector that contains the partial 
     * derivatives of the scalar field expressed at a given coordinate.
     * @param[in] coord The coordinate at which the gradient should be evaluated.
     *
     * @return The components of the gradient at a given coordinate, 
     * expressed on the contravariant basis.
     */
    KOKKOS_INLINE_FUNCTION DVectorType
    operator()(DVectorCov const& partial_derivatives, CoordArg const& coord) const
    {
        return tensor_mul(
                index<'i', 'j'>(m_metric_tensor.inverse(coord)),
                index<'j'>(partial_derivatives));
    }

    /**
     * @brief Compute the gradient of a scalar field using partial derivatives
     * of the field. The gradient is expressed on the contravariant basis.  
     *
     * @param[out] gradient A vector field that contains on output the 
     * value of the gradient components expressed on the contravariant
     * @param[in] partial_derivatives A vector field that contains the partial 
     * derivatives of the scalar field.
     * basis. 
     */
    template <class IdxRange>
    KOKKOS_FUNCTION void operator()(
            DVectorFieldType<IdxRange> const gradient,
            DVectorConstFieldCovType<IdxRange> const partial_derivatives) const
    {
        static_assert(is_accessible_v<Kokkos::DefaultExecutionSpace, MetricTensorType>);
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
