// SPDX-License-Identifier: MIT
#pragma once
#include <type_traits>

#include "ddc_aliases.hpp"
#include "indexed_tensor.hpp"
#include "inverse_jacobian_matrix.hpp"
#include "mapping_tools.hpp"
#include "tensor.hpp"
#include "vector_index_tools.hpp"
#include "view.hpp"

/**
 * @brief An operator for calculating the metric tensor.
 * @tparam Mapping The mapping providing the Jacobian operator.
 * @tparam PositionCoordinate The coordinate type where the metric tensor can be evaluated.
 */
template <class Mapping, class PositionCoordinate = typename Mapping::CoordArg>
class MetricTensorEvaluator
{
    static_assert(is_mapping_v<Mapping>);
    static_assert(has_2d_jacobian_v<Mapping, PositionCoordinate>);

    using Dims = ddc::to_type_seq_t<typename Mapping::CoordArg>;
    using Dims_cov = vector_index_set_dual_t<Dims>;

public:
    /// The type of a contravariant vector associated with this mapping.
    using ContravariantVectorType = DTensor<Dims>;

    /// The type of a covariant vector associated with this mapping.
    using CovariantVectorType = DTensor<vector_index_set_dual_t<Dims>>;

private:
    Mapping m_mapping;

public:
    /**
     * @brief A constructor for the metric tensor operator.
     *
     * @param[in] mapping The mapping which can be used to calculate the Jacobian.
     */
    explicit KOKKOS_FUNCTION MetricTensorEvaluator(Mapping mapping) : m_mapping(mapping) {}

    /**
     * @brief Compute the metric tensor associated with the mapping at a given position in space.
     *
     * The metric tensor matrix is defined as:
     * @f$ G = (J_{\mathcal{F}})^T J_{\mathcal{F}} @f$.
     * with @f$ J_{\mathcal{F}} @f$ the Jacobian matrix.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the metric tensor.
     * @return metric_tensor
     * 				A DTensor object containing the value of the metric tensor.
     */
    KOKKOS_FUNCTION DTensor<Dims_cov, Dims_cov> operator()(PositionCoordinate const& coord) const
    {
        Tensor J = m_mapping.jacobian_matrix(coord);
        return tensor_mul(index<'j', 'i'>(J), index<'j', 'k'>(J));
    }

    /**
     * @brief Compute the inverse metric tensor associated with the mapping at a given position in space.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the metric tensor.
     * @return inverse_metric_tensor
     * 				A DTensor object containing the value of the inverse of the metric tensor.
     */
    KOKKOS_FUNCTION DTensor<Dims, Dims> inverse(PositionCoordinate const& coord) const
    {
        InverseJacobianMatrix<Mapping, PositionCoordinate> get_inverse_jacobian(m_mapping);
        Tensor inv_J = get_inverse_jacobian(coord);
        return tensor_mul(index<'i', 'j'>(inv_J), index<'k', 'j'>(inv_J));
    }
};
