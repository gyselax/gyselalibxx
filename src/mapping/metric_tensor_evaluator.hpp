// SPDX-License-Identifier: MIT
#pragma once
#include <type_traits>

#include "ddc_aliases.hpp"
#include "mapping_tools.hpp"
#include "tensor.hpp"
#include "vector_index_tools.hpp"
#include "view.hpp"

/**
 * @brief An operator for calculating the metric tensor.
 * @tparam Mapping The mapping providing the Jacobian operator.
 * @tparam PositionCoordinate The coordinate type where the metric tensor can be evaluated.
 */
template <class Mapping, class PositionCoordinate>
class MetricTensorEvaluator
{
    static_assert(is_mapping_v<Mapping>);
    static_assert(has_2d_jacobian_v<Mapping, PositionCoordinate>);

    using Dims = ddc::to_type_seq_t<typename Mapping::CoordArg>;
    using Dims_cov = vector_index_set_dual_t<Dims>;
    using Dim0 = ddc::type_seq_element_t<0, Dims>;
    using Dim1 = ddc::type_seq_element_t<1, Dims>;
    using Dim0_cov = typename Dim0::Dual;
    using Dim1_cov = typename Dim1::Dual;

public:
    /// The type of the Jacobian matrix and its inverse
    using Matrix_2x2 = std::array<std::array<double, 2>, 2>;

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
    KOKKOS_FUNCTION DTensor<Dims, Dims> operator()(PositionCoordinate const& coord) const
    {
        const double J_11 = m_mapping.jacobian_11(coord);
        const double J_12 = m_mapping.jacobian_12(coord);
        const double J_21 = m_mapping.jacobian_21(coord);
        const double J_22 = m_mapping.jacobian_22(coord);
        DTensor<Dims, Dims> metric_tensor;
        ddcHelper::get<Dim0, Dim0>(metric_tensor) = (J_11 * J_11 + J_21 * J_21);
        ddcHelper::get<Dim0, Dim1>(metric_tensor) = (J_11 * J_12 + J_21 * J_22);
        ddcHelper::get<Dim1, Dim0>(metric_tensor) = (J_11 * J_12 + J_21 * J_22);
        ddcHelper::get<Dim1, Dim1>(metric_tensor) = (J_12 * J_12 + J_22 * J_22);

        return metric_tensor;
    }

    /**
     * @brief Compute the inverse metric tensor associated with the mapping at a given position in space.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the metric tensor.
     * @return inverse_metric_tensor
     * 				A DTensor object containing the value of the inverse of the metric tensor.
     */
    KOKKOS_FUNCTION DTensor<Dims_cov, Dims_cov> inverse(PositionCoordinate const& coord) const
    {
        const double J_11 = m_mapping.jacobian_11(coord);
        const double J_12 = m_mapping.jacobian_12(coord);
        const double J_21 = m_mapping.jacobian_21(coord);
        const double J_22 = m_mapping.jacobian_22(coord);
        const double jacob_2 = m_mapping.jacobian(coord) * m_mapping.jacobian(coord);

        DTensor<Dims_cov, Dims_cov> inverse_metric_tensor;
        ddcHelper::get<Dim0_cov, Dim0_cov>(inverse_metric_tensor)
                = (J_12 * J_12 + J_22 * J_22) / jacob_2;
        ddcHelper::get<Dim0_cov, Dim1_cov>(inverse_metric_tensor)
                = (-J_11 * J_12 - J_21 * J_22) / jacob_2;
        ddcHelper::get<Dim1_cov, Dim0_cov>(inverse_metric_tensor)
                = (-J_11 * J_12 - J_21 * J_22) / jacob_2;
        ddcHelper::get<Dim1_cov, Dim1_cov>(inverse_metric_tensor)
                = (J_11 * J_11 + J_21 * J_21) / jacob_2;

        return inverse_metric_tensor;
    }

    /**
     * @brief Compute the covariant vector from the contravariant vector
     *
     * @param[in] contravariant_vector
     * 				The metric tensor matrix.
     * @param[in] coord
     * 				The coordinate where we want to compute the convariant vector.
     *
     * @return A vector of the covariant
     */
    KOKKOS_FUNCTION CovariantVectorType to_covariant(
            ContravariantVectorType const& contravariant_vector,
            PositionCoordinate const& coord) const
    {
        DTensor<Dims_cov, Dims_cov> inv_metric_tensor = inverse(coord);
        CovariantVectorType covariant_vector;
        ddcHelper::get<Dim0_cov>(covariant_vector)
                = ddcHelper::get<Dim0_cov, Dim0_cov>(inv_metric_tensor)
                          * ddcHelper::get<Dim0>(contravariant_vector)
                  + ddcHelper::get<Dim0_cov, Dim1_cov>(inv_metric_tensor)
                            * ddcHelper::get<Dim1>(contravariant_vector);
        ddcHelper::get<Dim1_cov>(covariant_vector)
                = ddcHelper::get<Dim1_cov, Dim0_cov>(inv_metric_tensor)
                          * ddcHelper::get<Dim0>(contravariant_vector)
                  + ddcHelper::get<Dim1_cov, Dim1_cov>(inv_metric_tensor)
                            * ddcHelper::get<Dim1>(contravariant_vector);
        return covariant_vector;
    }
};
