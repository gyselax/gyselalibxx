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
class MetricTensor
{
    static_assert(is_mapping_v<Mapping>);
    static_assert(has_2d_jacobian_v<Mapping, PositionCoordinate>);

    using Dims = ddc::to_type_seq_t<typename Mapping::CoordArg>;
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
    explicit KOKKOS_FUNCTION MetricTensor(Mapping mapping) : m_mapping(mapping) {}

    /**
     * @brief Compute the metric tensor assignd to the mapping.
     *
     * The metric tensor matrix is defined as:
     * @f$ G = (J_{\mathcal{F}})^T J_{\mathcal{F}} @f$.
     * with @f$ J_{\mathcal{F}} @f$ the Jacobian matrix.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the metric tensor.
     * @param[out] matrix
     * 				The metric tensor matrix.
     */
    KOKKOS_FUNCTION void operator()(Matrix_2x2& matrix, PositionCoordinate const& coord) const
    {
        const double J_11 = m_mapping.jacobian_11(coord);
        const double J_12 = m_mapping.jacobian_12(coord);
        const double J_21 = m_mapping.jacobian_21(coord);
        const double J_22 = m_mapping.jacobian_22(coord);
        matrix[0][0] = (J_11 * J_11 + J_21 * J_21);
        matrix[0][1] = (J_11 * J_12 + J_21 * J_22);
        matrix[1][0] = (J_11 * J_12 + J_21 * J_22);
        matrix[1][1] = (J_12 * J_12 + J_22 * J_22);
    }

    /**
     * @brief Compute the inverse metric tensor associated to the mapping.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the metric tensor.
     * @param[out] matrix
     * 				The metric tensor matrix.
     */
    KOKKOS_FUNCTION void inverse(Matrix_2x2& matrix, PositionCoordinate const& coord) const
    {
        const double J_11 = m_mapping.jacobian_11(coord);
        const double J_12 = m_mapping.jacobian_12(coord);
        const double J_21 = m_mapping.jacobian_21(coord);
        const double J_22 = m_mapping.jacobian_22(coord);
        const double jacob_2 = m_mapping.jacobian(coord) * m_mapping.jacobian(coord);
        matrix[0][0] = (J_12 * J_12 + J_22 * J_22) / jacob_2;
        matrix[0][1] = (-J_11 * J_12 - J_21 * J_22) / jacob_2;
        matrix[1][0] = (-J_11 * J_12 - J_21 * J_22) / jacob_2;
        matrix[1][1] = (J_11 * J_11 + J_21 * J_21) / jacob_2;
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
        Matrix_2x2 inv_metric_tensor;
        inverse(inv_metric_tensor, coord);
        CovariantVectorType covariant_vector;
        ddcHelper::get<Dim0_cov>(covariant_vector)
                = inv_metric_tensor[0][0] * ddcHelper::get<Dim0>(contravariant_vector)
                  + inv_metric_tensor[0][1] * ddcHelper::get<Dim1>(contravariant_vector);
        ddcHelper::get<Dim1_cov>(covariant_vector)
                = inv_metric_tensor[1][0] * ddcHelper::get<Dim0>(contravariant_vector)
                  + inv_metric_tensor[1][1] * ddcHelper::get<Dim1>(contravariant_vector);
        return covariant_vector;
    }
};
