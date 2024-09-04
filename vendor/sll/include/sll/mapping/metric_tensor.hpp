#pragma once
#include <type_traits>

#include <sll/view.hpp>

/**
 * @brief An operator for calculating the metric tensor.
 * @tparam Mapping The mapping providing the Jacobian operator.
 * @tparam PositionCoordinate The coordinate type where the metric tensor can be evaluated.
 */
template <class Mapping, class PositionCoordinate>
class MetricTensor
{
public:
    /// The type of the Jacobian matrix and its inverse
    using Matrix_2x2 = std::array<std::array<double, 2>, 2>;

private:
    Mapping m_mapping;

public:
    /**
     * @brief A constructor for the metric tensor operator.
     *
     * @param[in] mapping The mapping which can be used to calculate the Jacobian.
     */
    MetricTensor(Mapping mapping) : m_mapping(mapping) {}

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
    void operator()(Matrix_2x2& matrix, PositionCoordinate const& coord) const
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
    void inverse(Matrix_2x2& matrix, PositionCoordinate const& coord) const
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
    std::array<double, 2> to_covariant(
            std::array<double, 2> const& contravariant_vector,
            PositionCoordinate const& coord) const
    {
        Matrix_2x2 inv_metric_tensor;
        inverse(inv_metric_tensor, coord);
        std::array<double, 2> covariant_vector;
        covariant_vector[0] = inv_metric_tensor[0][0] * contravariant_vector[0]
                              + inv_metric_tensor[0][1] * contravariant_vector[1];
        covariant_vector[1] = inv_metric_tensor[1][0] * contravariant_vector[0]
                              + inv_metric_tensor[1][1] * contravariant_vector[1];
        return covariant_vector;
    }
};
