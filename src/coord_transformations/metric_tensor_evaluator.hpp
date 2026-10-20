// SPDX-License-Identifier: MIT
#pragma once
#include <concepts>
#include <type_traits>

#include "coord_transformation_tools.hpp"
#include "ddc_aliases.hpp"
#include "indexed_tensor.hpp"
#include "inverse_jacobian_matrix.hpp"
#include "tensor.hpp"
#include "vector_index_tools.hpp"
#include "view.hpp"

namespace metric_tensor_detail {
/**
 * @brief A placeholder grid used to check that a metric tensor can be evaluated at
 * an index. It is never instantiated in a discrete space.
 */
template <class Dim>
struct PlaceholderGrid
{
    using continuous_dimension_type = Dim;
};
} // namespace metric_tensor_detail

namespace concepts {
/**
 * @brief A concept describing a type which can provide the metric tensor and its
 * inverse at a given coordinate or at a given index.
 *
 * The metric tensor must be expressed on the covariant basis (i.e. it acts on
 * contravariant vectors) and its inverse must be expressed on the contravariant
 * basis (i.e. it acts on covariant vectors).
 *
 * The operators must also accept any index IdxType such that
 * std::is_same_v<ddc::coordinate_of_t<IdxType>, CoordArg>. This is checked using
 * a placeholder index type.
 */
template <typename T>
concept MetricTensor = requires
{
    typename T::CoordArg;
    typename T::ContravariantVectorType;
    typename T::CovariantVectorType;
}
&&is_tensor_type_v<typename T::ContravariantVectorType>&&
        is_tensor_type_v<typename T::CovariantVectorType>&& requires(
                T const& t,
                typename T::CoordArg const& x,
                Idx<metric_tensor_detail::PlaceholderGrid<typename T::CoordArg>> const& idx)
{
    {
        t(x)
        } -> std::same_as<
                DTensor<ddc::type_seq_element_t<0, typename T::CovariantVectorType::index_set>,
                        ddc::type_seq_element_t<0, typename T::CovariantVectorType::index_set>>>;
    {
        t.inverse(x)
        } -> std::same_as<DTensor<
                ddc::type_seq_element_t<0, typename T::ContravariantVectorType::index_set>,
                ddc::type_seq_element_t<0, typename T::ContravariantVectorType::index_set>>>;
    {
        t(idx)
        } -> std::same_as<
                DTensor<ddc::type_seq_element_t<0, typename T::CovariantVectorType::index_set>,
                        ddc::type_seq_element_t<0, typename T::CovariantVectorType::index_set>>>;
    {
        t.inverse(idx)
        } -> std::same_as<DTensor<
                ddc::type_seq_element_t<0, typename T::ContravariantVectorType::index_set>,
                ddc::type_seq_element_t<0, typename T::ContravariantVectorType::index_set>>>;
};
} // namespace concepts

/**
 * @brief An operator for calculating the metric tensor.
 * @tparam Mapping The mapping providing the Jacobian operator.
 * @tparam PositionCoordinate The coordinate type where the metric tensor can be evaluated.
 */
template <
        concepts::MappingWithJacobian Mapping,
        class PositionCoordinate = typename Mapping::CoordJacobian>
class MetricTensorEvaluator
{
    static_assert(
            std::is_same_v<PositionCoordinate, typename Mapping::CoordJacobian>,
            "The metric tensor is calculated from the Jacobian matrix. In order to evaluate the "
            "metric tensor on a coordinate different to the one given as argument to the mapping, "
            "please define a specialisation of this class.");
    static_assert(
            (is_covariant_vector_index_set_v<ddc::to_type_seq_t<typename Mapping::CoordResult>>)&&(
                    is_contravariant_vector_index_set_v<
                            ddc::to_type_seq_t<typename Mapping::CoordResult>>),
            "The metric tensor can only be calculated from the Jacobian for coordinate "
            "transformations which map to the Cartesian space.");

    using Dims = ddc::to_type_seq_t<typename Mapping::CoordArg>;
    using Dims_cov = vector_index_set_dual_t<Dims>;

public:
    /// The type of a contravariant vector associated with this mapping.
    using ContravariantVectorType = DTensor<Dims>;

    /// The type of a covariant vector associated with this mapping.
    using CovariantVectorType = DTensor<vector_index_set_dual_t<Dims>>;

    /// The type of a coordinate associated with this mapping.
    using CoordArg = PositionCoordinate;

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
    KOKKOS_FUNCTION DTensor<Dims_cov, Dims_cov> operator()(CoordArg const& coord) const
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
    KOKKOS_FUNCTION DTensor<Dims, Dims> inverse(CoordArg const& coord) const
    {
        InverseJacobianMatrix get_inverse_jacobian(m_mapping);
        Tensor inv_J = get_inverse_jacobian(coord);
        return tensor_mul(index<'i', 'j'>(inv_J), index<'k', 'j'>(inv_J));
    }

    /**
     * @brief Compute the metric tensor associated with the mapping at the coordinate of an index.
     *
     * @param[in] idx
     * 				The index where we evaluate the metric tensor.
     * @return metric_tensor
     * 				A DTensor object containing the value of the metric tensor.
     */
    template <ddc::concepts::discrete_element IdxType>
    requires(std::is_same_v<ddc::coordinate_of_t<IdxType>, CoordArg>)
            KOKKOS_FUNCTION DTensor<Dims_cov, Dims_cov>
    operator()(IdxType const& idx) const
    {
        return (*this)(ddc::coordinate(idx));
    }

    /**
     * @brief Compute the inverse metric tensor associated with the mapping at the coordinate of an index.
     *
     * @param[in] idx
     * 				The index where we evaluate the metric tensor.
     * @return inverse_metric_tensor
     * 				A DTensor object containing the value of the inverse of the metric tensor.
     */
    template <ddc::concepts::discrete_element IdxType>
    requires(std::is_same_v<ddc::coordinate_of_t<IdxType>, CoordArg>)
            KOKKOS_FUNCTION DTensor<Dims, Dims> inverse(IdxType const& idx)
    const
    {
        return inverse(ddc::coordinate(idx));
    }
};

namespace mapping_detail {
template <class Mapping, class PositionCoordinate, class ExecSpace>
struct MappingAccessibility<ExecSpace, MetricTensorEvaluator<Mapping, PositionCoordinate>>
    : MappingAccessibility<ExecSpace, Mapping>
{
};
} // namespace mapping_detail
