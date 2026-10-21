

# File metric\_tensor\_evaluator.hpp

[**File List**](files.md) **>** [**coord\_transformations**](dir_67161c4ffadea73fddf46ea451c2f62c.md) **>** [**metric\_tensor\_evaluator.hpp**](metric__tensor__evaluator_8hpp.md)

[Go to the documentation of this file](metric__tensor__evaluator_8hpp.md)


```C++
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
template <class Dim>
struct PlaceholderGrid
{
    using continuous_dimension_type = Dim;
};

template <class CoordType>
struct IdxPlaceholder;

template <class... Dims>
struct IdxPlaceholder<Coord<Dims...>>
{
    using type = Idx<PlaceholderGrid<Dims>...>;
};

template <class CoordType>
using placeholder_idx_t = typename IdxPlaceholder<CoordType>::type;
} // namespace metric_tensor_detail

namespace concepts {
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
                metric_tensor_detail::placeholder_idx_t<typename T::CoordArg> const& idx)
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
    using ContravariantVectorType = DTensor<Dims>;

    using CovariantVectorType = DTensor<vector_index_set_dual_t<Dims>>;

    using CoordArg = PositionCoordinate;

private:
    Mapping m_mapping;

public:
    explicit KOKKOS_FUNCTION MetricTensorEvaluator(Mapping mapping) : m_mapping(mapping) {}

    KOKKOS_FUNCTION DTensor<Dims_cov, Dims_cov> operator()(CoordArg const& coord) const
    {
        Tensor J = m_mapping.jacobian_matrix(coord);
        return tensor_mul(index<'j', 'i'>(J), index<'j', 'k'>(J));
    }

    KOKKOS_FUNCTION DTensor<Dims, Dims> inverse(CoordArg const& coord) const
    {
        InverseJacobianMatrix get_inverse_jacobian(m_mapping);
        Tensor inv_J = get_inverse_jacobian(coord);
        return tensor_mul(index<'i', 'j'>(inv_J), index<'k', 'j'>(inv_J));
    }

    template <ddc::concepts::discrete_element IdxType>
    requires(std::is_same_v<ddc::coordinate_of_t<IdxType>, CoordArg>)
            KOKKOS_FUNCTION DTensor<Dims_cov, Dims_cov>
    operator()(IdxType const& idx) const
    {
        return (*this)(ddc::coordinate(idx));
    }

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
```


