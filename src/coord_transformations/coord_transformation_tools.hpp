// SPDX-License-Identifier: MIT
#pragma once

#include <array>
#include <concepts>
#include <type_traits>
#include <utility>

#include "tensor.hpp"
#include "vector_index_tools.hpp"
#include "view.hpp"

namespace mapping_detail {

/**
 * @brief A helper concept to determine if a type is a mapping.
 */
template <typename T>
concept IsMapping = requires
{
    typename T::CoordArg;
    typename T::CoordResult;
}
&&std::invocable<T, typename T::CoordArg>&&
        std::same_as<std::invoke_result_t<T, typename T::CoordArg>, typename T::CoordResult>;

template <typename T>
concept DefinesJacobian = IsMapping<T> && requires
{
    typename T::CoordJacobian;
} && requires(T const& t, typename T::CoordJacobian const& x)
{
    {
        t.jacobian_matrix(x)
        } -> std::same_as<
                DTensor<get_contravariant_dims_t<ddc::to_type_seq_t<typename T::CoordResult>>,
                        get_covariant_dims_t<ddc::to_type_seq_t<typename T::CoordArg>>>>;
    {
        t.template jacobian_component<int, int>(x)
        } -> std::same_as<double>;
    {
        t.jacobian(x)
        } -> std::same_as<double>;
};

template <typename T>
concept DefinesInvJacobian
        = DefinesJacobian<T> && requires(T const& t, typename T::CoordJacobian const& x)
{
    {
        t.inv_jacobian_matrix(x)
        } -> std::same_as<
                DTensor<get_contravariant_dims_t<ddc::to_type_seq_t<typename T::CoordArg>>,
                        get_covariant_dims_t<ddc::to_type_seq_t<typename T::CoordResult>>>>;
    {
        t.template inv_jacobian_component<int, int>(x)
        } -> std::same_as<double>;
};

template <typename T>
concept IsAnalyticalMapping = IsMapping<T> && requires(T const& t)
{
    {
    t.get_inverse_mapping()
    } -> IsMapping;
};


template <class ExecSpace, class Type>
struct MappingAccessibility : std::false_type
{
};


template <class Mapping>
struct HasOPoint : std::false_type
{
};

template <class Mapping>
struct SingularOPointInvJacobian : std::false_type
{
};

} // namespace mapping_detail

template <class ExecSpace, class Type>
static constexpr bool is_accessible_v = mapping_detail::
        MappingAccessibility<ExecSpace, std::remove_const_t<std::remove_reference_t<Type>>>::value;

template <class Mapping>
static constexpr bool is_mapping_v = mapping_detail::IsMapping<Mapping>;

template <class Mapping, bool RaiseError = true>
static constexpr bool has_jacobian_v = mapping_detail::DefinesJacobian<Mapping>;

template <class Mapping, bool RaiseError = true>
static constexpr bool has_inv_jacobian_v = mapping_detail::DefinesInvJacobian<Mapping>;

/// Indicates that a coordinate change operator is 2D with a curvilinear mapping showing an O-point.
template <class Mapping>
static constexpr bool is_coord_transform_with_o_point_v
        = mapping_detail::HasOPoint<std::remove_const_t<std::remove_reference_t<Mapping>>>::value;

template <class Mapping>
static constexpr bool is_analytical_mapping_v = mapping_detail::IsAnalyticalMapping<Mapping>;

template <class Mapping>
using inverse_mapping_t = decltype(std::declval<Mapping>().get_inverse_mapping());

template <class Mapping>
static constexpr bool has_singular_o_point_inv_jacobian_v
        = mapping_detail::SingularOPointInvJacobian<
                std::remove_const_t<std::remove_reference_t<Mapping>>>::value;
