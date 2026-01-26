// SPDX-License-Identifier: MIT
#pragma once

#include <array>
#include <concepts>
#include <type_traits>
#include <utility>

#include "tensor.hpp"
#include "vector_index_tools.hpp"
#include "view.hpp"

namespace concepts {
/**
 * @brief A helper concept to determine if a type is a mapping.
 */
template <typename T>
concept Mapping = requires
{
    typename T::CoordArg;
    typename T::CoordResult;
}
&&requires(T const& t, typename T::CoordArg const& x)
{
    {
        t(x)
        } -> std::same_as<typename T::CoordResult>;
}

template <typename T>
concept MappingWithJacobian = Mapping<T> && requires
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
concept MappingWithInvJacobian
        = MappingWithJacobian<T> && requires(T const& t, typename T::CoordJacobian const& x)
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
concept AnalyticalMapping = Mapping<T> && requires(T const& t)
{
    {
        t.get_inverse_mapping()
        } -> Mapping;
};
} // namespace concepts

namespace mapping_detail {

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
concept is_mapping_v = concepts::Mapping<Mapping>;

template <class Mapping>
concept has_jacobian_v = concepts::MappingWithJacobian<Mapping>;

template <class Mapping>
concept has_inv_jacobian_v = concepts::MappingWithInvJacobian<Mapping>;

/// Indicates that a coordinate change operator is 2D with a curvilinear mapping showing an O-point.
template <class Mapping>
static constexpr bool is_coord_transform_with_o_point_v
        = mapping_detail::HasOPoint<std::remove_const_t<std::remove_reference_t<Mapping>>>::value;

template <class Mapping>
concept is_analytical_mapping_v = concepts::AnalyticalMapping<Mapping>;

template <class Mapping>
using inverse_mapping_t = decltype(std::declval<Mapping>().get_inverse_mapping());

template <class Mapping>
static constexpr bool has_singular_o_point_inv_jacobian_v
        = mapping_detail::SingularOPointInvJacobian<
                std::remove_const_t<std::remove_reference_t<Mapping>>>::value;
