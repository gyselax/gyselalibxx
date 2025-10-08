// SPDX-License-Identifier: MIT
#pragma once

#include <array>
#include <type_traits>
#include <utility>

#include "tensor.hpp"
#include "vector_index_tools.hpp"
#include "view.hpp"

namespace mapping_detail {
template <class ExecSpace, class Type>
struct MappingAccessibility : std::false_type
{
};

template <typename Type, template <typename ClassType> typename Attribute>
class CheckClassAttributeExistence
{
    template <typename C, typename = Attribute<C>>
    static std::true_type test(int);
    template <typename C>
    static std::false_type test(...);

public:
    static constexpr bool value = decltype(test<Type>(0))::value;
};

/**
 * @brief A helper class to determine if a type is a mapping.
 */
template <typename Type>
class IsMapping
{
    template <typename ClassType>
    using coord_arg_type = typename ClassType::CoordArg;
    template <typename ClassType>
    using coord_result_type = typename ClassType::CoordResult;

    static bool constexpr is_mapping()
    {
        constexpr bool success = CheckClassAttributeExistence<Type, coord_arg_type>::value
                                 && CheckClassAttributeExistence<Type, coord_result_type>::value;
        if constexpr (success) {
            using CoordArg = typename Type::CoordArg;
            using CoordResult = typename Type::CoordResult;
            return std::is_invocable_r_v<CoordResult, Type, CoordArg>;
        }
        return success;
    }

public:
    /// True if the type describes a mapping, false otherwise
    static constexpr bool value = is_mapping();
};

template <typename Type, bool HideError>
class DefinesJacobian
{
    struct IdxTag;
    template <typename ClassType>
    using coord_arg_type = typename ClassType::CoordJacobian;
    template <typename ClassType>
    using jacobian_matrix = decltype(std::declval<ClassType>().jacobian_matrix(
            std::declval<typename ClassType::CoordJacobian>()));
    template <typename ClassType>
    using jacobian_component
            = decltype(std::declval<ClassType>().template jacobian_component<IdxTag, IdxTag>(
                    std::declval<typename ClassType::CoordJacobian>()));
    template <typename ClassType>
    using jacobian = decltype(std::declval<ClassType>().jacobian(
            std::declval<typename ClassType::CoordJacobian>()));

    static bool constexpr has_jacobian_methods()
    {
        if constexpr (!CheckClassAttributeExistence<Type, coord_arg_type>::value) {
            static_assert(HideError, "A Mapping must define the CoordJacobian alias");
            return false;
        };
        if constexpr (!CheckClassAttributeExistence<Type, jacobian_matrix>::value) {
            static_assert(HideError, "A Mapping must define the jacobian_matrix function");
            return false;
        }
        if constexpr (!CheckClassAttributeExistence<Type, jacobian_component>::value) {
            static_assert(HideError, "A Mapping must define the jacobian_component function");
            return false;
        }
        if constexpr (!CheckClassAttributeExistence<Type, jacobian>::value) {
            static_assert(HideError, "A Mapping must define the jacobian function");
            return false;
        }
        return true;
    }

    static bool constexpr has_jacobian()
    {
        static_assert(mapping_detail::IsMapping<Type>::value);
        constexpr bool success = has_jacobian_methods();
        if constexpr (success) {
            using ArgBasisCov = get_covariant_dims_t<ddc::to_type_seq_t<typename Type::CoordArg>>;
            using ResultBasis
                    = get_contravariant_dims_t<ddc::to_type_seq_t<typename Type::CoordResult>>;
            if constexpr (!std::is_same_v<
                                  DTensor<ResultBasis, ArgBasisCov>,
                                  jacobian_matrix<Type>>) {
                static_assert(
                        HideError,
                        "The jacobian_matrix method of a Mapping must take a Coordinate as an "
                        "argument and return a Tensor.");
                return false;
            }
            if constexpr (!std::is_same_v<double, jacobian_component<Type>>) {
                static_assert(
                        HideError,
                        "The jacobian_component method of a Mapping must take a Coordinate as an "
                        "argument and return a double.");
                return false;
            }
            if constexpr (!std::is_same_v<double, jacobian<Type>>) {
                static_assert(
                        HideError,
                        "The jacobian method of a Mapping must take a Coordinate as an argument "
                        "and return a double.");
                return false;
            }
            return true;
        }
        return success;
    }

public:
    /// True if the type describes a mapping with a Jacobian, false otherwise
    static constexpr bool value = has_jacobian();
};

template <typename Type, bool HideError>
class DefinesInvJacobian
{
    struct IdxTag;
    template <typename ClassType>
    using inv_jacobian_matrix = decltype(std::declval<ClassType>().inv_jacobian_matrix(
            std::declval<typename ClassType::CoordJacobian>()));
    template <typename ClassType>
    using inv_jacobian_component
            = decltype(std::declval<ClassType>().template inv_jacobian_component<IdxTag, IdxTag>(
                    std::declval<typename ClassType::CoordJacobian>()));

    static bool constexpr has_inv_jacobian_methods()
    {
        if constexpr (!CheckClassAttributeExistence<Type, inv_jacobian_matrix>::value) {
            static_assert(HideError, "A Mapping must define the inv_jacobian_matrix function");
            return false;
        }
        if constexpr (!CheckClassAttributeExistence<Type, inv_jacobian_component>::value) {
            static_assert(HideError, "A Mapping must define the inv_jacobian_component function");
            return false;
        }
        return true;
    }

    static bool constexpr has_inv_jacobian()
    {
        if constexpr (!DefinesJacobian<Type, HideError>::value) {
            static_assert(
                    HideError,
                    "A Mapping must define its Jacobian before defining its inverse");
            return false;
        }
        static_assert(mapping_detail::IsMapping<Type>::value);
        constexpr bool success = has_inv_jacobian_methods();
        if constexpr (success) {
            using ResultBasisCov
                    = get_covariant_dims_t<ddc::to_type_seq_t<typename Type::CoordResult>>;
            using ArgBasis = get_contravariant_dims_t<ddc::to_type_seq_t<typename Type::CoordArg>>;
            if constexpr (!std::is_same_v<
                                  DTensor<ArgBasis, ResultBasisCov>,
                                  inv_jacobian_matrix<Type>>) {
                static_assert(
                        HideError,
                        "The inv_jacobian_matrix method of a Mapping must take a Coordinate as an "
                        "argument and return a Tensor.");
                return false;
            }
            if constexpr (!std::is_same_v<double, inv_jacobian_component<Type>>) {
                static_assert(
                        HideError,
                        "The inv_jacobian_component method of a Mapping must take a Coordinate as "
                        "an argument and return a double.");
                return false;
            }
            return true;
        }
        return success;
    }

public:
    /// True if the type describes a mapping with an inverse jacobian, false otherwise
    static constexpr bool value = has_inv_jacobian();
};

template <class Mapping>
struct HasOPoint : std::false_type
{
};

template <typename Type>
class IsAnalyticalMapping
{
private:
    template <typename ClassType>
    using inverse_mapping = decltype(&ClassType::get_inverse_mapping);

    static bool constexpr is_analytical_mapping()
    {
        constexpr bool success = CheckClassAttributeExistence<Type, inverse_mapping>::value;
        if constexpr (success) {
            return std::is_invocable_v<inverse_mapping<Type>, Type>;
        }
        return false;
    }

public:
    /// True if the type describes an analytical mapping, false otherwise
    static constexpr bool value = is_analytical_mapping();
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
static constexpr bool is_mapping_v = mapping_detail::IsMapping<Mapping>::value;

template <class Mapping, bool RaiseError = true>
static constexpr bool has_jacobian_v = mapping_detail::DefinesJacobian<Mapping, !RaiseError>::value;

template <class Mapping, bool RaiseError = true>
static constexpr bool has_inv_jacobian_v
        = mapping_detail::DefinesInvJacobian<Mapping, !RaiseError>::value;

/// Indicates that a coordinate change operator is 2D with a curvilinear mapping showing an O-point.
template <class Mapping>
static constexpr bool has_o_point_v
        = mapping_detail::HasOPoint<std::remove_const_t<std::remove_reference_t<Mapping>>>::value;

template <class Mapping>
static constexpr bool is_analytical_mapping_v = mapping_detail::IsAnalyticalMapping<Mapping>::value;

template <class Mapping>
using inverse_mapping_t = decltype(std::declval<Mapping>().get_inverse_mapping());

template <class Mapping>
static constexpr bool has_singular_o_point_inv_jacobian_v
        = mapping_detail::SingularOPointInvJacobian<
                std::remove_const_t<std::remove_reference_t<Mapping>>>::value;
