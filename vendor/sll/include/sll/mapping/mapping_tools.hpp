// SPDX-License-Identifier: MIT
#pragma once

#include <array>
#include <type_traits>

#include <sll/view.hpp>

namespace mapping_detail {
template <class ExecSpace, class Type>
struct MappingAccessibility : std::false_type
{
};

/**
 * @brief A helper class to determine if a class has a given attribute (type alias or method).
 * @tparam Type The type whose attributes are being checked.
 * @tparam Attribute A template which returns the expected attribute when given a compatible class.
 */
template <typename Type, template <typename ClassType> typename Attribute>
class CheckClassAttributeExistence
{
private:
    // Class for SFINAE deduction
    template <typename U>
    class check
    {
    };

    // Function that will be chosen if ClassType has a attribute called idx_range with 0 arguments
    template <typename ClassType>
    static char attribute(check<Attribute<ClassType>>*);

    // Function that will be chosen by default
    template <typename ClassType>
    static long attribute(...);

public:
    /// True if the type has the expected attribute, false otherwise
    static constexpr bool has_attribute = (sizeof(attribute<Type>(0)) == sizeof(char));
};

/**
 * @brief A helper class to determine if a type is a 2D mapping.
 */
template <typename Type>
class Is2DMapping
{
    template <typename ClassType>
    using coord_arg_type = typename ClassType::CoordArg;
    template <typename ClassType>
    using coord_result_type = typename ClassType::CoordResult;

    static bool constexpr is_2d_mapping()
    {
        constexpr bool success
                = CheckClassAttributeExistence<Type, coord_arg_type>::has_attribute
                  && CheckClassAttributeExistence<Type, coord_result_type>::has_attribute;
        if constexpr (success) {
            using CoordArg = typename Type::CoordArg;
            using CoordResult = typename Type::CoordResult;
            return std::is_invocable_r_v<CoordResult, Type, CoordArg>;
        }
        return true;
    }

public:
    /// True if the type describes a 2d mapping, false otherwise
    static constexpr bool value = is_2d_mapping();
};

template <typename Type, typename CoordinateType>
class Defines2DJacobian
{
    template <typename ClassType>
    using jacobian_type = decltype(&ClassType::jacobian_matrix);
    template <typename ClassType>
    using jacobian_11 = decltype(&ClassType::jacobian_11);
    template <typename ClassType>
    using jacobian_12 = decltype(&ClassType::jacobian_12);
    template <typename ClassType>
    using jacobian_21 = decltype(&ClassType::jacobian_21);
    template <typename ClassType>
    using jacobian_22 = decltype(&ClassType::jacobian_22);
    template <typename ClassType>
    using jacobian = decltype(&ClassType::jacobian);

    static std::tuple<bool, const char*> constexpr has_2d_jacobian_methods()
    {
        if (!CheckClassAttributeExistence<Type, jacobian_type>::has_attribute) {
            return std::make_tuple(false, "A 2D Mapping must define the jacobian_matrix function");
        }
        if (!CheckClassAttributeExistence<Type, jacobian_11>::has_attribute) {
            return std::make_tuple(false, "A 2D Mapping must define the jacobian_11 function");
        }
        if (!CheckClassAttributeExistence<Type, jacobian_12>::has_attribute) {
            return std::make_tuple(false, "A 2D Mapping must define the jacobian_12 function");
        }
        if (!CheckClassAttributeExistence<Type, jacobian_21>::has_attribute) {
            return std::make_tuple(false, "A 2D Mapping must define the jacobian_21 function");
        }
        if (!CheckClassAttributeExistence<Type, jacobian_22>::has_attribute) {
            return std::make_tuple(false, "A 2D Mapping must define the jacobian_22 function");
        }
        if (!CheckClassAttributeExistence<Type, jacobian>::has_attribute) {
            return std::make_tuple(false, "A 2D Mapping must define the jacobian function");
        }
        return std::make_tuple(true, "");
    }

    static std::tuple<bool, const char*> constexpr has_2d_jacobian()
    {
        constexpr std::tuple<bool, const char*> success = has_2d_jacobian_methods();
        if constexpr (std::get<bool>(success)) {
            if (!std::is_invocable_r_v<
                        void,
                        decltype(&Type::jacobian_matrix),
                        Type,
                        CoordinateType,
                        Matrix_2x2&>) {
                return std::make_tuple(
                        false,
                        "The jacobian_matrix method of a 2D Mapping must take a Coordinate and a "
                        "Matrix_2x2& as an argument and return nothing.");
            }
            if (!std::is_invocable_r_v<
                        double,
                        decltype(&Type::jacobian_11),
                        Type,
                        CoordinateType>) {
                return std::make_tuple(
                        false,
                        "The jacobian_11 method of a 2D Mapping must take a Coordinate as an "
                        "argument and return a double.");
            }
            if (!std::is_invocable_r_v<
                        double,
                        decltype(&Type::jacobian_12),
                        Type,
                        CoordinateType>) {
                return std::make_tuple(
                        false,
                        "The jacobian_12 method of a 2D Mapping must take a Coordinate as an "
                        "argument and return a double.");
            }
            if (!std::is_invocable_r_v<
                        double,
                        decltype(&Type::jacobian_21),
                        Type,
                        CoordinateType>) {
                return std::make_tuple(
                        false,
                        "The jacobian_21 method of a 2D Mapping must take a Coordinate as an "
                        "argument and return a double.");
            }
            if (!std::is_invocable_r_v<
                        double,
                        decltype(&Type::jacobian_22),
                        Type,
                        CoordinateType>) {
                return std::make_tuple(
                        false,
                        "The jacobian_22 method of a 2D Mapping must take a Coordinate as an "
                        "argument and return a double.");
            }
            if (!std::is_invocable_r_v<double, decltype(&Type::jacobian), Type, CoordinateType>) {
                return std::make_tuple(
                        false,
                        "The jacobian method of a 2D Mapping must take a Coordinate as an argument "
                        "and return a double.");
            }
            return std::make_tuple(true, "");
        }
        return success;
    }

    static constexpr std::tuple<bool, const char*> has_2d_jacobian_output = has_2d_jacobian();

public:
    /// True if the type describes a 2d mapping, false otherwise
    static constexpr bool value = std::get<bool>(has_2d_jacobian_output);
    /// A string containing any explanation for why the type is not a jacobian (for debugging)
    static constexpr const char* error_msg = std::get<const char*>(has_2d_jacobian_output);
};

template <typename Type, typename CoordinateType>
class Defines2DInvJacobian
{
    template <typename ClassType>
    using inv_jacobian_type = decltype(&ClassType::inv_jacobian_matrix);
    template <typename ClassType>
    using inv_jacobian_11 = decltype(&ClassType::inv_jacobian_11);
    template <typename ClassType>
    using inv_jacobian_12 = decltype(&ClassType::inv_jacobian_12);
    template <typename ClassType>
    using inv_jacobian_21 = decltype(&ClassType::inv_jacobian_21);
    template <typename ClassType>
    using inv_jacobian_22 = decltype(&ClassType::inv_jacobian_22);

    static std::tuple<bool, const char*> constexpr has_2d_inv_jacobian_methods()
    {
        if (!CheckClassAttributeExistence<Type, inv_jacobian_type>::has_attribute) {
            return std::
                    make_tuple(false, "A 2D Mapping must define the inv_jacobian_matrix function");
        }
        if (!CheckClassAttributeExistence<Type, inv_jacobian_11>::has_attribute) {
            return std::make_tuple(false, "A 2D Mapping must define the inv_jacobian_11 function");
        }
        if (!CheckClassAttributeExistence<Type, inv_jacobian_12>::has_attribute) {
            return std::make_tuple(false, "A 2D Mapping must define the inv_jacobian_12 function");
        }
        if (!CheckClassAttributeExistence<Type, inv_jacobian_21>::has_attribute) {
            return std::make_tuple(false, "A 2D Mapping must define the inv_jacobian_21 function");
        }
        if (!CheckClassAttributeExistence<Type, inv_jacobian_22>::has_attribute) {
            return std::make_tuple(false, "A 2D Mapping must define the inv_jacobian_22 function");
        }
        return std::make_tuple(true, "");
    }

    static std::tuple<bool, const char*> constexpr has_2d_inv_jacobian()
    {
        constexpr std::tuple<bool, const char*> success = has_2d_inv_jacobian_methods();
        if constexpr (std::get<bool>(success)) {
            if (!std::is_invocable_r_v<
                        void,
                        decltype(&Type::inv_jacobian_matrix),
                        Type,
                        CoordinateType,
                        Matrix_2x2&>) {
                return std::make_tuple(
                        false,
                        "The inv_jacobian_matrix method of a 2D Mapping must take a Coordinate and "
                        "a Matrix_2x2& as an argument and return nothing.");
            }
            if (!std::is_invocable_r_v<
                        double,
                        decltype(&Type::inv_jacobian_11),
                        Type,
                        CoordinateType>) {
                return std::make_tuple(
                        false,
                        "The inv_jacobian_11 method of a 2D Mapping must take a Coordinate as an "
                        "argument and return a double.");
            }
            if (!std::is_invocable_r_v<
                        double,
                        decltype(&Type::inv_jacobian_12),
                        Type,
                        CoordinateType>) {
                return std::make_tuple(
                        false,
                        "The inv_jacobian_12 method of a 2D Mapping must take a Coordinate as an "
                        "argument and return a double.");
            }
            if (!std::is_invocable_r_v<
                        double,
                        decltype(&Type::inv_jacobian_21),
                        Type,
                        CoordinateType>) {
                return std::make_tuple(
                        false,
                        "The inv_jacobian_21 method of a 2D Mapping must take a Coordinate as an "
                        "argument and return a double.");
            }
            if (!std::is_invocable_r_v<
                        double,
                        decltype(&Type::inv_jacobian_22),
                        Type,
                        CoordinateType>) {
                return std::make_tuple(
                        false,
                        "The inv_jacobian_22 method of a 2D Mapping must take a Coordinate as an "
                        "argument and return a double.");
            }
            return std::make_tuple(true, "");
        }
        return success;
    }

    static constexpr std::tuple<bool, const char*> has_2d_inv_jacobian_output
            = has_2d_inv_jacobian();

public:
    /// True if the type describes a 2d mapping, false otherwise
    static constexpr bool value = std::get<bool>(has_2d_inv_jacobian_output);
    /// A string containing any explanation for why the type is not an inverse jacobian (for debugging)
    static constexpr const char* error_msg = std::get<const char*>(has_2d_inv_jacobian_output);
};

} // namespace mapping_detail

template <class ExecSpace, class Type>
static constexpr bool is_accessible_v
        = mapping_detail::MappingAccessibility<ExecSpace, Type>::value;

template <class Mapping>
static constexpr bool is_2d_mapping_v = mapping_detail::Is2DMapping<Mapping>::value;

template <class Mapping, class CoordinateType>
static constexpr bool has_2d_jacobian_v
        = mapping_detail::Defines2DJacobian<Mapping, CoordinateType>::value;

template <class Mapping, class CoordinateType>
static constexpr bool has_2d_inv_jacobian_v
        = mapping_detail::Defines2DInvJacobian<Mapping, CoordinateType>::value;
