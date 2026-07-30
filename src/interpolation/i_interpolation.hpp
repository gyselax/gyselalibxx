// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/kernels/splines.hpp>

#include "constant_identity_interpolation_extrapolation_rule.hpp"
#include "geometry_descriptors.hpp"
#include "i_interpolation_builder.hpp"
#include "i_interpolation_evaluator.hpp"
#include "lagrange_basis_non_uniform.hpp"
#include "lagrange_basis_uniform.hpp"
#include "type_seq_tools.hpp"

namespace concepts {

/**
 * @brief A concept describing an interpolation object that owns a matching builder–evaluator pair.
 *
 * An Interpolation bundles a builder (which converts function values on the
 * interpolation mesh into interpolation coefficients) and an evaluator (which
 * reconstructs function values from those coefficients) into a single owning object.
 *
 * Concrete models include SplineInterpolator and LagrangeInterpolator.
 *
 * @tparam T The interpolation type to check.
 */
template <typename T>
concept Interpolation = requires
{
    typename T::BuilderType;
    typename T::EvaluatorType;
    {
        T::rank()
        } -> std::same_as<std::size_t>;
}
&&concepts::InterpolationBuilder<typename T::BuilderType>&& concepts::InterpolationEvaluator<
        typename T::
                EvaluatorType> && (InterpolationBuilderTraits<typename T::BuilderType>::rank() == InterpolationEvaluatorTraits<typename T::EvaluatorType>::rank())
        && requires(T const& t)
{
    {
        t.get_builder()
        } -> std::same_as<typename T::BuilderType const&>;
    {
        t.get_evaluator()
        } -> std::same_as<typename T::EvaluatorType const&>;
};

/**
 * @brief A concept describing a 1D interpolation object.
 *
 * @tparam T The interpolation type to check.
 */
template <typename T>
concept Interpolation1D = Interpolation<T> && InterpolationBuilder1D<typename T::BuilderType>;

} // namespace concepts

/**
 * @brief A type trait that is true when Basis is a Lagrange basis type.
 *
 * Evaluates to true for both UniformLagrangeBasis and NonUniformLagrangeBasis
 * instantiations, false for all other types.
 *
 * @tparam Basis The basis type to test.
 */
template <class Basis>
constexpr bool is_lagrange_basis_v
        = is_uniform_lagrange_basis_v<Basis> || is_non_uniform_lagrange_basis_v<Basis>;

/**
 * @brief A type trait that is true when Basis is a B-spline basis type.
 *
 * Evaluates to true for both ddc::UniformBSplines and ddc::NonUniformBSplines
 * instantiations, false for all other types.
 *
 * @tparam Basis The basis type to test.
 */
template <class Basis>
constexpr bool is_spline_basis_v
        = ddc::is_uniform_bsplines_v<Basis> || ddc::is_non_uniform_bsplines_v<Basis>;

/**
 * @brief A namespace containing tag types describing how a function is extrapolated
 * outside the interpolation domain.
 *
 * Each tag exposes a type alias template that resolves, given the coefficient grid
 * and data type of the interpolator, to the concrete extrapolation rule class to use.
 * This keeps the set of extrapolation rules open for extension: user code may define
 * its own tag following the same pattern, or bypass tags entirely and pass an
 * already-concrete extrapolation rule class directly to SplineInterpolator or
 * LagrangeInterpolator (see extrapolation_rule_choice.hpp for the generic resolution
 * mechanism).
 */
namespace ExtrapolationRule {

/**
 * @brief Tag selecting periodic extrapolation.
 *
 * The value at a point outside the domain is taken as the value at the equivalent
 * point inside the domain.
 */
struct OneSidedPeriodic
{
    /// @brief The concrete extrapolation rule class for a given CoeffGrid/DataType.
    template <class Basis, class DataType>
    using type = ddc::PeriodicExtrapolationRule<typename CoeffGrid::continuous_dimension_type>;
};

/// @brief Convenience pairing of ExtrapolationRule::OneSidedPeriodic for both boundaries.
using Periodic = ddc::detail::TypeSeq<OneSidedPeriodic, OneSidedPeriodic>;

/// @brief Tag selecting null extrapolation: the function evaluates to zero outside the domain.
struct NullValue
{
    /// @brief The concrete extrapolation rule class for a given CoeffGrid/DataType.
    template <class Basis, class DataType>
    using type = ddc::NullExtrapolationRule;
};

/// @brief Convenience pairing of ExtrapolationRule::NullValue for both boundaries.
using Null_Null = ddc::detail::TypeSeq<NullValue, NullValue>;

/**
 * @brief Tag selecting constant extrapolation.
 *
 * The function is clamped to the value at the nearest boundary point.
 */
struct Constant
{
    /// @brief The concrete extrapolation rule class for a given CoeffGrid/DataType.
    template <class Basis, class DataType>
    using type = std::conditional_t<
            is_spline_basis_v<Basis>,
            ddc::ConstantExtrapolationRule<typename Basis::continuous_dimension_type>,
            ConstantIdentityInterpolationExtrapolationRule<Basis, DataType>>;
};

/// @brief Convenience pairing of ExtrapolationRule::Constant for both boundaries.
using Constant_Constant = ddc::detail::TypeSeq<Constant, Constant>;

} // namespace ExtrapolationRule
