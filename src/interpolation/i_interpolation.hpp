// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/kernels/splines.hpp>

#include "geometry_descriptors.hpp"
#include "i_interpolation_builder.hpp"
#include "i_interpolation_evaluator.hpp"
#include "lagrange_basis_non_uniform.hpp"
#include "lagrange_basis_uniform.hpp"
#include "type_seq_tools.hpp"

/**
 * @brief An enum describing how a function is extrapolated outside the interpolation domain.
 *
 * - PERIODIC : the function is assumed to be periodic. The value at a point outside
 *   the domain is taken as the value at the equivalent point inside the domain.
 * - NULL_VALUE : the function evaluates to zero outside the domain.
 * - CONSTANT : the function is clamped to the value at the nearest boundary point.
 */
enum ExtrapolationRule { PERIODIC, NULL_VALUE, CONSTANT };

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
