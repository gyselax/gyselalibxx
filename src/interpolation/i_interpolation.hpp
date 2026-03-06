// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/kernels/splines.hpp>

#include "geometry_descriptors.hpp"
#include "i_interpolation_builder.hpp"
#include "i_interpolation_evaluator.hpp"
#include "lagrange_basis_non_uniform.hpp"
#include "lagrange_basis_uniform.hpp"
#include "type_seq_tools.hpp"

enum ExtrapolationRule { PERIODIC, NULL_VALUE, CONSTANT };

namespace concepts {

template <typename T>
concept Interpolation = requires
{
    typename T::BuilderType;
    typename T::EvaluatorType;
}
&&concepts::InterpolationBuilder<typename T::BuilderType>&& concepts::InterpolationEvaluator<
        typename T::EvaluatorType>&& requires(T const& t)
{
    {
        t.get_builder()
        } -> std::same_as<typename T::BuilderType const&>;
    {
        t.get_evaluator()
        } -> std::same_as<typename T::EvaluatorType const&>;
};

} // namespace concepts

template <class Basis>
constexpr bool is_lagrange_basis_v
        = is_uniform_lagrange_basis_v<Basis> || is_non_uniform_lagrange_basis_v<Basis>;

template <class Basis>
constexpr bool is_spline_basis_v
        = ddc::is_uniform_bsplines_v<Basis> || ddc::is_non_uniform_bsplines_v<Basis>;
