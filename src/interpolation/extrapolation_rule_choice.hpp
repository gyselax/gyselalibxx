// SPDX-License-Identifier: MIT
#pragma once
#include "constant_identity_interpolation_extrapolation_rule.hpp"
#include "i_interpolation.hpp"

namespace details {

// Get the class describing the interpolation rule
template <ExtrapolationRule Rule, class CoeffGrid>
struct GetExtrapolationRuleClass;

template <class CoeffGrid>
struct GetExtrapolationRuleClass<PERIODIC, CoeffGrid>
{
    using type = ddc::PeriodicExtrapolationRule<typename CoeffGrid::continuous_dimension_type>;
};

template <class CoeffGrid>
struct GetExtrapolationRuleClass<NULL_VALUE, CoeffGrid>
{
    using type = ddc::NullExtrapolationRule;
};

template <class CoeffGrid>
struct GetExtrapolationRuleClass<CONSTANT, CoeffGrid>
{
    using type = std::conditional_t<
            is_spline_basis_v<CoeffGrid>,
            ddc::ConstantExtrapolationRule<typename CoeffGrid::continuous_dimension_type>,
            ConstantIdentityInterpolationExtrapolationRule<CoeffGrid>>;
};

} // namespace details

template <ExtrapolationRule Rule, class CoeffGrid>
using extrapolation_rule_t = details::GetExtrapolationRuleClass<Rule, CoeffGrid>::type;

template <ExtrapolationRule Rule, class CoeffGrid, class Basis = CoeffGrid>
extrapolation_rule_t<Rule, CoeffGrid> get_extrapolation(Extremity extremity)
{
    if constexpr (Rule == CONSTANT) {
        if constexpr (is_spline_basis_v<Basis>) {
            if (extremity == Extremity::FRONT) {
                return extrapolation_rule_t<Rule, CoeffGrid>(
                        ddc::discrete_space<CoeffGrid>().rmin());
            } else {
                return extrapolation_rule_t<Rule, CoeffGrid>(
                        ddc::discrete_space<CoeffGrid>().rmax());
            }
        } else {
            if (extremity == Extremity::FRONT) {
                return extrapolation_rule_t<Rule, CoeffGrid>(
                        ddc::discrete_space<Basis>().full_domain().front());
            } else {
                return extrapolation_rule_t<Rule, CoeffGrid>(
                        ddc::discrete_space<Basis>().full_domain().back());
            }
        }
    } else {
        return extrapolation_rule_t<Rule, CoeffGrid>();
    }
}
