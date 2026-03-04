#pragma once
#include "i_interpolation.hpp"

namespace details {

// Get the class describing the interpolation rule
template <ExtrapolationRule Rule, class Basis>
struct GetExtrapolationRuleClass;

template <class Basis>
struct GetExtrapolationRuleClass<PERIODIC, Basis>
{
    using type = ddc::PeriodicExtrapolationRule<typename Basis::continuous_dimension_type>;
};

template <class Basis>
struct GetExtrapolationRuleClass<NULL_VALUE, Basis>
{
    using type = ddc::NullExtrapolationRule;
};

template <class Basis>
struct GetExtrapolationRuleClass<CONSTANT, Basis>
{
    using type = std::conditional_t<
            is_spline_basis_v<Basis>,
            ddc::ConstantExtrapolationRule<typename Basis::continuous_dimension_type>,
            ConstantIdentityInterpolationExtrapolationRule<Basis>>;
};

} // namespace details

template <ExtrapolationRule Rule, class Basis>
using extrapolation_rule_t = GetExtrapolationRuleClass<Rule, Basis>::type;

template <ExtrapolationRule Rule, class Basis>
extrapolation_rule_t<Rule, Basis> get_extrapolation(Extremity extremity)
{
    if constexpr (Rule == CONSTANT) {
        if (extremity == Extremity::FRONT) {
            return extrapolation_rule_t<Rule, Basis>(ddc::discrete_space<Basis>().rmin());
        } else {
            return extrapolation_rule_t<Rule, Basis>(ddc::discrete_space<Basis>().rmax());
        }
    } else {
        return extrapolation_rule_t<Rule, Basis>();
    }
}
