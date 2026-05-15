// SPDX-License-Identifier: MIT
#pragma once
#include "constant_identity_interpolation_extrapolation_rule.hpp"
#include "i_interpolation.hpp"

namespace details {

// Get the class describing the interpolation rule
template <ExtrapolationRule Rule, class CoeffGrid, class DataType>
struct GetExtrapolationRuleClass;

template <class CoeffGrid, class DataType>
struct GetExtrapolationRuleClass<PERIODIC, CoeffGrid, DataType>
{
    using type = ddc::PeriodicExtrapolationRule<typename CoeffGrid::continuous_dimension_type>;
};

template <class CoeffGrid, class DataType>
struct GetExtrapolationRuleClass<NULL_VALUE, CoeffGrid, DataType>
{
    using type = ddc::NullExtrapolationRule;
};

template <class CoeffGrid, class DataType>
struct GetExtrapolationRuleClass<CONSTANT, CoeffGrid, DataType>
{
    using type = std::conditional_t<
            is_spline_basis_v<CoeffGrid>,
            ddc::ConstantExtrapolationRule<typename CoeffGrid::continuous_dimension_type>,
            ConstantIdentityInterpolationExtrapolationRule<CoeffGrid, DataType>>;
};

} // namespace details

template <ExtrapolationRule Rule, class CoeffGrid, class DataType>
using extrapolation_rule_t = details::GetExtrapolationRuleClass<Rule, CoeffGrid, DataType>::type;

template <ExtrapolationRule Rule, class CoeffGrid, class DataType, class Basis = CoeffGrid>
extrapolation_rule_t<Rule, CoeffGrid, DataType> get_extrapolation(Extremity extremity)
{
    if constexpr (Rule == CONSTANT) {
        if constexpr (is_spline_basis_v<Basis>) {
            if (extremity == Extremity::FRONT) {
                return extrapolation_rule_t<Rule, CoeffGrid, DataType>(
                        ddc::discrete_space<CoeffGrid>().rmin());
            } else {
                return extrapolation_rule_t<Rule, CoeffGrid, DataType>(
                        ddc::discrete_space<CoeffGrid>().rmax());
            }
        } else {
            if (extremity == Extremity::FRONT) {
                return extrapolation_rule_t<Rule, CoeffGrid, DataType>(
                        ddc::discrete_space<Basis>().full_domain().front());
            } else {
                return extrapolation_rule_t<Rule, CoeffGrid, DataType>(
                        ddc::discrete_space<Basis>().full_domain().back());
            }
        }
    } else {
        return extrapolation_rule_t<Rule, CoeffGrid, DataType>();
    }
}
