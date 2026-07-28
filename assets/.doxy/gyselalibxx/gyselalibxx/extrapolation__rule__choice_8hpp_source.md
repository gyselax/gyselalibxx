

# File extrapolation\_rule\_choice.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**extrapolation\_rule\_choice.hpp**](extrapolation__rule__choice_8hpp.md)

[Go to the documentation of this file](extrapolation__rule__choice_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include "constant_identity_interpolation_extrapolation_rule.hpp"
#include "i_interpolation.hpp"

namespace details {

template <class Rule, class CoeffGrid, class DataType, class = void>
struct ExtrapolationRuleResolver
{
    using type = Rule;
};

template <class Rule, class CoeffGrid, class DataType>
struct ExtrapolationRuleResolver<
        Rule,
        CoeffGrid,
        DataType,
        std::void_t<typename Rule::template type<CoeffGrid, DataType>>>
{
    using type = typename Rule::template type<CoeffGrid, DataType>;
};

} // namespace details

template <class Rule, class CoeffGrid, class DataType>
using extrapolation_rule_t =
        typename details::ExtrapolationRuleResolver<Rule, CoeffGrid, DataType>::type;

template <class Rule, class CoeffGrid, class DataType>
constexpr bool is_extrapolation_rule_auto_constructible_v
        = std::is_default_constructible_v<extrapolation_rule_t<
                  Rule,
                  CoeffGrid,
                  DataType>> || std::is_same_v<Rule, ExtrapolationRule::Constant>;

template <class Rule, class CoeffGrid, class DataType, class Basis = CoeffGrid>
extrapolation_rule_t<Rule, CoeffGrid, DataType> get_extrapolation(Extremity extremity)
{
    using RuleType = extrapolation_rule_t<Rule, CoeffGrid, DataType>;
    if constexpr (std::is_default_constructible_v<RuleType>) {
        return RuleType();
    } else if constexpr (std::is_same_v<Rule, ExtrapolationRule::Constant>) {
        if constexpr (is_spline_basis_v<Basis>) {
            if (extremity == Extremity::FRONT) {
                return RuleType(ddc::discrete_space<CoeffGrid>().rmin());
            } else {
                return RuleType(ddc::discrete_space<CoeffGrid>().rmax());
            }
        } else {
            if (extremity == Extremity::FRONT) {
                return RuleType(ddc::discrete_space<Basis>().full_domain().front());
            } else {
                return RuleType(ddc::discrete_space<Basis>().full_domain().back());
            }
        }
    } else {
        static_assert("Extrapolation rule initialisation cannot be deduced from type.");
    }
}
```


