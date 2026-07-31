

# File extrapolation\_rule\_choice.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**extrapolation\_rule\_choice.hpp**](extrapolation__rule__choice_8hpp.md)

[Go to the documentation of this file](extrapolation__rule__choice_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include "constant_identity_interpolation_extrapolation_rule.hpp"
#include "i_interpolation.hpp"

namespace details {

template <class Rule, class CoeffGrid, class DataType, class Basis, class = void>
struct ExtrapolationRuleResolver
{
    using type = Rule;
};

template <class Rule, class CoeffGrid, class DataType, class Basis>
struct ExtrapolationRuleResolver<
        Rule,
        CoeffGrid,
        DataType,
        Basis,
        std::void_t<typename Rule::template type<CoeffGrid, DataType, Basis>>>
{
    using type = typename Rule::template type<CoeffGrid, DataType, Basis>;
};

} // namespace details

template <class Rule, class CoeffGrid, class DataType, class Basis>
constexpr bool is_extrapolation_rule_auto_constructible_v
        = (std::is_default_constructible_v<Rule>)
          || (std::is_same_v<
                  Rule,
                  ddc::ConstantExtrapolationRule<typename CoeffGrid::continuous_dimension_type>>)
          || (std::is_same_v<
                  Rule,
                  ConstantIdentityInterpolationExtrapolationRule<CoeffGrid, DataType>>);

template <class Rule, class CoeffGrid, class DataType, class Basis>
using extrapolation_rule_t = ddc::detail::TypeSeq<
        typename details::ExtrapolationRuleResolver<
                ddc::type_seq_element_t<0, Rule>,
                CoeffGrid,
                DataType,
                Basis>::type,
        typename details::ExtrapolationRuleResolver<
                ddc::type_seq_element_t<1, Rule>,
                CoeffGrid,
                DataType,
                Basis>::type>;

template <class Rule, class CoeffGrid, class DataType, class Basis>
Rule get_extrapolation(Extremity extremity)
{
    if constexpr (std::is_default_constructible_v<Rule>) {
        return Rule();
    } else if constexpr (std::is_same_v<
                                 Rule,
                                 ddc::ConstantExtrapolationRule<
                                         typename Basis::continuous_dimension_type>>) {
        if (extremity == Extremity::FRONT) {
            return Rule(ddc::discrete_space<Basis>().rmin());
        } else {
            return Rule(ddc::discrete_space<Basis>().rmax());
        }
    } else if constexpr (
            std::is_same_v<
                    Rule,
                    ConstantIdentityInterpolationExtrapolationRule<CoeffGrid, DataType>>) {
        if (extremity == Extremity::FRONT) {
            return Rule(ddc::discrete_space<Basis>().full_domain().front());
        } else {
            return Rule(ddc::discrete_space<Basis>().full_domain().back());
        }
    } else {
        static_assert("Extrapolation rule initialisation cannot be deduced from type.");
    }
}
```


