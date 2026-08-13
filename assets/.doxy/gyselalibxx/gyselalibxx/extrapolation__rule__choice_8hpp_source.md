

# File extrapolation\_rule\_choice.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**extrapolation\_rule\_choice.hpp**](extrapolation__rule__choice_8hpp.md)

[Go to the documentation of this file](extrapolation__rule__choice_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include "constant_identity_interpolation_extrapolation_rule.hpp"
#include "i_interpolation.hpp"

namespace details {

template <class Rule, class DataType, class Basis, class IdxRangeCoeff, class = void>
struct ExtrapolationRuleResolver
{
    using type = Rule;
};

template <class Rule, class DataType, class Basis, class IdxRangeCoeff>
struct ExtrapolationRuleResolver<
        Rule,
        DataType,
        Basis,
        IdxRangeCoeff,
        std::void_t<typename Rule::template type<DataType, Basis, IdxRangeCoeff>>>
{
    using type = typename Rule::template type<DataType, Basis, IdxRangeCoeff>;
};

template <class Rule>
static constexpr bool is_ddc_constant_extrapolation_rule_v = false;

template <class... CDim>
static constexpr bool
        is_ddc_constant_extrapolation_rule_v<ddc::ConstantExtrapolationRule<CDim...>> = true;

template <class ExtrapRule>
struct ConstantExtrapRuleRank;

template <class... CDim>
struct ConstantExtrapRuleRank<ddc::ConstantExtrapolationRule<CDim...>>
{
    static constexpr std::size_t rank = sizeof...(CDim);
};

} // namespace details

template <class Rule, class CoeffGrid, class DataType, class Basis>
constexpr bool is_extrapolation_rule_auto_constructible_v
        = (std::is_default_constructible_v<Rule>)
          || (details::is_ddc_constant_extrapolation_rule_v<Rule>)
          || (std::is_same_v<
                  Rule,
                  ConstantIdentityInterpolationExtrapolationRule<CoeffGrid, DataType>>);

template <class Rule, class DataType, class Basis, class IdxRangeCoeff>
using extrapolation_rule_t = ddc::detail::TypeSeq<
        typename details::ExtrapolationRuleResolver<
                ddc::type_seq_element_t<0, Rule>,
                DataType,
                Basis,
                IdxRangeCoeff>::type,
        typename details::ExtrapolationRuleResolver<
                ddc::type_seq_element_t<1, Rule>,
                DataType,
                Basis,
                IdxRangeCoeff>::type>;

template <class Rule, class DataType, class CDim, class IdxRangeCoeff, class IdxRangeBasis>
Rule get_extrapolation(Extremity extremity)
{
    if constexpr (std::is_default_constructible_v<Rule>) {
        return Rule();
    }

    using Basis = find_grid_t<CDim, ddc::to_type_seq_t<IdxRangeBasis>>;
    if constexpr (details::is_ddc_constant_extrapolation_rule_v<Rule>) {
        static constexpr std::size_t rank = details::ConstantExtrapRuleRank<Rule>::rank;
        if constexpr (rank == 1) {
            if (extremity == Extremity::FRONT) {
                return Rule(ddc::discrete_space<Basis>().rmin());
            } else {
                return Rule(ddc::discrete_space<Basis>().rmax());
            }
        } else {
            static_assert(rank == 2);
            using NIBasis = ddc::type_seq_element_t<
                    0,
                    ddc::to_type_seq_t<ddc::remove_dims_of_t<IdxRangeBasis, Basis>>>;
            if (extremity == Extremity::FRONT) {
                return Rule(
                        ddc::discrete_space<Basis>().rmin(),
                        ddc::discrete_space<NIBasis>().rmin(),
                        ddc::discrete_space<NIBasis>().rmax());
            } else {
                return Rule(
                        ddc::discrete_space<Basis>().rmax(),
                        ddc::discrete_space<NIBasis>().rmin(),
                        ddc::discrete_space<NIBasis>().rmax());
            }
        }
    }
    using CoeffGrid = find_grid_t<CDim, ddc::to_type_seq_t<IdxRangeCoeff>>;
    if constexpr (std::is_same_v<
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


