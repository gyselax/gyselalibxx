// SPDX-License-Identifier: MIT
#pragma once
#include "constant_identity_interpolation_extrapolation_rule.hpp"
#include "i_interpolation.hpp"

namespace details {

/// Primary template: Rule has no nested `type` alias template, so it is assumed to
/// already be a concrete, usable extrapolation rule class - use it as-is.
template <class Rule, class DataType, class Basis, class IdxRangeCoeff, class = void>
struct ExtrapolationRuleResolver
{
    using type = Rule;
};

/// Specialisation selected via SFINAE when Rule exposes `Rule::type<Basis, DataType>`
/// (i.e. Rule is a self-describing tag such as ExtrapolationRule::Periodic) - resolve
/// through it to obtain the concrete extrapolation rule class.
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

} // namespace details

/// @brief True if get_extrapolation<Rule, CoeffGrid, DataType, ...> can build the rule
/// with no extra information: the resolved rule type is default-constructible, or Rule
/// is the ExtrapolationRule::Constant tag (which get_extrapolation knows how to build
/// from the boundary of the discrete space).
template <class Rule, class CoeffGrid, class DataType, class Basis>
constexpr bool is_extrapolation_rule_auto_constructible_v
        = (std::is_default_constructible_v<Rule>)
          || (details::is_ddc_constant_extrapolation_rule_v<Rule>)
          || (std::is_same_v<
                  Rule,
                  ConstantIdentityInterpolationExtrapolationRule<CoeffGrid, DataType>>);

/// @brief Resolve Rule (either a self-describing tag or a concrete extrapolation rule
/// class) to the concrete extrapolation rule class to use for a given Basis/DataType.
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

/**
 * @brief Initialise the extrapolation rule.
 *
 * Initialise the extrapolation rule using the default constructor if available. For a
 * constant extrapolation, initialise the extrapolation from the selected extremity.
 */
template <class Rule, class CoeffGrid, class DataType, class Basis>
Rule get_extrapolation(Extremity extremity)
{
    if constexpr (std::is_default_constructible_v<Rule>) {
        return Rule();
    } else if constexpr (details::is_ddc_constant_extrapolation_rule_v<Rule>) {
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
