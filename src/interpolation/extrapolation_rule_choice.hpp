// SPDX-License-Identifier: MIT
#pragma once
#include "constant_identity_interpolation_extrapolation_rule.hpp"
#include "i_interpolation.hpp"

namespace details {

/// Primary template: Rule has no nested `type` alias template, so it is assumed to
/// already be a concrete, usable extrapolation rule class - use it as-is.
template <class Rule, class CoeffGrid, class DataType, class = void>
struct ExtrapolationRuleResolver
{
    using type = Rule;
};

/// Specialisation selected via SFINAE when Rule exposes `Rule::type<CoeffGrid, DataType>`
/// (i.e. Rule is a self-describing tag such as ExtrapolationRule::Periodic) - resolve
/// through it to obtain the concrete extrapolation rule class.
template <class Rule, class Basis, class DataType>
struct ExtrapolationRuleResolver<
        Rule,
        CoeffGrid,
        DataType,
        std::void_t<typename Rule::template type<CoeffGrid, DataType>>>
{
    using type = typename Rule::template type<CoeffGrid, DataType>;
};

} // namespace details

/// @brief Resolve Rule (either a self-describing tag or a concrete extrapolation rule
/// class) to the concrete extrapolation rule class to use for a given CoeffGrid/DataType.
template <class Rule, class Basis, class DataType>
using extrapolation_rule_t =
        typename details::ExtrapolationRuleResolver<Rule, CoeffGrid, DataType>::type;

/// @brief True if get_extrapolation<Rule, CoeffGrid, DataType, ...> can build the rule
/// with no extra information: the resolved rule type is default-constructible, or Rule
/// is the ExtrapolationRule::Constant tag (which get_extrapolation knows how to build
/// from the boundary of the discrete space).
template <class Rule, class Basis, class DataType>
constexpr bool is_extrapolation_rule_auto_constructible_v
        = std::is_default_constructible_v<extrapolation_rule_t<
                  Rule,
                  CoeffGrid,
                  DataType>> || std::is_same_v<Rule, ExtrapolationRule::Constant>;

/**
 * @brief Initialise the extrapolation rule.
 *
 * Initialise the extrapolation rule using the default constructor if available. For a
 * constant extrapolation, initialise the extrapolation from the selected extremity.
 */
template <class Rule, class CoeffGrid, class DataType, class Basis = CoeffGrid>
extrapolation_rule_t<Rule, CoeffGrid, DataType> get_extrapolation(Extremity extremity)
{
    using RuleType = extrapolation_rule_t<Rule, CoeffGrid, DataType>;
    if constexpr (std::is_default_constructible_v<RuleType>) {
        return RuleType();
    } else if constexpr (std::is_same_v<Rule, ExtrapolationRule::Constant>) {
        if constexpr (is_spline_basis_v<Basis>) {
            if (extremity == Extremity::FRONT) {
                return RuleType(ddc::discrete_space<Basis>().rmin());
            } else {
                return RuleType(ddc::discrete_space<Basis>().rmax());
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
