// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "edge.hpp"
#include "edge_transformation.hpp"
#include "single_interface_derivatives_calculator.hpp"
#include "types.hpp"



template <class... DerivCalculatorType>
class SingleInterfaceDerivativesCalculatorCollection
{
    using DerivCalculatorTypeSeq = ddc::detail::TypeSeq<DerivCalculatorType...>;

    using InterfaceTypeSeq
            = ddc::detail::TypeSeq<typename DerivCalculatorType::associated_interface...>;

    static constexpr std::size_t n_calculators = sizeof...(DerivCalculatorType);


    template <class SingleDerivCalculatorType>
    struct is_same_as_associate_correct_type
    {
        using correct_type = SingleInterfaceDerivativesCalculator<
                typename SingleDerivCalculatorType::associated_interface,
                SingleDerivCalculatorType::boundary_condition_patch1,
                SingleDerivCalculatorType::boundary_condition_patch2>;

        static constexpr bool value = std::is_same_v<SingleDerivCalculatorType, correct_type>;
    };

    template <class SingleDerivCalculatorType>
    static constexpr bool is_same_as_associate_correct_type_v
            = is_same_as_associate_correct_type<SingleDerivCalculatorType>::value;


    static_assert(
            (is_same_as_associate_correct_type_v<DerivCalculatorType> && ...),
            "The input parameters should be SingleInterfaceDerivativesCalculator.");


    template <class Interface>
    struct get_deriv_calulator
    {
        using type_HH = SingleInterfaceDerivativesCalculator<
                Interface,
                ddc::BoundCond::HERMITE,
                ddc::BoundCond::HERMITE>;
        using type_HG = SingleInterfaceDerivativesCalculator<
                Interface,
                ddc::BoundCond::HERMITE,
                ddc::BoundCond::GREVILLE>;
        using type_GH = SingleInterfaceDerivativesCalculator<
                Interface,
                ddc::BoundCond::GREVILLE,
                ddc::BoundCond::HERMITE>;
        using type_GG = SingleInterfaceDerivativesCalculator<
                Interface,
                ddc::BoundCond::GREVILLE,
                ddc::BoundCond::GREVILLE>;

        using options_type_seq = ddc::detail::TypeSeq<type_HH, type_HG, type_GH, type_GG>;

        using options_type_seq_minus_correct_option
                = ddc::type_seq_remove_t<options_type_seq, DerivCalculatorTypeSeq>;

        using correct_option_type_seq
                = ddc::type_seq_remove_t<options_type_seq, options_type_seq_minus_correct_option>;

        static_assert(
                ddc::type_seq_size_v<correct_option_type_seq> == 1,
                "One and only one of the derivative calculators in the collection has to be define on the "
                "given interface.");
        using type = ddc::type_seq_element_t<0, correct_option_type_seq>;
    };

    template <class Interface>
    using get_deriv_calulator_t = typename get_deriv_calulator<Interface>::type;



    std::tuple<DerivCalculatorType const&...> m_derivative_calculator_collection;

public:
    SingleInterfaceDerivativesCalculatorCollection(
            DerivCalculatorType const&... derivative_calculators)
        : m_derivative_calculator_collection(derivative_calculators...)
    {
    }


    /**
     * @brief Get a derivative calculator of the collection. 
     * The output cannot be copied. This operator only allows to 
     * get a temporary reference to call one the operator of the 
     * SingleInterfaceDerivativesCalculator class. 
     */
    template <class Interface>
    get_deriv_calulator_t<Interface> const& get() const
    {
        static_assert(
                ddc::in_tags_v<Interface, InterfaceTypeSeq>,
                "No element defined on this Interface in this collection.");

        return std::get<get_deriv_calulator_t<Interface> const&>(
                m_derivative_calculator_collection);
    }
};