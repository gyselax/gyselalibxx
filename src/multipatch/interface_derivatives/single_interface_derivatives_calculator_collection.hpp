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

    static_assert(
            (is_single_derivative_calculator_v<DerivCalculatorType> && ...),
            "The input parameters should be SingleInterfaceDerivativesCalculator.");
    
    template <class Interface>
    using get_deriv_calulator_t =ddc::type_seq_element_t<ddc::type_seq_rank_v<Interface, InterfaceTypeSeq>, DerivCalculatorTypeSeq>;



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