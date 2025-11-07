// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "single_interface_derivatives_calculator.hpp"

template <class T>
inline constexpr bool enable_single_derivative_calculator_collection = false;

template <class T>
inline constexpr bool is_single_derivative_calculator_collection_v
        = enable_single_derivative_calculator_collection<
                std::remove_const_t<std::remove_reference_t<T>>>;


/**
 * @brief A class to store a collection of interface derivative calculators templated 
 * on the interfaces. 
 * 
 * The class stores a constant reference of interface derivative calculators. 
 * It should not be use to copy the elements outside of the class but it should be 
 * use to access the operators of the stored interface derivative calculators. 
 * 
 * @tparam Interfaces Types of interface that defined the interface derivative calculators. 
 * 
 * @warning For each interface, only one interface derivative calculator should be defined.
 * 
 * @see SingleInterfaceDerivativesCalculator. 
 */
template <class... Interfaces>
class SingleInterfaceDerivativesCalculatorCollection
{
    using DerivCalculatorTypeSeq
            = ddc::detail::TypeSeq<SingleInterfaceDerivativesCalculator<Interfaces>...>;

    using InterfaceTypeSeq = ddc::detail::TypeSeq<Interfaces...>;

    std::tuple<SingleInterfaceDerivativesCalculator<Interfaces> const&...>
            m_derivative_calculator_collection;

public:
    /**
     * @brief Instantiate a SingleInterfaceDerivativesCalculatorCollection 
     * from a list of interface derivative calculators. 
     *  
     * @param derivative_calculators Interface derivative calculators. 
     */
    explicit SingleInterfaceDerivativesCalculatorCollection(
            SingleInterfaceDerivativesCalculator<Interfaces> const&... derivative_calculators)
        : m_derivative_calculator_collection(derivative_calculators...)
    {
    }


    /**
     * @brief Get a derivative calculator of the collection. 
     * The output cannot be copied. This operator only allows to 
     * get a temporary reference to call one of the operators of the 
     * SingleInterfaceDerivativesCalculator class. 
     * 
     * @tparam Interface The interface where the required interface derivative 
     * calculator is defined. 
     * 
     * @return The required interface derivative calculator as a constant reference. 
     */
    template <class Interface>
    SingleInterfaceDerivativesCalculator<Interface> const& get() const
    {
        static_assert(
                ddc::in_tags_v<Interface, InterfaceTypeSeq>,
                "No element defined on this Interface in this collection.");

        return std::get<SingleInterfaceDerivativesCalculator<Interface> const&>(
                m_derivative_calculator_collection);
    }
};


// To help the template deduction.
template <class... DerivCalculatorType>
SingleInterfaceDerivativesCalculatorCollection(DerivCalculatorType const&... derivative_calculators)
        -> SingleInterfaceDerivativesCalculatorCollection<
                typename DerivCalculatorType::associated_interface...>;


template <class... DerivCalculatorType>
inline constexpr bool enable_single_derivative_calculator_collection<
        SingleInterfaceDerivativesCalculatorCollection<DerivCalculatorType...>> = true;
