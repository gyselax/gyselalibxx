// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "interface_derivative_coefficients.hpp"

template <class T>
inline constexpr bool enable_interface_derivative_coefficients_collection = false;

template <class T>
inline constexpr bool is_interface_derivative_coefficients_collection_v
        = enable_interface_derivative_coefficients_collection<
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
 * @see InterfaceDerivCoeffs. 
 */
template <class... Interfaces>
class InterfaceDerivCoeffsCollection
{
    using DerivCalculatorTypeSeq = ddc::detail::TypeSeq<InterfaceDerivCoeffs<Interfaces>...>;

    using InterfaceTypeSeq = ddc::detail::TypeSeq<Interfaces...>;

    std::tuple<InterfaceDerivCoeffs<Interfaces> const&...> m_deriv_coeffs_collection;

public:
    /**
     * @brief Instantiate a InterfaceDerivCoeffsCollection 
     * from a list of interface derivative calculators. 
     *  
     * @param deriv_coeffs Interface derivative calculators. 
     */
    explicit InterfaceDerivCoeffsCollection(InterfaceDerivCoeffs<Interfaces> const&... deriv_coeffs)
        : m_deriv_coeffs_collection(deriv_coeffs...)
    {
    }


    /**
     * @brief Get a derivative calculator of the collection. 
     * The output cannot be copied. This operator only allows to 
     * get a temporary reference to call one of the operators of the 
     * InterfaceDerivCoeffs class. 
     * 
     * @tparam Interface The interface where the required interface derivative 
     * calculator is defined. 
     * 
     * @return The required interface derivative calculator as a constant reference. 
     */
    template <class Interface>
    InterfaceDerivCoeffs<Interface> const& get() const
    {
        static_assert(
                ddc::in_tags_v<Interface, InterfaceTypeSeq>,
                "No element defined on this Interface in this collection.");

        return std::get<InterfaceDerivCoeffs<Interface> const&>(m_deriv_coeffs_collection);
    }
};


// To help the template deduction.
template <class... DerivCalculatorType>
InterfaceDerivCoeffsCollection(DerivCalculatorType const&... deriv_coeffs)
        -> InterfaceDerivCoeffsCollection<typename DerivCalculatorType::associated_interface...>;


template <class... DerivCalculatorType>
inline constexpr bool enable_interface_derivative_coefficients_collection<
        InterfaceDerivCoeffsCollection<DerivCalculatorType...>> = true;
