

# File interface\_derivative\_coefficients\_collection.hpp

[**File List**](files.md) **>** [**interface\_derivatives**](dir_d1bd52a3e76a422151eefdcc4e15c189.md) **>** [**interface\_derivative\_coefficients\_collection.hpp**](interface__derivative__coefficients__collection_8hpp.md)

[Go to the documentation of this file](interface__derivative__coefficients__collection_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "interface_derivative_coefficients.hpp"

template <class T>
inline constexpr bool enable_interface_derivative_coefficients_collection = false;

template <class T>
inline constexpr bool is_single_derivative_calculator_collection_v
        = enable_interface_derivative_coefficients_collection<
                std::remove_const_t<std::remove_reference_t<T>>>;


template <class... Interfaces>
class InterfaceDerivCoeffsCollection
{
    using DerivCalculatorTypeSeq = ddc::detail::TypeSeq<InterfaceDerivCoeffs<Interfaces>...>;

    using InterfaceTypeSeq = ddc::detail::TypeSeq<Interfaces...>;

    std::tuple<InterfaceDerivCoeffs<Interfaces> const&...> m_derivative_calculator_collection;

public:
    explicit InterfaceDerivCoeffsCollection(
            InterfaceDerivCoeffs<Interfaces> const&... derivative_calculators)
        : m_derivative_calculator_collection(derivative_calculators...)
    {
    }


    template <class Interface>
    InterfaceDerivCoeffs<Interface> const& get() const
    {
        static_assert(
                ddc::in_tags_v<Interface, InterfaceTypeSeq>,
                "No element defined on this Interface in this collection.");

        return std::get<InterfaceDerivCoeffs<Interface> const&>(m_derivative_calculator_collection);
    }
};


// To help the template deduction.
template <class... DerivCalculatorType>
InterfaceDerivCoeffsCollection(DerivCalculatorType const&... derivative_calculators)
        -> InterfaceDerivCoeffsCollection<typename DerivCalculatorType::associated_interface...>;


template <class... DerivCalculatorType>
inline constexpr bool enable_interface_derivative_coefficients_collection<
        InterfaceDerivCoeffsCollection<DerivCalculatorType...>> = true;
```


