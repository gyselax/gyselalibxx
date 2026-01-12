

# File single\_interface\_derivatives\_calculator\_collection.hpp

[**File List**](files.md) **>** [**interface\_derivatives**](dir_d1bd52a3e76a422151eefdcc4e15c189.md) **>** [**single\_interface\_derivatives\_calculator\_collection.hpp**](single__interface__derivatives__calculator__collection_8hpp.md)

[Go to the documentation of this file](single__interface__derivatives__calculator__collection_8hpp.md)


```C++
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


template <class... Interfaces>
class SingleInterfaceDerivativesCalculatorCollection
{
    using DerivCalculatorTypeSeq
            = ddc::detail::TypeSeq<SingleInterfaceDerivativesCalculator<Interfaces>...>;

    using InterfaceTypeSeq = ddc::detail::TypeSeq<Interfaces...>;

    std::tuple<SingleInterfaceDerivativesCalculator<Interfaces> const&...>
            m_derivative_calculator_collection;

public:
    explicit SingleInterfaceDerivativesCalculatorCollection(
            SingleInterfaceDerivativesCalculator<Interfaces> const&... derivative_calculators)
        : m_derivative_calculator_collection(derivative_calculators...)
    {
    }


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
```


