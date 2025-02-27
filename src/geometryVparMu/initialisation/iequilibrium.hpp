// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

/**
 * @brief An abstract class for initialising a distribution function in (species,vpar,mu).
 */
class IEquilibrium
{
public:
    virtual ~IEquilibrium() = default;

    /**
     * @brief Operator for initialising an equilibrium distribution function.
     * @param[in, out] allfequilibrium On input: the uninitialized distribution function.
     *                                 On output: the initialised distribution function.
     * @return The initialised equilibrium distribution function.
     */
    virtual DFieldSpVparMu operator()(DFieldSpVparMu allfequilibrium) const = 0;
};
