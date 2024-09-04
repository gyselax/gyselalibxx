// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

/**
 * @brief An abstract class for initializing a distribution function in (species,vpar,mu,r,theta).
 */
class IEquilibrium
{
public:
    virtual ~IEquilibrium() = default;

    /**
     * @brief Operator for initializing an equilibrium distribution function.
     * @param[in, out] allfequilibrium On input: the uninitialized distribution function.
     *                                 On output: the initialized distribution function.
     * @return The initialized equilibrium distribution function.
     */
    virtual DFieldSpV2DTor2D operator()(DFieldSpV2DTor2D allfequilibrium) const = 0;
};
