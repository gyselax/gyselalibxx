// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

/**
 * @brief An abstract class for initializing a distribution function that does not depend on space.
 */
class IEquilibrium
{
public:
    virtual ~IEquilibrium() = default;

    /**
     * @brief Operator for initializing a distribution function that does not depend on space.
     * @param[in, out] allfequilibrium On input: the uninitialized distribution function.
     *                                 On output: the initialized distribution function.
     * @return The initialized distribution function.
     */
    virtual device_t<DSpanSpVx> operator()(device_t<DSpanSpVx> allfequilibrium) const = 0;
};
