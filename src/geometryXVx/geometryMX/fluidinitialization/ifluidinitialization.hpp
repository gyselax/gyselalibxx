// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

/**
 * @brief An abstract class that allows for initializing a fluid species.
 */
class IFluidInitialization
{
public:
    virtual ~IFluidInitialization() = default;

    /**
     * @brief Operator for initializing a neutral species.
     * @param[in, out] fluid_moments On input: the uninitialized fluid species.
     *                                 On output: the initialized fluid species.
     * @return A span referencing the initialized fluid species.
     */
    virtual DFieldSpMomX operator()(DFieldSpMomX fluid_moments) const = 0;
};
