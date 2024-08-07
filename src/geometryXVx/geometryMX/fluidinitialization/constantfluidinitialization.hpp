// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

#include "ifluidinitialization.hpp"
/**
 * @brief A class that initializes a fluid species with constant moments.
 */
class ConstantFluidInitialization : public IFluidInitialization
{
    DFieldMemSpMom m_moments_alloc;

public:
    /**
     * @brief Creates an instance of the ConstantFluidInitialization class.
     * @param[in] moments The fluid moments the fluid species should be initialized with. 
     */
    ConstantFluidInitialization(host_t<DConstFieldSpMom> moments);

    ~ConstantFluidInitialization() override = default;

    /**
     * @brief Initializes the fluid species with a constant moments.
     * @param[inout] fluid_moments On input: a span referencing an uninitialized fluid species described through its moments.
     *                             On output: a span referencing a the fluid species initialized with constant moments.
     * @return A span referencing the initialized fluid species.
     */
    DFieldSpMomX operator()(DFieldSpMomX const fluid_moments) const override;
};
