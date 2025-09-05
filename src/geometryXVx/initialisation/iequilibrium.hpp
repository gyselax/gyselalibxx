// SPDX-License-Identifier: MIT

#pragma once

#include <memory>

#include <paraconf.h>

#include "geometry.hpp"
#include "paraconfpp.hpp"

/**
 * @brief An abstract class for initialising a distribution function that does not depend on space.
 */
class IEquilibrium
{
public:
    virtual ~IEquilibrium() = default;

    /**
     * @brief Operator for initialising a distribution function that does not depend on space.
     * @param[in, out] allfequilibrium On input: the uninitialized distribution function.
     *                                 On output: the initialised distribution function.
     * @return The initialised distribution function.
     */
    virtual DFieldSpVx operator()(DFieldSpVx allfequilibrium) const = 0;

    /**
     * @brief Determine the chosen equilibrium method from a YAML input file.
     *
     * Use the Algorithm.equilibrium key in a YAML input file to select
     * the chosen equilibrium method and call the init_from_input method
     * for that method.
     *
     * @param[in] idx_range_kinsp Index range for the kinetic species
     * @param[in] yaml_input_file YAML input file
     * @return A pointer to an equilibrium operator. A pointer is used for OOP.
     */
    static std::unique_ptr<IEquilibrium> init_from_input(
            IdxRangeSp idx_range_kinsp,
            PC_tree_t const& yaml_input_file);
};
