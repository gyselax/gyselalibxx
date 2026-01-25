// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

/**
 * @brief An abstract class for initialising the fields.
 * The fields only depend on spatial dimensions.
 */
class IFieldInitialisation
{
public:
    virtual ~IFieldInitialisation() = default;

    /**
     * @brief Operator for initialising the fields that only depend on space.
     *
     * @param[in, out] fields On input: the uninitialized fields.
     *                                 On output: the initialised fields.
     * @return The initialised fields.
     */

    virtual DFieldX operator()(DFieldX magnetic_field_x, DFieldX magnetic_field_y, DFieldX magnetic_field_z, DFieldX pressure) const = 0;
};
