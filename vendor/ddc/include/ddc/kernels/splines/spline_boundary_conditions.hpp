// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include <iostream>
#include <stdexcept>

namespace ddc {

/** @brief An enum representing a spline boundary condition. Please refer to
 * Emily Bourne's thesis (https://www.theses.fr/2022AIXM0412.pdf)
 */
enum class BoundCond {
    PERIODIC, ///< Periodic boundary condition u(1)=u(n)
    HERMITE, ///< Hermite boundary condition
    GREVILLE, ///< Use Greville points instead of conditions on derivative for B-Spline interpolation
};

/**
 * @brief Prints a boundary condition in a std::ostream.
 *
 * @param out The stream in which the boundary condition is printed.
 * @param degree The boundary condition.
 *
 * @return The stream in which the boundary condition is printed.
 **/
static inline std::ostream& operator<<(std::ostream& out, ddc::BoundCond const bc)
{
    switch (bc) {
    case ddc::BoundCond::PERIODIC:
        return out << "PERIODIC";
    case ddc::BoundCond::HERMITE:
        return out << "HERMITE";
    case ddc::BoundCond::GREVILLE:
        return out << "GREVILLE";
    default:
        throw std::runtime_error("ddc::BoundCond not handled");
    }
}

/**
 * @brief Return the number of equations needed to describe a given boundary condition.
 *
 * @param bc The boundary condition.
 * @param degree The degree of the spline.
 *
 * @return The number of equations.
 **/
constexpr int n_boundary_equations(ddc::BoundCond const bc, std::size_t const degree)
{
    if (bc == ddc::BoundCond::PERIODIC) {
        return 0;
    } else if (bc == ddc::BoundCond::HERMITE) {
        return degree / 2;
    } else if (bc == ddc::BoundCond::GREVILLE) {
        return 0;
    } else {
        throw std::runtime_error("ddc::BoundCond not handled");
    }
}
} // namespace ddc
