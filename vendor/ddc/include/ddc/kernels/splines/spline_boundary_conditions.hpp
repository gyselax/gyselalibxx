// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include <iostream>
#include <stdexcept>

namespace ddc {
enum class BoundCond {
    // Periodic boundary condition u(1)=u(n)
    PERIODIC,
    // Hermite boundary condition
    HERMITE,
    // Use Greville points instead of conditions on derivative for B-Spline
    // interpolation
    GREVILLE,
    // Natural boundary condition
    NATURAL
};

static inline std::ostream& operator<<(std::ostream& out, ddc::BoundCond const bc)
{
    switch (bc) {
    case ddc::BoundCond::PERIODIC:
        return out << "PERIODIC";
    case ddc::BoundCond::HERMITE:
        return out << "HERMITE";
    case ddc::BoundCond::GREVILLE:
        return out << "GREVILLE";
    case ddc::BoundCond::NATURAL:
        return out << "NATURAL";
    default:
        throw std::runtime_error("ddc::BoundCond not handled");
    }
}

constexpr int n_boundary_equations(ddc::BoundCond const bc, std::size_t const degree)
{
    if (bc == ddc::BoundCond::PERIODIC) {
        return 0;
    } else if (bc == ddc::BoundCond::HERMITE) {
        return degree / 2;
    } else if (bc == ddc::BoundCond::GREVILLE) {
        return 0;
    } else if (bc == ddc::BoundCond::NATURAL) {
        return degree / 2;
    } else {
        throw std::runtime_error("ddc::BoundCond not handled");
    }
}

constexpr int n_user_input(ddc::BoundCond const bc, std::size_t const degree)
{
    if (bc == ddc::BoundCond::PERIODIC) {
        return 0;
    } else if (bc == ddc::BoundCond::HERMITE) {
        return degree / 2;
    } else if (bc == ddc::BoundCond::GREVILLE) {
        return 0;
    } else if (bc == ddc::BoundCond::NATURAL) {
        return 0;
    } else {
        throw std::runtime_error("ddc::BoundCond not handled");
    }
}
} // namespace ddc
