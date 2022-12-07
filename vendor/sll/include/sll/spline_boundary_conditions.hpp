#pragma once

#include <iostream>
#include <stdexcept>

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

static inline std::ostream& operator<<(std::ostream& out, BoundCond const bc)
{
    switch (bc) {
    case BoundCond::PERIODIC:
        return out << "PERIODIC";
    case BoundCond::HERMITE:
        return out << "HERMITE";
    case BoundCond::GREVILLE:
        return out << "GREVILLE";
    case BoundCond::NATURAL:
        return out << "NATURAL";
    default:
        throw std::runtime_error("BoundCond not handled");
    }
}

constexpr int n_boundary_equations(BoundCond const bc, std::size_t const degree)
{
    if (bc == BoundCond::PERIODIC) {
        return 0;
    } else if (bc == BoundCond::HERMITE) {
        return degree / 2;
    } else if (bc == BoundCond::GREVILLE) {
        return 0;
    } else if (bc == BoundCond::NATURAL) {
        return degree / 2;
    } else {
        throw std::runtime_error("BoundCond not handled");
    }
}

constexpr int n_user_input(BoundCond const bc, std::size_t const degree)
{
    if (bc == BoundCond::PERIODIC) {
        return 0;
    } else if (bc == BoundCond::HERMITE) {
        return degree / 2;
    } else if (bc == BoundCond::GREVILLE) {
        return 0;
    } else if (bc == BoundCond::NATURAL) {
        return 0;
    } else {
        throw std::runtime_error("BoundCond not handled");
    }
}
