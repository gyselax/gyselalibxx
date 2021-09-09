#pragma once

enum class BoundCond {
    // Periodic boundary condition u(1)=u(n)
    PERIODIC,
    // Hermite boundary condition
    HERMITE,
    // Use Greville points instead of conditions on derivative for B-Spline
    // interpolation
    GREVILLE,
};

static std::ostream& operator<<(std::ostream& out, BoundCond bc)
{
    switch (bc) {
    case BoundCond::PERIODIC:
        return out << "PERIODIC";
    case BoundCond::HERMITE:
        return out << "HERMITE";
    case BoundCond::GREVILLE:
        return out << "GREVILLE";
    default:
        std::exit(1);
    }
}
