#pragma once

enum class BoundaryCondition {
    // Periodic boundary condition u(1)=u(n)
    sll_p_periodic,
    // Hermite boundary condition
    sll_p_hermite,
    // Use Greville points instead of conditions on derivative for B-Spline
    // interpolation
    sll_p_greville,
};
