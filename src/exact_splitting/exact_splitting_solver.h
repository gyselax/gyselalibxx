#ifndef EXACT_SPLITTING_SOLVER_H
#define EXACT_SPLITTING_SOLVER_H

#include <array>

void exact_splitting_solver(
    const std::array<double, 3>& magnetic_field,
    double t,
    std::array<double, 3>& yl,
    std::array<double, 3>& y2,
    std::array<double, 3>& y3,
    std::array<double, 3>& yr,
    const std::array<int, 3>& order
);

#endif // EXACT_SPLITTING_SOLVER_H
