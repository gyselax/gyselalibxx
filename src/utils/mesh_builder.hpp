// SPDX-License-Identifier: MIT
#pragma once
#include <cstdlib>
#include <ctime>
#include <vector>

#include "ddc_aliases.hpp"

template <class Dim, class Grid1D>
std::vector<Coord<Dim>> build_uniform_break_points(
        Coord<Dim> min,
        Coord<Dim> max,
        IdxStep<Grid1D> n_cells)
{
    static_assert(std::is_same_v<typename Grid1D::continuous_dimension_type, Dim>);
    std::vector<Coord<Dim>> break_points(n_cells + 1);
    double const delta((max - min) / double(n_cells));

    break_points[0] = min;
    for (int i(1); i < n_cells; ++i) {
        break_points[i] = min + i * delta;
    }
    break_points[n_cells] = max;
    return break_points;
}

template <class Dim, class Grid1D>
std::vector<Coord<Dim>> build_random_non_uniform_break_points(
        Coord<Dim> min,
        Coord<Dim> max,
        IdxStep<Grid1D> n_cells,
        double const non_uniformity = 1.)
{
    static_assert(std::is_same_v<typename Grid1D::continuous_dimension_type, Dim>);

    std::srand(std::time(nullptr)); // Seed with random value (the time)

    std::vector<Coord<Dim>> break_points(n_cells + 1);
    double const delta((max - min) / double(n_cells));

    break_points[0] = min;
    for (int i(1); i < n_cells; ++i) {
        double const random_perturbation = double(rand()) / RAND_MAX - 0.5;
        break_points[i] = min + (i + random_perturbation * non_uniformity) * delta;
    }
    break_points[n_cells] = max;
    return break_points;
}
