
// SPDX-License-Identifier: MIT
#pragma once 
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>


/**
 *  @brief Get interpolation points from the break points by placing 
 * the interpolation points on the break points and adding one on the
 * left boundary cell at 2/3 of the cell. 
 */
template <class CoordType>
std::vector<CoordType> get_interpolation_points_add_one_on_left(
        std::vector<CoordType> const& break_points)
{
    CoordType additional_point(break_points[0] * 2. / 3. + break_points[1] * 1. / 3.);
    std::vector<CoordType> interpolation_points(break_points);
    interpolation_points.insert(interpolation_points.begin() + 1, additional_point);
    return interpolation_points;
}

/**
 *  @brief Get interpolation points from the break points by placing 
 * the interpolation points on the break points and adding one on the
 * right boundary cell at 1/3 of the cell. 
 */
template <class CoordType>
std::vector<CoordType> get_interpolation_points_add_one_on_right(
        std::vector<CoordType> const& break_points)
{
    int n_bpoints = break_points.size();
    CoordType additional_point(
            break_points[n_bpoints - 1] * 2. / 3. + break_points[n_bpoints - 2] * 1. / 3.);
    std::vector<CoordType> interpolation_points(break_points);
    interpolation_points.insert(interpolation_points.end() - 1, additional_point);
    return interpolation_points;
}

/**
 * @brief Fill in a vector of points for the equivalent global mesh 
 * by conserving the same order of the given points.
 */
template <class CoordTypeG, class CoordTypeP>
void fill_in(std::vector<CoordTypeG>& points_global, std::vector<CoordTypeP> const& points_patch)
{
    for (CoordTypeP pt : points_patch) {
        points_global.push_back(CoordTypeG {double(pt)});
    }
}

/**
 * @brief Fill in a vector of points for the equivalent global mesh
 *  by reversing the order of the given points.
 */
template <class CoordTypeG, class CoordTypeP>
void fill_in_reverse(
        std::vector<CoordTypeG>& points_global,
        std::vector<CoordTypeP> const& points_patch)
{
    std::size_t const n_pt = points_patch.size();
    CoordTypeP const max = points_patch[n_pt - 1];
    CoordTypeP const min = points_patch[0];
    for (int i(0); i < n_pt; ++i) {
        points_global.push_back(CoordTypeG {double(min + max - points_patch[n_pt - 1 - i])});
    }
}
