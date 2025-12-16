
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

/** 
 * @brief Convert a coordinate on a dimension to another coordinate
 * on another dimension. The scalar value does not change. 
 */
template <class DimOut, class DimIn>
inline Coord<DimOut> constexpr convert_dim(Coord<DimIn> const& input_coord)
{
    return Coord<DimOut> {double(input_coord)};
}

/** 
 * @brief Convert a vector on coordinate on a dimension to a vector on a coordinate
 * on another dimension. The scalar values do not change. 
 */
template <class DimOut, class DimIn>
std::vector<Coord<DimOut>> const convert_dim(std::vector<Coord<DimIn>> const& input_vec)
{
    std::vector<Coord<DimOut>> output_vec;
    for (double pt : input_vec) {
        output_vec.push_back(Coord<DimOut>(pt));
    }
    return output_vec;
}


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
 *  @brief Get interpolation points from the break points by placing 
 * the interpolation points on the break points and adding one on the
 * right boundary cell at 1/3 of the cell and another on the left 
 * boundary cell at 2/3 of the cell. 
 */
template <class CoordType>
std::vector<CoordType> get_interpolation_points_add_one_on_left_and_one_on_right(
        std::vector<CoordType> const& break_points)
{
    return get_interpolation_points_add_one_on_left(
            get_interpolation_points_add_one_on_right(break_points));
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

/**
 * @brief Get an index range slice containing two indices: 
 *   * the first index of the given index range, 
 *   * another index at first index + number of elements in the given index range.
 */
template <class Grid1D>
IdxRangeSlice<Grid1D> get_bound_idx_range_slice(IdxRange<Grid1D> const idx_range)
{
    return IdxRangeSlice<Grid1D>(idx_range.front(), IdxStep<Grid1D>(2), idx_range.extents()-1);
}
