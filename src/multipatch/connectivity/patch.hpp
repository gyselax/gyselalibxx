// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

/**
 * @brief Base tag for a patch. 
 * @tparam T Variadic template parameter allowing us to use Patch struct 
 * for different numbers of dimensions. 
 */
template <class... T>
struct Patch;

/**
 * @brief Tag for a patch. 
 * 
 * @tparam grid1 Grid on the first logical dimension associated to the patch. 
 * @tparam grid2 Grid on the second logical dimension associated to the patch. 
 * @tparam bsplines_dim1 Bspline dimension defined on the first logical continuous dimension associated to the patch. 
 * @tparam bsplines_dim2 Bspline dimension defined on the second logical continuous dimension associated to the patch. 
 */
template <class grid1, class grid2, class bsplines_dim1, class bsplines_dim2>
struct Patch<grid1, grid2, bsplines_dim1, bsplines_dim2>
{
    /// @brief The number of dimensions of the patch.
    static int constexpr n_dims = 2;

    /// @brief Grid on the first logical dimension.
    using Grid1 = grid1;
    /// @brief Grid on the second logical dimension.
    using Grid2 = grid2;

    /// @brief First continuous dimension.
    using Dim1 = typename grid1::continuous_dimension_type;
    /// @brief Second continuous dimension.
    using Dim2 = typename grid2::continuous_dimension_type;

    /// @brief B-splines defined along the first dimension.
    using BSplines1 = bsplines_dim1;
    /// @brief B-splines defined along the second dimension.
    using BSplines2 = bsplines_dim2;

    /// @brief Coordinate type on the first continuous dimension.
    using Coord1 = ddc::Coordinate<Dim1>;
    /// @brief Coordinate type on the second continuous dimension.
    using Coord2 = ddc::Coordinate<Dim2>;
    /// @brief Coordinate type on the 2D continuous domain.
    using Coord12 = ddc::Coordinate<Dim1, Dim2>;

    /// @brief 1D index of a grid point along the first dimension.
    using Idx1 = ddc::DiscreteElement<Grid1>;
    /// @brief 1D index of a grid point along the second dimension.
    using Idx2 = ddc::DiscreteElement<Grid2>;
    /// @brief 2D index of a grid point along the first and second dimensions.
    using Idx12 = ddc::DiscreteElement<Grid1, Grid2>;

    /// @brief 1D index step between grid points along the first dimension.
    using IdxStep1 = ddc::DiscreteVector<Grid1>;
    /// @brief 1D index step between grid points along the second dimension.
    using IdxStep2 = ddc::DiscreteVector<Grid2>;
    /// @brief 2D index step between grid points along the first and second dimensions.
    using IdxStep12 = ddc::DiscreteVector<Grid1, Grid2>;

    /// @brief Index range of a grids over the first dimension.
    using IdxRange1 = ddc::DiscreteDomain<Grid1>;
    /// @brief Index range of a grids over the second dimension.
    using IdxRange2 = ddc::DiscreteDomain<Grid2>;
    /// @brief Index range of a grids over the first and second dimension.
    using IdxRange12 = ddc::DiscreteDomain<Grid1, Grid2>;

    /// @brief Index range of a grids over the first spline dimension.
    using BSIdxRange1 = ddc::DiscreteDomain<BSplines1>;
    /// @brief Index range of a grids over the second spline dimension.
    using BSIdxRange2 = ddc::DiscreteDomain<BSplines2>;
    /// @brief Index range of a grids over the first and second spline dimension.
    using BSIdxRange12 = ddc::DiscreteDomain<BSplines1, BSplines2>;
};
