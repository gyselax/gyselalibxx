// SPDX-License-Identifier: MIT
#pragma once
#include <type_traits>
#include <utility>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"

/**
 * @brief Base tag for a patch.
 * @tparam PatchId A unique identifier for the patch. This provides a compile-time
 * ordering of the patches.
 * @tparam T Variadic template parameter allowing us to use Patch struct
 * for different numbers of dimensions.
 */
template <std::size_t PatchId, class... T>
struct Patch;

/**
 * @brief Tag for a patch.
 *
 * @tparam PatchId A unique identifier for the patch. This provides a compile-time
 * ordering of the patches.
 * @tparam grid1 Grid on the first logical dimension associated to the patch.
 * @tparam grid2 Grid on the second logical dimension associated to the patch.
 * @tparam bsplines_dim1 Bspline dimension defined on the first logical continuous dimension associated to the patch.
 * @tparam bsplines_dim2 Bspline dimension defined on the second logical continuous dimension associated to the patch.
 */
template <std::size_t PatchId, class grid1, class grid2, class bsplines_dim1, class bsplines_dim2>
struct Patch<PatchId, grid1, grid2, bsplines_dim1, bsplines_dim2>
{
    /// @brief The unique identifier of the patch.
    static constexpr std::size_t id = PatchId;

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
    using Coord1 = Coord<Dim1>;
    /// @brief Coordinate type on the second continuous dimension.
    using Coord2 = Coord<Dim2>;
    /// @brief Coordinate type on the 2D continuous dimensions.
    using Coord12 = Coord<Dim1, Dim2>;

    /// @brief 1D index of a grid point along the first dimension.
    using Idx1 = Idx<Grid1>;
    /// @brief 1D index of a grid point along the second dimension.
    using Idx2 = Idx<Grid2>;
    /// @brief 2D index of a grid point along the first and second dimensions.
    using Idx12 = Idx<Grid1, Grid2>;

    /// @brief 1D index step between grid points along the first dimension.
    using IdxStep1 = IdxStep<Grid1>;
    /// @brief 1D index step between grid points along the second dimension.
    using IdxStep2 = IdxStep<Grid2>;
    /// @brief 2D index step between grid points along the first and second dimensions.
    using IdxStep12 = IdxStep<Grid1, Grid2>;

    /// @brief Index range of a grids over the first dimension.
    using IdxRange1 = IdxRange<Grid1>;
    /// @brief Index range of a grids over the second dimension.
    using IdxRange2 = IdxRange<Grid2>;
    /// @brief Index range of a grids over the first and second dimension.
    using IdxRange12 = IdxRange<Grid1, Grid2>;

    /// @brief Index range of a grids over the first spline dimension.
    using IdxRangeBS1 = IdxRange<BSplines1>;
    /// @brief Index range of a grids over the second spline dimension.
    using IdxRangeBS2 = IdxRange<BSplines2>;
    /// @brief Index range of a grids over the first and second spline dimension.
    using IdxRangeBS12 = IdxRange<BSplines1, BSplines2>;
};

namespace detail {
/**
 * @brief Overload selected if the argument points to a Patch or to a class
 * which derives (publicly and unambiguously) from a Patch.
 */
template <std::size_t PatchId, class... T>
std::true_type is_patch_impl(Patch<PatchId, T...> const*);

/// @brief Fallback overload selected for any type which is not a Patch.
std::false_type is_patch_impl(...);
} // namespace detail

/**
 * @brief A type trait to determine whether a type is a Patch, i.e. either an
 * instantiation of Patch or a class deriving from such an instantiation.
 * @tparam PatchType The type being examined.
 */
template <class PatchType>
struct is_patch : decltype(detail::is_patch_impl(std::declval<std::add_pointer_t<PatchType>>()))
{
};

template <class PatchType>
bool constexpr is_patch_v = is_patch<PatchType>::value;

namespace concepts {
template <class PatchType>
concept Patch = is_patch_v<PatchType>;
}; // namespace concepts
