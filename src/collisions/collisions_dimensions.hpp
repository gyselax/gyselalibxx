// SPDX-License-Identifier: MIT
#pragma once
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"

/**
 * @brief A namespace to collect classes which are necessary to create Fields with the
 * correct number of dimensions to be compatible with Koliop.
 */
namespace collisions_dimensions {
/**
 * Create dimensions to act as the radial/poloidal/toroidal tag for Koliop even if these dims don't
 * exist in the simulation
 */

/// Class from which fake dimensions inherit. These fake dimensions are inserted in order to match Koliop interface
struct InternalSpoofGrid
{
};

/// Fake radial dimension to be used if there is no radial dimension in the simulation
struct InternalSpoofGridR : InternalSpoofGrid
{
};

/// Fake poloidal dimension to be used if there is no poloidal dimension in the simulation
struct InternalSpoofGridTheta : InternalSpoofGrid
{
};

/// Check if a dimension is spoofed but is not present in the actual simulation
template <class Grid1D>
inline constexpr bool is_spoofed_dim_v = std::is_base_of_v<InternalSpoofGrid, Grid1D>;

/**
 * Class to get the type of the radial dimension from a field containing a radial profile.
 * @tparam FieldMem The type of the field containing the radial profile.
 */
template <class FieldMem>
struct ExtractRDim
{
    static_assert(!std::is_same_v<FieldMem, FieldMem>, "Unrecognised radial profile type");
};

/**
 * Class to get the type of the poloidal dimension from a field containing a profile on the poloidal plane.
 * @tparam FieldMem The type of the field containing the profile on the poloidal plane.
 * @tparam GridR The tag for the discrete radial dimension.
 */
template <class FieldMem, class GridR>
struct ExtractThetaDim
{
    static_assert(!std::is_same_v<FieldMem, FieldMem>, "Unrecognised poloidal profile type");
};

/**
 * @brief Get the index range for a specific grid from a multi-D index range.
 * If the dimension is spoofed and does not appear in the multi-D index range then an index
 * range which only iterates over the index(0) is returned.
 *
 * @tparam Grid The tag for the specific grid.
 * @param idx_range The multi-D index range.
 * @returns The index range for the specific grid.
 */
template <class Grid1D, class IdxRangeFDistrib>
inline IdxRange<Grid1D> get_1d_idx_range(IdxRangeFDistrib idx_range)
{
    if constexpr (is_spoofed_dim_v<Grid1D>) {
        return IdxRange<Grid1D>(Idx<Grid1D> {0}, IdxStep<Grid1D> {1});
    } else {
        return ddc::select<Grid1D>(idx_range);
    }
}

/**
 * @brief Get the index range for specific grid dimensions from a multi-D index range.
 *
 * @tparam Grid1D The tags for the specific grid dimensions.
 * @param idx_range The multi-D index range.
 * @returns The index range for the specific grid.
 */
template <class... Grid1D, class IdxRangeFDistrib>
inline IdxRange<Grid1D...> get_idx_range(IdxRangeFDistrib idx_range)
{
    return IdxRange<Grid1D...>(get_1d_idx_range<Grid1D>(idx_range)...);
}

/// If radial profile is stored in a double then the grid tag must be spoofed.
template <>
struct ExtractRDim<double>
{
    using type = InternalSpoofGridR;
};

/// If radial profile is stored in a 1D chunk then the grid tag is extracted.
template <class GridR, class Layout>
struct ExtractRDim<
        ConstField<double, IdxRange<GridR>, Layout, Kokkos::DefaultExecutionSpace::memory_space>>
{
    using type = GridR;
};

/// If radial profile is stored in a chunk then information about why the chunk has the wrong type is provided
template <class ElementType, class IdxRange, class Layout, class MemSpace>
struct ExtractRDim<Field<ElementType, IdxRange, Layout, MemSpace>>
{
    static_assert(
            std::is_same_v<ElementType, const double>,
            "The radial profile should be a double and should not be modifiable. Please use a "
            "constant view.");
    static_assert(
            std::is_same_v<MemSpace, Kokkos::DefaultExecutionSpace::memory_space>,
            "The radial profile should be provided on the GPU");
    static_assert(
            (ddc::type_seq_size_v<ddc::to_type_seq_t<IdxRange>>) > 1,
            "The radial profile should not be defined on more than 1 dimensions.");
};

/// If the profile on the poloidal plane is stored in a double then the grid tag must be spoofed.
template <>
struct ExtractThetaDim<double, InternalSpoofGridR>
{
    using type = InternalSpoofGridTheta;
};

/// If the profile on the poloidal plane is stored in a chunk on radial values then the grid tag must be spoofed.
template <class GridR, class Layout>
struct ExtractThetaDim<
        ConstField<double, IdxRange<GridR>, Layout, Kokkos::DefaultExecutionSpace::memory_space>,
        GridR>
{
    using type = InternalSpoofGridTheta;
};

/// If the profile on the poloidal plane is stored in a chunk on poloidal values then the grid tag is extracted.
template <class GridR, class GridTheta, class Layout>
struct ExtractThetaDim<
        ConstField<
                double,
                IdxRange<GridTheta>,
                Layout,
                Kokkos::DefaultExecutionSpace::memory_space>,
        GridR>
{
    using type = GridTheta;
};

/// If the profile on the poloidal plane is stored in a 2D chunk then the grid tag is extracted.
template <class GridR, class GridTheta, class Layout>
struct ExtractThetaDim<
        ConstField<
                double,
                IdxRange<GridTheta, GridR>,
                Layout,
                Kokkos::DefaultExecutionSpace::memory_space>,
        GridR>
{
    using type = GridTheta;
};

/// If poloidal profile is stored in a chunk then information about why the chunk has the wrong type is provided
template <class GridR, class ElementType, class IdxRange, class Layout, class MemSpace>
struct ExtractThetaDim<Field<ElementType, IdxRange, Layout, MemSpace>, GridR>
{
    static_assert(
            std::is_same_v<ElementType, const double>,
            "The poloidal profile should be a double and should not be modifiable. Please use a "
            "constant view.");
    static_assert(
            std::is_same_v<MemSpace, Kokkos::DefaultExecutionSpace::memory_space>,
            "The poloidal profile should be provided on the GPU");
    static_assert(
            (ddc::type_seq_size_v<ddc::to_type_seq_t<IdxRange>>) > 2,
            "The poloidal profile should not be defined on more than 2 dimensions.");
};

/**
 * A helper function to check if a grid is found at the correct location in a set of
 * ordered grids.
 *
 * @tparam Grid1D The grid we are checking.
 * @tparam OrderedGrids A ddc::TypeSeq containing an ordered set of grids.
 *
 * @param idx The index where the dimension is expected to be found
 *
 * @returns True if the grid is at the expected location or if the grid is spoofed.
 *          False otherwise.
 */
template <class Grid1D, class OrderedGrids>
constexpr bool check_dimension_location(int& idx)
{
    if constexpr (!is_spoofed_dim_v<Grid1D>) {
        bool success = ddc::type_seq_rank_v<Grid1D, OrderedGrids> == idx;
        idx += 1;
        return success;
    } else {
        return true;
    }
}

/**
 * A helper function to check if the final grids in an ordered set of grids match
 * the order of the provided expected grids.
 *
 * @tparam OrderedGrids A ddc::TypeSeq containing an ordered set of grids.
 * @tparam ExpectedGrids The expected grids in the order that they are expected.
 *
 * @returns True if the grids are ordered correctly, false otherwise.
 */
template <class OrderedGrids, class... ExpectedGrids>
constexpr bool order_of_last_grids()
{
    int n_real_dims = ((is_spoofed_dim_v<ExpectedGrids> ? 0 : 1) + ...);
    int idx(ddc::type_seq_size_v<OrderedGrids> - n_real_dims);
    return ((check_dimension_location<ExpectedGrids, OrderedGrids>(idx)) && ...);
}

/**
 * A helper function to check if the first grids in an ordered set of grids match
 * the order of the provided expected grids.
 *
 * @tparam OrderedGrids A ddc::TypeSeq containing an ordered set of grids.
 * @tparam ExpectedGrids The expected grids in the order that they are expected.
 *
 * @returns True if the grids are ordered correctly, false otherwise.
 */
template <class OrderedGrids, class... ExpectedGrids>
constexpr bool order_of_first_grids()
{
    int idx(0);
    return ((check_dimension_location<ExpectedGrids, OrderedGrids>(idx)) && ...);
}

} // namespace collisions_dimensions
