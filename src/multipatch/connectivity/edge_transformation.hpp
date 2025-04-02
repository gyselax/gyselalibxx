// SPDX-License-Identifier: MIT

#pragma once

#include <numeric>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "edge.hpp"
#include "interface.hpp"


/**
 * @brief Struct to identify the type of an index range on the current 
 * patch and the index range type on the target patch. 
 */
template <class Interface, class PatchOrDim>
struct IdxRangeCurrentAndTarget;

/// @brief Specification for current patch being the Patch1 of the given Interface.
template <class Interface>
struct IdxRangeCurrentAndTarget<Interface, typename Interface::Edge1::associated_patch>
{
    /// @brief Grid parallel to Edge1 of the Interface.
    using EdgeGrid1 = typename Interface::Edge1::parallel_grid;
    /// @brief Grid parallel to Edge2 of the Interface.
    using EdgeGrid2 = typename Interface::Edge2::parallel_grid;

    /// @brief IdxRange parallel to the Edge on the current patch.
    using IdxRangeCurrent = IdxRange<EdgeGrid1>;
    /// @brief IdxRange parallel to the Edge on the target patch.
    using IdxRangeTarget = IdxRange<EdgeGrid2>;
};

/// @brief Specification for current patch being the Patch2 of the given Interface.
template <class Interface>
struct IdxRangeCurrentAndTarget<Interface, typename Interface::Edge2::associated_patch>
{
    /// @brief Grid parallel to Edge1 of the Interface.
    using EdgeGrid1 = typename Interface::Edge1::parallel_grid;
    /// @brief Grid parallel to Edge2 of the Interface.
    using EdgeGrid2 = typename Interface::Edge2::parallel_grid;

    /// @brief IdxRange parallel to the Edge on the current patch.
    using IdxRangeCurrent = IdxRange<EdgeGrid2>;
    /// @brief IdxRange parallel to the Edge on the target patch.
    using IdxRangeTarget = IdxRange<EdgeGrid1>;
};

/**
 * @brief Specification for current continuous dimension being the 
 * continuous dimension of the grid on Edge1 of the given Interface. 
 */
template <class Interface>
struct IdxRangeCurrentAndTarget<
        Interface,
        typename Interface::Edge1::parallel_grid::continuous_dimension_type>
{
    /// @brief Grid parallel to Edge1 of the Interface.
    using EdgeGrid1 = typename Interface::Edge1::parallel_grid;
    /// @brief Grid parallel to Edge2 of the Interface.
    using EdgeGrid2 = typename Interface::Edge2::parallel_grid;

    /// @brief IdxRange parallel to the Edge on the current patch.
    using IdxRangeCurrent = IdxRange<EdgeGrid1>;
    /// @brief IdxRange parallel to the Edge on the target patch.
    using IdxRangeTarget = IdxRange<EdgeGrid2>;
};

/**
 * @brief Specification for current continuous dimension being the 
 * continuous dimension of the grid on Edge2 of the given Interface. 
 */
template <class Interface>
struct IdxRangeCurrentAndTarget<
        Interface,
        typename Interface::Edge2::parallel_grid::continuous_dimension_type>
{
    /// @brief Grid parallel to Edge1 of the Interface.
    using EdgeGrid1 = typename Interface::Edge1::parallel_grid;
    /// @brief Grid parallel to Edge2 of the Interface.
    using EdgeGrid2 = typename Interface::Edge2::parallel_grid;

    /// @brief IdxRange parallel to the Edge on the current patch.
    using IdxRangeCurrent = IdxRange<EdgeGrid2>;
    /// @brief IdxRange parallel to the Edge on the target patch.
    using IdxRangeTarget = IdxRange<EdgeGrid1>;
};


/**
 * @brief Get equivalent coordinate on the target patch. 
 * 
 * @tparam PatchOrDim Type for a continuous dimension or a patch. 
 * @tparam Interface Interface type given in EdgeTransformation class. 
 * @tparam CurrentDim The continuous dimension where the given coordinate is defined. 
 * 
 * @param current_coord Coordinate on the current patch. 
 * @param idx_range_patch_1 Index range on the Edge1 of the Interface. 
 * @param idx_range_patch_2 Index range on the Edge2 of the Interface. 
 * 
 * @return Equivalent coordinate on the target patch. 
 */
template <typename PatchOrDim, class Interface, class CurrentDim>
Coord<std::conditional_t<
        std::is_same_v<
                CurrentDim,
                typename Interface::Edge1::parallel_grid::continuous_dimension_type>,
        typename Interface::Edge2::parallel_grid::continuous_dimension_type,
        typename Interface::Edge1::parallel_grid::continuous_dimension_type>>
get_equivalent_target_coordinate(
        Coord<CurrentDim> const& current_coord,
        IdxRange<typename Interface::Edge1::parallel_grid> const& idx_range_patch_1,
        IdxRange<typename Interface::Edge2::parallel_grid> const& idx_range_patch_2)
{
    using EdgeGrid1 = typename Interface::Edge1::parallel_grid;
    using EdgeGrid2 = typename Interface::Edge2::parallel_grid;

    using EdgeDim1 = typename EdgeGrid1::continuous_dimension_type;
    using EdgeDim2 = typename EdgeGrid2::continuous_dimension_type;

    // The other continuous dimension
    using TargetDim = std::conditional_t<std::is_same_v<CurrentDim, EdgeDim1>, EdgeDim2, EdgeDim1>;

    using IdxRangeCurrent =
            typename IdxRangeCurrentAndTarget<Interface, PatchOrDim>::IdxRangeCurrent;
    using IdxRangeTarget = typename IdxRangeCurrentAndTarget<Interface, PatchOrDim>::IdxRangeTarget;

    // Gem_idx_range_patch_1min and length on each 1D index ranges:
    IdxRange<EdgeGrid1, EdgeGrid2> const combined_idx_range(idx_range_patch_1, idx_range_patch_2);
    IdxRangeCurrent const current_idx_range(combined_idx_range);
    IdxRangeTarget const target_idx_range(combined_idx_range);

    Coord<CurrentDim> const current_min = ddc::coordinate(current_idx_range.front());
    double const current_length = ddcHelper::total_interval_length(current_idx_range);

    Coord<TargetDim> const target_min = ddc::coordinate(target_idx_range.front());
    double target_length = ddcHelper::total_interval_length(target_idx_range);

    double rescale_x = (current_coord - current_min) / current_length * target_length;

    bool constexpr orientations_agree = Interface::orientations_agree;
    if constexpr (!orientations_agree) {
        rescale_x = target_length - rescale_x;
    }
    if constexpr (TargetDim::PERIODIC) {
        rescale_x = std::fmod(rescale_x, target_length);
    }
    return target_min + rescale_x;
}


/**
 * @brief Transform a coordinate or an index from one edge to the one on the other edge.
 * 
 * According to the orientation of the interface, we compute the equivalent coordinate
 * * if True,
 *  @f$ \\ t \mapsto min_2 + \frac{t - min_1}{max_1 - min_1}(max_2 - min_2) @f$
 * 
 * * if False, 
 *  @f$ \\ t \mapsto max_2 - \frac{t - min_1}{max_1 - min_1}(max_2 - min_2) @f$
 * 
 * where @f$ min_i @f$ and @f$ max_i @f$ are the minimum and maximum 
 * coordinates of the edge @f$ i @f$. 
 * 
 * For the indices, we look for an equivalent index corresponding to a coordinate 
 * equivalent to the coordinate of the initial index. 
 * 
 * @tparam Interface The Interface type where we want to compute the transformation.
*/
template <class Interface>
class EdgeTransformation
{
    static_assert(
            !std::is_same_v<
                    typename Interface::Edge1,
                    OutsideEdge> && !std::is_same_v<typename Interface::Edge2, OutsideEdge>,
            "The interface cannot be an interface with the outside domain.");

    using EdgeGrid1 = typename Interface::Edge1::parallel_grid;
    using EdgeGrid2 = typename Interface::Edge2::parallel_grid;

    using EdgeDim1 = typename EdgeGrid1::continuous_dimension_type;
    using EdgeDim2 = typename EdgeGrid2::continuous_dimension_type;

    using IdxRangeEdge1 = IdxRange<EdgeGrid1>;
    using IdxRangeEdge2 = IdxRange<EdgeGrid2>;

    using IdxEdge1 = Idx<EdgeGrid1>;
    using IdxEdge2 = Idx<EdgeGrid2>;

    using Patch1 = typename Interface::Edge1::associated_patch;
    using Patch2 = typename Interface::Edge2::associated_patch;

    IdxRangeEdge1 const& m_idx_range_patch_1;
    IdxRangeEdge2 const& m_idx_range_patch_2;


public:
    /**
     * @brief Instantiate an EdgeTransformation.
     * @param idx_range_patch_1 1D index range on the patch 1 of the interface.
     * @param idx_range_patch_2 1D index range on the patch 2 of the interface.
    */
    EdgeTransformation(
            IdxRangeEdge1 const& idx_range_patch_1,
            IdxRangeEdge2 const& idx_range_patch_2)
        : m_idx_range_patch_1(idx_range_patch_1)
        , m_idx_range_patch_2(idx_range_patch_2)
    {
    }

    ~EdgeTransformation() = default;


    /**
     * @brief Transform a coordinate on the edge in the dimension 
     * of the current patch to the analogous coordinate on the target patch. 
     * 
     * @param current_coord
     *      A coordinate on the edge of the current patch.
     * 
     * @tparam CurrentDim 
     *      The current continuous dimension of the given coordinate coord. 
     * 
     * @return The analogous coordinate on the target patch. 
     * @warning This operator is ill-defined when the two patches have the same
     * continuous dimension. 
    */
    template <class CurrentDim>
    Coord<std::conditional_t<std::is_same_v<CurrentDim, EdgeDim1>, EdgeDim2, EdgeDim1>> operator()(
            Coord<CurrentDim> const& current_coord) const
    {
        static_assert(
                !std::is_same_v<EdgeDim1, EdgeDim2>,
                "The two patches have the same conditinous dimension. "
                "Please specify the input Patch as template parameter "
                "using .template operator<CurrentPatch>()(current_coord).");
        return get_equivalent_target_coordinate<
                CurrentDim,
                Interface>(current_coord, m_idx_range_patch_1, m_idx_range_patch_2);
    }

    /**
     * @brief Transform a coordinate on the edge in the dimension 
     * of the current patch to the analogous coordinate on the target patch. 
     * 
     * @param current_coord
     *      A coordinate on the edge of the current patch.
     * 
     * @tparam CurrentPatch 
     *      The current patch of the given coordinate coord. 
     * 
     * @return The analogous coordinate on the target patch. 
    */
    template <class CurrentPatch>
    Coord<std::conditional_t<std::is_same_v<CurrentPatch, Patch1>, EdgeDim2, EdgeDim1>> operator()(
            Coord<std::conditional_t<
                    std::is_same_v<CurrentPatch, Patch1>,
                    EdgeDim1,
                    EdgeDim2>> const& current_coord) const
    {
        return get_equivalent_target_coordinate<
                CurrentPatch,
                Interface>(current_coord, m_idx_range_patch_1, m_idx_range_patch_2);
    }


    /**
     * @brief Transform an index on the edge in the dimension of the current patch 
     * to the analogous index on the target patch. 
     * 
     * If the grids are uniform, we can simplify the algorithm by using modulo. 
     * Otherwise, we need to check all the indices of the target grid. 
     * We suppose the coordinate transformation bijective, so we
     * can use a dichotomy method. 
     * 
     * This method mainly calls search_for_match.
     * 
     * @param current_idx
     *      A index on the edge of the current patch.
     * @tparam CurrentIdx 
     *      The current index type of the given coordinate index. 
     * 
     * @return The analogous index on the target patch. 
     */
    template <class CurrentIdx>
    auto operator()(CurrentIdx const& current_idx) const
    {
        using IdxTarget
                = std::conditional_t<std::is_same_v<CurrentIdx, IdxEdge1>, IdxEdge2, IdxEdge1>;
        IdxTarget target_idx;
        [[maybe_unused]] bool is_equivalent_idx_found = search_for_match(target_idx, current_idx);
        assert(is_equivalent_idx_found);
        return target_idx;
    }



    /**
     * @brief Check if a given index has an equivalent index on the other 
     * patch of an interface.
     * 
     * If the grids are uniform, we can simplify the algorithm by using modulo. 
     * Otherwise, we need to check all the indices of the target grid. 
     * We suppose the coordinate transformation bijective, so we
     * can use a dichotomy method.  
     * 
     * This method mainly calls search_for_match.
     * 
     * @param current_idx
     *      A index on the edge of the current patch.
     * @tparam CurrentIdx 
     *      The current index type of the given coordinate index. 
     * 
     * @return Boolean stating if there is an equivalent index. 
    */
    template <class CurrentIdx>
    bool is_match_available(CurrentIdx const& current_idx) const
    {
        using IdxTarget
                = std::conditional_t<std::is_same_v<CurrentIdx, IdxEdge1>, IdxEdge2, IdxEdge1>;
        IdxTarget potential_target_idx;
        return search_for_match(potential_target_idx, current_idx);
    }



    /**
     * @brief Check if a given index has an equivalent index and 
     * transform an index on the edge in the dimension of the current patch 
     * to the analogous index on the target patch. 
     * 
     * If the grids are uniform, we can simplify the algorithm by using modulo. 
     * Otherwise, we need to check all the indices of the target grid. 
     * We suppose the coordinate transformation bijective, so we
     * can use a dichotomy method. 
     * 
     * @warning target_idx is always replaced by the suspected index. 
     * If there is not equivalent index, the returned index is wrong but the 
     * closest that the algorithm found. 
     * 
     * @tparam CurrentGrid The grid where the input index is defined. 
     * @tparam TargetGrid The grid where the output index is defined.

     * @param[out] target_idx
     *      A index on the edge of the target patch.
     * @param[in] current_idx
     *      A index on the edge of the current patch.
     * @tparam CurrentIdx 
     *      The current index type of the given coordinate index. 
     * 
     * @return Boolean stating if there is an equivalent index. 
     * 
     */
    template <class CurrentGrid, class TargetGrid>
    bool search_for_match(Idx<TargetGrid>& target_idx, Idx<CurrentGrid> const& current_idx) const
    {
        static_assert(
                std::is_same_v<CurrentGrid, EdgeGrid1> || std::is_same_v<CurrentGrid, EdgeGrid2>,
                "Current index must be associated with one of the edges.");
        static_assert(
                std::is_same_v<TargetGrid, EdgeGrid1> || std::is_same_v<TargetGrid, EdgeGrid2>,
                "Target index must be associated with one of the edges.");
        static_assert(
                !std::is_same_v<TargetGrid, CurrentGrid>,
                "The types of the indices should be different");

        using TargetPatch
                = std::conditional_t<std::is_same_v<TargetGrid, EdgeGrid1>, Patch1, Patch2>;

        // Index range of CurrentIdx.
        using IdxRangeCurrent = IdxRange<CurrentGrid>;
        // Index range of other dimension.
        using IdxRangeTarget = IdxRange<TargetGrid>;
        // Coordinate on the discontinuous dimension of the current grid.
        using CurrentCoord = Coord<typename CurrentGrid::continuous_dimension_type>;
        // Index on the target grid
        using IdxTarget = typename IdxRangeTarget::discrete_element_type;
        // Index step on the target grid
        using IdxStepTarget = typename IdxRangeTarget::discrete_vector_type;


        bool is_equivalent_idx = false;

        // Get the 1D index range corresponding to the current and target domains.
        IdxRange<EdgeGrid1, EdgeGrid2> const
                combined_idx_range(m_idx_range_patch_1, m_idx_range_patch_2);
        IdxRangeCurrent const current_idx_range(combined_idx_range);
        IdxRangeTarget const target_idx_range(combined_idx_range);

        int n_cells_current
                = current_idx_range.size() - int(!CurrentGrid::continuous_dimension_type::PERIODIC);
        int n_cells_target
                = target_idx_range.size() - int(!TargetGrid::continuous_dimension_type::PERIODIC);

        if constexpr ( // Uniform case
                ddc::is_uniform_point_sampling_v<
                        EdgeGrid1> && ddc::is_uniform_point_sampling_v<EdgeGrid2>) {
            // Get greatest common divisor
            int const gcd_cells = std::gcd(n_cells_current, n_cells_target);

            int const current_idx_value = (current_idx - current_idx_range.front()).value();
            /*
                There is an equivalent index if the current index is a multiple of
                (number of current cells / gcd). If gcd == 1 then this condition is
                only true if the current index represents an extremity
                (current_idx_value == 0 or current_idx_value == n_cells_current).
            */
            is_equivalent_idx = (current_idx_value % (n_cells_current / gcd_cells) == 0);

            // If there is an equivalent index, update target_idx.
            if (is_equivalent_idx) {
                double const rescaling_factor = double(n_cells_current) / double(n_cells_target);
                IdxStepTarget target_idx_step_rescaled(int(current_idx_value / rescaling_factor));
                target_idx = target_idx_range.front() + target_idx_step_rescaled;

                if constexpr (!Interface::orientations_agree) {
                    IdxStepTarget target_idx_step = target_idx_range.back() - target_idx;
                    target_idx = IdxTarget(target_idx_range.front())
                                 + int(TargetGrid::continuous_dimension_type::PERIODIC)
                                 + target_idx_step;

                    if constexpr (TargetGrid::continuous_dimension_type::PERIODIC) {
                        // Apply periodicity property of the domain.
                        if (target_idx == target_idx_range.back() + 1) {
                            target_idx = target_idx_range.front();
                        }
                    }
                }
            }
        } else { // Non uniform case
            // Dichotomy method comparing the coordinates of indexes of the target edge.
            target_idx = get_target_idx(target_idx_range, current_idx);
            CurrentCoord const current_coord(ddc::coordinate(current_idx));
            CurrentCoord target_equivalent_coord
                    = (*this).template operator()<TargetPatch>(ddc::coordinate(target_idx));
            is_equivalent_idx = (abs(current_coord - target_equivalent_coord) < 1e-14);
        }

        return is_equivalent_idx;
    }



private:
    /// @brief Get index corresponding to the middle of two indexes.
    template <class IdxType>
    IdxType const get_mid(IdxType const& idx_1, IdxType const& idx_2) const
    {
        return idx_1 + (idx_2 - idx_1) / 2;
    }

    /// @brief Get the closest target index which is supposed to be the equivalent to the current index.
    template <class TargetGrid, class CurrentGrid>
    Idx<TargetGrid> const get_target_idx(
            IdxRange<TargetGrid> const& target_idx_range,
            Idx<CurrentGrid> const& current_idx) const
    {
        using TargetPatch
                = std::conditional_t<std::is_same_v<TargetGrid, EdgeGrid1>, Patch1, Patch2>;

        // Coordinate on the discontinuous dimension of the current grid.
        using CurrentCoord = Coord<typename CurrentGrid::continuous_dimension_type>;
        // Index on the target grid
        using IdxTarget = Idx<TargetGrid>;
        // Index step on the target grid
        using IdxStepTarget = IdxStep<TargetGrid>;

        // Dichotomy method
        CurrentCoord const current_coord(ddc::coordinate(current_idx));

        IdxTarget target_idx_min(target_idx_range.front());
        IdxTarget target_idx_max(
                target_idx_range.back() + int(TargetGrid::continuous_dimension_type::PERIODIC));
        IdxTarget target_idx_mid = get_mid(target_idx_min, target_idx_max);

        IdxStepTarget target_idx_step_diff = target_idx_max - target_idx_min;
        while (target_idx_step_diff != IdxStepTarget(1)) {
            CurrentCoord target_equivalent_coord_mid
                    = (*this).template operator()<TargetPatch>(ddc::coordinate(target_idx_mid));
            if constexpr (TargetGrid::continuous_dimension_type::PERIODIC) {
                if (target_idx_mid == IdxTarget(target_idx_range.back() + 1)) {
                    target_equivalent_coord_mid = (*this).template operator()<TargetPatch>(
                            ddc::coordinate(target_idx_range.front())
                            + ddcHelper::total_interval_length(target_idx_range));
                }
            }

            if ((current_coord > target_equivalent_coord_mid && Interface::orientations_agree)
                || (current_coord <= target_equivalent_coord_mid
                    && !Interface::orientations_agree)) {
                target_idx_min = target_idx_mid;
            } else {
                target_idx_max = target_idx_mid;
            }

            target_idx_mid = get_mid(target_idx_min, target_idx_max);
            target_idx_step_diff = target_idx_max - target_idx_min;
        }

        if constexpr (TargetGrid::continuous_dimension_type::PERIODIC) {
            if (target_idx_max == IdxTarget(target_idx_range.back() + 1)) {
                target_idx_max = target_idx_range.front();
            }
        }

        CurrentCoord target_coord_min
                = (*this).template operator()<TargetPatch>(ddc::coordinate(target_idx_min));
        CurrentCoord target_coord_max
                = (*this).template operator()<TargetPatch>(ddc::coordinate(target_idx_max));

        if (abs(current_coord - target_coord_min) < abs(current_coord - target_coord_max)) {
            return target_idx_min;
        } else {
            return target_idx_max;
        }
    }
};
