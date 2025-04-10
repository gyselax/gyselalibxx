// SPDX-License-Identifier: MIT

#pragma once

#include <numeric>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "edge.hpp"
#include "interface.hpp"


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
        using CurrentPatch
                = std::conditional_t<std::is_same_v<CurrentDim, EdgeDim1>, Patch1, Patch2>;
        return transform_edge_coord<CurrentPatch>(current_coord);
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
    Coord<std::conditional_t<std::is_same_v<CurrentPatch, Patch1>, EdgeDim2, EdgeDim1>>
    transform_edge_coord(Coord<std::conditional_t<
                                 std::is_same_v<CurrentPatch, Patch1>,
                                 EdgeDim1,
                                 EdgeDim2>> const& current_coord) const
    {
        using CurrentDim
                = std::conditional_t<std::is_same_v<CurrentPatch, Patch1>, EdgeDim1, EdgeDim2>;
        using TargetDim
                = std::conditional_t<std::is_same_v<CurrentPatch, Patch1>, EdgeDim2, EdgeDim1>;

        using IdxRangeCurrent = IdxRange<
                std::conditional_t<std::is_same_v<CurrentPatch, Patch1>, EdgeGrid1, EdgeGrid2>>;
        using IdxRangeTarget = IdxRange<
                std::conditional_t<std::is_same_v<CurrentPatch, Patch1>, EdgeGrid2, EdgeGrid1>>;

        // Get min and length on each 1D index range:
        IdxRange<EdgeGrid1, EdgeGrid2> const
                combined_idx_range(m_idx_range_patch_1, m_idx_range_patch_2);
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
     * @brief Transform a coordinate of an index on the edge of the current patch 
     * to the analogous coordinate on the target patch. 
     * 
     * @param current_idx
     *      A index on the edge of the current patch.
     * @tparam CurrentIdx 
     *      The current index type of the given coordinate index. 
     * 
     * @return The analogous coordinate on the target patch. 
    */
    template <class CurrentIdx>
    Coord<std::conditional_t<std::is_same_v<CurrentIdx, IdxEdge1>, EdgeDim2, EdgeDim1>>
    get_equivalent_coord_from_idx(CurrentIdx const& current_idx) const
    {
        using CurrentPatch
                = std::conditional_t<std::is_same_v<CurrentIdx, IdxEdge1>, Patch1, Patch2>;
        return transform_edge_coord<CurrentPatch>(ddc::coordinate(current_idx));
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
    bool search_for_match(Idx<TargetGrid>& target_idx, Idx<CurrentGrid> current_idx) const
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

        using TargetPatch = std::conditional<std::is_same_v<TargetGrid,EdgeGrid1>, Patch1, Patch2>; 

        // Index range of CurrentIdx.
        using IdxRangeCurrent = IdxRange<CurrentGrid>;
        // Index range of other dimension.
        using IdxRangeTarget = IdxRange<TargetGrid>;
        // Coordinate on the discontinuous dimension of the current grid.
        using CurrentCoord = Coord<typename CurrentGrid::continuous_dimension_type>;
        // Index step on the target grid
        using IdxStepTarget = typename IdxRangeTarget::discrete_vector_type;

        /*
            Add a point for the periodic case where the last interpolation point
            is usually not on the last break point. 
        */
        IdxRangeEdge1 const idx_range_patch_1_full(
                m_idx_range_patch_1.front(),
                m_idx_range_patch_1.extents() + int(EdgeDim1::PERIODIC));
        IdxRangeEdge2 const idx_range_patch_2_full(
                m_idx_range_patch_2.front(),
                m_idx_range_patch_2.extents() + int(EdgeDim2::PERIODIC));


        bool is_equivalent_idx = false;

        // Get the 1D index range corresponding to the current and target domains.
        IdxRange<EdgeGrid1, EdgeGrid2> const
                combined_idx_range(idx_range_patch_1_full, idx_range_patch_2_full);
        IdxRangeCurrent const current_idx_range(combined_idx_range);
        IdxRangeTarget const target_idx_range(combined_idx_range);

        int const n_cells_current = current_idx_range.size() - 1;
        int const n_cells_target = target_idx_range.size() - 1;

        // Periodicity property.
        if constexpr (CurrentGrid::continuous_dimension_type::PERIODIC) {
            current_idx = (current_idx == current_idx_range.back()) ? current_idx_range.front()
                                                                    : current_idx;
        }

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
                    target_idx = target_idx_range.front() + target_idx_step;
                }

                // Apply periodicity property of the domain.
                if constexpr (TargetGrid::continuous_dimension_type::PERIODIC) {
                    if (target_idx == target_idx_range.back()) {
                        target_idx = target_idx_range.front();
                    }
                }
            }
        } else { // Non uniform case
            // Dichotomy method comparing the coordinates of indexes of the target edge.
            target_idx = get_target_idx(target_idx_range, current_idx);
            CurrentCoord const current_coord(ddc::coordinate(current_idx));
            CurrentCoord target_equivalent_coord = transform_edge_coord<TargetPatch>(ddc::coordinate(target_idx));
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
        using TargetPatch = std::conditional<std::is_same_v<TargetGrid,EdgeGrid1>, Patch1, Patch2>; 

        // Coordinate on the discontinuous dimension of the current grid.
        using CurrentCoord = Coord<typename CurrentGrid::continuous_dimension_type>;
        // Index on the target grid
        using IdxTarget = Idx<TargetGrid>;
        // Index step on the target grid
        using IdxStepTarget = IdxStep<TargetGrid>;

        // Dichotomy method
        CurrentCoord const current_coord(ddc::coordinate(current_idx));

        IdxTarget target_idx_min(target_idx_range.front());
        IdxTarget target_idx_max(target_idx_range.back());
        IdxTarget target_idx_mid = get_mid(target_idx_min, target_idx_max);

        IdxStepTarget target_idx_step_diff = target_idx_max - target_idx_min;
        while (target_idx_step_diff != IdxStepTarget(1)) {
            CurrentCoord target_equivalent_coord_mid = transform_edge_coord<TargetPatch>(ddc::coordinate(target_idx_mid));

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

        // Periodicity property.
        if constexpr (TargetGrid::continuous_dimension_type::PERIODIC) {
            if (target_idx_max == target_idx_range.back()) {
                target_idx_max = target_idx_range.front();
            }
        }

        CurrentCoord target_coord_min = transform_edge_coord<TargetPatch>(ddc::coordinate(target_idx_min));
        CurrentCoord target_coord_max = transform_edge_coord<TargetPatch>(ddc::coordinate(target_idx_max));

        if (abs(current_coord - target_coord_min) < abs(current_coord - target_coord_max)) {
            return target_idx_min;
        } else {
            return target_idx_max;
        }
    }
};
