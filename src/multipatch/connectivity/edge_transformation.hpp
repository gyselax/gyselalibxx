// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <ddc_helper.hpp>

#include "interface.hpp"

/**
 * @brief Transform a coordinate or an index from one edge to the one on the other edge.
 * 
 * According to the orientation of the interface, we compute 
 * * if True,
 *  @f$ \\ t \mapsto min_2 + \frac{t - min_1}{max_1 - min_1}(max_2 - min_2) @f$
 * 
 * * if False, 
 *  @f$ \\ t \mapsto max_2 - \frac{t - min_1}{max_1 - min_1}(max_2 - min_2) @f$
 * 
 * where @f$ min_i @f$ and @f$ max_i @f$ are the minimum and maximum 
 * coordinates of the edge @f$ i @f$. 
 * 
 * @tparam Interface The Interface type where we want to compute the transformation.
*/
template <class Interface>
class EdgeTransformation
{
    using EdgeGrid1 = typename Interface::Edge1::grid;
    using EdgeGrid2 = typename Interface::Edge2::grid;

    using EdgeDim1 = typename EdgeGrid1::continuous_dimension_type;
    using EdgeDim2 = typename EdgeGrid2::continuous_dimension_type;

    using Patch1 = typename Interface::Patch1;
    using Patch2 = typename Interface::Patch2;

    using IdxRangeEdge1 = ddc::DiscreteDomain<EdgeGrid1>;
    using IdxRangeEdge2 = ddc::DiscreteDomain<EdgeGrid2>;

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
        , m_idx_range_patch_2(idx_range_patch_2) {};

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
    */
    template <class CurrentDim>
    ddc::Coordinate<
            std::conditional_t<std::is_same_v<CurrentDim, EdgeDim1>, EdgeDim2, EdgeDim1>> const
    operator()(ddc::Coordinate<CurrentDim> const& current_coord) const
    {
        // The other continuous dimension
        using ODim = std::conditional_t<std::is_same_v<CurrentDim, EdgeDim1>, EdgeDim2, EdgeDim1>;

        // The index range of CurrentDim
        using CurrentIdxRange = std::
                conditional_t<std::is_same_v<CurrentDim, EdgeDim1>, IdxRangeEdge1, IdxRangeEdge2>;
        // The index range of other dimension
        using OIdxRange = std::
                conditional_t<std::is_same_v<CurrentDim, EdgeDim1>, IdxRangeEdge2, IdxRangeEdge1>;

        // Get min and length on each 1D domains:
        CurrentIdxRange const current_idx_range(
                ddc::DiscreteDomain<
                        EdgeGrid1,
                        EdgeGrid2>(m_idx_range_patch_1, m_idx_range_patch_2));
        OIdxRange const target_idx_range(ddc::DiscreteDomain<
                                         EdgeGrid1,
                                         EdgeGrid2>(m_idx_range_patch_1, m_idx_range_patch_2));

        ddc::Coordinate<CurrentDim> const current_min = ddc::coordinate(current_idx_range.front());
        double const current_length = ddcHelper::total_interval_length(current_idx_range);

        ddc::Coordinate<ODim> const target_min = ddc::coordinate(target_idx_range.front());
        double const target_length = ddcHelper::total_interval_length(target_idx_range);

        double rescale_x = (current_coord - current_min) / current_length * target_length;

        bool constexpr orientations_agree = Interface::orientations_agree;
        if constexpr (!orientations_agree) {
            rescale_x = target_length - rescale_x;
        }
        return target_min + rescale_x; // This has type ddc::Coordinate<ODim>.
    };
};
