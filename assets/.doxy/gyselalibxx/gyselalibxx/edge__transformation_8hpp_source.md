

# File edge\_transformation.hpp

[**File List**](files.md) **>** [**connectivity**](dir_28b51abc9241105ab41b66c468e7d019.md) **>** [**edge\_transformation.hpp**](edge__transformation_8hpp.md)

[Go to the documentation of this file](edge__transformation_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <numeric>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "edge.hpp"
#include "interface.hpp"


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

    IdxRangeEdge1 const& m_idx_range_patch_1;
    IdxRangeEdge2 const& m_idx_range_patch_2;


public:
    EdgeTransformation(
            IdxRangeEdge1 const& idx_range_patch_1,
            IdxRangeEdge2 const& idx_range_patch_2)
        : m_idx_range_patch_1(idx_range_patch_1)
        , m_idx_range_patch_2(idx_range_patch_2)
    {
    }

    ~EdgeTransformation() = default;


    template <class CurrentDim>
    Coord<std::conditional_t<std::is_same_v<CurrentDim, EdgeDim1>, EdgeDim2, EdgeDim1>> operator()(
            Coord<CurrentDim> const& current_coord) const
    {
        // The other continuous dimension
        using TargetDim
                = std::conditional_t<std::is_same_v<CurrentDim, EdgeDim1>, EdgeDim2, EdgeDim1>;

        // The index range of CurrentDim
        using IdxRangeCurrent = std::
                conditional_t<std::is_same_v<CurrentDim, EdgeDim1>, IdxRangeEdge1, IdxRangeEdge2>;
        // The index range of other dimension
        using IdxRangeTarget = std::
                conditional_t<std::is_same_v<CurrentDim, EdgeDim1>, IdxRangeEdge2, IdxRangeEdge1>;

        // Get min and length on each 1D index ranges:
        IdxRange<EdgeGrid1, EdgeGrid2> const
                combined_idx_range(m_idx_range_patch_1, m_idx_range_patch_2);
        IdxRangeCurrent const current_idx_range(combined_idx_range);
        IdxRangeTarget const target_idx_range(combined_idx_range);

        Coord<CurrentDim> const current_min = ddc::coordinate(current_idx_range.front());
        double const current_length = ddcHelper::total_interval_length(current_idx_range);

        Coord<TargetDim> const target_min = ddc::coordinate(target_idx_range.front());
        double const target_length = ddcHelper::total_interval_length(target_idx_range);

        double rescale_x = (current_coord - current_min) / current_length * target_length;

        bool constexpr orientations_agree = Interface::orientations_agree;
        if constexpr (!orientations_agree) {
            rescale_x = target_length - rescale_x;
        }
        return target_min + rescale_x;
    }



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



    template <class CurrentIdx>
    bool is_match_available(CurrentIdx const& current_idx) const
    {
        using IdxTarget
                = std::conditional_t<std::is_same_v<CurrentIdx, IdxEdge1>, IdxEdge2, IdxEdge1>;
        IdxTarget potential_target_idx;
        return search_for_match(potential_target_idx, current_idx);
    }



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

        int const n_cells_current = current_idx_range.size() - 1;
        int const n_cells_target = target_idx_range.size() - 1;

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
                target_idx = IdxTarget(target_idx_range.front());
                target_idx += target_idx_step_rescaled;

                if constexpr (!Interface::orientations_agree) {
                    IdxStepTarget target_idx_step = IdxTarget(target_idx_range.back()) - target_idx;
                    target_idx = IdxTarget(target_idx_range.front());
                    target_idx += target_idx_step;
                }
            }
        } else { // Non uniform case
            // Dichotomy method comparing the coordinates of indexes of the target edge.
            target_idx = get_target_idx(target_idx_range, current_idx);
            CurrentCoord const current_coord(ddc::coordinate(current_idx));
            CurrentCoord const target_equivalent_coord((*this)(ddc::coordinate(target_idx)));
            is_equivalent_idx = (abs(current_coord - target_equivalent_coord) < 1e-14);
        }

        return is_equivalent_idx;
    }



private:
    template <class IdxType>
    IdxType const get_mid(IdxType const& idx_1, IdxType const& idx_2) const
    {
        return idx_1 + (idx_2 - idx_1) / 2;
    }

    template <class TargetGrid, class CurrentGrid>
    Idx<TargetGrid> const get_target_idx(
            IdxRange<TargetGrid> const& target_idx_range,
            Idx<CurrentGrid> const& current_idx) const
    {
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
            CurrentCoord target_equivalent_coord_mid((*this)(ddc::coordinate(target_idx_mid)));

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

        CurrentCoord target_coord_min((*this)(ddc::coordinate(target_idx_min)));
        CurrentCoord target_coord_max((*this)(ddc::coordinate(target_idx_max)));
        if (abs(current_coord - target_coord_min) < abs(current_coord - target_coord_max)) {
            return target_idx_min;
        } else {
            return target_idx_max;
        }
    }
};
```


