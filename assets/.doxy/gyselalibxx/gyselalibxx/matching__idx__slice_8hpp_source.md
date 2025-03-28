

# File matching\_idx\_slice.hpp

[**File List**](files.md) **>** [**connectivity**](dir_28b51abc9241105ab41b66c468e7d019.md) **>** [**matching\_idx\_slice.hpp**](matching__idx__slice_8hpp.md)

[Go to the documentation of this file](matching__idx__slice_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <numeric>
#include <stdexcept>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "edge.hpp"
#include "edge_transformation.hpp"
#include "idx_range_slice.hpp"

template <class Interface>
class MatchingIdxSlice
{
    static_assert(
            (!std::is_same_v<typename Interface::Edge1, OutsideEdge>)&&(
                    !std::is_same_v<typename Interface::Edge2, OutsideEdge>),
            "The interface cannot be an interface with an OutsideEdge.");

    using Patch1 = typename Interface::Edge1::associated_patch;
    using Patch2 = typename Interface::Edge2::associated_patch;

    using EdgeGrid1 = typename Interface::Edge1::parallel_grid;
    using EdgeGrid2 = typename Interface::Edge2::parallel_grid;

    using PerpEdgeGrid1 = typename Interface::Edge1::perpendicular_grid;
    using PerpEdgeGrid2 = typename Interface::Edge2::perpendicular_grid;

    using IdxRange1D_1 = IdxRange<EdgeGrid1>;
    using IdxRange1D_2 = IdxRange<EdgeGrid2>;

    using IdxRange2D_1 = typename Patch1::IdxRange12;
    using IdxRange2D_2 = typename Patch2::IdxRange12;

    using IdxRangeSlice1 = IdxRangeSlice<EdgeGrid1>;
    using IdxRangeSlice2 = IdxRangeSlice<EdgeGrid2>;

    static constexpr bool are_grids_uniform = (ddc::is_uniform_point_sampling_v<EdgeGrid1>)&&(
            ddc::is_uniform_point_sampling_v<EdgeGrid2>);

    IdxRange1D_1 const m_idx_range_edge_1;
    IdxRange1D_2 const m_idx_range_edge_2;

    IdxRangeSlice1 m_conforming_idx_1;
    IdxRangeSlice2 m_conforming_idx_2;

public:
    ~MatchingIdxSlice() = default;

    MatchingIdxSlice(IdxRange1D_1 const& idx_range_1, IdxRange1D_2 const& idx_range_2)
        : m_idx_range_edge_1(idx_range_1)
        , m_idx_range_edge_2(idx_range_2)
    {
        std::vector<Idx<EdgeGrid1>> conforming_indices_1;
        std::vector<Idx<EdgeGrid2>> conforming_indices_2;
        if constexpr (!are_grids_uniform) {
            get_conforming_idx_vector(conforming_indices_1, conforming_indices_2);
        }
        m_conforming_idx_1 = get_idx_range_slice(m_idx_range_edge_1, conforming_indices_1);
        m_conforming_idx_2 = get_idx_range_slice(m_idx_range_edge_2, conforming_indices_2);
    }


    MatchingIdxSlice(IdxRange2D_1 const& idx_range_1, IdxRange2D_2 const& idx_range_2)
        : MatchingIdxSlice(ddc::select<EdgeGrid1>(idx_range_1), ddc::select<EdgeGrid2>(idx_range_2))
    {
    }


    template <
            class ParallelGrid,
            std::enable_if_t<
                    (std::is_same_v<ParallelGrid, EdgeGrid1>)
                            || (std::is_same_v<ParallelGrid, EdgeGrid2>),
                    bool> = true>
    IdxRangeSlice<ParallelGrid> get() const
    {
        if constexpr (std::is_same_v<ParallelGrid, EdgeGrid1>) {
            return m_conforming_idx_1;
        } else {
            return m_conforming_idx_2;
        }
    }


    template <
            class PerpendicularGrid,
            std::enable_if_t<
                    (std::is_same_v<PerpendicularGrid, PerpEdgeGrid1>)
                            || (std::is_same_v<PerpendicularGrid, PerpEdgeGrid2>),
                    bool> = true>
    auto get_from_perp() const
    {
        if constexpr (std::is_same_v<PerpendicularGrid, PerpEdgeGrid1>) {
            return m_conforming_idx_1;
        } else {
            return m_conforming_idx_2;
        }
    }


private:
    void get_conforming_idx_vector(
            std::vector<Idx<EdgeGrid1>>& conforming_idx_vec_1,
            std::vector<Idx<EdgeGrid2>>& conforming_idx_vec_2)
    {
        EdgeTransformation<Interface> edge_transformation(m_idx_range_edge_1, m_idx_range_edge_2);
        Idx<EdgeGrid2> idx_2;
        ddc::for_each(m_idx_range_edge_1, [&](Idx<EdgeGrid1> const& idx_1) {
            if (edge_transformation.search_for_match(idx_2, idx_1)) {
                conforming_idx_vec_1.push_back(idx_1);
                conforming_idx_vec_2.push_back(idx_2);
            }
        });
    }


    template <class Grid1D>
    IdxStep<Grid1D> get_idx_step(std::vector<Idx<Grid1D>> const& conforming_idx_vec) const
    {
        if constexpr (are_grids_uniform) {
            std::size_t const n_cells_1 = m_idx_range_edge_1.size() - 1;
            std::size_t const n_cells_2 = m_idx_range_edge_2.size() - 1;

            // Get greatest common divisor
            int const gcd_cells = std::gcd(n_cells_1, n_cells_2);
            // Select right number of cells
            std::size_t const n_cells = std::is_same_v<Grid1D, EdgeGrid1> ? n_cells_1 : n_cells_2;
            int const step = n_cells / gcd_cells;
            return IdxStep<Grid1D>(step);
        } else {
            IdxStep<Grid1D> const first_idx_step = conforming_idx_vec[1] - conforming_idx_vec[0];
            for (std::size_t i(2); i < conforming_idx_vec.size(); i++) {
                Idx<Grid1D> const idx = conforming_idx_vec[i];
                Idx<Grid1D> const previous_idx = conforming_idx_vec[i - 1];

                IdxStep<Grid1D> idx_step = idx - previous_idx;
                if (idx_step != first_idx_step) {
                    throw std::invalid_argument(
                            "The steps between conforming indexes have to be uniform.");
                }
            };
            return first_idx_step;
        }
    }


    template <class Grid1D>
    IdxRangeSlice<Grid1D> get_idx_range_slice(
            IdxRange<Grid1D> const& idx_range,
            std::vector<Idx<Grid1D>> const& conforming_idx_vec) const
    {
        Idx<Grid1D> idx_front = idx_range.front();
        IdxStep<Grid1D> stride = get_idx_step(conforming_idx_vec);
        IdxStep<Grid1D> size((idx_range.size() - 1) / stride.value() + 1);
        return IdxRangeSlice<Grid1D>(idx_front, size, stride);
    }
};
```


