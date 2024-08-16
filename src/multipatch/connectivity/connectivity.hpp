// SPDX-License-Identifier: MIT
#pragma once

#include "ddc/ddc.hpp"

#include "connectivity_details.hpp"
#include "ddc_aliases.hpp"
#include "multipatch_type.hpp"

/**
 * @brief A helper class which provides functionalities to recognise how different patches are
 * connected.
 *
 * @tparam Interfaces The interfaces which describe the connectivity in this geometry.
 */
template <class... Interfaces>
class MultipatchConnectivity
{
public:
    /// A type sequence of all the interfaces handled by this class.
    using interface_collection = ddc::detail::TypeSeq<Interfaces...>;

    /// A type sequence of all the edges handled by this class except the OuterEdge types.
    using inner_edges
            = strip_outside_edges_t<typename Interfaces::Edge1..., typename Interfaces::Edge2...>;

    /// A type sequence of all the patches handled by this class.
    using all_patches = extract_patches_t<inner_edges>;

    static_assert(
            ddc::type_seq_size_v<inner_edges> == 4 * ddc::type_seq_size_v<all_patches>,
            "Missing edges. There should be 4 edges for each patch.");

    /**
     * @brief A tool to find a type sequence of all patches directly connected (via an interface)
     * to the patch.
     *
     * @tparam QueryPatch The patch whose connections we are interested in.
     */
    template <class QueryPatch>
    using get_type_seq_connections_t = direct_patch_connections_t<QueryPatch, interface_collection>;

    /**
     * @brief A tool to find a tuple of all patches directly connected (via an interface) to
     * the patch.
     *
     * @tparam QueryPatch The patch whose connections we are interested in.
     */
    template <class QueryPatch>
    using get_connections_t = to_tuple_t<get_type_seq_connections_t<QueryPatch>>;

    /**
     * @brief A function to return all index ranges which can be used to obtain coordinates
     * along a line which passes through the requested grid.
     *
     * @tparam Grid1D The grid indicating the direction of interest.
     * @param all_idx_ranges A tuple containing all available index ranges.
     *
     * @return A tuple of index ranges along the line of interest.
     */
    template <class Grid1D, class... IdxRanges>
    static auto get_all_idx_ranges_along_direction(std::tuple<IdxRanges...> all_idx_ranges)
    {
        using StartPatch = find_patch_t<Grid1D, all_patches>;
        using RelevantGrids = collect_grids_on_dim_t<StartPatch, Grid1D, interface_collection>;
        static_assert(ddc::type_seq_size_v<RelevantGrids> > 0);
        return get_idx_range<RelevantGrids>(
                all_idx_ranges,
                std::make_index_sequence<ddc::type_seq_size_v<RelevantGrids>> {});
    }

    /**
     * @brief A function to return all index ranges which can be used to obtain coordinates
     * along a line which passes through the requested grid.
     *
     * @tparam Grid1D The grid indicating the direction of interest.
     * @param all_idx_ranges A MultipatchType containing all available index ranges.
     *
     * @return A tuple of index ranges along the line of interest.
     */
    template <class Grid1D, template <typename P> typename T, class... Patches>
    static auto get_all_idx_ranges_along_direction(MultipatchType<T, Patches...> all_idx_ranges)
    {
        static_assert(ddc::type_seq_same_v<all_patches, ddc::detail::TypeSeq<Patches...>>);
        return get_all_idx_ranges_along_direction<Grid1D>(all_idx_ranges.get_tuple());
    }

private:
    /// A compile-time loop to get all the relevant index ranges identified by the associated grid.
    template <class RelevantGrids, class IdxRangeTuple, std::size_t... PatchIndex>
    static auto get_idx_range(IdxRangeTuple all_idx_ranges, std::index_sequence<PatchIndex...>)
    {
        return std::make_tuple(
                get_idx_range<PatchIndex, RelevantGrids, IdxRangeTuple>(all_idx_ranges)...);
    }

    /// Get the relevant index range from the input index ranges using the grid to identify it.
    template <std::size_t PatchIndex, class RelevantGrids, class IdxRangeTuple>
    static auto get_idx_range(IdxRangeTuple all_idx_ranges)
    {
        using GridToLocate = ddc::type_seq_element_t<PatchIndex, RelevantGrids>;
        using GridLocation = find_relevant_idx_range_t<GridToLocate, IdxRangeTuple>;
        return ddc::select<GridToLocate>(std::get<GridLocation>(all_idx_ranges));
    }
};
