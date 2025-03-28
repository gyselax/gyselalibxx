

# File connectivity.hpp

[**File List**](files.md) **>** [**connectivity**](dir_28b51abc9241105ab41b66c468e7d019.md) **>** [**connectivity.hpp**](connectivity_8hpp.md)

[Go to the documentation of this file](connectivity_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "ddc/ddc.hpp"

#include "connectivity_details.hpp"
#include "ddc_aliases.hpp"
#include "multipatch_type.hpp"

template <class... Interfaces>
class MultipatchConnectivity
{
public:
    using interface_collection = ddc::detail::TypeSeq<Interfaces...>;

    using inner_edges
            = strip_outside_edges_t<typename Interfaces::Edge1..., typename Interfaces::Edge2...>;

    using all_patches = extract_patches_t<inner_edges>;

    static_assert(
            ddc::type_seq_size_v<inner_edges> == 4 * ddc::type_seq_size_v<all_patches>,
            "Missing edges. There should be 4 edges for each patch.");

    template <class QueryPatch>
    using get_type_seq_connections_t = interfaces_of_patch_t<QueryPatch, interface_collection>;

    template <class QueryPatch>
    using get_connections_t = to_tuple_t<get_type_seq_connections_t<QueryPatch>>;

    template <class QueryPatch1, class QueryPatch2>
    using find_connections_t = ddcHelper::type_seq_intersection_t<
            get_type_seq_connections_t<QueryPatch1>,
            get_type_seq_connections_t<QueryPatch2>>;

    template <class Grid1D>
    using get_all_interfaces_along_direction_t = collect_interfaces_on_dim_t<
            find_patch_t<Grid1D, all_patches>,
            Grid1D,
            interface_collection>;

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

    template <class Grid1D, template <typename P> typename T, class... Patches>
    static auto get_all_idx_ranges_along_direction(MultipatchType<T, Patches...> all_idx_ranges)
    {
        static_assert(ddc::type_seq_same_v<all_patches, ddc::detail::TypeSeq<Patches...>>);
        return get_all_idx_ranges_along_direction<Grid1D>(all_idx_ranges.get_tuple());
    }

private:
    template <class RelevantGrids, class IdxRangeTuple, std::size_t... PatchIndex>
    static auto get_idx_range(IdxRangeTuple all_idx_ranges, std::index_sequence<PatchIndex...>)
    {
        return std::make_tuple(
                get_idx_range<PatchIndex, RelevantGrids, IdxRangeTuple>(all_idx_ranges)...);
    }

    template <std::size_t PatchIndex, class RelevantGrids, class IdxRangeTuple>
    static auto get_idx_range(IdxRangeTuple all_idx_ranges)
    {
        using GridToLocate = ddc::type_seq_element_t<PatchIndex, RelevantGrids>;
        using GridLocation = find_relevant_idx_range_t<GridToLocate, IdxRangeTuple>;
        return ddc::select<GridToLocate>(std::get<GridLocation>(all_idx_ranges));
    }
};
```


