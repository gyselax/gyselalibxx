// SPDX-License-Identifier: MIT
#pragma once

#include "connectivity_details.hpp"
#include "ddc_helper.hpp"

template <class... Interfaces>
class MultipatchConnectivity
{
public:
    using interface_collection = ddc::detail::TypeSeq<Interfaces...>;

    using inner_edges
            = filter_edges_t<typename Interfaces::Edge1..., typename Interfaces::Edge2...>;

    using all_patches = extract_patches_t<inner_edges>;

    static_assert(ddcHelper::type_seq_length_v<inner_edges> == 4*ddcHelper::type_seq_length_v<all_patches>,
            "Missing edges. There should be 4 edges for each patch.");

    template <class Patch>
    using get_type_seq_connections_t = direct_patch_connections_t<Patch, interface_collection>;

    template <class Patch>
    using get_connections_t = to_tuple_t<get_type_seq_connections_t<Patch>>;

    template <class Grid1D, class... Domains>
    static auto get_all_idx_ranges_along_direction(std::tuple<Domains...> all_domains)
    {
        using StartPatch = find_patch_t<Grid1D, all_patches>;
        using RelevantGrids = collect_grids_on_dim_t<StartPatch, Grid1D, interface_collection>;
        static_assert(ddcHelper::type_seq_length_v<RelevantGrids> > 0);
        return get_idx_range<RelevantGrids>(
                all_domains,
                std::make_index_sequence<ddcHelper::type_seq_length_v<RelevantGrids>>{});
    }

private:
    template <class RelevantGrids, class DomainTuple, std::size_t... PatchIndex>
    static auto get_idx_range(DomainTuple all_domains, std::index_sequence<PatchIndex...>)
    {
        return std::make_tuple(get_idx_range<PatchIndex, RelevantGrids, DomainTuple>(all_domains)...);
    }

    template <std::size_t PatchIndex, class RelevantGrids, class DomainTuple>
    static auto get_idx_range(DomainTuple all_domains)
    {
        using GridToLocate = ddc::type_seq_element_t<PatchIndex, RelevantGrids>;
        using GridLocation
                = find_relevant_idx_range_t<GridToLocate, DomainTuple>;
        return ddc::select<GridToLocate>(std::get<GridLocation>(all_domains));
    }
};
