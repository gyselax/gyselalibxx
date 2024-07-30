// SPDX-License-Identifier: MIT
#pragma once

#include "connectivity_details.hpp"

template <class... Interfaces>
struct MultipatchConnectivity
{
    using interface_collection = ddc::detail::TypeSeq<Interfaces...>;

    using inner_edges = filter_edges_t<typename Interfaces::Edge1..., typename Interfaces::Edge2...>;

    using all_patches = extract_patches_t<inner_edges>;

    template <class Patch>
    using get_type_seq_connections_t = direct_patch_connections_t<Patch, interface_collection>;

    template <class Patch>
    using get_connections_t = to_tuple_t<get_type_seq_connections_t<Patch>>;

    template <class Grid1D, class... Domains>
    auto get_all_idx_ranges_along_direction(std::tuple<Domains..>)
    {
        using StartPatch = find_patch_t<Grid1D>;
        using FrontEdge = Edge<StartPatch, Grid1D, FRONT>;
        using BackEdge = Edge<StartPatch, Grid1D, BACK>;
    }

    auto get_idx_range() {
    }
};
