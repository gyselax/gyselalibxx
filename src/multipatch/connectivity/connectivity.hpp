// SPDX-License-Identifier: MIT
#pragma once

#include "connectivity_details.hpp"

template <class... Interfaces>
struct MultipatchConnectivity
{
    using interface_collection = ddc::detail::TypeSeq<Interfaces...>;

    template <class Patch>
    using get_type_seq_connections_t =
            typename connectivity_details::PatchConnection<Patch, interface_collection>::type;

    template <class Patch>
    using get_connections_t = connectivity_details::to_tuple_t<get_type_seq_connections_t<Patch>>;
};
