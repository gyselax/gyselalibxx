

# File interface.hpp

[**File List**](files.md) **>** [**connectivity**](dir_28b51abc9241105ab41b66c468e7d019.md) **>** [**interface.hpp**](interface_8hpp.md)

[Go to the documentation of this file](interface_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "edge.hpp"

template <class Edge1Type, class Edge2Type, bool orientations_agree_bool>
struct Interface
{
    static_assert(
            !std::is_same_v<Edge1Type, OutsideEdge> || !std::is_same_v<Edge2Type, OutsideEdge>,
            "At least one edge should belong to a patch.");

    using Edge1 = Edge1Type;
    using Edge2 = Edge2Type;

    template <
            class Edge,
            class = std::enable_if_t<std::is_same_v<Edge, Edge1> || std::is_same_v<Edge, Edge2>>>
    using OtherEdge = std::conditional_t<std::is_same_v<Edge, Edge1>, Edge2, Edge1>;

    static constexpr bool orientations_agree = orientations_agree_bool;

    template <class Patch>
    static constexpr bool connected_to_patch()
    {
        if constexpr (std::is_same_v<Edge1, OutsideEdge>) {
            return std::is_same_v<typename Edge2::associated_patch, Patch>;
        } else if constexpr (std::is_same_v<Edge2, OutsideEdge>) {
            return std::is_same_v<typename Edge1::associated_patch, Patch>;
        } else {
            return (std::is_same_v<typename Edge1::associated_patch, Patch>)
                   || (std::is_same_v<typename Edge2::associated_patch, Patch>);
        }
    }
};
```


