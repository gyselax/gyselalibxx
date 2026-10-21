

# File coord\_or\_idx.hpp

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**utils**](dir_313caf1132e152dd9b58bea13a4052ca.md) **>** [**coord\_or\_idx.hpp**](coord__or__idx_8hpp.md)

[Go to the documentation of this file](coord__or__idx_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <type_traits>

#include <ddc/ddc.hpp>

namespace concepts {
template <class T, class CoordType>
concept CoordOrIdx = std::is_same_v<T, CoordType> ||(
        ddc::is_discrete_element_v<T>&& std::is_same_v<ddc::coordinate_of_t<T>, CoordType>);
} // namespace concepts

namespace ddcHelper {
template <class CoordType, concepts::CoordOrIdx<CoordType> CoordOrIdxType>
KOKKOS_INLINE_FUNCTION CoordType get_coordinate(CoordOrIdxType const& coord_or_idx)
{
    if constexpr (std::is_same_v<CoordOrIdxType, CoordType>) {
        return coord_or_idx;
    } else {
        return ddc::coordinate(coord_or_idx);
    }
}
} // namespace ddcHelper
```


