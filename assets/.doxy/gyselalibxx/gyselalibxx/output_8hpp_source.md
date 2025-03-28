

# File output.hpp

[**File List**](files.md) **>** [**io**](dir_c184e51c84f2c3f0345bbc8a0d75d3e1.md) **>** [**output.hpp**](output_8hpp.md)

[Go to the documentation of this file](output_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include "ddc_aliases.hpp"

template <class Mesh>
void expose_mesh_to_pdi(std::string pdi_name, IdxRange<Mesh> idx_range)
{
    using Dim = typename Mesh::continuous_dimension_type;
    using FieldMem = host_t<FieldMem<Coord<Dim>, IdxRange<Mesh>>>;
    using Idx = Idx<Mesh>;

    FieldMem mesh_coord(idx_range);
    for (Idx const idx : idx_range) {
        mesh_coord(idx) = ddc::coordinate(idx);
    }
    ddc::expose_to_pdi(pdi_name, mesh_coord);
}
```


