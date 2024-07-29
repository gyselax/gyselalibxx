// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"

/**
 * @brief Expose a IdxRange to PDI.
 *
 * @param pdi_name The name given in PDI.
 * @param dom IdxRange that is exposed
 *
 */
template <class Mesh>
void expose_mesh_to_pdi(std::string pdi_name, IdxRange<Mesh> dom)
{
    using Dim = typename Mesh::continuous_dimension_type;
    using FieldMem = host_t<FieldMem<Coord<Dim>, IdxRange<Mesh>>>;
    using Idx = Idx<Mesh>;

    FieldMem mesh_coord(dom);
    for (Idx const idx : dom) {
        mesh_coord(idx) = ddc::coordinate(idx);
    }
    ddc::expose_to_pdi(pdi_name, mesh_coord);
}
