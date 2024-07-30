// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

/**
 * @brief Expose a ddc::DiscreteDomain to PDI.
 *
 * @param pdi_name The name given in PDI.
 * @param dom ddc::DiscreteDomain that is exposed
 *
 */
template <class Mesh>
void expose_mesh_to_pdi(std::string pdi_name, ddc::DiscreteDomain<Mesh> dom)
{
    using Dim = typename Mesh::continuous_dimension_type;
    using Field = ddc::Chunk<ddc::Coordinate<Dim>, ddc::DiscreteDomain<Mesh>>;
    using Index = ddc::DiscreteElement<Mesh>;

    Field mesh_coord(dom);
    for (Index const idx : dom) {
        mesh_coord(idx) = ddc::coordinate(idx);
    }
    ddc::expose_to_pdi(pdi_name, mesh_coord);
}
