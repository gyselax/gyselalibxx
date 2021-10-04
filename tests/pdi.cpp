#include <array>
#include <string>
#include <string_view>

#include <ddc/MCoord>
#include <ddc/ProductMDomain>
#include <ddc/TaggedVector>
#include <ddc/expose_to_pdi>

#include <paraconf.h>
#include <pdi.h>

#include "geometry.h"

std::string_view path;

int main(int argc, char** argv)
{
    PC_tree_t conf;
    if (!path.empty()) {
        conf = PC_parse_path(path.data());
    } else {
        return 0;
    }

    PDI_init(conf);

    MeshX mesh_x(RCoordX(0.), RCoordX(2.));
    MeshVx mesh_vx(RCoordVx(0.), RCoordVx(2.));
    ProductMesh mesh_x_vx(mesh_x, mesh_vx);
    MDomainXVx const dom(mesh_x_vx, MCoordXVx(0, 0), MCoordXVx(9, 9));
    DBlockXVx fdistribu(dom);
    auto const fdistribu_s = fdistribu.cview();

    for (const auto& x : fdistribu.domain<MeshX>()) {
        for (const auto& vx : fdistribu.domain<MeshVx>()) {
            fdistribu(x, vx) = x + vx * 2;
        }
    }

    expose_to_pdi("fdistribu", fdistribu_s);

    // finalize PDI
    PDI_finalize();

    // destroy the paraconf configuration tree
    PC_tree_destroy(&conf);

    return 0;
}
