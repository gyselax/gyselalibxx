#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>

#include <ddc/MCoord>
#include <ddc/ProductMDomain>
#include <ddc/ProductMesh>
#include <ddc/pdi.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "geometry.h"
#include "nulladvectionvx.h"
#include "nullefieldsolver.h"
#include "pdi_out.yml.h"
#include "predcorr.h"
#include "splineadvectionx.h"
#include "splitvlasovsolver.h"

using std::cerr;
using std::cout;
using std::endl;
using std::exp;
namespace fs = std::filesystem;

void gaussian_initialization(DSpanXVx fdistribu)
{
    auto gridx = fdistribu.domain<MeshX>();
    auto gridvx = fdistribu.domain<MeshVx>();
    for (MCoordX ix : gridx) {
        for (MCoordVx iv : gridvx) {
            const RCoordVx v = gridvx.to_real(iv);
            fdistribu(ix, iv) = exp(-v * v / 2.);
        }
    }
}

int main(int argc, char** argv)
{
    PC_tree_t conf;
    if (argc > 1) {
        conf = PC_parse_path(fs::path(argv[1]).c_str());
    } else {
        cerr << "usage: " << argv[0] << " <config_file.yml>" << endl;
        return EXIT_FAILURE;
    }

    // Reading config

    long steps;
    PC_int(PC_get(conf, ".steps"), &steps);

    double time_step;
    PC_double(PC_get(conf, ".time_step"), &time_step);

    double mass_ratio;
    PC_double(PC_get(conf, ".mass_ratio"), &mass_ratio);

    RCoordX x_min = [&]() {
        double x_min;
        PC_double(PC_get(conf, ".MeshX.min"), &x_min);
        return x_min;
    }();
    RCoordX x_max = [&]() {
        double x_max;
        PC_double(PC_get(conf, ".MeshX.max"), &x_max);
        return x_max;
    }();
    MLengthX x_size = [&]() {
        long x_size;
        PC_int(PC_get(conf, ".MeshX.size"), &x_size);
        return x_size;
    }();

    RCoordVx vx_min = [&]() {
        double vx_min;
        PC_double(PC_get(conf, ".MeshVx.min"), &vx_min);
        return vx_min;
    }();
    RCoordVx vx_max = [&]() {
        double vx_max;
        PC_double(PC_get(conf, ".MeshVx.max"), &vx_max);
        return vx_max;
    }();
    MLengthVx vx_size = [&]() {
        double vx_size;
        PC_double(PC_get(conf, ".MeshVx.size"), &vx_size);
        return vx_size;
    }();

    PDI_init(PC_parse_string(PDI_CFG));

    // Creating mesh & supports

    BSplinesX const bsplines_x(x_min, x_max, x_size);

    SplineXBuilder const builder_x(bsplines_x);

    MeshVx const mesh_vx(vx_min, vx_max, vx_size);

    MeshXVx const mesh_x_vx(builder_x.interpolation_domain().mesh().get<MeshX>(), mesh_vx);

    MDomainXVx const
            dom2d(mesh_x_vx, MCoordXVx(builder_x.interpolation_domain().extents(), vx_size));

    // Creating operators

    SplineAdvectionX const advection_x(bsplines_x, builder_x);

    NullAdvectionVx const advection_vx;

    SplitVlasovSolver const vlasov(advection_x, advection_vx);

    NullEfieldSolver const efield;

    PredCorr const predcorr(vlasov, efield, time_step);

    // Creating data and initialization

    DBlockXVx fdistribu(dom2d);
    gaussian_initialization(fdistribu);

    BlockX<RCoordX> meshX_coord(fdistribu.domain<MeshX>());
    MDomainX gridx = fdistribu.domain<MeshX>();
    for (MCoordX ix : gridx) {
        meshX_coord(ix) = gridx.to_real(ix);
    }
    BlockVx<RCoordVx> meshVx_coord(fdistribu.domain<MeshVx>());
    MDomainVx gridvx = fdistribu.domain<MeshVx>();
    for (MCoordVx ivx : gridvx) {
        meshVx_coord(ivx) = gridvx.to_real(ivx);
    }

    // Starting the code

    expose_to_pdi("Nx", x_size);
    expose_to_pdi("Nvx", vx_size);
    expose_to_pdi("MeshX", meshX_coord);
    expose_to_pdi("MeshVx", meshVx_coord);

    predcorr(fdistribu, mass_ratio, steps);

    PDI_finalize();

    return EXIT_SUCCESS;
}
