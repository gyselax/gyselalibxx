#include <cstdlib>
#include <iostream>

#include <ddc/MCoord>
#include <ddc/ProductMDomain>
#include <ddc/ProductMesh>

#include <paraconf.h>

#include "geometry.h"
#include "nulladvectionvx.h"
#include "nullefieldsolver.h"
#include "predcorr.h"
#include "splineadvectionx.h"
#include "splitvlasovsolver.h"

using std::cerr;
using std::endl;

int main(int argc, char** argv)
{
    PC_tree_t conf;
    if (argc > 1) {
        conf = PC_parse_path(argv[1]);
    } else {
        cerr << "usage: " << argv[0] << " <config_file.yml>" << endl;
        return EXIT_FAILURE;
    }

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
    MLengthElement x_size = [&]() {
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
    MLengthElement vx_size = [&]() {
        double vx_size;
        PC_double(PC_get(conf, ".MeshVx.size"), &vx_size);
        return vx_size;
    }();

    BSplinesX bsplines_x(x_min, x_max, x_size);

    SplineXBuilder builder_x(bsplines_x);

    MeshVx const mesh_vx(vx_min, vx_max, vx_size);

    MeshXVx const mesh_x_vx(builder_x.interpolation_domain().mesh().get<MeshX>(), mesh_vx);

    MDomainXVx const
            dom2d(mesh_x_vx,
                  MCoordXVx(builder_x.interpolation_domain().extents() - 1, vx_size - 1));

    SplineAdvectionX const advection_x(bsplines_x, builder_x);

    NullAdvectionVx const advection_vx;

    SplitVlasovSolver const vlasov(advection_x, advection_vx);

    NullEfieldSolver const efield;

    PredCorr const predcorr(vlasov, efield, time_step);

    DBlockXVx fdistribu(dom2d);

    predcorr(fdistribu, mass_ratio, steps);

    return EXIT_SUCCESS;
}
